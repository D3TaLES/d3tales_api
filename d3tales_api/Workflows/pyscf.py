"""
PySCF counterpart of Gaussian.py (d3tales_api.Workflows).

Every FireTask in Gaussian.py has a PySCF twin here with the same task keywords, the same FWAction update_spec
keys, and the same working-directory / file-name conventions:

    Gaussian.py                  Psyscf.py
    ---------------------------  ------------------------------------------------
    GaussianBase                 PySCFBase
    RunGaussianEnergy            RunPySCFEnergy
    RunGaussianOpt               RunPySCFOpt            (opt -> freq -> re-opt loop, max 5 runs)
    RunWtuning                   RunPySCFWtuning        (uses PySCFWtuningJob, a WtuningJob subclass)
    RunGaussianTDDFT             RunPySCFTDDFT          (+ NTO analysis when run_nto is set)

What changes because the engine is PySCF (there is no Gaussian input/output text to write or parse):
    * <prefix>.com        -> <prefix>.xyz       (geometry input)
    * <prefix>.chk        -> <prefix>.chk       (PySCF HDF5 checkpoint)
    * <prefix>.fchk       -> <prefix>.molden    (what `formchk` produced; PySCF has no fchk writer)
    * <prefix>.log        -> <prefix>.log       (PySCF log + a Gaussian-style summary block at the end so the
                                                 existing helpers in utils.py -- parse_tddft, get_tddft_excitations,
                                                 runtime_from_log -- keep working on it)
    * Gaussian route keywords are translated: IOp(3/107,3/108) -> omega, SCRF=(PCM,Solvent=) -> IEF-PCM,
      SCF=(xqc,vshift,conver) -> SCF ladder, Opt -> geomeTRIC, TD(NStates,50-50) -> PySCF TDDFT.
    * Results are stored as local JSON records (see utils.save_local_result). PySCF output cannot be
      uploaded to the D3TaLES website (its parser reads Gaussian logs), so nothing is ever uploaded;
      the D3TaLES database is still used for READING (ground state, SMILES, geometry, omega).

Requires: pyscf, geometric (geomeTRIC; for geometry optimization), scipy (omega tuning).
"""
import abc
import multiprocessing
import time
import uuid
import warnings

import numpy as np
from six import add_metaclass
from atomate.utils.utils import get_logger, env_chk
from fireworks import FiretaskBase, explicit_serialize, FWAction
from pyscf import gto, scf, dft, lib
import pyscf

from d3tales_api.Workflows.utils import *
from d3tales_api.Workflows.wtuning import WtuningJob

logger = get_logger(__name__)
cpus = [multiprocessing.cpu_count() if multiprocessing.cpu_count() < 16 else 16]
nprocs = str(cpus[0])

# Copyright 2021, University of Kentucky

HARTREE_TO_EV = 27.211396132  # same constant wtuning.omega_extract uses
PYSCF_VERBOSE = 4
PYSCF_WTUNING_KEY = "wtuning_pyscf"  # local-record key for PySCF-tuned omega (kept apart from Gaussian's "wtuning")
DEFAULT_MAX_MEMORY_MB = 48000  # Gaussian.py asks for %mem=48GB

# ----------------------------------------------------------------------------------------------------------------------
# Gaussian name -> PySCF/libxc translation tables
# ----------------------------------------------------------------------------------------------------------------------
# Gaussian functional keyword (lower case, without a U/R/RO prefix) -> PySCF xc string.
# NOTE: Gaussian's LC-wHPBE has no identical libxc twin; LC_WPBE (long-range-corrected wPBE, 100% long-range HF) is the
# closest, and its omega is overridden by the tuned value exactly as IOp(3/107,3/108) does in Gaussian.
FUNCTIONAL_MAP = {
    'hf': 'HF',
    'b3lyp': 'B3LYPG',  # B3LYPG = Gaussian's B3LYP (VWN3); PySCF's plain 'B3LYP' uses VWN5
    'pbe0': 'PBE0', 'pbe1pbe': 'PBE0',
    'pbe': 'PBE', 'pbepbe': 'PBE',
    'blyp': 'BLYP',
    'bp86': 'B88,P86',
    'b3pw91': 'B3PW91',
    'tpss': 'TPSS', 'tpsstpss': 'TPSS',
    'm06': 'M06', 'm062x': 'M062X', 'm06-2x': 'M062X',
    'cam-b3lyp': 'CAM_B3LYP', 'camb3lyp': 'CAM_B3LYP',
    'wb97xd': 'WB97X_D', 'wb97x-d': 'WB97X_D',
    'hse06': 'HSE06', 'hseh1pbe': 'HSE06',
    'lc-wpbe': 'LC_WPBE', 'lc-whpbe': 'LC_WPBE', 'lcwpbe': 'LC_WPBE', 'lcwhpbe': 'LC_WPBE',
}

# Dielectric constants (Gaussian's values) for SCRF=(PCM,Solvent=...)
SOLVENT_EPS = {
    'water': 78.3553, 'acetonitrile': 35.688, 'dichloromethane': 8.93, 'methylenechloride': 8.93,
    'tetrahydrofuran': 7.4257, 'thf': 7.4257, 'toluene': 2.3741, 'dimethylsulfoxide': 46.826, 'dmso': 46.826,
    'methanol': 32.613, 'ethanol': 24.852, 'chloroform': 4.7113, 'acetone': 20.493, 'benzene': 2.2706,
    'n,n-dimethylformamide': 37.219, 'dimethylformamide': 37.219, 'dmf': 37.219, 'n-hexane': 1.8819,
    'hexane': 1.8819, 'diethylether': 4.24, 'ether': 4.24,
}


def split_functional(func):
    """
    Split a Gaussian functional keyword into (PySCF xc string, reference prefix).
    Prefix is 'U', 'R', 'RO' or None. The outermost prefix wins ("R" + "uLC-wHPBE" -> 'R', as RunWtuning builds it).
    """
    name = str(func).strip()
    low, prefix = name.lower(), None
    while low not in FUNCTIONAL_MAP:
        for p in ('ro', 'u', 'r'):
            rest = low[len(p):]
            if low.startswith(p) and rest and (rest in FUNCTIONAL_MAP or rest[0] in 'ur'):
                prefix, low = prefix or p.upper(), rest
                break
        else:
            return name, prefix  # unknown keyword: hand it to PySCF/libxc untouched
    return FUNCTIONAL_MAP[low], prefix


def to_pyscf_basis(basis):
    """Gaussian basis keyword -> PySCF basis name ('6-31G(d,p)' -> '6-31g**', 'def2tzvp' stays)."""
    return str(basis).strip().lower().replace(' ', '').replace('(d,p)', '**').replace('(d)', '*')


def iop_to_omega(iop_str):
    """Inverse of utils.omega_to_iop: '0123400000' -> 0.1234 (the value IOp(3/107,3/108) carries)."""
    return int(str(iop_str)) / 1e9


def _route_text(route_parameters, prefix):
    """Flatten every route keyword starting with `prefix` (keys + values, dict values included) into one lower string."""
    chunks = []
    for key, val in (route_parameters or {}).items():
        if str(key).lower().startswith(prefix):
            if isinstance(val, dict):
                val = ' '.join(str(k) if v in ('', None) else '{}={}'.format(k, v) for k, v in val.items())
            chunks.append('{} {}'.format(key, val or ''))
    return ' '.join(chunks).lower()


def parse_scf_options(route_parameters):
    """SCF=(xqc,fermi,noincfock,ndamp=35,conver=6,vshift=500,novaracc) -> PySCF SCF ladder settings."""
    txt = _route_text(route_parameters, 'scf')
    conver = re.search(r'conver\s*=\s*(\d+)', txt)
    vshift = re.search(r'vshift\s*=\s*(\d+)', txt)
    conver = int(conver.group(1)) if conver else 8
    return {
        'conv_tol': 10.0 ** -(conver + 1),
        'conv_tol_grad': 10.0 ** -(max(conver - 1, 3)),
        'level_shift': int(vshift.group(1)) / 1000.0 if vshift else 0.3,  # Gaussian VShift is in milliHartree
        'xqc': 'noxqc' not in txt,  # quadratically convergent fallback (default in Gaussian)
    }


def parse_opt_options(route_parameters):
    txt = _route_text(route_parameters, 'opt')
    maxc = re.search(r'maxcycles?\s*=\s*(\d+)', txt)
    conv_set = 'GAU'
    for word, name in (('verytight', 'GAU_VERYTIGHT'), ('tight', 'GAU_TIGHT'), ('loose', 'GAU_LOOSE')):
        if word in txt:
            conv_set = name
            break
    return {'max_cycles': int(maxc.group(1)) if maxc else 100, 'convergence_set': conv_set}


def parse_td_options(route_parameters):
    """TD(NStates=n, 50-50 | Singlets | Triplets, TDA) -> number of singlet / triplet roots."""
    txt = _route_text(route_parameters, 'td')
    n = re.search(r'nstates?\s*=\s*(\d+)', txt)
    n = int(n.group(1)) if n else 3  # Gaussian default
    if '50-50' in txt:  # half singlets, half triplets (NStates in total)
        n_singlet = (n + 1) // 2
        n_triplet = n - n_singlet
    elif 'triplet' in txt:
        n_singlet, n_triplet = 0, n
    else:
        n_singlet, n_triplet = n, 0
    return {'n_singlet': n_singlet, 'n_triplet': n_triplet, 'tda': bool(re.search(r'tamm|\btda\b', txt))}


# ----------------------------------------------------------------------------------------------------------------------
# Hash ids / conditions (PySCF versions of utils.orig_conditions / orig_hash_id so results never collide with Gaussian's)
# ----------------------------------------------------------------------------------------------------------------------
def pyscf_conditions(functional, basis_set, tuning_parameter=None, solvent=None):
    conditions = orig_conditions(functional, basis_set, tuning_parameter=tuning_parameter, solvent=solvent)
    conditions.update({"code_name": "PySCF", "code_version": pyscf.__version__.split('.')[0]})
    return conditions


def pyscf_hash_id(_id, calculation_type, functional, basis_set, tuning_parameter=None, solvent=None):
    hash_dict = {"_id": _id, "calculation_type": calculation_type,
                 "conditions": pyscf_conditions(functional, basis_set, tuning_parameter, solvent)}
    return hashlib.md5(json.dumps(hash_dict, sort_keys=True).encode()).hexdigest()


# ----------------------------------------------------------------------------------------------------------------------
# PySCF engine
# ----------------------------------------------------------------------------------------------------------------------
class PySCFSettings:
    """Everything needed to turn a (structure, charge, multiplicity) into a converged PySCF calculation."""

    def __init__(self, functional, basis, omega=None, solvent=None, solvent_eps=None, route_parameters=None,
                 density_fit=False, xc=None, max_memory=DEFAULT_MAX_MEMORY_MB, nthreads=None):
        parsed_xc, self.prefix = split_functional(functional)
        self.functional = functional
        self.xc = xc or parsed_xc
        self.basis = to_pyscf_basis(basis)
        self.omega = omega
        self.solvent = solvent
        if solvent:
            self.solvent_eps = solvent_eps or SOLVENT_EPS.get(str(solvent).lower())
            if not self.solvent_eps:
                raise ValueError("No dielectric constant known for solvent '{}'. Pass `solvent_eps` to the "
                                 "task or add it to SOLVENT_EPS.".format(solvent))
        self.scf_opts = parse_scf_options(route_parameters)
        self.opt_opts = parse_opt_options(route_parameters)
        self.td_opts = parse_td_options(route_parameters)
        self.density_fit = density_fit
        self.ladder_stage = 0
        self.max_memory = max_memory
        lib.num_threads(int(nthreads or nprocs))

    @property
    def label(self):
        return self.functional

    def build_mol(self, structure, charge, multiplicity, output=None):
        mol = gto.Mole()
        mol.atom = [(site.specie.symbol, tuple(site.coords)) for site in structure]
        mol.unit = 'Angstrom'
        mol.basis = self.basis
        if self.basis.startswith('def2'):  # PySCF does not attach def2 ECPs on its own (needed beyond Kr)
            heavy = {site.specie.symbol for site in structure if site.specie.Z > 36}
            if heavy:
                mol.ecp = {sym: self.basis for sym in heavy}
        mol.charge = int(charge)
        mol.spin = int(multiplicity) - 1
        mol.max_memory = self.max_memory
        mol.verbose = PYSCF_VERBOSE
        mol.symmetry = False
        if output:
            mol.output = output
        mol.build()
        return mol

    def build_mf(self, mol, chkfile=None):
        pre = self.prefix
        hf_only = self.xc.upper() == 'HF'
        if mol.spin == 0 and pre != 'U':
            mf = scf.RHF(mol) if hf_only else dft.RKS(mol)
        elif pre in ('R', 'RO'):
            mf = scf.ROHF(mol) if hf_only else dft.ROKS(mol)
        else:
            mf = scf.UHF(mol) if hf_only else dft.UKS(mol)
        if not hf_only:
            mf.xc = self.xc
            mf.grids.atom_grid = (99, 590)  # Gaussian's default UltraFine grid
            mf.grids.prune = None
            if self.omega is not None:
                mf.omega = self.omega  # == IOp(3/107, 3/108)
        if self.density_fit:
            mf = mf.density_fit()
        if self.solvent:  # == SCRF=(PCM,Solvent=...)  (Gaussian's PCM is IEF-PCM)
            mf = mf.PCM()
            mf.with_solvent.method = 'IEF-PCM'
            mf.with_solvent.eps = self.solvent_eps
        mf.conv_tol = self.scf_opts['conv_tol']
        mf.conv_tol_grad = self.scf_opts['conv_tol_grad']
        mf.max_cycle = 128
        if chkfile:
            mf.chkfile = chkfile
        return mf

    def run_scf(self, mf, dm0=None):
        """Plain DIIS -> level shift (VShift) -> second-order (XQC). Returns the converged mean-field object."""
        self.ladder_stage = 0
        mf.kernel(dm0=dm0)
        if not mf.converged:
            self.ladder_stage = 1
            mf.level_shift = self.scf_opts['level_shift']
            mf.kernel(dm0=mf.make_rdm1())
            mf.level_shift = 0
        if not mf.converged and self.scf_opts['xqc']:
            self.ladder_stage = 2
            mf2 = mf.newton()
            mf2.kernel(mf.mo_coeff, mf.mo_occ)
            mf = mf2
        if not mf.converged:
            raise RuntimeError("PySCF SCF did not converge ({})".format(self.label))
        return mf


class PySCFResult(dict):
    """Light container for what a finished PySCF stage produced (kept JSON friendly)."""
    __getattr__ = dict.get


def structure_from_mol(mol):
    return Molecule([mol.atom_symbol(i) for i in range(mol.natm)], mol.atom_coords(unit='Angstrom'))


def sites_of(structure):
    return structure.as_dict()['sites']


def frontier_orbitals(mf):
    """HOMO / LUMO in eV as wtuning.omega_extract defines them: HOMO = max(alpha, beta), LUMO = min(alpha, beta)."""
    mo_e, occ = np.asarray(mf.mo_energy), np.asarray(mf.mo_occ)
    if mo_e.ndim == 2:  # unrestricted
        channels = [(mo_e[0], occ[0] > 0), (mo_e[1], occ[1] > 0)]
    else:  # restricted / restricted-open
        channels = [(mo_e, occ >= 1), (mo_e, occ >= 2)]
    homos = [e[o].max() for e, o in channels if o.any()]
    lumos = [e[~o].min() for e, o in channels if (~o).any()]
    return max(homos) * HARTREE_TO_EV, min(lumos) * HARTREE_TO_EV


def append_gaussian_style_summary(logfile, mol, energy, label, cycles, elapsed, excitations=None):
    """
    Append a block in Gaussian's log format so utils.get_tddft_excitations / runtime_from_log / get_scfs-style parsers
    can read PySCF logs. `excitations` = [(kind, energy_eV, f, s2), ...] with kind 'Singlet' | 'Triplet' | '<S2>-A'.
    """
    if mol.stdout is not None and hasattr(mol.stdout, 'flush'):
        mol.stdout.flush()
    days, rem = divmod(elapsed, 86400)
    hours, rem = divmod(rem, 3600)
    minutes, seconds = divmod(rem, 60)
    lines = ["\n", " SCF Done:  E({}) = {:19.10f}     A.U. after {:4d} cycles\n".format(label, energy, int(cycles or 0))]
    if excitations:
        lines.append(" Excitation energies and oscillator strengths:\n")
        for i, (kind, e_ev, f, s2) in enumerate(excitations, start=1):
            lines.append(" Excited State {:3d}:      {}-A     {:8.4f} eV {:9.2f} nm  f={:.4f}  <S**2>={:.3f}\n".format(
                i, kind, e_ev, 1239.84198 / e_ev, f, s2))
    lines.append(" Job cpu time:  {:d} days {:2d} hours {:2d} minutes {:4.1f} seconds.\n".format(
        int(days), int(hours), int(minutes), seconds))
    lines.append(" Normal termination of PySCF job.\n")
    with open(logfile, 'a') as fn:
        fn.writelines(lines)


def write_xyz(filename, structure, comment=""):
    with open(filename, 'w') as fn:
        fn.write("{}\n{}\n".format(len(structure), comment))
        for site in structure:
            fn.write("{:<3s} {:16.10f} {:16.10f} {:16.10f}\n".format(site.specie.symbol, *site.coords))
    return filename


def write_molden(mf, filename):
    """What `formchk` produced for Gaussian: a portable orbital file. Raises like the original on failure."""
    from pyscf.tools import molden
    try:
        molden.from_scf(mf, filename)
    except Exception as err:
        raise RuntimeError("molden file {} could not be created: {}".format(filename, err))
    return filename


def run_energy(settings, structure, charge, multiplicity, output=None, chkfile=None, dm0=None):
    """Single point. Returns (PySCFResult, converged mean-field object)."""
    start = time.time()
    mol = settings.build_mol(structure, charge, multiplicity, output=output)
    mf = settings.run_scf(settings.build_mf(mol, chkfile=chkfile), dm0=dm0)
    homo, lumo = frontier_orbitals(mf)
    res = PySCFResult(properly_terminated=True, energy=float(mf.e_tot), energy_ev=float(mf.e_tot) * HARTREE_TO_EV,
                      homo_ev=homo, lumo_ev=lumo, sites=sites_of(structure_from_mol(mol)), elapsed=0.0)
    res['elapsed'] = time.time() - start
    if output:
        append_gaussian_style_summary(output, mol, mf.e_tot, settings.label, getattr(mf, 'cycles', 0), res['elapsed'])
    return res, mf


def run_optimization(settings, structure, charge, multiplicity, output=None, chkfile=None):
    """
    Geometry optimization with geomeTRIC (Gaussian-like convergence criteria), followed by a single point at the
    optimized geometry. `properly_terminated` is False when the optimizer ran out of steps (Gaussian: no normal
    termination); `sites` then holds the last geometry so the caller can restart from it.
    """
    from pyscf.geomopt import geometric_solver
    start = time.time()
    mol = settings.build_mol(structure, charge, multiplicity, output=output)
    # converge the starting point with the SCF ladder, then hand a plain mean-field object (with those orbitals and,
    # if the ladder needed it, the level shift) to geomeTRIC, which re-solves the SCF at every step
    probe = settings.run_scf(settings.build_mf(mol, chkfile=chkfile))
    mf = settings.build_mf(mol, chkfile=chkfile)
    mf.mo_coeff, mf.mo_occ, mf.mo_energy = probe.mo_coeff, probe.mo_occ, probe.mo_energy
    if settings.ladder_stage >= 1:
        mf.level_shift = settings.scf_opts['level_shift']
    conv, mol_eq = geometric_solver.kernel(
        mf, assert_convergence=False, maxsteps=settings.opt_opts['max_cycles'],
        convergence_set=settings.opt_opts['convergence_set'])
    final_structure = structure_from_mol(mol_eq)
    if not conv:
        return PySCFResult(properly_terminated=False, sites=sites_of(final_structure)), None
    res, mf_final = run_energy(settings, final_structure, charge, multiplicity, output=output, chkfile=chkfile)
    res['elapsed'] = time.time() - start
    res['converged_opt'] = True
    return res, mf_final


def run_frequencies(settings, structure, charge, multiplicity, output=None, chkfile=None, guess_chk=None):
    """
    Harmonic frequencies + thermochemistry (298.15 K, 1 atm) at `structure`. Returns (PySCFResult, mf).
    `frequencies` are in cm^-1 with imaginary modes reported as negative numbers (as Gaussian prints them);
    `gibbs_correction` is the thermal correction to the Gibbs free energy in Hartree (G - E_elec).
    """
    from pyscf.hessian import thermo
    start = time.time()
    mol = settings.build_mol(structure, charge, multiplicity, output=output)
    mf = settings.build_mf(mol, chkfile=chkfile)
    dm0 = None
    if guess_chk and os.path.isfile(guess_chk):  # Guess=TCheck: start from the optimization orbitals
        try:
            dm0 = mf.from_chk(guess_chk)
        except Exception:
            dm0 = None
    mf = settings.run_scf(mf, dm0=dm0)
    hess = mf.Hessian().kernel()
    freq_info = thermo.harmonic_analysis(mol, hess)
    wavenumbers = np.asarray(freq_info['freq_wavenumber'])
    freqs = np.where(np.iscomplex(wavenumbers), -np.abs(wavenumbers.imag), wavenumbers.real)
    freqs = np.sort(freqs)
    thermo_data = thermo.thermo(mf, freq_info['freq_au'], 298.15, 101325)
    e_elec = float(mf.e_tot)
    gibbs_corr = float(thermo_data['G_tot'][0]) - e_elec
    res = PySCFResult(properly_terminated=True, energy=e_elec, frequencies=[float(f) for f in freqs],
                      gibbs_correction=gibbs_corr, zpe=float(thermo_data['ZPE'][0]),
                      sites=sites_of(structure_from_mol(mol)), elapsed=time.time() - start)
    if output:
        append_gaussian_style_summary(output, mol, mf.e_tot, settings.label, getattr(mf, 'cycles', 0), res['elapsed'])
    return res, mf


def run_tddft(settings, structure, charge, multiplicity, output=None, chkfile=None):
    """
    Ground-state SCF + TDDFT. Returns (PySCFResult, mf, [td objects]). Singlet and/or triplet roots are requested
    from the TD(NStates=..., 50-50) route keyword (closed-shell only for triplets, as in Gaussian).
    """
    start = time.time()
    opts = settings.td_opts
    mol = settings.build_mol(structure, charge, multiplicity, output=output)
    mf = settings.build_mf(mol, chkfile=chkfile)
    if isinstance(mf, (scf.rohf.ROHF,)):
        raise NotImplementedError("TDDFT is not available for restricted open-shell (ROKS) references.")
    mf = settings.run_scf(mf)
    closed = np.asarray(mf.mo_energy).ndim == 1
    tds, excitations = [], []
    jobs = []
    if closed:
        jobs += [('Singlet', True, opts['n_singlet']), ('Triplet', False, opts['n_triplet'])]
    else:  # unrestricted: no singlet/triplet distinction
        jobs += [('{:.3f}'.format(multiplicity), None, opts['n_singlet'] + opts['n_triplet'])]
    for kind, singlet, nstates in jobs:
        if not nstates:
            continue
        td = mf.TDA() if opts['tda'] else mf.TDDFT()
        td.nstates = nstates
        if singlet is not None:
            td.singlet = singlet
        td.kernel()
        if not np.all(td.converged):
            raise RuntimeError("PySCF TDDFT ({}) did not converge".format(kind))
        f = td.oscillator_strength() if singlet is not False else np.zeros(len(td.e))
        s2 = {'Singlet': 0.0, 'Triplet': 2.0}.get(kind, (multiplicity ** 2 - 1) / 4.0)
        for e, osc in zip(td.e, f):
            excitations.append((kind, float(e) * HARTREE_TO_EV, float(osc), s2))
        tds.append((kind, td))
    res = PySCFResult(properly_terminated=True, energy=float(mf.e_tot), excitations=excitations,
                      sites=sites_of(structure_from_mol(mol)), elapsed=time.time() - start)
    if output:
        append_gaussian_style_summary(output, mol, mf.e_tot, settings.label, getattr(mf, 'cycles', 0), res['elapsed'],
                                      excitations=excitations)
    return res, mf, tds


def write_nto(td, mf, filename_stem, state=1, weight_cutoff=1e-3):
    """
    NTO analysis of excited state `state` (1-based) -- what write_nto/Pop=NTO did in Gaussian. Transition density is
    X+Y; SVD gives hole/particle NTOs. Written as Molden file(s) (occupation = NTO weight).
    :return: list of molden file paths (one per spin channel)
    """
    from pyscf.tools import molden
    x, y = td.xy[state - 1]
    mo_coeff, mo_occ = mf.mo_coeff, mf.mo_occ
    if isinstance(x, (tuple, list)):  # unrestricted
        channels = [('alpha', mo_coeff[0], mo_occ[0], x[0] + y[0]), ('beta', mo_coeff[1], mo_occ[1], x[1] + y[1])]
    else:
        channels = [('', mo_coeff, mo_occ, x + y)]
    files = []
    for spin, coeff, occ, t in channels:
        u, s, vt = np.linalg.svd(t)
        weights = s ** 2 / np.sum(s ** 2)
        keep = weights > weight_cutoff
        hole = coeff[:, occ > 0] @ u[:, :len(s)][:, keep]
        particle = coeff[:, occ == 0] @ vt[:len(s)][keep].T
        nto_occ = np.concatenate([weights[keep], weights[keep]])
        name = "{}{}_nto.molden".format(filename_stem, "_" + spin if spin else "")
        molden.from_mo(mf.mol, name, np.hstack([hole, particle]), occ=nto_occ)
        files.append(name)
    return files


# ----------------------------------------------------------------------------------------------------------------------
# FireTasks
# ----------------------------------------------------------------------------------------------------------------------
@add_metaclass(abc.ABCMeta)
class PySCFBase(FiretaskBase):
    _fw_name = "PySCFBase"

    def setup_files(self, fw_spec, calc_type='opt'):
        # get parameters
        self.data_store = get_data_store(self, fw_spec)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.full_name = name_tag + self['name']
        self.calc_name = self['name']
        self.paramset = self["paramset"]
        self.identifier = fw_spec.get("identifier", ) or self.get("identifier", uuid.uuid4().__str__())
        self.smiles = fw_spec.get("smiles", ) or self.get("smiles")
        self.runfile_log = env_chk(self.get('runfile_log'), fw_spec)
        self.skip_freq = fw_spec.get("skip_freq", ) or self.get("skip_freq") or False
        self.skip_freq = True if self.calc_name == "opt_mol" else self.skip_freq
        self.check_if_already_run = spec_or_task(self, fw_spec, "check_if_already_run", True)
        self.submit = spec_or_task(self, fw_spec, "submit", True)  # kept for API parity; PySCF results are never uploaded
        self.restricted = spec_or_task(self, fw_spec, "restricted", True)
        self.run_nto = fw_spec.get("run_nto") or self.get("run_nto") or False
        self.density_fit = spec_or_task(self, fw_spec, "density_fit", False)

        self.pyscf_file_name = fw_spec.get("pyscf_file_name") or self.get("pyscf_file_name") or "pyscf"

        # set prefix for the file names
        if self.get("prefix"):
            prefix = name_tag + self["prefix"]
        elif self.get("subtype", ):
            prefix = name_tag + calc_type + "_" + self["subtype"]
        else:
            prefix = name_tag + calc_type + "_mol"

        # get calc directory
        if self.get("path", ):
            self.calc_dir = "{}/{}/{}".format(env_chk(self.get('path'), fw_spec), self.identifier,
                                              self.pyscf_file_name)
        else:
            self.calc_dir = self.get("calc_dir", os.getcwd())

        # folder for local JSON records
        self.store_dir = get_store_dir(self.identifier, get_local_output_dir(self, fw_spec), self.calc_dir)

        # get the type of pyscf calculation
        if self.get("type", ):
            self.type = self["type"]
            if self.type == "coupling":
                self.working_dir = '{}/{}/{}/{}/{}'.format(self.calc_dir, self.type, calc_type, self["label"],
                                                           self["subtype"])
            elif self.get("label"):
                self.working_dir = '{}/{}/{}/{}'.format(self.calc_dir, self.type, calc_type, self["label"])
            elif self.type.startswith("solv_"):
                self.working_dir = '{}/{}/{}/{}'.format(self.calc_dir, self.type, calc_type, self["subtype"])
            else:
                self.working_dir = '{}/{}/{}'.format(self.calc_dir, self.type, calc_type)
        else:
            self.type = 'mol_{}'.format(calc_type)
            self.working_dir = "{}/{}".format(self.calc_dir, self.type)

        # Set default files (.com -> .xyz, .fchk -> .molden; see module docstring)
        self.file_xyz = "{}/{}.xyz".format(self.working_dir, prefix)
        self.file_chk = "{}/{}.chk".format(self.working_dir, prefix)
        self.file_molden = "{}/{}.molden".format(self.working_dir, prefix)
        self.file_log = "{}/{}.log".format(self.working_dir, prefix)
        self.file_result = "{}/{}_result.json".format(self.working_dir, prefix)  # lets a restart reuse finished stages
        self.freq_xyz = "{}/freq_{}.xyz".format(self.working_dir, prefix)
        self.freq_chk = "{}/freq_{}.chk".format(self.working_dir, prefix)
        self.freq_molden = "{}/freq_{}.molden".format(self.working_dir, prefix)
        self.freq_log = "{}/freq_{}.log".format(self.working_dir, prefix)
        self.freq_result = "{}/freq_{}_result.json".format(self.working_dir, prefix)

        # create folders
        os.makedirs(self.working_dir, exist_ok=True)
        os.chdir(self.working_dir)

    def setup_calc(self, fw_spec, calc_type='opt'):
        self.setup_files(fw_spec, calc_type=calc_type)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.solvent = fw_spec.get("solvent", ) or self.get("solvent", )
        self.gs_charge = fw_spec.get("gs_charge") or self.get("gs_charge") or get_groundState(
            self.identifier, self.smiles, data_store=self.data_store) or 0
        self.gs_spin = fw_spec.get("gs_spin") or self.get("gs_spin") or get_groundState(
            self.identifier, self.smiles, prop='spin', data_store=self.data_store) or 1
        self.use_iop = fw_spec.get("use_iop", True) if self.get("use_iop", True) else self.get("use_iop")
        run_from_com = fw_spec.get("run_from_com") or self.get("run_from_com") or \
                       fw_spec.get("run_from_xyz") or self.get("run_from_xyz") or False

        # get iop string (omega value) for the calculation: PySCF-tuned omega first, then the data store's
        self.iop_str = fw_spec.get("iop_str", ) or self.get("iop_str", )
        if not self.iop_str and self.use_iop:
            tuned_w = get_pyscf_omega(self.identifier, data_store=self.data_store, store_dir=self.store_dir)
            if tuned_w:
                self.iop_str = omega_to_iop(tuned_w)
        self.result_iop_str = self.iop_str  # fixed copy used for local record hashes

        # get starting geometry for calculation
        geometry = fw_spec.get("geometry", ) or self.get("geometry", )
        geometry_sites = fw_spec.get("{}_geom".format(geometry), )
        if run_from_com:
            self.mol = Molecule.from_file(self.file_xyz)
        else:
            if geometry_sites:
                self.mol = Molecule.from_sites([Site.from_dict(sd) for sd in geometry_sites])
            else:
                geometry_hash = fw_spec.get("{}_hash".format(geometry),
                                            pyscf_hash_id(self.identifier, self.calc_name, self.paramset.functional,
                                                          self.paramset.basis_set, tuning_parameter=self.iop_str,
                                                          solvent=self.solvent))
                self.mol = get_db_geom(geometry_hash, data_store="local", store_dir=self.store_dir)
                if self.mol is None and self.data_store == "d3tales":
                    self.mol = get_db_geom(geometry_hash, data_store=self.data_store, store_dir=self.store_dir)
                self.mol = self.mol or start_from_smiles(self.identifier, self.smiles, data_store=self.data_store)

            # End job if the total number of atoms is greater than 200
            num_atoms = len(self.mol.sites)
            if num_atoms > 200:
                return FWAction(defuse_children=True,
                                update_spec={"defuser_reason": "The molecule has {} atoms".format(num_atoms)})

        # charge / multiplicity bookkeeping (identical to Gaussian.py)
        radical_electrons = (self.gs_spin - 1 - self.paramset.charge) % 2
        self.paramset.multiplicity = radical_electrons + 1  # calculate spin multiplicity with Hand's rule
        self.paramset.charge += self.gs_charge
        # Update functional for unrestricted calculation
        self.functional = self.paramset.functional if self.restricted else "R{}".format(self.paramset.functional)
        print("Charge: ", self.paramset.charge, "\t Multiplicity: ", self.paramset.multiplicity)

        # translate Gaussian route options (omega IOp, PCM solvent, SCF/Opt/TD keywords) into PySCF settings
        omega = iop_to_omega(self.iop_str) if (self.iop_str and self.use_iop) else None
        self.settings = PySCFSettings(
            self.functional, self.paramset.basis_set, omega=omega, solvent=self.solvent,
            solvent_eps=self.get("solvent_eps"), route_parameters=self.paramset.route_parameters,
            density_fit=self.density_fit, xc=self.get("xc"))
        self.charge, self.multiplicity = self.paramset.charge, self.paramset.multiplicity

    def local_conditions(self):
        return pyscf_conditions(self.paramset.functional, self.paramset.basis_set,
                                tuning_parameter=self.result_iop_str, solvent=self.solvent)

    def local_key(self, calc_name=None):
        # key for local records; same hash scheme as the D3TaLES backend, but with PySCF conditions
        return pyscf_hash_id(self.identifier, calc_name or self.calc_name, self.paramset.functional,
                             self.paramset.basis_set, tuning_parameter=self.result_iop_str, solvent=self.solvent)

    def existing_data(self, _hash=None, calc_name=None):
        # check if this job has already run and been stored (PySCF results only live in the local store)
        return get_existing_data(_hash or self.local_key(calc_name), data_store="local", store_dir=self.store_dir)

    def post_job(self, data=None, upload_files=None, delete_files=None, calc_type=None, name_tag=''):
        # Store data: a local JSON record (PySCF output is not uploadable to the D3TaLES website)
        delete_files = []  # delete_files or [self.file_chk, self.file_xyz, self.file_molden, self.file_log]
        calc_type = calc_type or self.calc_name
        _hash = self.local_key(calc_type) if os.path.isfile(self.file_log) else None
        data = dict(data or {})
        data.setdefault("conditions", self.local_conditions())
        data.setdefault("mol_info", {"groundState_charge": self.gs_charge})
        record = local_record(self.local_key(calc_type), self.identifier, calc_type, data)
        print("Record saved to", save_local_result(self.store_dir, record, key=self.local_key(calc_type)))
        # Write runfile to runfile_log so the runfile can be deleted after calculation
        if self.runfile_log:
            with open(self.runfile_log, 'a') as fn:
                fn.write("{}, {}\n".format(_hash, self.working_dir))
        # TODO write script that checks process status, approves, and deletes run dir
        # Remove excess files
        os.chdir(self.working_dir)
        if delete_files:
            for file in delete_files:
                os.system("rm -fr {}".format(file))

    @staticmethod
    def energy_data(res):
        """The `data` payload shared by energy / opt records (mirrors the fields Gaussian.py reads back)."""
        # `value` is in Hartree so a cached hit returns the same units as a fresh run (Gaussian.py's fresh runs return
        # pymatgen's final_energy, in Hartree)
        return {"scf_total_energy": {"value": res.get("energy"), "unit": "Eh"},
                "scf_total_energy_ev": res.get("energy_ev"), "geometry": res.get("sites"),
                "homo": res.get("homo_ev"), "lumo": res.get("lumo_ev"),
                "timings": {"pyscf_seconds": res.get("elapsed")}}

    @staticmethod
    def save_stage(path, res):
        with open(path, 'w') as fn:
            json.dump(dict(res), fn)

    @staticmethod
    def load_stage(path, *needed_files):
        """Reuse a finished stage on restart (the Gaussian code skipped re-running when log + chk already existed)."""
        if os.path.isfile(path) and all(os.path.isfile(f) for f in needed_files):
            with open(path) as fn:
                return PySCFResult(json.load(fn))
        return None


def get_pyscf_omega(identifier, data_store="d3tales", store_dir=None):
    """Tuned omega: a PySCF-tuned local record wins; otherwise whatever the data store has (as utils.get_omega)."""
    return read_local_result(store_dir, PYSCF_WTUNING_KEY).get("omega") or get_omega(
        identifier, data_store=data_store, store_dir=store_dir)


@explicit_serialize
class RunPySCFEnergy(PySCFBase):

    def run_task(self, fw_spec):
        setup_obj = self.setup_calc(fw_spec, calc_type='energy')
        if isinstance(setup_obj, FWAction):
            return setup_obj

        # generate the input for the pyscf energy run
        write_xyz(self.file_xyz, self.mol)
        if self.check_if_already_run:
            existing_data = self.existing_data()
            if existing_data:
                return FWAction(
                    update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier,
                                 "gs_charge": self.gs_charge, "gs_spin": self.gs_spin,
                                 "{}_eng".format(self.full_name): existing_data.get("scf_total_energy", {}).get(
                                     "value"),
                                 "{}_hash".format(self.full_name): self.local_key(),
                                 "iop_str": self.iop_str})

        # run pyscf
        print('SUBMITTING PYSCF ENERGY JOB {} for {}'.format(self.full_name, self.identifier))
        res, mf = run_energy(self.settings, self.mol, self.charge, self.multiplicity,
                             output=self.file_log, chkfile=self.file_chk)
        final_energy = res.energy

        # check for normal termination of the job
        if not res.properly_terminated:
            raise RuntimeError(self.file_xyz + " not terminated. calc_dir " + self.working_dir)
        else:
            write_molden(mf, self.file_molden)  # == formchk

            # Clean up files and store results
            self.post_job(data=self.energy_data(res),
                          delete_files=[self.file_log, self.file_chk, self.file_molden, self.file_xyz])
        return FWAction(
            update_spec={"gaussrun_dir": self.working_dir, "identifier": self.identifier, "iop_str": self.iop_str,
                         "gs_charge": self.gs_charge, "gs_spin": self.gs_spin,
                         "{}_eng".format(self.full_name): final_energy})


@explicit_serialize
class RunPySCFOpt(PySCFBase):

    def run_task(self, fw_spec):
        setup_obj = self.setup_calc(fw_spec, calc_type='opt')
        if isinstance(setup_obj, FWAction):
            return setup_obj

        # run pyscf optimize until normal termination is achieved. max runs is 5
        runs = 0
        converged = False
        final_structure, final_energy, gibbs_correction = None, None, None
        opt_res, freq_res = None, None
        freq_name = "freq_{}".format(self.full_name.split('_')[-1])

        while not converged and runs < 5:
            runs += 1

            # write input file for the pyscf optimization calculation
            write_xyz(self.file_xyz, self.mol)
            opt_hash = self.local_key()
            freq_hash = self.local_key(freq_name)
            if self.check_if_already_run:
                existing_data = self.existing_data()
                if existing_data:
                    return FWAction(
                        update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier,
                                     "gs_charge": self.gs_charge, "gs_spin": self.gs_spin,
                                     "{}_hash".format(self.full_name): opt_hash,
                                     "{}_hash".format(freq_name): freq_hash,
                                     "{}_geom".format(self.full_name): existing_data.get("geometry"),
                                     "{}_eng".format(self.full_name): existing_data.get("scf_total_energy", {}).get(
                                         "value"),
                                     "{}_gibb".format(self.full_name): self.existing_data(
                                         _hash=freq_hash, calc_name=freq_name).get(
                                         "gibbs_correction", {}).get("value"),
                                     "iop_str": self.iop_str})

            # run pyscf optimization (a finished first run is reused, like the Gaussian log/chk check)
            opt_res = self.load_stage(self.file_result, self.file_log, self.file_chk) if runs == 1 else None
            mf = None
            if opt_res is None:
                print('SUBMITTING PYSCF OPT JOB {} for {}'.format(self.full_name, self.identifier))
                opt_res, mf = run_optimization(self.settings, self.mol, self.charge, self.multiplicity,
                                               output=self.file_log, chkfile=self.file_chk)
                self.save_stage(self.file_result, opt_res)

            # check for normal termination
            if not opt_res.properly_terminated:
                self.mol = Molecule.from_sites([Site.from_dict(sd) for sd in opt_res.sites])  # restart from last geometry
                continue
            if mf is not None:
                write_molden(mf, self.file_molden)  # == formchk
            opt_structure = Molecule.from_sites([Site.from_dict(sd) for sd in opt_res.sites])

            if self.skip_freq:
                final_structure = opt_res.sites
                final_energy = opt_res.energy  # was left as None in Gaussian.py's skip_freq branch (bug fix)
                converged = True
                break

            # run pyscf frequency at the optimized geometry (== Geom=AllCheck Guess=TCheck SCRF=Check)
            write_xyz(self.freq_xyz, opt_structure)
            freq_res = self.load_stage(self.freq_result, self.freq_log, self.freq_chk) if runs == 1 else None
            fmf = None
            if freq_res is None:
                print('SUBMITTING PYSCF FREQUENCY JOB {} for {}'.format(self.full_name, self.identifier))
                try:
                    freq_res, fmf = run_frequencies(self.settings, opt_structure, self.charge, self.multiplicity,
                                                    output=self.freq_log, chkfile=self.freq_chk,
                                                    guess_chk=self.file_chk)
                except Exception as err:
                    raise RuntimeError("Frequency calculation failed Normal termination: {}".format(err))
                self.save_stage(self.freq_result, freq_res)

            if min(freq_res.frequencies) < 0:
                self.mol.perturb(0.1)  # perturb structure to find structure with no negative frequencies
                continue
            else:
                converged = True
                if fmf is not None:
                    write_molden(fmf, self.freq_molden)

            final_structure = opt_res.sites
            final_energy = opt_res.energy
            gibbs_correction = freq_res.gibbs_correction

        if not converged:
            raise RuntimeError("Structure not converged. calc_dir " + self.working_dir)

        # Store results (frequency record first, then the optimization record)
        if not self.skip_freq:
            self.post_job(calc_type=freq_name, name_tag="freq_",
                          data={"gibbs_correction": {"value": freq_res.gibbs_correction, "unit": "Eh"},
                                "zero_point_energy": {"value": freq_res.zpe, "unit": "Eh"},
                                "frequency": freq_res.frequencies, "geometry": freq_res.sites,
                                "scf_total_energy": {"value": freq_res.energy, "unit": "Eh"},
                                "timings": {"pyscf_seconds": freq_res.elapsed}},
                          delete_files=[self.freq_log, self.freq_chk, self.freq_molden, self.freq_xyz])
        self.post_job(data=self.energy_data(opt_res))

        return FWAction(
            update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier, "gs_charge": self.gs_charge,
                         "gs_spin": self.gs_spin,
                         "{}_hash".format(self.full_name): opt_hash,
                         "{}_hash".format(freq_name): freq_hash,
                         "{}_geom".format(self.full_name): final_structure,
                         "{}_eng".format(self.full_name): final_energy,
                         "{}_gibb".format(self.full_name): gibbs_correction,
                         "iop_str": self.get('iop_str')})


class PySCFWtuningJob(WtuningJob):
    """
    WtuningJob (omega tuning with scipy's fminbound, scheme Jh/Jl/J2/Jn2/O2) with Gaussian replaced by PySCF.
    The tuning logic, log format (output.log), directory layout and cycle bookkeeping are inherited unchanged; only the
    "run a calculation and read E/HOMO/LUMO" steps are overridden.
    """

    def __init__(self, *args, solvent=None, solvent_eps=None, density_fit=False, xc=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.solvent, self.solvent_eps, self.density_fit, self.xc = solvent, solvent_eps, density_fit, xc

    def pyscf_settings(self):
        # omega is quantized exactly like the Gaussian IOp string (4 decimals) so the tuning objective is identical
        return PySCFSettings(self.func, self.basis, omega=iop_to_omega(self.omega_iopstr), solvent=self.solvent,
                             solvent_eps=self.solvent_eps, route_parameters=self.route_params,
                             density_fit=self.density_fit, xc=self.xc, max_memory=int(self.mem) * 1000,
                             nthreads=self.nproc)

    def pyscf_extract(self, charge, spin, label):
        """Run one tuning calculation; returns [E (eV), HOMO (eV), LUMO (eV)] like WtuningJob.omega_extract."""
        res, mf = run_energy(self.pyscf_settings(), self.mol, charge, spin,
                             output='{}_tune_{}{}.log'.format(self.name, label, self.cycle))
        return [res.energy * HARTREE_TO_EV, res.homo_ev, res.lumo_ev]

    def geo_opt(self, rps=None):
        res, mf = run_optimization(self.pyscf_settings(), self.mol, self.n_charge, self.n_spin,
                                   output='{}_opt_{}.log'.format(self.name, self.ocycle))
        if not res.properly_terminated:
            warnings.warn('PySCF geometry optimization inside tuning did not converge; using last geometry')
        return Molecule.from_sites([Site.from_dict(sd) for sd in res.sites])

    def omega_gauss_do(self):
        """
        Run the neutral / cation / anion PySCF calculations the scheme needs and extract data from them.
        (Name kept from WtuningJob so wtuning_cycle / omega_wtune / alpha_atune call it unchanged.)
        """
        self.n_e, self.n_homo, self.n_lumo = self.pyscf_extract(self.n_charge, self.n_spin, 'n')
        if self.scheme == 'Jh':
            self.c_e, self.c_homo, self.c_lumo = self.pyscf_extract(self.c_charge, self.c_spin, 'c')
        elif self.scheme == 'Jl':
            self.a_e, self.a_homo, self.a_lumo = self.pyscf_extract(self.a_charge, self.a_spin, 'a')
        else:
            self.c_e, self.c_homo, self.c_lumo = self.pyscf_extract(self.c_charge, self.c_spin, 'c')
            self.a_e, self.a_homo, self.a_lumo = self.pyscf_extract(self.a_charge, self.a_spin, 'a')
        self.log_add({'c_e': getattr(self, 'c_e', None), 'c_homo': getattr(self, 'c_homo', None),
                      'c_lumo': getattr(self, 'c_lumo', None),
                      'n_e': self.n_e, 'n_homo': self.n_homo, 'n_lumo': self.n_lumo})
        self.cycle += 1


@explicit_serialize
class RunPySCFWtuning(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameters
        data_store = get_data_store(self, fw_spec)
        path = env_chk(self.get('path'), fw_spec)
        runfile_log = env_chk(self.get('runfile_log'), fw_spec)
        paramset = self["paramset"]
        identifier = fw_spec.get("identifier", ) or self.get("identifier")
        smiles = fw_spec.get("smiles", ) or self.get("smiles")
        gs_charge = fw_spec.get("gs_charge") or self.get("gs_charge") or get_groundState(
            identifier, smiles, data_store=data_store)
        gs_spin = fw_spec.get("gs_spin") or self.get("gs_spin") or get_groundState(
            identifier, smiles, prop='spin', data_store=data_store)
        pyscf_file_name = fw_spec.get("pyscf_file_name") or self.get("pyscf_file_name") or "pyscf"
        calc_dir = "{}/{}/{}".format(path, identifier or smiles, pyscf_file_name)
        store_dir = get_store_dir(identifier or smiles, get_local_output_dir(self, fw_spec), calc_dir)
        check_if_already_run = spec_or_task(self, fw_spec, "check_if_already_run", True)
        submit = spec_or_task(self, fw_spec, "submit", True)  # API parity; PySCF results are never uploaded
        restricted = spec_or_task(self, fw_spec, "restricted", True)
        solvent = fw_spec.get("solvent", ) or self.get("solvent", )

        radical_electrons = (gs_spin - 1 - paramset.charge) % 2
        paramset.multiplicity = radical_electrons + 1  # calculate spin multiplicity with Hand's rule
        paramset.charge += gs_charge
        geometry = fw_spec.get("geometry", ) or self.get("geometry", )
        geometry_sites = fw_spec.get("{}_geom".format(geometry), )
        if geometry_sites:
            mol = Molecule.from_sites([Site.from_dict(sd) for sd in geometry_sites])
        else:
            geometry_hash = fw_spec.get("{}_hash".format(geometry), )
            mol = get_db_geom(geometry_hash, data_store="local", store_dir=store_dir)
            if mol is None and data_store == "d3tales":
                mol = get_db_geom(geometry_hash, data_store=data_store, store_dir=store_dir)
            mol = mol or start_from_smiles(identifier, smiles, data_store=data_store)

        # create folder
        working_dir = "{}/{}".format(calc_dir, "wtuning")
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)
        write_xyz('wtuning.xyz', mol)

        # check if job has already run (only a PySCF-tuned omega counts: this task exists to tune with PySCF)
        if check_if_already_run:
            tuned_w = read_local_result(store_dir, PYSCF_WTUNING_KEY).get("omega")
            if tuned_w:
                iop_str = omega_to_iop(tuned_w)
                return FWAction(
                    update_spec={"iop_str": iop_str, "gaussrun_dir": calc_dir, "identifier": identifier,
                                 "gs_charge": gs_charge, "gs_spin": gs_spin})

        # run the tuning in pyscf
        functional = paramset.functional if restricted else "R{}".format(paramset.functional)
        wt_mol = PySCFWtuningJob(func=functional, basis=paramset.basis_set, name="wtuning", nproc=nprocs, mem=48,
                                 n_charge=paramset.charge, n_spin=paramset.multiplicity, wdir='./', scheme='Jh',
                                 wbmin=0.05, wbmax=0.5, solvent=solvent, solvent_eps=self.get("solvent_eps"),
                                 density_fit=spec_or_task(self, fw_spec, "density_fit", False), xc=self.get("xc"))
        wt_mol.mol = mol
        wt_mol.wtuning_cycle(max_cycles=0)

        # get the tuned w value from the run
        os.chdir(working_dir)
        with open("output.log", "r") as f:
            output = f.read()
        tuned_w = float(output.strip().split()[-3])
        iop_str = omega_to_iop(tuned_w)

        # Store data (local JSON record; never uploaded)
        conditions = pyscf_conditions(paramset.functional, paramset.basis_set)  # no tuning_parameter, as in back2front
        record = local_record(PYSCF_WTUNING_KEY, identifier, self['name'],
                              {"omega": tuned_w, "conditions": conditions,
                               "output_file": os.path.join(working_dir, "output.log")})
        print("Record saved to", save_local_result(store_dir, record))
        data_hash = pyscf_hash_id(identifier, self['name'], paramset.functional, paramset.basis_set)
        if runfile_log:
            with open(runfile_log, 'a') as fn:
                fn.write("{}, {}\n".format(data_hash, working_dir))
        os.system("rm -rf *.chk *.molden tuning_wtuning_ocycle*")
        return FWAction(
            update_spec={"iop_str": iop_str, "gaussrun_dir": calc_dir, "identifier": identifier, "gs_charge": gs_charge,
                         "gs_spin": gs_spin})


@explicit_serialize
class RunPySCFTDDFT(PySCFBase):

    # Solvent is applied through PySCF's PCM if the installed PySCF supports TDDFT with it.
    def run_task(self, fw_spec):
        setup_obj = self.setup_calc(fw_spec, calc_type='energy')
        if isinstance(setup_obj, FWAction):
            return setup_obj

        # generate the input for the pyscf TDDFT run
        write_xyz(self.file_xyz, self.mol)
        if self.check_if_already_run:
            existing_data = self.existing_data()
            if existing_data:
                return FWAction(
                    update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier, "iop_str": self.iop_str,
                                 "gs_charge": self.gs_charge, "gs_spin": self.gs_spin})

        # run pyscf TDDFT
        print('SUBMITTING PYSCF TDDFT JOB {} for {}'.format(self.full_name, self.identifier))
        res, mf, tds = run_tddft(self.settings, self.mol, self.charge, self.multiplicity,
                                 output=self.file_log, chkfile=self.file_chk)

        # check for normal termination of the job
        if not res.properly_terminated:
            raise RuntimeError(self.file_xyz + " not terminated. calc_dir " + self.working_dir)
        else:
            write_molden(mf, self.file_molden)  # == formchk

            if self.run_nto:
                # NTOs of excited state 1 (first root of the first TDDFT block, as write_nto(log, 1) in Gaussian.py)
                nto_files = write_nto(tds[0][1], mf, self.file_log.rsplit(".", 1)[0] + "_state01", state=1)
                if not all(os.path.isfile(f) for f in nto_files):
                    raise RuntimeError("NTO analysis molden file could not be created")

            data = self.energy_data(res)
            try:
                data.update(parse_tddft(self.file_log))  # descriptors, absorption spectrum, timings (utils.py)
            except Exception as err:  # e.g. no triplet roots requested -> vs0t1 cannot be formed
                warnings.warn("Could not build TDDFT descriptors/spectrum from the log: {}".format(err))
                data["excitations"] = res.excitations
            self.post_job(data=data)
        return FWAction(
            update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier, "iop_str": self.iop_str,
                         "gs_charge": self.gs_charge, "gs_spin": self.gs_spin})
