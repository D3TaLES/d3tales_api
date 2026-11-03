"""ORCA implementation of the D3TaLES electronic-structure FireTasks.

This module mirrors the workflow interface used by Gaussian.py and Psyscf.py,
while executing calculations with the ORCA command-line program. It preserves
D3TaLES/FireWorks naming, storage, geometry, hash, and update_spec conventions.

Runtime requirements:
  * ORCA installed and available as ``orca`` or configured with ``orca_cmd``
  * Open MPI compatible with the installed ORCA release
  * NumPy, FireWorks, atomate, pymatgen, and the D3TaLES API package

ORCA output is parsed directly. GBW files are retained as native checkpoints;
JSON records replace Gaussian formatted-checkpoint files.
"""

import abc
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import uuid
from copy import deepcopy

import numpy as np
from six import add_metaclass

from atomate.utils.utils import env_chk, get_logger
from fireworks import FWAction, FiretaskBase, explicit_serialize

from d3tales_api.Workflows.utils import (
    ElectronicStructureTaskMixin,
    Molecule,
    Site,
    WTUNING_KEY,
    get_data_store,
    get_db_geom,
    get_groundState,
    get_hash_id,
    get_local_output_dir,
    get_omega,
    get_store_dir,
    local_record,
    omega_to_iop,
    orig_conditions,
    orig_hash_id,
    save_local_result,
    spec_or_task,
    start_from_smiles,
    zip_files
)
from d3tales_api.Workflows.d3tales_utils import d3t_submit_calculation

logger = get_logger(__name__)
NPROCS = min(multiprocessing.cpu_count(), 16)
HARTREE_TO_EV = 27.211386245988


def _value(obj, name, default=None):
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def _json_dump(path, data):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x))


def _orca_method(functional):
    """Translate common Gaussian/PySCF functional spellings to ORCA keywords."""
    name = str(functional or "B3LYP").strip()
    aliases = {
        "b3lyp": "B3LYP", "pbe0": "PBE0", "m062x": "M06-2X",
        "m06-2x": "M06-2X", "cam-b3lyp": "CAM-B3LYP",
        "camb3lyp": "CAM-B3LYP", "wb97xd": "wB97X-D3",
        "wb97x-d": "wB97X-D3", "lc-wpbe": "LC-PBE", "hf": "HF",
        "rhf": "HF", "uhf": "HF", "rohf": "HF",
    }
    return aliases.get(name.lower(), name)


def _orca_basis(basis):
    aliases = {
        "def2-svp": "def2-SVP", "def2-tzvp": "def2-TZVP",
        "def2-tzvpp": "def2-TZVPP", "6-31g(d)": "6-31G(d)",
        "6-311+g(d,p)": "6-311+G(d,p)",
    }
    name = str(basis or "def2-SVP").strip()
    return aliases.get(name.lower(), name)


def _omega_from_iop(iop_str):
    if iop_str in (None, ""):
        return None
    if isinstance(iop_str, (int, float)):
        return float(iop_str)
    try:
        value = float(str(iop_str).split(",")[0].strip())
    except ValueError:
        return None
    return value / 10000.0 if abs(value) > 2 else value


def _coords_from_molecule(mol):
    rows = []
    for site in mol.sites:
        symbol = getattr(site.specie, "symbol", str(site.specie))
        x, y, z = (float(v) for v in site.coords)
        rows.append((symbol, x, y, z))
    return rows


def _sites_from_coords(coords):
    try:
        molecule = Molecule([row[0] for row in coords], [row[1:] for row in coords])
        return molecule.as_dict()["sites"]
    except Exception:
        return [{"name": row[0], "xyz": [float(v) for v in row[1:]]} for row in coords]


def _read_xyz(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    natoms = int(lines[0].strip())
    coords = []
    for line in lines[2:2 + natoms]:
        fields = line.split()
        coords.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
    return coords


def _terminated_normally(text):
    return "ORCA TERMINATED NORMALLY" in text


def _parse_energy(text):
    matches = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)", text)
    if not matches:
        raise RuntimeError("Could not find FINAL SINGLE POINT ENERGY in ORCA output")
    return float(matches[-1])


def _parse_frequencies(text):
    values = []
    for match in re.finditer(r"^\s*\d+:\s+(-?\d+(?:\.\d+)?)\s+cm\*\*-1", text, re.MULTILINE):
        values.append(float(match.group(1)))
    return values


def _parse_gibbs_correction(text, electronic_energy=None):
    patterns = [
        r"G-E\(el\)\s+\.\.\.\s+(-?\d+\.\d+)",
        r"Gibbs free energy correction\s+\.\.\.\s+(-?\d+\.\d+)",
    ]
    for pattern in patterns:
        matches = re.findall(pattern, text, re.IGNORECASE)
        if matches:
            return float(matches[-1])
    totals = re.findall(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)", text, re.IGNORECASE)
    if totals and electronic_energy is not None:
        return float(totals[-1]) - electronic_energy
    return None


def _parse_tddft(text):
    states = []
    pattern = re.compile(
        r"STATE\s+(\d+):\s+E=\s*([\d.]+)\s+au\s+([\d.]+)\s+eV.*?f=\s*([\d.Ee+-]+)",
        re.IGNORECASE,
    )
    for state, au, ev, osc in pattern.findall(text):
        states.append({"state": int(state), "energy_ha": float(au), "energy_ev": float(ev), "oscillator_strength": float(osc)})
    if states:
        return states
    # ORCA versions that print oscillator strengths in a separate spectrum table.
    energies = re.findall(r"STATE\s+(\d+):\s+E=\s*([\d.]+)\s+au\s+([\d.]+)\s+eV", text, re.IGNORECASE)
    spectrum = {int(s): float(f) for s, f in re.findall(r"^\s*(\d+)\s+[\d.]+\s+[\d.]+\s+([\d.Ee+-]+)", text, re.MULTILINE)}
    return [{"state": int(s), "energy_ha": float(au), "energy_ev": float(ev), "oscillator_strength": spectrum.get(int(s))} for s, au, ev in energies]


@add_metaclass(abc.ABCMeta)
class OrcaBase(ElectronicStructureTaskMixin, FiretaskBase):
    engine_name = "Orca"
    _fw_name = "OrcaBase"

    def setup_files(self, fw_spec, calc_type="opt"):
        self.data_store = get_data_store(self, fw_spec)
        name_tag = fw_spec.get("name_tag") or self.get("name_tag") or ""
        self.full_name = name_tag + self["name"]
        self.calc_name = self["name"]
        self.paramset = deepcopy(self["paramset"])
        self.identifier = fw_spec.get("identifier") or self.get("identifier", str(uuid.uuid4()))
        self.smiles = fw_spec.get("smiles") or self.get("smiles")
        self.runfile_log = env_chk(self.get("runfile_log"), fw_spec)
        self.skip_freq = fw_spec.get("skip_freq") or self.get("skip_freq") or False
        self.skip_freq = True if self.calc_name == "opt_mol" else self.skip_freq
        self.check_if_already_run = spec_or_task(self, fw_spec, "check_if_already_run", True)
        self.submit = spec_or_task(self, fw_spec, "submit", True)
        self.restricted = spec_or_task(self, fw_spec, "restricted", True)
        self.run_nto = fw_spec.get("run_nto") or self.get("run_nto") or False
        self.orca_cmd = env_chk(self.get("orca_cmd", "orca"), fw_spec)
        self.orca_file_name = fw_spec.get("orca_file_name") or self.get("orca_file_name") or "orca"

        if self.get("prefix"):
            prefix = name_tag + self["prefix"]
        elif self.get("subtype"):
            prefix = name_tag + calc_type + "_" + self["subtype"]
        else:
            prefix = name_tag + calc_type + "_mol"

        if self.get("path"):
            self.calc_dir = "{}/{}/{}".format(env_chk(self.get("path"), fw_spec), self.identifier, self.orca_file_name)
        else:
            self.calc_dir = self.get("calc_dir", os.getcwd())
        self.store_dir = get_store_dir(self.identifier, get_local_output_dir(self, fw_spec), self.calc_dir)

        if self.get("type"):
            self.type = self["type"]
            if self.type == "coupling":
                self.working_dir = "{}/{}/{}/{}/{}".format(self.calc_dir, self.type, calc_type, self["label"], self["subtype"])
            elif self.get("label"):
                self.working_dir = "{}/{}/{}/{}".format(self.calc_dir, self.type, calc_type, self["label"])
            elif self.type.startswith("solv_"):
                self.working_dir = "{}/{}/{}/{}".format(self.calc_dir, self.type, calc_type, self["subtype"])
            else:
                self.working_dir = "{}/{}/{}".format(self.calc_dir, self.type, calc_type)
        else:
            self.type = "mol_{}".format(calc_type)
            self.working_dir = "{}/{}".format(self.calc_dir, self.type)

        os.makedirs(self.working_dir, exist_ok=True)
        os.chdir(self.working_dir)
        self.prefix = prefix
        self.file_input = os.path.join(self.working_dir, prefix + ".inp")
        self.file_log = os.path.join(self.working_dir, prefix + ".out")
        self.file_gbw = os.path.join(self.working_dir, prefix + ".gbw")
        self.file_xyz = os.path.join(self.working_dir, prefix + ".xyz")
        self.file_result = os.path.join(self.working_dir, prefix + ".result.json")
        self.freq_result = os.path.join(self.working_dir, "freq_" + prefix + ".result.json")
        # Compatibility names expected by legacy workflow code.
        self.file_com = self.file_input
        self.file_chk = self.file_gbw
        self.file_fchk = self.file_result

    def setup_calc(self, fw_spec, calc_type="opt"):
        self.setup_files(fw_spec, calc_type)
        self.solvent = fw_spec.get("solvent") or self.get("solvent")
        self.gs_charge = fw_spec.get("gs_charge") or self.get("gs_charge") or get_groundState(self.identifier, self.smiles, data_store=self.data_store) or 0
        self.gs_spin = fw_spec.get("gs_spin") or self.get("gs_spin") or get_groundState(self.identifier, self.smiles, prop="spin", data_store=self.data_store) or 1
        use_iop = fw_spec.get("use_iop", True) if self.get("use_iop", True) else self.get("use_iop")
        self.iop_str = fw_spec.get("iop_str") or self.get("iop_str")
        if not self.iop_str and use_iop:
            tuned_w = get_omega(self.identifier, data_store=self.data_store, store_dir=self.store_dir)
            if tuned_w:
                self.iop_str = omega_to_iop(tuned_w)
        self.omega = _omega_from_iop(self.iop_str) if use_iop else None

        geometry = fw_spec.get("geometry") or self.get("geometry")
        geometry_sites = fw_spec.get("{}_geom".format(geometry))
        if geometry_sites:
            self.mol = Molecule.from_sites([Site.from_dict(sd) for sd in geometry_sites])
        else:
            geometry_hash = fw_spec.get("{}_hash".format(geometry), orig_hash_id(
                self.identifier, self.calc_name, _value(self.paramset, "functional"),
                _value(self.paramset, "basis_set"), tuning_parameter=self.iop_str, solvent=self.solvent
            ))
            self.mol = get_db_geom(geometry_hash, data_store=self.data_store, store_dir=self.store_dir) or start_from_smiles(
                self.identifier, self.smiles, data_store=self.data_store
            )
        if len(self.mol.sites) > 200:
            return FWAction(defuse_children=True, update_spec={"defuser_reason": "The molecule has {} atoms".format(len(self.mol.sites))})

        self.charge = int(_value(self.paramset, "charge", 0)) + int(self.gs_charge)
        radical_electrons = (int(self.gs_spin) - 1 - int(_value(self.paramset, "charge", 0))) % 2
        self.multiplicity = int(radical_electrons + 1)
        self.functional = _orca_method(_value(self.paramset, "functional", "B3LYP"))
        self.basis = _orca_basis(_value(self.paramset, "basis_set", "def2-SVP"))
        print("Charge:", self.charge, "\t Multiplicity:", self.multiplicity)




    def input_text(self, job, coords=None, extra_keywords=None, blocks=None, charge=None, multiplicity=None, prefix=None):
        keywords = [self.functional, self.basis, "TightSCF"]
        keywords.extend(extra_keywords or [])
        if job:
            keywords.append(job)
        if self.solvent:
            keywords.append("CPCM({})".format(self.solvent))
        if not self.restricted and (charge if charge is not None else self.charge) == 0:
            keywords.append("UHF")
        lines = ["! " + " ".join(keywords), "%pal nprocs {} end".format(int(self.get("nprocs", NPROCS)))]
        lines.append("%maxcore {}".format(int(self.get("maxcore_mb", 3000))))
        if self.omega is not None:
            lines.extend(["%method", "  RangeSepMu {:.10f}".format(self.omega), "end"])
        for block in blocks or []:
            lines.append(block.rstrip())
        rows = coords or _coords_from_molecule(self.mol)
        lines.append("* xyz {} {}".format(self.charge if charge is None else charge, self.multiplicity if multiplicity is None else multiplicity))
        for symbol, x, y, z in rows:
            lines.append("  {:<3s} {: .12f} {: .12f} {: .12f}".format(symbol, x, y, z))
        lines.append("*")
        return "\n".join(lines) + "\n"

    def run_orca(self, input_text, input_path=None, output_path=None):
        input_path = input_path or self.file_input
        output_path = output_path or self.file_log
        with open(input_path, "w", encoding="utf-8") as handle:
            handle.write(input_text)
        executable = shutil.which(self.orca_cmd) if os.path.sep not in str(self.orca_cmd) else self.orca_cmd
        if not executable or not os.path.exists(executable):
            raise RuntimeError("ORCA executable was not found. Set orca_cmd to the ORCA binary path.")
        with open(output_path, "w", encoding="utf-8") as output:
            proc = subprocess.run([executable, input_path], cwd=self.working_dir, stdout=output, stderr=subprocess.STDOUT, check=False)
        with open(output_path, encoding="utf-8", errors="replace") as handle:
            text = handle.read()
        if proc.returncode != 0 or not _terminated_normally(text):
            raise RuntimeError("ORCA failed for {}. See {}".format(input_path, output_path))
        return text

    def post_job(self, record, upload_files=None, calc_type=None, name_tag=""):
        calc_type = calc_type or self.calc_name
        result_path = self.file_result if calc_type == self.calc_name else self.freq_result
        _json_dump(result_path, record)
        files = [f for f in (upload_files or [self.file_input, self.file_log, self.file_gbw, result_path]) if os.path.isfile(f)]
        data_hash = get_hash_id(self.identifier, self.file_input, calc_type, output_file=self.file_log)
        if self.data_store == "d3tales":
            archive = zip_files(files, zip_name="{}_{}.zip".format(self.identifier, name_tag + self.full_name))
            if self.submit:
                d3t_submit_calculation(archive, self.identifier, calc_type)
        else:
            print("Record saved to", save_local_result(self.store_dir, record, key=self.local_key(calc_type)))
        with open(self.runfile_log, "a", encoding="utf-8") as handle:
            handle.write("{}, {}\n".format(data_hash, self.working_dir))


@explicit_serialize
class RunOrcaEnergy(OrcaBase):
    def run_task(self, fw_spec):
        action = self.setup_calc(fw_spec, "energy")
        if isinstance(action, FWAction):
            return action
        if self.check_if_already_run:
            existing = self.existing_data()
            if existing:
                return FWAction(update_spec={"orcarun_dir": self.calc_dir, "identifier": self.identifier, "gs_charge": self.gs_charge, "gs_spin": self.gs_spin, "iop_str": self.iop_str, "{}_eng".format(self.full_name): existing.get("scf_total_energy", {}).get("value")})
        text = self.run_orca(self.input_text("SP"))
        energy = _parse_energy(text)
        record = local_record(self.local_key(), self.identifier, self.calc_name, {
            "engine": "ORCA", "geometry": self.mol.as_dict()["sites"],
            "scf_total_energy": {"value": energy, "unit": "Ha"},
            "conditions": self.local_conditions(), "output_file": self.file_log,
            "checkpoint_file": self.file_gbw,
        })
        self.post_job(record)
        return FWAction(update_spec={"orcarun_dir": self.working_dir, "gaussrun_dir": self.working_dir, "identifier": self.identifier, "iop_str": self.iop_str, "gs_charge": self.gs_charge, "gs_spin": self.gs_spin, "{}_eng".format(self.full_name): energy})


@explicit_serialize
class RunOrcaOpt(OrcaBase):
    def run_task(self, fw_spec):
        action = self.setup_calc(fw_spec, "opt")
        if isinstance(action, FWAction):
            return action
        freq_name = "freq_{}".format(self.full_name.split("_")[-1])
        if self.check_if_already_run:
            existing = self.existing_data()
            if existing:
                freq_existing = self.existing_data(freq_name) or {}
                return FWAction(update_spec={"orcarun_dir": self.calc_dir, "gaussrun_dir": self.calc_dir, "identifier": self.identifier, "gs_charge": self.gs_charge, "gs_spin": self.gs_spin, "{}_hash".format(self.full_name): self.local_key(), "{}_hash".format(freq_name): self.local_key(freq_name), "{}_geom".format(self.full_name): existing.get("geometry"), "{}_eng".format(self.full_name): existing.get("scf_total_energy", {}).get("value"), "{}_gibb".format(self.full_name): freq_existing.get("gibbs_correction", {}).get("value"), "iop_str": self.iop_str})

        job = "TightOpt" if self.skip_freq else "TightOpt Freq"
        blocks = ["%geom\n  MaxIter {}\nend".format(int(self.get("max_opt_steps", 100)))]
        text = self.run_orca(self.input_text(job, blocks=blocks))
        energy = _parse_energy(text)
        xyz_candidates = [self.file_xyz, os.path.join(self.working_dir, self.prefix + ".xyz")]
        coords = _coords_from_molecule(self.mol)
        for path in xyz_candidates:
            if os.path.isfile(path):
                coords = _read_xyz(path)
                break
        sites = _sites_from_coords(coords)
        gibbs = None
        if not self.skip_freq:
            frequencies = _parse_frequencies(text)
            imaginary = [f for f in frequencies if f < -1.0]
            if imaginary:
                raise RuntimeError("ORCA optimization retained imaginary frequencies: {}".format(imaginary))
            gibbs = _parse_gibbs_correction(text, energy)
            freq_record = local_record(self.local_key(freq_name), self.identifier, freq_name, {
                "engine": "ORCA", "geometry": sites,
                "frequencies": {"value": frequencies, "unit": "cm^-1"},
                "gibbs_correction": {"value": gibbs, "unit": "Ha"},
                "conditions": self.local_conditions(), "output_file": self.file_log,
            })
            self.post_job(freq_record, calc_type=freq_name, name_tag="freq_")
        opt_record = local_record(self.local_key(), self.identifier, self.calc_name, {
            "engine": "ORCA", "geometry": sites,
            "scf_total_energy": {"value": energy, "unit": "Ha"},
            "conditions": self.local_conditions(), "output_file": self.file_log,
            "checkpoint_file": self.file_gbw,
        })
        self.post_job(opt_record)
        return FWAction(update_spec={"orcarun_dir": self.calc_dir, "gaussrun_dir": self.calc_dir, "identifier": self.identifier, "gs_charge": self.gs_charge, "gs_spin": self.gs_spin, "{}_hash".format(self.full_name): self.local_key(), "{}_hash".format(freq_name): self.local_key(freq_name), "{}_geom".format(self.full_name): sites, "{}_eng".format(self.full_name): energy, "{}_gibb".format(self.full_name): gibbs, "iop_str": self.iop_str})


@explicit_serialize
class RunOrcaTDDFT(OrcaBase):
    def run_task(self, fw_spec):
        action = self.setup_calc(fw_spec, "energy")
        if isinstance(action, FWAction):
            return action
        if self.check_if_already_run and self.existing_data():
            return FWAction(update_spec={"orcarun_dir": self.calc_dir, "gaussrun_dir": self.calc_dir, "identifier": self.identifier, "iop_str": self.iop_str, "gs_charge": self.gs_charge, "gs_spin": self.gs_spin})
        nroots = int(self.get("nstates", self.get("n_states", 10)))
        tddft_block = "%tddft\n  nroots {}\n  tda {}\n  DoNTO {}\nend".format(nroots, "true" if self.get("tda", False) else "false", "true" if self.run_nto else "false")
        text = self.run_orca(self.input_text("SP", blocks=[tddft_block]))
        energy = _parse_energy(text)
        states = _parse_tddft(text)
        nto_files = []
        if self.run_nto:
            nto_files = [os.path.join(self.working_dir, f) for f in os.listdir(self.working_dir) if "nto" in f.lower()]
        record = local_record(self.local_key(), self.identifier, self.calc_name, {
            "engine": "ORCA", "geometry": self.mol.as_dict()["sites"],
            "scf_total_energy": {"value": energy, "unit": "Ha"},
            "excited_states": states,
            "excitation_energies": {"value": [s["energy_ha"] for s in states], "unit": "Ha"},
            "excitation_energies_ev": {"value": [s["energy_ev"] for s in states], "unit": "eV"},
            "oscillator_strengths": [s["oscillator_strength"] for s in states],
            "nto_files": nto_files, "conditions": self.local_conditions(), "output_file": self.file_log,
        })
        self.post_job(record, upload_files=[self.file_input, self.file_log, self.file_gbw, self.file_result] + nto_files)
        return FWAction(update_spec={"orcarun_dir": self.calc_dir, "gaussrun_dir": self.calc_dir, "identifier": self.identifier, "iop_str": self.iop_str, "gs_charge": self.gs_charge, "gs_spin": self.gs_spin})


@explicit_serialize
class RunWtuning(FiretaskBase):
    """Tune ORCA's range-separation parameter using an IP-theorem objective."""

    def run_task(self, fw_spec):
        # Reuse OrcaBase machinery through a temporary energy task, then perform
        # bounded scalar minimization over neutral/cation/anion single points.
        task = RunOrcaEnergy(**dict(self))
        action = task.setup_calc(fw_spec, "wtuning")
        if isinstance(action, FWAction):
            return action
        if task.check_if_already_run:
            tuned = get_omega(task.identifier, data_store=task.data_store, store_dir=task.store_dir)
            if tuned:
                return FWAction(update_spec={"iop_str": omega_to_iop(tuned), "orcarun_dir": task.calc_dir, "gaussrun_dir": task.calc_dir, "identifier": task.identifier, "gs_charge": task.gs_charge, "gs_spin": task.gs_spin})
        from scipy.optimize import minimize_scalar

        def run_state(charge, multiplicity, omega, label):
            task.omega = float(omega)
            inp = os.path.join(task.working_dir, "{}_{}.inp".format(label, str(omega).replace(".", "_")))
            out = os.path.splitext(inp)[0] + ".out"
            text = task.run_orca(task.input_text("SP", charge=charge, multiplicity=multiplicity), inp, out)
            energy = _parse_energy(text)
            homo_matches = re.findall(r"^\s*\d+\s+[12]\.0000\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", text, re.MULTILINE)
            if not homo_matches:
                raise RuntimeError("Could not parse occupied orbital energies during ORCA omega tuning")
            homo_ev = max(float(row[1]) for row in homo_matches)
            return energy, homo_ev / HARTREE_TO_EV

        def objective(omega):
            spin = task.multiplicity - 1
            e_n, h_n = run_state(task.charge, task.multiplicity, omega, "neutral")
            e_c, _ = run_state(task.charge + 1, 2 if spin == 0 else max(1, task.multiplicity - 1), omega, "cation")
            e_a, h_a = run_state(task.charge - 1, 2 if spin == 0 else task.multiplicity + 1, omega, "anion")
            return (h_n + e_c - e_n) ** 2 + (h_a + e_n - e_a) ** 2

        result = minimize_scalar(objective, bounds=(float(self.get("wbmin", 0.05)), float(self.get("wbmax", 0.5))), method="bounded", options={"xatol": float(self.get("omega_tolerance", 1e-3))})
        if not result.success:
            raise RuntimeError("ORCA omega tuning failed: {}".format(result.message))
        tuned = float(result.x)
        record = local_record(WTUNING_KEY, task.identifier, self["name"], {"engine": "ORCA", "omega": tuned, "objective": float(result.fun), "conditions": orig_conditions(task.functional, task.basis)})
        _json_dump(os.path.join(task.working_dir, "wtuning.result.json"), record)
        if task.data_store == "local":
            save_local_result(task.store_dir, record)
        elif task.submit:
            archive = zip_files([f for f in os.listdir(task.working_dir) if f.endswith((".inp", ".out", ".json"))], zip_name="{}_wtuning.zip".format(task.identifier))
            d3t_submit_calculation(archive, task.identifier, self["name"])
        return FWAction(update_spec={"iop_str": omega_to_iop(tuned), "orcarun_dir": task.calc_dir, "gaussrun_dir": task.calc_dir, "identifier": task.identifier, "gs_charge": task.gs_charge, "gs_spin": task.gs_spin})


# Compatibility aliases let existing workflow definitions switch modules with
# minimal changes while exposing ORCA-native class names for new workflows.
GaussianBase = OrcaBase
RunGaussianEnergy = RunOrcaEnergy
RunGaussianOpt = RunOrcaOpt
RunGaussianTDDFT = RunOrcaTDDFT
