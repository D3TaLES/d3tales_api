import os
import re
import json
import hashlib
import glob
import zipfile
from monty.json import MontyEncoder

from pymatgen.core.sites import Site
from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianInput, GaussianOutput

from rdkit.Chem import MolFromSmiles, AllChem, MolToInchiKey
from rdkit.Chem.rdmolops import GetFormalCharge, AddHs
from rdkit.Chem.Descriptors import NumRadicalElectrons
from d3tales_api.Workflows.d3tales_utils import (d3t_get_smiles, d3t_get_init_structure, d3t_get_ground_state,
                                                 d3t_get_omega, d3t_get_raw_computation)
file_path = os.path.dirname(os.path.realpath(__file__))

# Copyright 2021, University of Kentucky

# ENVIRONMENT VARIABLES

ZINDOLIB = ">>ZINDOLIB<<"
ZINDOBIN = ">>ZINDOBIN<<"
ZINDOCTBIN = ">>ZINDOCTBIN<<"
DB_FILE = ">>db_file<<"
AFLOW = ">>aflow<<"
VASP_CMD = ">>vasp_cmd<<"
GAMMA_VASP_CMD = ">>gamma_vasp_cmd<<"
G16_CMD = ">>g16_cmd<<"
MA_DEF2 = ">>ma_def2<<"
PATH = ">>path<<"
RUNFILE_LOG = ">>runfile_log<<"
PROC_VM_KEY = ">>processing_vm_key<<"
LOCAL_OUTPUT_DIR = ">>local_output_dir<<"  # optional; set under `env:` in my_fworker.yaml

# DATA STORES
# "d3tales": read/write through the D3TaLES website + databases (original behavior, default)
# "local":   no database; each calculation's record is written as a JSON file
DATA_STORES = ("d3tales", "local")
DEFAULT_DATA_STORE = "d3tales"
LOCAL_RESULTS_DIR = "results"  # used under the calculation directory when no local_output_dir is set
WTUNING_KEY = "wtuning"


# DATA STORE FUNCTIONS
def get_data_store(task=None, fw_spec=None):
    """
    Get the data store for a FireTask. Precedence: fw_spec, then task parameter, then DEFAULT_DATA_STORE.
    """
    fw_spec = fw_spec or {}
    store = fw_spec.get("data_store") or (task.get("data_store") if task else None) or DEFAULT_DATA_STORE
    store = str(store).lower()
    if store not in DATA_STORES:
        raise ValueError("data_store must be one of {}, not {}".format(DATA_STORES, store))
    return store


def get_local_output_dir(task=None, fw_spec=None):
    """
    Folder for local JSON records. Precedence:
        1. `local_output_dir` in the fw_spec or task kwargs (a path, or an env key such as ">>local_output_dir<<")
        2. `local_output_dir` in the FireWorker env (my_fworker.yaml), same mechanism as g16_cmd
        3. None -> records go to <calc_dir>/results
    Unlike env_chk(strict=True), a missing FireWorker key returns None instead of raising.
    """
    fw_spec = fw_spec or {}
    val = fw_spec.get("local_output_dir") or (task.get("local_output_dir") if task else None) or LOCAL_OUTPUT_DIR
    if isinstance(val, str) and val.startswith(">>") and val.endswith("<<"):
        val = fw_spec.get("_fw_env", {}).get(val[2:-2])
    return val or None


def get_store_dir(identifier, local_output_dir=None, calc_dir=None):
    """
    Directory holding the local JSON records for one molecule.
        local_output_dir set: <local_output_dir>/<identifier>
        otherwise:            <calc_dir>/results
    Always absolute, because tasks chdir into their working directories.
    """
    if local_output_dir:
        return os.path.join(os.path.abspath(os.path.expandvars(os.path.expanduser(str(local_output_dir)))),
                            str(identifier))
    if calc_dir:
        return os.path.join(os.path.abspath(calc_dir), LOCAL_RESULTS_DIR)
    return None


def local_identifier(smiles):
    """Database-free molecule identifier (InChIKey; filesystem safe)."""
    return MolToInchiKey(MolFromSmiles(smiles))


def omega_to_iop(tuned_w):
    return str(int(tuned_w * 1e4)).zfill(5) + "00000"


def local_record(_hash, identifier, calc_type, data):
    """Envelope that mirrors a D3TaLES backend computation document."""
    return {"_id": _hash, "mol_id": identifier, "calculation_type": calc_type, "data": data}


def write_json(filename, data):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    tmp_filename = filename + ".tmp"
    with open(tmp_filename, 'w') as fn:
        json.dump(data, fn, cls=MontyEncoder, indent=2)
    os.replace(tmp_filename, filename)  # atomic: safe when several Fireworks write to the same folder


def save_local_result(store_dir, record):
    """Write a record to <store_dir>/<calculation_type>__<_id>.json (one file per calculation)."""
    filename = os.path.join(store_dir, "{}__{}.json".format(record["calculation_type"], record["_id"]))
    write_json(filename, record)
    return filename


def read_local_result(store_dir, key):
    """Return the `data` dict of the record whose _id is `key` (newest if several); {} if none."""
    if not store_dir or not key:
        return {}
    matches = sorted(glob.glob(os.path.join(store_dir, "*__{}.json".format(key))), key=os.path.getmtime)
    if not matches:
        return {}
    with open(matches[-1], 'r') as fn:
        return json.load(fn).get("data", {})


def parse_gaussian_log(logfile, conditions=None):
    """
    Minimal, database-free summary of a Gaussian log file. Keys mirror the D3TaLES backend `data` fields
    used by the workflow and by back2front (geometry, scf_total_energy, gibbs_correction, charge, ...).
    Energies are in Hartree.
    """
    gout = GaussianOutput(logfile)
    data = {"conditions": conditions or {}, "log_file": os.path.abspath(logfile),
            "properly_terminated": gout.properly_terminated,
            "charge": gout.charge, "spin_multiplicity": gout.spin_multiplicity}
    try:
        data["geometry"] = gout.final_structure.as_dict()['sites']
    except Exception:
        pass
    if gout.final_energy is not None:
        data["scf_total_energy"] = {"value": gout.final_energy, "unit": "Hartree"}
    gibbs = (gout.corrections or {}).get("Gibbs Free Energy")
    if gibbs is not None:
        data["gibbs_correction"] = {"value": gibbs, "unit": "Hartree"}
    return data


def get_omega(identifier, data_store=DEFAULT_DATA_STORE, store_dir=None):
    """Return a previously tuned omega value (float) or None."""
    if data_store == "d3tales":
        return d3t_get_omega(identifier)
    return read_local_result(store_dir, WTUNING_KEY).get("omega")


def get_existing_data(_hash, data_store=DEFAULT_DATA_STORE, store_dir=None):
    """Return the stored data dict for a calculation, or {} if it has not been run."""
    if data_store == "d3tales":
        return d3t_get_raw_computation(_hash) or {}
    return read_local_result(store_dir, _hash)


# FUNCTIONS

def runtime_from_log(logfile):
    import re

    time_patt = re.compile(r"\d+\.d+|\d+")
    time_data = []

    with open(logfile, "r") as f:
        line = f.readline()

        while line != "":
            if re.match("Job cpu time", line.strip()):
                time_data.extend(time_patt.findall(line))
            line = f.readline()

    if time_data:
        time_data = [float(time) for time in time_data]
        runtime = (time_data[0] * 86400 + time_data[1] * 3600 + time_data[2] * 60 + time_data[3]) / 3600
    else:
        runtime = 0

    return round(runtime, 3)


def generate_gaussian_input(paramset=None, mol=None, dieze_tag="#P"):
    route_parameters = paramset.route_parameters
    input_parameters = paramset.input_parameters
    link0_parameters = paramset.link0_parameters
    charge = paramset.charge
    multiplicity = paramset.multiplicity
    functional = paramset.functional
    basis_set = paramset.basis_set

    ginput = GaussianInput(mol=mol, charge=charge, spin_multiplicity=multiplicity,
                           title=None, functional=functional, basis_set=basis_set, link0_parameters=link0_parameters,
                           route_parameters=route_parameters, input_parameters=input_parameters, dieze_tag=dieze_tag)
    return ginput


def get_central_dihed(rdk_mol):
    """
    function for recording central dihedral atom numbers
    """
    # find all potential center bonds (single bonds between carbons)
    potential_cntbond = []
    for bond in rdk_mol.GetBonds():
        if str(bond.GetBondType()) == 'SINGLE':
            atom_a = bond.GetBeginAtom()
            atom_b = bond.GetEndAtom()
            bond_atoms = [atom_a, atom_b]
            if atom_a.GetAtomicNum() == 6 and atom_b.GetAtomicNum() == 6:
                potential_cntbond.append(bond_atoms)

    # find central bond
    num_cnt_bond = int((len(potential_cntbond) - 1) / 2)
    cnt_bond = potential_cntbond[num_cnt_bond]

    cntatom_a = cnt_bond[0]
    cntatom_b = cnt_bond[1]
    dihed1 = []
    dihed2 = []

    # assemble list of atoms in first dihedral angle
    for n in cntatom_a.GetNeighbors():
        if n.GetIdx() != cntatom_b.GetIdx():
            dihed1.append(n.GetIdx())
            dihed1.extend((cntatom_a.GetIdx(), cntatom_b.GetIdx()))
            break
    for n in cntatom_b.GetNeighbors():
        if n.GetIdx() != cntatom_a.GetIdx():
            dihed1.append(n.GetIdx())
            break
    # assemble list of atoms in second dihedral angle
    for n in cntatom_a.GetNeighbors():
        dihed2.append(n.GetIdx()) if n.GetIdx() not in dihed1 else dihed1
    dihed2.extend((cntatom_a.GetIdx(), cntatom_b.GetIdx()))
    for n in cntatom_b.GetNeighbors():
        dihed2.append(n.GetIdx()) if n.GetIdx() not in dihed1 else dihed1

    return [dihed1, dihed2]


def get_scfs(logfile):
    energies = []
    with open(logfile, "r") as f:
        line = f.readline()
        while line != "":
            if re.search(r"^\sSCF Done:", line):
                energies.append(line.split(" ")[4])
    return energies


def parse_tddft(tddft_logfile):
    # read the log files
    excitations = get_tddft_excitations(tddft_logfile)

    # get descriptors
    singlet_excitation_energy = excitations["Singlet"]
    triplet_excitation_energy = excitations["Triplet"]

    # get absorption spectrum (pymatgen.io.gaussian Gaussian.get_spectre_plot truncated)
    import scipy.constants as cst
    import numpy as np
    from scipy.stats import norm

    transitions = singlet_excitation_energy
    sigma = 0.50
    step = 0.01

    minval = min([val[0] for val in transitions]) - 5.0 * sigma
    maxval = max([val[0] for val in transitions]) + 5.0 * sigma
    npts = int((maxval - minval) / step) + 1

    eneval = np.linspace(minval, maxval, npts)  # in eV
    lambdaval = [cst.h * cst.c / (val * cst.e) * 1.e9
                 for val in eneval]  # in nm

    # sum of gaussian functions
    spectre = np.zeros(npts)
    for trans in transitions:
        spectre += trans[2] * norm.pdf(eneval, trans[0], sigma)
    spectre /= spectre.max()
    singlet_spectra_data = {"energies": eneval, "lambda": lambdaval, "xas": spectre}

    # for plotting in frontend
    abs_plot = [{"x": [round(x, 0) for x in singlet_spectra_data["lambda"]],
                 "y": singlet_spectra_data["xas"], "mode": 'lines',
                 "name": "absorption",
                 "line": {
                     "color": "#003396",
                     "width": 3
                 }}]
    singlet_spectra_data.update({"abs_plot": abs_plot})

    # add results to DB
    descriptors = {
        "vs0s1": singlet_excitation_energy[0][0],
        "vs0t1": triplet_excitation_energy[0][0],
    }

    visuals = {"exec_singlet": singlet_spectra_data}

    timings = {"tddft": runtime_from_log(tddft_logfile)}

    d_results = {"descriptors": descriptors,
                 "visuals": visuals,
                 "timings": timings,
                 }

    return d_results


def get_hash_id(identifier, filepath, calc_type, output_file=None, **kwargs):
    from d3tales_api.Processors.parser_dft import ProcessGausLog  # lazy: Processors is the next module to decouple
    metadata = {"calculation_type": calc_type, "mol_file": filepath, "wtuning_output": output_file, **kwargs}
    process_data = ProcessGausLog(_id=identifier, metadata=metadata)
    return process_data.hash_id


def orig_hash_id(_id, calculation_type, functional, basis_set, tuning_parameter=None, solvent=None):
    """
    Hash ID
    """
    conditions = orig_conditions(functional=functional, basis_set=basis_set, tuning_parameter=tuning_parameter, solvent=solvent)
    hash_dict = {
        "_id": _id,
        "calculation_type": calculation_type,
        "conditions": conditions,
    }
    print("CONDITIONS", conditions, hash_dict)
    dhash = hashlib.md5()
    encoded = json.dumps(hash_dict, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()


def orig_conditions(functional, basis_set, tuning_parameter=None, solvent=None, ):
    """
    Dictionary of conditions (in accordance with D3TaLES backend schema)
    """
    data_dict = {
        "data_source": 'dft',
        "code_name": 'Gaussian',
        "code_version": '16',
        "functional": functional,
        "basis_set": basis_set,
    }
    if tuning_parameter:
        if "." not in str(tuning_parameter):
            tuning_parameter = float('0.{}'.format(tuning_parameter[1:5]))
        data_dict['tuning_parameter'] = tuning_parameter
    if solvent:
        dielectric_constant = {"acetonitrile": 35.688}
        data_dict['solvent'] = {
            'name': solvent.capitalize(),
            'model': 'implicit_solvent',
            'dielectric_constant': dielectric_constant.get(solvent.lower())
        }
    return data_dict


def start_from_smiles(identifier, smiles, data_store=DEFAULT_DATA_STORE):
    """
    get input geometry (from the frontend DB if data_store is d3tales, else from SMILES)
    :return: pymatgen mol object
    """
    if data_store == "d3tales":
        try:
            return d3t_get_init_structure(identifier)
        except Exception:
            smiles = smiles or d3t_get_smiles(identifier)
    if not smiles:
        raise ValueError("Cannot generate a starting geometry for {}: no SMILES provided.".format(identifier))
    from d3tales_api.Processors.info_from_smiles import find_lowest_e_conf  # lazy: see get_hash_id
    structure = find_lowest_e_conf(smiles)
    return Molecule.from_str(structure, 'xyz')


def get_groundState(identifier, smiles=None, prop='charge', data_store=DEFAULT_DATA_STORE):
    """
    get ground state charge or spin (from the frontend DB if data_store is d3tales, else from SMILES)
    :return: int
    """
    if data_store == "d3tales":
        value = d3t_get_ground_state(identifier, prop)
        if value is not None:
            return value
    if smiles:
        rdkmol = MolFromSmiles(smiles)
        rdkmol_hs = AddHs(rdkmol)
        AllChem.EmbedMolecule(rdkmol_hs)
        return GetFormalCharge(rdkmol) if prop == 'charge' else NumRadicalElectrons(rdkmol) + 1  # calculate spin multiplicity with Hand's rule
    raise ValueError(f"No ground state found for {identifier} (data_store={data_store}). A SMILES must be specified")


def get_db_geom(input_geometry_hash, error_raise=False, data_store=DEFAULT_DATA_STORE, store_dir=None):
    """
    get input geometry from the data store
    :param store_dir: str, folder holding local JSON records (local data store only)
    :return: pymatgen mol object
    """
    if data_store == "d3tales":
        record = d3t_get_raw_computation(input_geometry_hash, fields="data.geometry")
    else:
        record = read_local_result(store_dir, input_geometry_hash) or None
    if record is None:
        if error_raise:
            raise LookupError("Starting geometry with hash {} not found.".format(input_geometry_hash))
        return None
    dict_sites = record.get('geometry')
    if not dict_sites:
        if error_raise:
            raise LookupError("There is not starting gemetry for computation item {}".format(input_geometry_hash))
        return None
    return Molecule.from_sites([Site.from_dict(sd) for sd in dict_sites])


def zip_files(file_list, zip_name='upload_zip.zip'):
    zip_path = os.path.join(os.getcwd(), zip_name)
    with zipfile.ZipFile(zip_path, 'w') as zipF:
        for file in file_list:
            zipF.write(file, compress_type=zipfile.ZIP_DEFLATED)
    return zip_path


def get_tddft_excitations(logfile):
    """
    Read a excitation energies after a TD-DFT calculation.

    Returns:

        A list: A list of tuple for each transition such as
                [(energie (eV), lambda (nm), oscillatory strength), ... ]
    """

    float_patt = re.compile(r"\s*([+-]?\d+\.\d+)")
    state_patt = re.compile(r"[a-zA-Z]*let")
    transitions = {"Singlet": [],
                   "Triplet": []
                   }

    # read in file
    with open(logfile, "r") as f:
        line = f.readline()
        td = False
        while line != "":
            if re.search(r"^\sExcitation energies and oscillator strengths:", line):
                td = True

            if td:
                if re.search(r"^\sExcited State\s*\d", line):
                    val = [float(v) for v in float_patt.findall(line)]
                    try:
                        state = state_patt.findall(line)[0]
                    except:
                        state_val = val.pop(0)
                        if state_val == 1:
                            state = "Singlet"
                        else:
                            state = "Triplet"
                    # transitions.append(tuple(val[0:3]))
                    transitions[state].append(tuple(val))
            line = f.readline()
    return transitions


def write_nto(tddft_log, excitation_num, tddft_chk=None):
    """
    Run NTO analysis after a TDDFT calculation. https://gaussian.com/faq4/
    :param tddft_log: str, path to TDDFT calculation log file
    :param excitation_num: int, number for excitation on which to run the NTO analysis
    :param tddft_chk: str, path to TDDFT calculation checkpoint file (assumed to be the same name as log if not specified)

    """
    tddft_chk = tddft_chk or tddft_log.split(".")[0] + ".chk"
    nto_com = "{}_state{:02d}.com".format(tddft_log.split(".")[0], excitation_num)
    nto_chk = "{}.chk".format(nto_com.split(".")[0])

    gout = GaussianOutput(tddft_log)
    ginp = gout.to_input()
    ginp.route_parameters = ginp.route_parameters or {}
    if ginp.route_parameters.get("iop(3107", ):
        iop_str = ginp.route_parameters["iop(3107"].split(",")[0]
        ginp.route_parameters = {"iop(3/107={}, 3/108={})".format(iop_str, iop_str): ""}
    ginp.route_parameters.update({"Geom": "AllCheck", "Guess": "(Read,Only)", "Density": "(Check,Transition={})".format(excitation_num),
                                  "Pop": "(Minimal,NTO,SaveNTO)"})
    ginp.link0_parameters.update({"%chk": nto_chk, "%oldchk": tddft_chk})
    ginp.write_file(nto_com, cart_coords=True)
    return nto_com, nto_chk


if __name__ == "__main__":
    get_groundState("05AATS")