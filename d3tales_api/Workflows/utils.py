import os
import zipfile
from d3tales_api.Processors.d3tales_parser import *
from d3tales_api.D3database.restapi import RESTAPI
from d3tales_api.D3database.info_from_smiles import find_lowest_e_conf

from pymatgen.core.sites import Site
from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianInput, GaussianOutput

from rdkit.Chem.rdmolops import AddHs
from rdkit.Chem import MolFromSmiles, AllChem

file_path = os.path.dirname(os.path.realpath(__file__))

# Copyright 2021, University of Kentucky


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


def get_hash_id(identifier, filepath, calc_type, output_file=None):
    metadata = {"calculation_type": calc_type, "mol_file": filepath, "wtuning_output": output_file}
    process_data = ProcessDFT(_id=identifier, metadata=metadata)
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
        data_dict['tuning_parameter'] = tuning_parameter
    if solvent:
        dielectric_constant = {"acetonitrile": 35.688}
        data_dict['solvent'] = {
            'name': solvent,
            'model': 'implicit_solvent',
            'dielectric_constant': dielectric_constant.get(solvent)
        }
    return data_dict


def start_from_smiles(identifier):
    """
    get input geometry frontend DB smiles
    :param identifier: str
    :return: pymatgen mol object
    """
    try:
        response = RESTAPI(method='get', endpoint="restapi/molecules/_id={}/mol_info.init_structure=1".format(identifier),
                           url="https://d3tales.as.uky.edu", return_json=True).response[0]
        return Molecule.from_str(response.get("mol_info", {}).get('init_structure'), 'xyz')
    except:
        response = RESTAPI(method='get', endpoint="restapi/molecules/_id={}/mol_info.smiles=1".format(identifier),
                           url="https://d3tales.as.uky.edu", return_json=True).response[0]
        smiles = response.get("mol_info", {}).get('smiles', '')
        structure = find_lowest_e_conf(smiles)
        return Molecule.from_str(structure, 'xyz')


def get_groundState(identifier, prop='charge'):
    """
    get input geometry frontend DB smiles
    :param identifier: str
    :param prop: str
    :return: int, ground state charge
    """
    response = RESTAPI(
        method='get', endpoint="restapi/molecules/_id={}/mol_info.groundState_{}=1".format(identifier, prop),
        url="https://d3tales.as.uky.edu", return_json=True
    ).response[0]
    return int(response.get("mol_info", {}).get('groundState_'+str(prop)))


def get_db_geom(input_geometry_hash, error_raise=False):
    """
    get input geometry from DB
    :param input_geometry_hash: str
    :return: pymatgen mol object
    """
    response = RESTAPI(
        method='get', endpoint="restapi/rawdata/computation/_id={}/data.geometry=1".format(input_geometry_hash),
        url="https://d3tales.as.uky.edu", return_json=True
    ).response
    if not response:
        if error_raise:
            raise LookupError("Starting geometry with hash {} not found.".format(input_geometry_hash))
        return None
    dict_sites = response[0].get("data", {}).get('geometry')
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


