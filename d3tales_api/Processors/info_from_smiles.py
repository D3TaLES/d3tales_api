import base64
import math
import sys, os
import selfies as sf
from PIL import Image
from io import BytesIO
import pubchempy as pcp

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import GetFormalCharge, AddHs
from rdkit.Chem.Descriptors import NumRadicalElectrons
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey
from rdkit.Chem import Draw, MolFromSmiles, MolToSmiles, AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, 'SA_Score'))
import sascorer


def find_lowest_e_conf(smiles, num_conf=50):
    """
    Find the lowest-energy RDKit conformer after force-field minimization.
    MMFF is preferred, with UFF used when MMFF parameters are unavailable.
    :param smiles: str, SMILES string
    :param num_conf: number of conformers to search
    :return: str, xyz coordinates
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES string: {!r}".format(smiles))
    molecule = Chem.AddHs(molecule)

    params = AllChem.ETKDGv3()
    params.randomSeed = 0xD3
    params.numThreads = 1
    conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=num_conf, params=params))
    if not conformer_ids:
        raise ValueError("RDKit could not generate a conformer for {!r}".format(smiles))

    if AllChem.MMFFHasAllMoleculeParams(molecule):
        optimized = AllChem.MMFFOptimizeMoleculeConfs(molecule, maxIters=5000, numThreads=1)
    elif AllChem.UFFHasAllMoleculeParams(molecule):
        optimized = AllChem.UFFOptimizeMoleculeConfs(molecule, maxIters=5000, numThreads=1)
    else:
        raise ValueError("No RDKit MMFF or UFF parameters are available for {!r}".format(smiles))

    converged = [
        (conf_id, result[1]) for conf_id, result in zip(conformer_ids, optimized)
        if result[0] == 0 and math.isfinite(result[1])
    ]
    candidates = converged or [
        (conf_id, result[1]) for conf_id, result in zip(conformer_ids, optimized)
        if math.isfinite(result[1])
    ]
    if not candidates:
        raise ValueError("RDKit force-field optimization produced no finite conformer energies")
    best_conf_id, _ = min(candidates, key=lambda candidate: candidate[1])
    return Chem.rdmolfiles.MolToXYZBlock(molecule, confId=best_conf_id)


def image_to_base64(img):
    """
    Produce base 64 string representation of image

    :param img: python Image object
    :return: str representing the image
    """
    output_buffer = BytesIO()
    img.save(output_buffer, format='PNG')
    byte_data = output_buffer.getvalue()
    base64_str = base64.b64encode(byte_data)
    return base64_str.decode('ascii')


def base64_to_image(base64_str, image_path=None):
    """
    Produce python Image object from base 64 string representation of image

    :param base64_str: base64 string representing image
    :param image_path: path for image to be saved (optional)
    :return: python Image object
    """
    byte_data = base64.b64decode(base64_str)
    image_data = BytesIO(byte_data)
    img = Image.open(image_data)
    if image_path:
        img.save(image_path)
    return img


class GenerateMolInfo:
    """
    Generate json object for insertion from smiles string
    Copyright 2021, University of Kentucky

    :param names: list of names for molecule
    :param smiles: smiles string
    :param origin_group: which group the molecule comes from
    :return: mol_info class object
    """
    def __init__(self, smiles, origin_group="", names=None, extra_info=True, ground_spin=None, ground_charge=None):
        self.smiles = smiles
        self.origin_group = origin_group
        self.names = names or []
        self.extra_info = extra_info
        self.ground_spin = ground_spin
        self.ground_charge = ground_charge
        self.mol_info_dict = self.get_mol_info()

    def get_mol_info(self):
        """
        Get molecule information

        :return: mol_info as dict
        """
        mol_info = {}
        # Generate rdkit mol and final (cleaned) smiles
        rdkmol = MolFromSmiles(self.smiles)
        clean_smile = MolToSmiles(rdkmol)
        rdkmol_hs = AddHs(rdkmol)
        AllChem.EmbedMolecule(rdkmol_hs)
        try:
            pcpmol = pcp.get_compounds(clean_smile, namespace="smiles")[0]
            mol_info["iupac_name"] = str(pcpmol.iupac_name)
            synonyms = pcpmol.synonyms
        except pcp.BadRequestError:
            synonyms = []

        # Populate class
        mol_info["smiles"] = clean_smile
        mol_info["selfies"] = sf.encoder(clean_smile)
        if self.origin_group:
            mol_info["source_group"] = self.origin_group
        mol_info["inchi"] = MolToInchi(rdkmol)
        mol_info["inchi_key"] = MolToInchiKey(rdkmol)
        mol_info["molecular_formula"] = CalcMolFormula(rdkmol)
        mol_info["number_of_atoms"] = Mol.GetNumAtoms(rdkmol)
        mol_info["molecular_weight"] = CalcExactMolWt(rdkmol)
        mol_info["groundState_charge"] = self.ground_charge or GetFormalCharge(rdkmol)
        mol_info["groundState_spin"] = self.ground_spin or NumRadicalElectrons(rdkmol) + 1  # calculate spin multiplicity with Hand's rule
        mol_info["sa_score"] = round(sascorer.calculateScore(rdkmol), 2)
        if self.extra_info:
            mol_info["d2_image"] = image_to_base64(Draw.MolToImage(rdkmol))
            mol_info["init_structure"] = find_lowest_e_conf(clean_smile)
        try:
            mol_info["synonyms"] = self.names + synonyms
        except TypeError:
            mol_info["synonyms"] = self.names

        return mol_info
