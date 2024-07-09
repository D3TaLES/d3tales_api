import base64
import sys, os
import selfies as sf
from PIL import Image
from tqdm import tqdm
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
    Find the lowest energy conformer for a molecule with RDKit
    :param smiles: str, SMILES string
    :param num_conf: float, number of conformers to serach
    :return: str, xyz coordinates
    """
    rdkmol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    results = {}
    AllChem.EmbedMultipleConfs(rdkmol, numConfs=num_conf, params=AllChem.ETKDG())
    results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(rdkmol, maxIters=5000)
    for i, result in tqdm(enumerate(results_MMFF)):
        results[i] = result[1]
    best_idx = min(results, key=results.get)
    structure = Chem.rdmolfiles.MolToXYZBlock(rdkmol, confId=best_idx)
    return structure


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

