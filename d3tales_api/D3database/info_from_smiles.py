import base64
import selfies as sf
from PIL import Image
from io import BytesIO
import pubchempy as pcp
from d3tales_api.D3database.schema2class import Schema2Class

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import GetFormalCharge, AddHs
from rdkit.Chem.Descriptors import NumRadicalElectrons
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey
from rdkit.Chem import Draw, MolFromSmiles, MolToSmiles, AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt


class GenerateMolInfo:
    """
        Generate json object for insertion from smiles string
        :param names: list of names for molecule
        :param smiles: smiles string
        :param origin_group: which group the molecule comes from
        :return: mol_info class object
        Copyright 2021, University of Kentucky
        """
    def __init__(self, smiles, origin_group="", names=[], extra_info=True, database='frontend', schema_name="mol_info"):
        self.smiles = smiles
        self.origin_group = origin_group
        self.names = names
        self.database = database
        self.schema_name = schema_name
        self.extra_info = extra_info
        self.mol_info_dict = self.get_mol_info()

    def get_mol_info(self):
        # Fetch schema and build class
        s2c = Schema2Class(schema_name=self.schema_name, database=self.database)
        mol_info = s2c.MolInfo()
        # Generate rdkit mol and final (cleaned) smiles
        rdkmol = MolFromSmiles(self.smiles)
        clean_smile = MolToSmiles(rdkmol)
        rdkmol_hs = AddHs(rdkmol)
        AllChem.EmbedMolecule(rdkmol_hs)
        pcpmol = pcp.get_compounds(clean_smile, namespace="smiles")[0]

        # Populate class
        mol_info.smiles = clean_smile
        mol_info.selfies = sf.encoder(clean_smile)
        if self.origin_group:
            mol_info.source_group = self.origin_group
        mol_info.inchi = MolToInchi(rdkmol)
        mol_info.inchi_key = MolToInchiKey(rdkmol)
        mol_info.iupac_name = str(pcpmol.iupac_name)
        mol_info.molecular_formula = CalcMolFormula(rdkmol)
        mol_info.number_of_atoms = Mol.GetNumAtoms(rdkmol)
        mol_info.molecular_weight = CalcExactMolWt(rdkmol)
        mol_info.groundState_charge = GetFormalCharge(rdkmol)
        mol_info.groundState_spin = NumRadicalElectrons(rdkmol) + 1  # calculate spin multiplicity with Hand's rule
        if self.extra_info:
            mol_info.d2_image = self.image_to_base64(Draw.MolToImage(rdkmol))
            mol_info.init_structure = AllChem.rdmolfiles.MolToXYZBlock(rdkmol_hs)
        try:
            mol_info.synonyms = self.names + pcpmol.synonyms
        except TypeError:
            mol_info.synonyms = self.names

        return mol_info.as_dict()

    @staticmethod
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

    @staticmethod
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
