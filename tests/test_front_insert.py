import json
from sys import argv
from d3tales_api.D3database.d3database import FrontDB
from d3tales_api.D3database.info_from_smiles import GenerateMolInfo
from rdkit.Chem import MolFromSmiles, MolToSmiles

if __name__ == "__main__":
    if len(argv) > 1:
        smiles_dict_file = argv[1]
        with open(smiles_dict_file) as fn:
            smiles_dict = json.load(fn)
    else:
        smiles_dict = {"water": "CC"}
    try:
        group = argv[2]
    except IndexError:
        group = "Risko"

    for name, smi in smiles_dict.items():
        rdkmol = MolFromSmiles(smi)
        clean_smiles = MolToSmiles(rdkmol)
        instance = GenerateMolInfo(clean_smiles, names=[name], origin_group=group, database="frontend").mol_info_dict
        print(instance)
        insertion = FrontDB(schema_layer='mol_info', instance=instance, smiles=clean_smiles, group=group)
        print("{} inserted into frontend database".format(insertion.id))
