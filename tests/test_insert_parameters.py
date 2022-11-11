from d3tales_api.D3database.d3database import ParamsDB
from d3tales_api.D3database.info_from_smiles import GenerateMolInfo

electrode_dict = {
    "_id": "standard hydrogen electrode",
    "smiles": "[H]",
    "size": {},
    "roughness": {},
    "absolute_potential": {"value": 4.42, "unit": "eV"}
}
ag_electrode_dict = {
    "_id": "silver electrode",
    "smiles": "[Ag]",
    "size": {"value": 5, "unit": "nm"},
    "roughness": {},
    "absolute_potential": {}
}

il_dict = {
    "_id": "1-Benzyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide",
    "purity":  {"value": 99, "unit": "%"},
    "cas_number": "433337-24-7",
    "supplier": "IoLiTec",
    "catalog_num": "IL-0241-HP",
    "viscosity": {"value": 153, "unit": "cP @25C"},
    "density": {"value": 1.49, "unit": "g/cm3"},
    "smiles": "C[N+]1=CN(C=C1)CC2=CC=CC=C2.C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F",
    "conductivity": {"value": 1.77, "unit": "mScm-1 @30C"},
    "electrochemical_window": {"value": 4.5, "unit": "V"}
}
ftir_dict = {
    "_id": "Thermo is50",
    "detectors": ["DTGS/MCT-A"],
    "frequency_range": "650-4000",
    "resolution": {"value": 1.77, "unit": "mScm-1 @30C"},
    "beamsplitter": "KBr",
    "windows": [{"_id": "CaF2"}],
    "attachments": ["Pike VeeMax 3", "Jackfish SEC w/ Si wafers +/- Au Coating", "External bench for IRRAS/PMIRRAS"]
}

data_dict = electrode_dict
data_dict["_id"] = "_".join(data_dict.get('_id').lower().split(' '))
print(data_dict.get("_id"))
mol_info = GenerateMolInfo(smiles=data_dict["smiles"], extra_info=False).mol_info_dict
data_dict.update(mol_info)

ParamsDB(instance=data_dict, collection_name='electrode', schema_directory='materials', update=False)
