import os
from d3tales_api.Processors.d3tales_parser import ProcessCV

cv_file = os.path.join(os.getcwd(), 'raw_data', 'cv_test.csv')
metadata = {'experiment_run_id': '71236a03-5560-4d2d-a6f7-d0e4c45c89bd',
            'working_electrode_surface_area': {'value': 4, 'unit': 'cm^2'}, 'molecule_id': '06TNKR',
            'instrument': 'electrochemistry__chi_660d,_chi_1100b,_pine_wavenow',
            'electrode_counter': 'Platinum wire electrode', 'electrode_reference': 'Silver/Silver electrode',
            'electrode_working': 'Polycrystalline platinum electrode',
            'solvent': [{'name': 'Acetonitrile', 'purity': ''}], 'electrolyte': [{'name': 'TBAP', 'purity': ''}],
            'ionic_liquid': []}
metadata = {"data_type": "cv",
            "electrode_counter": "standard_hydrogen_electrode",
            "electrode_reference": "standard_hydrogen_electrode",
            "electrode_working": "standard_hydrogen_electrode",
            "electrolyte": [
                "test"
            ],
            "instrument": "chi_660d,_chi_1100b,_pine_wavenow",
            "ionic_liquid": [
                "test"
            ],
            "redox_mol_concentration": "0.1 M",
            "solvent": [
                "test"
            ],
            "temperature": "273 K",
            "working_electrode_surface_area": "0.05 cm^2"}

_id = "06TNKR"
instance = ProcessCV(cv_file, _id=_id, metadata=metadata).data_dict
# print(instance.get("data"))
# BackDB(collection_name='experimentation', instance=instance)
