from d3tales_api.Processors.d3tales_parser import *

cv_file = os.path.join(os.getcwd(), '../raw_data', 'pot_data', 'CA.csv')
ca_metadata = {"data_type": "cv",
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
instance = ProcessCA(cv_file, _id=_id, metadata=ca_metadata)
print(instance.data_dict)
# print(instance.data_dict.get("data"))


connector = {
    "A": "data.conditions.working_electrode_surface_area",
    "v": "data.conditions.scan_rate",
    "C": "data.conditions.redox_mol_concentration",
    "n": "num_electrodes",
    "T": "data.conditions.temperature",
    "D": "diffusion",

    "scan_data": "data.scan_data",
    "variable_prop": "data.conditions.scan_rate.value"
}

# dirty_calc = DirtyElectrodeDetector(connector=connector)
# dirty = dirty_calc.calculate(instance.data_dict)

print("CV TESTING SUCCESSFUL")
