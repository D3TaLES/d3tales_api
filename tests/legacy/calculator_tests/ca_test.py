import os
from d3tales_api.Processors.parser_echem import ProcessChiCA


ca_file = os.path.join(os.getcwd(), '../raw_data', 'ca_data', 'CA.csv')
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

# _id = "06TNKR"
# instance = ProcessCA(ca_file, _id=_id, metadata=ca_metadata)
# print(instance.data_dict)

for fn in os.listdir(os.path.join(os.getcwd(), '../raw_data', 'ca_data')):
    f_path = os.path.join(os.path.join(os.getcwd(), '../raw_data', 'ca_data'), fn)
    if not f_path.endswith(".csv"):
        continue
    resistance = ProcessChiCA(f_path, metadata=ca_metadata).measured_resistance
    print("{}: \t{:.2f}".format(fn.strip(".csv"), resistance))


print("CA TESTING SUCCESSFUL")
