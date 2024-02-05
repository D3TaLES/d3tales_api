from d3tales_api.Processors.d3tales_parser import *
from d3tales_api.Processors.back2front import CV2Front

ca_file = os.path.join(os.getcwd(), '../raw_data', 'cv_data', 'CV_microelectrode.csv')
ca_metadata = {"data_type": "cv",
               "electrode_counter": "pt_working_electrode",
               "electrode_reference": "ag_reference_electrode",
               "electrode_working": "pt_counter_electrode",
               "electrolyte": ["0.1M TBAPF6"],
               "instrument": "robot_potentiostat",
               "redox_mol_concentration": "0.1 M",
               "solvent": ["ACN"],
               "temperature": "273 K",
               "working_electrode_radius": "0.001cm"}

_id = "06TNKR"
instance = ProcessCVMicro(ca_file, _id=_id, metadata=ca_metadata).data_dict

instance.update(dict(num_electrodes=1))
connector = {
    "i_ss": "data.i_ss",
    "r": "data.conditions.working_electrode_radius",
    "C": "data.conditions.redox_mol_concentration",
    "n": "num_electrodes",
}

diff_calc = CVDiffusionCalculatorMicro(connector=connector)
diff = diff_calc.calculate(instance)
print(diff)

metadata_dict = CV2Front(backend_data=[instance], run_anodic=False, insert=False,
                         micro_electrodes=True).meta_dict
[print(k, v, '\n') for k, v in metadata_dict.items()]

print("CV TESTING SUCCESSFUL")
