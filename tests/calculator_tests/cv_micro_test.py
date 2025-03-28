import json
import os
from d3tales_api.Processors.parser_echem import ProcessChiCV
from d3tales_api.Calculators.calculators import CVChargeTransferCalculatorMicro, CVDiffusionCalculatorMicro


cv_file = os.path.join(os.getcwd(), '../raw_data', 'cv_data', 'chamini_CV_microelectrode.csv')
metadata = {"data_type": "cv",
               "electrode_counter": "pt_working_electrode",
               "electrode_reference": "ag_reference_electrode",
               "electrode_working": "pt_counter_electrode",
               "electrolyte": ["0.1M TBAPF6"],
               "instrument": "robot_potentiostat",
               "redox_mol_concentration": "0.01 M",
               "solvent": ["ACN"],
               "e_ref": "0.3 V",
               "temperature": "298 K",
               "working_electrode_radius": "0.0005cm"}

_id = "06TNKR"
instance = ProcessChiCV(cv_file, _id=_id, metadata=metadata, micro_electrodes=True).data_dict

instance.update(dict(num_electrodes=1))
connector = {
    "i_ss": "data.i_ss",
    "r": "data.conditions.working_electrode_radius",
    "C": "data.conditions.redox_mol_concentration",
    "n": "num_electrodes",
    "e_half": "data.e_half",
    "e_ref": "data.e_ref",
    "T": "data.conditions.temperature",
    "D": "D",
}

# diff_calc = CVDiffusionCalculatorMicro(connector=connector)
# diff = diff_calc.calculate(instance)
# print(diff)
#
# metadata_dict = CV2Front(backend_data=[instance], run_anodic=False, insert=False,
#                          micro_electrodes=True).meta_dict
# [print(k, v, '\n') for k, v in metadata_dict.items()]

# instance = {
#           "i_ss": "0.39e-6 A",
#           "r": "5e-6 m",
#           "C": "50 mM",
#           "n": 1,
#
#           "e_half": "0.26 V",
#           "e_ref": "0.3 V",
#           "T": "298 K",
#       }
# with open("parsed_micro_data.json", 'w') as fn:
#     json.dump(instance, fn)

print(instance["data"].keys())
print(instance["data"]["i_ss"])
diffusion_coef = CVDiffusionCalculatorMicro(connector=connector).calculate(instance)

instance.update(dict(D=diffusion_coef))
transfer_rate = CVChargeTransferCalculatorMicro(connector=connector).calculate(instance, sci_notation=True)
#
print(diffusion_coef, transfer_rate)

print("CV TESTING SUCCESSFUL")
