"""
Testing of diffusion calculator for ferricynide

"""


import os
from d3tales_api.Calculators.calculators import *
from d3tales_api.Processors.parser_echem import ProcessChiCV

num_electrons = 1
temperature = "298 K"
redox_mol_concentration = "0.011 M"
working_electrode_surface_area = "0.070685835 cm^2"
curve_type = "anodic"

# Iterate through each file in the current directory
multi_data = []
cv_dir = os.path.abspath(os.path.join(os.path.join(os.getcwd(), os.pardir), 'raw_data', 'cv_data'))
for f in [f for f in os.listdir(cv_dir)]:
    print("Parsing file ", f, "...")
    cv_i_data = ProcessChiCV(filepath=os.path.join(cv_dir, f), _id='test').data_dict
    multi_data.append(cv_i_data)

# Iterate through each data item, setting global variables for each one
processed_data = []
for i in multi_data:
    data_dict = i.get('data')
    data_dict["n"] = num_electrons
    data_dict["temperature"] = temperature
    data_dict["redox_mol_concentration"] = redox_mol_concentration
    data_dict["working_electrode_surface_area"] = working_electrode_surface_area
    data_dict.get("reverse", {})[num_electrons - 1][1]
    processed_data.append(data_dict)

# Establish the D3TaLES API Calculator connector
connector = {
    "middle_scan": "middle_sweep",
    "v": "conditions.scan_rate",
    "X": "peak_splittings.{}".format(num_electrons - 1),
    "n": "n",
    "D": "diffusion",
    "T": "temperature",
    "C": "redox_mol_concentration",
    "A": "working_electrode_surface_area",
}

# Calculate diffusion constant
diffusion_cal = CVDiffusionCalculator(connector=connector)
diffusion_coef = diffusion_cal.calculate(processed_data, sci_notation=True)
assert float(diffusion_coef[0]) == 2.357e-06
print("Diffusion Coefficent (Randles–Ševčík): \t", diffusion_coef[0], ' cm^2/s')
assert float(diffusion_coef[1]) == 1.885e-06
print("Diffusion Coefficent (mean method): \t", diffusion_coef[1], ' cm^2/s')

# Calculate charge transfer rates
for data_dict in processed_data:
    data_dict["diffusion"] = float(diffusion_coef[1])
transfer_cal = CVChargeTransferCalculator(connector=connector)
transfer_rate = transfer_cal.calculate(processed_data, sci_notation=True)
assert float(transfer_rate) == 7.048e-04
print("Charge Transfer Rate: \t\t\t", transfer_rate, ' cm/s')

