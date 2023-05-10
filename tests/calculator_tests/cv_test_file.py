import os
from d3tales_api.Calculators.plotters import *
from d3tales_api.Calculators.calculators import *
from d3tales_api.Processors.d3tales_parser import *

SELF_STD = True
NUM_ELECTRONS = 1

ROBOT_DATA = True
DEFAULT_CONCENTRATION = "0.01M"
DEFAULT_TEMPERATURE = "293K"
DEFAULT_WORKING_ELECTRODE_AREA = "0.03142cm^2"

if ROBOT_DATA:
    cv_dir = "../raw_data/robotic_data"
    cv_locations = [os.path.join(cv_dir, f) for f in os.listdir(cv_dir) if f.endswith(".csv")]
    cv_entries = [ProcessCV(loc, _id="test", submission_info={}, metadata={"redox_mol_concentration": DEFAULT_CONCENTRATION, "temperature": DEFAULT_TEMPERATURE, "working_electrode_surface_area": DEFAULT_WORKING_ELECTRODE_AREA}, parsing_class=ParseChiCV).data_dict for loc in cv_locations]
else:
    cv_data = ProcessCV("../raw_data/aman_cv2.csv", _id='test', parsing_class=ParseChiCV).data_dict
    cv_entries = [cv_data]

[d.update({"num_electrons": NUM_ELECTRONS}) for d in cv_entries]
[d.update({"current_cathodic": d.get("data", {}).get("forward", [[]])[NUM_ELECTRONS - 1][1]}) for d in cv_entries]


connector = {
    "n": "num_electrons",
    "i_p": "current_cathodic",
    "A": "data.conditions.working_electrode_surface_area",
    "v": "data.conditions.scan_rate",
    "C": "data.conditions.redox_mol_concentration",
    "X": "data.peak_splittings.{}".format(NUM_ELECTRONS - 1),
    "T": "data.conditions.temperature",
    "D": "diffusion",

    "scan_data": "data.scan_data",
    "variable_prop": "data.conditions.scan_rate.value",
}



diffusion_cal = CVDiffusionCalculator(connector=connector)
diffusion = diffusion_cal.calculate(cv_entries)
print("Average diffusion", diffusion[0])
print("Fitted diffusion", diffusion[1])
[d.update({"diffusion": diffusion[1]}) for d in cv_entries]
#
# charge_transfer_cal = CVChargeTransferCalculator(connector=connector)
# charge_transfer = charge_transfer_cal.calculate(cv_entries)
# print("Charge Transfer", charge_transfer)

# e_half_transfer_cal = AvgEHalfCalculator(connector=connector)
# e_half = e_half_transfer_cal.calculate(cv_entries)
# print("Avg E half", e_half)

# descriptor_cal = CVDescriptorCalculator(connector=connector)
# print(descriptor_cal.peaks(cv_entries[0]))
# descriptor_cal.reversibility(cv_entries[0])
# descriptor_cal.e_half(cv_entries[0])[0]
# descriptor_cal.peak_splittings(cv_entries[0])
# len(descriptor_cal.middle_sweep(cv_entries[0]))

cv_plotter = CVPlotter(connector=connector).live_plot(cv_entries[0], fig_path="cv_test.png", self_standard=SELF_STD,
                                                                  title="CV Plot", xlabel="x", ylabel='y')
cv_plotter_multi = CVPlotter(connector=connector).live_plot_multi(cv_entries, fig_path="cv_test_multi.png", self_standard=SELF_STD,
                                                                  title="CV Plot", xlabel="Potential (V)", ylabel='Current (A)', legend_title="Scan Rate (V/s)")

print("CV TESTING SUCCESSFUL")
