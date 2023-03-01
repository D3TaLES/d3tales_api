from d3tales_api.Calculators.plotters import *
from d3tales_api.Calculators.calculators import *
from d3tales_api.Processors.d3tales_parser import *


NUM_ELECTRONS = 1
cv_data = ProcessCV("../raw_data/aman_cv2.csv", _id='test', parsing_class=ParseChiCV).data_dict
cv_entries = [cv_data]
[d.update({"num_electrodes": NUM_ELECTRONS}) for d in cv_entries]

connector = {
    "A": "data.conditions.working_electrode_surface_area",
    "v": "data.conditions.scan_rate",
    "C": "data.conditions.redox_mol_concentration",
    "n": "num_electrodes",

    "scan_data": "data.scan_data",
    "variable_prop": "data.conditions.scan_rate.value"
}



# diffusion_cal = CVDiffusionCalculator(connector=connector)
# diffusion = diffusion_cal.calculate(cv_entries)
# print("Average diffusion", diffusion[0])
# print("Fitted diffusion", diffusion[1])

descriptor_cal = CVDescriptorCalculator(connector=connector)
print(descriptor_cal.peaks(cv_entries[0]))
# descriptor_cal.reversibility(cv_entries[0])
# descriptor_cal.e_half(cv_entries[0])[0]
# descriptor_cal.peak_splittings(cv_entries[0])
# len(descriptor_cal.middle_sweep(cv_entries[0]))

cv_plotter = CVPlotter(connector=connector).live_plot(cv_entries[0], fig_path="cv_test.png")
cv_plotter_multi = CVPlotter(connector=connector).live_plot_multi(cv_entries, fig_path="cv_test_multi.png")

print("CV TESTING SUCCESSFUL")
