"""
this example shows the usages of calculators with schema2class
Easy to work with for objects with different schema
The code of calculator hardly changes

"""

from d3tales_api.Calculators.calculators import *
from d3tales_api.D3database.d3database import Expflow


NUM_ELECTRONS = 1

cv_entries = Expflow().make_query({"mol_id": "06QMKG"}, output='json')
[d.update({"num_electrodes": NUM_ELECTRONS}) for d in cv_entries]

connector = {
    "A": "data.conditions.working_electrode_surface_area",
    "v": "data.conditions.scan_rate",
    "C": "data.conditions.redox_mol_concentration",
    "n": "num_electrodes",

    "scan_data": "data.scan_data"
}


diffusion_cal = CVDiffusionCalculator(connector=connector)
diffusion = diffusion_cal.calculate(cv_entries)
assert diffusion == [2.969e-06, 2.52e-06]
print("Average diffusion", diffusion[0])
print("Fitted diffusion", diffusion[1])

descriptor_cal = CVDescriptorCalculator(connector=connector)
assert descriptor_cal.peaks(cv_entries[0]) == {'forward': [[0.091, 4.796e-05]], 'reverse': [[-0.05, -6.478e-05]]}
assert descriptor_cal.reversibility(cv_entries[0]) == ['quasi-reversible']
assert descriptor_cal.e_half(cv_entries[0])[0] == [0.02]
assert descriptor_cal.peak_splittings(cv_entries[0]) ==[0.141]
assert len(descriptor_cal.middle_sweep(cv_entries[0])) == 2

print("CV TESTING SUCCESSFUL")
