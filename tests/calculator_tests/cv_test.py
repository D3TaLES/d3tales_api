"""
this example shows the usages of calculators with schema2class
Easy to work with for objects with different schema
The code of calculator hardly changes

"""

import python_jsonschema_objects as pjs
from urllib import request
from d3tales_api.Calculators.calculators import *
from d3tales_api.D3database.d3database import Expflow


from pymongo import MongoClient

cv_entries = Expflow().make_query({"mol_id": "06QMKG"}, output='json')

# can be defined for a schema this is a demo
connector = {
    "i_p": "data.conditions.high_e.value",
    "A": "data.conditions.high_e.value",
    "v": "data.conditions.scan_rate.value",
    "C": "data.data_points_per_scan",
    "n": "data.conditions.num_scans",

    "scan_data": "data.scan_data"
}

diffusion_cal = CVDiffusionCalculator(connector=connector)
diffusion = diffusion_cal.calculate(cv_entries)
print("Average diffusion", diffusion[0])
print("Fitted diffusion", diffusion[1])

descriptor_cal = CVDescriptorCalculator(connector=connector)
print("Peaks: ", descriptor_cal.peaks(cv_entries[0]))
print("Reversibility: ", descriptor_cal.reversibility(cv_entries[0]))
print("E 1/2: ", descriptor_cal.e_half(cv_entries[0]))
print("Peak Splittings: ", descriptor_cal.peak_splittings(cv_entries[0]))
print("Numer of middle sweeps: ", len(descriptor_cal.middle_sweep(cv_entries[0])))
