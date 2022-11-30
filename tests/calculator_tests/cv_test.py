"""
this example shows the usages of calculators with schema2class
Easy to work with for objects with different schema
The code of calculator hardly changes

"""

import python_jsonschema_objects as pjs
from urllib import request
import json, time
import jsonschema
from d3tales_api.Calculators.calculators import *


class Schema2Class():
    def __init__(self, schema_name=None, database=None, named_only=False):

        # fetch schema
        self.schema_name = schema_name
        self.database = database
        schema_url = "https://raw.githubusercontent.com/D3TaLES/schema/main/schema_{}/{}.schema.json".format(
            self.database, self.schema_name)
        response = request.urlopen(schema_url)
        self.schema = json.loads(response.read().decode())

        # generating classes
        builder = pjs.ObjectBuilder(self.schema)
        ns = builder.build_classes(named_only=named_only)

        # get all name space
        for name_space in dir(ns):
            setattr(self, name_space, getattr(ns, name_space))

        # highest-level name space for validation
        # self.high_level = getattr(ns,schema_name.title())

        # required values
        if self.schema.get("required", ):
            self.required = self.schema.get("required")
        else:
            self.required = None


from pymongo import MongoClient

entry = MongoClient("mongodb://10.33.30.17:23771/backend?retryWrites=true&w=majority")["backend"][
    "experimentation"].find_one({"_id": "37dde7f6-204d-4e04-a451-c99a58485f69"})

s2c = Schema2Class(schema_name="experimentation", database="backend")

exp_obj = s2c.Experimentation()

# schema errors
entry["data"].pop("date_recorded")
entry["data"].pop("plot_data")

# added for demo
# entry["data"].update({"ip_cathode": -3.071e-05})
# entry["data"].update({"conc": {"value":1e-5, "unit":"mol/cm3"} })
# entry["data"].update({"area": 0.0725})
# entry["data"].update({"n_elect": 1})

# can be defined for a schema this is a demo
connector = {
    "i_p": "data.conditions.high_e.value",
    "A": "data.conditions.high_e.value",
    "v": "data.conditions.scan_rate.value",
    "C": "data.data_points_per_scan",
    "n": "data.conditions.num_scans",

    "scan_data": "data.scan_data"
}

exp_obj.data = entry["data"]
# diffusion_cal = CVDiffusionCalculator(connector=connector)
# diffusion = diffusion_cal.calculate([exp_obj])
# print("Average diffusion", diffusion[0])
# print("Fitted diffusion", diffusion[1])

descriptor_cal = CVDescriptorCalculator(connector=connector)
print("Peaks: ", descriptor_cal.peaks(exp_obj))
print("Reversibility: ", descriptor_cal.reversibility(exp_obj))
print("E 1/2: ", descriptor_cal.e_half(exp_obj))
print("Peak Splittings: ", descriptor_cal.peak_splittings(exp_obj))
print("Numer of middle sweeps: ", len(descriptor_cal.middle_sweep(exp_obj)))
