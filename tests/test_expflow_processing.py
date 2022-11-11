import json
from d3tales_api.Processors.expflow_parser import *

with open("raw_data/wf_test_3_workflow.json") as fn:
    data = json.load(fn)

expflow_obj = data.get("experiments")[0]

x = ProcessExpFlowObj(expflow_obj, 'Risko')
print(x)
