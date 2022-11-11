"""
Testing of diffusion calculator for ferricynide

"""


import pandas as pd
from ..calculators import *
from d3tales_api.D3database.restapi import RESTAPI

CSV_TEST = False

if CSV_TEST:
    # test with csv
    df = pd.read_csv("fecn.csv")
    data = df.to_dict('records')
    connector = {
        "i_p": "Current_anode",
        "A": "Anodic",
        "v": "scan rate",
        "C": "conc",
        "n": "n"
    }
else:
    # test with REST API
    processing_ids = [
        "b5b51e31-8fce-42f2-a22e-827f80124b5d",
        "d478d3f2-7ee2-40fb-af91-5e2d781940b4",
        "2c544b8e-581a-4a6b-897c-415cde354cd4",
        "59606dd9-f960-48fb-89b2-10cd9e1a7d8d",
        "575c5005-5acb-4366-955c-ea3ae412ca65"
    ]
    num_electrons = 1
    data = []
    for i in processing_ids:
        query = RESTAPI(method='get', endpoint="restapi/rawdata/experimentation/_id="+i, url="https://d3tales.as.uky.edu",
                        login_endpoint='login', params=dict(approved='on'), return_json=True).response
        if query:
            data.append(query[0].get("data"))
    for data_dict in data:
        data_dict["current"] = data_dict.get("forward", {})[num_electrons-1][1]
        data_dict["n"] = num_electrons
    connector = {
        "i_p": "current",
        "A": "conditions.working_electrode_surface_area",
        "v": "conditions.scan_rate",
        "C": "conditions.redox_mol_concentration",
        "n": "n"
    }
# 7.7*10^-6 (anodic)

diffusion_cal = CVDiffusionCalculator(connector=connector)
diffusion = diffusion_cal.calculate(data)
print("Average diffusion", diffusion[0])
print("Fitted diffusion", diffusion[1])

if not CSV_TEST:
    for data_dict in data:
        data_dict["diffusion"] = diffusion[1]
        data_dict["temperature"] = 293

    t_connector = {
            "T": "temperature",
            "D": "diffusion",
            "v": "conditions.scan_rate",
            "n": "n",
            "X": "peak_splittings.{}".format(num_electrons-1)
        }
    # 0.00127 cm/s

    transfer_cal = CVChargeTransferCalculator(connector=t_connector)
    transfer = transfer_cal.calculate(data)
    print("Charge Transfer Rate", transfer)
