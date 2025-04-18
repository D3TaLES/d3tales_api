"""
Testing of diffusion calculator for ferricynide

"""

from d3tales_api.Calculators.calculators import *
from d3tales_api.D3database.restapi import RESTAPI

REST_API = False

data = []
num_electrons = 1

if REST_API:
    # test with REST API
    processing_ids = [
        "b5b51e31-8fce-42f2-a22e-827f80124b5d",
        "d478d3f2-7ee2-40fb-af91-5e2d781940b4",
        "2c544b8e-581a-4a6b-897c-415cde354cd4",
        "59606dd9-f960-48fb-89b2-10cd9e1a7d8d",
        "575c5005-5acb-4366-955c-ea3ae412ca65"
    ]


    for i in processing_ids:
        query = RESTAPI(method='get', endpoint="restapi/rawdata/experimentation/_id="+i, url="https://d3tales.as.uky.edu",
                        login_endpoint='login', params=dict(approved='on'), return_json=True).response
        if query:
            data.append(query[0].get("data"))
    for data_dict in data:
        data_dict["current"] = data_dict.get("forward", {})[num_electrons-1][1]
        data_dict["n"] = num_electrons
    connector = {
        "A": "conditions.working_electrode_surface_area",
        "v": "conditions.scan_rate",
        "C": "conditions.redox_mol_concentration",
        "n": "n"
    }
else:
    # test with std pre-extracted cv data
    df = pd.read_csv("fecn.csv")
    data = [df.to_dict('records')]
    connector = {
        "A": "Anodic",
        "v": "scan rate",
        "C": "conc",
        "n": "n"
    }


diffusion_cal = CVDiffusionCalculator(connector=connector)
diffusion = diffusion_cal.calculate(data)
print("Average diffusion", diffusion[0])
print("Fitted diffusion", diffusion[1])

if REST_API:
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
