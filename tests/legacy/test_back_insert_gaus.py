import os
import zipfile
from d3tales_api.D3database.d3database import *
from d3tales_api.Processors.parser_dft import ProcessGausLog

ID='05XICU'
home = os.path.join(os.getcwd(), 'raw_data', "computation_data", ID, "zips")

# for mol in os.listdir(home):
for zip_file in os.listdir(home):
    calculation_type = "_".join(zip_file.split("_")[1:-1])
    print("----------------"+calculation_type+"------------------")
    metadata = {'calculation_type': calculation_type}
    # if calculation_type != 'opt_cation1':
    #     continue
    with zipfile.ZipFile(os.path.join(home, zip_file), "r") as target_zip:
        names = target_zip.namelist()
        if calculation_type == 'wtuning':
            output_file = target_zip.extract(member="output.log")
            mol_file = target_zip.extract(member="wtuning.com")
            metadata["wtuning_output"] = output_file
        else:
            logfile = [f for f in names if ".log" in f][0]
            mol_file = target_zip.extract(member=logfile)
        metadata["mol_file"] = mol_file
        data = ProcessGausLog(_id=ID, metadata=metadata).data_dict
        BackDB(collection_name='computation', instance=data)
        print(data.get("_id"))
