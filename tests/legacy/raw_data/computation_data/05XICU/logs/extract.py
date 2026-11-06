import os
import zipfile

for zip_file in [f for f in os.listdir(os.getcwd()) if f.endswith('.zip')]:
    zipdata = zipfile.ZipFile(zip_file)
    zipinfos = zipdata.infolist()

    # iterate through each file
    for zipinfo in zipinfos:
        if zipinfo.filename.endswith(".log"):
            # This will do the renaming
            zipinfo.filename = zip_file.replace(".zip", ".log")
            zipdata.extract(zipinfo)
