# Processors

This module contains data parsers to process both computational and experimental 
raw data files. The D<sup>3</sup>TaLES parser then converts the processed data to the 
D<sup>3</sup>TaLES schema. It also contains a submodule for converting data from the backend
D<sup>3</sup>TaLES database to frontend data and pushing the new data to the frontend 
database. 

Full documentation can be found [here](d3tales_api.Processors.html).

### Current Parsers: 
* Molecular DFT
  * Gaussian logfiles
  * Psi4 logfiles (in development)
* Cyclic Voltammetry
  * Chi CV files
* UV/Vis Spectroscopy
  * Excel output from a ___ instrument 

## Processing Molecular DFT

```python
import shutil
from urllib import request

# Pull example file form GitHub
dft_url = "https://raw.githubusercontent.com/D3TaLES/d3tales_api/main/tests/raw_data/05XICU/logs/05XICU_opt_groundState_3H0.log"
with request.urlopen(dft_url) as response, open("gaussian_ex.log", 'wb') as out_file:
    shutil.copyfileobj(response, out_file)
    
```
```python    
from d3tales_api.Processors.d3tales_parser import *

dft_data = ProcessDFT(filepath="gaussian_ex.log", parsing_class=ProcessGausLog)
print(dft_data.data_dict)
```

A user can also include `submission_info` and `metadata` in the processing. 
```python    
from d3tales_api.Processors.d3tales_parser import *

submission_info = {}

metadata = {}

dft_data = ProcessDFT(filepath="gaussian_ex.log", submission_info=submission_info, metadata=metadata, parsing_class=ProcessGausLog)
print(dft_data.data_dict)
```

## Processing Cyclic Voltmeter 

For this example, we must get an example CV file. Here, we pull a CV text file from GitHub 
[D<sup>3</sup>TaLES API repo](https://github.com/D3TaLES/d3tales_api). 
```python
import shutil
from urllib import request

# Pull example file form GitHub
cv_url = "https://raw.githubusercontent.com/D3TaLES/d3tales_api/main/tests/raw_data/cv_test.csv"
with request.urlopen(cv_url) as response, open("cv_test.txt", 'wb') as out_file:
    shutil.copyfileobj(response, out_file)
    
```

To parse the CV, file, simply import the `d3tales_parser` module and use teh `ProcessCV`
module to parse the data file (in this case, `cv_test.txt`).  Here we chose the `ParseChiCV`
parsing class because of the format of our data file. 
```python    
from d3tales_api.Processors.d3tales_parser import *

cv_data = ProcessCV(filepath="cv_test.txt", _id='test', parsing_class=ParseChiCV)
print(cv_data.data_dict)
```

A user can also include `submission_info` and `metadata` in the processing. 
```python    
from d3tales_api.Processors.d3tales_parser import *

submission_info = {
    "source" : "Risko",
    "author" : "d3tales@gmail.com",
    "author_email" : "d3tales@gmail.com",
    "upload_time" : "2022-04-14T22:11:23.490652+00:00",
    "file_type" : "csv",
    "data_category" : "experimentation",
    "data_type" : "cv",
    "approved" : False
}

metadata = {
    "electrode_counter" : "standard_hydrogen_electrode",
    "electrode_working" : "standard_hydrogen_electrode",
    "electrode_reference" : "standard_hydrogen_electrode",
    "solvent" : ["acetonitrile"],
    "electrolyte" : ["tetrabutylammonium hexafluorophosphate"],
    "ionic_liquid" : [],
    "instrument" : "chi_660d,_chi_1100b,_pine_wavenow",
    "working_electrode_surface_area" : "0.05 cm^2",
    "temperature" : "273 K",
    "redox_mol_concentration" : "0.1 M",
    "data_type" : "cv"
}

cv_data = ProcessCV(filepath="cv_test.txt", _id='test', submission_info=submission_info, metadata=metadata, parsing_class=ParseChiCV)
print(cv_data.data_dict)
```


