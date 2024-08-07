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
  * Output file from `chi_660d, chi_1100b, pine_wavenow` instruments
* UV/Vis Spectroscopy (in development) 
  * Excell files
* Literature Articles 
  * Generate article metadata from DOI


## Basic Data Generation

The `GenerateMolInfo` module can be used to generate several 2D molecular descriptors 
from a SMILES string input. The descriptors can be accessed as attributes of the resulting
class or accessed as a dictionary with the `mol_info_dict` attribute. These descriptors include
`smiles`, `selfies`, `inchi`, `inchi_key`, `iupac_name`, `synonyms`, `init_structure` (coordinates 
for estimated 3D geometry), `molecular_formula`, `groundState_charge`, `number_of_atoms`, 
`molecular_weight`, `d2_image` (bit string for molecule image), `source_group`, and `groundState_spin`.

```python
from d3tales_api.Processors.info_from_smiles import GenerateMolInfo

# Generate basic molecule information using the GenerateMolInfo module
smiles = "CC"
instance = GenerateMolInfo(smiles)
print(instance.mol_info_dict)
```

## Processing Molecular DFT

For this example, we must get an example log file. Here, we pull a [Gaussian](https://gaussian.com/)
log file from GitHub [D<sup>3</sup>TaLES API repo](https://github.com/D3TaLES/d3tales_api).
```python
import shutil
from urllib import request

# Pull example file form GitHub
dft_url = "https://raw.githubusercontent.com/D3TaLES/d3tales_api/main/tests/raw_data/05XICU/logs/05XICU_opt_groundState_3H0.log"
with request.urlopen(dft_url) as response, open("gaussian_ex.log", 'wb') as out_file:
    shutil.copyfileobj(response, out_file)
    
```

To parse the DFT file, simply import the `parser_dft` module and use the `ProcessGausLog`
module to parse the data file (in this case, `gaussian_ex.log`).  Here we chose the `ProcessGausLog`
parsing class because our data file is Gaussian format. 

```python    
from d3tales_api.Processors.parser_dft import ProcessGausLog

dft_data = ProcessGausLog(filepath="gaussian_ex.log")
print(dft_data.data_dict)
```

A user can also include `submission_info` and `metadata` in the processing. 
```python    
from d3tales_api.Processors.parser_dft import ProcessGausLog

submission_info = {
    "source" : "Risko",
    "author" : "d3tales@gmail.com",
    "author_email" : "d3tales@gmail.com",
    "upload_time" : "2021-10-01T21:07:29.546377+00:00",
    "file_type" : "zip",
    "data_category" : "computation",
    "data_type" : "gaussian",
    "all_files_in_zip" : [ 
        "opt_groundState.log", 
        "opt_groundState.fchk"
    ],
    "approved" : True
}

metadata = {
    "id" : "opt_groundState",
    "calculation_type" : "tddft_cation1"
}

dft_data = ProcessGausLog(filepath="gaussian_ex.log", submission_info=submission_info, metadata=metadata)
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

To parse the CV file, simply import the `parser_echem` module and use the `ProcessCV`
module to parse the data file (in this case, `cv_test.txt`).  Here we chose the `ProcessChiCV`
parsing class because our data file comes from a CHI instrument. 
```python    
from d3tales_api.Processors.parser_echem import ProcessChiCV

cv_data = ProcessChiCV(filepath="cv_test.txt", _id='test')
print(cv_data.data_dict)
```

A user can also include `submission_info` and `metadata` in the processing. 
```python    
from d3tales_api.Processors.parser_echem import ProcessChiCV

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

cv_data = ProcessChiCV(filepath="cv_test.txt", _id='test', submission_info=submission_info, metadata=metadata)
print(cv_data.data_dict)
```




## Processing Literature Articles

Here we generate metadata for an article based on its DOI. To do this, we 
simply import the `parser_nlp` module and use the `ProcessNlp`
module to generate the article metadata. 
```python    
from d3tales_api.Processors.parser_nlp import ProcessNlp

doi = "10.1039/c7cs00569e"
nlp_data = ProcessNlp(doi)
print(nlp_data.data_dict)
```

We can also instruct the parser to download the article PDF with the `article_download` kwarg and 
store it in a directory specified with the `download_dir` kwarg. The parser will record the PDF 
location in the output data. 
```python    
from d3tales_api.Processors.parser_nlp import ProcessNlp

doi = "10.1039/c7cs00569e"
nlp_data = ProcessNlp(doi, article_download=True, download_dir="temp/")
print(nlp_data.data_dict)
```
