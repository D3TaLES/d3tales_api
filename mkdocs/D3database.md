# D3database

The D3database module provides a python interface for users to interact with the D<sup>3</sup>TaLES
database. D<sup>3</sup>TaLES uses a [MongoDB](https://www.mongodb.com/) No-SQL structure, so some  base modules can be adapted to any MongoDB database. This module includes applications to
for data insertion, schema validation, basic data generation, database queries, and 
Python-based interaction with the [D<sup>3</sup>TaLES REST API](https://d3tales.as.uky.edu/docs/restapi.html). 

Full documentation can be found <a href="d3tales_api.D3database.html">here</a>.


## Data Insertion

Use `D3Database` with the unified `d3tales` database and an explicit collection.
The database connection is loaded from the `d3tales` entry in `DB_INFO_FILE`.

```python
from d3tales_api.D3database.d3database import D3Database
from d3tales_api.Processors.info_from_smiles import GenerateMolInfo

smiles = "CC"
mol_info = GenerateMolInfo(smiles).mol_info_dict
database = D3Database(
    database="d3tales",
    collection_name="base",
    instance={"_id": "molecule-id", "mol_info": mol_info},
    validate_schema=False,
)
database.insert("molecule-id")
```

`FrontDB` and `BackDB` remain as deprecated compatibility adapters. They use
the same `d3tales` database, with `base` and `computation` as their respective
default collections, and emit a `DeprecationWarning`. New code should use
`D3Database`.

### Schema Validation

`D3Database` can validate a document against the D3TaLES schema before insertion.
For example, the following invalid string value for `groundState_charge` fails
schema validation:

```python
from d3tales_api.D3database.d3database import D3Database

instance = {"_id": "molecule-id", "mol_info": {"groundState_charge": "string"}}
database = D3Database(
    database="d3tales", collection_name="base", instance=instance, schema_db="frontend"
)
```

### Database Queries

Query molecule and computation collections from the same database:

```python
from d3tales_api.D3database.d3database import D3Database

base = D3Database(database="d3tales", collection_name="base", validate_schema=False)
molecule = base.coll.find_one({"mol_info.smiles": "CC"})

computation = D3Database(
    database="d3tales", collection_name="computation", validate_schema=False
)
records = computation.make_query({"mol_id": "05XICU"})
```

### D<sup>3</sup>TaLES REST API 

The `RESTAPI` class may be used to Pythonically access the D<sup>3</sup>TaLES REST API. This is 
useful for searching and scraping the D<sup>3</sup>TaLES database on a larger scale. 

Documentation for the REST API syntax and URL interaction with the D<sup>3</sup>TaLES REST API can be found
[here](https://d3tales.as.uky.edu/docs/restapi.html). 

The following examples show how one might (1) get 200 of the SMILES from the D<sup>3</sup>TaLES
database or (2) fnd the ground state charge for the molecule with ID `05XICU`.  
```python
from d3tales_api.D3database.restapi import RESTAPI

# Get 200 of the SMILES from the D<sup>3</sup>TaLES
endpoint="restapi/molecules/{}/mol_info.smiles=1/limit=200/"
response_1 = RESTAPI(method='get', endpoint=endpoint, url="https://d3tales.as.uky.edu", return_json=True).response

# Fnd the ground state charge for the molecule with ID 05XICU
endpoint="restapi/molecules/_id=05XICU/mol_info.groundState_charge=1"
response_2 = RESTAPI(method='get', endpoint=endpoint, url="https://d3tales.as.uky.edu", return_json=True).response
```

Alternatively, you can use the `D3talesData` class, which produces a pandas DataFrame. 
```python
from d3tales_api.D3database.restapi import D3talesData

# Create data collector object 
data_collector = D3talesData(username='USERNAME', password='PASSWORD')
# Get 10 oxidation potentials 
data_collector.get_prop_data('mol_characterization.oxidation_potential.0.value', limit=10)
```

This class will also gather all `oxidation_potential` properties in the database and produce a one-dimensional 
histogram of the resulting values. Note `USERNAME` and `PASSWORD` should be the user's username and password. 
```python
from d3tales_api.D3database.restapi import D3talesData

# Create data collector object 
data_collector = D3talesData(username='USERNAME', password='PASSWORD')
# Plot 1 dimensional histogram
data_collector.hist_1d('mol_characterization.oxidation_potential.0.value', min_cutoff=-10, max_cutoff=10)
```
Likewise, this example will gather all `molecular_weight` and associated `globular_volume` properties in the 
database and produce a two-dimensional histogram of the resulting values. 
```python
from d3tales_api.D3database.restapi import D3talesData

# Create data collector object 
data_collector = D3talesData(username='USERNAME', password='PASSWORD')
# Plot 1 dimensional histogram
data_collector.hist_2d("species_characterization.groundState.globular_volume.0.value", "mol_info.molecular_weight")
```

There is also a Colab notebook with examples of how to access the D3TaLES data with these models: 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1oK1hOZs0rTpc_SoSFg54qQA5U4Qekqu8?usp=sharing)

Note that the D<sup>3</sup>TaLES REST API is user restricted. This means that you must have a
D<sup>3</sup>TaLES website account and user permission to access this tool. When using the 
D<sup>3</sup>TaLES API `RESTAPI` class, you must include your D<sup>3</sup>TaLES website 
username and password as keyword arguments. Alternatively, you may establish your 
username and password as the environment variables `UPLOAD_USER` `UPLOAD_USER`, respectively. 
The `RESTAPI` class will automatically check these environment variables if no keyword 
arguments are provided. 

If you do not have a D<sup>3</sup>TaLES website account, you may create one 
[here](https://d3tales.as.uky.edu/register/). 

If you have an account but do not have permission to access the REST API, you 
may request permission [here](https://d3tales.as.uky.edu/request-permission/). 