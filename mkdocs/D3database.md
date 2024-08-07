# D3database

The D3database module provides a python interface for users to interact with the D<sup>3</sup>TaLES
database. D<sup>3</sup>TaLES uses a [MongoDB](https://www.mongodb.com/) No-SQL structure, so some  base modules can be adapted to any MongoDB database. This module includes applications to
for data insertion, schema validation, basic data generation, database queries, and 
Python-based interaction with the [D<sup>3</sup>TaLES REST API](https://d3tales.as.uky.edu/docs/restapi.html). 

Full documentation can be found [here](d3tales_api.D3database.html).


## Data Insertion


The following example show the insertion of generated molecule information (see above)
into the frontend database. Note that information will be inserted into the database 
specified in the `DB_INFO_FILE`.

```python
from d3tales_api.D3database.d3database import FrontDB
from d3tales_api.Processors.info_from_smiles import GenerateMolInfo

# Generate basic molecule information using the GenerateMolInfo module
smiles = "CC"
instance = GenerateMolInfo(smiles, database="frontend").mol_info_dict

# Insert basic molecule information into the Frontend database
db_insertion = FrontDB(schema_layer='mol_info', instance=instance, smiles=smiles, group="Non-D3TaLES")
```


### Schema Validation

The D3database modules automatically validate instances with the [appropriate schema](https://github.com/D3TaLES/schema). 
```python
from d3tales_api.D3database.d3database import FrontDB

# Create instance with an erroneous field
instance = {"groundState_charge": "string"}

# Insert basic molecule information into the Frontend database
db_insertion = FrontDB(_id="05XICU", schema_layer='mol_info', instance=instance)
```
Here, the code tries to insert a `string` instead of an `integer` for the attribute `groundState_charge`. 
However, this violates the frontend database schema. Because the module automatically validates 
the schema when and `_id` is present, this should result in an error something like this: 
```
Failed validating 'type' in schema['properties']['mol_info']['properties']['groundState_charge']:
    {'description': 'Charge of the ground state molecule', 'type': 'number'}

On instance['mol_info']['groundState_charge']:
    'string'
```
Alternatively, this should not produce an error: 
```python
from d3tales_api.D3database.d3database import FrontDB

# Create instance with an erroneous field
instance = {"groundState_charge": 0}

# Insert basic molecule information into the Frontend database
db_insertion = FrontDB(_id="05XICU", schema_layer='mol_info', instance=instance)

```

## Database Queries

One useful function for database queries is the `FrondDB` method `check_if_in_db`. This 
will check if the SMILES string associated with the `FrondDB` already exists in the database. 
If it exists, the method will return the molecule ID, if not it will return `False`. 
```python
from d3tales_api.D3database.d3database import FrontDB

smiles = "CC"
FrontDB(smiles=smiles).check_if_in_db()
```
This should return `06PCFL` if the `DB_INFO_FILE` contains the database information for the 
D<sup>3</sup>TaLES database. If not, the result will vary depending on the database. 


This API can also be used for more general MongoDB queries. The following example show
how a user might query the backend computational database for calculations with the mol_id
of `05XICU`. The query can also be filtered to return only calculation types. 
```python
from d3tales_api.D3database.d3database import BackDB

# Query all computational entries with the mol_id 05XICU
BackDB(collection_name="computation").make_query({"mol_id": "05XICU"})

# Query all computational entries with the mol_id 05XICU and return only calculation_types
BackDB(collection_name="computation").make_query({"mol_id": "05XICU"}, {"calculation_type": 1})
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