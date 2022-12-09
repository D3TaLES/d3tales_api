# D3database

The D3database module provides a python interface for users to interact with the D<sup>3</sup>TaLES
database. D<sup>3</sup>TaLES uses a [MongoDB](https://www.mongodb.com/) No-SQL structure, so some  base modules can be adapted to any MongoDB database. This module includes applications to
for data insertion, schema validation, basic data generation, database queries, and 
Python-based interaction with the [D<sup>3</sup>TaLES REST API](https://d3tales.as.uky.edu/docs/restapi.html). 

Full documentation can be found [here](d3tales_api.D3database.html).

## Data Insertion

```python
from rdkit.Chem import MolFromSmiles, MolToSmiles
from d3tales_api.D3database.d3database import FrontDB

smiles = "CC"
clean_smiles = MolToSmiles(MolFromSmiles(smiles))
FrontDB(smiles=clean_smiles).check_if_in_db()
```
Should return ``. 


inserting...
```python
from rdkit.Chem import MolFromSmiles, MolToSmiles
from d3tales_api.D3database.d3database import FrontDB
from d3tales_api.D3database.info_from_smiles import GenerateMolInfo

smiles = "CC"
clean_smiles = MolToSmiles(MolFromSmiles(smiles))

instance = GenerateMolInfo(clean_smiles, database="frontend").mol_info_dict
db_insertion = FrontDB(schema_layer='mol_info', instance=instance, smiles=clean_smiles, group="Non-D3TaLES")
```

### Basic Data Generation

```python

```

### Schema Validation

```python

```

## Database Queries

```python

```

### D<sup>3</sup>TaLES REST API 

```python

```
