# D<sup>3</sup>TaLES API 
![](D3TaLES_api_logo.png)


API designed for the [D<sup>3</sup>TaLES project](https://d3tales.as.uky.edu/). Includes data processing, database access functions, and property calculators.

## Installation 

### PIP Install
This module can be installed with pip install by running the following code: 

```bash
pip install git+https://github.com/d3tales/d3tales_api.git
```

### Environment Variables
Note that many operation require a `DB_INFO_FILE` environment variable. This environment variable
should be a path to a JSON file containing connection information for the databases. The keys should 
be the database names such as `frontend`, `backend`, `expflow`, and `fireworks`.
```json
{
  "frontend": {
    "host": "mongodb://USERNAME:PASSWORD@DATABASE_IP:DATABASE_PORT/frontend",
    "database": "frontend"
  }
}
```
To set this variable run:
```bash
export DB_INFO_FILE=/home/USER/db_infos.json
```

REST API (`d3tales_api.restapi.RESTAPI`) will use `UPLOAD_USER` and `UPLOAD_PASS` environment 
variables if `username` and `password` arguments are not provided. These environment variables 
can be set in the same was as the `DB_INFO_FILE` variable. 

## Modules 
This API contains three modules: 
* `D3database`: Interface for interacting with the D<sup>3</sup>TaLES databases. Includes 
queries, data insertion with validation, and a REST API interface.
* `Processors`: Processing code for extracting data from computational and experimental
data files as well as [ExpFlow](https://d3tales.as.uky.edu/expflow) run objects. 
* `Calculators`: Property calculators and plotters for cyclic voltammetry
