# Installation 

## PIP Install
This module can be installed with pip install by running the following code: 

```bash
pip install git+https://github.com/d3tales/d3tales_api.git
```

**Note**: If your python version is greater than 3.8, then you will need to install 
[pymatgen](https://pymatgen.org/) directly into the environment you are working in. It is best to use 
[conda](https://anaconda.org/conda-forge/pymatgen) for this.

## Environment Variables
Note that many operation require a `DB_INFO_FILE` environment variable. This environment variable
should be a path to a JSON file containing connection information for the databases. The keys should 
be the database names such as `frontend`, `backend`, `expflow`, and `fireworks`. An example 
is shown below

```json
{
  "frontend": {
    "host": "mongodb://USERNAME:PASSWORD@DATABASE_IP:DATABASE_PORT/frontend",
    "database": "ui"
  }
}
```

To set this variable run:
```bash
export DB_INFO_FILE=$PWD/db_infos_ex.json
```

REST API (`d3tales_api.restapi.RESTAPI`) will use `UPLOAD_USER` and `UPLOAD_PASS` environment 
variables if `username` and `password` arguments are not provided. These environment variables 
can be set in the same was as the `DB_INFO_FILE` variable. 

