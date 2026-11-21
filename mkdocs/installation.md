# Installation

Install the extra that matches the modules you plan to use:

```bash
pip install "d3tales_api[workflows]"
pip install "d3tales_api[workflows,database]"
pip install "d3tales_api[calculators,processors,database]"
```

The `workflows` extra supports local execution. Add `database` when workflows
use the D3TaLES REST API or MongoDB.

For the full environment and tests:

```bash
pip install "d3tales_api[all,test]"
```

PySCF is optional and can be added to a workflow installation with
`pip install "d3tales_api[workflows,pyscf]"`. Gaussian and ORCA executables are
installed separately and made available on the FireWorks workers. The
`requirements.txt` file is retained for existing environments.

NLP support is optional: install `pip install "d3tales_api[nlp]"`. The
ChemDataExtractor2 parser dependency may need to be installed separately
because of its dependency constraints.

## Database configuration

Database operations require `DB_INFO_FILE`, a path to a JSON file containing
connection information for databases such as `frontend`, `backend`, `expflow`,
and `fireworks`. An example is available in `db_info_ex.json`.

```bash
export DB_INFO_FILE="$PWD/db_infos_ex.json"
```

REST API operations may also use `UPLOAD_USER` and `UPLOAD_PASS`. NLP
processing uses `ELS_APIKEY`; installing ChemDataExtractor2 may require a
separate step because of its dependency constraints.
