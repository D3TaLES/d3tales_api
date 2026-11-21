# Installation

Install the extra that matches the modules you plan to use:

```bash
pip install "d3tales_api[workflows]"
pip install "d3tales_api[workflows,database]"
pip install "d3tales_api[calculators,processors,database]"
```

The `workflows` extra supports local execution. Add `database` when workflows
use a D3TaLES REST API or direct MongoDB connection.

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

Direct database operations require `DB_INFO_FILE`, a path to a JSON file
containing a `d3tales` connection entry (database name `d3tales`). Other
application databases may have their own entries. An example is available in
`db_info_ex.json`.

```bash
export DB_INFO_FILE="$PWD/db_info_ex.json"
```

Use `d3tales_direct` only on workers inside the database network firewall with
write access to the D3TaLES database. `d3tales_web` uses the REST API and
website upload flow instead; it does not need MongoDB credentials. REST API
operations may also use `UPLOAD_USER` and `UPLOAD_PASS`. NLP
processing uses `ELS_APIKEY`; installing ChemDataExtractor2 may require a
separate step because of its dependency constraints.
