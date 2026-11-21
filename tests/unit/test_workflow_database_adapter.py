import subprocess
import sys

import pytest

from d3tales_api.D3database import workflow_adapter
from d3tales_api.D3database.d3database import DBconnector
from d3tales_api.D3database.restapi import RESTAPI


@pytest.mark.unit
def test_workflow_adapter_reexports_canonical_database_clients():
    assert workflow_adapter.DBconnector is DBconnector
    assert workflow_adapter.RESTAPI is RESTAPI


@pytest.mark.unit
def test_workflow_adapter_import_does_not_load_database_clients():
    code = (
        "import sys; "
        "import d3tales_api.D3database.workflow_adapter; "
        "assert 'pymongo' not in sys.modules; "
        "assert 'requests' not in sys.modules"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
