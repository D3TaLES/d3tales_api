import pytest
pytestmark=[pytest.mark.integration, pytest.mark.database]

def test_insert_parameters_contract():
    pytest.skip("Requires legacy parameter payload plus exact insertion API")
