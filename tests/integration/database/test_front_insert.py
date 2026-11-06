import pytest
pytestmark=[pytest.mark.integration, pytest.mark.database]

def test_front_insert_contract():
    pytest.skip("Requires isolated frontend collection and exact insertion API")
