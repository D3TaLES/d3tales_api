import pytest
pytestmark = [pytest.mark.database, pytest.mark.integration]

def test_database_tests_are_opt_in():
    pytest.skip("Provide an isolated MongoDB test URI and implement adapter fixtures before enabling")
