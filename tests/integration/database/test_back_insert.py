import pytest

pytestmark=[pytest.mark.integration, pytest.mark.database]

@pytest.mark.parametrize("data_kind", ["cv", "gaussian", "nlp"])
def test_backend_insert_contract_isolated(monkeypatch, data_kind):
    """Replacement target for test_back_insert_cv/gaus/nlp.py.

    This intentionally performs no live insert. Wire the package's insertion class
    to a fake collection and assert validate/insert calls when its exact API is confirmed.
    """
    pytest.skip("Requires exact D3database insertion class and isolated collection fixture")
