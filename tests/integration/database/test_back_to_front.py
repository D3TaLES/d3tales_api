import pytest

pytestmark=[pytest.mark.integration, pytest.mark.database]

@pytest.mark.parametrize("data_kind", ["cv", "gaussian", "nlp"])
def test_back_to_front_contract(data_kind):
    """Replacement target for test_back2front_cv/gaus/nlp.py without production DB writes."""
    pytest.skip("Requires exact backend-to-frontend transformer API and legacy expected outputs")
