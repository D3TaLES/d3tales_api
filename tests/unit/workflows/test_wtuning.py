import pytest

from d3tales_api.Workflows.wtuning import WtuningJob


@pytest.mark.unit
@pytest.mark.parametrize("scheme", ["__import__('os').system('false')", "", None])
def test_unknown_tuning_scheme_is_rejected(scheme):
    with pytest.raises(ValueError, match="scheme must be one of"):
        WtuningJob(scheme=scheme)
