import pytest
calculators = pytest.importorskip("d3tales_api.Calculators.calculators")

@pytest.mark.electrochemistry
def test_cv_diffusion_returns_two_positive_values(ferricyanide_records, cv_connector):
    result = calculators.CVDiffusionCalculator(connector=cv_connector).calculate(ferricyanide_records)
    assert len(result) == 2
    assert all(float(value) > 0 for value in result)

@pytest.mark.electrochemistry
def test_cv_diffusion_does_not_mutate_input(ferricyanide_records, cv_connector):
    import copy
    original = copy.deepcopy(ferricyanide_records)
    calculators.CVDiffusionCalculator(connector=cv_connector).calculate(ferricyanide_records)
    assert ferricyanide_records == original
