import copy
import pytest
from d3tales_api.Calculators.calculators import CVDiffusionCalculator

@pytest.mark.unit
@pytest.mark.electrochemistry
def test_diffusion_legacy_series(ferricyanide_records, cv_connector):
    value = CVDiffusionCalculator(connector=cv_connector).calculate(ferricyanide_records)
    values = value if isinstance(value, (list, tuple)) else [value]
    assert values
    assert all(float(v) > 0 for v in values)

@pytest.mark.unit
@pytest.mark.electrochemistry
def test_diffusion_is_deterministic(ferricyanide_records, cv_connector):
    calc = CVDiffusionCalculator(connector=cv_connector)
    assert calc.calculate(copy.deepcopy(ferricyanide_records)) == calc.calculate(copy.deepcopy(ferricyanide_records))

@pytest.mark.unit
@pytest.mark.electrochemistry
def test_diffusion_does_not_mutate_records(ferricyanide_records, cv_connector):
    original = copy.deepcopy(ferricyanide_records)
    CVDiffusionCalculator(connector=cv_connector).calculate(ferricyanide_records)
    assert ferricyanide_records == original
