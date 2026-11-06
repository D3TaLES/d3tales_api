import pytest
from d3tales_api.Calculators.calculators import ConcentrationCalculator

@pytest.mark.unit
def test_concentration_from_legacy_example():
    data = {
        "concentration_smiles":"CC1(CCCC(N1[O])(C)C)C",
        "concentration_mass":"17.19mg",
        "solv_density":"1g/mL",
        "concentration_volume":"10mL",
    }
    connector = {
        "smiles":"concentration_smiles", "weight":"concentration_mass",
        "solv_density":"solv_density", "volume":"concentration_volume",
    }
    result = ConcentrationCalculator(connector=connector).calculate(data)
    assert float(result) == pytest.approx(0.01101, rel=5e-3)

@pytest.mark.unit
def test_concentration_does_not_mutate_input():
    data = {"smiles":"O", "weight":"18.015g", "solv_density":"1g/mL", "volume":"1L"}
    original = dict(data)
    ConcentrationCalculator(connector={k:k for k in data}).calculate(data)
    assert data == original
