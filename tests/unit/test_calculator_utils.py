import pint
import pytest

from d3tales_api.Calculators.utils import (
    float_from_str,
    get_electrode_potential,
    rgetkeys,
    unit_conversion,
)


@pytest.mark.unit
@pytest.mark.parametrize(
    "measurement, default_unit, expected",
    [
        ("1e-3", "M", 0.001),
        ("2.5V", "V", 2.5),
        ({"value": "1e-3", "unit": "M"}, "M", 0.001),
        ({"value": 500, "unit": "mV"}, "V", 0.5),
    ],
)
def test_unit_conversion_handles_scientific_and_dict_values(measurement, default_unit, expected):
    assert unit_conversion(measurement, default_unit) == pytest.approx(expected)


@pytest.mark.unit
def test_unit_conversion_accepts_pint_quantity():
    ureg = pint.UnitRegistry()
    assert unit_conversion(ureg.Quantity(250, "mV"), "V") == pytest.approx(0.25)


@pytest.mark.unit
def test_unit_conversion_return_dict():
    assert unit_conversion("500mV", "V", return_dict=True) == {
        "value": pytest.approx(0.5),
        "unit": "V",
    }


@pytest.mark.unit
def test_float_from_str_normalizes_unicode_and_scientific_notation():
    assert float_from_str("1.25 ×10⁻³") == pytest.approx(0.00125)


@pytest.mark.unit
def test_rgetkeys_reads_list_and_nested_dict_indexes():
    data = {"scans": [{"peak": 0.25}]}
    assert rgetkeys(data, "scans.0.peak") == 0.25


@pytest.mark.unit
def test_get_electrode_potential_rejects_unknown_name():
    with pytest.raises(ValueError, match="not found"):
        get_electrode_potential("unknown-electrode")
