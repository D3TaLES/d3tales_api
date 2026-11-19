import pytest

from d3tales_api.Calculators.calculators import (
    AvgEHalfCalculator,
    CAResistanceCalculator,
    CVDescriptorCalculator,
    CVDiffusionCalculatorMicro,
    CVChargeTransferCalculatorMicro,
)


@pytest.mark.unit
def test_cv_peaks_for_analysis_rejects_mismatched_counts():
    cal = CVDescriptorCalculator({"scan_data": "scan_data"})
    cal.peaks = lambda data, **kwargs: {"forward": [[0.1, 1.0]], "reverse": []}
    with pytest.raises(ValueError, match="forward peaks"):
        cal.peaks_for_analysis({"scan_data": []})


@pytest.mark.unit
def test_cv_peak_splitting_and_e_half():
    cal = CVDescriptorCalculator({"scan_data": "scan_data"})
    cal.peaks_for_analysis = lambda data, **kwargs: ([0.2, 0.4], [0.1, 0.3])
    assert cal.peak_splittings({}) == [0.1, 0.1]
    assert cal.e_half({}) == [0.15, 0.35]


@pytest.mark.unit
@pytest.mark.parametrize(
    "delta, upper, quasi, expected",
    [(0.050, 63, 200, "reversible"), (0.080, 63, 200, "quasi-reversible"), (0.250, 63, 200, "irreversible")],
)
def test_cv_reversibility_boundaries(delta, upper, quasi, expected):
    cal = CVDescriptorCalculator({"scan_data": "scan_data"})
    cal.peaks_for_analysis = lambda data, **kwargs: ([0.1], [0.1 - delta])
    assert cal.reversibility({}, rev_upperbound=upper, quasi_rev_upperbound=quasi) == [expected]


@pytest.mark.unit
def test_avg_e_half_accepts_explicit_values():
    cal = AvgEHalfCalculator({"e": "e"})
    assert cal.calculate([{"e": "0.1V"}, {"e": "0.3V"}]) == pytest.approx(0.2)


@pytest.mark.unit
def test_microelectrode_diffusion_positive_and_scientific_output():
    data = {"i_ss": "1e-6A", "n": 1, "C": "1e-3mol/cm^3", "r": "0.002cm"}
    cal = CVDiffusionCalculatorMicro({k: k for k in data})
    assert cal.calculate(data) > 0
    assert isinstance(cal.calculate(data, sci_notation=True), str)


@pytest.mark.unit
def test_microelectrode_charge_transfer_positive_or_nan_boundary():
    data = {"e_half": "0.2V", "e_ref": "0.1V", "n": 1, "T": "298K", "D": "1e-5cm^2/s", "r": "0.002cm"}
    cal = CVChargeTransferCalculatorMicro({k: k for k in data})
    assert isinstance(cal.calculate(data), float)


@pytest.mark.unit
def test_ca_resistance_with_error():
    data = {
        "i_s": [1.0, 1.1, 0.9, 1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3],
        "t_s": list(range(10)),
        "pulse_width": "1s",
        "steps": 4,
        "low_e": "1V",
    }
    result = CAResistanceCalculator({k: k for k in data}).calculate(data, return_error=True)
    assert len(result) == 2
    assert result[0] > 0
