import pytest

from d3tales_api.Calculators.plotters import CAPlotter, CVPlotter, DFTSpecPlotter


@pytest.mark.unit
def test_cv_plot_data_milliamps():
    data = {"scan_data": [[[0.0, 0.001], [0.1, 0.002]]]}
    result = CVPlotter({"scan_data": "scan_data"}).plot_data(data, a_to_ma=True)
    assert result["y"] == [1.0, 2.0]


@pytest.mark.unit
def test_ca_plot_data_milliamps_and_contract():
    data = {"t_s": [0, 1], "i_s": [0.001, 0.002]}
    result = CAPlotter({"t_s": "t_s", "i_s": "i_s"}).plot_data(data, a_to_ma=True)
    assert result["x"] == [0, 1]
    assert result["y"] == [1.0, 2.0]


@pytest.mark.unit
def test_ca_plot_rejects_mismatched_arrays():
    data = {"t_s": [0, 1], "i_s": [0.001]}
    with pytest.raises(ValueError, match="different lengths"):
        CAPlotter({"t_s": "t_s", "i_s": "i_s"}).plot_data(data)


@pytest.mark.unit
def test_dft_spectrum_does_not_mutate_transitions():
    transitions = [(2.0, 620.0, 0.1), (3.0, 413.0, 0.2)]
    data = {"transitions": transitions, "sigma": 0.1, "step": 0.01}
    original = list(transitions)
    result = DFTSpecPlotter({"transitions": "transitions", "sigma": "sigma", "step": "step"}).plot_data(data)
    assert data["transitions"] == original
    assert len(result["energies"]) > 1
