import pytest
import copy
import json

from d3tales_api.Calculators.calculators import (
    AvgEHalfCalculator,
    CVChargeTransferCalculator,
    D3Calculator,
    DeltaGSolvCalc,
    DirtyElectrodeDetector,
    EnergyDiffCalc,
)



@pytest.fixture
def cv_entry(data_dir):
    path = (
        data_dir
        / "electrochemistry"
        / "cv"
        / "extracted_cv_data.json"
    )

    with path.open("r", encoding="utf-8") as f:
        entry = json.load(f)

    entry = copy.deepcopy(entry)
    entry["num_electrodes"] = 1

    return entry

@pytest.mark.unit
def test_make_connections_reads_nested_values():
    calc = EnergyDiffCalc(connector={
        "energy_final": "results.final",
        "energy_initial": "results.initial",
    })
    connected = calc.make_connections({"results": {"final": "2eV", "initial": "1eV"}})
    assert set(connected) >= {"energy_final", "energy_initial"}


@pytest.mark.unit
def test_all_connections_exposes_connector_keys():
    connector = {"energy_final": "final", "energy_initial": "initial"}
    calc = EnergyDiffCalc(connector=connector)
    assert set(calc.all_connections) >= set(connector)


@pytest.mark.unit
def test_energy_difference_handles_units_and_precision():
    calc = EnergyDiffCalc(connector={"energy_final": "final", "energy_initial": "initial"})
    assert calc.calculate({"final": "2.125eV", "initial": "1eV"}, precision=3) == pytest.approx(1.12, abs=0.01)


@pytest.mark.unit
def test_energy_difference_negative_result():
    calc = EnergyDiffCalc(connector={"energy_final": "final", "energy_initial": "initial"})
    assert calc.calculate({"final": "1eV", "initial": "2eV"}) == pytest.approx(-1.0)


@pytest.mark.unit
def test_delta_g_solv_zero_cycle():
    connector = {
        "init_eng": "init_eng", "init_corr": "init_corr",
        "init_eng_solv": "init_eng_solv", "fin_eng": "fin_eng",
        "fin_corr": "fin_corr", "fin_eng_solv": "fin_eng_solv",
    }
    data = {
        "init_eng": "0eV", "init_corr": "0eV", "init_eng_solv": "0eV",
        "fin_eng": "0eV", "fin_corr": "0eV", "fin_eng_solv": "0eV",
    }
    assert DeltaGSolvCalc(connector=connector).calculate(data) == pytest.approx(0.0)


@pytest.mark.unit
def test_average_e_half_from_cv_scan_data(cv_entry):
    calc = AvgEHalfCalculator({"scan_data": "data.scan_data"})
    result = calc.calculate([cv_entry])
    assert float(result) == pytest.approx(-0.156, abs=0.01)

@pytest.mark.unit
def test_dirty_electrode_detects_small_current_range(cv_entry):
    calc = DirtyElectrodeDetector(connector={"scan_data": "data.scan_data"})
    result = calc.calculate(cv_entry)
    assert result == True


@pytest.mark.unit
def test_charge_transfer_returns_one_result_per_entry():
    connector = {"T":"T", "D":"D", "v":"v", "n":"n", "X":"X"}
    entries = [
        {"T":"298K", "D":"1e-5cm^2/s", "v":"0.1V/s", "n":1, "X":"0.08V"},
        {"T":"298K", "D":"1e-5cm^2/s", "v":"0.2V/s", "n":1, "X":"0.09V"},
    ]
    result = CVChargeTransferCalculator(connector=connector).calculate(entries)
    assert float(result) == pytest.approx(0.017, abs=0.01)
