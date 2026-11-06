import copy
import json

import pytest

from d3tales_api.Calculators.calculators import CVDescriptorCalculator


NUM_ELECTRONS = 1


@pytest.fixture
def cv_connector():
    return {
        "A": "data.conditions.working_electrode_surface_area",
        "v": "data.conditions.scan_rate",
        "C": "data.conditions.redox_mol_concentration",
        "n": "num_electrodes",
        "X": f"data.peak_splittings.{NUM_ELECTRONS - 1}",
        "T": "data.conditions.temperature",
        "D": "diffusion",
        "scan_data": "data.scan_data",
        "variable_prop": "data.conditions.scan_rate.value",
    }


@pytest.fixture
def cv_entry(data_dir):
    """
    Load a single processed CV experiment for regression testing.

    Using a local fixture avoids requiring a live D3TaLES database
    connection during unit testing.
    """
    path = (
        data_dir
        / "electrochemistry"
        / "cv"
        / "extracted_cv_data.json"
    )

    if not path.exists():
        pytest.skip(
            "Missing CV regression fixture: "
            "tests/data/electrochemistry/cv/extracted_cv_data.json"
        )

    with path.open("r", encoding="utf-8") as f:
        entry = json.load(f)

    assert isinstance(entry, dict)
    assert "data" in entry
    assert "scan_data" in entry["data"]

    entry = copy.deepcopy(entry)

    # Required by the calculator connector.
    entry["num_electrodes"] = NUM_ELECTRONS

    return entry


@pytest.fixture
def descriptor_calculator(cv_connector):
    return CVDescriptorCalculator(
        connector=cv_connector
    )


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_peaks(
    descriptor_calculator,
    cv_entry,
):
    result = descriptor_calculator.peaks(cv_entry)

    assert isinstance(result, dict)
    assert "forward" in result
    assert "reverse" in result

    assert result["forward"]
    assert result["reverse"]


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_reversibility(
    descriptor_calculator,
    cv_entry,
):
    result = descriptor_calculator.reversibility(cv_entry)

    assert result == [
        "quasi-reversible",
        "quasi-reversible",
    ]


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_e_half(
    descriptor_calculator,
    cv_entry,
):
    result = descriptor_calculator.e_half(cv_entry)

    assert len(result) == 2

    assert result[0] == pytest.approx(-0.156, abs=0.001)
    assert result[1] == pytest.approx(0.566, abs=0.001)


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_peak_splittings(
    descriptor_calculator,
    cv_entry,
):
    result = descriptor_calculator.peak_splittings(cv_entry)

    assert result == pytest.approx(
        [0.088, 0.088],
        abs=0.001,
    )


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_middle_sweep(
    descriptor_calculator,
    cv_entry,
):
    result = descriptor_calculator.middle_sweep(cv_entry)

    assert len(result) == 2

    # Each branch should contain CV points.
    assert len(result[0]) > 0
    assert len(result[1]) > 0


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_descriptor_methods_are_deterministic(
    descriptor_calculator,
    cv_entry,
):
    assert (
        descriptor_calculator.peaks(cv_entry)
        == descriptor_calculator.peaks(cv_entry)
    )

    assert (
        descriptor_calculator.reversibility(cv_entry)
        == descriptor_calculator.reversibility(cv_entry)
    )

    assert (
        descriptor_calculator.e_half(cv_entry)
        == descriptor_calculator.e_half(cv_entry)
    )

    assert (
        descriptor_calculator.peak_splittings(cv_entry)
        == descriptor_calculator.peak_splittings(cv_entry)
    )

    assert (
        descriptor_calculator.middle_sweep(cv_entry)
        == descriptor_calculator.middle_sweep(cv_entry)
    )


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_descriptor_methods_do_not_mutate_input(
    descriptor_calculator,
    cv_entry,
):
    original = copy.deepcopy(cv_entry)

    descriptor_calculator.peaks(cv_entry)
    descriptor_calculator.reversibility(cv_entry)
    descriptor_calculator.e_half(cv_entry)
    descriptor_calculator.peak_splittings(cv_entry)
    descriptor_calculator.middle_sweep(cv_entry)

    assert cv_entry == original