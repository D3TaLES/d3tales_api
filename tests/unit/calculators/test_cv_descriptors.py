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
def cv_entries(data_dir):
    """
    Load the CV entries that were used by the legacy test for mol_id 06QMKG.

    The fixture replaces:

        Expflow().make_query({"mol_id": "06QMKG"}, output="json")

    so the unit test does not require a live database connection.
    """
    path = data_dir / "electrochemistry" / "cv" / "extracted_cv_data.json"

    if not path.exists():
        pytest.skip(
            "Missing CV regression fixture: "
            "tests/data/electrochemistry/cv/extracted_cv_data.json"
        )

    with path.open("r", encoding="utf-8") as f:
        entries = json.load(f)

    entries = copy.deepcopy(entries)

    for entry in entries:
        entry["num_electrodes"] = NUM_ELECTRONS

    return entries


@pytest.fixture
def descriptor_calculator(cv_connector):
    return CVDescriptorCalculator(connector=cv_connector)


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_peaks(descriptor_calculator, cv_entries):
    result = descriptor_calculator.peaks(cv_entries[0])

    assert result == {
        "forward": [[0.091, 4.796e-05]],
        "reverse": [[-0.05, -6.478e-05]],
    }


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_reversibility(descriptor_calculator, cv_entries):
    result = descriptor_calculator.reversibility(cv_entries[0])

    assert result == ["quasi-reversible"]


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_e_half(descriptor_calculator, cv_entries):
    result = descriptor_calculator.e_half(cv_entries[0])

    assert result[0] == [0.02]


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_peak_splittings(descriptor_calculator, cv_entries):
    result = descriptor_calculator.peak_splittings(cv_entries[0])

    assert result == [0.141]


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_middle_sweep(descriptor_calculator, cv_entries):
    result = descriptor_calculator.middle_sweep(cv_entries[0])

    assert len(result) == 2


@pytest.mark.unit
@pytest.mark.electrochemistry
def test_descriptor_methods_are_deterministic(
    descriptor_calculator,
    cv_entries,
):
    entry = cv_entries[0]

    assert descriptor_calculator.peaks(entry) == descriptor_calculator.peaks(entry)
    assert (
        descriptor_calculator.reversibility(entry)
        == descriptor_calculator.reversibility(entry)
    )
    assert (
        descriptor_calculator.e_half(entry)
        == descriptor_calculator.e_half(entry)
    )
    assert (
        descriptor_calculator.peak_splittings(entry)
        == descriptor_calculator.peak_splittings(entry)
    )