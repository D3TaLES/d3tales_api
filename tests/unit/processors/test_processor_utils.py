import pytest

from d3tales_api.Processors.utils import break_list, break_scan_data, convert_list_units


@pytest.mark.unit
@pytest.mark.parametrize(
    "values, expected_breaks",
    [
        ([1, 3, 5, 4, 2, 1], [3]),
        ([1, 2, 3, 3, 2], [4]),
        ([1, 1, 1], []),
    ],
)
def test_break_list_detects_turning_points(values, expected_breaks):
    assert break_list(values, return_brks=True) == expected_breaks


@pytest.mark.unit
def test_break_scan_data_splits_forward_and_reverse_sweeps():
    scan = [[0.0, 1.0], [0.5, 2.0], [1.0, 1.0], [0.5, 0.5], [0.0, 0.2]]
    breaks = break_scan_data(scan)
    assert len(breaks) == 2
