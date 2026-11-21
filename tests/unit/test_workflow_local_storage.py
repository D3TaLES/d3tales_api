import json
import pytest

from d3tales_api.Workflows import utils


@pytest.mark.unit
def test_local_record_round_trip(tmp_path):
    store = tmp_path / "M1"
    record = utils.local_record("abc", "M1", "opt_groundState", {"value": 1})
    path = utils.save_local_result(str(store), record)
    assert utils.read_local_result(str(store), "abc") == {"value": 1}
    with open(path, "r", encoding="utf-8") as handle:
        assert json.load(handle)["mol_id"] == "M1"


@pytest.mark.unit
def test_get_store_dir_is_absolute(tmp_path):
    path = utils.get_store_dir("M1", local_output_dir=str(tmp_path))
    assert path == str((tmp_path / "M1").resolve())


@pytest.mark.unit
@pytest.mark.parametrize("task_store, fw_store, expected", [
    ("d3tales_web", "local", "local"),
    (None, "d3tales_web", "d3tales_web"),
])
def test_get_data_store_precedence(task_store, fw_store, expected):
    task = {"data_store": task_store} if task_store else None
    assert utils.get_data_store(task, {"data_store": fw_store}) == expected


@pytest.mark.unit
def test_get_data_store_rejects_unknown():
    with pytest.raises(ValueError, match="data_store"):
        utils.get_data_store({"data_store": "remote"})
