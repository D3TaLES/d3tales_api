import json
import pytest
from pathlib import Path
from d3tales_api.Processors.expflow_parser import ProcessExpFlowObj

@pytest.mark.integration
def test_expflow_legacy_fixture(data_dir):
    folder_path=data_dir / "expflow"
    path=next(Path(folder_path).glob("*workflow.json"), None)
    if not path.exists(): pytest.skip("Move *workflow.json to tests/data/expflow/*workflow.json")
    obj=json.loads(path.read_text(encoding="utf-8"))
    result=ProcessExpFlowObj(obj, "Risko")
    assert result is not None
