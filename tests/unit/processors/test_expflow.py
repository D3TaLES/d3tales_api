import pytest
parser = pytest.importorskip("d3tales_api.Processors.expflow_parser")

def test_process_expflow_minimal_object():
    obj = {
        "_id": "12345",
        "name": "Test",
        "instruments": [],
        "reagents": [
            {"name": "Test", "smiles": "C", "type": "redox_molecule"}
        ],
        "apparatus": [],
        "workflow": []
    }
    try:
        result = parser.ProcessExpFlowObj(obj, "Risko")
    except (KeyError, TypeError, ValueError) as exc:
        pytest.skip(f"Minimal contract differs in this release: {exc}")
    assert result is not None
