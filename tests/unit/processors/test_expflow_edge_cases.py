import pytest

from d3tales_api.Processors.expflow_parser import ProcessExpFlowObj, get_id


def molecule(redox_uuid="r1"):
    return {
        "_id": "run1",
        "name": "Run",
        "instruments": [{"name": "CHI"}],
        "reagents": [
            {"uuid": redox_uuid, "name": "Ferricyanide", "smiles": "C1=CC=CC=C1", "type": "redox_molecule", "purity": 1.0},
            {"uuid": "s1", "name": "Water", "smiles": "O", "type": "solvent", "purity": 1.0},
        ],
        "apparatus": [{"uuid": "b1", "name": "Beaker", "type": "beaker"}],
        "workflow": [],
    }


@pytest.mark.unit
def test_expflow_id_helper_prefers_id_then_uuid():
    assert get_id(type("A", (), {"_id": "M1"})()) == "M1"
    assert get_id(type("B", (), {"uuid": "U1"})()) == "U1"


@pytest.mark.unit
def test_expflow_properties_from_dict():
    obj = ProcessExpFlowObj(molecule(), "Risko")
    assert obj.object_id == "run1"
    assert obj.instrument_name == "CHI"
    assert obj.redox_mol.name == "Ferricyanide"
    assert obj.solvent.name == "Water"
    assert obj.get_apparatus("beaker") == "Beaker"
    assert obj.get_reagents("solvent") == [{"name": "Water", "purity": 1.0}]
