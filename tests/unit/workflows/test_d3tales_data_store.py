"""D3TaLES-backed workflow state and data lookups, isolated from live services."""
from pathlib import Path

import pytest
from fireworks import Workflow

from d3tales_api.Workflows import Initialize, utils
from d3tales_api.Workflows.wf_writer import d3tales_wf
from tests.helpers.workflow_helpers import only_task

LOCAL_DB_INFO = Path(__file__).resolve().parents[3] / "db_info.json"
pytestmark = [
    pytest.mark.unit,
    pytest.mark.gaussian,
    pytest.mark.skipif(
        not LOCAL_DB_INFO.is_file(),
        reason="D3TaLES workflow tests require a local db_info.json configuration",
    ),
]


def test_d3tales_store_is_forwarded_to_all_workflow_tasks(paramset):
    workflow = d3tales_wf(
        paramset,
        identifier="TEST123",
        data_store="d3tales",
        wtune=False,
    )
    restored = Workflow.from_dict(workflow.to_dict())

    assert all(only_task(fw)["data_store"] == "d3tales" for fw in restored.fws)


def test_initialize_existing_d3tales_molecule_populates_downstream_spec(
    monkeypatch, water_geom
):
    looked_up = []
    monkeypatch.setattr(Initialize, "d3t_get_smiles", lambda identifier: "O")

    def ground_state(identifier, smiles, prop, data_store):
        looked_up.append((identifier, smiles, prop, data_store))
        return {"charge": 0, "spin": 1}[prop]

    monkeypatch.setattr(Initialize, "get_groundState", ground_state)
    update = Initialize.MoleculeInit(
        identifier="TEST123",
        data_store="d3tales",
        mol_opt_geom=water_geom,
    ).run_task({}).update_spec

    assert update == {
        "identifier": "TEST123",
        "smiles": "O",
        "gs_charge": 0,
        "gs_spin": 1,
        "data_store": "d3tales",
        "mol_opt_geom": water_geom,
        "opt_mol_geom": water_geom,
    }
    assert looked_up == [
        ("TEST123", "O", "charge", "d3tales"),
        ("TEST123", "O", "spin", "d3tales"),
    ]


def test_d3tales_ground_state_uses_remote_value_then_smiles_fallback(monkeypatch):
    calls = []

    def remote_ground_state(identifier, prop):
        calls.append((identifier, prop))
        return -1 if prop == "charge" else None

    monkeypatch.setattr(utils, "d3t_get_ground_state", remote_ground_state)

    assert utils.get_groundState("TEST123", smiles="O", data_store="d3tales") == -1
    assert utils.get_groundState(
        "TEST123", smiles="[CH3]", prop="spin", data_store="d3tales"
    ) == 2
    assert calls == [("TEST123", "charge"), ("TEST123", "spin")]


def test_d3tales_starting_and_optimized_geometries_use_remote_data(monkeypatch):
    initial = utils.Molecule.from_str("1\nH\nH 0 0 0", "xyz")
    optimized_sites = [
        {
            "name": "H",
            "species": [{"element": "H", "occu": 1}],
            "xyz": [0.0, 0.0, 0.0],
            "properties": {},
        }
    ]
    calls = []
    monkeypatch.setattr(
        utils,
        "d3t_get_init_structure",
        lambda identifier: calls.append(("initial", identifier)) or initial,
    )
    monkeypatch.setattr(
        utils,
        "d3t_get_raw_computation",
        lambda geometry_hash, fields: calls.append((geometry_hash, fields))
        or {"geometry": optimized_sites},
    )

    assert utils.start_from_smiles("TEST123", None, data_store="d3tales") is initial
    optimized = utils.get_db_geom("HASH123", data_store="d3tales")
    assert len(optimized) == 1
    assert optimized[0].specie.symbol == "H"
    assert calls == [
        ("initial", "TEST123"),
        ("HASH123", "data.geometry"),
    ]


def test_initialize_d3tales_registration_receives_store_and_molecule_options(monkeypatch):
    calls = {}
    monkeypatch.setattr(Initialize, "d3t_get_smiles", lambda identifier: "")

    def register(smiles, names, origin_group, public):
        calls.update(
            smiles=smiles,
            names=names,
            origin_group=origin_group,
            public=public,
        )
        return "NEW123", 0, 1

    monkeypatch.setattr(Initialize, "d3t_register_molecule", register)
    update = Initialize.MoleculeInit(
        smiles="O",
        mol_name="water",
        data_store="d3tales",
        origin_group="Risko",
        public=True,
    ).run_task({}).update_spec

    assert calls == {
        "smiles": "O",
        "names": ["water"],
        "origin_group": "Risko",
        "public": True,
    }
    assert update == {
        "identifier": "NEW123",
        "smiles": "O",
        "gs_charge": 0,
        "gs_spin": 1,
        "data_store": "d3tales",
        "mol_opt_geom": None,
    }


def test_initialize_rejects_unknown_d3tales_identifier_without_smiles(monkeypatch):
    monkeypatch.setattr(Initialize, "d3t_get_smiles", lambda identifier: "")
    with pytest.raises(ValueError, match="SMILES"):
        Initialize.MoleculeInit(
            identifier="MISSING",
            data_store="d3tales",
        ).run_task({})
