import pytest

from d3tales_api.Workflows import d3tales_utils as helpers
from d3tales_api.Workflows import utils as workflow_utils


@pytest.mark.unit
@pytest.mark.parametrize(
    "function, response, expected",
    [
        ("d3t_get_smiles", [{"mol_info": {"smiles": "CC"}}], "CC"),
        ("d3t_get_ground_state", [{"mol_info": {"groundState_charge": "-1"}}], -1),
        ("d3t_get_omega", [{"mol_characterization": {"omega": {"value": 0.2}}}], 0.2),
        ("d3t_get_raw_computation", [{"data": {"energy": 3}}], {"energy": 3}),
    ],
)
def test_rest_read_helpers_decode_d3tales_responses(monkeypatch, function, response, expected):
    monkeypatch.setattr(helpers, "_rest_get", lambda endpoint: response)
    if function == "d3t_get_raw_computation":
        result = getattr(helpers, function)("hash-1", fields="data")
    else:
        result = getattr(helpers, function)("M1")
    assert result == expected


@pytest.mark.unit
def test_rest_read_helpers_return_empty_values_for_missing_fields(monkeypatch):
    monkeypatch.setattr(helpers, "_rest_get", lambda endpoint: [])
    assert helpers.d3t_get_smiles("M1") == ""
    assert helpers.d3t_get_ground_state("M1") is None
    assert helpers.d3t_get_omega("M1") is None
    assert helpers.d3t_get_raw_computation("hash") is None


@pytest.mark.unit
def test_get_initial_structure_converts_xyz_and_raises_when_absent(monkeypatch):
    monkeypatch.setattr(helpers, "_mol_info_field", lambda identifier, field: "1\nH\nH 0 0 0")
    structure = helpers.d3t_get_init_structure("M1")
    assert structure[0].specie.symbol == "H"

    monkeypatch.setattr(helpers, "_mol_info_field", lambda identifier, field: None)
    with pytest.raises(LookupError, match="init_structure"):
        helpers.d3t_get_init_structure("M1")


@pytest.mark.unit
def test_local_result_store_selection_and_omega_lookup(tmp_path):
    assert workflow_utils.get_data_store() == "d3tales_web"
    assert workflow_utils.get_local_output_dir(fw_spec={"_fw_env": {"local_output_dir": "/tmp/out"}}) == "/tmp/out"
    assert workflow_utils.get_local_output_dir(
        {"local_output_dir": ">>local_output_dir<<"},
        {"_fw_env": {"local_output_dir": str(tmp_path)}},
    ) == str(tmp_path)

    store = str(tmp_path / "M1")
    workflow_utils.save_local_result(
        store,
        workflow_utils.local_record("wtuning", "M1", "wtuning", {"omega": 0.15}),
        key="wtuning",
    )
    assert workflow_utils.get_omega("M1", data_store="local", store_dir=store) == 0.15
    assert workflow_utils.get_existing_data("missing", data_store="local", store_dir=store) == {}
    assert workflow_utils.get_store_dir("M1") is None
    assert workflow_utils.get_store_dir("M1", calc_dir=str(tmp_path)) == str(tmp_path / "results")


@pytest.mark.unit
def test_local_identifier_validates_smiles_and_is_stable():
    identifier = workflow_utils.local_identifier("CCO")
    assert identifier == workflow_utils.local_identifier("CCO")
    assert "-" in identifier
    with pytest.raises(ValueError, match="required"):
        workflow_utils.local_identifier("")
    with pytest.raises(ValueError, match="Invalid SMILES"):
        workflow_utils.local_identifier("not a smiles")


@pytest.mark.unit
def test_register_molecule_delegates_to_canonical_front_db(monkeypatch):
    calls = {}

    class MoleculeInfo:
        def __init__(self, smiles, names, origin_group):
            self.mol_info_dict = {"smiles": smiles, "names": names, "groundState_charge": 0}

    class FrontDB:
        def __init__(self, **kwargs):
            calls.update(kwargs)
            self.id = "07ABCD"

    monkeypatch.setattr("d3tales_api.D3database.d3database.FrontDB", FrontDB)
    monkeypatch.setattr("d3tales_api.Processors.info_from_smiles.GenerateMolInfo", MoleculeInfo)
    identifier, charge, spin = helpers.d3t_register_molecule("CC", names=["ethane"], public=True)

    assert (identifier, charge, spin) == ("07ABCD", 0, None)
    assert calls == {
        "schema_layer": "mol_info",
        "instance": {"smiles": "CC", "names": ["ethane"], "groundState_charge": 0},
        "smiles": "CC",
        "group": "Risko",
        "public": True,
        "collection_name": "base",
    }


@pytest.mark.unit
def test_front_db_generates_unique_molecule_id(monkeypatch):
    from d3tales_api.D3database.d3database import FrontDB

    class Collection:
        def find_one(self, query, projection=None):
            return None

        def count_documents(self, query):
            return int(query == {"_id": "07ABCD"})

    monkeypatch.setattr("d3tales_api.D3database.d3database.source_groups", lambda: {"Risko": "07"})
    monkeypatch.setattr(
        "d3tales_api.D3database.d3database.generate",
        lambda alphabet, size: "ABCD" if size == 4 else "ABCDEF",
    )
    front_db = FrontDB.__new__(FrontDB)
    front_db.smiles = "CC"
    front_db.group = "Risko"
    front_db.coll = Collection()
    front_db.new_molecule = False

    assert front_db.generate_id() == "07ABCDEF"
