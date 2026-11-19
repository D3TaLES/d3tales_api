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
    assert workflow_utils.get_data_store() == "d3tales"
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
def test_nested_set_and_molecule_id_generation(monkeypatch):
    class Collection:
        def __init__(self):
            self.updates = []

        def update_one(self, query, update, upsert=False):
            self.updates.append((query, update, upsert))

        def find_one(self, query, projection=None):
            return None

        def count_documents(self, query):
            return 1 if query == {"_id": "07ABCD"} else 0

    collection = Collection()
    helpers._nested_set(collection, "M1", "mol_info", {"smiles": "CC", "names": ["ethane"]})
    assert collection.updates == [
        ({"_id": "M1"}, {"$set": {"mol_info.smiles": "CC"}}, True),
        ({"_id": "M1"}, {"$set": {"mol_info.names": ["ethane"]}}, True),
    ]

    monkeypatch.setattr("nanoid.generate", lambda alphabet, size: "ABCD" if size == 4 else "ABCDEF")
    assert helpers._generate_mol_id(collection, "CC", "Risko", {"Risko": "07"}) == "07ABCDEF"
    existing = type("Collection", (), {
        "find_one": lambda self, query, projection=None: {"_id": "existing"},
    })()
    assert helpers._generate_mol_id(existing, "CC", "Risko", {}) == "existing"


@pytest.mark.unit
def test_register_molecule_generates_id_and_stores_nested_info(monkeypatch):
    updates = []

    class Collection:
        def find_one(self, query, projection=None):
            return None

        def count_documents(self, query):
            return 0

        def update_one(self, query, update, upsert=False):
            updates.append((query, update, upsert))

    collection = Collection()
    monkeypatch.setattr(helpers, "DBconnector", lambda info: type(
        "Connector", (), {"get_collection": lambda self, name: collection})())
    monkeypatch.setattr("d3tales_api.database_info.db_info", lambda: {"frontend": {}})
    monkeypatch.setattr("d3tales_api.database_info.source_groups", lambda: {"Risko": "07"})
    monkeypatch.setattr("nanoid.generate", lambda alphabet, size: "ABCD")

    class MoleculeInfo:
        def __init__(self, smiles, names, origin_group):
            self.mol_info_dict = {"smiles": smiles, "names": names, "groundState_charge": 0}

    monkeypatch.setattr("d3tales_api.Processors.info_from_smiles.GenerateMolInfo", MoleculeInfo)
    identifier, charge, spin = helpers.d3t_register_molecule("CC", names=["ethane"], public=True)

    assert (identifier, charge, spin) == ("07ABCD", 0, None)
    assert updates[0] == ({"_id": "07ABCD"}, {"$set": {"public": True}}, True)
    assert updates[1:] == [
        ({"_id": "07ABCD"}, {"$set": {"mol_info.smiles": "CC"}}, True),
        ({"_id": "07ABCD"}, {"$set": {"mol_info.names": ["ethane"]}}, True),
        ({"_id": "07ABCD"}, {"$set": {"mol_info.groundState_charge": 0}}, True),
    ]


@pytest.mark.unit
def test_mail_sends_to_each_recipient(monkeypatch):
    sent = []

    class SMTP:
        def __init__(self, server, port, context):
            assert server == "mail.example"
            assert port == 465

        def login(self, sender, password):
            assert (sender, password) == ("sender@example.com", "secret")

        def sendmail(self, sender, recipient, message):
            sent.append((recipient, message))

        def quit(self):
            sent.append(("quit", ""))

    monkeypatch.setattr(helpers.smtplib, "SMTP_SSL", SMTP)
    helpers.Mail(sender="sender@example.com", password="secret", smtp_server="mail.example").send(
        ["one@example.com", "two@example.com"], "subject", "body"
    )
    assert [recipient for recipient, _ in sent] == ["one@example.com", "two@example.com", "quit"]
    with pytest.raises(EnvironmentError, match="credentials"):
        helpers.Mail(sender="", password="")
