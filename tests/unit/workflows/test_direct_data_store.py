"""Tests for the direct MongoDB workflow store without contacting live services."""
import sys
import types
import zipfile

import pytest

from d3tales_api.D3database import workflow_adapter
from d3tales_api.Workflows import Initialize, utils
from d3tales_api.Workflows import Gaussian, Orca, pyscf
from d3tales_api.Workflows.ParamSet import GaussianParameters


@pytest.mark.unit
def test_direct_store_reads_molecule_and_computation_documents(monkeypatch):
    molecule = {
        "_id": "M1",
        "mol_info": {
            "smiles": "O",
            "groundState_charge": 0,
            "groundState_spin": 1,
            "init_structure": "1\nH\nH 0 0 0",
        },
        "mol_characterization": {"omega": {"value": 0.18}},
    }

    class Collection:
        def find_one(self, query, projection=None, sort=None):
            if query == {"_id": "M1"}:
                return molecule
            if query == {"_id": "HASH"}:
                return {"_id": "HASH", "data": {"geometry": [{"name": "H", "xyz": [0, 0, 0]}]}}
            return None

    monkeypatch.setattr(
        workflow_adapter,
        "_direct_collection",
        lambda name: Collection() if name == "base" else Collection(),
    )
    assert workflow_adapter.d3t_direct_get_smiles("M1") == "O"
    assert workflow_adapter.d3t_direct_get_ground_state("M1", "charge") == 0
    assert workflow_adapter.d3t_direct_get_omega("M1") == 0.18
    assert workflow_adapter.d3t_direct_get_init_structure("M1")[0].specie.symbol == "H"
    assert workflow_adapter.d3t_direct_get_raw_computation("HASH", "data.geometry") == {
        "geometry": [{"name": "H", "xyz": [0, 0, 0]}]
    }


@pytest.mark.unit
def test_direct_geometry_and_ground_state_helpers_use_direct_adapter(monkeypatch):
    calls = []
    monkeypatch.setattr(utils, "d3t_direct_get_ground_state", lambda identifier, prop: calls.append(("gs", prop)) or 0)
    monkeypatch.setattr(
        utils,
        "d3t_direct_get_raw_computation",
        lambda key, fields=None: calls.append(("geometry", key, fields))
        or {"geometry": [{
            "name": "H",
            "species": [{"element": "H", "occu": 1}],
            "xyz": [0, 0, 0],
            "properties": {},
        }]},
    )
    assert utils.get_groundState("M1", data_store="d3tales_direct") == 0
    geometry = utils.get_db_geom("HASH", data_store="d3tales_direct")
    assert geometry[0].specie.symbol == "H"
    assert calls == [
        ("gs", "charge"),
        ("geometry", "HASH", "data.geometry"),
    ]


@pytest.mark.unit
def test_direct_molecule_init_never_uses_website_adapter(monkeypatch):
    monkeypatch.setattr(Initialize, "d3t_get_smiles", lambda identifier: pytest.fail("unexpected REST lookup"))
    monkeypatch.setattr(Initialize, "d3t_direct_get_smiles", lambda identifier: "O")
    monkeypatch.setattr(Initialize, "get_groundState", lambda *args, **kwargs: 0 if kwargs["prop"] == "charge" else 1)

    update = Initialize.MoleculeInit(
        identifier="M1",
        data_store="d3tales_direct",
    ).run_task({}).update_spec
    assert update["identifier"] == "M1"
    assert update["data_store"] == "d3tales_direct"
    assert update["gs_charge"] == 0
    assert update["gs_spin"] == 1


@pytest.mark.unit
def test_direct_molecule_registration_uses_base_collection_without_legacy_wrapper(monkeypatch):
    inserted = []

    class Collection:
        def find_one(self, query):
            return None

        def count_documents(self, query):
            return 0

        def insert_one(self, document):
            inserted.append(document)

    monkeypatch.setattr(workflow_adapter, "_direct_collection", lambda name: Collection())
    monkeypatch.setattr("d3tales_api.database_info.source_groups", lambda: {"Risko": "06"})

    class MoleculeInfo:
        mol_info_dict = {"smiles": "O", "groundState_charge": 0, "groundState_spin": 1}

        def __init__(self, smiles, names, origin_group):
            assert (smiles, names, origin_group) == ("O", ["water"], "Risko")

    monkeypatch.setattr("d3tales_api.Processors.info_from_smiles.GenerateMolInfo", MoleculeInfo)
    monkeypatch.setattr("nanoid.generate", lambda alphabet, size: "ABCD")
    identifier, charge, spin = workflow_adapter.d3t_direct_register_molecule(
        "O", names=["water"], origin_group="Risko", public=True
    )
    assert (identifier, charge, spin) == ("06ABCD", 0, 1)
    assert inserted == [{
        "_id": "06ABCD",
        "mol_info": MoleculeInfo.mol_info_dict,
        "public": True,
    }]


@pytest.mark.unit
def test_direct_result_storage_inserts_record_and_archives_existing_worker_files(monkeypatch, tmp_path):
    calls = []

    class FakeDatabase:
        def __init__(self, **kwargs):
            calls.append(("init", kwargs))

        def insert(self, identifier):
            calls.append(("insert", identifier))

    monkeypatch.setattr("d3tales_api.D3database.d3database.D3Database", FakeDatabase)
    logfile = tmp_path / "calc.log"
    fchk = tmp_path / "calc.fchk"
    logfile.write_text("log", encoding="utf-8")
    fchk.write_text("checkpoint", encoding="utf-8")
    archive_path = tmp_path / "calc.zip"
    record = {"_id": "HASH", "mol_id": "M1", "data": {"geometry": []}}

    assert workflow_adapter.d3t_direct_store_calculation(
        record, files=[str(logfile), str(fchk), str(tmp_path / "missing")], archive_path=str(archive_path)
    ) == str(archive_path)
    assert calls[0][1]["database"] == "d3tales"
    assert calls[0][1]["collection_name"] == "computation"
    assert calls[-1] == ("insert", "HASH")
    with zipfile.ZipFile(archive_path) as archive:
        assert sorted(archive.namelist()) == ["calc.fchk", "calc.log"]


@pytest.mark.unit
def test_direct_computation_runs_gaussian_back_to_front_after_computation_insert(monkeypatch):
    calls = []

    class FakeDatabase:
        def __init__(self, **kwargs):
            calls.append(("database_init", kwargs["collection_name"]))

        def insert(self, identifier):
            calls.append(("database_insert", identifier))

    monkeypatch.setattr("d3tales_api.D3database.d3database.D3Database", FakeDatabase)
    monkeypatch.setattr(
        workflow_adapter,
        "d3t_direct_update_frontend",
        lambda record: calls.append(("frontend_update", record["_id"])),
    )
    workflow_adapter.d3t_direct_store_calculation({
        "_id": "HASH",
        "mol_id": "M1",
        "calculation_type": "opt_groundState",
        "data": {"conditions": {"code_name": "Gaussian"}},
    })

    assert calls == [
        ("database_init", "computation"),
        ("database_insert", "HASH"),
        ("frontend_update", "HASH"),
    ]


@pytest.mark.unit
def test_direct_frontend_update_runs_gaussian_characterization(monkeypatch):
    calls = []

    class Characterization:
        def __init__(self, identifier, **kwargs):
            calls.append((identifier, kwargs))

    monkeypatch.setattr(
        "d3tales_api.D3database.back2front.Gaus2FrontCharacterization",
        Characterization,
    )
    workflow_adapter.d3t_direct_update_frontend({
        "_id": "HASH",
        "mol_id": "M1",
        "calculation_type": "opt_groundState",
        "data": {"conditions": {"code_name": "Gaussian", "functional": "B3LYP"}},
    })

    assert calls == [("M1", {
        "calculation_type": "opt_groundState",
        "conditions": {"code_name": "Gaussian", "functional": "B3LYP"},
        "rmsd": False,
        "insert": True,
        "verbose": 0,
    })]


@pytest.mark.unit
def test_direct_wtuning_updates_molecule_omega_without_gaussian_characterization(monkeypatch):
    calls = []

    class FakeDatabase:
        def __init__(self, **kwargs):
            calls.append(("init", kwargs))

        def insert(self, identifier, nested=False):
            calls.append(("insert", identifier, nested))

    monkeypatch.setattr("d3tales_api.D3database.d3database.D3Database", FakeDatabase)
    workflow_adapter.d3t_direct_update_frontend({
        "_id": "W-HASH",
        "mol_id": "M1",
        "calculation_type": "wtuning",
        "data": {
            "conditions": {"code_name": "Gaussian", "tuning_parameter": 0.1234},
            "omega": 0.12349,
        },
    })

    init_kwargs = calls[0][1]
    assert init_kwargs["database"] == "d3tales"
    assert init_kwargs["collection_name"] == "base"
    assert init_kwargs["instance"]["mol_characterization"]["omega"] == {
        "source_hash_ids": ["W-HASH"],
        "conditions": {"code_name": "Gaussian"},
        "value": 0.1234,
        "unit": "",
    }
    assert calls[1] == ("insert", "M1", True)


@pytest.mark.unit
def test_engine_post_jobs_route_direct_results_to_database_not_website(monkeypatch, tmp_path):
    stored = []
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        Gaussian,
        "d3t_direct_store_calculation",
        lambda record, **kwargs: stored.append(("gaussian", record, kwargs)),
    )
    monkeypatch.setattr(Gaussian, "d3t_submit_calculation", lambda *args, **kwargs: pytest.fail("website upload"))
    parser = types.ModuleType("d3tales_api.Processors.parser_dft")

    class ProcessGausLog:
        def __init__(self, _id, metadata):
            self.data_dict = {"_id": "G-HASH", "mol_id": _id, "data": {"geometry": []}}

    parser.ProcessGausLog = ProcessGausLog
    monkeypatch.setitem(sys.modules, "d3tales_api.Processors.parser_dft", parser)
    monkeypatch.setattr(Gaussian, "get_hash_id", lambda *args, **kwargs: "G-HASH")
    gaussian = Gaussian.RunGaussianEnergy.__new__(Gaussian.RunGaussianEnergy)
    gaussian.data_store = "d3tales_direct"
    gaussian.file_log = str(tmp_path / "calc.log")
    gaussian.file_fchk = str(tmp_path / "calc.fchk")
    gaussian.file_chk = str(tmp_path / "calc.chk")
    gaussian.file_com = str(tmp_path / "calc.com")
    gaussian.calc_name = gaussian.full_name = "opt_mol"
    gaussian.identifier = "M1"
    gaussian.gs_charge = 0
    gaussian.working_dir = str(tmp_path)
    gaussian.store_dir = str(tmp_path / "results")
    gaussian.runfile_log = None
    (tmp_path / "calc.log").write_text("log", encoding="utf-8")
    (tmp_path / "calc.fchk").write_text("fchk", encoding="utf-8")
    gaussian.post_job(delete_files=[str(tmp_path / "calc.chk")])

    monkeypatch.setattr(
        Orca,
        "d3t_direct_store_calculation",
        lambda record, **kwargs: stored.append(("orca", record, kwargs)),
    )
    monkeypatch.setattr(Orca, "d3t_submit_calculation", lambda *args, **kwargs: pytest.fail("website upload"))
    orca = Orca.RunOrcaEnergy.__new__(Orca.RunOrcaEnergy)
    orca.data_store = "d3tales_direct"
    orca.file_result = str(tmp_path / "orca.result.json")
    orca.file_log = str(tmp_path / "orca.out")
    orca.file_input = str(tmp_path / "orca.inp")
    orca.file_gbw = str(tmp_path / "orca.gbw")
    orca.calc_name = "energy"
    orca.identifier = "M1"
    orca.full_name = "energy"
    orca.working_dir = str(tmp_path)
    orca.store_dir = str(tmp_path / "results")
    orca.local_key = lambda calc_name: "O-HASH"
    orca.runfile_log = str(tmp_path / "orca-runfiles.log")
    (tmp_path / "orca.out").write_text("out", encoding="utf-8")
    orca.post_job({"_id": "O-HASH", "mol_id": "M1", "data": {}})

    monkeypatch.setattr(
        pyscf,
        "save_local_result",
        lambda *args, **kwargs: pytest.fail("direct result should not be stored only locally"),
    )
    monkeypatch.setattr(
        "d3tales_api.D3database.workflow_adapter.d3t_direct_store_calculation",
        lambda record, **kwargs: stored.append(("pyscf", record, kwargs)),
    )
    pyscf_task = pyscf.RunPySCFEnergy.__new__(pyscf.RunPySCFEnergy)
    pyscf_task.data_store = "d3tales_direct"
    pyscf_task.identifier = "M1"
    pyscf_task.calc_name = "energy"
    pyscf_task.store_dir = str(tmp_path)
    pyscf_task.working_dir = str(tmp_path)
    pyscf_task.file_log = str(tmp_path / "pyscf.log")
    pyscf_task.file_chk = str(tmp_path / "pyscf.chk")
    pyscf_task.file_molden = str(tmp_path / "pyscf.molden")
    pyscf_task.file_xyz = str(tmp_path / "pyscf.xyz")
    pyscf_task.runfile_log = str(tmp_path / "pyscf-runfiles.log")
    pyscf_task.gs_charge = 0
    pyscf_task.local_key = lambda calc_name: "P-HASH"
    pyscf_task.local_conditions = lambda: {"code_name": "PySCF"}
    (tmp_path / "pyscf.log").write_text("log", encoding="utf-8")
    pyscf_task.post_job(data={"energy": -1.0})

    assert [entry[0] for entry in stored] == ["gaussian", "orca", "pyscf"]
    assert all(entry[2]["archive_path"].endswith(".zip") for entry in stored)


@pytest.mark.unit
@pytest.mark.parametrize(
    "data_store, result_store",
    [("local", "local"), ("d3tales_web", "local"), ("d3tales_direct", "d3tales_direct")],
)
def test_pyscf_cache_uses_the_store_where_results_are_persisted(data_store, result_store, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    paramset = GaussianParameters(
        route_parameters={}, charge=0, multiplicity=1, basis_set="sto-3g", functional="HF"
    )
    task = pyscf.RunPySCFEnergy(
        name="energy",
        paramset=paramset,
        path=str(tmp_path),
        data_store=data_store,
    )
    task.setup_files({}, calc_type="energy")
    assert task.result_data_store == result_store


@pytest.mark.unit
def test_old_data_store_name_is_rejected():
    with pytest.raises(ValueError, match="data_store"):
        utils.get_data_store({"data_store": "d3tales"})
    with pytest.raises(ValueError, match="data_store"):
        utils.get_groundState("M1", smiles="O", data_store="d3tales")
    with pytest.raises(ValueError, match="data_store"):
        utils.get_existing_data("HASH", data_store="d3tales")
