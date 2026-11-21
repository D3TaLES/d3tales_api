import pytest

from d3tales_api.Workflows import utils


@pytest.mark.unit
def test_original_conditions_hash_and_task_precedence():
    conditions = utils.orig_conditions(
        "B3LYP",
        "6-31G",
        tuning_parameter="0130000",
        solvent="acetonitrile",
    )
    assert conditions["tuning_parameter"] == pytest.approx(0.13)
    assert conditions["solvent"] == {
        "name": "Acetonitrile",
        "model": "implicit_solvent",
        "dielectric_constant": 35.688,
    }
    assert utils.orig_hash_id("M1", "opt", "B3LYP", "6-31G") == utils.orig_hash_id(
        "M1", "opt", "B3LYP", "6-31G"
    )
    assert utils.spec_or_task({"value": 2}, {"value": 0}, "value", default=3) == 0
    assert utils.spec_or_task({"value": 2}, {}, "value", default=3) == 2
    assert utils.spec_or_task({}, {}, "missing", default=3) == 3


@pytest.mark.unit
def test_configured_spin_is_authoritative_and_checked_for_parity():
    from types import SimpleNamespace

    triplet = SimpleNamespace(charge=0, multiplicity=3)
    with pytest.warns(RuntimeWarning, match="preserving the configured value"):
        charge, multiplicity = utils.resolve_charge_multiplicity({}, triplet, {}, 0, 1)
    assert (charge, multiplicity) == (0, 3)
    assert triplet.charge == 0

    unusual_ground_state = SimpleNamespace(charge=0, multiplicity=3)
    with pytest.warns(RuntimeWarning, match="parity-derived value 1"):
        assert utils.resolve_charge_multiplicity({}, unusual_ground_state, {}, 0, 3) == (0, 3)

    override = SimpleNamespace(charge=0, multiplicity=1)
    assert utils.resolve_charge_multiplicity(
        {}, override, {"charge": 1, "multiplicity": 2}, 0, 1
    ) == (1, 2)
    with pytest.raises(ValueError, match="incompatible with the electron-count parity"):
        utils.resolve_charge_multiplicity({}, override, {"multiplicity": 2}, 0, 1)


@pytest.mark.unit
def test_geometry_key_candidates_prefer_tagged_producer_output():
    assert utils.geometry_key_candidates({"name_tag": "batch_"}, {}, "opt_groundState") == [
        "batch_opt_groundState",
        "opt_groundState",
    ]
    assert utils.geometry_key_candidates({"name_tag": "batch_"}, {}, "batch_opt_groundState") == [
        "batch_opt_groundState"
    ]


@pytest.mark.unit
def test_start_from_smiles_and_ground_state_local_fallback(monkeypatch):
    monkeypatch.setattr(utils, "d3t_get_init_structure", lambda identifier: (_ for _ in ()).throw(LookupError()))
    monkeypatch.setattr(utils, "d3t_get_smiles", lambda identifier: "O")
    monkeypatch.setattr("d3tales_api.Processors.info_from_smiles.find_lowest_e_conf",
                        lambda smiles: "1\nwater\nO 0 0 0")
    structure = utils.start_from_smiles("M1", None, data_store="d3tales_web")
    assert len(structure) == 1
    assert utils.get_groundState("M1", smiles="CC", data_store="local", prop="charge") == 0
    assert utils.get_groundState("M1", smiles="[CH3]", data_store="local", prop="spin") == 2

    monkeypatch.setattr(utils, "d3t_get_ground_state", lambda identifier, prop: -1)
    assert utils.get_groundState("M1", data_store="d3tales_web") == -1
    with pytest.raises(ValueError, match="SMILES must be specified"):
        utils.get_groundState("M1", data_store="local")
    with pytest.raises(ValueError, match="prop must"):
        utils.get_groundState("M1", smiles="CC", prop="other", data_store="local")
    with pytest.raises(ValueError, match="Invalid SMILES"):
        utils.get_groundState("M1", smiles="not smiles", data_store="local")


@pytest.mark.unit
def test_start_from_smiles_uses_database_structure_or_rejects_missing_smiles(monkeypatch):
    structure = utils.Molecule.from_str("1\nwater\nO 0 0 0", "xyz")
    monkeypatch.setattr(utils, "d3t_get_init_structure", lambda identifier: structure)
    assert utils.start_from_smiles("M1", None, data_store="d3tales_web") is structure
    monkeypatch.setattr(
        utils,
        "d3t_get_init_structure",
        lambda identifier: (_ for _ in ()).throw(LookupError()),
    )
    monkeypatch.setattr(utils, "d3t_get_smiles", lambda identifier: "")
    with pytest.raises(ValueError, match="no SMILES provided"):
        utils.start_from_smiles("M1", None, data_store="d3tales_web")


@pytest.mark.unit
def test_get_db_geometry_reads_local_records_and_raises_for_missing(tmp_path):
    store = str(tmp_path / "M1")
    sites = [{"name": "H", "species": [{"element": "H", "occu": 1}], "xyz": [0, 0, 0], "properties": {}}]
    utils.save_local_result(
        store,
        utils.local_record("geometry-hash", "M1", "opt", {"geometry": sites}),
    )
    structure = utils.get_db_geom("geometry-hash", data_store="local", store_dir=store)
    assert len(structure) == 1
    assert utils.get_db_geom("missing", data_store="local", store_dir=store) is None
    with pytest.raises(LookupError, match="not found"):
        utils.get_db_geom("missing", error_raise=True, data_store="local", store_dir=store)
    utils.save_local_result(store, utils.local_record("empty", "M1", "opt", {"value": 1}))
    with pytest.raises(LookupError, match="starting gemetry"):
        utils.get_db_geom("empty", error_raise=True, data_store="local", store_dir=store)


@pytest.mark.unit
def test_zip_files_and_tddft_extraction_and_summary(tmp_path, monkeypatch):
    log = tmp_path / "excited.log"
    log.write_text(
        " Excitation energies and oscillator strengths:\n"
        " Excited State 1: Singlet-A 4.00 eV 310.0 nm f=0.25\n"
        " Excited State 2: Triplet-A 5.00 eV 248.0 nm f=0.05\n",
        encoding="utf-8",
    )
    excitations = utils.get_tddft_excitations(str(log))
    assert excitations["Singlet"] == [(4.0, 310.0, 0.25)]
    monkeypatch.setattr(utils, "runtime_from_log", lambda path: 0.5)
    result = utils.parse_tddft(str(log))
    assert result["descriptors"] == {"vs0s1": 4.0, "vs0t1": 5.0}
    assert result["timings"] == {"tddft": 0.5}

    input_file = tmp_path / "input.txt"
    input_file.write_text("data", encoding="utf-8")
    archive_path = utils.zip_files([str(input_file)], str(tmp_path / "input.zip"))
    assert utils.zipfile.is_zipfile(archive_path)


@pytest.mark.unit
def test_write_nto_configures_checkpoint_and_transition(monkeypatch, tmp_path):
    observed = {}

    class Input:
        route_parameters = {"iop(3107": "1300000,1"}
        link0_parameters = {}

        def write_file(self, path, cart_coords):
            observed["path"] = path
            observed["cart_coords"] = cart_coords
            observed["route"] = self.route_parameters
            observed["link0"] = self.link0_parameters

    class Output:
        def __init__(self, path):
            pass

        def to_input(self):
            return Input()

    monkeypatch.setattr(utils, "GaussianOutput", Output)
    logfile = str(tmp_path / "calc.log")
    com_path, chk_path = utils.write_nto(logfile, 2)
    assert com_path.endswith("_state02.com")
    assert chk_path.endswith("_state02.chk")
    assert observed["route"]["Density"] == "(Check,Transition=2)"
    assert observed["route"]["iop(3/107=1300000, 3/108=1300000)"] == ""
    assert observed["link0"]["%oldchk"].endswith("calc.chk")


@pytest.mark.unit
def test_electronic_structure_task_mixin_helpers(tmp_path, monkeypatch):
    class Task(utils.ElectronicStructureTaskMixin):
        identifier = "M1"
        calc_name = "opt_groundState"
        data_store = "local"
        result_data_store = None
        store_dir = str(tmp_path)
        functional = "B3LYP"
        basis = "6-31G"
        gs_charge = 0
        gs_spin = 1
        smiles = "CC"
        working_dir = "/work/M1"
        runfile_log = str(tmp_path / "run.log")

    task = Task()
    task.append_runfile_log("hash-1")
    assert (tmp_path / "run.log").read_text(encoding="utf-8") == "hash-1, /work/M1\n"
    update = task.common_update_spec(extra="value")
    assert update == {
        "identifier": "M1",
        "gs_charge": 0,
        "gs_spin": 1,
        "iop_str": None,
        "data_store": "local",
        "smiles": "CC",
        "extra": "value",
    }
    monkeypatch.setattr(utils, "get_existing_data", lambda key, **kwargs: {"hash": key})
    assert task.existing_data("known") == {"hash": "known"}
    assert task.local_key().startswith("")


@pytest.mark.unit
def test_tddft_parser_handles_triplet_fallback(tmp_path):
    logfile = tmp_path / "spin.log"
    logfile.write_text(
        " Excitation energies and oscillator strengths:\n"
        " Excited State 1: 1.00 2.00 3.00 4.00\n"
        " Excited State 2: 3.00 2.00 1.00 0.00\n",
        encoding="utf-8",
    )
    result = utils.get_tddft_excitations(str(logfile))
    assert result["Singlet"]
    assert result["Triplet"]
