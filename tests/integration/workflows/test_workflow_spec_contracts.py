"""
Data hand-off between steps. FireWorks passes data only through update_spec, and a task with
``geometry="X"`` reads ``fw_spec["X_geom"]``; these tests check that producers and consumers agree.
"""
import pytest

from tests.helpers.workflow_helpers import example_dft_workflow, only_task

pytestmark = [pytest.mark.integration, pytest.mark.gaussian]

WATER_INCHIKEY = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"


@pytest.fixture
def wf(paramset, water_geom):
    return example_dft_workflow(paramset, smiles="O", mol_opt_geom=water_geom, data_store="local",
                                check_if_already_run=False)


@pytest.fixture
def init_update_spec(water_geom):
    """Run MoleculeInit for real in local mode (RDKit only; no database, no files)."""
    from d3tales_api.Workflows.Initialize import MoleculeInit
    return MoleculeInit(smiles="O", data_store="local", mol_opt_geom=water_geom).run_task({}).update_spec


def test_initialization_local_mode(init_update_spec, water_geom):
    assert init_update_spec["identifier"] == WATER_INCHIKEY
    assert (init_update_spec["smiles"], init_update_spec["gs_charge"], init_update_spec["gs_spin"]) == ("O", 0, 1)
    assert init_update_spec["data_store"] == "local"


def test_initialization_publishes_geometry_for_first_calculations(wf, init_update_spec, water_geom):
    """Steps directly after initialization must find their input geometry in the init update_spec."""
    for fw in wf.fws:
        if fw.name == "wtuning":
            assert init_update_spec["{}_geom".format(only_task(fw)["geometry"])] == water_geom


def test_initialization_without_smiles_requires_identifier():
    from d3tales_api.Workflows.Initialize import MoleculeInit
    with pytest.raises(ValueError):
        MoleculeInit(data_store="local").run_task({})
    assert MoleculeInit(identifier="ABC", data_store="local").run_task({}).update_spec["identifier"] == "ABC"


def test_downstream_geometries_come_from_a_parent_optimization(wf):
    parents = {fw.fw_id: [] for fw in wf.fws}
    for pid, cids in wf.links.items():
        for cid in cids:
            parents[cid].append(wf.id_fw[pid])
    checked = 0
    for fw in wf.fws:
        t = only_task(fw)
        if type(t).__name__ in ("RunGaussianEnergy", "RunGaussianTDDFT"):
            produced = {"{}_geom".format(only_task(p)["name"]) for p in parents[fw.fw_id]
                        if type(only_task(p)).__name__ == "RunGaussianOpt"}
            assert "{}_geom".format(t["geometry"]) in produced, fw.name
            checked += 1
    assert checked == 7


def test_wtuning_forwards_molecule_inputs(tmp_path, monkeypatch, paramset, water_geom):
    """W-tuning sits between initialization and the optimizations, so it must pass the inputs on.
    Uses the 'already tuned' path (a local w-tuning record) so Gaussian is never called."""
    from d3tales_api.Workflows.Gaussian import RunWtuning
    from d3tales_api.Workflows.utils import WTUNING_KEY, local_record, save_local_result

    monkeypatch.chdir(tmp_path)
    results_dir = tmp_path / WATER_INCHIKEY / "gaussian" / "results"
    save_local_result(str(results_dir), local_record(WTUNING_KEY, WATER_INCHIKEY, "wtuning", {"omega": 0.2}))

    task = RunWtuning(paramset=paramset.wtuning, name="wtuning", path=str(tmp_path), geometry="opt_mol",
                      data_store="local", check_if_already_run=True)
    update = task.run_task({"identifier": WATER_INCHIKEY, "smiles": "O", "gs_charge": 0, "gs_spin": 1,
                            "opt_mol_geom": water_geom}).update_spec

    assert update["iop_str"] == "0200000000"
    assert update["opt_mol_geom"] == water_geom
    assert update["smiles"] == "O"
    assert update["data_store"] == "local"
