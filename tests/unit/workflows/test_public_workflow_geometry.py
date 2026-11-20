"""Geometry and completion dependencies in the public workflow builders."""
import pytest

from d3tales_api.Workflows.wf_writer import d3tales_wf, hf_wf, huckaba_wf, just_anion, just_nmr
from tests.helpers.workflow_helpers import only_task

pytestmark = [pytest.mark.unit, pytest.mark.gaussian]


def _parents_by_name(workflow):
    parents = {fw.fw_id: [] for fw in workflow.fws}
    for parent_id, children in workflow.links.items():
        for child_id in children:
            parents[child_id].append(workflow.id_fw[parent_id].name)
    return {fw.name: set(parents[fw.fw_id]) for fw in workflow.fws}


def _fw_by_task(workflow, task_name):
    return next(fw for fw in workflow.fws if only_task(fw).get("name") == task_name)


def test_just_anion_vertical_energy_uses_ground_state_geometry(paramset):
    workflow = just_anion(paramset, smiles="O", data_store="local")
    parents = _parents_by_name(workflow)
    assert "opt_groundState" in parents
    assert parents["energy_gsa1"] == {"opt_groundState"}
    assert only_task(_fw_by_task(workflow, "energy_gsa1"))["geometry"] == "opt_groundState"


def test_huckaba_tddft_uses_its_solvent_optimized_parent(paramset):
    workflow = huckaba_wf(paramset, smiles="O", data_store="local", solvent="acetonitrile")
    tddft = _fw_by_task(workflow, "tddft_groundState")
    task = only_task(tddft)
    assert task["geometry"] == "solv_opt_groundState"
    assert task["solvent"] == "acetonitrile"
    assert _parents_by_name(workflow)[tddft.name] == {"solv_opt_groundState"}


def test_nmr_workflow_optimizes_and_uses_matching_solvent_geometry(paramset):
    paramset.nmr_groundState = paramset.energy_groundState
    workflow = just_nmr(paramset, smiles="O", data_store="local", name_tag="batch_")
    nmr = next(fw for fw in workflow.fws if only_task(fw).get("paramset") is paramset.nmr_groundState)
    task = only_task(nmr)
    assert task["geometry"] == "batch_solv_opt_groundState"
    assert _parents_by_name(workflow)[nmr.name] == {"batch_solv_opt_groundState"}


def test_name_tagged_geometry_key_is_used_by_consumers(paramset):
    workflow = d3tales_wf(paramset, smiles="O", data_store="local", wtune=False, name_tag="batch_")
    task = only_task(_fw_by_task(workflow, "energy_gsa1"))
    from d3tales_api.Workflows.utils import geometry_key_candidates

    assert geometry_key_candidates(task, {}, task["geometry"]) == [
        "batch_opt_groundState",
        "opt_groundState",
    ]


def test_hf_workflow_custom_tags_match_geometry_producers(paramset):
    workflow = hf_wf(paramset, smiles="O", data_store="local", name_tag="batch_")
    for fw in workflow.fws:
        task = only_task(fw)
        if type(task).__name__ != "RunGaussianEnergy":
            continue
        parent_names = _parents_by_name(workflow)[fw.name]
        parent_geometry_keys = {
            only_task(parent)["name_tag"] + only_task(parent)["name"]
            for parent in workflow.fws
            if parent.name in parent_names and type(only_task(parent)).__name__ == "RunGaussianOpt"
        }
        assert task["geometry"] in parent_geometry_keys, fw.name
