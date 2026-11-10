"""Composing Fireworks into a multi-step Workflow: dependencies, naming, parameter routing, serialization."""
import pytest
from fireworks import Workflow

from tests.helpers.workflow_helpers import (INIT_NAME, WTUNE_NAME, example_dft_workflow, fws_by_name,
                                            only_task, parents_by_name)

pytestmark = [pytest.mark.unit, pytest.mark.gaussian]

OPTS = ("opt_groundState", "opt_anion1", "opt_cation1")
ENERGIES = ("energy_a1gs", "energy_gsa1", "energy_c1gs", "energy_gsc1")
TDDFTS = ("tddft_groundState", "tddft_anion1", "tddft_cation1")
EXPECTED_PARENTS = {
    INIT_NAME: set(), WTUNE_NAME: {INIT_NAME},
    **{o: {WTUNE_NAME} for o in OPTS},
    "energy_a1gs": {"opt_anion1"}, "energy_gsa1": {"opt_groundState"},
    "energy_c1gs": {"opt_cation1"}, "energy_gsc1": {"opt_groundState"},
    "tddft_groundState": {"opt_groundState"}, "tddft_anion1": {"opt_anion1"}, "tddft_cation1": {"opt_cation1"},
}


@pytest.fixture
def build(paramset, runtime_kwargs, water_geom):
    def _build(**overrides):
        kw = dict(runtime_kwargs, paramset=paramset, smiles="O", mol_opt_geom=water_geom)
        kw.update(overrides)
        return example_dft_workflow(**kw)
    return _build


def test_dependency_graph(build):
    wf = build()
    assert isinstance(wf, Workflow)
    assert parents_by_name(wf) == EXPECTED_PARENTS
    assert [wf.id_fw[i].name for i in wf.root_fw_ids] == [INIT_NAME]
    assert {wf.id_fw[i].name for i in wf.leaf_fw_ids} == {*ENERGIES, *TDDFTS}


def test_optional_step_can_be_removed(build):
    parents = parents_by_name(build(wtune=False))
    assert WTUNE_NAME not in parents
    assert all(parents[o] == {INIT_NAME} for o in OPTS)
    assert all(parents[n] == EXPECTED_PARENTS[n] for n in (*ENERGIES, *TDDFTS))


def test_species_subset(build):
    names = set(fws_by_name(build(species=("groundState", "cation1"))))
    assert names == {INIT_NAME, WTUNE_NAME, "opt_groundState", "opt_cation1", "energy_c1gs", "energy_gsc1",
                     "tddft_groundState", "tddft_cation1"}


@pytest.mark.parametrize("overrides, expected", [
    ({}, "gaus_O"), ({"wf_tag": "batch_"}, "batch_gaus_O"), ({"identifier": "05XICU"}, "gaus_05XICU"),
])
def test_workflow_name(build, overrides, expected):
    assert build(**overrides).name == expected


def test_priority_and_categories(build):
    for name, fw in fws_by_name(build(priority=7)).items():
        assert fw.spec["_priority"] == 7
        assert fw.spec["_category"] == ("processing" if name == INIT_NAME else "gaussian")


def test_software_switch_applies_to_every_calculation(build):
    for name, fw in fws_by_name(build(software="orca")).items():
        if name != INIT_NAME:
            assert fw.spec["_category"] == "orca"
            assert only_task(fw)["software"] == "orca"


@pytest.mark.parametrize("name, block", [
    (WTUNE_NAME, "wtuning"),
    ("opt_groundState", "opt_groundState"), ("opt_anion1", "opt_anion1"), ("opt_cation1", "opt_cation1"),
    ("energy_a1gs", "energy_groundState"), ("energy_gsa1", "energy_anion1"),
    ("energy_c1gs", "energy_groundState"), ("energy_gsc1", "energy_cation1"),
    ("tddft_groundState", "tddft_groundState"), ("tddft_anion1", "tddft_anion1"), ("tddft_cation1", "tddft_cation1"),
])
def test_each_step_gets_matching_parameter_block(build, paramset, name, block):
    assert only_task(fws_by_name(build())[name])["paramset"] is getattr(paramset, block)


def test_runtime_kwargs_reach_every_step(build, runtime_kwargs):
    for name, fw in fws_by_name(build()).items():
        t = only_task(fw)
        for key in ("data_store", "origin_group", "public"):
            assert t[key] == runtime_kwargs[key], (name, key)
        if name != INIT_NAME:
            for key in ("check_if_already_run", "restricted", "gaussian_file_name", "submit", "run_nto"):
                assert t[key] == runtime_kwargs[key], (name, key)


def test_initialization_receives_molecule(build, water_geom):
    t = only_task(fws_by_name(build())[INIT_NAME])
    assert (t["smiles"], t["mol_opt_geom"], t["identifier"]) == ("O", water_geom, None)


def test_dict_round_trip_preserves_graph_and_parameters(build):
    from d3tales_api.Workflows.ParamSet import GaussianParameters
    restored = Workflow.from_dict(build().to_dict())
    assert parents_by_name(restored) == EXPECTED_PARENTS
    for fw in restored.fws:
        if fw.name != INIT_NAME:
            ps = only_task(fw)["paramset"]
            assert isinstance(ps, GaussianParameters) and ps.functional == "LC-wHPBE", fw.name
