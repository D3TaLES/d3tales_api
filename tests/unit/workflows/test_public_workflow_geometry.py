"""Geometry and completion dependencies in the public workflow builders."""
import pytest

from d3tales_api.Workflows.wf_writer import (
    d3tales_wf,
    hf_wf,
    huckaba_wf,
    just_anion,
    just_initialize,
    just_nmr,
    just_tddft,
    molDisc_wf,
    sams_wf,
    solv_wf,
    specialized_job,
)
from d3tales_api.Workflows.utils import forcefield_conformer_for_task, geometry_key_candidates
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


def _fw_by_name(workflow, name):
    return next(fw for fw in workflow.fws if fw.name == name)


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


def test_huckaba_tddft_uses_gas_geometry_when_no_solvent_is_requested(paramset):
    workflow = huckaba_wf(paramset, smiles="O", data_store="local", solvent=None)
    tddft = _fw_by_task(workflow, "tddft_groundState")

    assert only_task(tddft)["geometry"] == "opt_groundState"
    assert _parents_by_name(workflow)[tddft.name] == {"opt_groundState"}


def test_solvent_workflow_requires_a_solvent(paramset):
    with pytest.raises(ValueError, match="requires a solvent"):
        solv_wf(paramset, smiles="O", data_store="local", solvent=None)


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


@pytest.mark.parametrize("wtune", [False, True])
def test_d3tales_workflow_geometry_and_wtuning_contract(paramset, water_geom, wtune):
    workflow = d3tales_wf(
        paramset,
        smiles="O",
        data_store="local",
        mol_opt_geom=water_geom,
        wtune=wtune,
    )
    parents = _parents_by_name(workflow)
    init = _fw_by_name(workflow, "initial_molecule_data")
    assert only_task(init)["mol_opt_geom"] == water_geom
    if wtune:
        opt_mol = _fw_by_task(workflow, "opt_mol")
        tuning = _fw_by_task(workflow, "wtuning")
        assert only_task(opt_mol)["geometry"] == "opt_mol"
        assert only_task(tuning)["geometry"] == "opt_mol"
        assert parents[opt_mol.name] == {init.name}
        assert parents[tuning.name] == {opt_mol.name}
        assert all(parents[name] == {tuning.name} for name in (
            "opt_groundState", "opt_anion1", "opt_anion2", "opt_cation1", "opt_cation2",
        ))
    else:
        assert all(
            only_task(fw).get("name") not in {"opt_mol", "wtuning"}
            for fw in workflow.fws
        )
        assert {
            *parents["opt_groundState"],
            *parents["opt_anion1"],
            *parents["opt_anion2"],
            *parents["opt_cation1"],
            *parents["opt_cation2"],
        } == {init.name}


@pytest.mark.parametrize(
    "task_name, geometry, parent_name, solvent",
    [
        ("energy_a1gs", "opt_anion1", "opt_anion1", None),
        ("energy_gsa1", "opt_groundState", "opt_groundState", None),
        ("energy_c1gs", "opt_cation1", "opt_cation1", None),
        ("energy_gsc1", "opt_groundState", "opt_groundState", None),
        ("solv_energy_gsgs", "opt_groundState", "opt_groundState", "acetonitrile"),
        ("solv_energy_a1a1", "opt_anion1", "opt_anion1", "acetonitrile"),
        ("solv_energy_a2a2", "opt_anion2", "opt_anion2", "acetonitrile"),
        ("solv_energy_c1c1", "opt_cation1", "opt_cation1", "acetonitrile"),
        ("solv_energy_c2c2", "opt_cation2", "opt_cation2", "acetonitrile"),
        ("tddft_groundState", "opt_groundState", "opt_groundState", None),
        ("tddft_anion1", "opt_anion1", "opt_anion1", None),
        ("tddft_anion2", "opt_anion2", "opt_anion2", None),
        ("tddft_cation1", "opt_cation1", "opt_cation1", None),
        ("tddft_cation2", "opt_cation2", "opt_cation2", None),
    ],
)
def test_d3tales_calculations_use_expected_geometry_and_solvent(
    paramset, task_name, geometry, parent_name, solvent
):
    workflow = d3tales_wf(paramset, smiles="O", data_store="local", wtune=False)
    fw = _fw_by_task(workflow, task_name)
    task = only_task(fw)

    assert task["geometry"] == geometry
    assert task.get("solvent") == solvent
    assert _parents_by_name(workflow)[fw.name] == {parent_name}


@pytest.mark.parametrize("wtune", [False, True])
def test_sams_workflow_has_valid_initial_geometry_and_wtuning_contract(paramset, water_geom, wtune):
    workflow = sams_wf(
        paramset,
        smiles="O",
        mol_opt_geom=water_geom,
        data_store="local",
        wtune=wtune,
    )
    parents = _parents_by_name(workflow)
    init = _fw_by_name(workflow, "initial_molecule_data")

    assert only_task(init)["mol_opt_geom"] == water_geom
    if wtune:
        tuning = _fw_by_task(workflow, "wtuning")
        assert tuning in workflow.fws
        assert only_task(tuning)["geometry"] == "opt_mol"
        assert parents[tuning.name] == {init.name}
        assert all(parents[name] == {tuning.name} for name in (
            "opt_groundState", "opt_anion1", "opt_cation1",
        ))
    else:
        assert all(only_task(fw).get("name") != "wtuning" for fw in workflow.fws)
        assert all(parents[name] == {init.name} for name in (
            "opt_groundState", "opt_anion1", "opt_cation1",
        ))


def test_sams_workflow_requires_geometry_for_wtuning(paramset):
    with pytest.raises(ValueError, match="requires mol_opt_geom"):
        sams_wf(paramset, smiles="O", wtune=True)


def test_sams_workflow_can_explicitly_fallback_to_smiles_for_wtuning(paramset):
    workflow = sams_wf(
        paramset,
        smiles="O",
        wtune=True,
        geom_from_smiles=True,
        data_store="local",
    )
    tuning = _fw_by_task(workflow, "wtuning")

    assert only_task(tuning)["geom_from_smiles"] is True
    assert only_task(tuning)["geometry"] == "opt_mol"


def test_sams_workflow_without_tuning_can_start_from_smiles(paramset):
    workflow = sams_wf(paramset, smiles="O", wtune=False, data_store="local")
    opt = _fw_by_task(workflow, "opt_groundState")

    assert only_task(opt)["geometry"] == "opt_mol"
    assert only_task(opt)["geom_from_smiles"] is True
    assert _parents_by_name(workflow)[opt.name] == {"initial_molecule_data"}


def test_forcefield_conformer_for_molopt_ignores_fallback_flag(monkeypatch, water_geom):
    expected = object()
    calls = []

    def generate(smiles):
        calls.append(smiles)
        return expected

    monkeypatch.setattr(
        "d3tales_api.Workflows.utils.generate_forcefield_conformer",
        generate,
    )
    task = {"name": "opt_mol", "geom_from_smiles": False, "smiles": "O"}

    assert forcefield_conformer_for_task(task, {"opt_mol_geom": water_geom}) is expected
    assert calls == ["O"]


def test_forcefield_conformer_requires_smiles_for_opt_mol():
    with pytest.raises(ValueError, match="SMILES string is required"):
        forcefield_conformer_for_task({"name": "opt_mol"}, {})


def test_rdkit_forcefield_conformer_generation_produces_3d_geometry():
    from d3tales_api.Workflows.utils import generate_forcefield_conformer

    molecule = generate_forcefield_conformer("O")

    assert [site.specie.symbol for site in molecule] == ["O", "H", "H"]
    assert all(site.coords.shape == (3,) for site in molecule)
    assert all(0.7 < molecule.get_distance(0, index) < 1.2 for index in (1, 2))


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


@pytest.mark.parametrize("solvent, expected_geometry", [
    (None, "opt_groundState"),
    ("acetonitrile", "solv_opt_groundState"),
])
def test_tddft_geometry_matches_optimization_solvent(paramset, solvent, expected_geometry):
    workflow = just_tddft(paramset, smiles="O", data_store="local", solvent=solvent)
    tddft = _fw_by_task(workflow, "tddft_groundState")

    assert only_task(tddft)["geometry"] == expected_geometry
    assert _parents_by_name(workflow)[tddft.name] == {expected_geometry}


@pytest.mark.parametrize("builder", [
    d3tales_wf,
    just_anion,
    solv_wf,
    hf_wf,
    huckaba_wf,
    just_nmr,
    just_tddft,
    just_initialize,
    specialized_job,
    molDisc_wf,
    sams_wf,
])
def test_public_workflow_geometry_consumers_have_an_ancestor_source(builder, paramset, water_geom):
    paramset.nmr_groundState = paramset.energy_groundState
    workflow = builder(
        paramset,
        smiles="O",
        data_store="local",
        mol_opt_geom=water_geom,
        name_tag="batch_",
    )
    direct_parents = _parents_by_name(workflow)
    fw_by_name = {fw.name: fw for fw in workflow.fws}

    def ancestors(name):
        found = set()
        pending = list(direct_parents[name])
        while pending:
            parent = pending.pop()
            if parent not in found:
                found.add(parent)
                pending.extend(direct_parents[parent])
        return found

    for fw in workflow.fws:
        task = only_task(fw)
        geometry = task.get("geometry")
        if not geometry:
            continue

        candidates = set(geometry_key_candidates(task, {}, geometry))
        ancestor_names = ancestors(fw.name)
        source_keys = set()
        has_initial_geometry = False
        for parent_name in ancestor_names:
            source_task = only_task(fw_by_name[parent_name])
            if type(source_task).__name__ == "MoleculeInit":
                has_initial_geometry = True
            elif type(source_task).__name__ in {
                "RunGaussianOpt", "RunOrcaOpt", "RunPySCFOpt",
            }:
                source_keys.add(source_task.get("name_tag", "") + source_task.get("name", ""))

        assert candidates & source_keys or (has_initial_geometry and "opt_mol" in candidates), (
            builder.__name__,
            fw.name,
            geometry,
            ancestor_names,
            source_keys,
        )
