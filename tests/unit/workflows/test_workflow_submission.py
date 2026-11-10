"""Submitting workflows for a molecule file to a LaunchPad (in-memory LaunchPad; no MongoDB)."""
import pytest

from d3tales_api.Workflows.ParamSet import GausParamSet
from tests.helpers.workflow_helpers import INIT_NAME, MOLECULE_FILE, PARAM_FILE, submit_workflows

pytestmark = [pytest.mark.unit, pytest.mark.gaussian]


@pytest.fixture
def spy_workflow():
    from fireworks import Firework, ScriptTask, Workflow
    calls = []

    def _wf(**kwargs):
        calls.append(kwargs)
        return Workflow([Firework(ScriptTask.from_str("echo test"))], name="spy")
    _wf.calls = calls
    return _wf


@pytest.mark.parametrize("filename", [None, "does_not_exist.json"])
def test_no_molecules_submits_nothing(recording_launchpad, filename):
    assert submit_workflows(recording_launchpad, filename, PARAM_FILE) == []
    assert recording_launchpad.workflows == []


def test_one_workflow_per_molecule(recording_launchpad):
    ids = submit_workflows(recording_launchpad, MOLECULE_FILE, PARAM_FILE)
    assert ids == [1]
    wf = recording_launchpad.workflows[0]
    assert wf.name == "gaus_O" and len(wf.fws) == 12
    assert wf.fws[0].name == INIT_NAME  # the returned id maps to the root step


def test_multiple_molecules_keep_order(recording_launchpad, two_molecule_file):
    assert submit_workflows(recording_launchpad, two_molecule_file, PARAM_FILE) == [1, 13]
    assert [wf.name for wf in recording_launchpad.workflows] == ["gaus_O", "gaus_C"]


def test_bad_parameter_file_raises(recording_launchpad, tmp_path):
    with pytest.raises((FileNotFoundError, OSError)):
        submit_workflows(recording_launchpad, MOLECULE_FILE, tmp_path / "missing.json")


def test_default_kwargs_passed_to_workflow_function(recording_launchpad, spy_workflow, water_geom):
    submit_workflows(recording_launchpad, MOLECULE_FILE, PARAM_FILE, workflow_function=spy_workflow)
    kw = dict(spy_workflow.calls[0])
    assert isinstance(kw.pop("paramset"), GausParamSet)
    assert kw == dict(public=False, submit=True, data_store="local", check_if_already_run=False, priority=4,
                      gaussian_file_name="gaussian", restricted=True, wtune=True, run_nto=False,
                      origin_group="Risko", wf_tag="", smiles="O", mol_opt_geom=water_geom)


@pytest.mark.parametrize("restricted, file_name", [(True, "gaussian"), (False, "gaussian_ur")])
def test_restricted_selects_output_folder(recording_launchpad, spy_workflow, restricted, file_name):
    submit_workflows(recording_launchpad, MOLECULE_FILE, PARAM_FILE, workflow_function=spy_workflow,
                     restricted=restricted)
    assert spy_workflow.calls[0]["gaussian_file_name"] == file_name


def test_options_flow_into_submitted_fireworks(recording_launchpad):
    submit_workflows(recording_launchpad, MOLECULE_FILE, PARAM_FILE, priority=9, wtune=False, wf_tag="b1_")
    wf = recording_launchpad.workflows[0]
    assert wf.name == "b1_gaus_O" and len(wf.fws) == 11
    assert {fw.spec["_priority"] for fw in wf.fws} == {9}
