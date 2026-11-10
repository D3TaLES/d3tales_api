"""Shared pytest fixtures for workflow tests (imported by the workflows conftest files)."""
import json

import pytest

from tests.helpers.workflow_helpers import MOLECULE_FILE, PARAM_FILE, RecordingLaunchPad


@pytest.fixture(scope="session")
def molecule_data():
    with open(MOLECULE_FILE, "r") as f:
        return json.load(f)


@pytest.fixture
def water_geom(molecule_data):
    return molecule_data["O"]


@pytest.fixture
def paramset():
    """Fresh GausParamSet per test (tasks mutate charge/multiplicity in place)."""
    from d3tales_api.Workflows.ParamSet import GausParamSet
    return GausParamSet.from_json(str(PARAM_FILE))


@pytest.fixture
def runtime_kwargs():
    """Typical kwargs a launcher forwards to every Firework."""
    return dict(public=False, submit=True, data_store="local", check_if_already_run=False, priority=4,
                gaussian_file_name="gaussian", restricted=True, run_nto=False, origin_group="Risko")


@pytest.fixture
def recording_launchpad():
    return RecordingLaunchPad()


@pytest.fixture
def two_molecule_file(tmp_path, molecule_data):
    path = tmp_path / "two_mols.json"
    path.write_text(json.dumps({"O": molecule_data["O"], "C": molecule_data["O"]}))
    return path
