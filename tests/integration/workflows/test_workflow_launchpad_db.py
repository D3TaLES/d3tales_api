"""
Workflow submission against a real, throw-away FireWorks LaunchPad (MongoDB).
Configure with FW_TEST_MONGO_HOST / FW_TEST_MONGO_PORT (default localhost:27017).
"""
import os
import uuid

import pytest

from tests.helpers.workflow_helpers import INIT_NAME, MOLECULE_FILE, PARAM_FILE, submit_workflows

pymongo = pytest.importorskip("pymongo")
pytestmark = [pytest.mark.integration, pytest.mark.database, pytest.mark.slow, pytest.mark.gaussian]

HOST = os.environ.get("FW_TEST_MONGO_HOST", "localhost")
PORT = int(os.environ.get("FW_TEST_MONGO_PORT", "27017"))


@pytest.fixture
def launchpad():
    from fireworks import LaunchPad
    client = pymongo.MongoClient(HOST, PORT, serverSelectionTimeoutMS=1500)
    try:
        client.admin.command("ping")
    except Exception as exc:  # pragma: no cover
        client.close()
        pytest.skip("MongoDB not reachable at {}:{} ({})".format(HOST, PORT, exc))
    db_name = "d3tales_pytest_{}".format(uuid.uuid4().hex[:10])
    try:
        yield LaunchPad(host=HOST, port=PORT, name=db_name)
    finally:
        client.drop_database(db_name)
        client.close()


def test_workflow_is_stored(launchpad):
    ids = submit_workflows(launchpad, MOLECULE_FILE, PARAM_FILE, priority=5)
    wf = launchpad.get_wf_by_fw_id(ids[0])
    assert wf.name == "gaus_O" and len(wf.fws) == 12
    assert {fw.spec["_priority"] for fw in wf.fws} == {5}


def test_only_root_step_is_ready(launchpad):
    fw_id = submit_workflows(launchpad, MOLECULE_FILE, PARAM_FILE)[0]
    assert launchpad.get_fw_by_id(fw_id).name == INIT_NAME
    states = {fw.name: fw.state for fw in launchpad.get_wf_by_fw_id(fw_id).fws}
    assert states.pop(INIT_NAME) == "READY"
    assert set(states.values()) == {"WAITING"}


def test_parameters_survive_database_round_trip(launchpad):
    from d3tales_api.Workflows.ParamSet import GaussianParameters
    fw_id = submit_workflows(launchpad, MOLECULE_FILE, PARAM_FILE)[0]
    stored = {fw.name: fw.tasks[0] for fw in launchpad.get_wf_by_fw_id(fw_id).fws}
    for name, (charge, mult) in {"opt_anion1": (-1, 2), "energy_gsc1": (1, 2), "tddft_cation1": (1, 2)}.items():
        ps = stored[name]["paramset"]
        assert isinstance(ps, GaussianParameters) and (ps.charge, ps.multiplicity) == (charge, mult), name
