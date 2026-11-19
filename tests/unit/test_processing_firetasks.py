import zipfile

import pytest

from d3tales_api.Processors import processing_firetasks as tasks


@pytest.mark.unit
def test_process_file_cv_zip_dispatches_parser_and_returns_spec(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    archive = tmp_path / "input.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("scan.txt", "fixture")
        zipped.writestr("readme.md", "metadata")

    captured = {}

    class ParsedCV:
        def __init__(self, path, **kwargs):
            captured.update(path=path, kwargs=kwargs)
            self.data_dict = {"_id": "hash-1", "data": {"conditions": {}}}

    monkeypatch.setattr(tasks, "ProcessChiCV", ParsedCV)
    task = tasks.ProcessFile({
        "filepath": str(archive),
        "mol_id": "M1",
        "submission_data": {"data_type": "cv", "data_category": "experimentation"},
        "metadata": {"source": "test"},
    })

    action = task.run_task({})

    assert action.update_spec["hash_id"] == "hash-1"
    assert action.update_spec["mol_id"] == "M1"
    assert action.update_spec["submission_data"]["all_files_in_zip"] == ["scan.txt", "readme.md"]
    assert captured["kwargs"] == {
        "_id": "M1",
        "submission_info": action.update_spec["submission_data"],
        "metadata": {"source": "test"},
    }
    assert not (tmp_path / "scan.txt").exists()


@pytest.mark.unit
def test_process_file_gaussian_automatic_approval(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    archive = tmp_path / "gaussian.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("output.log", "log")

    parsed = {"_id": "gaussian-hash", "data": {}}
    approved = []

    class ParsedGaussian:
        def __init__(self, **kwargs):
            assert kwargs["metadata"]["mol_file"].endswith("output.log")
            self.data_dict = parsed

    class Characterization:
        @classmethod
        def from_data(cls, data, **kwargs):
            approved.append((data, kwargs))

    monkeypatch.setattr(tasks, "ProcessGausLog", ParsedGaussian)
    monkeypatch.setattr(tasks, "Gaus2FrontCharacterization", Characterization)
    task = tasks.ProcessFile({
        "filepath": str(archive),
        "mol_id": "M2",
        "submission_data": {"data_type": "gaussian"},
        "metadata": {"calculation_type": "opt_groundState"},
        "automatic_approval": True,
    })

    action = task.run_task({})

    assert action.update_spec["processed_data"] == parsed
    assert action.update_spec["automatic_approval"] is True
    assert approved == [(parsed, {"rmsd": tasks.RMSD_DEFAULT, "insert": True})]
    assert not (tmp_path / "output.log").exists()


@pytest.mark.unit
def test_file_transfer_can_skip_external_transfer():
    task = tasks.FileTransferTask({
        "files": ["/remote/file"],
        "dest": "/local/file",
        "remote_server": {"hostname": "example"},
        "skip_transfer": True,
    })
    assert task.run_task({"_fw_env": {"processing_vm_key": "/tmp/key"}}) is None


@pytest.mark.unit
def test_file_transfer_downloads_files_and_closes_sftp(monkeypatch):
    events = []

    class SFTP:
        def get(self, source, destination):
            events.append(("get", source, destination))

        def remove(self, path):
            events.append(("remove", path))

        def close(self):
            events.append(("sftp-close",))

    class SSH:
        def load_system_host_keys(self):
            events.append(("load-keys",))

        def set_missing_host_key_policy(self, policy):
            events.append(("policy",))

        def connect(self, **kwargs):
            events.append(("connect", kwargs))

        def open_sftp(self):
            return SFTP()

        def close(self):
            events.append(("ssh-close",))

    monkeypatch.setattr(tasks.paramiko, "SSHClient", SSH)
    task = tasks.FileTransferTask({
        "files": ["/remote/file"],
        "dest": "/local/file",
        "remote_server": {"hostname": "example"},
    })
    task.run_task({"_fw_env": {"processing_vm_key": "/tmp/key"}})
    assert ("get", "/remote/file", "/local/file") in events
    assert ("remove", "/remote/file") in events
    assert events[-2:] == [("sftp-close",), ("ssh-close",)]


@pytest.mark.unit
def test_send_to_storage_copies_locally_and_updates_submission(tmp_path, monkeypatch):
    source = tmp_path / "original.txt"
    source.write_text("data", encoding="utf-8")
    backend_collection = type("Collection", (), {"update_one": lambda *args, **kwargs: None})()
    monkeypatch.setattr(tasks, "LOCAL_STORAGE", True)
    monkeypatch.setattr(tasks, "TESTING", True)
    monkeypatch.setattr(tasks, "DBconnector", lambda info: type(
        "Connector", (), {"get_collection": lambda self, name: backend_collection})())

    task = tasks.SendToStorage({"filepath": str(source), "hash_id": "hash-1"})
    spec = {
        "mol_id": "M1",
        "_fw_env": {"storage_base_loc": str(tmp_path / "storage")},
        "submission_data": {"data_category": "computation", "data_type": "gaussian"},
        "processed_data": {"submission_info": {}},
    }

    action = task.run_task(spec)

    destination = tmp_path / "storage" / "d3tales" / "M1" / "computation" / "gaussian" / "original.txt"
    assert destination.read_text(encoding="utf-8") == "data"
    assert action.update_spec["processed_data"]["submission_info"]["stored_location"] == str(destination).replace("\\", "/")


@pytest.mark.unit
def test_send_to_storage_transfers_remotely(monkeypatch, tmp_path):
    calls = []

    class SFTP:
        def put(self, source, destination):
            calls.append(("put", source, destination))

        def close(self):
            calls.append(("sftp-close",))

    class SSH:
        def load_system_host_keys(self):
            pass

        def set_missing_host_key_policy(self, policy):
            pass

        def connect(self, **kwargs):
            calls.append(("connect", kwargs))

        def open_sftp(self):
            return SFTP()

        def close(self):
            calls.append(("ssh-close",))

    collection = type("Collection", (), {"update_one": lambda *args, **kwargs: None})()
    monkeypatch.setattr(tasks, "LOCAL_STORAGE", False)
    monkeypatch.setattr(tasks, "TESTING", True)
    monkeypatch.setattr(tasks.paramiko, "SSHClient", SSH)
    monkeypatch.setattr(tasks, "mkdir_p", lambda sftp, path: calls.append(("mkdir", path)))
    monkeypatch.setattr(tasks, "DBconnector", lambda info: type(
        "Connector", (), {"get_collection": lambda self, name: collection})())
    source = str(tmp_path / "input.txt")
    task = tasks.SendToStorage({
        "filepath": source,
        "remote_server": {"hostname": "storage.example"},
    })
    action = task.run_task({
        "mol_id": "M1",
        "_fw_env": {"storage_base_loc": "/data", "storage_key": "/key"},
        "submission_data": {"data_category": "computation", "data_type": "gaussian"},
        "processed_data": {"submission_info": {}},
    })

    assert ("connect", {"hostname": "storage.example", "key_filename": "/key"}) in calls
    assert ("mkdir", "/data/d3tales/M1/computation/gaussian") in calls
    assert any(call[0] == "put" and call[1] == source for call in calls)
    assert action.update_spec["processed_data"]["submission_info"]["stored_location"].endswith("input.txt")


@pytest.mark.unit
def test_update_backend_and_approval_status_dispatch(monkeypatch):
    calls = []
    monkeypatch.setattr(tasks, "BackDB", lambda **kwargs: calls.append(("back", kwargs)))

    tasks.UpdateBackendDB({}).run_task({
        "processed_data": {"_id": "hash"},
        "submission_data": {"data_category": "computation"},
    })
    assert calls == [("back", {"collection_name": "computation", "instance": {"_id": "hash"}})]

    collection = type("Collection", (), {"update_one": lambda self, *args, **kwargs: calls.append((args, kwargs))})()
    monkeypatch.setattr(tasks, "DBconnector", lambda info: type(
        "Connector", (), {"get_collection": lambda self, name: collection})())
    task = tasks.UpdateApprovalStatus({
        "hash_id": "hash",
        "approved_status": True,
        "data_category": "computation",
    })
    task.run_task({})
    assert calls[-1] == (
        ({"_id": "hash"}, {"$set": {"submission_info.approved": True}}),
        {"upsert": True},
    )
    assert tasks.UpdateApprovalStatus({"already_updated": True}).run_task({}) is None


@pytest.mark.unit
def test_update_user_db_handles_already_updated_and_stores_approval(monkeypatch):
    calls = []

    class Collection:
        def find_one(self, query):
            return {"_id": "submission-id"}

        def update_one(self, query, update, upsert=False):
            calls.append(update)

    collection = Collection()
    monkeypatch.setattr(tasks, "DBconnector", lambda info: type(
        "Connector", (), {"get_collection": lambda self, name: collection})())
    task = tasks.UpdateUserDB({})
    task.run_task({
        "hash_id": "processed-hash",
        "automatic_approval": True,
        "processed_data": {"data": {"value": 2}},
    })
    assert calls == [
        {"$set": {"processing_id": "processed-hash"}},
        {"$set": {"processed_data": {"value": 2}}},
        {"$set": {"approved": True}},
    ]

    collection.find_one = lambda query: None
    action = task.run_task({"hash_id": "missing", "processed_data": {"data": {}}})
    assert action.update_spec == {"already_updated": True}


@pytest.mark.unit
def test_update_frontend_db_only_characterizes_gaussian_data(monkeypatch):
    data = {
        "_id": "hash",
        "mol_id": "M1",
        "calculation_type": "opt_groundState",
        "submission_info": {"data_type": "gaussian"},
        "data": {"conditions": {"functional": "B3LYP"}},
    }
    collection = type("Collection", (), {"find_one": lambda self, query: data})()
    monkeypatch.setattr(tasks, "DBconnector", lambda info: type(
        "Connector", (), {"get_collection": lambda self, name: collection})())
    calls = []
    monkeypatch.setattr(tasks, "Gaus2FrontCharacterization", lambda *args, **kwargs: calls.append((args, kwargs)))

    tasks.UpdateFrontendDB({"hash_id": "hash", "data_category": "computation"}).run_task({})
    assert calls == [
        (("M1",), {
            "calculation_type": "opt_groundState",
            "conditions": {"functional": "B3LYP"},
            "rmsd": tasks.RMSD_DEFAULT,
            "insert": True,
        })
    ]

    data["submission_info"]["data_type"] = "cv"
    tasks.UpdateFrontendDB({"hash_id": "hash", "data_category": "experimentation"}).run_task({})
    assert len(calls) == 1


@pytest.mark.unit
def test_approve_jobs_deduplicates_and_calls_rest_api(monkeypatch):
    calls = []

    class Collection:
        def find(self, query):
            return [{"processing_id": "hash"}]

        def count_documents(self, query):
            return 2 if not calls else 1

        def find_one_and_delete(self, query, sort):
            calls.append(("delete", query, sort))

    collection = Collection()
    monkeypatch.setattr(tasks, "D3Database", lambda **kwargs: type(
        "Database", (), {"coll": collection})())
    monkeypatch.setattr(tasks, "RESTAPI", lambda **kwargs: calls.append(("approve", kwargs)))

    tasks.ApproveJobs({"data_category": "computation"}).run_task({})

    assert calls[0][0] == "delete"
    assert calls[1][0] == "approve"
    assert calls[1][1]["endpoint"] == "tools/upload/hash"
