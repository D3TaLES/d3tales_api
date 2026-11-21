import pytest

from d3tales_api.Workflows.runtime import (
    GAUSSIAN_GPU_CPU,
    G16_CMD,
    LOCAL_OUTPUT_DIR,
    PATH,
    RUNFILE_LOG,
    WorkerRuntime,
    spec_or_task,
)


@pytest.mark.unit
def test_worker_runtime_resolves_worker_values_and_preserves_task_boundaries():
    task = {
        "data_store": "local",
        "path": PATH,
        "runfile_log": RUNFILE_LOG,
        "local_output_dir": LOCAL_OUTPUT_DIR,
        "g16_cmd": G16_CMD,
        "gpu_cpu": GAUSSIAN_GPU_CPU,
        "paramset": {"charge": 0, "multiplicity": 3},
    }
    spec = {
        "data_store": "LOCAL",
        "path": "/spec/path",
        "_fw_env": {
            "runfile_log": "/worker/runfile.log",
            "local_output_dir": "/worker/results",
            "g16_cmd": "/opt/gaussian/g16",
            "gaussian_gpu_cpu": "0=1",
        },
    }

    runtime = WorkerRuntime.from_task(task, spec, software="gaussian")

    assert runtime.data_store == "local"
    assert runtime.path == "/spec/path"
    assert runtime.runfile_log == "/worker/runfile.log"
    assert runtime.local_output_dir == "/worker/results"
    assert runtime.engine_options["g16_cmd"] == "/opt/gaussian/g16"
    assert runtime.engine_options["gpu_cpu"] == "0=1"
    assert task["paramset"] == {"charge": 0, "multiplicity": 3}


@pytest.mark.unit
def test_spec_values_override_task_values_without_losing_false():
    assert spec_or_task({"restricted": True}, {"restricted": False}, "restricted", True) is False
    assert spec_or_task({"nthreads": 4}, {"nthreads": 0}, "nthreads", 2) == 0


@pytest.mark.unit
def test_worker_runtime_requires_path_environment_when_placeholder_is_used():
    with pytest.raises(KeyError, match="path"):
        WorkerRuntime.from_task({"path": PATH}, {}, software="orca")


@pytest.mark.unit
def test_worker_runtime_defaults_optional_engine_settings():
    runtime = WorkerRuntime.from_task(
        {"orca_cmd": ">>orca_cmd<<", "path": "/scratch"},
        {"_fw_env": {}},
        software="orca",
    )
    assert runtime.engine_options["orca_cmd"] == "orca"

    gaussian_runtime = WorkerRuntime.from_task(
        {"g16_cmd": G16_CMD, "gpu_cpu": GAUSSIAN_GPU_CPU},
        {"_fw_env": {"g16_cmd": "/opt/gaussian/g16"}},
        software="gaussian",
    )
    assert gaussian_runtime.engine_options["gpu_cpu"] is None


@pytest.mark.unit
def test_worker_runtime_requires_gaussian_executable_when_placeholder_is_used():
    with pytest.raises(KeyError, match="g16_cmd"):
        WorkerRuntime.from_task(
            {"g16_cmd": G16_CMD}, {}, software="gaussian", require_executables=True
        )


@pytest.mark.unit
@pytest.mark.parametrize("data_store", ["remote", ""])
def test_worker_runtime_rejects_unknown_data_store(data_store):
    with pytest.raises(ValueError, match="data_store"):
        WorkerRuntime.from_task({"data_store": data_store}, {}, software="pyscf")


@pytest.mark.unit
def test_worker_runtime_uses_default_data_store_when_unspecified():
    assert WorkerRuntime.from_task({"data_store": None}, {}).data_store == "d3tales"
