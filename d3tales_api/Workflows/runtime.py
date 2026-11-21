"""Resolve FireWorks task settings separately from worker environment settings."""

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional

DATA_STORES = ("d3tales", "local")
DEFAULT_DATA_STORE = "d3tales"
LOCAL_RESULTS_DIR = "results"
WTUNING_KEY = "wtuning"

ZINDOLIB = ">>ZINDOLIB<<"
ZINDOBIN = ">>ZINDOBIN<<"
ZINDOCTBIN = ">>ZINDOCTBIN<<"
DB_FILE = ">>db_file<<"
AFLOW = ">>aflow<<"
VASP_CMD = ">>vasp_cmd<<"
GAMMA_VASP_CMD = ">>gamma_vasp_cmd<<"
G16_CMD = ">>g16_cmd<<"
GAUSSIAN_GPU_CPU = ">>gaussian_gpu_cpu<<"
GAUSSIAN_CPU_LIST = ">>gaussian_cpu_list<<"
GAUSSIAN_MEMORY = ">>gaussian_memory<<"
GAUSSIAN_NPROCSHARED = ">>gaussian_nprocshared<<"
MA_DEF2 = ">>ma_def2<<"
PATH = ">>path<<"
RUNFILE_LOG = ">>runfile_log<<"
PROC_VM_KEY = ">>processing_vm_key<<"
LOCAL_OUTPUT_DIR = ">>local_output_dir<<"
ORCA_CMD = ">>orca_cmd<<"

_ENGINE_ENVIRONMENT = {
    "gaussian": {
        "g16_cmd": G16_CMD,
        "proc_vm_key": PROC_VM_KEY,
        "gpu_cpu": GAUSSIAN_GPU_CPU,
        "cpu_list": GAUSSIAN_CPU_LIST,
        "memory": GAUSSIAN_MEMORY,
        "nprocshared": GAUSSIAN_NPROCSHARED,
    },
    "orca": {"orca_cmd": ORCA_CMD},
    "pyscf": {},
}


def spec_or_task(task, fw_spec, key, default=None):
    """Resolve explicit FireWorks spec values before FireTask kwargs.

    ``False`` and ``0`` are valid overrides; only ``None`` means unspecified.
    """
    fw_spec = fw_spec or {}
    value = fw_spec.get(key)
    if value is None and task is not None:
        value = task.get(key)
    return default if value is None else value


def resolve_worker_value(value, fw_spec, default=None, required=False):
    """Resolve a ``>>worker_env_key<<`` placeholder from FireWorks' ``_fw_env``."""
    if not (isinstance(value, str) and value.startswith(">>") and value.endswith("<<")):
        return default if value is None else value

    name = value[2:-2]
    worker_env = (fw_spec or {}).get("_fw_env") or {}
    if not isinstance(worker_env, Mapping):
        raise TypeError("FireWorks _fw_env must be a mapping")
    resolved = worker_env.get(name)
    if resolved is None:
        if required:
            raise KeyError("FireWorker environment variable {!r} is required".format(name))
        return default
    return resolved


def resolve_runtime_field(task, fw_spec, key, default=None, required_env=False):
    value = spec_or_task(task, fw_spec, key, default)
    return resolve_worker_value(value, fw_spec, default=None, required=required_env)


def resolve_data_store(task=None, fw_spec=None):
    value = str(spec_or_task(task, fw_spec, "data_store", DEFAULT_DATA_STORE)).lower()
    if value not in DATA_STORES:
        raise ValueError("data_store must be one of {}, not {}".format(DATA_STORES, value))
    return value


@dataclass(frozen=True)
class WorkerRuntime:
    """Validated runtime settings resolved at a FireTask boundary.

    Scientific inputs remain in the task's ``paramset`` and FireWorks spec.
    Worker-specific paths, executables, and resources are resolved from the
    FireWorker's ``_fw_env`` using FireWorks' ``>>key<<`` convention.
    """

    software: Optional[str]
    data_store: str
    local_output_dir: Optional[str]
    path: Optional[str]
    runfile_log: Optional[str]
    engine_options: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_task(cls, task, fw_spec, software=None, require_executables=False):
        data_store = resolve_data_store(task, fw_spec)

        local_output_dir = resolve_runtime_field(
            task, fw_spec, "local_output_dir", LOCAL_OUTPUT_DIR
        )
        path = resolve_runtime_field(task, fw_spec, "path", required_env=True)
        runfile_log = resolve_runtime_field(task, fw_spec, "runfile_log")

        engine_options = {}
        for key, placeholder in _ENGINE_ENVIRONMENT.get(software, {}).items():
            default = "orca" if key == "orca_cmd" else None
            engine_options[key] = resolve_worker_value(
                spec_or_task(task, fw_spec, key, placeholder),
                fw_spec,
                default=default,
                required=require_executables and key in {"g16_cmd"},
            )

        return cls(
            software=software,
            data_store=data_store,
            local_output_dir=local_output_dir,
            path=path,
            runfile_log=runfile_log,
            engine_options=engine_options,
        )


def worker_task_kwargs(software):
    """Return FireTask references for common and selected-engine worker settings."""
    kwargs = {"path": PATH, "runfile_log": RUNFILE_LOG}
    for key, placeholder in _ENGINE_ENVIRONMENT.get(software, {}).items():
        kwargs[key] = placeholder
    return kwargs
