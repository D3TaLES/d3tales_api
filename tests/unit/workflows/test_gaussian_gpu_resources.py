import pytest

pytest.importorskip("fireworks")

from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

from d3tales_api.Workflows.Gaussian import RunGaussianEnergy
from d3tales_api.Workflows.ParamSet import GaussianParameters
from d3tales_api.Workflows.utils import (
    GAUSSIAN_CPU_LIST,
    GAUSSIAN_GPU_CPU,
    GAUSSIAN_MEMORY,
    GAUSSIAN_NPROCSHARED,
)


def _task(**kwargs):
    kwargs.setdefault("paramset", GaussianParameters(link0_parameters={"%mem": "12GB"}))
    task = RunGaussianEnergy(**kwargs)
    task.paramset = task["paramset"]
    return task


def _configured_parameters(task, env):
    task.configure_link0_resources({"_fw_env": env})
    return task.paramset.link0_parameters


def test_worker_gpu_settings_build_matching_link0_directives():
    task = _task(
        gpu_cpu=GAUSSIAN_GPU_CPU,
        cpu_list=GAUSSIAN_CPU_LIST,
        memory=GAUSSIAN_MEMORY,
        nprocshared=GAUSSIAN_NPROCSHARED,
    )

    directives = _configured_parameters(task, {
        "gaussian_gpu_cpu": "0-1=0,2",
        "gaussian_cpu_list": "0-3",
        "gaussian_memory": "24GB",
        "gaussian_nprocshared": "8",
    })

    assert directives == {"%mem": "24GB", "%CPU": "0-3", "%GPUCPU": "0-1=0,2"}

    gaussian_input = GaussianInput(
        Molecule(["H"], [[0, 0, 0]], charge=0, spin_multiplicity=2),
        route_parameters={"#p": "HF/3-21G"},
        link0_parameters=directives,
    )
    rendered = str(gaussian_input)
    assert "%CPU=0-3" in rendered
    assert "%GPUCPU=0-1=0,2" in rendered


def test_missing_gpu_settings_keep_cpu_only_directives():
    task = _task(
        gpu_cpu=GAUSSIAN_GPU_CPU,
        cpu_list=GAUSSIAN_CPU_LIST,
        memory=GAUSSIAN_MEMORY,
        nprocshared=GAUSSIAN_NPROCSHARED,
    )

    directives = _configured_parameters(task, {})

    assert directives["%mem"] == "12GB"
    assert directives["%nprocshared"].isdigit()
    assert "%CPU" not in directives
    assert "%GPUCPU" not in directives


@pytest.mark.parametrize("env, message", [
    ({"gaussian_gpu_cpu": "0-1=0", "gaussian_cpu_list": "0-3"}, "one control CPU per GPU"),
    ({"gaussian_gpu_cpu": "0=2", "gaussian_cpu_list": "0-1"}, "include every GPU control CPU"),
    ({"gaussian_gpu_cpu": "0=0", "gaussian_cpu_list": "0,0"}, "duplicate indices"),
    ({"gaussian_gpu_cpu": "0=0"}, "cpu_list must be set"),
])
def test_invalid_worker_gpu_assignments_are_rejected(env, message):
    task = _task(gpu_cpu=GAUSSIAN_GPU_CPU, cpu_list=GAUSSIAN_CPU_LIST)

    with pytest.raises(ValueError, match=message):
        _configured_parameters(task, env)


def test_worker_cpu_list_can_configure_cpu_only_job():
    task = _task(cpu_list=GAUSSIAN_CPU_LIST)

    directives = _configured_parameters(task, {"gaussian_cpu_list": "0-7"})

    assert directives["%CPU"] == "0-7"
    assert "%GPUCPU" not in directives
    assert "%nprocshared" not in directives
