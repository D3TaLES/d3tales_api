"""Engine setup preserves configured charge, spin, and reference overrides."""
import pytest

from d3tales_api.Workflows.ParamSet import GaussianParameters
from d3tales_api.Workflows.Gaussian import RunGaussianEnergy
from d3tales_api.Workflows.Orca import RunOrcaEnergy, get_orca_omega, orca_conditions, orca_hash_id
from d3tales_api.Workflows.pyscf import RunPySCFEnergy, pyscf_conditions, split_functional
from d3tales_api.Workflows.utils import WTUNING_KEY, local_record, orig_hash_id, save_local_result

pytestmark = [pytest.mark.unit, pytest.mark.gaussian]


@pytest.mark.parametrize("engine", ["gaussian", "orca", "pyscf"])
def test_engine_setup_keeps_configured_triplet_and_unrestricted_override(
    engine, tmp_path, monkeypatch, water_geom
):
    monkeypatch.chdir(tmp_path)
    paramset = GaussianParameters(
        route_parameters={}, charge=0, multiplicity=3,
        basis_set="def2-SVP", functional="LC-wHPBE",
    )
    task_kwargs = {
        "name": "energy_groundState",
        "paramset": paramset,
        "path": str(tmp_path),
        "geometry": "opt_groundState",
        "data_store": "local",
        "restricted": True,
        "skip_freq": True,
        "use_iop": True,
    }
    if engine == "gaussian":
        task_kwargs["g16_cmd"] = "g16"
        task = RunGaussianEnergy(**task_kwargs)
    elif engine == "orca":
        task = RunOrcaEnergy(**task_kwargs)
    else:
        task_kwargs["pyscf_file_name"] = "pyscf"
        task = RunPySCFEnergy(**task_kwargs)

    fw_spec = {
        "identifier": "M1",
        "smiles": "O",
        "gs_charge": 0,
        "gs_spin": 1,
        "data_store": "local",
        "geometry": "opt_groundState",
        "opt_groundState_geom": water_geom,
        "restricted": False,
        "skip_freq": False,
        "use_iop": False,
    }
    with pytest.warns(RuntimeWarning, match="preserving the configured value"):
        action = task.setup_calc(fw_spec, calc_type="energy")

    assert action is None
    assert task.restricted is False
    assert task.skip_freq is False
    if engine == "gaussian":
        assert (task.paramset.charge, task.paramset.multiplicity) == (0, 3)
    else:
        assert (task.charge, task.multiplicity) == (0, 3)
    assert paramset.charge == 0 and paramset.multiplicity == 3

    if engine == "gaussian":
        assert task.functional == "ULC-wHPBE"
        assert task.paramset.functional == "ULC-wHPBE"
    elif engine == "orca":
        assert "UKS" in task.input_text("SP")
    else:
        assert task.settings.prefix == "U"


@pytest.mark.unit
def test_orca_records_use_orca_specific_conditions_and_hashes():
    conditions = orca_conditions("B3LYP", "def2-SVP")
    assert conditions["code_name"] == "ORCA"
    assert orca_hash_id("M1", "energy", "B3LYP", "def2-SVP") != orig_hash_id(
        "M1", "energy", "B3LYP", "def2-SVP"
    )
    assert pyscf_conditions("B3LYP", "def2-SVP")["code_name"] == "PySCF"
    assert split_functional("ULC-wHPBE")[1] == "U"


@pytest.mark.unit
def test_orca_uses_only_orca_wtuning_results(tmp_path):
    save_local_result(
        str(tmp_path),
        local_record(WTUNING_KEY, "M1", "wtuning", {"omega": 0.1}),
    )
    assert get_orca_omega(str(tmp_path)) is None

    save_local_result(
        str(tmp_path),
        local_record(WTUNING_KEY, "M1", "wtuning", {"engine": "ORCA", "omega": 0.2}),
    )
    assert get_orca_omega(str(tmp_path)) == pytest.approx(0.2)
