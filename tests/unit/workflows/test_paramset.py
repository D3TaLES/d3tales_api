"""GausParamSet / GaussianParameters loading and serialization."""
import pytest

from d3tales_api.Workflows.ParamSet import GausParamSet, GaussianParameters, sett
from tests.helpers.workflow_helpers import PARAM_FILE

pytestmark = [pytest.mark.unit, pytest.mark.gaussian]


def test_from_json_builds_parameter_blocks(paramset):
    assert paramset.name == "D3TaLES_v01"
    assert isinstance(paramset.opt_groundState, GaussianParameters)
    assert paramset.tddft_groundState.route_parameters["td"]["nstates"] == 15


@pytest.mark.parametrize("block, charge, multiplicity", [
    ("opt_groundState", 0, 1), ("opt_anion1", -1, 2), ("opt_cation1", 1, 2),
    ("opt_anion2", -2, 1), ("opt_cation2", 2, 1), ("opt_triplet", 0, 3),
    ("energy_groundState", 0, 1), ("energy_anion1", -1, 2), ("energy_cation1", 1, 2),
    ("tddft_groundState", 0, 1), ("tddft_anion1", -1, 2), ("tddft_cation1", 1, 2),
])
def test_charge_and_multiplicity(paramset, block, charge, multiplicity):
    p = getattr(paramset, block)
    assert (p.charge, p.multiplicity) == (charge, multiplicity)


def test_missing_values_get_defaults(paramset):
    # "wtuning" has no multiplicity in the file
    assert paramset.wtuning.multiplicity == 1
    p = GaussianParameters()
    assert (p.charge, p.multiplicity, p.route_parameters, p.input_parameters, p.link0_parameters) == (0, 1, {}, {}, {})


def test_dict_round_trip(paramset):
    restored = GausParamSet.from_dict(paramset.as_dict())
    for block in ("opt_groundState", "wtuning", "tddft_cation1"):
        assert getattr(restored, block).as_dict() == getattr(paramset, block).as_dict()


def test_sett_accepts_none_dict_or_instance():
    proto = GaussianParameters()
    assert isinstance(sett(None, proto), GaussianParameters)
    assert sett({"charge": 1}, proto).charge == 1
    inst = GaussianParameters(charge=2)
    assert sett(inst, proto) is inst
    with pytest.raises(TypeError):
        sett("bad", proto)


def test_missing_file_raises(tmp_path):
    with pytest.raises((FileNotFoundError, OSError)):
        GausParamSet.from_json(str(tmp_path / "nope.json"))


def test_test_data_file_exists():
    assert PARAM_FILE.is_file()
