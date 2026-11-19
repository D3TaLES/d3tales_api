from types import SimpleNamespace

import pytest

from d3tales_api.Processors.parser_dft import ProcessDFTBase, ProcessGausLog, ProcessCCLIBMixin


class StubDFT(ProcessDFTBase):
    @property
    def runtime(self):
        return 0.25

    def get_radical_stability_score(self, spin_type=""):
        return {
            "radical_buried_vol": {"value": 1.0, "unit": "A^3"},
            "radical_spin": {"value": 0.5},
            "radical_stability_score": {"value": 0.75},
        }


def _stub(calc_type, metadata=None):
    obj = StubDFT(
        _id="M1",
        submission_info={"data_type": "gaussian"},
        metadata={"calculation_type": calc_type, "mol_info": {"groundState_charge": 0}, **(metadata or {})},
    )
    obj.code_name = "Gaussian"
    obj.code_version = "16"
    obj.functional = "B3LYP"
    obj.basis_set = "6-31G"
    obj.charge = 0
    obj.spin_multiplicity = 1
    obj.electrons = [2, 2]
    obj.final_energy = -10.0
    obj.dipole_moments = [(1.0, 2.0, 2.0)]
    obj.homo = -5.0
    obj.lumo = -1.0
    obj.homo_1 = -6.0
    obj.lumo_1 = 0.0
    obj.final_structure = [{"species": [{"element": "H"}]}]
    obj.gibbs_correction = 0.1
    obj.frequency_dicts = [{"frequency": 100.0}]
    obj.tddft_excitations = {"Singlet": [(2.0, 600.0, 0.2)], "Doublet": [], "Triplet": []}
    obj.tddft_spectra_data = {"energies": [1.0, 2.0]}
    obj.omega = 0.15
    return obj


@pytest.mark.unit
@pytest.mark.parametrize(
    "calc_type, expected_key",
    [
        ("opt_groundState", "scf_total_energy"),
        ("energy_gsc1", "scf_total_energy"),
        ("freq_groundState", "gibbs_correction"),
        ("tddft_groundState", "excitations"),
    ],
)
def test_dft_data_dict_shapes_calculation_payloads(calc_type, expected_key):
    record = _stub(calc_type).data_dict
    assert record["mol_id"] == "M1"
    assert record["runtime"] == 0.25
    assert record["data"][expected_key]
    assert record["data"]["number_of_electrons"] == 4
    assert record["data"]["is_groundState"] is True
    if calc_type == "freq_groundState":
        assert record["data"]["gibbs_correction"]["value"] == pytest.approx(2.72114)
    if calc_type == "opt_groundState":
        assert record["data"]["geometry"] == [{"species": [{"element": "H"}]}]
    if calc_type == "tddft_groundState":
        assert record["data"]["singlet_plotting"] == {"energies": [1.0, 2.0]}


@pytest.mark.unit
def test_dft_wtuning_payload_omits_standard_calculation_fields():
    record = _stub("wtuning").data_dict
    assert record["data"] == {
        "conditions": record["data"]["conditions"],
        "omega": 0.15,
    }
    assert "runtime" not in record


@pytest.mark.unit
def test_dft_conditions_include_optional_tuning_and_solvent():
    obj = _stub("opt_groundState", {"solvent": "Water", "dielectric_constant": 78.4})
    obj.tuning_parameter = 0.2
    assert obj.conditions == {
        "data_source": "dft",
        "code_name": "Gaussian",
        "code_version": "16",
        "functional": "B3LYP",
        "basis_set": "6-31G",
        "tuning_parameter": 0.2,
        "solvent": {"name": "Water", "model": "implicit_solvent", "dielectric_constant": 78.4},
    }


@pytest.mark.unit
def test_dft_ground_state_requires_ground_state_charge():
    obj = _stub("opt_groundState", {"mol_info": {}})
    with pytest.raises(IOError, match="groundState_charge"):
        _ = obj.is_groundState


@pytest.mark.unit
def test_cclib_frequency_extraction_and_missing_data_fallback():
    parsed = SimpleNamespace(
        vibfreqs=[100.0, 200.0],
        vibrmasses=[1.0, 2.0],
        vibfconsts=[3.0, 4.0],
        vibirs=[5.0, 6.0],
        vibsyms=["A", "B"],
    )
    assert ProcessCCLIBMixin.get_freq_dicts(parsed) == [
        {"frequency": 100.0, "r_mass": 1.0, "f_constant": 3.0, "IR_intensity": 5.0, "symmetry": "A"},
        {"frequency": 200.0, "r_mass": 2.0, "f_constant": 4.0, "IR_intensity": 6.0, "symmetry": "B"},
    ]
    assert ProcessCCLIBMixin.get_freq_dicts(SimpleNamespace()) == []


@pytest.mark.unit
def test_gaussian_runtime_tddft_dipole_and_pcm_parsing(tmp_path):
    logfile = tmp_path / "sample.log"
    logfile.write_text(
        " Job cpu time: 0 days 1 hours 30 minutes 0.0 seconds.\n"
        " Excitation energies and oscillator strengths:\n"
        " Excited State 1: Singlet-A 2.00 eV 620.0 nm f=0.25\n"
        " Excited State 2: Triplet-A 1.00 eV 1240.0 nm f=0.05\n"
        " Dipole moment (field-independent basis, Debye):\n"
        " X= 1.000 Y= -2.000 Z= 3.000 Tot= 3.742\n"
        " Polarizable Continuum Model\n"
        " Solvent: Water, Eps= 78.4 EpsInf= 1.7\n",
        encoding="utf-8",
    )
    obj = ProcessGausLog.__new__(ProcessGausLog)
    obj.log_path = str(logfile)

    assert obj.runtime == 1.5
    excitations = obj.get_tddft_excitations(str(logfile))
    assert excitations["Singlet"] == [(2.0, 620.0, 0.25)]
    assert excitations["Triplet"] == [(1.0, 1240.0, 0.05)]
    assert obj.dipole_moments == [(1.0, -2.0, 3.0)]
    assert obj.pcm_info == ["Water", 78.4]


@pytest.mark.unit
def test_gaussian_tddft_filters_negative_absorptions(monkeypatch, tmp_path):
    obj = ProcessGausLog.__new__(ProcessGausLog)
    obj.log_path = str(tmp_path / "unused.log")
    obj.calculation_type = "tddft_groundState"
    monkeypatch.setattr(
        obj,
        "get_tddft_excitations",
        lambda path: {"Singlet": [(-1.0, 100.0, 0.0), (2.0, 500.0, 0.1)], "Doublet": [], "Triplet": []},
    )
    obj.parse_tddft()
    assert obj.tddft_spectra_data
    assert obj.tddft_excitations["Singlet"] == [(-1.0, 100.0, 0.0), (2.0, 500.0, 0.1)]


@pytest.mark.unit
def test_gaussian_tuning_parameter_parses_both_route_parameter_forms():
    obj = ProcessGausLog.__new__(ProcessGausLog)
    obj.pmgmol = SimpleNamespace(route_parameters={"3108": "1300000"})
    assert obj.tuning_parameter == pytest.approx(0.3)
    obj.pmgmol.route_parameters = {"IOP(3107": "1200000,1"}
    assert obj.tuning_parameter == pytest.approx(0.2)
    obj.pmgmol.route_parameters = {}
    assert obj.tuning_parameter is None
