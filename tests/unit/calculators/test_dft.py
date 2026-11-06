import copy
import pytest
from d3tales_api.Calculators.calculators import (
    DeltaGSolvCalc, EnergyDiffCalc, RedoxPotentialCalc,
    RelaxationCalc, ReorganizationCalc, RMSDCalc,
)

CONNECTOR = {
 "energy_final":"solv_energy_gsgs.scf_total_energy.value",
 "energy_initial":"opt_groundState.scf_total_energy.value",
 "gs_opt":"opt_groundState.scf_total_energy.value",
 "ion_opt":"opt_cation1.scf_total_energy.value",
 "gs_energy":"energy_c1gs.scf_total_energy.value",
 "ion_energy":"energy_gsc1.scf_total_energy.value",
 "opt_energy":"opt_cation1.scf_total_energy.value",
 "energy":"energy_gsc1.scf_total_energy.value",
 "geom_initial":"opt_groundState.geometry",
 "geom_final":"opt_cation1.geometry",
 "init_eng":"opt_groundState.scf_total_energy.value",
 "init_corr":"freq_groundState.gibbs_correction.value",
 "init_eng_solv":"solv_energy_gsgs.scf_total_energy.value",
 "fin_eng":"opt_anion1.scf_total_energy.value",
 "fin_corr":"freq_anion1.gibbs_correction.value",
 "fin_eng_solv":"solv_energy_a1a1.scf_total_energy.value",
 "electrode":"electrode", "num_electrons":"num_electrons",
}

@pytest.fixture
def prepared_gaussian_data(gaussian_data):
    data=copy.deepcopy(gaussian_data)
    data.update({"electrode":"standard_hydrogen_electrode", "num_electrons":1})
    return data

@pytest.mark.unit
@pytest.mark.gaussian
@pytest.mark.parametrize("calculator, expected, tolerance", [
    (EnergyDiffCalc, -0.1879, 1e-4),
    (ReorganizationCalc, 0.2643, 1e-4),
    (RelaxationCalc, 0.1248, 1e-4),
    (RMSDCalc, 0.05498, 1e-2),
    (RedoxPotentialCalc, 3.77, 1e-2),
])
def test_legacy_dft_regressions(prepared_gaussian_data, calculator, expected, tolerance):
    result=calculator(connector=CONNECTOR).calculate(prepared_gaussian_data)
    assert float(result) == pytest.approx(expected, abs=tolerance)
