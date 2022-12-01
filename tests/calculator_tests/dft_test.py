"""
this example shows the usages of calculators with schema2class
Easy to work with for objects with different schema
The code of calculator hardly changes

"""

from d3tales_api.Calculators.plotters import *
from d3tales_api.Calculators.calculators import *

# Get DFT data
with open('../raw_data/gaus_ex_data.json', 'r') as fn:
    dft_data = json.load(fn)
# Update with redox potential parameters
dft_data.update({"electrode": "standard_hydrogen_electrode", "num_electrons": 1})
# Update with TDDFT parameters
dft_data.update({"sigma": 0.10, "step": 0.01})

connector = {
    # Energy difference connections
    "energy_final": "solv_energy_gsgs.scf_total_energy.value",
    "energy_initial": "opt_groundState.scf_total_energy.value",

    # Reorganization connections
    "gs_opt": "opt_groundState.scf_total_energy.value",
    "ion_opt": "opt_cation1.scf_total_energy.value",
    "gs_energy": "energy_c1gs.scf_total_energy.value",
    "ion_energy": "energy_gsc1.scf_total_energy.value",

    # Relaxation connections
    "opt_energy": "opt_cation1.scf_total_energy.value",
    "energy": "energy_gsc1.scf_total_energy.value",

    # RMSD connections
    "geom_initial": "opt_groundState.geometry",
    "geom_final": "opt_cation1.geometry",

    # Reduction potential connections
    "init_eng": "opt_groundState.scf_total_energy.value",
    "init_corr": "freq_groundState.gibbs_correction.value",
    "init_eng_solv": "solv_energy_gsgs.scf_total_energy.value",
    "fin_eng": "opt_anion1.scf_total_energy.value",
    "fin_corr": "freq_anion1.gibbs_correction.value",
    "fin_eng_solv": "solv_energy_a1a1.scf_total_energy.value",
    "electrode": "electrode", "num_electrons": "num_electrons",

    # TDDFT plotting connections
    "transitions": "tddft_groundState.excitations.Singlet",
    "sigma": "sigma", "step": "step"
}
solv_eng = EnergyDiffCalc(connector=connector).calculate(dft_data)
assert solv_eng == -0.1879

reorg_energy = ReorganizationCalc(connector=connector).calculate(dft_data)
assert reorg_energy == 0.2643

relax_energy = RelaxationCalc(connector=connector).calculate(dft_data)
assert relax_energy == 0.1248

dist = RMSDCalc(connector=connector).calculate(dft_data)
assert dist == 0.05498

reduction = RedoxPotentialCalc(connector=connector).calculate(dft_data)
assert reduction == 3.77

plotting = DFTSpecPlotter(connector=connector).plot_data(dft_data)
assert list(plotting.keys()) == ['energies', 'lambda', 'xas', 'abs_plot']

# Stability score testing
c_data = {"log_file": '../raw_data/energy_cation1.log', "spin_type": "mulliken"}
connector = {"log_file": "log_file", "spin_type": "spin_type"}
rss = RSSCalc(connector=connector).calculate(c_data)
assert rss == 78.55

print("DFT TESTING SUCCESSFUL")
