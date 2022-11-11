from d3tales_api.Calculators.calculators import *

data = {
    "concentration_smiles": "CC1(CCCC(N1[O])(C)C)C",
    "concentration_mass": "17.19mg",
    "solv_density": "1 g/mL",
    "concentration_volume": "10mL"
}

connector = {
            "smiles": "concentration_smiles",
            "weight": "concentration_mass",
            "solv_density": "solv_density",
            "volume": "concentration_volume"
        }
concentration = ConcentrationCalculator(connector=connector).calculate(data)
print(concentration)  # correct concentration = 0.01101 M
