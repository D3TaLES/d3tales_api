# Calculators

D<sup>3</sup>TaLES API Calculators allow users to calculate useful experimental and 
computational properties from nested data. Additionally, all calculators contain unit conversion
features. 

Full documentation can be found [here](d3tales_api.Calculators.html).

### Calculator Submodules

The Calculators module contains two primary submodules (`calculators` and `plotters`) and a minor submodules
(`utils`). Most users will only the primary submodules, specifically the `calculators` primary submodule. 

Submodules: 
* `calculators`: This contains many calculator classes. Each calculator class must receive 
`data` and a data `connector` that connects the input data with the connection points. 
Generally, each calculator calculates a specific property with a `calculate` method. 
* `plotters`: The submodule is adapted from the calculators submodule and focuses on producing 
plotting data. Most plotters have a `plot_data` method that produces data for a [plotly](https://plotly.com/)
plot and a `live_plot` method that produces a [Matplotlib](https://matplotlib.org/) plot.  
* `utils`: This contains several useful functions used in the primary submodules, most notable
of which is the `unit_conversion` function. 

## Basic Use

The following examples show the basic use for the D<sup>3</sup>TaLES Calculators, including 
examples for both computational and experimental data. In general, calculating a property requires: 
1. Importing the appropriate calculator class, `CALCULATOR_CLASS`
2. Defining the raw data, `DATA`
3. Defining a connector from the raw data structure to the calculator connection points (see documentation for the imported 
calculator class), `CONNECTOR`
4. Establishing the calculator instance and calculating the property. The general format 
for this final step is: `CALCULATOR_CLASS(connector=CONNECTOR).calculate(DATA)`. 

Generally, calculating the property from a calculators class instance requires calling the 
`calculate` method. However, there are some cases (EX: `CVDescriptorCalculator`) where several
properties can be calculated from one class, so more specific methods must be called(EX: 
`peaks`, `reversibility`, `e_half`). 

Note the `plotter` submodule has a very similar use to the `calculators` submodule. However, 
instead of calling the `calculate` method, the user must call the `plot_data` or `live_plot` method. 

### Molecular DFT

This example shows the calculation of the adiabatic ionization energy from `calc_data` extracted
from a calculation output file and transformed to the D<sup>3</sup>TaLES backend schema. 
```python
from d3tales_api.Calculators.calculators import EnergyDiffCalc

calc_data = {
  "opt_groundState": {
    "scf_total_energy": {
      "value": -0.60, 
      "unit": "eV"
    }
  },
  "opt_cation1": {
    "scf_total_energy": "-0.25eV"
  }
}

connector = {   
  "energy_final":   "opt_cation1.scf_total_energy",
  "energy_initial":  "opt_groundState.scf_total_energy"
}
energy = EnergyDiffCalc(connector=connector).calculate(calc_data)
```
The resulting `energy` should be `0.35`. 

Useful molecular DFT calculations include: 
* `ReorganizationCalc`: Reorganization energy.
* `RMSDCalc`: Root mean squared error between atoms in two geometries.  
* `DeltaGSolvCalc`: Change in Gibbs energy for solvation.
* `RedoxPotentialCalc`: Redox potential.
* `RadBuriedVolCalc`: Radical buried volume for the atom with the atom with the highest portion
        of spin. Uses [DBSTEP](https://github.com/patonlab/DBSTEP).
* `RadicalSpinCalc`: Radical spin density. Uses [DBSTEP](https://github.com/patonlab/DBSTEP).

For further examples of the D<sup>3</sup>TaLES API Calculators used for molecular DFT calculations, 
see this Colab notebook: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1q2O6-76WaGRU81g0QwwJqRQIDtfiYhD3?usp=sharing)

### Cyclic Voltammetry

CV scan data is assumed to be in the format produced by the D<sup>3</sup>TaLES CV parsers, 
a list of scans where each list contains a list of points and each point is a tuple containing
the voltage the potential. EX: `[[[volt_1a, pot_1a], [volt_1b, pot1b], ...], [[volt_2a, pot_2a], 
[volt_2b, pot2b], ...], ...]`
```python
SCAN_DATA1 = [] # Some extracted CV scan data
SCAN_DATA2 = [] # Some extracted CV scan data
SCAN_DATA3 = [] # Some extracted CV scan data

cv_data = [
    {
        "conditions": {
            "scan_rate": "0.15V/s",
            "working_electrode_surface_area": {
                "value": "0.007",
                "unit": "cm^3"
            },
            "redox_mol_concentration": "0.01M"
        }, 
        "num_electrons": 1,
        "scan_data": SCAN_DATA1
    },
    {
        "conditions": {
            "scan_rate": "0.15V/s",
            "working_electrode_surface_area": "0.007cm^3",
            "redox_mol_concentration": "0.01M"
        }, 
        "num_electrons": 1,
        "scan_data": SCAN_DATA2
    },
    {
        "conditions": {
            "scan_rate": "0.15V/s",
            "working_electrode_surface_area": "0.007cm^3",
            "redox_mol_concentration": "0.01M"
        }, 
        "num_electrons": 1,
        "scan_data": SCAN_DATA3
    },
]

connector = {
    "A": "conditions.working_electrode_surface_area",
    "v": "conditions.scan_rate",
    "C": "conditions.redox_mol_concentration",
    "n": "num_electrons",
    "scan_data": "scan_data"
}
```

#### General CV Analysis

Extracting basic cyclic voltammagram properties from the first
CV in the `cv_data` list. 
```python
from d3tales_api.Calculators.calculators import CVDescriptorCalculator

descriptor_cal = CVDescriptorCalculator(connector=connector)
print("Peaks: ", descriptor_cal.peaks(cv_entries[0]))
print("Reversibility: ", descriptor_cal.reversibility(cv_entries[0]))
print("E 1/2: ", descriptor_cal.e_half(cv_entries[0]))
print("Peak Splittings: ", descriptor_cal.peak_splittings(cv_entries[0]))
print("Numer of middle sweeps: ", len(descriptor_cal.middle_sweep(cv_entries[0])))
```

#### Diffusion Coefficient Calculation

Using the CVs in the `cv_data` list to calculate diffusion coefficient 
```python
from d3tales_api.Calculators.calculators import CVDiffusionCalculator

diffusion_cal = CVDiffusionCalculator(connector=connector)
diffusion = diffusion_cal.calculate(cv_data)
print("Average diffusion", diffusion[0])
print("Fitted diffusion", diffusion[1])
```

Useful CV calculations include: 
* `ConcentrationCalculator`: Solution concentration.
* `CVDescriptorCalculator`: Various CV properties including peaks, reversibility, E<sub>1/2</sub>,
peak splitting, etc.
* `CVDiffusionCalculator`: Diffusion constant using Randles-Scidwick equation.
* `CVChargeTransferCalculator`: Charge transfer rate.

For further examples of the D<sup>3</sup>TaLES API Calculators used for CV analysis, 
see this Colab notebook: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/12E5Go9KujNtify6eUTx97WabZL7axU61?usp=share_link)

## Considering Units

The D<sup>3</sup>TaLES Calculators use [pint](https://pint.readthedocs.io/en/stable/) to 
manage unit conversions. Therefore, parameter units can be specified in any manner 
compatible with pint (EX: `1 molar`, `1.0M`, pint objects). Additionally, these calculators 
accept the D<sup>3</sup>TaLES `data_dict` format: `{"value": 1, "unit": "M"}`. 