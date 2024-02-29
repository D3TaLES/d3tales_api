# Workflows

The Workflows module provides the [Fireworks](https://materialsproject.github.io/fireworks/)-based workflows 
for the high-throughput [computational workflow](https://d3tales.as.uky.edu/docs/methodology/computational.html) 
that populated the [D<sup>3</sup>TaLES database](https://d3tales.as.uky.edu/database). 

## Modules 
Important Fireworks-based modules include: 

* Calculation-specific modules for **FireTasks**: 
  * `Initialize`: FireTasks for initialize molecules and getting basic molecular information 
  * `Gaussian`: FireTasks for Density Functional Theory calculations with [Gaussian](https://gaussian.com/)
  * _In development_ `Gromacs`: FireTasks for Molecular Dynamics simulations with [Gromacs](https://www.gromacs.org/)
* `D3TaLES_FW`: defines all **Fireworks** 
* `wf_writer`: defines **Workflows**

The other modules are not of much importance for upper-level use because their functionalities are simply imported
by the Fireworks based modules. But, to summarize, these modules include: 
* `ParamSet`: translates dictionary-based calculation parameters (e.g., functional, basis set, etc.) for use in workflow
* `wtuning`: module adapted from the [OCELOT API](https://caer200.github.io/ocelot_api/) for performing IP fitting 
* `utils`: utility functions and environment variables 


## Setup
### Config Directory 

The directory where jobs are launched should also contain a `config`
directory. This should contain [fireworks configuration files](https://materialsproject.github.io/fireworks/config_tutorial.html), 
including `FW_config.yaml`, `my_fireworker.yaml`, and `my_launchpad.yaml`. The directory
`ex_config` provides an example for the config directory. Remember, the config 
directory must be in the same parent directory as `d3tales_fw`.

### Parameter file
The `GausParamSet` class generates object-based parameters for the Workflows module to use throughout the design and launch
of calculation workflows. The object parameters, `GausParamSet` reads a JSON file containing job parameters. The first-level
keys in this parameter file should be the job names and the values should be the associated parameters. For example, 
a JSON file with the following data would establish DFT parameters for `opt_groundState` and `opt_cation1` jobs. 
```json
{
   "opt_groundState": {
        "route_parameters":{
            "opt": "",
            "scf(maxcyc=512)" : ""
        },
        "functional": "LC-wHPBE",
        "basis_set": "Def2SVP",
        "charge": 0,
        "multiplicity": 1
    },

  "opt_cation1": {
      "route_parameters":{
          "opt": "",
          "scf(maxcyc=512)" : ""
      },
      "functional": "LC-wHPBE",
      "basis_set": "Def2SVP",
      "charge": 1,
      "multiplicity": 2
  }
}
```

**Note on charge and spin state**: The charge and multiplicity set in the parameters file assumes a non-radical, 
neutrally-charged ground state. If a molecule has a charged or radical ground state that is represented in the 
associated SMILES string, the Workflow code will automatically adjust the calculation charge and multiplicity. 


## Adding a workflow

The following code submits a basic workflow. 

```python
import os
from fireworks.core.launchpad import LaunchPad
from d3tales_api.Workflows.wf_writer import d3tales_wf
from d3tales_api.Workflows.ParamSet import GausParamSet


BASE_DIR = "" # BASE_DIR must be defined 
lpad_file = os.path.join(BASE_DIR, 'config', 'my_launchpad.yaml')  # See notes about config directory 
param_file = os.path.join(BASE_DIR, 'parameter_file.json')

# Establish calculation parameters from parm_file json file
paramset = GausParamSet().from_json(param_file)

# Define Fireworks Workflow 
wf = d3tales_wf(paramset, smiles="CCCC", wtune=True, solvent='acetonitrile')

# Add the workflow to the Fireworks launchpad
info = LaunchPad().from_file(lpad_file).add_wf(wf)
```

## Launch jobs

Once jobs are added to the launchpad, you can launch jobs from the commandline of any server with the Fireworks
`rlaunch` commands. See the  [Fireworks documentation](https://materialsproject.github.io/fireworks/) for more detail. 
