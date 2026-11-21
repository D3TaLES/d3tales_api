# Workflows

The Workflows module provides the [Fireworks](https://materialsproject.github.io/fireworks/)-based workflows 
for the high-throughput [computational workflow](https://d3tales.as.uky.edu/docs/methodology/computational.html) 
that populated the [D<sup>3</sup>TaLES database](https://d3tales.as.uky.edu/database). 

## Modules 
Important FireWorks-based modules include:

* Calculation-specific modules for **FireTasks**: 
  * `Initialize`: FireTasks for initialize molecules and getting basic molecular information 
  * `Gaussian`: FireTasks for calculations with [Gaussian](https://gaussian.com/)
  * `Orca`: FireTasks for calculations with [ORCA](https://www.orcasoftware.de/)
  * `pyscf`: FireTasks for calculations with [PySCF](https://pyscf.org/)
* `Fireworks`: defines engine-aware **Fireworks** and routes them with FireWorks worker categories
* `wf_writer`: defines **Workflows**

The other modules are not of much importance for upper-level use because their functionalities are simply imported
by the Fireworks based modules. But, to summarize, these modules include: 
* `ParamSet`: translates dictionary-based calculation parameters (e.g., functional, basis set, etc.) for use in workflow
* `wtuning`: module adapted from the [OCELOT API](https://caer200.github.io/ocelot_api/) for performing IP fitting 
* `utils`: shared chemistry, geometry, storage, and result helpers

Gromacs is an empty placeholder and is not a supported workflow engine. VASP
parameter classes are deprecated; neither VASP nor ZINDO is registered in the
current workflow factory.

## Installing and configuring workers

Install the engine-neutral workflow dependencies with:

```bash
pip install "d3tales_api[workflows]"
```

For D3TaLES-backed data stores, add the database dependencies with
`pip install "d3tales_api[workflows,database]"`. Use `d3tales_web` on workers
outside the MongoDB firewall; it reads through the REST API and uploads
calculation files to the website. Use `d3tales_direct` only on workers with
network access and write permissions to the configured MongoDB database; it
processes results locally, writes computation documents directly, and archives
calculation files on the worker. The `local` store writes JSON results locally
and does not require database credentials or network access.

PySCF is optional: install it with `pip install "d3tales_api[workflows,pyscf]"`.
Gaussian and ORCA executables are provided by the compute environment rather
than pip.

The workflow uses FireWorks' existing separation of concerns:

* Chemistry and calculation inputs are stored in the parameter set and task
  arguments, then carried between Fireworks in `fw_spec` (for example charge,
  spin, solvent, and geometry).
* Engine selection sets each calculation Firework's `_category`; FireWorks
  schedules it on a worker with that category.
* Worker-specific paths, executable names, and Gaussian CPU/GPU resources are
  referenced with `>>environment_key<<` placeholders and resolved from the
  worker's `_fw_env`. Do not put node-specific resources in the chemistry
  parameter file.

For example:

```yaml
env:
  path: /scratch/calculations
  runfile_log: /scratch/runfiles.log
  g16_cmd: /path/to/g16
  gaussian_gpu_cpu: "0-1=0,2"
  gaussian_cpu_list: "0-3"
  gaussian_memory: "48GB"
  gaussian_nprocshared: "4"
  orca_cmd: /path/to/orca
```

Worker environment values can differ by worker without changing the saved
workflow or scientific parameters. Explicit workflow-spec values override
task arguments; missing optional resource values retain engine defaults.

The database connection is configured under the `d3tales` key in
`db_info.json`; the database name must be `d3tales`. Set `DB_INFO_FILE` on
workers using `d3tales_direct`. `d3tales_web` needs website API access and
upload credentials, but does not require a MongoDB connection. New molecules
must already be registered for `d3tales_web`; new-molecule registration is
available in `d3tales_direct`.

| `data_store` | Reads | Result handling |
| --- | --- | --- |
| `d3tales_web` | D3TaLES website REST API | Uploads Gaussian files for website-side processing; ORCA website uploads are unsupported. |
| `d3tales_direct` | Configured `d3tales` MongoDB database | Processes results on the worker, inserts calculation records directly, updates frontend characterizations for Gaussian computations (including omega from w-tuning), and stores local ZIP archives. ORCA/PySCF records are stored in `computation`; Gaussian-specific frontend characterization is not applied to their different record formats. |
| `local` | Local result files (with SMILES-derived fallback state) | Writes JSON records in the configured local results directory. |


## Setup
### Config Directory 

The directory where jobs are launched should also contain a `config`
directory. This should contain [fireworks configuration files](https://materialsproject.github.io/fireworks/config_tutorial.html), 
including `FW_config.yaml`, `my_fireworker.yaml`, and `my_launchpad.yaml`. The directory
`ex_config` provides an example for the config directory. Remember, the config 
directory must be in the same parent directory as `d3tales_fw`.

### Gaussian GPU jobs

Gaussian resource settings are optional FireWorker environment variables in
`my_fireworker.yaml`, so each HPC worker can provide values that match its
scheduler allocation. For example:

```yaml
env:
  g16_cmd: /path/to/g16
  gaussian_gpu_cpu: "0-1=0,2"
  gaussian_cpu_list: "0-3"
  gaussian_memory: "48GB"
  gaussian_nprocshared: "4"
```

`gaussian_gpu_cpu` supplies Gaussian's `%GPUCPU` value: the GPU indices,
followed by `=`, then one distinct control CPU per GPU. `gaussian_cpu_list`
supplies the complete allocated CPU list and must include those control CPUs.
Both lists must match the resources assigned by that worker's scheduler; the
indices above are examples only. When GPU variables are omitted, tasks use
CPU-only `%NProcShared` settings. `gaussian_memory` and `gaussian_nprocshared`
are optional; absent values retain the task's configured memory and default
CPU count. The separate w-tuning task does not use these GPU settings.

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

**Charge and spin state**: A task's `charge` is a relative charge, added to the molecule's ground-state charge.
The configured `multiplicity` is authoritative; it is not replaced by an automatic singlet/doublet estimate.
The workflow reports when its parity-derived estimate differs from the configured value and rejects a
multiplicity that is incompatible with the electron-count parity. This permits valid higher-spin states such as
triplets. Set `charge` and `multiplicity` in the task parameters, or override them with task kwargs or the
FireWorks spec. `restricted=False` selects an unrestricted reference; it does not change the configured spin.

**Geometry handoff**: A task configured with `geometry="X"` consumes `X_geom` (or `X_hash`) from an upstream
Firework. If that requested geometry is absent, the task fails rather than silently calculating from the SMILES
geometry. In builders that include a pre-w-tuning `opt_mol` job, that job always searches SMILES-derived RDKit
conformers, minimizes them with MMFF (or UFF when MMFF parameters are unavailable), and uses the lowest-energy
conformer as its starting geometry. This is independent of `geom_from_smiles` and any `mol_opt_geom` supplied to
initialization.


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
paramset = GausParamSet.from_json(param_file)

# Define Fireworks Workflow 
wf = d3tales_wf(paramset, smiles="CCCC", wtune=True, solvent='acetonitrile',
                data_store="local")

# Add the workflow to the Fireworks launchpad
info = LaunchPad().from_file(lpad_file).add_wf(wf)
```

## Launch jobs

Once jobs are added to the launchpad, you can launch jobs from the commandline of any server with the Fireworks
`rlaunch` commands. See the  [Fireworks documentation](https://materialsproject.github.io/fireworks/) for more detail. 
