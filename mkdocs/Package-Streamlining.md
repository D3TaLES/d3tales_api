# Package streamlining changes

This note summarizes the changes that affect chemists preparing and launching
computational workflows with `d3tales_api`. The package retains its
`Workflows`, `Processors`, `Calculators`, and `D3database` areas, while workflow
chemistry, compute-worker configuration, and data storage have clearer
boundaries.

## Building a workflow

`GausParamSet` reads a JSON parameter file whose entries define each calculation
step (for example `opt_groundState`, `energy_anion1`, and `tddft_groundState`).
Set the charge and multiplicity for each species/step in that file. A workflow
builder assembles those calculations into a FireWorks `Workflow`; it does not
run calculations by itself.

```python
from fireworks import LaunchPad
from d3tales_api.Workflows.ParamSet import GausParamSet
from d3tales_api.Workflows.wf_writer import d3tales_wf

paramset = GausParamSet.from_json("parameter_file.json")

workflow = d3tales_wf(
    paramset,
    smiles="CCCC",
    wtune=True,
    solvent="acetonitrile",
    software="gaussian",
    data_store="local",
    name_tag="batch_",
    wf_tag="screen1_",
    check_if_already_run=True,
)

launchpad = LaunchPad.from_file("config/my_launchpad.yaml")
launchpad.add_wf(workflow)
```

The workflow is stored on the configured FireWorks LaunchPad. A FireWorks
launcher (`rlaunch`) on a configured worker then runs the ready jobs. The
worker's FireWorks category must match the engine selected by `software`.
For example, launch ready jobs from a configured worker with `rlaunch
rapidfire`.
Install the basic workflow dependencies with `pip install
"d3tales_api[workflows]"`; use `pip install
"d3tales_api[workflows,database]"` for either D3TaLES-backed store.

The principal `d3tales_wf` steps optimize the ground state, two anion states,
and two cation states; calculate vertical and solvent single-point energies;
and run TDDFT steps. `sams_wf` (the SAMS builder's actual function name) builds
a smaller ground-state/first-anion/first-cation workflow. It accepts a
precomputed molecular geometry for its optional w-tuning step:

```python
from d3tales_api.Workflows.wf_writer import sams_wf

workflow = sams_wf(
    paramset,
    smiles="O",
    mol_opt_geom=optimized_sites,
    wtune=True,
    software="gaussian",
    data_store="local",
)
```

`optimized_sites` is a serializable list of molecular site dictionaries, such
as the geometry format used in the workflow `*_geom` entries. If `wtune=True`,
`sams_wf` requires `mol_opt_geom` by default: w-tuning consumes `opt_mol`, and
this builder intentionally does not insert a geometry optimization to create
that input. You may explicitly pass `geom_from_smiles=True` to permit an
unoptimized SMILES-derived starting structure instead. With `wtune=False`, the
per-species optimizations can start from the provided geometry or generate a
starting structure from SMILES.

## Important workflow arguments

| Argument | Purpose |
| --- | --- |
| `paramset` | Calculation parameters loaded with `GausParamSet.from_json(...)`. |
| `smiles` / `identifier` | Identify the molecule. `smiles` is required to generate initial local structures and local ground-state data; use `identifier` when resolving a molecule already in the configured data store. |
| `mol_opt_geom` | Optional starting geometry in the serialized site-dictionary format. It seeds the `opt_mol` geometry; `sams_wf` requires it when `wtune=True` unless `geom_from_smiles=True` is explicitly selected. |
| `wtune` | Add the molecular optimization and omega-tuning precursor (`d3tales_wf`), or the omega-tuning precursor (`sams_wf`). Set false to omit w-tuning. |
| `hf_mol_opt` | Use the configured HF molecular-optimization parameter set where the builder supports this option. |
| `solvent` | Solvent model/name for solvent optimizations or single points. `d3tales_wf` defaults to acetonitrile; `None` requests gas-phase behavior where supported. |
| `software` | Selects a supported engine: `gaussian`, `orca`, or `pyscf`. PySCF must be installed separately. |
| `data_store` | Choose `local`, `d3tales_web`, or `d3tales_direct` according to the worker's network access and desired persistence. |
| `name_tag` | Prefix Firework names and geometry keys. Use it to distinguish repeated/batched workflows. |
| `wf_tag` | Prefix the workflow's display name; it does not change calculation names or geometry keys. |
| `check_if_already_run` | Reuse existing calculations/results when the engine task supports the cache check. |
| `restricted` | Select restricted (`True`) or unrestricted (`False`) reference behavior where supported. It does not set the spin multiplicity. |
| `priority` | Set FireWorks priority for the generated calculation jobs. |
| `geom_from_smiles` | Explicitly allow a calculation to fall back to a SMILES-derived starting structure when its requested geometry is absent. Use only when that fallback is scientifically intended. |

Other task options, including `run_nto`, `origin_group`, and `public`, can be
passed as keyword arguments where applicable. FireWorks spec values override
task arguments, including explicit `False` and `0` values.

## Chemistry and geometry correctness fixes

- **Configured charge and spin are preserved.** The configured relative charge
  and multiplicity are authoritative; the workflow no longer replaces an
  explicitly configured triplet or other higher-spin state with a
  singlet/doublet-only parity estimate. It checks electron-count parity,
  rejects incompatible combinations, and warns when the basic parity-based
  estimate differs from the configured spin. `restricted=False` controls the
  reference method and is not a substitute for setting multiplicity.
- **Requested structures must be available.** A task requesting geometry
  `X` consumes `X_geom` (or its stored hash) from an upstream Firework. It no
  longer silently performs an unrelated calculation on a SMILES-derived
  structure when that requested geometry is missing. SMILES-based geometry
  generation remains an explicit starting-geometry path for initial
  optimizations; an intentional fallback can be enabled with
  `geom_from_smiles=True`.
- **Geometry flow follows the intended species and solvent.** Vertical
  anion/cation energies consume the matching optimized species geometry;
  solvent TDDFT requests the solvent-optimized geometry, while gas-phase
  TDDFT requests the gas-phase geometry. This also fixes mismatched parents
  and consumer keys in public builders such as `just_anion`, `huckaba_wf`,
  and `just_nmr`.
- **Name-tagged workflows keep producer and consumer keys aligned.** Geometry
  lookups prefer the matching `name_tag` key, so a batched workflow does not
  accidentally consume another calculation's untagged structure.
- **W-tuning forwards required molecular state.** Its downstream calculations
  receive the geometry and molecule inputs needed for consistent execution.
  The SAMS builder now fails at construction with a clear error if w-tuning is
  requested without its required `mol_opt_geom` or an explicit SMILES fallback,
  rather than building a
  workflow that fails later when the tuning task looks up `opt_mol`.

These contracts are covered by tests of the public workflow graphs, including
geometry producers, consumers, parent dependencies, solvent choices, and
w-tuning enabled/disabled cases.

## FireWorks and worker configuration

Scientific inputs (parameter sets, charge, multiplicity, solvent, and geometry)
belong to calculation tasks and `fw_spec`. The engine is selected independently
with `software`; the calculation Firework's `_category` routes it to a
compatible worker. Paths, executable locations, and hardware resources belong
in each worker's `my_fireworker.yaml`, not in chemistry parameter files:

```yaml
env:
  path: /scratch/calculations
  runfile_log: /scratch/runfiles.log
  g16_cmd: /opt/gaussian/g16
  gaussian_gpu_cpu: "0-1=0,2"
  gaussian_cpu_list: "0-3"
  gaussian_memory: "48GB"
  gaussian_nprocshared: "4"
  orca_cmd: /opt/orca/orca
```

The `gaussian_gpu_cpu` and `gaussian_cpu_list` examples must be changed to
match the resources allocated to that worker. GPU values are optional; without
them, Gaussian uses CPU-only settings. Worker environment values are resolved
at task runtime through FireWorks' `_fw_env`, so the same workflow can be
launched on workers with different resources.

## Database and result storage modes

Configure the single MongoDB connection under the `d3tales` key in
`db_info.json` (with `"database": "d3tales"`). The old `frontend` and `backend`
connection entries are no longer used. `FrontDB` and `BackDB` remain
deprecated compatibility wrappers over the unified database and emit a
`DeprecationWarning`; use `D3Database(database="d3tales", collection_name=...)`
for new code.

| Store | Network location | Molecule reads | Calculation result handling |
| --- | --- | --- | --- |
| `d3tales_web` | Workers can reach the D3TaLES website but are outside the DB firewall. | Website REST API. Molecules must be pre-registered; this API does not expose a workflow registration endpoint. | Gaussian files are zipped and uploaded for website-side processing. ORCA is rejected because the website upload endpoint expects Gaussian output. PySCF remains local because website processing expects Gaussian logs. |
| `d3tales_direct` | Workers are inside the DB firewall and have direct MongoDB write access. | `base` and `computation` collections in the configured `d3tales` database; no REST calls. New molecule registration is supported. | Results are processed on the worker and inserted into `computation`. Engine output/checkpoint files are zipped and retained locally; no website upload is performed. |
| `local` | No database or website access required. | Local molecule data/results, with explicit SMILES-based starting geometry fallback. | Processed results are stored as per-calculation JSON files. |

Website uploads in `d3tales_web` are always enabled; selecting
`d3tales_direct` means the workflow writes results to MongoDB. A local
`db_info.json` is intentionally not required by tests: live credentials and
database writes are never used by the direct-store unit tests.
The former `data_store="d3tales"` value is rejected; choose `d3tales_web` or
`d3tales_direct` explicitly.

## Package boundaries and dependencies

Database client classes have one canonical implementation in `D3database`;
the workflow-facing adapter loads database support lazily. PySCF and
ChemDataExtractor2 are optional dependencies. Packaging metadata, optional
dependency groups, and pytest configuration are consolidated in
`pyproject.toml`; `requirements.txt` remains available for existing
environments. Wildcard imports and dynamic `eval`-based workflow parameter
selection were replaced with explicit imports and lookups. VASP parameter
classes are deprecated, and Gromacs, VASP, and ZINDO are not supported workflow
engines.
