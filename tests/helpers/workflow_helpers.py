"""
Test-only helpers for building and submitting D3TaLES DFT workflows.

These are self-contained stand-ins for the kind of launcher scripts users write on top of
``d3tales_api.Workflows``. They are not imported from any script; they exist so the tests can
exercise the package's real Firework factories, ParamSet loading, and LaunchPad submission.
"""
import json
import os
from pathlib import Path

from fireworks import Workflow

from d3tales_api.Workflows.Fireworks import InitializeMolecule, WTuning, Optimization, Energy, TDDFT
from d3tales_api.Workflows.ParamSet import GausParamSet

TESTS_DIR = Path(__file__).resolve().parents[1]
GAUSSIAN_DATA_DIR = TESTS_DIR / "data" / "gaussian"
PARAM_FILE = GAUSSIAN_DATA_DIR / "gaus_parameter_file.json"
MOLECULE_FILE = GAUSSIAN_DATA_DIR / "workflow_populate_data.json"

INIT_NAME = "initial_molecule_data"
WTUNE_NAME = "wtuning"


def example_dft_workflow(paramset, identifier=None, smiles=None, mol_opt_geom=None, wtune=True, wf_tag="",
                         species=("groundState", "anion1", "cation1"), **kwargs):
    """
    Example multi-step workflow:

        init -> [wtuning] -> opt_<species> -> tddft_<species>
                                           -> energy at the ground-state geometry (charged species)
                             opt_<charged> -> ground-state energy at the charged geometry
    """
    init = InitializeMolecule(identifier=identifier, smiles=smiles, mol_opt_geom=mol_opt_geom, **kwargs)
    fws = [init]
    opt_parent = init
    if wtune:
        opt_parent = WTuning(paramset=paramset.wtuning, parents=init, **kwargs)
        fws.append(opt_parent)

    opts = {}
    for sp in species:
        opts[sp] = Optimization(paramset=getattr(paramset, "opt_" + sp), species=sp, parents=opt_parent, **kwargs)
        fws.append(opts[sp])

    for sp in species:
        if sp == "groundState":
            continue
        fws.append(Energy(paramset=paramset.energy_groundState, species="groundState", geometry="opt_" + sp,
                          parents=[opts[sp]], **kwargs))
        if "groundState" in opts:
            fws.append(Energy(paramset=getattr(paramset, "energy_" + sp), species=sp, geometry="opt_groundState",
                              parents=[opts["groundState"]], **kwargs))

    for sp in species:
        fws.append(TDDFT(paramset=getattr(paramset, "tddft_" + sp), species=sp, parents=[opts[sp]], **kwargs))

    return Workflow(fws, name="{}gaus_{}".format(wf_tag, identifier or smiles))


def load_molecules(filename):
    """{smiles: geometry sites}; empty when the file is missing."""
    if filename and os.path.isfile(filename):
        with open(filename, "r") as f:
            return json.load(f)
    return {}


def submit_workflows(launchpad, molecule_file, param_file, workflow_function=example_dft_workflow, priority=4,
                     public=False, origin_group="Risko", submit=True, restricted=True, wtune=True, run_nto=False,
                     wf_tag="", data_store="local"):
    """Build one workflow per molecule and add it to ``launchpad``. Returns the first fw_id of each."""
    paramset = GausParamSet.from_json(str(param_file))
    kwarg_dict = dict(paramset=paramset, public=public, submit=submit, data_store=data_store,
                      check_if_already_run=False, priority=priority,
                      gaussian_file_name="gaussian" if restricted else "gaussian_ur", restricted=restricted,
                      wtune=wtune, run_nto=run_nto, origin_group=origin_group, wf_tag=wf_tag)
    fw_ids = []
    for smiles, geom in load_molecules(molecule_file).items():
        wf = workflow_function(**dict(kwarg_dict, smiles=smiles, mol_opt_geom=geom))
        fw_ids.append(list(launchpad.add_wf(wf).values())[0])
    return fw_ids


class RecordingLaunchPad:
    """In-memory LaunchPad stand-in: records workflows and assigns sequential fw_ids like add_wf."""

    def __init__(self):
        self.workflows = []
        self._next_id = 1

    def add_wf(self, wf, *args, **kwargs):
        self.workflows.append(wf)
        old_new = {fw.fw_id: self._next_id + i for i, fw in enumerate(wf.fws)}
        self._next_id += len(wf.fws)
        return old_new


# ----------------------------------------------------------------- graph helpers used by assertions
def fws_by_name(wf):
    names = [fw.name for fw in wf.fws]
    assert len(names) == len(set(names)), "Firework names must be unique: {}".format(names)
    return {fw.name: fw for fw in wf.fws}


def parents_by_name(wf):
    id_name = {fw.fw_id: fw.name for fw in wf.fws}
    parents = {name: set() for name in id_name.values()}
    for parent_id, child_ids in wf.links.items():
        for child_id in child_ids:
            parents[id_name[child_id]].add(id_name[parent_id])
    return parents


def only_task(fw):
    assert len(fw.tasks) == 1, fw.name
    return fw.tasks[0]
