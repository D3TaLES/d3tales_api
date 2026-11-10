from fireworks import Firework

from d3tales_api.Workflows.Gaussian import (
    RunGaussianEnergy,
    RunGaussianOpt,
    RunGaussianTDDFT,
    RunWtuning as RunGaussianWtuning,
)
from d3tales_api.Workflows.Initialize import MoleculeInit
from d3tales_api.Workflows.Orca import (
    RunOrcaEnergy,
    RunOrcaOpt,
    RunOrcaTDDFT,
    RunWtuning as RunOrcaWtuning,
)
from d3tales_api.Workflows.pyscf import (
    RunPySCFEnergy,
    RunPySCFOpt,
    RunPySCFTDDFT,
    RunPySCFWtuning,
)
from d3tales_api.Workflows.utils import G16_CMD, PATH, PROC_VM_KEY, RUNFILE_LOG
from d3tales_api.Workflows.d3tales_utils import EmailStart, EmailEnd  # noqa: F401

# Copyright 2021-2022, University of Kentucky

DEFAULT_SOFTWARE = "gaussian"
SUPPORTED_SOFTWARE = ("gaussian", "orca", "pyscf")

_FIRETASKS = {
    "gaussian": {
        "opt": RunGaussianOpt,
        "energy": RunGaussianEnergy,
        "tddft": RunGaussianTDDFT,
        "wtuning": RunGaussianWtuning,
    },
    "orca": {
        "opt": RunOrcaOpt,
        "energy": RunOrcaEnergy,
        "tddft": RunOrcaTDDFT,
        "wtuning": RunOrcaWtuning,
    },
    "pyscf": {
        "opt": RunPySCFOpt,
        "energy": RunPySCFEnergy,
        "tddft": RunPySCFTDDFT,
        "wtuning": RunPySCFWtuning,
    },
}


def _normalize_software(software):
    software = (software or DEFAULT_SOFTWARE).lower()
    if software not in SUPPORTED_SOFTWARE:
        raise ValueError(
            "Unsupported software {!r}. Choose one of: {}".format(
                software, ", ".join(SUPPORTED_SOFTWARE)
            )
        )
    return software


def _calc_spec(software, priority=None):
    """Use the selected engine as the FireWorks category."""
    spec = {"_category": software}
    if priority is not None:
        spec["_priority"] = priority
    return spec


def _task_kwargs(software, **kwargs):
    """Add only executable/runtime kwargs relevant to the selected engine."""
    kwargs["software"] = software
    if software == "gaussian":
        kwargs.update(g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY)
    return kwargs


class InitializeMolecule(Firework):
    def __init__(self, name="initial_molecule_data", parents=None, priority=None,
                 name_tag="", software=DEFAULT_SOFTWARE, **kwargs):
        # software is accepted for a uniform Firework API, although MoleculeInit
        # itself is independent of the electronic-structure engine.
        software = _normalize_software(software)
        spec = {"_category": "processing"}
        if priority is not None:
            spec["_priority"] = priority
        t = [MoleculeInit(software=software, **kwargs)]
        super().__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class MolOpt(Firework):
    def __init__(self, smiles=None, name="opt_mol", parents=None, priority=None,
                 name_tag="", software=DEFAULT_SOFTWARE, **kwargs):
        software = _normalize_software(software)
        kwargs.pop("use_iop", None)
        task_cls = _FIRETASKS[software]["opt"]
        task_kwargs = _task_kwargs(
            software,
            smiles=smiles,
            name=name,
            name_tag=name_tag,
            path=PATH,
            runfile_log=RUNFILE_LOG,
            skip_freq=True,
            use_iop=False,
            **kwargs,
        )
        t = [task_cls(**task_kwargs)]
        super().__init__(t, parents=parents, spec=_calc_spec(software, priority),
                         name="{}{}".format(name_tag, name))


class WTuning(Firework):
    def __init__(self, smiles=None, name="wtuning", parents=None, priority=None,
                 name_tag="", software=DEFAULT_SOFTWARE, **kwargs):
        software = _normalize_software(software)
        task_cls = _FIRETASKS[software]["wtuning"]
        task_kwargs = _task_kwargs(
            software,
            smiles=smiles,
            path=PATH,
            runfile_log=RUNFILE_LOG,
            name=name,
            name_tag=name_tag,
            geometry="opt_mol",
            **kwargs,
        )
        t = [task_cls(**task_kwargs)]
        super().__init__(t, parents=parents, spec=_calc_spec(software, priority),
                         name="{}{}".format(name_tag, name))


class Optimization(Firework):
    def __init__(self, smiles=None, species="groundState", geometry="opt_mol",
                 solvent=None, parents=None, priority=None, name_tag="",
                 software=DEFAULT_SOFTWARE, **kwargs):
        software = _normalize_software(software)
        _type = "solv_" + solvent if solvent else "gas_phase"
        name = "solv_opt_" + species if solvent else "opt_" + species
        task_cls = _FIRETASKS[software]["opt"]
        task_kwargs = _task_kwargs(
            software,
            smiles=smiles,
            name=name,
            name_tag=name_tag,
            path=PATH,
            runfile_log=RUNFILE_LOG,
            geometry=geometry,
            type=_type,
            subtype=species,
            solvent=solvent,
            **kwargs,
        )
        t = [task_cls(**task_kwargs)]
        super().__init__(t, parents=parents, spec=_calc_spec(software, priority),
                         name="{}{}".format(name_tag, name))


class Energy(Firework):
    def __init__(self, smiles=None, species="groundState", geometry="groundState",
                 parents=None, priority=None, name_tag="", solvent=None,
                 software=DEFAULT_SOFTWARE, **kwargs):
        software = _normalize_software(software)
        abbrev_dict = {
            "groundState": "gs", "cation1": "c1", "anion1": "a1",
            "cation2": "c2", "anion2": "a2",
        }
        _type = "solv_" + solvent if solvent else "gas_phase"
        prefix = "solv_energy" if solvent else "energy"
        name = "{}_{}{}".format(
            prefix, abbrev_dict[geometry.split("_")[-1]], abbrev_dict[species]
        )
        task_cls = _FIRETASKS[software]["energy"]
        task_kwargs = _task_kwargs(
            software,
            smiles=smiles,
            name=name,
            name_tag=name_tag,
            path=PATH,
            runfile_log=RUNFILE_LOG,
            geometry=geometry,
            type=_type,
            subtype=species,
            solvent=solvent,
            **kwargs,
        )
        t = [task_cls(**task_kwargs)]
        super().__init__(t, parents=parents, spec=_calc_spec(software, priority),
                         name="{}{}".format(name_tag, name))


class TDDFT(Firework):
    def __init__(self, smiles=None, species="groundState", solvent=None,
                 parents=None, priority=None, name_tag="", solv_geom=False,
                 software=DEFAULT_SOFTWARE, **kwargs):
        software = _normalize_software(software)
        name = "tddft_{}".format(species)
        geom_prefix = "solv_opt_" if solv_geom else "opt_"

        paramset = kwargs.get("paramset")
        route_parameters = getattr(paramset, "route_parameters", {}) if paramset is not None else {}
        route_keys = list(map(str.lower, route_parameters.get("td", [])))
        tddft_prefix = (
            "singlet_" if "singlets" in route_keys
            else "triplet_" if "triplets" in route_keys
            else "singlet_triplet_" if "50-50" in route_keys
            else "tddft_"
        )

        task_cls = _FIRETASKS[software]["tddft"]
        task_kwargs = _task_kwargs(
            software,
            smiles=smiles,
            name=name,
            name_tag=name_tag,
            path=PATH,
            runfile_log=RUNFILE_LOG,
            type="tddft",
            geometry="{}{}".format(geom_prefix, species),
            prefix=tddft_prefix + species,
            solvent=solvent,
            **kwargs,
        )
        t = [task_cls(**task_kwargs)]
        super().__init__(t, parents=parents, spec=_calc_spec(software, priority),
                         name="{}tddft_{}".format(name_tag, species))
