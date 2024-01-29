from fireworks import Firework
from d3tales_api.Workflows.Gaussian import *
from d3tales_api.Workflows.Initialize import *
from d3tales_api.Workflows.utils import G16_CMD, PATH, PROC_VM_KEY, RUNFILE_LOG

# Copyright 2021-2022, University of Kentucky

CALC_CATEGORY = 'gaussian'  # gaussian gaussian_mcc


class InitializeMolecule(Firework):
    def __init__(self, name="initial_molecule_data", parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': 'processing', '_priority': priority} if priority else {'_category': 'processing'}
        t = [MoleculeInit(**kwargs)]
        super(InitializeMolecule, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class MolOpt(Firework):
    def __init__(self, name="opt_mol", parents=None, priority=None, name_tag='', **kwargs):
        kwargs.pop("use_iop") if "use_iop" in kwargs.keys() else None
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        t = [RunGaussianOpt(name=name, name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY,
                            path=PATH, runfile_log=RUNFILE_LOG, skip_freq=True, use_iop=False, **kwargs)]
        super(MolOpt, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class WTuning(Firework):
    def __init__(self, name="wtuning", parents=None,  priority=None, name_tag='', **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        t = [RunWtuning(path=PATH, runfile_log=RUNFILE_LOG, name=name, name_tag=name_tag, geometry="opt_mol", **kwargs)]
        super(WTuning, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class Optimization(Firework):
    def __init__(self, species="groundState", geometry="opt_mol", solvent=None, parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        _type = "solv_" + solvent if solvent else "gas_phase"
        name = "solv_opt_" + species if solvent else "opt_" + species
        t = [RunGaussianOpt(name=name, name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH, runfile_log=RUNFILE_LOG,
                            geometry=geometry, type=_type, subtype=species, solvent=solvent, **kwargs)]
        super(Optimization, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class Energy(Firework):
    def __init__(self, species="groundState", geometry="groundState", parents=None, priority=None, name_tag='', solvent=None, **kwargs):
        abbrev_dict = {"groundState": "gs", "cation1": "c1", "anion1": "a1", "cation2": "c2", "anion2": "a2", }
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        _type = "solv_" + solvent if solvent else "gas_phase"
        prefix = "solv_energy" if solvent else "energy"
        name = "{}_{}{}".format(prefix, abbrev_dict[geometry.split("_")[-1]], abbrev_dict[species])
        t = [RunGaussianEnergy(name=name, name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH, runfile_log=RUNFILE_LOG,
                               geometry=geometry, type=_type, subtype=species, solvent=solvent, **kwargs)]
        super(Energy, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class TDDFT(Firework):
    def __init__(self, species='groundState', solvent=None, parents=None, priority=None, name_tag='', solv_geom=False, **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        name = "tddft_{}".format(species)
        geom_prefix = "solv_opt_" if solv_geom else "opt_"
        route_keys = list(map(str.lower, kwargs.get("paramset").route_parameters.get("td", [])))
        tddft_prefix = "singlet_" if "singlets" in route_keys else "triplet_" if "triplets" in route_keys else "singlet_triplet_" if "50-50" in route_keys else "tddft_"
        t = [RunGaussianTDDFT(name=name, name_tag=name_tag, path=PATH, runfile_log=RUNFILE_LOG, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, type='tddft',
                              geometry="{}{}".format(geom_prefix, species), prefix=tddft_prefix+species, solvent=solvent, **kwargs)]
        super(TDDFT, self).__init__(t, parents=parents, spec=spec, name="{}tddft_{}".format(name_tag, species))


# RETIRED CLASS. This class needs the OCELOT API to be installed.

# class LowestEConformer(Firework):
#     def __init__(self, parents=None, priority=None, name_tag='', **kwargs):
#         spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
#         t = [GetLowestEConformer(name="conformers", name_tag=name_tag, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH,
#                                  runfile_log=RUNFILE_LOG, type="gas_phase", **kwargs)]
#         super(LowestEConformer, self).__init__(t, parents=parents, spec=spec, name="{}conformers".format(name_tag))


class DihedRot(Firework):
    def __init__(self, geometry="groundState", dihed_degree=0, parents=None, priority=None, name_tag='', **kwargs):
        spec = {'_category': CALC_CATEGORY, '_priority': priority} if priority else {'_category': CALC_CATEGORY}
        name = "dihedrot_" + str(dihed_degree).zfill(3)
        t = [RunGaussianDihedRot(name=name, name_tag=name_tag, dihed_degree=dihed_degree, g16_cmd=G16_CMD, proc_vm_key=PROC_VM_KEY, path=PATH, runfile_log=RUNFILE_LOG,
                                 geometry=geometry, type="gas_phase", subtype=str(dihed_degree).zfill(3), **kwargs)]
        super(DihedRot, self).__init__(t, parents=parents, spec=spec, name="{}{}".format(name_tag, name))


class EmailStart(Firework):
    def __init__(self, parents=None, priority=None, name_tag='', identifier="", email="", username="", **kwargs):
        spec = {'_category': 'processing', '_priority': priority} if priority else {'_category': 'processing'}
        t = [EmailStarting(identifier=identifier, email=email, username=username)]
        super(EmailStart, self).__init__(t, parents=parents, spec=spec, name="{}email_starting".format(name_tag))


class EmailEnd(Firework):
    def __init__(self, parents=None, priority=None, name_tag='', identifier="", email="", username="", **kwargs):
        spec = {'_category': 'processing', '_priority': priority} if priority else {'_category': 'processing'}
        t = [EmailFinished(identifier=identifier, email=email, username=username)]
        super(EmailEnd, self).__init__(t, parents=parents, spec=spec, name="{}email_finished".format(name_tag))


