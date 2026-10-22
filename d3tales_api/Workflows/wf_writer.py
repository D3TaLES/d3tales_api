from fireworks import Workflow
from d3tales_api.Workflows.D3TaLES_FW import *


# Copyright 2021, University of Kentucky


def d3tales_wf(paramset, identifier=None, smiles=None, wtune=True, solvent='acetonitrile',
               hf_mol_opt=False, email=None, username=None, wf_tag="", **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    opt_parents = f30 if wtune else f10
    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=opt_parents, **kwargs)
    f35 = Optimization(paramset=paramset.opt_anion1, species="anion1", parents=opt_parents, **kwargs)
    f36 = Optimization(paramset=paramset.opt_anion2, species="anion2", parents=opt_parents, **kwargs)
    f37 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=opt_parents, **kwargs)
    f38 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=opt_parents, **kwargs)

    f45 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="opt_anion1", parents=[f35],
                 **kwargs)
    f46 = Energy(paramset=paramset.energy_anion1, species="anion1", geometry="opt_groundState", parents=[f31], **kwargs)
    f47 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="opt_cation1", parents=[f37],
                 **kwargs)
    f48 = Energy(paramset=paramset.energy_cation1, species="cation1", geometry="opt_groundState", parents=[f31],
                 **kwargs)

    f51 = Energy(paramset=paramset.energy_groundState, parents=[f31], species='groundState', geometry="opt_groundState",
                 solvent=solvent, **kwargs)
    f55 = Energy(paramset=paramset.energy_anion1, parents=[f35], species='anion1', geometry="opt_anion1",
                 solvent=solvent, **kwargs)
    f56 = Energy(paramset=paramset.energy_anion2, parents=[f36], species='anion2', geometry="opt_anion2",
                 solvent=solvent, **kwargs)
    f57 = Energy(paramset=paramset.energy_cation1, parents=[f37], species='cation1', geometry="opt_cation1",
                 solvent=solvent, **kwargs)
    f58 = Energy(paramset=paramset.energy_cation2, parents=[f38], species='cation2', geometry="opt_cation2",
                 solvent=solvent, **kwargs)

    f61 = TDDFT(paramset=paramset.tddft_groundState, parents=[f31], **kwargs)
    f65 = TDDFT(paramset=paramset.tddft_anion1, species='anion1', parents=[f35], **kwargs)
    f66 = TDDFT(paramset=paramset.tddft_anion2, species='anion2', parents=[f36], **kwargs)
    f67 = TDDFT(paramset=paramset.tddft_cation1, species='cation1', parents=[f37], **kwargs)
    f68 = TDDFT(paramset=paramset.tddft_cation2, species='cation2', parents=[f38], **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f31, f35, f36, f37, f38, f45, f46, f47, f48, f51, f55, f56, f57, f58, f61, f65, f66, f67, f68, ]
    if wtune:
        fws.extend([f26, f30])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(
            EmailEnd(identifier=identifier, email=email, username=username, parents=[f51, f55, f57, f61, f65, f67]))

    wf = Workflow(fws, name="{}gaus_{}".format(wf_tag, identifier or smiles))

    return wf


def just_anion(paramset, identifier=None, smiles=None, wtune=False, solvent='acetonitrile',
               hf_mol_opt=False, email=None, username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    opt_parents = f30 if wtune else f10

    f35 = Optimization(paramset=paramset.opt_anion1, species="anion1", parents=opt_parents, **kwargs)
    f36 = Optimization(paramset=paramset.opt_anion2, species="anion2", parents=opt_parents, **kwargs)

    f45 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="opt_anion1", parents=[f35],
                 **kwargs)
    f46 = Energy(paramset=paramset.energy_anion1, species="anion1", geometry="opt_groundState", parents=[opt_parents],
                 **kwargs)

    f55 = Energy(paramset=paramset.energy_anion1, parents=[f35], species='anion1', geometry="opt_anion1",
                 solvent=solvent, **kwargs)
    f56 = Energy(paramset=paramset.energy_anion2, parents=[f36], species='anion2', geometry="opt_anion2",
                 solvent=solvent, **kwargs)

    f65 = TDDFT(paramset=paramset.tddft_anion1, species='anion1', parents=[f35], **kwargs)
    f66 = TDDFT(paramset=paramset.tddft_anion2, species='anion2', parents=[f36], **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f35, f36, f45, f46, f55, f56, f65, f66]

    if wtune:
        fws.extend([f26, f30])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(EmailEnd(identifier=identifier, email=email, username=username, parents=[f55, f65, f66]))

    wf = Workflow(fws, name="anion_{}".format(identifier or smiles))

    return wf


def solv_wf(paramset, identifier=None, smiles=None, wtune=True, solvent='acetonitrile',
            hf_mol_opt=False, email=None, username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles,
                             **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)

    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    opt_parents = f30 if wtune else f10
    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=opt_parents, solvent=solvent,
                       **kwargs)
    f37 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=opt_parents, solvent=solvent, **kwargs)
    f38 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=opt_parents, solvent=solvent, **kwargs)

    f47 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="solv_opt_cation1",
                 parents=[f37], solvent=solvent, **kwargs)
    f48 = Energy(paramset=paramset.energy_cation1, species="cation1", geometry="solv_opt_groundState", parents=[f31],
                 solvent=solvent, **kwargs)

    f61 = TDDFT(paramset=paramset.tddft_groundState, parents=[f31], solvent=solvent, solv_geom=True, **kwargs)
    f67 = TDDFT(paramset=paramset.tddft_cation1, species='cation1', parents=[f37], solvent=solvent, solv_geom=True,
                **kwargs)
    f68 = TDDFT(paramset=paramset.tddft_cation2, species='cation2', parents=[f38], solvent=solvent, solv_geom=True,
                **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f31, f37, f38, f47, f48, f61, f67, f68, ]
    if wtune:
        fws.extend([f26, f30])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(EmailEnd(identifier=identifier, email=email, username=username, parents=[f61, f67]))

    wf = Workflow(fws, name="solv_{}".format(identifier or smiles))

    return wf


def hf_wf(paramset, identifier=None, smiles=None, hf_mol_opt=False, email=None, username=None, solvent=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    f36 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=f30, **kwargs)
    f37 = Optimization(paramset=paramset.hf_opt_groundState, species="groundState", parents=f30, name_tag='hf_',
                       **kwargs)
    f46 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=f30, **kwargs)
    f47 = Optimization(paramset=paramset.hf_opt_cation1, species="cation1", parents=f30, name_tag='hf_', **kwargs)
    f56 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=f30, **kwargs)
    f57 = Optimization(paramset=paramset.hf_opt_cation2, species="cation2", parents=f30, name_tag='hf_', **kwargs)

    f38 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="hf_opt_groundState",
                 parents=[f37], name_tag='hfgeom_', **kwargs)
    f39 = Energy(paramset=paramset.hfgeom_energy_groundState, species="groundState", geometry="opt_groundState",
                 parents=[f36], name_tag='hf_dftgeom_', **kwargs)
    f48 = Energy(paramset=paramset.energy_cation1, species="cation1", geometry="hf_opt_cation1", parents=[f47],
                 name_tag='hfgeom_', **kwargs)
    f49 = Energy(paramset=paramset.hfgeom_energy_cation1, species="cation1", geometry="opt_cation1", parents=[f46],
                 name_tag='hf_dftgeom_', **kwargs)
    f58 = Energy(paramset=paramset.energy_cation2, species="cation2", geometry="hf_opt_cation2", parents=[f57],
                 name_tag='hfgeom_', **kwargs)
    f59 = Energy(paramset=paramset.hfgeom_energy_cation2, species="cation2", geometry="opt_cation2", parents=[f56],
                 name_tag='hf_dftgeom_', **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f26, f30, f36, f37, f38, f39, f46, f47, f48, f49, f56, f57, f58, f59]
    wf = Workflow(fws, name="gaus_hf_{}".format(identifier or smiles))

    return wf


def huckaba_wf(paramset, identifier=None, smiles=None, solvent='DiMethylSulfoxide', wtune=True, hf_mol_opt=False,
               email=None, username=None, **kwargs):
    kwargs.update(dict(iop_str="0170700000"))
    # name = kwargs.get('mol_name')
    # species = "groundState" if "_p0" in name else "cation1" if "_p1" in name else "cation2"
    species = "groundState"

    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    f26 = MolOpt(paramset=paramset.opt_groundState, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)
    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=f30, **kwargs)

    opt_parents = f30 if wtune else f10

    f37 = Optimization(paramset=eval("paramset.opt_" + species), species=species, parents=opt_parents, solvent=solvent,
                       **kwargs)
    f63 = TDDFT(paramset=eval("paramset.tddft_" + species), parents=[f37], **kwargs)

    fws = [f10, f37, f63]
    if wtune:
        fws.extend([f26, f30, f31])

    wf = Workflow(fws, name="special_{}".format(identifier or smiles))
    return wf


def just_nmr(paramset, identifier=None, smiles=None, solvent='acetonitrile', hf_mol_opt=False, email=None,
             username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    f76 = Energy(paramset=paramset.nmr_groundState, parents=[f10], species='groundState', geometry="opt_groundState",
                 name_tag='nmr_', solvent=solvent, **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f76, ]
    wf = Workflow(fws, name="gaus_nmr_{}".format(identifier or smiles))

    return wf


def just_tddft(paramset, identifier=None, smiles=None, wtune=False, solvent=None,
               hf_mol_opt=False, email=None, username=None, check_if_already_run=True, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState

    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)
    opt_parents = f30 if wtune else f10

    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=opt_parents, solvent=solvent,
                       check_if_already_run=check_if_already_run, **kwargs)
    f37 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=opt_parents, solvent=solvent,
                       check_if_already_run=check_if_already_run, **kwargs)
    f38 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=opt_parents, solvent=solvent,
                       check_if_already_run=check_if_already_run, **kwargs)

    f61 = TDDFT(paramset=paramset.tddft_groundState, parents=[f31], solvent=solvent, solv_geom=True, **kwargs)
    f67 = TDDFT(paramset=paramset.tddft_cation1, species='cation1', parents=[f37], solvent=solvent, solv_geom=True,
                **kwargs)
    f68 = TDDFT(paramset=paramset.tddft_cation2, species='cation2', parents=[f38], solvent=solvent, solv_geom=True,
                **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f31, f37, f38, f61, f67, f68, ]
    if wtune:
        fws.extend([f26, f30])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(EmailEnd(identifier=identifier, email=email, username=username, parents=[f61, f67]))

    wf = Workflow(fws, name="tddft_{}".format(identifier or smiles))

    return wf


def just_initialize(paramset, identifier=None, smiles=None, solvent='acetonitrile', hf_mol_opt=False, email=None,
                    username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    fws = [f10, ]
    wf = Workflow(fws, name="init_{}".format(identifier or smiles))
    return wf


# RETIRED WORKFLOW. This class needs the OCELOT API to be installed.
# def dihed_rot(paramset, identifier=None, smiles=None, solvent='acetonitrile', hf_mol_opt=False, email=None,
#               username=None, **kwargs):
#     f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
#     # f12 = LowestEConformer(smiles=smiles, paramset=paramset.semi_empirical, parents=f10, **kwargs)
#     f26 = Optimization(paramset=paramset.opt_groundState, geometry="conformer", species='groundState', parents=f12,
#                        **kwargs)
#     f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)
#     f36 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=f30, **kwargs)
#
#     # fws = [f10, f36, ]  # exclude omega tuning
#     fws = [f10, f12, f26, f30, f36, ]  # include omega tuning
#
#     for degree in np.linspace(start=0, stop=180, num=19):
#         fw = DihedRot(paramset=paramset.opt_groundState, dihed_degree=degree, parents=f36, **kwargs)
#         fws.append(fw)
#
#     wf = Workflow(fws, name="dihedrot_{}".format(identifier or smiles))
#     return wf


def specialized_job(paramset, identifier=None, smiles=None, solvent='acetonitrile', wtune=True, hf_mol_opt=False,
                    email=None, username=None, check_if_already_run=False, **kwargs):
    kwargs.update(dict(iop_str="0170700000"))
    # name = kwargs.get('mol_name')
    # species = "groundState" if "_p0" in name else "cation1" if "_p1" in name else "cation2"
    species = "groundState"

    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    f26 = MolOpt(paramset=paramset.opt_groundState, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    opt_parents = f30 if wtune else f10

    f36 = Optimization(paramset=eval("paramset.opt_" + species), species=species, parents=opt_parents, **kwargs)
    f37 = Optimization(paramset=eval("paramset.opt_" + species), species=species, parents=opt_parents, solvent=solvent,
                       **kwargs)
    f63 = TDDFT(paramset=eval("paramset.tddft_" + species), parents=[f36], **kwargs)

    fws = [f10, f36, f37, f63]
    if wtune:
        fws.extend([f26, f30])

    wf = Workflow(fws, name="special_{}".format(identifier or smiles))
    return wf
