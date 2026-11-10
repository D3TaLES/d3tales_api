"""The Firework factory classes in d3tales_api.Workflows.Fireworks, tested one at a time."""
import pytest
from fireworks import Firework

from d3tales_api.Workflows import Fireworks as fw_mod
from d3tales_api.Workflows.Fireworks import InitializeMolecule, MolOpt, WTuning, Optimization, Energy, TDDFT
from tests.helpers.workflow_helpers import only_task

pytestmark = [pytest.mark.unit, pytest.mark.gaussian]


# ----------------------------------------------------------------- shared behaviour
@pytest.fixture
def factories(paramset):
    """name -> zero-arg builder for each calculation factory."""
    return {
        "molopt": lambda **kw: MolOpt(paramset=paramset.opt_groundState, **kw),
        "wtuning": lambda **kw: WTuning(paramset=paramset.wtuning, **kw),
        "opt": lambda **kw: Optimization(paramset=paramset.opt_groundState, **kw),
        "energy": lambda **kw: Energy(paramset=paramset.energy_groundState, **kw),
        "tddft": lambda **kw: TDDFT(paramset=paramset.tddft_groundState, **kw),
    }


CALC_KINDS = ["molopt", "wtuning", "opt", "energy", "tddft"]
TASK_KIND = {"molopt": "opt", "wtuning": "wtuning", "opt": "opt", "energy": "energy", "tddft": "tddft"}


@pytest.mark.parametrize("kind", CALC_KINDS)
def test_default_software_is_gaussian(factories, kind):
    fw = factories[kind]()
    assert isinstance(fw, Firework)
    assert fw.spec == {"_category": "gaussian"}
    t = only_task(fw)
    assert type(t) is fw_mod._FIRETASKS["gaussian"][TASK_KIND[kind]]
    assert t["software"] == "gaussian"
    assert t["g16_cmd"] == ">>g16_cmd<<"
    assert t["proc_vm_key"] == ">>processing_vm_key<<"
    assert t["path"] == ">>path<<"
    assert t["runfile_log"] == ">>runfile_log<<"


@pytest.mark.parametrize("kind", CALC_KINDS)
@pytest.mark.parametrize("software", ["orca", "pyscf", "ORCA"])
def test_other_software_routes_task_and_category(factories, kind, software):
    fw = factories[kind](software=software)
    sw = software.lower()
    assert fw.spec["_category"] == sw
    t = only_task(fw)
    assert type(t) is fw_mod._FIRETASKS[sw][TASK_KIND[kind]]
    assert "g16_cmd" not in t and "proc_vm_key" not in t


@pytest.mark.parametrize("kind", CALC_KINDS + ["init"])
def test_unsupported_software_raises(factories, kind):
    build = factories.get(kind, InitializeMolecule)
    with pytest.raises(ValueError, match="Unsupported software"):
        build(software="qchem")


@pytest.mark.parametrize("kind", CALC_KINDS)
def test_priority_and_name_tag(factories, kind):
    fw = factories[kind](priority=3, name_tag="t1_")
    assert fw.spec["_priority"] == 3
    assert fw.name.startswith("t1_")
    assert only_task(fw)["name_tag"] == "t1_"


@pytest.mark.parametrize("kind", CALC_KINDS)
def test_extra_kwargs_reach_task(factories, kind):
    t = only_task(factories[kind](data_store="local", restricted=False, custom_flag=42))
    assert (t["data_store"], t["restricted"], t["custom_flag"]) == ("local", False, 42)


# ----------------------------------------------------------------- InitializeMolecule
def test_initialize_molecule():
    fw = InitializeMolecule(smiles="O", identifier=None, mol_opt_geom=[{"x": 1}], priority=2, origin_group="Risko")
    assert fw.name == "initial_molecule_data"
    assert fw.spec == {"_category": "processing", "_priority": 2}
    t = only_task(fw)
    assert type(t).__name__ == "MoleculeInit"
    assert (t["smiles"], t["mol_opt_geom"], t["origin_group"], t["software"]) == ("O", [{"x": 1}], "Risko", "gaussian")


# ----------------------------------------------------------------- per-factory naming / geometry rules
def test_molopt_disables_iop_and_freq(paramset):
    t = only_task(MolOpt(paramset=paramset.opt_groundState, use_iop=True))
    assert t["use_iop"] is False and t["skip_freq"] is True
    assert t["name"] == "opt_mol"


def test_wtuning_reads_initial_geometry(paramset):
    fw = WTuning(paramset=paramset.wtuning)
    assert fw.name == "wtuning"
    assert only_task(fw)["geometry"] == "opt_mol"


@pytest.mark.parametrize("species, solvent, name, type_", [
    ("groundState", None, "opt_groundState", "gas_phase"),
    ("anion1", None, "opt_anion1", "gas_phase"),
    ("cation1", "acetonitrile", "solv_opt_cation1", "solv_acetonitrile"),
])
def test_optimization_naming(paramset, species, solvent, name, type_):
    fw = Optimization(paramset=paramset.opt_groundState, species=species, solvent=solvent)
    t = only_task(fw)
    assert fw.name == t["name"] == name
    assert (t["type"], t["subtype"], t["solvent"], t["geometry"]) == (type_, species, solvent, "opt_mol")


@pytest.mark.parametrize("species, geometry, solvent, name", [
    ("groundState", "opt_anion1", None, "energy_a1gs"),
    ("anion1", "opt_groundState", None, "energy_gsa1"),
    ("groundState", "opt_cation1", None, "energy_c1gs"),
    ("cation1", "opt_groundState", None, "energy_gsc1"),
    ("cation2", "opt_groundState", None, "energy_gsc2"),
    ("anion1", "opt_anion1", "acetonitrile", "solv_energy_a1a1"),
])
def test_energy_naming(paramset, species, geometry, solvent, name):
    fw = Energy(paramset=paramset.energy_groundState, species=species, geometry=geometry, solvent=solvent)
    t = only_task(fw)
    assert fw.name == t["name"] == name
    assert (t["geometry"], t["subtype"]) == (geometry, species)
    assert t["type"] == ("solv_" + solvent if solvent else "gas_phase")


def test_energy_unknown_species_raises(paramset):
    with pytest.raises(KeyError):
        Energy(paramset=paramset.energy_groundState, species="triplet", geometry="opt_groundState")


@pytest.mark.parametrize("td, prefix", [
    ({"50-50": "", "nstates": 15}, "singlet_triplet_"),
    ({"singlets": "", "nstates": 10}, "singlet_"),
    ({"Triplets": ""}, "triplet_"),
    ({"nstates": 10}, "tddft_"),
])
def test_tddft_prefix_from_route_parameters(td, prefix):
    from d3tales_api.Workflows.ParamSet import GaussianParameters
    ps = GaussianParameters(route_parameters={"td": td})
    assert only_task(TDDFT(paramset=ps, species="anion1"))["prefix"] == prefix + "anion1"


def test_tddft_without_paramset_uses_generic_prefix():
    assert only_task(TDDFT(species="cation1"))["prefix"] == "tddft_cation1"


@pytest.mark.parametrize("solv_geom, geometry", [(False, "opt_cation1"), (True, "solv_opt_cation1")])
def test_tddft_geometry(paramset, solv_geom, geometry):
    fw = TDDFT(paramset=paramset.tddft_cation1, species="cation1", solv_geom=solv_geom)
    t = only_task(fw)
    assert fw.name == t["name"] == "tddft_cation1"
    assert (t["geometry"], t["type"]) == (geometry, "tddft")


def test_parents_are_linked(paramset):
    parent = InitializeMolecule(smiles="O")
    child = Optimization(paramset=paramset.opt_groundState, parents=parent)
    assert child.parents == [parent]
