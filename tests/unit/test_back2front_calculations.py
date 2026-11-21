import pytest

from d3tales_api.D3database.back2front import (
    CV2Front,
    DOI2Front,
    Gaus2FrontCharacterization,
)


def _gaus_characterization():
    obj = Gaus2FrontCharacterization.__new__(Gaus2FrontCharacterization)
    obj.id = "M1"
    obj.calculation_type = None
    obj.schema_version = "v2"
    obj.bare_conditions = {"functional": "B3LYP", "basis_set": "6-31G"}
    obj.solvent = {"name": "Acetonitrile", "model": "implicit_solvent", "dielectric_constant": 35.688}
    obj.groundState_charge = 0
    obj.species_dict = {"anion2": -2, "anion1": -1, "groundState": 0, "cation1": 1, "cation2": 2}
    obj.verbose = 0
    obj.all_props = True
    data = {}
    values = {
        "opt_groundState": -10.0,
        "opt_cation1": -8.0,
        "opt_cation2": -5.0,
        "opt_anion1": -12.0,
        "opt_anion2": -13.0,
        "energy_gsc1": -9.5,
        "energy_c1gs": -9.0,
        "energy_gsa1": -10.5,
        "energy_a1gs": -11.0,
        "energy_c1c2": -7.0,
        "energy_a1a2": -11.5,
        "solv_energy_gsgs": -10.2,
        "solv_energy_c1c1": -8.4,
        "solv_energy_a1a1": -12.1,
    }
    for key, energy in values.items():
        data[key] = {"scf_total_energy": {"value": energy, "unit": "eV"}}
    for calc in ("freq_groundState", "freq_cation1", "freq_anion1"):
        data[calc] = {"gibbs_correction": {"value": 0.1, "unit": "eV"}}
    for calc, calc_data in data.items():
        setattr(obj, calc, (calc_data, f"hash-{calc}"))
    return obj


@pytest.mark.unit
@pytest.mark.parametrize(
    "method, name, expected",
    [
        ("hole_reorganization_energy", "hole_reorganization_energy", -0.5),
        ("electron_reorganization_energy", "electron_reorganization_energy", 0.5),
        ("relaxation_groundState_cation1", "relaxation_groundState_cation1", -1.5),
        ("relaxation_cation1_groundState", "relaxation_cation1_groundState", 1.0),
        ("relaxation_groundState_anion1", "relaxation_groundState_anion1", 1.5),
        ("relaxation_anion1_groundState", "relaxation_anion1_groundState", -1.0),
        ("vertical_ionization_energy", "vertical_ionization_energy", 0.5),
        ("vertical_ionization_energy_2", "vertical_ionization_energy_2", 1.0),
        ("vertical_electron_affinity", "vertical_electron_affinity", -0.5),
        ("vertical_electron_affinity_2", "vertical_electron_affinity_2", 0.5),
        ("adiabatic_ionization_energy", "adiabatic_ionization_energy", 2.0),
        ("adiabatic_ionization_energy_2", "adiabatic_ionization_energy_2", 3.0),
        ("adiabatic_electron_affinity", "adiabatic_electron_affinity", -2.0),
        ("adiabatic_electron_affinity_2", "adiabatic_electron_affinity_2", -1.0),
    ],
)
def test_gaussian_descriptor_calculations(method, name, expected):
    obj = _gaus_characterization()
    result = getattr(obj, method)()
    assert result[name]["value"] == pytest.approx(expected)
    assert result[name]["conditions"] == obj.bare_conditions
    assert result[name]["source_hash_ids"]


@pytest.mark.unit
def test_gaussian_solvation_energy_adds_solvent_condition():
    obj = _gaus_characterization()
    result = obj.solvation_energy("groundState")["solvation_energy"]
    assert result["value"] == pytest.approx(-0.2)
    assert result["conditions"]["solvent"] == obj.solvent


@pytest.mark.unit
def test_gaussian_omega_uses_truncated_tuning_value():
    obj = _gaus_characterization()
    obj.find_data = lambda calc, conditions=None: ({"omega": 0.123456}, "omega-hash")
    obj.bare_conditions["tuning_parameter"] = 0.2
    value = obj.omega()["omega"]
    assert value["value"] == pytest.approx(0.1234)
    assert "tuning_parameter" not in value["conditions"]


@pytest.mark.unit
def test_gaussian_get_data_rejects_unrelated_calculation_type():
    obj = _gaus_characterization()
    obj.calculation_type = "opt_groundState"
    obj.all_props = False
    with pytest.raises(ValueError, match="does not use"):
        obj.get_data(["opt_cation1"])


@pytest.mark.unit
def test_gaussian_get_data_reports_missing_required_calculations():
    obj = _gaus_characterization()
    obj.calculation_type = "opt_groundState"
    obj.all_props = False
    obj.opt_cation1 = (None, "missing")
    with pytest.raises(ValueError, match="missing"):
        obj.get_data(["opt_groundState", "opt_cation1"])


@pytest.mark.unit
def test_gaussian_constructor_loads_all_expected_calculation_slots(monkeypatch):
    class Collection:
        def find_one(self, query):
            return {"mol_info": {"groundState_charge": 0}}

        def find(self, query):
            return [{"_id": "existing", "data": {"scf_total_energy": {"value": -1}}}]

    collection = Collection()

    class Connector:
        def __init__(self, info):
            pass

        def get_collection(self, name):
            return collection

    monkeypatch.setattr("d3tales_api.D3database.back2front.DBconnector", Connector)
    monkeypatch.setattr("d3tales_api.D3database.back2front.db_info", lambda: {"d3tales": {}})
    obj = Gaus2FrontCharacterization(
        "M1",
        conditions={"functional": "B3LYP", "basis_set": "6-31G"},
        verbose=0,
    )
    assert obj.groundState_charge == 0
    assert len(obj.all_calcs) == 26
    assert len(obj.mol_properties) == 17
    assert obj.opt_groundState[0] is None


@pytest.mark.unit
def test_gaussian_oxidation_and_reduction_potentials_use_solvent_conditions():
    obj = _gaus_characterization()
    oxidation = obj.oxidation_potential()["oxidation_potential"]
    reduction = obj.reduction_potential()["reduction_potential"]
    assert oxidation["conditions"]["reference_electrode"] == "standard_hydrogen_electrode"
    assert reduction["conditions"]["solvent"] == obj.solvent
    assert oxidation["unit"] == reduction["unit"] == "eV"


@pytest.mark.unit
def test_gaussian_insert_methods_dispatch_frontend_records(monkeypatch):
    insertions = []
    monkeypatch.setattr(
        Gaus2FrontCharacterization,
        "_insert_frontend",
        lambda self, schema_layer, instance: insertions.append((schema_layer, instance)),
    )
    obj = _gaus_characterization()
    obj.calculation_type = None
    obj.verbose = 0
    assert obj.insert() is None

    obj.calculation_type = "opt_groundState"
    obj.species_descriptors = lambda calc, **kwargs: {"groundState": {"homo": {"value": -5}}}
    obj.mol_properties = [lambda: {"energy": {"value": -1}}]
    obj.all_calcs = ["opt_groundState"]
    obj.insert()
    obj.insert_all_species()
    obj.insert_all_mol()
    assert [entry[0] for entry in insertions] == [
        "species_characterization",
        "mol_characterization",
        "species_characterization",
        "mol_characterization",
    ]


@pytest.mark.unit
def test_gaussian_frontend_insertion_uses_unified_database(monkeypatch):
    calls = []

    class Database:
        def __init__(self, **kwargs):
            calls.append(("init", kwargs))

        def insert(self, identifier, nested):
            calls.append(("insert", identifier, nested))

    monkeypatch.setattr("d3tales_api.D3database.back2front.D3Database", Database)
    obj = _gaus_characterization()
    obj.id = "M1"
    obj.verbose = 0
    obj._insert_frontend("species_characterization", {"groundState": {"charge": 0}})

    kwargs = calls[0][1]
    assert kwargs["database"] == "d3tales"
    assert kwargs["collection_name"] == "base"
    assert kwargs["schema_db"] == "frontend"
    assert kwargs["schema_layer"] == "species_characterization"
    assert kwargs["default_id"] == "M1"
    assert calls[1] == ("insert", "M1", True)


@pytest.mark.unit
def test_doi2front_construction_and_molecule_lookup_warning(monkeypatch):
    collection = object()

    class Connector:
        def __init__(self, info):
            pass

        def get_collection(self, name):
            return collection

    monkeypatch.setattr("d3tales_api.D3database.back2front.DBconnector", Connector)
    monkeypatch.setattr("d3tales_api.D3database.back2front.db_info", lambda: {"d3tales": {}})
    original_get_mol_id = DOI2Front.get_mol_id
    monkeypatch.setattr(DOI2Front, "get_mol_id", lambda self, mol_data, nlp_group: "M1")
    backend = {
        "_id": "10/example",
        "extracted_molecules": [{
            "molecule_name": "sample",
            "extracted_properties": {"solubility": [{"conditions": {}, "value": 1}]},
        }],
    }
    doi = DOI2Front(backend_data=backend, insert=False)
    assert doi.mol_ids == {"M1"}
    assert doi.extracted_mol_data["M1"]["literature_data"]["extracted_properties"]["solubility"][0]["conditions"]["doi"] == "10/example"

    monkeypatch.setattr(DOI2Front, "get_mol_id", original_get_mol_id)
    monkeypatch.setattr("d3tales_api.D3database.back2front.pcp.get_compounds", lambda name, kind: [])
    doi.cir_convert = lambda name: None
    with pytest.warns(UserWarning, match="No PubChem or CIR"):
        assert doi.get_mol_id({"molecule_name": "unlisted"}, "Risko") is None

    compounds = [type("Compound", (), {"iupac_name": "A"})(), type("Compound", (), {"iupac_name": "B"})()]
    monkeypatch.setattr("d3tales_api.D3database.back2front.pcp.get_compounds", lambda name, kind: compounds)
    with pytest.warns(UserWarning, match="Multiple PubChem"):
        assert doi.get_mol_id({"molecule_name": "ambiguous"}, "Risko") is None


@pytest.mark.unit
def test_species_descriptors_for_opt_calculation_without_geometry():
    obj = _gaus_characterization()
    obj.opt_groundState = (
        {
            "conditions": {"functional": "B3LYP"},
            "charge": 0,
            "spin_multiplicity": 1,
            "is_groundState": True,
            "homo": {"value": -5.0, "unit": "eV"},
            "lumo": {"value": -1.0, "unit": "eV"},
            "scf_dipole_moment": {"value": [1.0, 2.0, 2.0], "unit": "Debye"},
        },
        "hash-opt",
    )
    result = obj.species_descriptors("opt_groundState", skip_geom=True)["groundState"]
    assert result["homo_lumo_gap"]["value"] == pytest.approx(4.0)
    assert result["dipole_moment"]["value"] == pytest.approx(3.0)
    assert "geometry" not in result


@pytest.mark.unit
def test_species_descriptors_for_tddft():
    obj = _gaus_characterization()
    obj.tddft_groundState = (
        {
            "conditions": {},
            "charge": 0,
            "spin_multiplicity": 1,
            "is_groundState": True,
            "excitations": {"Singlet": [(2.0, 600, 0.2)], "Triplet": [(1.0, 1200, 0.1)]},
            "singlet_plotting": {"energies": [1.0]},
        },
        "hash-td",
    )
    result = obj.species_descriptors("tddft_groundState")["groundState"]
    assert result["singlet_states"]["excitations"] == [(2.0, 600, 0.2)]
    assert result["triplet_states"]["excitations"] == [(1.0, 1200, 0.1)]
    assert "spectra" in result


@pytest.mark.unit
def test_doi2front_adds_doi_to_extracted_conditions():
    obj = DOI2Front.__new__(DOI2Front)
    obj.doi = "10.0000/example"
    props = {"energy": [{"value": 1.2, "conditions": {}}, {"value": 1.3, "conditions": {"doi": "existing"}}]}
    result = obj.check_extr_props(props)
    assert result["energy"][0]["conditions"]["doi"] == obj.doi
    assert result["energy"][1]["conditions"]["doi"] == "existing"
    assert obj.get_literature_data({"extracted_properties": props})["literature_data"]["related_literature"] == [obj.doi]


@pytest.mark.unit
def test_cv2front_shared_conditions_and_e_half_entries():
    data = [{
        "_id": "cv-1",
        "mol_id": "M1",
        "data": {
            "conditions": {"scan_rate": {"value": 0.1, "unit": "V/s"}, "temperature": {"value": 298, "unit": "K"}},
            "e_half": [{"value": 0.25, "unit": "V"}, {"value": 0.5, "unit": "V"}],
        },
    }]
    cv = CV2Front(backend_data=data, run_processing=False, insert=False, verbose=0, e_half_scan_rate=0.1)
    assert cv.conditions == data[0]["data"]["conditions"]
    assert [entry["value"] for entry in cv.e_halfs] == [0.25, 0.5]
    assert cv.prop_dict(1, "eV", notes="test")["notes"] == "test"


@pytest.mark.unit
def test_cv2front_shared_conditions_excludes_varying_fields():
    base = {"scan_rate": {"value": 0.1}, "temperature": {"value": 298}}
    data = [
        {"_id": "cv-1", "mol_id": "M1", "data": {"conditions": base}},
        {"_id": "cv-2", "mol_id": "M1", "data": {"conditions": {"scan_rate": {"value": 0.2}, "temperature": {"value": 298}}}},
    ]
    cv = CV2Front(backend_data=data, run_processing=False, insert=False, verbose=0, e_half_scan_rate=0.1)
    assert cv.conditions == {"temperature": {"value": 298}}


@pytest.mark.unit
def test_cv2front_processing_builds_anodic_and_cathodic_descriptors():
    data = [{
        "_id": "cv-1",
        "mol_id": "M1",
        "data": {
            "conditions": {"scan_rate": {"value": 0.1, "unit": "V/s"}},
            "e_half": [{"value": 0.25, "unit": "V"}],
        },
    }]
    cv = CV2Front(backend_data=data, run_processing=False, insert=False, verbose=0, run_anodic=True)
    calls = []

    def calculate(electron_num, curve_type):
        calls.append((electron_num, curve_type))
        return (1e-5, 2e-6) if curve_type == "cathodic" else (3e-5, 4e-6)

    cv.meta_calcs_func = calculate
    cv.process()

    assert calls == [(1, "cathodic"), (1, "anodic")]
    assert len(cv.meta_dict["diffusion_coefficient"]) == 2
    assert [entry["notes"] for entry in cv.meta_dict["charge_transfer_rate"]] == ["cathodic", "anodic"]


@pytest.mark.unit
def test_cv2front_cv_calculators_filter_scan_rates_and_calculate(monkeypatch):
    seen = {}

    class Diffusion:
        def __init__(self, connector):
            seen["diffusion_connector"] = connector

        def calculate(self, data, **kwargs):
            seen["processed"] = data
            seen["diffusion_options"] = kwargs
            return ("1.0e-5", 1e-5)

    class Transfer:
        def __init__(self, connector):
            seen["transfer_connector"] = connector

        def calculate(self, data, **kwargs):
            seen["transfer_data"] = data
            return 2e-6

    monkeypatch.setattr("d3tales_api.D3database.back2front.CVDiffusionCalculator", Diffusion)
    monkeypatch.setattr("d3tales_api.D3database.back2front.CVChargeTransferCalculator", Transfer)
    records = [
        {
            "_id": "cv-1",
            "mol_id": "M1",
            "data": {
                "conditions": {
                    "scan_rate": {"value": 0.1, "unit": "V/s"},
                    "temperature": {"value": 298, "unit": "K"},
                },
                "peak_splittings": [0.2],
                "middle_sweep": [],
            },
        },
        {
            "_id": "cv-2",
            "mol_id": "M1",
            "data": {
                "conditions": {
                    "scan_rate": {"value": 1.0, "unit": "V/s"},
                    "temperature": {"value": 298, "unit": "K"},
                },
            },
        },
    ]
    cv = CV2Front(backend_data=records, run_processing=False, insert=False, verbose=0, max_scan_rate=0.5)
    diffusion, transfer = cv.cv_meta_calcs(electron_num=1, curve_type="cathodic")

    assert (diffusion, transfer) == (1e-5, 2e-6)
    assert len(seen["processed"]) == 1
    assert seen["processed"][0]["n"] == 1
    assert seen["processed"][0]["diffusion"] == 1e-5
    assert seen["diffusion_options"] == {"sci_notation": True, "cathodic_anodic": "cathodic"}


@pytest.mark.unit
def test_cv2front_microelectrode_calculators_use_single_scan(monkeypatch):
    seen = {}

    class Diffusion:
        def __init__(self, connector):
            seen["diffusion_connector"] = connector

        def calculate(self, data):
            seen["input"] = data.copy()
            return 1e-6

    class Transfer:
        def __init__(self, connector):
            seen["transfer_connector"] = connector

        def calculate(self, data, **kwargs):
            seen["transfer_input"] = data.copy()
            return 2e-7

    monkeypatch.setattr("d3tales_api.D3database.back2front.CVDiffusionCalculatorMicro", Diffusion)
    monkeypatch.setattr("d3tales_api.D3database.back2front.CVChargeTransferCalculatorMicro", Transfer)
    data = [{
        "_id": "micro-1",
        "mol_id": "M1",
        "data": {"conditions": {"temperature": 298}, "e_half": {"value": 0.2}},
    }]
    cv = CV2Front(backend_data=data, run_processing=False, insert=False, verbose=0, micro_electrodes=True)
    assert cv.micro_meta_calcs(curve_type="cathodic") == (1e-6, 2e-7)
    assert seen["input"]["num_electrodes"] == 1
    assert seen["transfer_input"]["diffusion"] == 1e-6
