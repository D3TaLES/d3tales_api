import copy
import urllib
import hashlib
import pubchempy as pcp

from d3tales_api.database_info import db_info
from rdkit.Chem.AllChem import ComputeMolVolume
from rdkit.Chem import MolFromSmiles, MolToSmiles
from d3tales_api.Calculators.calculators import *
from d3tales_api.D3database.d3database import DBconnector, FrontDB
from d3tales_api.Processors.info_from_smiles import GenerateMolInfo

DEFAULT_SOLV = {"name": "Acetonitrile", "model": "implicit_solvent", "dielectric_constant": 35.688}
RMSD_DEFAULT = False
FRONT_COLL = "base"
SCHEMA_VERSION = "v2"


class Gaus2FrontCharacterization:
    """
    Update frontend db with backend db data for a particular set of data.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, _id, conditions, calculation_type=None, insert=False, data=None, rmsd=RMSD_DEFAULT,
                 gs_charge=None, all_props=False, verbose=2, schema_version=SCHEMA_VERSION):
        """

        :param _id: str, molecule ID
        :param calculation_type: str, calculation type
        :param conditions: dict, calculation conditions
        :param insert: bool, insert generated data to the frontend D3TaLES database if True
        """
        # connect to databases
        self.front_dbc = DBconnector(db_info.get("frontend"))
        self.front_coll = self.front_dbc.get_collection(FRONT_COLL)
        self.back_dbc = DBconnector(db_info.get("backend"))
        self.back_coll = self.back_dbc.get_collection("computation")
        self.schema_version = schema_version
        self.all_props = all_props
        self.verbose = verbose

        # basic variables
        self.id = _id
        self.calculation_type = calculation_type
        self.bare_conditions = copy.deepcopy(conditions)
        self.solvent = self.bare_conditions.pop("solvent") if self.bare_conditions.get("solvent") else DEFAULT_SOLV
        self.groundState_charge = self.front_coll.find_one({"_id": self.id})["mol_info"][
            "groundState_charge"] if gs_charge is None else gs_charge
        print("Starting query for all backend {} data...".format(self.id)) if verbose else None
        self.mol_data = list(self.back_coll.find({"mol_id": self.id}) or [])
        if data:
            self.mol_data.append(data)
        self.mol_hashes = [d.get("_id") for d in self.mol_data]
        print("Found {} backend data documents for {}...".format(len(self.mol_data), self.id)) if verbose else None

        # all calculations for this molecule
        self.all_calcs = ['opt_groundState', 'freq_groundState', 'solv_energy_gsgs', 'tddft_groundState',
                          'opt_cation1', 'freq_cation1', 'solv_energy_c1c1', 'tddft_cation1',
                          'opt_cation2', 'freq_cation2', 'solv_energy_c2c2', 'tddft_cation2',
                          'opt_anion1', 'freq_anion1', 'solv_energy_a1a1', 'tddft_anion1',
                          'opt_anion2', 'freq_anion2', 'solv_energy_a2a2', 'tddft_anion2',
                          'energy_gsc1', 'energy_gsa1', 'energy_c1gs', 'energy_a1gs', 'energy_c1c2', 'energy_a1a2']
        self.species_dict = {
            "anion2": -2,
            "anion1": -1,
            "groundState": 0,
            "cation1": 1,
            "cation2": 2,
        }
        # Get data for all calculations for this molecule
        for calc_type in self.all_calcs:
            if "solv_" in calc_type:
                setattr(self, calc_type, self.find_data(calc_type=calc_type, solvent=self.solvent))
            else:
                setattr(self, calc_type, self.find_data(calc_type=calc_type))

        # Set up mol properties to calculate
        self.mol_properties = [
            self.omega,
            self.hole_reorganization_energy, self.electron_reorganization_energy,
            self.relaxation_groundState_cation1,
            self.relaxation_cation1_groundState, self.relaxation_groundState_anion1,
            self.relaxation_anion1_groundState,
            self.vertical_ionization_energy, self.vertical_ionization_energy_2,
            self.vertical_electron_affinity, self.vertical_electron_affinity_2,
            self.adiabatic_ionization_energy, self.adiabatic_ionization_energy_2,
            self.adiabatic_electron_affinity, self.adiabatic_electron_affinity_2,
            self.oxidation_potential, self.reduction_potential
        ]
        if rmsd:
            self.mol_properties.extend(
                [self.rmsd_groundState_cation1, self.rmsd_cation1_cation2, self.rmsd_groundState_anion1,
                 self.rmsd_anion1_anion2])

        if insert:
            self.insert()

    def insert(self, **kwargs):
        # Get and insert species data
        if not self.calculation_type:
            print(
                f"No data for {self.id} inserted because no calculation type specified.") if self.verbose > 1 else None
            return None
        species_dict = self.species_descriptors(self.calculation_type, **kwargs)
        FrontDB(schema_layer="species_characterization", instance=species_dict, _id=self.id,
                verbose=self.verbose, schema_version=self.schema_version)

        # Get and insert mol data
        mol_char_dict = {}
        for prop in self.mol_properties:
            try:
                mol_char_dict.update(prop())
            except ValueError as e:
                print(e) if self.verbose > 1 else None
        FrontDB(schema_layer="mol_characterization", instance=mol_char_dict, _id=self.id,
                verbose=self.verbose, schema_version=self.schema_version)

    def insert_all_species(self, **kwargs):
        # Get and insert all species data
        print(f"Inserting all data for molecule {self.id}...") if self.verbose else None
        all_species_data = {s: {} for s in self.species_dict.keys()}
        for calc_type in self.all_calcs:
            species_dict = self.species_descriptors(calc_type, **kwargs)
            if species_dict:
                for species, s_data in species_dict.items():
                    all_species_data[species].update(s_data)

        FrontDB(schema_layer="species_characterization", instance=all_species_data, _id=self.id,
                verbose=self.verbose, schema_version=self.schema_version)
        [print(s, ':\t', list(s_data.keys())) for s, s_data in all_species_data.items()] if self.verbose > 1 else None

    def insert_all_mol(self):
        # Get and insert mol data
        self.all_props = True
        mol_char_dict = {}
        for prop in self.mol_properties:
            try:
                mol_char_dict.update(prop())
            except ValueError as e:
                print(e) if self.verbose > 1 else None
        FrontDB(schema_layer="mol_characterization", instance=mol_char_dict, _id=self.id,
                verbose=self.verbose, schema_version=self.schema_version)

    def species_descriptors(self, calculation_type, skip_geom=False):
        """Descriptor dict for species descriptors"""
        species_data, _hash = getattr(self, calculation_type, (None, None))
        if not species_data:
            print(f"No {calculation_type} data for {self.id}") if self.verbose > 1 else None
            return {}
        conditions = species_data.get("conditions")
        charge = species_data.get("charge")

        species = {v: k for k, v in self.species_dict.items()}[charge - self.groundState_charge]
        data_dict = {
            species: {
                "charge": charge,
                "spin_multiplicity": species_data["spin_multiplicity"],
                "is_groundState": species_data["is_groundState"],
            }
        }

        if "solv_energy" in calculation_type:
            try:
                data_dict[species].update(self.solvation_energy(species))
            except ValueError as e:
                print(e) if self.verbose > 1 else None
        elif "opt" in calculation_type:
            data_dict[species].update({"homo": self.prop_entry(hashes=[_hash], conditions=conditions,
                                                               value=species_data["homo"]["value"],
                                                               unit=species_data["homo"]["unit"])})
            data_dict[species].update({"lumo": self.prop_entry(hashes=[_hash], conditions=conditions,
                                                               value=species_data["lumo"]["value"],
                                                               unit=species_data["lumo"]["unit"])})
            if species_data["lumo"]["unit"] == species_data["homo"]["unit"]:
                gap = species_data["lumo"]["value"] - species_data["homo"]["value"]
                unit = species_data["lumo"]["unit"]
                data_dict[species].update({"homo_lumo_gap": self.prop_entry(hashes=[_hash], conditions=conditions,
                                                                            value=gap, unit=unit)})
            dipole_moment = self.prop_entry(hashes=[_hash], conditions=conditions,
                                            value=np.linalg.norm(
                                                np.asarray(species_data["scf_dipole_moment"]["value"])),
                                            unit=species_data["scf_dipole_moment"]["unit"])
            data_dict[species].update({"dipole_moment": dipole_moment})
            for opt_descriptor in ["radical_buried_vol", "radical_spin", "radical_stability_score"]:
                if species_data.get(opt_descriptor):
                    species_data[opt_descriptor].update({"conditions": conditions})
                    data_dict[species].update(
                        {opt_descriptor: self.prop_entry(hashes=[_hash], **species_data[opt_descriptor])})
            if not skip_geom:
                pmg_mol = Molecule.from_sites([Site.from_dict(sd) for sd in species_data["geometry"]])
                try:
                    rdmol = pmgmol_to_rdmol(pmg_mol)[0]
                    vol = ComputeMolVolume(rdmol)
                    globular_volume = self.prop_entry(hashes=[_hash], conditions=conditions, value=round(vol, 3),
                                                      unit="A^3")
                    data_dict[species].update({"globular_volume": globular_volume})
                except Chem.rdchem.AtomValenceException:
                    pass
                geom_data = self.prop_entry(hashes=[_hash], conditions=conditions, sites=species_data["geometry"])
                data_dict[species].update({"geometry": geom_data})
        elif "tddft" in calculation_type:
            if species_data.get("dipole_moment"):
                dipole_moment = self.prop_entry(hashes=[_hash], conditions=conditions,
                                                value=np.linalg.norm(
                                                    np.asarray(species_data["scf_dipole_moment"]["value"])),
                                                unit=species_data["scf_dipole_moment"]["unit"])
                data_dict[species].update({"dipole_moment": dipole_moment})
            if species_data["excitations"].get("Singlet"):
                singlet_data = self.prop_entry(hashes=[_hash], conditions=conditions,
                                               excitations=species_data["excitations"]["Singlet"])
                data_dict[species].update({"singlet_states": singlet_data})
            if species_data["excitations"].get("Triplet"):
                triplet_data = self.prop_entry(hashes=[_hash], conditions=conditions,
                                               excitations=species_data["excitations"]["Triplet"])
                data_dict[species].update({"triplet_states": triplet_data})
            if species_data.get("singlet_plotting", ):
                plotting_data = self.prop_entry(hashes=[_hash], conditions=conditions,
                                                plotting_data=species_data["singlet_plotting"])
                data_dict[species].update({"spectra": plotting_data})
        else:
            return {}
        return data_dict

    @classmethod
    def from_data(cls, processing_data, **kwargs):
        """
        Generate data class from data dict

        :param processing_data: dict, data dict
        :return: data class
        """
        m_id = processing_data.get("mol_id")
        calculation_type = processing_data.get("calculation_type")
        conditions = processing_data.get("data", {}).get("conditions")
        return cls(_id=m_id, calculation_type=calculation_type, conditions=conditions, data=processing_data, **kwargs)

    @staticmethod
    def prop_entry(hashes=None, conditions=None, **kwargs):
        prop_dict = {
            "source_hash_ids": hashes or [],
            "conditions": conditions or {},
        }
        prop_dict.update(kwargs)
        return prop_dict

    def mol_descriptor_dict(self, value, unit="", hashes=None, name="", order=1, condition_addition=None):
        """
        Generate descriptor dictionary in accordance with D3TaLES schema

        :param value: data value
        :param unit: str, unit
        :param hashes: list, hash ids
        :param name: str, property name
        :param order: int, property order
        :param condition_addition: dict, additional data to add to conditions
        :return: descriptor dictionary
        """
        conditions = copy.deepcopy(self.bare_conditions)
        if condition_addition:
            conditions.update(condition_addition)
        return {
            name: self.prop_entry(hashes=hashes, conditions=conditions, value=value, unit=unit, order=order)
        }

    def find_data(self, calc_type, solvent=None, conditions=None):
        """
        Find calculation data

        :param calc_type: str, calculation type
        :param solvent: str, solvent
        :return: [data dict, hash ID]
        """
        data_conditions = conditions or copy.deepcopy(self.bare_conditions)
        if solvent:
            data_conditions["solvent"] = solvent
        _hash = self.hash_id(self.id, calc_type, data_conditions)
        try:
            db_request = [i for i in self.mol_data if i.get("_id") == _hash][0]
            h_data = db_request.get("data", )
            return [h_data, _hash]
        except IndexError:
            return [None, _hash]

    def get_data(self, calculations, disregard_missing=False):
        """
        Get all data from a list of calculations

        :param calculations: list, list of calculation names
        :param disregard_missing: bool, disregard missing jobs if True
        :return: [list of hash IDs, dict of calculation name with associated data]
        """
        if (self.calculation_type not in calculations) and not self.all_props:
            raise ValueError("Calculation does not use this calculation type")
        calcs_data = [getattr(self, c) for c in calculations]
        jobs, h_ids = [c[0] for c in calcs_data], [c[1] for c in calcs_data]
        c_data = {calculations[i]: d for i, d in enumerate(jobs)}
        if None in jobs and not disregard_missing:
            missing_calcs = [k for k, v in c_data.items() if v is None]
            raise ValueError(
                "The following calculations needed for this property are missing: {}".format(", ".join(missing_calcs)))
        return h_ids, c_data

    @staticmethod
    def hash_id(_id, calc_type, gaus_conditions):
        """
        Generate hash ID

        :param _id: str, molecule ID
        :param calc_type: str, calculation type
        :param gaus_conditions: dict, calculation conditions
        :return:
        """
        hash_dict = {
            "_id": _id,
            "calculation_type": calc_type,
            "conditions": gaus_conditions,
        }
        dhash = hashlib.md5()
        encoded = json.dumps(hash_dict, sort_keys=True).encode()
        dhash.update(encoded)
        return dhash.hexdigest()

    def solvation_energy(self, species):
        """Descriptor dict for the solvation energy"""
        species_abbrev_dict = {"anion2": "a2", "anion1": "a1", "groundState": "gs", "cation1": "c1", "cation2": "c2"}
        species_abbrev = species_abbrev_dict.get(species)
        gas_phase = "opt_{}".format(species)
        solv = "solv_energy_{}{}".format(species_abbrev, species_abbrev)

        h_ids, c_data = self.get_data([gas_phase, solv])
        connector = {"energy_final": f"{solv}.scf_total_energy.value",
                     "energy_initial": f"{gas_phase}.scf_total_energy.value"}
        solv_eng = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(solv_eng, unit='eV', hashes=h_ids, name="solvation_energy",
                                        condition_addition={"solvent": self.solvent})

    def omega(self):
        """Descriptor dict for the hole reorganization energy"""
        w_conditions = copy.deepcopy(self.bare_conditions)
        w_conditions.pop("tuning_parameter") if w_conditions.get("tuning_parameter") else None
        c_data, h_id = self.find_data("wtuning", conditions=w_conditions)
        if c_data:
            omega = math.floor(c_data.get("omega", 0) * 10000) / 10000
            return {"omega": self.prop_entry(hashes=[h_id], conditions=w_conditions, value=omega, unit="")}

        raise ValueError

    def hole_reorganization_energy(self):
        """Descriptor dict for the hole reorganization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsc1", "opt_cation1", "energy_c1gs"])
        connector = {"gs_opt": "opt_groundState.scf_total_energy.value",
                     "ion_opt": "opt_cation1.scf_total_energy.value",
                     "gs_energy": "energy_c1gs.scf_total_energy.value",
                     "ion_energy": "energy_gsc1.scf_total_energy.value"}
        energy = ReorganizationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="hole_reorganization_energy")

    def electron_reorganization_energy(self):
        """Descriptor dict for the electron reorganization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsa1", "opt_anion1", "energy_a1gs"])
        connector = {"gs_opt": "opt_groundState.scf_total_energy.value",
                     "ion_opt": "opt_anion1.scf_total_energy.value",
                     "gs_energy": "energy_a1gs.scf_total_energy.value",
                     "ion_energy": "energy_gsa1.scf_total_energy.value"}
        energy = ReorganizationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="electron_reorganization_energy")

    def relaxation_groundState_cation1(self):
        """Descriptor dict for the relaxation energy of the ground state geometry to the +1 cation geometry"""
        h_ids, c_data = self.get_data(["energy_gsc1", "opt_cation1"])
        connector = {"opt_energy": "opt_cation1.scf_total_energy.value",
                     "energy": "energy_gsc1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_groundState_cation1")

    def relaxation_cation1_groundState(self):
        """Descriptor dict for the relaxation energy of the +1 cation geometry to the ground state geometry"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_c1gs"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_c1gs.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_cation1_groundState")

    def relaxation_groundState_anion1(self):
        """Descriptor dict for the relaxation energy of the ground state geometry to the -1 anion geometry"""
        h_ids, c_data = self.get_data(["energy_gsa1", "opt_anion1"])
        connector = {"opt_energy": "opt_anion1.scf_total_energy.value",
                     "energy": "energy_gsa1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_groundState_anion1")

    def relaxation_anion1_groundState(self):
        """Descriptor dict for the relaxation energy of the ground -1 anion to the ground state geometry"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_a1gs"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_a1gs.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_anion1_groundState")

    def vertical_ionization_energy(self):
        """Descriptor dict for the vertical ionization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsc1"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_gsc1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_ionization_energy")

    def vertical_ionization_energy_2(self):
        """Descriptor dict for the vertical ionization energy"""
        h_ids, c_data = self.get_data(["opt_cation1", "energy_c1c2"])
        connector = {"opt_energy": "opt_cation1.scf_total_energy.value",
                     "energy": "energy_c1c2.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_ionization_energy_2")

    def vertical_electron_affinity(self):
        """Descriptor dict for the """
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsa1"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_gsa1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_electron_affinity")

    def vertical_electron_affinity_2(self):
        """Descriptor dict for the """
        h_ids, c_data = self.get_data(["opt_anion1", "energy_a1a2"])
        connector = {"opt_energy": "opt_anion1.scf_total_energy.value",
                     "energy": "energy_a1a2.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_electron_affinity_2")

    def adiabatic_ionization_energy(self):
        """Descriptor dict for the adiabatic ionization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_cation1"])
        connector = {"energy_final": "opt_cation1.scf_total_energy.value",
                     "energy_initial": "opt_groundState.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_ionization_energy")

    def adiabatic_ionization_energy_2(self):
        """Descriptor dict for the second adiabatic ionization energy"""
        h_ids, c_data = self.get_data(["opt_cation1", "opt_cation2"])
        connector = {"energy_final": "opt_cation2.scf_total_energy.value",
                     "energy_initial": "opt_cation1.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_ionization_energy_2")

    def adiabatic_electron_affinity(self):
        """Descriptor dict for the adiabatic electron affinity"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_anion1"])
        connector = {"energy_final": "opt_anion1.scf_total_energy.value",
                     "energy_initial": "opt_groundState.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_electron_affinity")

    def adiabatic_electron_affinity_2(self):
        """Descriptor dict for the second adiabatic electron affinity"""
        h_ids, c_data = self.get_data(["opt_anion1", "opt_anion2"])
        connector = {"energy_final": "opt_anion2.scf_total_energy.value",
                     "energy_initial": "opt_anion1.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_electron_affinity_2")

    def rmsd_groundState_cation1(self):
        """Descriptor dict for the RMSD between the ground state and +1 cation geometries"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_cation1"])
        connector = {"geom_initial": "opt_groundState.geometry",
                     "geom_final": "opt_cation1.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_groundState_cation1")

    def rmsd_cation1_cation2(self):
        """Descriptor dict for the RMSD between the +1 cation and +2 cation geometries"""
        h_ids, c_data = self.get_data(["opt_cation2", "opt_cation1"])
        connector = {"geom_initial": "opt_cation1.geometry",
                     "geom_final": "opt_cation2.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_cation1_cation2")

    def rmsd_groundState_anion1(self):
        """Descriptor dict for the RMSD between the ground state and -1 anion geometries"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_anion1"])
        connector = {"geom_initial": "opt_groundState.geometry",
                     "geom_final": "opt_anion1.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_groundState_anion1")

    def rmsd_anion1_anion2(self):
        """Descriptor dict for the RMSD between the -1 anion and -2 anion geometries"""
        h_ids, c_data = self.get_data(["opt_anion1", "opt_anion2"])
        connector = {"geom_initial": "opt_anion1.geometry",
                     "geom_final": "opt_anion2.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_cation1_anion2")

    def oxidation_potential(self, electrode="standard_hydrogen_electrode", num_electrons=1):
        """Descriptor dict for the oxidation potential"""
        h_ids, c_data = self.get_data(
            ["opt_groundState", "freq_groundState", "solv_energy_gsgs", "opt_cation1", "freq_cation1",
             "solv_energy_c1c1"])
        c_data.update({"electrode": electrode, "num_electrons": num_electrons})
        connector = {"init_eng": "opt_groundState.scf_total_energy.value",
                     "init_corr": "freq_groundState.gibbs_correction.value",
                     "init_eng_solv": "solv_energy_gsgs.scf_total_energy.value",
                     "fin_eng": "opt_cation1.scf_total_energy.value",
                     "fin_corr": "freq_cation1.gibbs_correction.value",
                     "fin_eng_solv": "solv_energy_c1c1.scf_total_energy.value",
                     "electrode": "electrode", "num_electrons": "num_electrons"}
        energy = -RedoxPotentialCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="oxidation_potential",
                                        condition_addition={"reference_electrode": electrode, "solvent": self.solvent})

    def reduction_potential(self, electrode="standard_hydrogen_electrode", num_electrons=1):
        """Descriptor dict for the reduction potential"""
        h_ids, c_data = self.get_data(
            ["opt_groundState", "freq_groundState", "solv_energy_gsgs", "opt_anion1", "freq_anion1",
             "solv_energy_a1a1"])
        c_data.update({"electrode": electrode, "num_electrons": num_electrons})
        connector = {"init_eng": "opt_groundState.scf_total_energy.value",
                     "init_corr": "freq_groundState.gibbs_correction.value",
                     "init_eng_solv": "solv_energy_gsgs.scf_total_energy.value",
                     "fin_eng": "opt_anion1.scf_total_energy.value",
                     "fin_corr": "freq_anion1.gibbs_correction.value",
                     "fin_eng_solv": "solv_energy_a1a1.scf_total_energy.value",
                     "electrode": "electrode", "num_electrons": "num_electrons"}
        energy = RedoxPotentialCalc(connector=connector).calculate(c_data)
        return self.mol_descriptor_dict(energy, unit='eV', hashes=h_ids, name="reduction_potential",
                                        condition_addition={"reference_electrode": electrode, "solvent": self.solvent})

    def get_all_data(self):
        """
        Get all calculations data for object molecule
        :return: dictionary containing all backend data for the object molecule
        """
        h_ids, c_data = self.get_data(self.all_calcs, disregard_missing=True)
        c_data.update({"h_ids": h_ids})
        return c_data


class DOI2Front:
    """
    Update frontend db with backend db data for a particular set of data.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, doi=None, backend_data=None, insert=True, nlp_group="Sarkar"):
        """

        :param doi: str, molecule ID
        :param backend_data: dict, calculation data
        :param insert: bool, insert generated data to the frontend D3TaLES database if True
        :param nlp_group: str, name of NLP group
        """

        # connect to databases
        self.front_coll = DBconnector(db_info.get("frontend")).get_collection(FRONT_COLL)
        self.back_coll = DBconnector(db_info.get("backend")).get_collection("nlp")

        # Basic variables
        if not doi and not backend_data:
            raise ValueError(
                "The DOI2Front class requires either the 'doi' kwarg or 'backend_data' kwarg. Neither were provided. ")
        self.doi = doi or backend_data.get("_id")
        self.backend_data = backend_data if backend_data else self.back_coll.find_one({"_id": self.doi})
        if not self.backend_data:
            raise ValueError("No backend data found for DOI {}.".format(self.doi))
        self.raw_mol_data = self.backend_data.get("extracted_molecules", [])

        # Set data for this NLP extraction
        self.extracted_mol_data = {self.get_mol_id(d, nlp_group): self.get_literature_data(d) for d in
                                   self.raw_mol_data}

        self.mol_ids = set([i for i in self.extracted_mol_data.keys() if i])
        if insert:
            for i in self.mol_ids:
                FrontDB(instance=self.extracted_mol_data[i], _id=i)

    def get_literature_data(self, mol_data):
        extracted_properties = mol_data.get("extracted_properties", "")
        return {
            "literature_data": {
                "related_literature": [self.doi],
                "extracted_properties": self.check_extr_props(extracted_properties)
            }
        }

    def check_extr_props(self, extracted_props):
        for _, p_list in extracted_props.items():
            for p in p_list:
                if not p.get("conditions", {}).get("doi"):
                    p["conditions"]["doi"] = self.doi
        return extracted_props

    @staticmethod
    def cir_convert(mol_name):
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + urllib.parse.quote(mol_name) + '/smiles'
        try:
            return urllib.request.urlopen(url).read().decode('utf8')
        except urllib.error.HTTPError:
            return None

    def get_mol_id(self, mol_data, nlp_group):
        mol_name = mol_data.get("molecule_name", "")
        pub_mols = pcp.get_compounds(mol_name, 'name')
        if len(pub_mols) == 0:
            smiles = self.cir_convert(mol_name)
            if not smiles:
                warnings.warn("No PubChem or CIR molecules found for molecule {}".format(mol_name))
        elif len(pub_mols) > 1:
            warnings.warn("Multiple PubChem molecules found for molecule {}: {}".format(mol_name, ", ".join(
                [p.iupac_name for p in pub_mols])))
        else:
            smiles = pub_mols[0].isomeric_smiles
            rdkmol = MolFromSmiles(smiles)
            clean_smiles = MolToSmiles(rdkmol)
            db_check = FrontDB(smiles=clean_smiles).check_if_in_db()
            if db_check:
                return db_check
            instance = GenerateMolInfo(clean_smiles, origin_group=nlp_group).mol_info_dict
            db_insertion = FrontDB(schema_layer='mol_info', instance=instance, smiles=clean_smiles, group=nlp_group)
            return db_insertion.id

    @classmethod
    def from_json(cls, json_path, **kwargs):
        """
        Generate data class from data dict

        :param json_path: str, path to JSON file
        :return: data class
        """
        with open(json_path) as fn:
            processing_data = json.load(fn)
        doi = processing_data.get("_id")
        return cls(doi=doi, backend_data=processing_data, **kwargs)


class CV2Front:
    """
    Update frontend db with backend db data for a particular set of data.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, id_list=None, backend_data=None, metadata_dict=None, mol_id=None, e_half_scan_rate=0.1, max_scan_rate=None,
                 run_anodic=False, run_processing=True, insert=True, redox_event="oxidation", micro_electrodes=False,
                 verbose=1, backend_db="backend"):
        """
        :param id_list: list, list of backend ids
        :param backend_data: list, list of JSON containing processed backend data
        :param metadata_dict: dict, default metadata dictionary
        :param mol_id: str, default molecule ID for data if none found
        :param e_half_scan_rate: float, scan rate at which to extract e half
        :param run_anodic: bool, run processing for anodic sweep if True
        :param run_processing: bool, run processing (takes a few minutes) if True
        :param insert: bool, insert generated data to the frontend D3TaLES database if True
        :param micro_electrodes: bool,
        :param verbose: int, level of verbosity
        :param backend_db: str, backend DB key word
        """

        # connect to databases
        self.front_coll = DBconnector(db_info.get("frontend")).get_collection(FRONT_COLL)
        self.back_coll = DBconnector(db_info.get(backend_db)).get_collection("experimentation")

        # Basic variables
        self.verbose = verbose
        if not id_list and not backend_data:
            raise ValueError(
                "The CV2Front class requires either the 'id_list' kwarg or 'backend_data' kwarg. Neither were provided. ")
        self.multi_data = backend_data or [self.back_coll.find_one({"_id": i}) for i in id_list]
        if not self.multi_data or not self.multi_data[0]:
            raise ValueError("No backend data found for IDs {} and no backend data provided.".format(id_list))
        self.conditions = self.get_conditions()
        self.p_ids = [d.get("_id") for d in self.multi_data]
        self.e_half_scan_rate = e_half_scan_rate
        self.max_scan_rate = max_scan_rate
        self.run_anodic = run_anodic
        self.redox_event = redox_event

        # Check Mol ID
        mol_ids = list(set(filter(lambda i: i is not None, [d.get("mol_id") for d in self.multi_data])))
        if len(mol_ids) > 1:
            raise ValueError(
                "Not all submitted data is associated with the same molecule! Submitted data are associated"
                "with molecules {}.".format(", ".join(mol_ids)))
        self.mol_id = mol_ids[0] if mol_ids else mol_id

        # Start Metadict
        self.meta_dict = metadata_dict or {}

        # Set meta calcs function
        self.meta_calcs_func = self.micro_meta_calcs if micro_electrodes else self.cv_meta_calcs

        if run_processing:
            self.process()

        [setattr(self, k, v) for k, v in self.meta_dict.items()]

        if insert:
            FrontDB(instance={"experiment_data.mol_characterization": self.meta_dict}, _id=self.mol_id)
            FrontDB(instance={"experiment_data.experiment_ids": self.p_ids}, _id=self.mol_id)

    def process(self):
        diffusions, charge_transfers = [], []
        if self.verbose:
            print("STARTING META PROPERTY PROCESSING...This could take a few minutes.")

        self.meta_dict.update({f"{self.redox_event}_potential": self.e_halfs})

        for i, e_half in enumerate(self.e_halfs):
            if self.verbose:
                print("E 1/2s: ", ", ".join([str(e.get("value")) for e in self.e_halfs]))
            c_diffusion_coef, c_transfer_rate = self.meta_calcs_func(electron_num=i + 1, curve_type="cathodic")
            diffusions.append(self.prop_dict(c_diffusion_coef, unit="cm^2/s", order=i + 1, notes="cathodic"))
            charge_transfers.append(self.prop_dict(c_transfer_rate, unit="cm/s", order=i + 1, notes="cathodic"))

            if self.run_anodic:
                a_diffusion_coef, a_transfer_rate = self.meta_calcs_func(electron_num=i + 1, curve_type="anodic")
                diffusions.append(self.prop_dict(a_diffusion_coef, unit="cm^2/s", order=i + 1, notes="anodic"))
                charge_transfers.append(self.prop_dict(a_transfer_rate, unit="cm/s", order=i + 1, notes="anodic"))

        self.meta_dict.update({"diffusion_coefficient": diffusions, "charge_transfer_rate": charge_transfers})

    def get_conditions(self):
        all_cond = [d.get("data", {}).get("conditions") for d in self.multi_data]
        if not all_cond:
            raise ValueError("No experiments have associated conditions.")
        return {k: v for k, v in all_cond[0].items() if len(set([self.generate_hash(c[k]) for c in all_cond])) == 1}

    @staticmethod
    def generate_hash(hash_dict):
        dhash = hashlib.md5()
        encoded = json.dumps(hash_dict, sort_keys=True).encode()
        dhash.update(encoded)
        return dhash.hexdigest()

    def prop_dict(self, value, unit="", order=1, conditions=None, hashes=None, notes=None):
        """
        Generate descriptor dictionary in accordance with D3TaLES schema

        :param value: data value
        :param unit: str, unit
        :param order: int, property order
        :param conditions: list, property conditions
        :param hashes: list, property hash ids
        :param notes: str, property notes
        :return: descriptor dictionary
        """
        cond = copy.deepcopy(conditions or self.conditions)
        prop_dict = {
            "source_hash_ids": hashes or self.p_ids,
            "conditions": cond,
            "order": order,
            "value": value,
            "unit": unit
        }
        if notes:
            prop_dict.update({"notes": notes})
        return prop_dict

    def cv_meta_calcs(self, electron_num=1, curve_type="anodic"):
        processed_data = []
        for i in self.multi_data:
            data_dict = i.get('data', {})
            scan_rate = unit_conversion(data_dict.get("conditions", {}).get("scan_rate", 0), default_unit='V/s')
            if self.max_scan_rate:
                if scan_rate > self.max_scan_rate:
                    print(f"WARNING! Data with scan rate {scan_rate} skipped because max scan rate is {self.max_scan_rate}")
                    continue
            data_dict["n"] = electron_num
            processed_data.append(data_dict)
        connector = {
            "n": "n",
            "X": "peak_splittings.{}".format(electron_num - 1),
            "v": "conditions.scan_rate",
            "T": "conditions.temperature",
            "C": "conditions.redox_mol_concentration",
            "A": "conditions.working_electrode_surface_area",

            "D": "diffusion",
            "scan_data": "middle_sweep",
        }

        # Calculate diffusion constant
        if self.verbose:
            print("Calculating {} diffusion coefficient for oxidation {}...".format(curve_type, electron_num))
        diffusion_cal = CVDiffusionCalculator(connector=connector)
        diffusion_coef = diffusion_cal.calculate(processed_data, sci_notation=True, cathodic_anodic=curve_type)

        # Calculate charge transfer rates
        if self.verbose:
            print("Calculating {} charge transfer rate for oxidation {}...".format(curve_type, electron_num))
        [d.update({"diffusion": float(diffusion_coef[1])}) for d in processed_data]
        transfer_cal = CVChargeTransferCalculator(connector=connector)
        transfer_rate = transfer_cal.calculate(processed_data, sci_notation=True)
        return diffusion_coef[1], transfer_rate

    def micro_meta_calcs(self, electron_num=1, curve_type="anodic"):
        if len(self.multi_data) > 1:
            raise warnings.warn("Micro electrode calculators only use one CV scan. While {} scans were inputted, "
                                "only the first one is being used for calculations. ".format(len(self.multi_data)))
        instance = self.multi_data[0]
        instance.update(dict(num_electrodes=1))
        connector = {
            "i_ss": "data.i_ss",
            "r": "data.conditions.working_electrode_radius",
            "C": "data.conditions.redox_mol_concentration",
            "n": "num_electrodes",

            "e_half": "data.e_half",
            "e_ref": "data.e_ref",
            "T": "data.conditions.temperature",
            "D": "diffusion",
        }

        # Calculate diffusion constant
        if self.verbose:
            print("Calculating {} diffusion coefficient for oxidation {}...".format(curve_type, electron_num))
        diff_calc = CVDiffusionCalculatorMicro(connector=connector)
        diffusion_coef = diff_calc.calculate(instance)

        # Calculate charge transfer rates
        if self.verbose:
            print("Calculating {} charge transfer rate for oxidation {}...".format(curve_type, electron_num))
        transfer_cal = CVChargeTransferCalculatorMicro(connector=connector)
        instance.update(dict(diffusion=diffusion_coef))
        # raise Exception(instance)
        transfer_rate = transfer_cal.calculate(instance, sci_notation=True)
        return float(diffusion_coef), float(transfer_rate)

    @property
    def e_halfs(self):
        e_half_data = {}
        for d in self.multi_data:
            d_sr = d.get("data", {}).get("conditions", {}).get("scan_rate", {}).get("value", None)
            if self.e_half_scan_rate == d_sr or len(self.multi_data) == 1:
                e_half_data = d
                break
        self.e_halfs_id = e_half_data.get("_id")
        self.e_halfs_conditions = e_half_data.get("data", {}).get("conditions")
        e_halfs = e_half_data.get("data", {}).get("e_half")
        print("E Half: ", e_halfs)
        if isinstance(e_halfs, list):
            return [self.prop_dict(e_half.get("value"), unit=e_half.get("unit"), order=i + 1, hashes=[self.e_halfs_id],
                                   conditions=self.e_halfs_conditions) for i, e_half in enumerate(e_halfs)]
        elif isinstance(e_halfs, dict):
            return [self.prop_dict(e_halfs.get("value"), unit=e_halfs.get("unit"), order=1, hashes=[self.e_halfs_id],
                                   conditions=self.e_halfs_conditions)]

    @classmethod
    def from_json(cls, json_path, **kwargs):
        """
        Generate data class from data dict

        :param json_path: str, path to JSON file
        :return: data class
        """
        with open(json_path) as fn:
            processing_data = json.load(fn)
        multi_data = processing_data if isinstance(list, processing_data) else [processing_data]
        return cls(backend_data=multi_data, **kwargs)


if __name__ == "__main__":
    front_coll = DBconnector(db_info.get("frontend")).get_collection(FRONT_COLL)
    back_dbc = DBconnector(db_info.get("backend"))
    back_coll = back_dbc.get_collection("computation")
    # _id = "28ef46bbcbb793fc0a8a4656e01aba3e"  # 502dc467c4780db94cc0c324a12c2b6b
    # d = back_coll.find_one({"_id": _id})
    # mol_id = d["mol_id"]
    # c_type = d["calculation_type"]
    # cond = d["data"]["conditions"]
    # chg = d["data"]["charge"]

    mol_id = "80MJUY"
    omega_q = \
        front_coll.find_one({"_id": mol_id}, {"mol_characterization.omega": 1}).get("mol_characterization", {}).get(
            "omega",
            [None])[
            0]
    cond = omega_q["conditions"]
    cond.update({"tuning_parameter": omega_q["value"]})
    g2c = Gaus2FrontCharacterization(
        _id=mol_id,
        conditions=cond,
    )
    g2c.insert_all_mol()

    # print(CV2Front(backend_data=ex_data, run_anodic=False, insert=False).meta_dict)
