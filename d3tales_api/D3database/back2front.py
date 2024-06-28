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
RMSD_DEFAULT = True


class Gaus2FrontCharacterization:
    """
    Update frontend db with backend db data for a particular set of data.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, _id, calculation_type, conditions, charge, data=None, insert=True, all_props=False,
                 rmsd=RMSD_DEFAULT):
        """

        :param _id: str, molecule ID
        :param calculation_type: str, calculation type
        :param conditions: dict, calculation conditions
        :param charge: int, charge
        :param data: dict, calculation data
        :param insert: bool, insert generated data to the frontend D3TaLES database if True
        :param all_props: bool, calculate all properties for the molecule if True
        """
        # connect to databases
        self.front_dbc = DBconnector(db_info.get("frontend"))
        self.front_coll = self.front_dbc.get_collection("base")
        self.back_dbc = DBconnector(db_info.get("backend"))
        self.back_coll = self.back_dbc.get_collection("computation")
        self.all_props = all_props

        # basic variables
        self.id = _id
        self.calculation_type = calculation_type
        self.conditions = copy.deepcopy(conditions)
        self.bare_conditions = copy.deepcopy(conditions)
        self.solvent = {"solvent": self.bare_conditions.pop("solvent")} if self.bare_conditions.get("solvent") else None
        self.charge = charge
        self.groundState_charge = self.front_coll.find_one({"_id": self.id})["mol_info"]["groundState_charge"]
        self.mol_data = list(self.back_coll.find({"mol_id": self.id}) or [])
        temp_data, self._hash = self.find_data(calc_type=calculation_type, solvent=self.solvent)
        self.data = data or temp_data
        if not self.data:
            raise ValueError("No data found for {} calculation {}.".format(_id, calculation_type))

        # all calculations for this molecule
        self.all_calcs = ['opt_groundState', 'freq_groundState', 'solv_energy_gsgs', 'tddft_groundState',
                          'opt_cation1', 'freq_cation1', 'solv_energy_c1c1', 'tddft_cation1',
                          'opt_cation2', 'freq_cation2', 'solv_energy_c2c2', 'tddft_cation2',
                          'opt_anion1', 'freq_anion1', 'solv_energy_a1a1', 'tddft_anion1',
                          'opt_anion2', 'freq_anion2', 'solv_energy_a2a2', 'tddft_anion2',
                          'energy_gsc1', 'energy_gsa1', 'energy_c1gs', 'energy_a1gs', 'energy_c1c2', 'energy_a1a2']

        # Get data for all calculations for this molecule
        for calc_type in self.all_calcs:
            if "solv_" in calc_type:
                setattr(self, calc_type, self.find_data(calc_type=calc_type, solvent=self.solvent or DEFAULT_SOLV))
            else:
                setattr(self, calc_type, self.find_data(calc_type=calc_type))

        # Set data for this calculation
        setattr(self, self.calculation_type, [self.data, self._hash])

        if insert:
            # generate data dictionary
            self.character_dict = {}
            properties = [
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
                properties.extend(
                    [self.rmsd_groundState_cation1, self.rmsd_cation1_cation2, self.rmsd_groundState_anion1,
                     self.rmsd_anion1_anion2])
            for prop in properties:
                try:
                    self.character_dict.update(prop())
                except ValueError as e:
                    print(e)
            FrontDB(schema_layer="species_characterization", instance=self.species_descriptors, _id=self.id)
            FrontDB(schema_layer="mol_characterization", instance=self.character_dict, _id=self.id)

    @property
    def species_descriptors(self):
        """Descriptor dict for species descriptors"""
        if "wtuning" in self.calculation_type:
            self.character_dict.update(
                self.return_descriptor_dict(self.data.get("omega"), name="omega", hashes=[self._hash]))
            return None
        species_dict = {
            -2: "anion2",
            -1: "anion1",
            0: "groundState",
            1: "cation1",
            2: "cation2",
        }
        species = species_dict[self.charge - self.groundState_charge]
        data_dict = {
            species: {
                "charge": self.charge,
                "spin_multiplicity": self.data["spin_multiplicity"],
                "is_groundState": self.data["is_groundState"],
            }
        }
        if "opt" in self.calculation_type or "energy" in self.calculation_type:
            gap = self.data["lumo"]["value"] - self.data["homo"]["value"]
            if self.data["lumo"]["unit"] == self.data["homo"]["unit"]:
                unit = self.data["lumo"]["unit"]
                homo_lumo_data = [{
                    "source_hash_ids": [self._hash],
                    "conditions": self.conditions,
                    "value": gap,
                    "unit": unit
                }]
                data_dict[species].update({"homo_lumo_gap": homo_lumo_data})
            dipole_moment = [{
                "source_hash_ids": [self._hash],
                "conditions": self.conditions,
                "value": np.linalg.norm(np.asarray(self.data["scf_dipole_moment"]["value"])),
                "unit": self.data["scf_dipole_moment"]["unit"]
            }]
            data_dict[species].update({"dipole_moment": dipole_moment})
            if "opt" in self.calculation_type:
                pmg_mol = Molecule.from_sites([Site.from_dict(sd) for sd in self.data["geometry"]])
                rdmol = pmgmol_to_rdmol(pmg_mol)[0]
                vol = ComputeMolVolume(rdmol)
                globular_volume = [{
                    "source_hash_ids": [self._hash],
                    "conditions": self.conditions,
                    "value": round(vol, 3),
                    "unit": "A^3"
                }]
                geom_data = [{
                    "source_hash_ids": [self._hash],
                    "conditions": self.conditions,
                    "sites": self.data["geometry"]
                }]
                data_dict[species].update({"globular_volume": globular_volume, "geometry": geom_data})
                for opt_descriptor in ["radical_buried_vol", "radical_spin", "radical_stability_score"]:
                    if self.data.get(opt_descriptor):
                        self.data[opt_descriptor].update(
                            {"source_hash_ids": [self._hash], "conditions": self.conditions})
                        data_dict[species].update({opt_descriptor: [self.data[opt_descriptor]]})
        if "tddft" in self.calculation_type:
            if self.data.get("dipole_moment"):
                dipole_moment = [{
                    "source_hash_ids": [self._hash],
                    "conditions": self.conditions,
                    "value": np.linalg.norm(np.asarray(self.data["scf_dipole_moment"]["value"])),
                    "unit": self.data["scf_dipole_moment"]["unit"]
                }]
                data_dict[species].update({"dipole_moment": dipole_moment})
            if self.data["excitations"].get("Singlet"):
                singlet_data = [{
                    "source_hash_ids": [self._hash],
                    "conditions": self.conditions,
                    "excitations": self.data["excitations"]["Singlet"]
                }]
                data_dict[species].update({"singlet_states": singlet_data})
            if self.data["excitations"].get("Triplet"):
                triplet_data = [{
                    "source_hash_ids": [self._hash],
                    "conditions": self.conditions,
                    "excitations": self.data["excitations"]["Triplet"]
                }]
                data_dict[species].update({"triplet_states": triplet_data})
            if self.data.get("singlet_plotting", ):
                plotting_data = [{
                    "source_hash_ids": [self._hash],
                    "conditions": self.conditions,
                    "plotting_data": self.data["singlet_plotting"]
                }]
                data_dict[species].update({"spectra": plotting_data})
        try:
            data_dict[species].update(self.solvation_energy(species))
        except ValueError:
            pass
        return data_dict

    @classmethod
    def from_data(cls, processing_data, **kwargs):
        """
        Generate data class from data dict

        :param processing_data: dict, data dict
        :return: data class
        """
        _id = processing_data.get("mol_id")
        calculation_type = processing_data.get("calculation_type")
        data = processing_data.get("data", {})
        conditions = data.get("conditions")
        charge = data.get("charge")
        return cls(_id=_id, calculation_type=calculation_type, conditions=conditions, charge=charge, data=data,
                   **kwargs)

    def return_descriptor_dict(self, value, unit="", hashes=None, name="", order=1, condition_addition=None, **kwargs):
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
        cond = copy.deepcopy(self.conditions)
        if condition_addition:
            cond.update(condition_addition)
        return {
            name:
                [{
                    "source_hash_ids": hashes or [],
                    "conditions": cond,
                    "order": order,
                    "value": value,
                    "unit": unit
                }]
        }

    def find_data(self, calc_type, solvent=None):
        """
        Find calculation data

        :param calc_type: str, calculation type
        :param solvent: str, solvent
        :return: [data dict, hash ID]
        """
        if solvent:
            data_conditions = copy.deepcopy(self.conditions)
            data_conditions["solvent"] = solvent.get("solvent", )
        else:
            data_conditions = self.bare_conditions
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
        if self.calculation_type not in calculations and not self.all_props:
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
        return self.return_descriptor_dict(solv_eng, unit='eV', hashes=h_ids, name="solvation_energy")

    def hole_reorganization_energy(self):
        """Descriptor dict for the hole reorganization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsc1", "opt_cation1", "energy_c1gs"])
        connector = {"gs_opt": "opt_groundState.scf_total_energy.value",
                     "ion_opt": "opt_cation1.scf_total_energy.value",
                     "gs_energy": "energy_c1gs.scf_total_energy.value",
                     "ion_energy": "energy_gsc1.scf_total_energy.value"}
        energy = ReorganizationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="hole_reorganization_energy")

    def electron_reorganization_energy(self):
        """Descriptor dict for the electron reorganization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsa1", "opt_anion1", "energy_a1gs"])
        connector = {"gs_opt": "opt_groundState.scf_total_energy.value",
                     "ion_opt": "opt_anion1.scf_total_energy.value",
                     "gs_energy": "energy_a1gs.scf_total_energy.value",
                     "ion_energy": "energy_gsa1.scf_total_energy.value"}
        energy = ReorganizationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="electron_reorganization_energy")

    def relaxation_groundState_cation1(self):
        """Descriptor dict for the relaxation energy of the ground state geometry to the +1 cation geometry"""
        h_ids, c_data = self.get_data(["energy_gsc1", "opt_cation1"])
        connector = {"opt_energy": "opt_cation1.scf_total_energy.value",
                     "energy": "energy_gsc1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_groundState_cation1")

    def relaxation_cation1_groundState(self):
        """Descriptor dict for the relaxation energy of the +1 cation geometry to the ground state geometry"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_c1gs"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_c1gs.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_cation1_groundState")

    def relaxation_groundState_anion1(self):
        """Descriptor dict for the relaxation energy of the ground state geometry to the -1 anion geometry"""
        h_ids, c_data = self.get_data(["energy_gsa1", "opt_anion1"])
        connector = {"opt_energy": "opt_anion1.scf_total_energy.value",
                     "energy": "energy_gsa1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_groundState_anion1")

    def relaxation_anion1_groundState(self):
        """Descriptor dict for the relaxation energy of the ground -1 anion to the ground state geometry"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_a1gs"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_a1gs.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="relaxation_anion1_groundState")

    def vertical_ionization_energy(self):
        """Descriptor dict for the vertical ionization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsc1"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_gsc1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_ionization_energy")

    def vertical_ionization_energy_2(self):
        """Descriptor dict for the vertical ionization energy"""
        h_ids, c_data = self.get_data(["opt_cation1", "energy_c1c2"])
        connector = {"opt_energy": "opt_cation1.scf_total_energy.value",
                     "energy": "energy_c1c2.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_ionization_energy_2")

    def vertical_electron_affinity(self):
        """Descriptor dict for the """
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsa1"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_gsa1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_electron_affinity")

    def vertical_electron_affinity_2(self):
        """Descriptor dict for the """
        h_ids, c_data = self.get_data(["opt_anion1", "energy_a1a2"])
        connector = {"opt_energy": "opt_anion1.scf_total_energy.value",
                     "energy": "energy_a1a2.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_electron_affinity_2")

    def adiabatic_ionization_energy(self):
        """Descriptor dict for the adiabatic ionization energy"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_cation1"])
        connector = {"energy_final": "opt_cation1.scf_total_energy.value",
                     "energy_initial": "opt_groundState.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_ionization_energy")

    def adiabatic_ionization_energy_2(self):
        """Descriptor dict for the second adiabatic ionization energy"""
        h_ids, c_data = self.get_data(["opt_cation1", "opt_cation2"])
        connector = {"energy_final": "opt_cation2.scf_total_energy.value",
                     "energy_initial": "opt_cation1.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_ionization_energy_2")

    def adiabatic_electron_affinity(self):
        """Descriptor dict for the adiabatic electron affinity"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_anion1"])
        connector = {"energy_final": "opt_anion1.scf_total_energy.value",
                     "energy_initial": "opt_groundState.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_electron_affinity")

    def adiabatic_electron_affinity_2(self):
        """Descriptor dict for the second adiabatic electron affinity"""
        h_ids, c_data = self.get_data(["opt_anion1", "opt_anion2"])
        connector = {"energy_final": "opt_anion2.scf_total_energy.value",
                     "energy_initial": "opt_anion1.scf_total_energy.value"}
        energy = EnergyDiffCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="adiabatic_electron_affinity_2")

    def rmsd_groundState_cation1(self):
        """Descriptor dict for the RMSD between the ground state and +1 cation geometries"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_cation1"])
        connector = {"geom_initial": "opt_groundState.geometry",
                     "geom_final": "opt_cation1.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_groundState_cation1")

    def rmsd_cation1_cation2(self):
        """Descriptor dict for the RMSD between the +1 cation and +2 cation geometries"""
        h_ids, c_data = self.get_data(["opt_cation2", "opt_cation1"])
        connector = {"geom_initial": "opt_cation1.geometry",
                     "geom_final": "opt_cation2.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_groundState_cation1")

    def rmsd_groundState_anion1(self):
        """Descriptor dict for the RMSD between the ground state and -1 anion geometries"""
        h_ids, c_data = self.get_data(["opt_groundState", "opt_anion1"])
        connector = {"geom_initial": "opt_groundState.geometry",
                     "geom_final": "opt_anion1.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_groundState_anion1")

    def rmsd_anion1_anion2(self):
        """Descriptor dict for the RMSD between the -1 anion and -2 anion geometries"""
        h_ids, c_data = self.get_data(["opt_anion1", "opt_anion2"])
        connector = {"geom_initial": "opt_anion1.geometry",
                     "geom_final": "opt_anion2.geometry"}
        dist = RMSDCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(dist, unit='A', hashes=h_ids, name="rmsd_groundState_anion1")

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
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="oxidation_potential",
                                           condition_addition={"reference_electrode": electrode})

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
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="reduction_potential",
                                           condition_addition={"reference_electrode": electrode})

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
        self.front_coll = DBconnector(db_info.get("frontend")).get_collection("base")
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
            instance = GenerateMolInfo(clean_smiles, origin_group=nlp_group, database="frontend").mol_info_dict
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

    def __init__(self, id_list=None, backend_data=None, metadata_dict=None, mol_id=None, e_half_scan_rate=0.1,
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
        self.front_coll = DBconnector(db_info.get("frontend")).get_collection("base")
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
            print("E 1/2s: ", ", ".join([str(e.get("value")) for e in self.e_halfs]))

        self.meta_dict.update({f"{self.redox_event}_potential": self.e_halfs})

        for i, e_half in enumerate(self.e_halfs):
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
            try:
                data_dict = i.get('data')
                data_dict["n"] = electron_num
                data_dict["current_cathodic"] = data_dict.get("forward", {})[electron_num - 1][1]
                data_dict["current_anodic"] = data_dict.get("reverse", {})[electron_num - 1][1]
                processed_data.append(data_dict)
            except:
                pass
        connector = {
            "n": "n",
            "i_p": "current_{}".format(curve_type),
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
        diffusion_coef = diffusion_cal.calculate(processed_data, sci_notation=True)

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

            "e_half": "data.e_half.0",
            "e_rev": "data.e_ref",
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
        return [self.prop_dict(e_half.get("value"), unit=e_half.get("unit"), order=i + 1, hashes=[self.e_halfs_id],
                               conditions=self.e_halfs_conditions) for i, e_half in enumerate(e_halfs)]

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
    back_dbc = DBconnector(db_info.get("backend"))
    back_coll = back_dbc.get_collection("computation")
    _id = "28ef46bbcbb793fc0a8a4656e01aba3e"  # 502dc467c4780db94cc0c324a12c2b6b
    data = back_coll.find_one({"_id": _id})
    mol_id = data["mol_id"]
    calculation_type = data["calculation_type"]
    conditions = data["data"]["conditions"]
    charge = data["data"]["charge"]
    ex_data = [
        {
            "_id": "796c0b80-9598-477e-bed8-787de66377f4",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:44.181483",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv01_12_20_53.bin",
                "header": "CV cycle01_cv01_12_20_53",
                "note": "",
                "date_recorded": "2024-01-31T12:21:02",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.5,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            3.86e-08
                        ],
                        [
                            0.45,
                            5.488e-08
                        ],
                        [
                            0.46,
                            1.128e-08
                        ],
                        [
                            0.47,
                            -8.676e-09
                        ],
                        [
                            0.48,
                            2.756e-08
                        ],
                        [
                            0.49,
                            5.673e-08
                        ],
                        [
                            0.5,
                            6.838e-08
                        ],
                        [
                            0.51,
                            2.909e-08
                        ],
                        [
                            0.52,
                            1.189e-08
                        ],
                        [
                            0.53,
                            4.352e-08
                        ],
                        [
                            0.54,
                            8.282e-08
                        ],
                        [
                            0.55,
                            1.003e-07
                        ],
                        [
                            0.56,
                            7.145e-08
                        ],
                        [
                            0.57,
                            8.098e-08
                        ],
                        [
                            0.58,
                            1.442e-07
                        ],
                        [
                            0.59,
                            2.201e-07
                        ],
                        [
                            0.6,
                            3.014e-07
                        ],
                        [
                            0.61,
                            3.635e-07
                        ],
                        [
                            0.62,
                            5.096e-07
                        ],
                        [
                            0.63,
                            7.58e-07
                        ],
                        [
                            0.64,
                            1.06e-06
                        ],
                        [
                            0.65,
                            1.413e-06
                        ],
                        [
                            0.66,
                            1.89e-06
                        ],
                        [
                            0.67,
                            2.514e-06
                        ],
                        [
                            0.68,
                            3.202e-06
                        ],
                        [
                            0.69,
                            3.843e-06
                        ],
                        [
                            0.7,
                            4.344e-06
                        ],
                        [
                            0.71,
                            4.81e-06
                        ],
                        [
                            0.72,
                            5.08e-06
                        ],
                        [
                            0.73,
                            5.1e-06
                        ],
                        [
                            0.74,
                            4.941e-06
                        ],
                        [
                            0.75,
                            4.672e-06
                        ],
                        [
                            0.76,
                            4.376e-06
                        ],
                        [
                            0.77,
                            4.07e-06
                        ],
                        [
                            0.78,
                            3.801e-06
                        ],
                        [
                            0.79,
                            3.56e-06
                        ],
                        [
                            0.8,
                            3.333e-06
                        ],
                        [
                            0.81,
                            3.11e-06
                        ],
                        [
                            0.82,
                            2.926e-06
                        ],
                        [
                            0.83,
                            2.819e-06
                        ],
                        [
                            0.84,
                            2.73e-06
                        ],
                        [
                            0.85,
                            2.624e-06
                        ],
                        [
                            0.86,
                            2.508e-06
                        ],
                        [
                            0.87,
                            2.418e-06
                        ],
                        [
                            0.88,
                            2.385e-06
                        ],
                        [
                            0.89,
                            2.363e-06
                        ],
                        [
                            0.9,
                            2.311e-06
                        ],
                        [
                            0.91,
                            2.234e-06
                        ],
                        [
                            0.92,
                            2.18e-06
                        ],
                        [
                            0.93,
                            2.181e-06
                        ],
                        [
                            0.94,
                            2.176e-06
                        ],
                        [
                            0.95,
                            2.146e-06
                        ],
                        [
                            0.96,
                            2.086e-06
                        ],
                        [
                            0.97,
                            2.038e-06
                        ],
                        [
                            0.98,
                            2.055e-06
                        ],
                        [
                            0.98,
                            1.993e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.947e-06
                        ],
                        [
                            0.96,
                            1.893e-06
                        ],
                        [
                            0.95,
                            1.85e-06
                        ],
                        [
                            0.94,
                            1.868e-06
                        ],
                        [
                            0.93,
                            1.883e-06
                        ],
                        [
                            0.92,
                            1.858e-06
                        ],
                        [
                            0.91,
                            1.812e-06
                        ],
                        [
                            0.9,
                            1.782e-06
                        ],
                        [
                            0.89,
                            1.805e-06
                        ],
                        [
                            0.88,
                            1.819e-06
                        ],
                        [
                            0.87,
                            1.799e-06
                        ],
                        [
                            0.86,
                            1.748e-06
                        ],
                        [
                            0.85,
                            1.71e-06
                        ],
                        [
                            0.84,
                            1.728e-06
                        ],
                        [
                            0.83,
                            1.715e-06
                        ],
                        [
                            0.82,
                            1.659e-06
                        ],
                        [
                            0.81,
                            1.566e-06
                        ],
                        [
                            0.8,
                            1.461e-06
                        ],
                        [
                            0.79,
                            1.377e-06
                        ],
                        [
                            0.78,
                            1.214e-06
                        ],
                        [
                            0.77,
                            9.327e-07
                        ],
                        [
                            0.76,
                            5.861e-07
                        ],
                        [
                            0.75,
                            1.648e-07
                        ],
                        [
                            0.74,
                            -3.504e-07
                        ],
                        [
                            0.73,
                            -1.014e-06
                        ],
                        [
                            0.72,
                            -1.771e-06
                        ],
                        [
                            0.71,
                            -2.389e-06
                        ],
                        [
                            0.7,
                            -2.875e-06
                        ],
                        [
                            0.69,
                            -3.225e-06
                        ],
                        [
                            0.68,
                            -3.411e-06
                        ],
                        [
                            0.67,
                            -3.382e-06
                        ],
                        [
                            0.66,
                            -3.154e-06
                        ],
                        [
                            0.65,
                            -2.879e-06
                        ],
                        [
                            0.64,
                            -2.549e-06
                        ],
                        [
                            0.63,
                            -2.242e-06
                        ],
                        [
                            0.62,
                            -1.978e-06
                        ],
                        [
                            0.61,
                            -1.741e-06
                        ],
                        [
                            0.6,
                            -1.541e-06
                        ],
                        [
                            0.59,
                            -1.318e-06
                        ],
                        [
                            0.58,
                            -1.138e-06
                        ],
                        [
                            0.57,
                            -1.019e-06
                        ],
                        [
                            0.56,
                            -9.375e-07
                        ],
                        [
                            0.55,
                            -8.515e-07
                        ],
                        [
                            0.54,
                            -7.29e-07
                        ],
                        [
                            0.53,
                            -6.363e-07
                        ],
                        [
                            0.52,
                            -5.896e-07
                        ],
                        [
                            0.51,
                            -5.699e-07
                        ],
                        [
                            0.5,
                            -5.402e-07
                        ],
                        [
                            0.49,
                            -4.557e-07
                        ],
                        [
                            0.48,
                            -3.998e-07
                        ],
                        [
                            0.47,
                            -3.885e-07
                        ],
                        [
                            0.46,
                            -3.891e-07
                        ],
                        [
                            0.45,
                            -3.808e-07
                        ],
                        [
                            0.44,
                            -3.145e-07
                        ],
                        [
                            0.43,
                            -2.7e-07
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.73,
                        5.1e-06
                    ]
                ],
                "reverse": [
                    [
                        0.68,
                        -3.411e-06
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            3.86e-08,
                            5.488e-08,
                            1.128e-08,
                            -8.676e-09,
                            2.756e-08,
                            5.673e-08,
                            6.838e-08,
                            2.909e-08,
                            1.189e-08,
                            4.352e-08,
                            8.282e-08,
                            1.003e-07,
                            7.145e-08,
                            8.098e-08,
                            1.442e-07,
                            2.201e-07,
                            3.014e-07,
                            3.635e-07,
                            5.096e-07,
                            7.58e-07,
                            1.06e-06,
                            1.413e-06,
                            1.89e-06,
                            2.514e-06,
                            3.202e-06,
                            3.843e-06,
                            4.344e-06,
                            4.81e-06,
                            5.08e-06,
                            5.1e-06,
                            4.941e-06,
                            4.672e-06,
                            4.376e-06,
                            4.07e-06,
                            3.801e-06,
                            3.56e-06,
                            3.333e-06,
                            3.11e-06,
                            2.926e-06,
                            2.819e-06,
                            2.73e-06,
                            2.624e-06,
                            2.508e-06,
                            2.418e-06,
                            2.385e-06,
                            2.363e-06,
                            2.311e-06,
                            2.234e-06,
                            2.18e-06,
                            2.181e-06,
                            2.176e-06,
                            2.146e-06,
                            2.086e-06,
                            2.038e-06,
                            2.055e-06,
                            1.993e-06,
                            1.947e-06,
                            1.893e-06,
                            1.85e-06,
                            1.868e-06,
                            1.883e-06,
                            1.858e-06,
                            1.812e-06,
                            1.782e-06,
                            1.805e-06,
                            1.819e-06,
                            1.799e-06,
                            1.748e-06,
                            1.71e-06,
                            1.728e-06,
                            1.715e-06,
                            1.659e-06,
                            1.566e-06,
                            1.461e-06,
                            1.377e-06,
                            1.214e-06,
                            9.327e-07,
                            5.861e-07,
                            1.648e-07,
                            -3.504e-07,
                            -1.014e-06,
                            -1.771e-06,
                            -2.389e-06,
                            -2.875e-06,
                            -3.225e-06,
                            -3.411e-06,
                            -3.382e-06,
                            -3.154e-06,
                            -2.879e-06,
                            -2.549e-06,
                            -2.242e-06,
                            -1.978e-06,
                            -1.741e-06,
                            -1.541e-06,
                            -1.318e-06,
                            -1.138e-06,
                            -1.019e-06,
                            -9.375e-07,
                            -8.515e-07,
                            -7.29e-07,
                            -6.363e-07,
                            -5.896e-07,
                            -5.699e-07,
                            -5.402e-07,
                            -4.557e-07,
                            -3.998e-07,
                            -3.885e-07,
                            -3.891e-07,
                            -3.808e-07,
                            -3.145e-07,
                            -2.7e-07
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "reversible"
                ],
                "e_half": [
                    0.705
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            3.86e-08
                        ],
                        [
                            0.45,
                            5.488e-08
                        ],
                        [
                            0.46,
                            1.128e-08
                        ],
                        [
                            0.47,
                            -8.676e-09
                        ],
                        [
                            0.48,
                            2.756e-08
                        ],
                        [
                            0.49,
                            5.673e-08
                        ],
                        [
                            0.5,
                            6.838e-08
                        ],
                        [
                            0.51,
                            2.909e-08
                        ],
                        [
                            0.52,
                            1.189e-08
                        ],
                        [
                            0.53,
                            4.352e-08
                        ],
                        [
                            0.54,
                            8.282e-08
                        ],
                        [
                            0.55,
                            1.003e-07
                        ],
                        [
                            0.56,
                            7.145e-08
                        ],
                        [
                            0.57,
                            8.098e-08
                        ],
                        [
                            0.58,
                            1.442e-07
                        ],
                        [
                            0.59,
                            2.201e-07
                        ],
                        [
                            0.6,
                            3.014e-07
                        ],
                        [
                            0.61,
                            3.635e-07
                        ],
                        [
                            0.62,
                            5.096e-07
                        ],
                        [
                            0.63,
                            7.58e-07
                        ],
                        [
                            0.64,
                            1.06e-06
                        ],
                        [
                            0.65,
                            1.413e-06
                        ],
                        [
                            0.66,
                            1.89e-06
                        ],
                        [
                            0.67,
                            2.514e-06
                        ],
                        [
                            0.68,
                            3.202e-06
                        ],
                        [
                            0.69,
                            3.843e-06
                        ],
                        [
                            0.7,
                            4.344e-06
                        ],
                        [
                            0.71,
                            4.81e-06
                        ],
                        [
                            0.72,
                            5.08e-06
                        ],
                        [
                            0.73,
                            5.1e-06
                        ],
                        [
                            0.74,
                            4.941e-06
                        ],
                        [
                            0.75,
                            4.672e-06
                        ],
                        [
                            0.76,
                            4.376e-06
                        ],
                        [
                            0.77,
                            4.07e-06
                        ],
                        [
                            0.78,
                            3.801e-06
                        ],
                        [
                            0.79,
                            3.56e-06
                        ],
                        [
                            0.8,
                            3.333e-06
                        ],
                        [
                            0.81,
                            3.11e-06
                        ],
                        [
                            0.82,
                            2.926e-06
                        ],
                        [
                            0.83,
                            2.819e-06
                        ],
                        [
                            0.84,
                            2.73e-06
                        ],
                        [
                            0.85,
                            2.624e-06
                        ],
                        [
                            0.86,
                            2.508e-06
                        ],
                        [
                            0.87,
                            2.418e-06
                        ],
                        [
                            0.88,
                            2.385e-06
                        ],
                        [
                            0.89,
                            2.363e-06
                        ],
                        [
                            0.9,
                            2.311e-06
                        ],
                        [
                            0.91,
                            2.234e-06
                        ],
                        [
                            0.92,
                            2.18e-06
                        ],
                        [
                            0.93,
                            2.181e-06
                        ],
                        [
                            0.94,
                            2.176e-06
                        ],
                        [
                            0.95,
                            2.146e-06
                        ],
                        [
                            0.96,
                            2.086e-06
                        ],
                        [
                            0.97,
                            2.038e-06
                        ],
                        [
                            0.98,
                            2.055e-06
                        ],
                        [
                            0.98,
                            1.993e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.947e-06
                        ],
                        [
                            0.96,
                            1.893e-06
                        ],
                        [
                            0.95,
                            1.85e-06
                        ],
                        [
                            0.94,
                            1.868e-06
                        ],
                        [
                            0.93,
                            1.883e-06
                        ],
                        [
                            0.92,
                            1.858e-06
                        ],
                        [
                            0.91,
                            1.812e-06
                        ],
                        [
                            0.9,
                            1.782e-06
                        ],
                        [
                            0.89,
                            1.805e-06
                        ],
                        [
                            0.88,
                            1.819e-06
                        ],
                        [
                            0.87,
                            1.799e-06
                        ],
                        [
                            0.86,
                            1.748e-06
                        ],
                        [
                            0.85,
                            1.71e-06
                        ],
                        [
                            0.84,
                            1.728e-06
                        ],
                        [
                            0.83,
                            1.715e-06
                        ],
                        [
                            0.82,
                            1.659e-06
                        ],
                        [
                            0.81,
                            1.566e-06
                        ],
                        [
                            0.8,
                            1.461e-06
                        ],
                        [
                            0.79,
                            1.377e-06
                        ],
                        [
                            0.78,
                            1.214e-06
                        ],
                        [
                            0.77,
                            9.327e-07
                        ],
                        [
                            0.76,
                            5.861e-07
                        ],
                        [
                            0.75,
                            1.648e-07
                        ],
                        [
                            0.74,
                            -3.504e-07
                        ],
                        [
                            0.73,
                            -1.014e-06
                        ],
                        [
                            0.72,
                            -1.771e-06
                        ],
                        [
                            0.71,
                            -2.389e-06
                        ],
                        [
                            0.7,
                            -2.875e-06
                        ],
                        [
                            0.69,
                            -3.225e-06
                        ],
                        [
                            0.68,
                            -3.411e-06
                        ],
                        [
                            0.67,
                            -3.382e-06
                        ],
                        [
                            0.66,
                            -3.154e-06
                        ],
                        [
                            0.65,
                            -2.879e-06
                        ],
                        [
                            0.64,
                            -2.549e-06
                        ],
                        [
                            0.63,
                            -2.242e-06
                        ],
                        [
                            0.62,
                            -1.978e-06
                        ],
                        [
                            0.61,
                            -1.741e-06
                        ],
                        [
                            0.6,
                            -1.541e-06
                        ],
                        [
                            0.59,
                            -1.318e-06
                        ],
                        [
                            0.58,
                            -1.138e-06
                        ],
                        [
                            0.57,
                            -1.019e-06
                        ],
                        [
                            0.56,
                            -9.375e-07
                        ],
                        [
                            0.55,
                            -8.515e-07
                        ],
                        [
                            0.54,
                            -7.29e-07
                        ],
                        [
                            0.53,
                            -6.363e-07
                        ],
                        [
                            0.52,
                            -5.896e-07
                        ],
                        [
                            0.51,
                            -5.699e-07
                        ],
                        [
                            0.5,
                            -5.402e-07
                        ],
                        [
                            0.49,
                            -4.557e-07
                        ],
                        [
                            0.48,
                            -3.998e-07
                        ],
                        [
                            0.47,
                            -3.885e-07
                        ],
                        [
                            0.46,
                            -3.891e-07
                        ],
                        [
                            0.45,
                            -3.808e-07
                        ],
                        [
                            0.44,
                            -3.145e-07
                        ],
                        [
                            0.43,
                            -2.7e-07
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 5.1e-06,
                "current_anodic": -3.411e-06,
                "diffusion": 1.451e-09
            }
        },
        {
            "_id": "a67777ff-0447-4df1-879c-bd7916ce8c22",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:46.066759",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv02_12_21_08.bin",
                "header": "CV cycle01_cv02_12_21_08",
                "note": "",
                "date_recorded": "2024-01-31T12:21:18",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.4,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            -2.495e-08
                        ],
                        [
                            0.45,
                            4.106e-08
                        ],
                        [
                            0.46,
                            -2.003e-08
                        ],
                        [
                            0.47,
                            4.291e-08
                        ],
                        [
                            0.48,
                            -1.513e-08
                        ],
                        [
                            0.49,
                            4.751e-08
                        ],
                        [
                            0.5,
                            -1.113e-08
                        ],
                        [
                            0.51,
                            5.212e-08
                        ],
                        [
                            0.52,
                            -5.918e-09
                        ],
                        [
                            0.53,
                            5.948e-08
                        ],
                        [
                            0.54,
                            4.217e-09
                        ],
                        [
                            0.55,
                            7.452e-08
                        ],
                        [
                            0.56,
                            2.663e-08
                        ],
                        [
                            0.57,
                            1.086e-07
                        ],
                        [
                            0.58,
                            7.33e-08
                        ],
                        [
                            0.59,
                            1.789e-07
                        ],
                        [
                            0.6,
                            1.722e-07
                        ],
                        [
                            0.61,
                            3.26e-07
                        ],
                        [
                            0.62,
                            3.797e-07
                        ],
                        [
                            0.63,
                            6.204e-07
                        ],
                        [
                            0.64,
                            7.906e-07
                        ],
                        [
                            0.65,
                            1.178e-06
                        ],
                        [
                            0.66,
                            1.527e-06
                        ],
                        [
                            0.67,
                            2.105e-06
                        ],
                        [
                            0.68,
                            2.607e-06
                        ],
                        [
                            0.69,
                            3.268e-06
                        ],
                        [
                            0.7,
                            3.728e-06
                        ],
                        [
                            0.71,
                            4.207e-06
                        ],
                        [
                            0.72,
                            4.377e-06
                        ],
                        [
                            0.73,
                            4.496e-06
                        ],
                        [
                            0.74,
                            4.334e-06
                        ],
                        [
                            0.75,
                            4.199e-06
                        ],
                        [
                            0.76,
                            3.882e-06
                        ],
                        [
                            0.77,
                            3.685e-06
                        ],
                        [
                            0.78,
                            3.368e-06
                        ],
                        [
                            0.79,
                            3.206e-06
                        ],
                        [
                            0.8,
                            2.943e-06
                        ],
                        [
                            0.81,
                            2.838e-06
                        ],
                        [
                            0.82,
                            2.63e-06
                        ],
                        [
                            0.83,
                            2.578e-06
                        ],
                        [
                            0.84,
                            2.409e-06
                        ],
                        [
                            0.85,
                            2.397e-06
                        ],
                        [
                            0.86,
                            2.255e-06
                        ],
                        [
                            0.87,
                            2.27e-06
                        ],
                        [
                            0.88,
                            2.144e-06
                        ],
                        [
                            0.89,
                            2.18e-06
                        ],
                        [
                            0.9,
                            2.066e-06
                        ],
                        [
                            0.91,
                            2.113e-06
                        ],
                        [
                            0.92,
                            2.007e-06
                        ],
                        [
                            0.93,
                            2.063e-06
                        ],
                        [
                            0.94,
                            1.963e-06
                        ],
                        [
                            0.95,
                            2.022e-06
                        ],
                        [
                            0.96,
                            1.925e-06
                        ],
                        [
                            0.97,
                            1.988e-06
                        ],
                        [
                            0.98,
                            1.894e-06
                        ],
                        [
                            0.98,
                            1.921e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.821e-06
                        ],
                        [
                            0.96,
                            1.878e-06
                        ],
                        [
                            0.95,
                            1.789e-06
                        ],
                        [
                            0.94,
                            1.853e-06
                        ],
                        [
                            0.93,
                            1.765e-06
                        ],
                        [
                            0.92,
                            1.832e-06
                        ],
                        [
                            0.91,
                            1.747e-06
                        ],
                        [
                            0.9,
                            1.815e-06
                        ],
                        [
                            0.89,
                            1.729e-06
                        ],
                        [
                            0.88,
                            1.798e-06
                        ],
                        [
                            0.87,
                            1.71e-06
                        ],
                        [
                            0.86,
                            1.776e-06
                        ],
                        [
                            0.85,
                            1.681e-06
                        ],
                        [
                            0.84,
                            1.739e-06
                        ],
                        [
                            0.83,
                            1.632e-06
                        ],
                        [
                            0.82,
                            1.665e-06
                        ],
                        [
                            0.81,
                            1.533e-06
                        ],
                        [
                            0.8,
                            1.52e-06
                        ],
                        [
                            0.79,
                            1.332e-06
                        ],
                        [
                            0.78,
                            1.227e-06
                        ],
                        [
                            0.77,
                            9.345e-07
                        ],
                        [
                            0.76,
                            6.711e-07
                        ],
                        [
                            0.75,
                            2.173e-07
                        ],
                        [
                            0.74,
                            -2.528e-07
                        ],
                        [
                            0.73,
                            -8.641e-07
                        ],
                        [
                            0.72,
                            -1.457e-06
                        ],
                        [
                            0.71,
                            -2.04e-06
                        ],
                        [
                            0.7,
                            -2.474e-06
                        ],
                        [
                            0.69,
                            -2.764e-06
                        ],
                        [
                            0.68,
                            -2.839e-06
                        ],
                        [
                            0.67,
                            -2.791e-06
                        ],
                        [
                            0.66,
                            -2.565e-06
                        ],
                        [
                            0.65,
                            -2.357e-06
                        ],
                        [
                            0.64,
                            -2.034e-06
                        ],
                        [
                            0.63,
                            -1.83e-06
                        ],
                        [
                            0.62,
                            -1.523e-06
                        ],
                        [
                            0.61,
                            -1.381e-06
                        ],
                        [
                            0.6,
                            -1.13e-06
                        ],
                        [
                            0.59,
                            -1.048e-06
                        ],
                        [
                            0.58,
                            -8.432e-07
                        ],
                        [
                            0.57,
                            -8.097e-07
                        ],
                        [
                            0.56,
                            -6.393e-07
                        ],
                        [
                            0.55,
                            -6.415e-07
                        ],
                        [
                            0.54,
                            -4.96e-07
                        ],
                        [
                            0.53,
                            -5.214e-07
                        ],
                        [
                            0.52,
                            -3.934e-07
                        ],
                        [
                            0.51,
                            -4.361e-07
                        ],
                        [
                            0.5,
                            -3.185e-07
                        ],
                        [
                            0.49,
                            -3.716e-07
                        ],
                        [
                            0.48,
                            -2.638e-07
                        ],
                        [
                            0.47,
                            -3.249e-07
                        ],
                        [
                            0.46,
                            -2.224e-07
                        ],
                        [
                            0.45,
                            -2.856e-07
                        ],
                        [
                            0.44,
                            -1.895e-07
                        ],
                        [
                            0.43,
                            -2.546e-07
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.73,
                        4.496e-06
                    ]
                ],
                "reverse": [
                    [
                        0.68,
                        -2.839e-06
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            -2.495e-08,
                            4.106e-08,
                            -2.003e-08,
                            4.291e-08,
                            -1.513e-08,
                            4.751e-08,
                            -1.113e-08,
                            5.212e-08,
                            -5.918e-09,
                            5.948e-08,
                            4.217e-09,
                            7.452e-08,
                            2.663e-08,
                            1.086e-07,
                            7.33e-08,
                            1.789e-07,
                            1.722e-07,
                            3.26e-07,
                            3.797e-07,
                            6.204e-07,
                            7.906e-07,
                            1.178e-06,
                            1.527e-06,
                            2.105e-06,
                            2.607e-06,
                            3.268e-06,
                            3.728e-06,
                            4.207e-06,
                            4.377e-06,
                            4.496e-06,
                            4.334e-06,
                            4.199e-06,
                            3.882e-06,
                            3.685e-06,
                            3.368e-06,
                            3.206e-06,
                            2.943e-06,
                            2.838e-06,
                            2.63e-06,
                            2.578e-06,
                            2.409e-06,
                            2.397e-06,
                            2.255e-06,
                            2.27e-06,
                            2.144e-06,
                            2.18e-06,
                            2.066e-06,
                            2.113e-06,
                            2.007e-06,
                            2.063e-06,
                            1.963e-06,
                            2.022e-06,
                            1.925e-06,
                            1.988e-06,
                            1.894e-06,
                            1.921e-06,
                            1.821e-06,
                            1.878e-06,
                            1.789e-06,
                            1.853e-06,
                            1.765e-06,
                            1.832e-06,
                            1.747e-06,
                            1.815e-06,
                            1.729e-06,
                            1.798e-06,
                            1.71e-06,
                            1.776e-06,
                            1.681e-06,
                            1.739e-06,
                            1.632e-06,
                            1.665e-06,
                            1.533e-06,
                            1.52e-06,
                            1.332e-06,
                            1.227e-06,
                            9.345e-07,
                            6.711e-07,
                            2.173e-07,
                            -2.528e-07,
                            -8.641e-07,
                            -1.457e-06,
                            -2.04e-06,
                            -2.474e-06,
                            -2.764e-06,
                            -2.839e-06,
                            -2.791e-06,
                            -2.565e-06,
                            -2.357e-06,
                            -2.034e-06,
                            -1.83e-06,
                            -1.523e-06,
                            -1.381e-06,
                            -1.13e-06,
                            -1.048e-06,
                            -8.432e-07,
                            -8.097e-07,
                            -6.393e-07,
                            -6.415e-07,
                            -4.96e-07,
                            -5.214e-07,
                            -3.934e-07,
                            -4.361e-07,
                            -3.185e-07,
                            -3.716e-07,
                            -2.638e-07,
                            -3.249e-07,
                            -2.224e-07,
                            -2.856e-07,
                            -1.895e-07,
                            -2.546e-07
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "reversible"
                ],
                "e_half": [
                    0.705
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            -2.495e-08
                        ],
                        [
                            0.45,
                            4.106e-08
                        ],
                        [
                            0.46,
                            -2.003e-08
                        ],
                        [
                            0.47,
                            4.291e-08
                        ],
                        [
                            0.48,
                            -1.513e-08
                        ],
                        [
                            0.49,
                            4.751e-08
                        ],
                        [
                            0.5,
                            -1.113e-08
                        ],
                        [
                            0.51,
                            5.212e-08
                        ],
                        [
                            0.52,
                            -5.918e-09
                        ],
                        [
                            0.53,
                            5.948e-08
                        ],
                        [
                            0.54,
                            4.217e-09
                        ],
                        [
                            0.55,
                            7.452e-08
                        ],
                        [
                            0.56,
                            2.663e-08
                        ],
                        [
                            0.57,
                            1.086e-07
                        ],
                        [
                            0.58,
                            7.33e-08
                        ],
                        [
                            0.59,
                            1.789e-07
                        ],
                        [
                            0.6,
                            1.722e-07
                        ],
                        [
                            0.61,
                            3.26e-07
                        ],
                        [
                            0.62,
                            3.797e-07
                        ],
                        [
                            0.63,
                            6.204e-07
                        ],
                        [
                            0.64,
                            7.906e-07
                        ],
                        [
                            0.65,
                            1.178e-06
                        ],
                        [
                            0.66,
                            1.527e-06
                        ],
                        [
                            0.67,
                            2.105e-06
                        ],
                        [
                            0.68,
                            2.607e-06
                        ],
                        [
                            0.69,
                            3.268e-06
                        ],
                        [
                            0.7,
                            3.728e-06
                        ],
                        [
                            0.71,
                            4.207e-06
                        ],
                        [
                            0.72,
                            4.377e-06
                        ],
                        [
                            0.73,
                            4.496e-06
                        ],
                        [
                            0.74,
                            4.334e-06
                        ],
                        [
                            0.75,
                            4.199e-06
                        ],
                        [
                            0.76,
                            3.882e-06
                        ],
                        [
                            0.77,
                            3.685e-06
                        ],
                        [
                            0.78,
                            3.368e-06
                        ],
                        [
                            0.79,
                            3.206e-06
                        ],
                        [
                            0.8,
                            2.943e-06
                        ],
                        [
                            0.81,
                            2.838e-06
                        ],
                        [
                            0.82,
                            2.63e-06
                        ],
                        [
                            0.83,
                            2.578e-06
                        ],
                        [
                            0.84,
                            2.409e-06
                        ],
                        [
                            0.85,
                            2.397e-06
                        ],
                        [
                            0.86,
                            2.255e-06
                        ],
                        [
                            0.87,
                            2.27e-06
                        ],
                        [
                            0.88,
                            2.144e-06
                        ],
                        [
                            0.89,
                            2.18e-06
                        ],
                        [
                            0.9,
                            2.066e-06
                        ],
                        [
                            0.91,
                            2.113e-06
                        ],
                        [
                            0.92,
                            2.007e-06
                        ],
                        [
                            0.93,
                            2.063e-06
                        ],
                        [
                            0.94,
                            1.963e-06
                        ],
                        [
                            0.95,
                            2.022e-06
                        ],
                        [
                            0.96,
                            1.925e-06
                        ],
                        [
                            0.97,
                            1.988e-06
                        ],
                        [
                            0.98,
                            1.894e-06
                        ],
                        [
                            0.98,
                            1.921e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.821e-06
                        ],
                        [
                            0.96,
                            1.878e-06
                        ],
                        [
                            0.95,
                            1.789e-06
                        ],
                        [
                            0.94,
                            1.853e-06
                        ],
                        [
                            0.93,
                            1.765e-06
                        ],
                        [
                            0.92,
                            1.832e-06
                        ],
                        [
                            0.91,
                            1.747e-06
                        ],
                        [
                            0.9,
                            1.815e-06
                        ],
                        [
                            0.89,
                            1.729e-06
                        ],
                        [
                            0.88,
                            1.798e-06
                        ],
                        [
                            0.87,
                            1.71e-06
                        ],
                        [
                            0.86,
                            1.776e-06
                        ],
                        [
                            0.85,
                            1.681e-06
                        ],
                        [
                            0.84,
                            1.739e-06
                        ],
                        [
                            0.83,
                            1.632e-06
                        ],
                        [
                            0.82,
                            1.665e-06
                        ],
                        [
                            0.81,
                            1.533e-06
                        ],
                        [
                            0.8,
                            1.52e-06
                        ],
                        [
                            0.79,
                            1.332e-06
                        ],
                        [
                            0.78,
                            1.227e-06
                        ],
                        [
                            0.77,
                            9.345e-07
                        ],
                        [
                            0.76,
                            6.711e-07
                        ],
                        [
                            0.75,
                            2.173e-07
                        ],
                        [
                            0.74,
                            -2.528e-07
                        ],
                        [
                            0.73,
                            -8.641e-07
                        ],
                        [
                            0.72,
                            -1.457e-06
                        ],
                        [
                            0.71,
                            -2.04e-06
                        ],
                        [
                            0.7,
                            -2.474e-06
                        ],
                        [
                            0.69,
                            -2.764e-06
                        ],
                        [
                            0.68,
                            -2.839e-06
                        ],
                        [
                            0.67,
                            -2.791e-06
                        ],
                        [
                            0.66,
                            -2.565e-06
                        ],
                        [
                            0.65,
                            -2.357e-06
                        ],
                        [
                            0.64,
                            -2.034e-06
                        ],
                        [
                            0.63,
                            -1.83e-06
                        ],
                        [
                            0.62,
                            -1.523e-06
                        ],
                        [
                            0.61,
                            -1.381e-06
                        ],
                        [
                            0.6,
                            -1.13e-06
                        ],
                        [
                            0.59,
                            -1.048e-06
                        ],
                        [
                            0.58,
                            -8.432e-07
                        ],
                        [
                            0.57,
                            -8.097e-07
                        ],
                        [
                            0.56,
                            -6.393e-07
                        ],
                        [
                            0.55,
                            -6.415e-07
                        ],
                        [
                            0.54,
                            -4.96e-07
                        ],
                        [
                            0.53,
                            -5.214e-07
                        ],
                        [
                            0.52,
                            -3.934e-07
                        ],
                        [
                            0.51,
                            -4.361e-07
                        ],
                        [
                            0.5,
                            -3.185e-07
                        ],
                        [
                            0.49,
                            -3.716e-07
                        ],
                        [
                            0.48,
                            -2.638e-07
                        ],
                        [
                            0.47,
                            -3.249e-07
                        ],
                        [
                            0.46,
                            -2.224e-07
                        ],
                        [
                            0.45,
                            -2.856e-07
                        ],
                        [
                            0.44,
                            -1.895e-07
                        ],
                        [
                            0.43,
                            -2.546e-07
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 4.496e-06,
                "current_anodic": -2.839e-06,
                "diffusion": 1.451e-09
            }
        },
        {
            "_id": "7ac8b1bf-4dfe-4af8-afe8-041814d3424e",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:48.007120",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv03_12_21_23.bin",
                "header": "CV cycle01_cv03_12_21_23",
                "note": "",
                "date_recorded": "2024-01-31T12:21:33",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.3,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            -2.236e-09
                        ],
                        [
                            0.45,
                            -1e-09
                        ],
                        [
                            0.46,
                            8.348e-10
                        ],
                        [
                            0.47,
                            2.37e-09
                        ],
                        [
                            0.48,
                            5.752e-09
                        ],
                        [
                            0.49,
                            6.976e-09
                        ],
                        [
                            0.5,
                            1.036e-08
                        ],
                        [
                            0.51,
                            1.282e-08
                        ],
                        [
                            0.52,
                            1.558e-08
                        ],
                        [
                            0.53,
                            2.049e-08
                        ],
                        [
                            0.54,
                            2.571e-08
                        ],
                        [
                            0.55,
                            3.216e-08
                        ],
                        [
                            0.56,
                            4.321e-08
                        ],
                        [
                            0.57,
                            5.764e-08
                        ],
                        [
                            0.58,
                            8.005e-08
                        ],
                        [
                            0.59,
                            1.114e-07
                        ],
                        [
                            0.6,
                            1.59e-07
                        ],
                        [
                            0.61,
                            2.268e-07
                        ],
                        [
                            0.62,
                            3.235e-07
                        ],
                        [
                            0.63,
                            4.617e-07
                        ],
                        [
                            0.64,
                            6.567e-07
                        ],
                        [
                            0.65,
                            9.195e-07
                        ],
                        [
                            0.66,
                            1.262e-06
                        ],
                        [
                            0.67,
                            1.687e-06
                        ],
                        [
                            0.68,
                            2.17e-06
                        ],
                        [
                            0.69,
                            2.674e-06
                        ],
                        [
                            0.7,
                            3.138e-06
                        ],
                        [
                            0.71,
                            3.502e-06
                        ],
                        [
                            0.72,
                            3.721e-06
                        ],
                        [
                            0.73,
                            3.783e-06
                        ],
                        [
                            0.74,
                            3.718e-06
                        ],
                        [
                            0.75,
                            3.562e-06
                        ],
                        [
                            0.76,
                            3.363e-06
                        ],
                        [
                            0.77,
                            3.151e-06
                        ],
                        [
                            0.78,
                            2.947e-06
                        ],
                        [
                            0.79,
                            2.767e-06
                        ],
                        [
                            0.8,
                            2.607e-06
                        ],
                        [
                            0.81,
                            2.473e-06
                        ],
                        [
                            0.82,
                            2.362e-06
                        ],
                        [
                            0.83,
                            2.269e-06
                        ],
                        [
                            0.84,
                            2.194e-06
                        ],
                        [
                            0.85,
                            2.131e-06
                        ],
                        [
                            0.86,
                            2.08e-06
                        ],
                        [
                            0.87,
                            2.037e-06
                        ],
                        [
                            0.88,
                            2.004e-06
                        ],
                        [
                            0.89,
                            1.973e-06
                        ],
                        [
                            0.9,
                            1.949e-06
                        ],
                        [
                            0.91,
                            1.929e-06
                        ],
                        [
                            0.92,
                            1.91e-06
                        ],
                        [
                            0.93,
                            1.894e-06
                        ],
                        [
                            0.94,
                            1.881e-06
                        ],
                        [
                            0.95,
                            1.869e-06
                        ],
                        [
                            0.96,
                            1.856e-06
                        ],
                        [
                            0.97,
                            1.847e-06
                        ],
                        [
                            0.98,
                            1.837e-06
                        ],
                        [
                            0.98,
                            1.812e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.796e-06
                        ],
                        [
                            0.96,
                            1.783e-06
                        ],
                        [
                            0.95,
                            1.774e-06
                        ],
                        [
                            0.94,
                            1.766e-06
                        ],
                        [
                            0.93,
                            1.758e-06
                        ],
                        [
                            0.92,
                            1.751e-06
                        ],
                        [
                            0.91,
                            1.746e-06
                        ],
                        [
                            0.9,
                            1.739e-06
                        ],
                        [
                            0.89,
                            1.737e-06
                        ],
                        [
                            0.88,
                            1.729e-06
                        ],
                        [
                            0.87,
                            1.724e-06
                        ],
                        [
                            0.86,
                            1.714e-06
                        ],
                        [
                            0.85,
                            1.704e-06
                        ],
                        [
                            0.84,
                            1.686e-06
                        ],
                        [
                            0.83,
                            1.663e-06
                        ],
                        [
                            0.82,
                            1.628e-06
                        ],
                        [
                            0.81,
                            1.577e-06
                        ],
                        [
                            0.8,
                            1.505e-06
                        ],
                        [
                            0.79,
                            1.399e-06
                        ],
                        [
                            0.78,
                            1.253e-06
                        ],
                        [
                            0.77,
                            1.045e-06
                        ],
                        [
                            0.76,
                            7.687e-07
                        ],
                        [
                            0.75,
                            4.043e-07
                        ],
                        [
                            0.74,
                            -3.755e-08
                        ],
                        [
                            0.73,
                            -5.386e-07
                        ],
                        [
                            0.72,
                            -1.066e-06
                        ],
                        [
                            0.71,
                            -1.546e-06
                        ],
                        [
                            0.7,
                            -1.922e-06
                        ],
                        [
                            0.69,
                            -2.143e-06
                        ],
                        [
                            0.68,
                            -2.201e-06
                        ],
                        [
                            0.67,
                            -2.124e-06
                        ],
                        [
                            0.66,
                            -1.961e-06
                        ],
                        [
                            0.65,
                            -1.746e-06
                        ],
                        [
                            0.64,
                            -1.525e-06
                        ],
                        [
                            0.63,
                            -1.309e-06
                        ],
                        [
                            0.62,
                            -1.12e-06
                        ],
                        [
                            0.61,
                            -9.55e-07
                        ],
                        [
                            0.6,
                            -8.144e-07
                        ],
                        [
                            0.59,
                            -6.977e-07
                        ],
                        [
                            0.58,
                            -6.003e-07
                        ],
                        [
                            0.57,
                            -5.211e-07
                        ],
                        [
                            0.56,
                            -4.536e-07
                        ],
                        [
                            0.55,
                            -3.986e-07
                        ],
                        [
                            0.54,
                            -3.532e-07
                        ],
                        [
                            0.53,
                            -3.148e-07
                        ],
                        [
                            0.52,
                            -2.835e-07
                        ],
                        [
                            0.51,
                            -2.558e-07
                        ],
                        [
                            0.5,
                            -2.328e-07
                        ],
                        [
                            0.49,
                            -2.147e-07
                        ],
                        [
                            0.48,
                            -1.975e-07
                        ],
                        [
                            0.47,
                            -1.822e-07
                        ],
                        [
                            0.46,
                            -1.696e-07
                        ],
                        [
                            0.45,
                            -1.585e-07
                        ],
                        [
                            0.44,
                            -1.478e-07
                        ],
                        [
                            0.43,
                            -1.386e-07
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.73,
                        3.783e-06
                    ]
                ],
                "reverse": [
                    [
                        0.68,
                        -2.201e-06
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            -2.236e-09,
                            -1e-09,
                            8.348e-10,
                            2.37e-09,
                            5.752e-09,
                            6.976e-09,
                            1.036e-08,
                            1.282e-08,
                            1.558e-08,
                            2.049e-08,
                            2.571e-08,
                            3.216e-08,
                            4.321e-08,
                            5.764e-08,
                            8.005e-08,
                            1.114e-07,
                            1.59e-07,
                            2.268e-07,
                            3.235e-07,
                            4.617e-07,
                            6.567e-07,
                            9.195e-07,
                            1.262e-06,
                            1.687e-06,
                            2.17e-06,
                            2.674e-06,
                            3.138e-06,
                            3.502e-06,
                            3.721e-06,
                            3.783e-06,
                            3.718e-06,
                            3.562e-06,
                            3.363e-06,
                            3.151e-06,
                            2.947e-06,
                            2.767e-06,
                            2.607e-06,
                            2.473e-06,
                            2.362e-06,
                            2.269e-06,
                            2.194e-06,
                            2.131e-06,
                            2.08e-06,
                            2.037e-06,
                            2.004e-06,
                            1.973e-06,
                            1.949e-06,
                            1.929e-06,
                            1.91e-06,
                            1.894e-06,
                            1.881e-06,
                            1.869e-06,
                            1.856e-06,
                            1.847e-06,
                            1.837e-06,
                            1.812e-06,
                            1.796e-06,
                            1.783e-06,
                            1.774e-06,
                            1.766e-06,
                            1.758e-06,
                            1.751e-06,
                            1.746e-06,
                            1.739e-06,
                            1.737e-06,
                            1.729e-06,
                            1.724e-06,
                            1.714e-06,
                            1.704e-06,
                            1.686e-06,
                            1.663e-06,
                            1.628e-06,
                            1.577e-06,
                            1.505e-06,
                            1.399e-06,
                            1.253e-06,
                            1.045e-06,
                            7.687e-07,
                            4.043e-07,
                            -3.755e-08,
                            -5.386e-07,
                            -1.066e-06,
                            -1.546e-06,
                            -1.922e-06,
                            -2.143e-06,
                            -2.201e-06,
                            -2.124e-06,
                            -1.961e-06,
                            -1.746e-06,
                            -1.525e-06,
                            -1.309e-06,
                            -1.12e-06,
                            -9.55e-07,
                            -8.144e-07,
                            -6.977e-07,
                            -6.003e-07,
                            -5.211e-07,
                            -4.536e-07,
                            -3.986e-07,
                            -3.532e-07,
                            -3.148e-07,
                            -2.835e-07,
                            -2.558e-07,
                            -2.328e-07,
                            -2.147e-07,
                            -1.975e-07,
                            -1.822e-07,
                            -1.696e-07,
                            -1.585e-07,
                            -1.478e-07,
                            -1.386e-07
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "reversible"
                ],
                "e_half": [
                    0.705
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            -2.236e-09
                        ],
                        [
                            0.45,
                            -1e-09
                        ],
                        [
                            0.46,
                            8.348e-10
                        ],
                        [
                            0.47,
                            2.37e-09
                        ],
                        [
                            0.48,
                            5.752e-09
                        ],
                        [
                            0.49,
                            6.976e-09
                        ],
                        [
                            0.5,
                            1.036e-08
                        ],
                        [
                            0.51,
                            1.282e-08
                        ],
                        [
                            0.52,
                            1.558e-08
                        ],
                        [
                            0.53,
                            2.049e-08
                        ],
                        [
                            0.54,
                            2.571e-08
                        ],
                        [
                            0.55,
                            3.216e-08
                        ],
                        [
                            0.56,
                            4.321e-08
                        ],
                        [
                            0.57,
                            5.764e-08
                        ],
                        [
                            0.58,
                            8.005e-08
                        ],
                        [
                            0.59,
                            1.114e-07
                        ],
                        [
                            0.6,
                            1.59e-07
                        ],
                        [
                            0.61,
                            2.268e-07
                        ],
                        [
                            0.62,
                            3.235e-07
                        ],
                        [
                            0.63,
                            4.617e-07
                        ],
                        [
                            0.64,
                            6.567e-07
                        ],
                        [
                            0.65,
                            9.195e-07
                        ],
                        [
                            0.66,
                            1.262e-06
                        ],
                        [
                            0.67,
                            1.687e-06
                        ],
                        [
                            0.68,
                            2.17e-06
                        ],
                        [
                            0.69,
                            2.674e-06
                        ],
                        [
                            0.7,
                            3.138e-06
                        ],
                        [
                            0.71,
                            3.502e-06
                        ],
                        [
                            0.72,
                            3.721e-06
                        ],
                        [
                            0.73,
                            3.783e-06
                        ],
                        [
                            0.74,
                            3.718e-06
                        ],
                        [
                            0.75,
                            3.562e-06
                        ],
                        [
                            0.76,
                            3.363e-06
                        ],
                        [
                            0.77,
                            3.151e-06
                        ],
                        [
                            0.78,
                            2.947e-06
                        ],
                        [
                            0.79,
                            2.767e-06
                        ],
                        [
                            0.8,
                            2.607e-06
                        ],
                        [
                            0.81,
                            2.473e-06
                        ],
                        [
                            0.82,
                            2.362e-06
                        ],
                        [
                            0.83,
                            2.269e-06
                        ],
                        [
                            0.84,
                            2.194e-06
                        ],
                        [
                            0.85,
                            2.131e-06
                        ],
                        [
                            0.86,
                            2.08e-06
                        ],
                        [
                            0.87,
                            2.037e-06
                        ],
                        [
                            0.88,
                            2.004e-06
                        ],
                        [
                            0.89,
                            1.973e-06
                        ],
                        [
                            0.9,
                            1.949e-06
                        ],
                        [
                            0.91,
                            1.929e-06
                        ],
                        [
                            0.92,
                            1.91e-06
                        ],
                        [
                            0.93,
                            1.894e-06
                        ],
                        [
                            0.94,
                            1.881e-06
                        ],
                        [
                            0.95,
                            1.869e-06
                        ],
                        [
                            0.96,
                            1.856e-06
                        ],
                        [
                            0.97,
                            1.847e-06
                        ],
                        [
                            0.98,
                            1.837e-06
                        ],
                        [
                            0.98,
                            1.812e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.796e-06
                        ],
                        [
                            0.96,
                            1.783e-06
                        ],
                        [
                            0.95,
                            1.774e-06
                        ],
                        [
                            0.94,
                            1.766e-06
                        ],
                        [
                            0.93,
                            1.758e-06
                        ],
                        [
                            0.92,
                            1.751e-06
                        ],
                        [
                            0.91,
                            1.746e-06
                        ],
                        [
                            0.9,
                            1.739e-06
                        ],
                        [
                            0.89,
                            1.737e-06
                        ],
                        [
                            0.88,
                            1.729e-06
                        ],
                        [
                            0.87,
                            1.724e-06
                        ],
                        [
                            0.86,
                            1.714e-06
                        ],
                        [
                            0.85,
                            1.704e-06
                        ],
                        [
                            0.84,
                            1.686e-06
                        ],
                        [
                            0.83,
                            1.663e-06
                        ],
                        [
                            0.82,
                            1.628e-06
                        ],
                        [
                            0.81,
                            1.577e-06
                        ],
                        [
                            0.8,
                            1.505e-06
                        ],
                        [
                            0.79,
                            1.399e-06
                        ],
                        [
                            0.78,
                            1.253e-06
                        ],
                        [
                            0.77,
                            1.045e-06
                        ],
                        [
                            0.76,
                            7.687e-07
                        ],
                        [
                            0.75,
                            4.043e-07
                        ],
                        [
                            0.74,
                            -3.755e-08
                        ],
                        [
                            0.73,
                            -5.386e-07
                        ],
                        [
                            0.72,
                            -1.066e-06
                        ],
                        [
                            0.71,
                            -1.546e-06
                        ],
                        [
                            0.7,
                            -1.922e-06
                        ],
                        [
                            0.69,
                            -2.143e-06
                        ],
                        [
                            0.68,
                            -2.201e-06
                        ],
                        [
                            0.67,
                            -2.124e-06
                        ],
                        [
                            0.66,
                            -1.961e-06
                        ],
                        [
                            0.65,
                            -1.746e-06
                        ],
                        [
                            0.64,
                            -1.525e-06
                        ],
                        [
                            0.63,
                            -1.309e-06
                        ],
                        [
                            0.62,
                            -1.12e-06
                        ],
                        [
                            0.61,
                            -9.55e-07
                        ],
                        [
                            0.6,
                            -8.144e-07
                        ],
                        [
                            0.59,
                            -6.977e-07
                        ],
                        [
                            0.58,
                            -6.003e-07
                        ],
                        [
                            0.57,
                            -5.211e-07
                        ],
                        [
                            0.56,
                            -4.536e-07
                        ],
                        [
                            0.55,
                            -3.986e-07
                        ],
                        [
                            0.54,
                            -3.532e-07
                        ],
                        [
                            0.53,
                            -3.148e-07
                        ],
                        [
                            0.52,
                            -2.835e-07
                        ],
                        [
                            0.51,
                            -2.558e-07
                        ],
                        [
                            0.5,
                            -2.328e-07
                        ],
                        [
                            0.49,
                            -2.147e-07
                        ],
                        [
                            0.48,
                            -1.975e-07
                        ],
                        [
                            0.47,
                            -1.822e-07
                        ],
                        [
                            0.46,
                            -1.696e-07
                        ],
                        [
                            0.45,
                            -1.585e-07
                        ],
                        [
                            0.44,
                            -1.478e-07
                        ],
                        [
                            0.43,
                            -1.386e-07
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 3.783e-06,
                "current_anodic": -2.201e-06,
                "diffusion": 1.451e-09
            }
        },
        {
            "_id": "fb783838-345f-41f9-a6c5-27a454dd5882",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:49.928173",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv04_12_21_39.bin",
                "header": "CV cycle01_cv04_12_21_39",
                "note": "",
                "date_recorded": "2024-01-31T12:21:51",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.2,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            9.709e-09
                        ],
                        [
                            0.45,
                            1.128e-08
                        ],
                        [
                            0.46,
                            1.288e-08
                        ],
                        [
                            0.47,
                            1.44e-08
                        ],
                        [
                            0.48,
                            1.61e-08
                        ],
                        [
                            0.49,
                            1.765e-08
                        ],
                        [
                            0.5,
                            1.922e-08
                        ],
                        [
                            0.51,
                            2.084e-08
                        ],
                        [
                            0.52,
                            2.285e-08
                        ],
                        [
                            0.53,
                            2.556e-08
                        ],
                        [
                            0.54,
                            2.945e-08
                        ],
                        [
                            0.55,
                            3.504e-08
                        ],
                        [
                            0.56,
                            4.333e-08
                        ],
                        [
                            0.57,
                            5.541e-08
                        ],
                        [
                            0.58,
                            7.302e-08
                        ],
                        [
                            0.59,
                            9.83e-08
                        ],
                        [
                            0.6,
                            1.347e-07
                        ],
                        [
                            0.61,
                            1.869e-07
                        ],
                        [
                            0.62,
                            2.615e-07
                        ],
                        [
                            0.63,
                            3.662e-07
                        ],
                        [
                            0.64,
                            5.122e-07
                        ],
                        [
                            0.65,
                            7.108e-07
                        ],
                        [
                            0.66,
                            9.709e-07
                        ],
                        [
                            0.67,
                            1.296e-06
                        ],
                        [
                            0.68,
                            1.674e-06
                        ],
                        [
                            0.69,
                            2.079e-06
                        ],
                        [
                            0.7,
                            2.467e-06
                        ],
                        [
                            0.71,
                            2.787e-06
                        ],
                        [
                            0.72,
                            3.001e-06
                        ],
                        [
                            0.73,
                            3.091e-06
                        ],
                        [
                            0.74,
                            3.075e-06
                        ],
                        [
                            0.75,
                            2.979e-06
                        ],
                        [
                            0.76,
                            2.84e-06
                        ],
                        [
                            0.77,
                            2.686e-06
                        ],
                        [
                            0.78,
                            2.537e-06
                        ],
                        [
                            0.79,
                            2.403e-06
                        ],
                        [
                            0.8,
                            2.287e-06
                        ],
                        [
                            0.81,
                            2.191e-06
                        ],
                        [
                            0.82,
                            2.112e-06
                        ],
                        [
                            0.83,
                            2.049e-06
                        ],
                        [
                            0.84,
                            1.998e-06
                        ],
                        [
                            0.85,
                            1.957e-06
                        ],
                        [
                            0.86,
                            1.925e-06
                        ],
                        [
                            0.87,
                            1.898e-06
                        ],
                        [
                            0.88,
                            1.877e-06
                        ],
                        [
                            0.89,
                            1.86e-06
                        ],
                        [
                            0.9,
                            1.846e-06
                        ],
                        [
                            0.91,
                            1.834e-06
                        ],
                        [
                            0.92,
                            1.824e-06
                        ],
                        [
                            0.93,
                            1.816e-06
                        ],
                        [
                            0.94,
                            1.808e-06
                        ],
                        [
                            0.95,
                            1.801e-06
                        ],
                        [
                            0.96,
                            1.795e-06
                        ],
                        [
                            0.97,
                            1.79e-06
                        ],
                        [
                            0.98,
                            1.785e-06
                        ],
                        [
                            0.98,
                            1.751e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.74e-06
                        ],
                        [
                            0.96,
                            1.733e-06
                        ],
                        [
                            0.95,
                            1.727e-06
                        ],
                        [
                            0.94,
                            1.721e-06
                        ],
                        [
                            0.93,
                            1.717e-06
                        ],
                        [
                            0.92,
                            1.713e-06
                        ],
                        [
                            0.91,
                            1.709e-06
                        ],
                        [
                            0.9,
                            1.705e-06
                        ],
                        [
                            0.89,
                            1.7e-06
                        ],
                        [
                            0.88,
                            1.695e-06
                        ],
                        [
                            0.87,
                            1.688e-06
                        ],
                        [
                            0.86,
                            1.68e-06
                        ],
                        [
                            0.85,
                            1.669e-06
                        ],
                        [
                            0.84,
                            1.653e-06
                        ],
                        [
                            0.83,
                            1.631e-06
                        ],
                        [
                            0.82,
                            1.599e-06
                        ],
                        [
                            0.81,
                            1.553e-06
                        ],
                        [
                            0.8,
                            1.488e-06
                        ],
                        [
                            0.79,
                            1.395e-06
                        ],
                        [
                            0.78,
                            1.264e-06
                        ],
                        [
                            0.77,
                            1.083e-06
                        ],
                        [
                            0.76,
                            8.421e-07
                        ],
                        [
                            0.75,
                            5.339e-07
                        ],
                        [
                            0.74,
                            1.616e-07
                        ],
                        [
                            0.73,
                            -2.536e-07
                        ],
                        [
                            0.72,
                            -6.743e-07
                        ],
                        [
                            0.71,
                            -1.047e-06
                        ],
                        [
                            0.7,
                            -1.321e-06
                        ],
                        [
                            0.69,
                            -1.47e-06
                        ],
                        [
                            0.68,
                            -1.496e-06
                        ],
                        [
                            0.67,
                            -1.425e-06
                        ],
                        [
                            0.66,
                            -1.294e-06
                        ],
                        [
                            0.65,
                            -1.135e-06
                        ],
                        [
                            0.64,
                            -9.733e-07
                        ],
                        [
                            0.63,
                            -8.235e-07
                        ],
                        [
                            0.62,
                            -6.923e-07
                        ],
                        [
                            0.61,
                            -5.818e-07
                        ],
                        [
                            0.6,
                            -4.9e-07
                        ],
                        [
                            0.59,
                            -4.151e-07
                        ],
                        [
                            0.58,
                            -3.545e-07
                        ],
                        [
                            0.57,
                            -3.056e-07
                        ],
                        [
                            0.56,
                            -2.663e-07
                        ],
                        [
                            0.55,
                            -2.345e-07
                        ],
                        [
                            0.54,
                            -2.087e-07
                        ],
                        [
                            0.53,
                            -1.88e-07
                        ],
                        [
                            0.52,
                            -1.711e-07
                        ],
                        [
                            0.51,
                            -1.573e-07
                        ],
                        [
                            0.5,
                            -1.458e-07
                        ],
                        [
                            0.49,
                            -1.364e-07
                        ],
                        [
                            0.48,
                            -1.281e-07
                        ],
                        [
                            0.47,
                            -1.209e-07
                        ],
                        [
                            0.46,
                            -1.147e-07
                        ],
                        [
                            0.45,
                            -1.091e-07
                        ],
                        [
                            0.44,
                            -1.04e-07
                        ],
                        [
                            0.43,
                            -9.948e-08
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.73,
                        3.091e-06
                    ]
                ],
                "reverse": [
                    [
                        0.68,
                        -1.496e-06
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            9.709e-09,
                            1.128e-08,
                            1.288e-08,
                            1.44e-08,
                            1.61e-08,
                            1.765e-08,
                            1.922e-08,
                            2.084e-08,
                            2.285e-08,
                            2.556e-08,
                            2.945e-08,
                            3.504e-08,
                            4.333e-08,
                            5.541e-08,
                            7.302e-08,
                            9.83e-08,
                            1.347e-07,
                            1.869e-07,
                            2.615e-07,
                            3.662e-07,
                            5.122e-07,
                            7.108e-07,
                            9.709e-07,
                            1.296e-06,
                            1.674e-06,
                            2.079e-06,
                            2.467e-06,
                            2.787e-06,
                            3.001e-06,
                            3.091e-06,
                            3.075e-06,
                            2.979e-06,
                            2.84e-06,
                            2.686e-06,
                            2.537e-06,
                            2.403e-06,
                            2.287e-06,
                            2.191e-06,
                            2.112e-06,
                            2.049e-06,
                            1.998e-06,
                            1.957e-06,
                            1.925e-06,
                            1.898e-06,
                            1.877e-06,
                            1.86e-06,
                            1.846e-06,
                            1.834e-06,
                            1.824e-06,
                            1.816e-06,
                            1.808e-06,
                            1.801e-06,
                            1.795e-06,
                            1.79e-06,
                            1.785e-06,
                            1.751e-06,
                            1.74e-06,
                            1.733e-06,
                            1.727e-06,
                            1.721e-06,
                            1.717e-06,
                            1.713e-06,
                            1.709e-06,
                            1.705e-06,
                            1.7e-06,
                            1.695e-06,
                            1.688e-06,
                            1.68e-06,
                            1.669e-06,
                            1.653e-06,
                            1.631e-06,
                            1.599e-06,
                            1.553e-06,
                            1.488e-06,
                            1.395e-06,
                            1.264e-06,
                            1.083e-06,
                            8.421e-07,
                            5.339e-07,
                            1.616e-07,
                            -2.536e-07,
                            -6.743e-07,
                            -1.047e-06,
                            -1.321e-06,
                            -1.47e-06,
                            -1.496e-06,
                            -1.425e-06,
                            -1.294e-06,
                            -1.135e-06,
                            -9.733e-07,
                            -8.235e-07,
                            -6.923e-07,
                            -5.818e-07,
                            -4.9e-07,
                            -4.151e-07,
                            -3.545e-07,
                            -3.056e-07,
                            -2.663e-07,
                            -2.345e-07,
                            -2.087e-07,
                            -1.88e-07,
                            -1.711e-07,
                            -1.573e-07,
                            -1.458e-07,
                            -1.364e-07,
                            -1.281e-07,
                            -1.209e-07,
                            -1.147e-07,
                            -1.091e-07,
                            -1.04e-07,
                            -9.948e-08
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "reversible"
                ],
                "e_half": [
                    0.705
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            9.709e-09
                        ],
                        [
                            0.45,
                            1.128e-08
                        ],
                        [
                            0.46,
                            1.288e-08
                        ],
                        [
                            0.47,
                            1.44e-08
                        ],
                        [
                            0.48,
                            1.61e-08
                        ],
                        [
                            0.49,
                            1.765e-08
                        ],
                        [
                            0.5,
                            1.922e-08
                        ],
                        [
                            0.51,
                            2.084e-08
                        ],
                        [
                            0.52,
                            2.285e-08
                        ],
                        [
                            0.53,
                            2.556e-08
                        ],
                        [
                            0.54,
                            2.945e-08
                        ],
                        [
                            0.55,
                            3.504e-08
                        ],
                        [
                            0.56,
                            4.333e-08
                        ],
                        [
                            0.57,
                            5.541e-08
                        ],
                        [
                            0.58,
                            7.302e-08
                        ],
                        [
                            0.59,
                            9.83e-08
                        ],
                        [
                            0.6,
                            1.347e-07
                        ],
                        [
                            0.61,
                            1.869e-07
                        ],
                        [
                            0.62,
                            2.615e-07
                        ],
                        [
                            0.63,
                            3.662e-07
                        ],
                        [
                            0.64,
                            5.122e-07
                        ],
                        [
                            0.65,
                            7.108e-07
                        ],
                        [
                            0.66,
                            9.709e-07
                        ],
                        [
                            0.67,
                            1.296e-06
                        ],
                        [
                            0.68,
                            1.674e-06
                        ],
                        [
                            0.69,
                            2.079e-06
                        ],
                        [
                            0.7,
                            2.467e-06
                        ],
                        [
                            0.71,
                            2.787e-06
                        ],
                        [
                            0.72,
                            3.001e-06
                        ],
                        [
                            0.73,
                            3.091e-06
                        ],
                        [
                            0.74,
                            3.075e-06
                        ],
                        [
                            0.75,
                            2.979e-06
                        ],
                        [
                            0.76,
                            2.84e-06
                        ],
                        [
                            0.77,
                            2.686e-06
                        ],
                        [
                            0.78,
                            2.537e-06
                        ],
                        [
                            0.79,
                            2.403e-06
                        ],
                        [
                            0.8,
                            2.287e-06
                        ],
                        [
                            0.81,
                            2.191e-06
                        ],
                        [
                            0.82,
                            2.112e-06
                        ],
                        [
                            0.83,
                            2.049e-06
                        ],
                        [
                            0.84,
                            1.998e-06
                        ],
                        [
                            0.85,
                            1.957e-06
                        ],
                        [
                            0.86,
                            1.925e-06
                        ],
                        [
                            0.87,
                            1.898e-06
                        ],
                        [
                            0.88,
                            1.877e-06
                        ],
                        [
                            0.89,
                            1.86e-06
                        ],
                        [
                            0.9,
                            1.846e-06
                        ],
                        [
                            0.91,
                            1.834e-06
                        ],
                        [
                            0.92,
                            1.824e-06
                        ],
                        [
                            0.93,
                            1.816e-06
                        ],
                        [
                            0.94,
                            1.808e-06
                        ],
                        [
                            0.95,
                            1.801e-06
                        ],
                        [
                            0.96,
                            1.795e-06
                        ],
                        [
                            0.97,
                            1.79e-06
                        ],
                        [
                            0.98,
                            1.785e-06
                        ],
                        [
                            0.98,
                            1.751e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.74e-06
                        ],
                        [
                            0.96,
                            1.733e-06
                        ],
                        [
                            0.95,
                            1.727e-06
                        ],
                        [
                            0.94,
                            1.721e-06
                        ],
                        [
                            0.93,
                            1.717e-06
                        ],
                        [
                            0.92,
                            1.713e-06
                        ],
                        [
                            0.91,
                            1.709e-06
                        ],
                        [
                            0.9,
                            1.705e-06
                        ],
                        [
                            0.89,
                            1.7e-06
                        ],
                        [
                            0.88,
                            1.695e-06
                        ],
                        [
                            0.87,
                            1.688e-06
                        ],
                        [
                            0.86,
                            1.68e-06
                        ],
                        [
                            0.85,
                            1.669e-06
                        ],
                        [
                            0.84,
                            1.653e-06
                        ],
                        [
                            0.83,
                            1.631e-06
                        ],
                        [
                            0.82,
                            1.599e-06
                        ],
                        [
                            0.81,
                            1.553e-06
                        ],
                        [
                            0.8,
                            1.488e-06
                        ],
                        [
                            0.79,
                            1.395e-06
                        ],
                        [
                            0.78,
                            1.264e-06
                        ],
                        [
                            0.77,
                            1.083e-06
                        ],
                        [
                            0.76,
                            8.421e-07
                        ],
                        [
                            0.75,
                            5.339e-07
                        ],
                        [
                            0.74,
                            1.616e-07
                        ],
                        [
                            0.73,
                            -2.536e-07
                        ],
                        [
                            0.72,
                            -6.743e-07
                        ],
                        [
                            0.71,
                            -1.047e-06
                        ],
                        [
                            0.7,
                            -1.321e-06
                        ],
                        [
                            0.69,
                            -1.47e-06
                        ],
                        [
                            0.68,
                            -1.496e-06
                        ],
                        [
                            0.67,
                            -1.425e-06
                        ],
                        [
                            0.66,
                            -1.294e-06
                        ],
                        [
                            0.65,
                            -1.135e-06
                        ],
                        [
                            0.64,
                            -9.733e-07
                        ],
                        [
                            0.63,
                            -8.235e-07
                        ],
                        [
                            0.62,
                            -6.923e-07
                        ],
                        [
                            0.61,
                            -5.818e-07
                        ],
                        [
                            0.6,
                            -4.9e-07
                        ],
                        [
                            0.59,
                            -4.151e-07
                        ],
                        [
                            0.58,
                            -3.545e-07
                        ],
                        [
                            0.57,
                            -3.056e-07
                        ],
                        [
                            0.56,
                            -2.663e-07
                        ],
                        [
                            0.55,
                            -2.345e-07
                        ],
                        [
                            0.54,
                            -2.087e-07
                        ],
                        [
                            0.53,
                            -1.88e-07
                        ],
                        [
                            0.52,
                            -1.711e-07
                        ],
                        [
                            0.51,
                            -1.573e-07
                        ],
                        [
                            0.5,
                            -1.458e-07
                        ],
                        [
                            0.49,
                            -1.364e-07
                        ],
                        [
                            0.48,
                            -1.281e-07
                        ],
                        [
                            0.47,
                            -1.209e-07
                        ],
                        [
                            0.46,
                            -1.147e-07
                        ],
                        [
                            0.45,
                            -1.091e-07
                        ],
                        [
                            0.44,
                            -1.04e-07
                        ],
                        [
                            0.43,
                            -9.948e-08
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 3.091e-06,
                "current_anodic": -1.496e-06,
                "diffusion": 1.451e-09
            }
        },
        {
            "_id": "0692b884-b880-4a6e-9e61-7903454b38cb",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:51.788339",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv05_12_21_57.bin",
                "header": "CV cycle01_cv05_12_21_57",
                "note": "",
                "date_recorded": "2024-01-31T12:22:15",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.1,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            2.773e-09
                        ],
                        [
                            0.45,
                            3.971e-09
                        ],
                        [
                            0.46,
                            5.082e-09
                        ],
                        [
                            0.47,
                            6.221e-09
                        ],
                        [
                            0.48,
                            7.479e-09
                        ],
                        [
                            0.49,
                            8.599e-09
                        ],
                        [
                            0.5,
                            9.719e-09
                        ],
                        [
                            0.51,
                            1.105e-08
                        ],
                        [
                            0.52,
                            1.242e-08
                        ],
                        [
                            0.53,
                            1.42e-08
                        ],
                        [
                            0.54,
                            1.671e-08
                        ],
                        [
                            0.55,
                            2.029e-08
                        ],
                        [
                            0.56,
                            2.55e-08
                        ],
                        [
                            0.57,
                            3.304e-08
                        ],
                        [
                            0.58,
                            4.383e-08
                        ],
                        [
                            0.59,
                            5.96e-08
                        ],
                        [
                            0.6,
                            8.22e-08
                        ],
                        [
                            0.61,
                            1.147e-07
                        ],
                        [
                            0.62,
                            1.614e-07
                        ],
                        [
                            0.63,
                            2.27e-07
                        ],
                        [
                            0.64,
                            3.193e-07
                        ],
                        [
                            0.65,
                            4.455e-07
                        ],
                        [
                            0.66,
                            6.138e-07
                        ],
                        [
                            0.67,
                            8.283e-07
                        ],
                        [
                            0.68,
                            1.085e-06
                        ],
                        [
                            0.69,
                            1.37e-06
                        ],
                        [
                            0.7,
                            1.655e-06
                        ],
                        [
                            0.71,
                            1.907e-06
                        ],
                        [
                            0.72,
                            2.098e-06
                        ],
                        [
                            0.73,
                            2.21e-06
                        ],
                        [
                            0.74,
                            2.249e-06
                        ],
                        [
                            0.75,
                            2.232e-06
                        ],
                        [
                            0.76,
                            2.181e-06
                        ],
                        [
                            0.77,
                            2.114e-06
                        ],
                        [
                            0.78,
                            2.046e-06
                        ],
                        [
                            0.79,
                            1.984e-06
                        ],
                        [
                            0.8,
                            1.931e-06
                        ],
                        [
                            0.81,
                            1.887e-06
                        ],
                        [
                            0.82,
                            1.853e-06
                        ],
                        [
                            0.83,
                            1.825e-06
                        ],
                        [
                            0.84,
                            1.803e-06
                        ],
                        [
                            0.85,
                            1.786e-06
                        ],
                        [
                            0.86,
                            1.773e-06
                        ],
                        [
                            0.87,
                            1.762e-06
                        ],
                        [
                            0.88,
                            1.754e-06
                        ],
                        [
                            0.89,
                            1.747e-06
                        ],
                        [
                            0.9,
                            1.741e-06
                        ],
                        [
                            0.91,
                            1.736e-06
                        ],
                        [
                            0.92,
                            1.732e-06
                        ],
                        [
                            0.93,
                            1.728e-06
                        ],
                        [
                            0.94,
                            1.725e-06
                        ],
                        [
                            0.95,
                            1.722e-06
                        ],
                        [
                            0.96,
                            1.72e-06
                        ],
                        [
                            0.97,
                            1.717e-06
                        ],
                        [
                            0.98,
                            1.715e-06
                        ],
                        [
                            0.98,
                            1.698e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.692e-06
                        ],
                        [
                            0.96,
                            1.688e-06
                        ],
                        [
                            0.95,
                            1.684e-06
                        ],
                        [
                            0.94,
                            1.681e-06
                        ],
                        [
                            0.93,
                            1.678e-06
                        ],
                        [
                            0.92,
                            1.675e-06
                        ],
                        [
                            0.91,
                            1.672e-06
                        ],
                        [
                            0.9,
                            1.669e-06
                        ],
                        [
                            0.89,
                            1.666e-06
                        ],
                        [
                            0.88,
                            1.662e-06
                        ],
                        [
                            0.87,
                            1.657e-06
                        ],
                        [
                            0.86,
                            1.651e-06
                        ],
                        [
                            0.85,
                            1.643e-06
                        ],
                        [
                            0.84,
                            1.631e-06
                        ],
                        [
                            0.83,
                            1.614e-06
                        ],
                        [
                            0.82,
                            1.59e-06
                        ],
                        [
                            0.81,
                            1.556e-06
                        ],
                        [
                            0.8,
                            1.507e-06
                        ],
                        [
                            0.79,
                            1.437e-06
                        ],
                        [
                            0.78,
                            1.34e-06
                        ],
                        [
                            0.77,
                            1.206e-06
                        ],
                        [
                            0.76,
                            1.029e-06
                        ],
                        [
                            0.75,
                            8.041e-07
                        ],
                        [
                            0.74,
                            5.357e-07
                        ],
                        [
                            0.73,
                            2.405e-07
                        ],
                        [
                            0.72,
                            -5.454e-08
                        ],
                        [
                            0.71,
                            -3.133e-07
                        ],
                        [
                            0.7,
                            -5.05e-07
                        ],
                        [
                            0.69,
                            -6.154e-07
                        ],
                        [
                            0.68,
                            -6.493e-07
                        ],
                        [
                            0.67,
                            -6.255e-07
                        ],
                        [
                            0.66,
                            -5.677e-07
                        ],
                        [
                            0.65,
                            -4.944e-07
                        ],
                        [
                            0.64,
                            -4.209e-07
                        ],
                        [
                            0.63,
                            -3.539e-07
                        ],
                        [
                            0.62,
                            -2.97e-07
                        ],
                        [
                            0.61,
                            -2.504e-07
                        ],
                        [
                            0.6,
                            -2.127e-07
                        ],
                        [
                            0.59,
                            -1.828e-07
                        ],
                        [
                            0.58,
                            -1.593e-07
                        ],
                        [
                            0.57,
                            -1.406e-07
                        ],
                        [
                            0.56,
                            -1.26e-07
                        ],
                        [
                            0.55,
                            -1.142e-07
                        ],
                        [
                            0.54,
                            -1.05e-07
                        ],
                        [
                            0.53,
                            -9.757e-08
                        ],
                        [
                            0.52,
                            -9.163e-08
                        ],
                        [
                            0.51,
                            -8.677e-08
                        ],
                        [
                            0.5,
                            -8.28e-08
                        ],
                        [
                            0.49,
                            -7.937e-08
                        ],
                        [
                            0.48,
                            -7.635e-08
                        ],
                        [
                            0.47,
                            -7.366e-08
                        ],
                        [
                            0.46,
                            -7.119e-08
                        ],
                        [
                            0.45,
                            -6.887e-08
                        ],
                        [
                            0.44,
                            -6.672e-08
                        ],
                        [
                            0.43,
                            -6.483e-08
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.74,
                        2.249e-06
                    ]
                ],
                "reverse": [
                    [
                        0.68,
                        -6.493e-07
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            2.773e-09,
                            3.971e-09,
                            5.082e-09,
                            6.221e-09,
                            7.479e-09,
                            8.599e-09,
                            9.719e-09,
                            1.105e-08,
                            1.242e-08,
                            1.42e-08,
                            1.671e-08,
                            2.029e-08,
                            2.55e-08,
                            3.304e-08,
                            4.383e-08,
                            5.96e-08,
                            8.22e-08,
                            1.147e-07,
                            1.614e-07,
                            2.27e-07,
                            3.193e-07,
                            4.455e-07,
                            6.138e-07,
                            8.283e-07,
                            1.085e-06,
                            1.37e-06,
                            1.655e-06,
                            1.907e-06,
                            2.098e-06,
                            2.21e-06,
                            2.249e-06,
                            2.232e-06,
                            2.181e-06,
                            2.114e-06,
                            2.046e-06,
                            1.984e-06,
                            1.931e-06,
                            1.887e-06,
                            1.853e-06,
                            1.825e-06,
                            1.803e-06,
                            1.786e-06,
                            1.773e-06,
                            1.762e-06,
                            1.754e-06,
                            1.747e-06,
                            1.741e-06,
                            1.736e-06,
                            1.732e-06,
                            1.728e-06,
                            1.725e-06,
                            1.722e-06,
                            1.72e-06,
                            1.717e-06,
                            1.715e-06,
                            1.698e-06,
                            1.692e-06,
                            1.688e-06,
                            1.684e-06,
                            1.681e-06,
                            1.678e-06,
                            1.675e-06,
                            1.672e-06,
                            1.669e-06,
                            1.666e-06,
                            1.662e-06,
                            1.657e-06,
                            1.651e-06,
                            1.643e-06,
                            1.631e-06,
                            1.614e-06,
                            1.59e-06,
                            1.556e-06,
                            1.507e-06,
                            1.437e-06,
                            1.34e-06,
                            1.206e-06,
                            1.029e-06,
                            8.041e-07,
                            5.357e-07,
                            2.405e-07,
                            -5.454e-08,
                            -3.133e-07,
                            -5.05e-07,
                            -6.154e-07,
                            -6.493e-07,
                            -6.255e-07,
                            -5.677e-07,
                            -4.944e-07,
                            -4.209e-07,
                            -3.539e-07,
                            -2.97e-07,
                            -2.504e-07,
                            -2.127e-07,
                            -1.828e-07,
                            -1.593e-07,
                            -1.406e-07,
                            -1.26e-07,
                            -1.142e-07,
                            -1.05e-07,
                            -9.757e-08,
                            -9.163e-08,
                            -8.677e-08,
                            -8.28e-08,
                            -7.937e-08,
                            -7.635e-08,
                            -7.366e-08,
                            -7.119e-08,
                            -6.887e-08,
                            -6.672e-08,
                            -6.483e-08
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "reversible"
                ],
                "e_half": [
                    0.71
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            2.773e-09
                        ],
                        [
                            0.45,
                            3.971e-09
                        ],
                        [
                            0.46,
                            5.082e-09
                        ],
                        [
                            0.47,
                            6.221e-09
                        ],
                        [
                            0.48,
                            7.479e-09
                        ],
                        [
                            0.49,
                            8.599e-09
                        ],
                        [
                            0.5,
                            9.719e-09
                        ],
                        [
                            0.51,
                            1.105e-08
                        ],
                        [
                            0.52,
                            1.242e-08
                        ],
                        [
                            0.53,
                            1.42e-08
                        ],
                        [
                            0.54,
                            1.671e-08
                        ],
                        [
                            0.55,
                            2.029e-08
                        ],
                        [
                            0.56,
                            2.55e-08
                        ],
                        [
                            0.57,
                            3.304e-08
                        ],
                        [
                            0.58,
                            4.383e-08
                        ],
                        [
                            0.59,
                            5.96e-08
                        ],
                        [
                            0.6,
                            8.22e-08
                        ],
                        [
                            0.61,
                            1.147e-07
                        ],
                        [
                            0.62,
                            1.614e-07
                        ],
                        [
                            0.63,
                            2.27e-07
                        ],
                        [
                            0.64,
                            3.193e-07
                        ],
                        [
                            0.65,
                            4.455e-07
                        ],
                        [
                            0.66,
                            6.138e-07
                        ],
                        [
                            0.67,
                            8.283e-07
                        ],
                        [
                            0.68,
                            1.085e-06
                        ],
                        [
                            0.69,
                            1.37e-06
                        ],
                        [
                            0.7,
                            1.655e-06
                        ],
                        [
                            0.71,
                            1.907e-06
                        ],
                        [
                            0.72,
                            2.098e-06
                        ],
                        [
                            0.73,
                            2.21e-06
                        ],
                        [
                            0.74,
                            2.249e-06
                        ],
                        [
                            0.75,
                            2.232e-06
                        ],
                        [
                            0.76,
                            2.181e-06
                        ],
                        [
                            0.77,
                            2.114e-06
                        ],
                        [
                            0.78,
                            2.046e-06
                        ],
                        [
                            0.79,
                            1.984e-06
                        ],
                        [
                            0.8,
                            1.931e-06
                        ],
                        [
                            0.81,
                            1.887e-06
                        ],
                        [
                            0.82,
                            1.853e-06
                        ],
                        [
                            0.83,
                            1.825e-06
                        ],
                        [
                            0.84,
                            1.803e-06
                        ],
                        [
                            0.85,
                            1.786e-06
                        ],
                        [
                            0.86,
                            1.773e-06
                        ],
                        [
                            0.87,
                            1.762e-06
                        ],
                        [
                            0.88,
                            1.754e-06
                        ],
                        [
                            0.89,
                            1.747e-06
                        ],
                        [
                            0.9,
                            1.741e-06
                        ],
                        [
                            0.91,
                            1.736e-06
                        ],
                        [
                            0.92,
                            1.732e-06
                        ],
                        [
                            0.93,
                            1.728e-06
                        ],
                        [
                            0.94,
                            1.725e-06
                        ],
                        [
                            0.95,
                            1.722e-06
                        ],
                        [
                            0.96,
                            1.72e-06
                        ],
                        [
                            0.97,
                            1.717e-06
                        ],
                        [
                            0.98,
                            1.715e-06
                        ],
                        [
                            0.98,
                            1.698e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.692e-06
                        ],
                        [
                            0.96,
                            1.688e-06
                        ],
                        [
                            0.95,
                            1.684e-06
                        ],
                        [
                            0.94,
                            1.681e-06
                        ],
                        [
                            0.93,
                            1.678e-06
                        ],
                        [
                            0.92,
                            1.675e-06
                        ],
                        [
                            0.91,
                            1.672e-06
                        ],
                        [
                            0.9,
                            1.669e-06
                        ],
                        [
                            0.89,
                            1.666e-06
                        ],
                        [
                            0.88,
                            1.662e-06
                        ],
                        [
                            0.87,
                            1.657e-06
                        ],
                        [
                            0.86,
                            1.651e-06
                        ],
                        [
                            0.85,
                            1.643e-06
                        ],
                        [
                            0.84,
                            1.631e-06
                        ],
                        [
                            0.83,
                            1.614e-06
                        ],
                        [
                            0.82,
                            1.59e-06
                        ],
                        [
                            0.81,
                            1.556e-06
                        ],
                        [
                            0.8,
                            1.507e-06
                        ],
                        [
                            0.79,
                            1.437e-06
                        ],
                        [
                            0.78,
                            1.34e-06
                        ],
                        [
                            0.77,
                            1.206e-06
                        ],
                        [
                            0.76,
                            1.029e-06
                        ],
                        [
                            0.75,
                            8.041e-07
                        ],
                        [
                            0.74,
                            5.357e-07
                        ],
                        [
                            0.73,
                            2.405e-07
                        ],
                        [
                            0.72,
                            -5.454e-08
                        ],
                        [
                            0.71,
                            -3.133e-07
                        ],
                        [
                            0.7,
                            -5.05e-07
                        ],
                        [
                            0.69,
                            -6.154e-07
                        ],
                        [
                            0.68,
                            -6.493e-07
                        ],
                        [
                            0.67,
                            -6.255e-07
                        ],
                        [
                            0.66,
                            -5.677e-07
                        ],
                        [
                            0.65,
                            -4.944e-07
                        ],
                        [
                            0.64,
                            -4.209e-07
                        ],
                        [
                            0.63,
                            -3.539e-07
                        ],
                        [
                            0.62,
                            -2.97e-07
                        ],
                        [
                            0.61,
                            -2.504e-07
                        ],
                        [
                            0.6,
                            -2.127e-07
                        ],
                        [
                            0.59,
                            -1.828e-07
                        ],
                        [
                            0.58,
                            -1.593e-07
                        ],
                        [
                            0.57,
                            -1.406e-07
                        ],
                        [
                            0.56,
                            -1.26e-07
                        ],
                        [
                            0.55,
                            -1.142e-07
                        ],
                        [
                            0.54,
                            -1.05e-07
                        ],
                        [
                            0.53,
                            -9.757e-08
                        ],
                        [
                            0.52,
                            -9.163e-08
                        ],
                        [
                            0.51,
                            -8.677e-08
                        ],
                        [
                            0.5,
                            -8.28e-08
                        ],
                        [
                            0.49,
                            -7.937e-08
                        ],
                        [
                            0.48,
                            -7.635e-08
                        ],
                        [
                            0.47,
                            -7.366e-08
                        ],
                        [
                            0.46,
                            -7.119e-08
                        ],
                        [
                            0.45,
                            -6.887e-08
                        ],
                        [
                            0.44,
                            -6.672e-08
                        ],
                        [
                            0.43,
                            -6.483e-08
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 2.249e-06,
                "current_anodic": -6.493e-07,
                "diffusion": 1.451e-09
            }
        },
        {
            "_id": "7a6661d4-0d6f-4ff8-af27-89ffe035c250",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:53.770269",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv06_12_22_20.bin",
                "header": "CV cycle01_cv06_12_22_20",
                "note": "",
                "date_recorded": "2024-01-31T12:22:42",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.075,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            -1.75e-09
                        ],
                        [
                            0.45,
                            -7.178e-10
                        ],
                        [
                            0.46,
                            3.149e-10
                        ],
                        [
                            0.47,
                            1.407e-09
                        ],
                        [
                            0.48,
                            2.528e-09
                        ],
                        [
                            0.49,
                            3.679e-09
                        ],
                        [
                            0.5,
                            4.741e-09
                        ],
                        [
                            0.51,
                            5.98e-09
                        ],
                        [
                            0.52,
                            7.337e-09
                        ],
                        [
                            0.53,
                            8.93e-09
                        ],
                        [
                            0.54,
                            1.108e-08
                        ],
                        [
                            0.55,
                            1.413e-08
                        ],
                        [
                            0.56,
                            1.858e-08
                        ],
                        [
                            0.57,
                            2.491e-08
                        ],
                        [
                            0.58,
                            3.401e-08
                        ],
                        [
                            0.59,
                            4.718e-08
                        ],
                        [
                            0.6,
                            6.622e-08
                        ],
                        [
                            0.61,
                            9.351e-08
                        ],
                        [
                            0.62,
                            1.326e-07
                        ],
                        [
                            0.63,
                            1.877e-07
                        ],
                        [
                            0.64,
                            2.651e-07
                        ],
                        [
                            0.65,
                            3.715e-07
                        ],
                        [
                            0.66,
                            5.135e-07
                        ],
                        [
                            0.67,
                            6.959e-07
                        ],
                        [
                            0.68,
                            9.158e-07
                        ],
                        [
                            0.69,
                            1.164e-06
                        ],
                        [
                            0.7,
                            1.417e-06
                        ],
                        [
                            0.71,
                            1.648e-06
                        ],
                        [
                            0.72,
                            1.832e-06
                        ],
                        [
                            0.73,
                            1.952e-06
                        ],
                        [
                            0.74,
                            2.011e-06
                        ],
                        [
                            0.75,
                            2.021e-06
                        ],
                        [
                            0.76,
                            1.999e-06
                        ],
                        [
                            0.77,
                            1.962e-06
                        ],
                        [
                            0.78,
                            1.92e-06
                        ],
                        [
                            0.79,
                            1.88e-06
                        ],
                        [
                            0.8,
                            1.845e-06
                        ],
                        [
                            0.81,
                            1.816e-06
                        ],
                        [
                            0.82,
                            1.792e-06
                        ],
                        [
                            0.83,
                            1.774e-06
                        ],
                        [
                            0.84,
                            1.759e-06
                        ],
                        [
                            0.85,
                            1.747e-06
                        ],
                        [
                            0.86,
                            1.737e-06
                        ],
                        [
                            0.87,
                            1.73e-06
                        ],
                        [
                            0.88,
                            1.724e-06
                        ],
                        [
                            0.89,
                            1.718e-06
                        ],
                        [
                            0.9,
                            1.714e-06
                        ],
                        [
                            0.91,
                            1.71e-06
                        ],
                        [
                            0.92,
                            1.707e-06
                        ],
                        [
                            0.93,
                            1.704e-06
                        ],
                        [
                            0.94,
                            1.702e-06
                        ],
                        [
                            0.95,
                            1.699e-06
                        ],
                        [
                            0.96,
                            1.697e-06
                        ],
                        [
                            0.97,
                            1.695e-06
                        ],
                        [
                            0.98,
                            1.694e-06
                        ],
                        [
                            0.98,
                            1.68e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.675e-06
                        ],
                        [
                            0.96,
                            1.672e-06
                        ],
                        [
                            0.95,
                            1.669e-06
                        ],
                        [
                            0.94,
                            1.666e-06
                        ],
                        [
                            0.93,
                            1.664e-06
                        ],
                        [
                            0.92,
                            1.661e-06
                        ],
                        [
                            0.91,
                            1.659e-06
                        ],
                        [
                            0.9,
                            1.656e-06
                        ],
                        [
                            0.89,
                            1.654e-06
                        ],
                        [
                            0.88,
                            1.65e-06
                        ],
                        [
                            0.87,
                            1.645e-06
                        ],
                        [
                            0.86,
                            1.639e-06
                        ],
                        [
                            0.85,
                            1.632e-06
                        ],
                        [
                            0.84,
                            1.62e-06
                        ],
                        [
                            0.83,
                            1.605e-06
                        ],
                        [
                            0.82,
                            1.583e-06
                        ],
                        [
                            0.81,
                            1.551e-06
                        ],
                        [
                            0.8,
                            1.506e-06
                        ],
                        [
                            0.79,
                            1.443e-06
                        ],
                        [
                            0.78,
                            1.354e-06
                        ],
                        [
                            0.77,
                            1.233e-06
                        ],
                        [
                            0.76,
                            1.072e-06
                        ],
                        [
                            0.75,
                            8.696e-07
                        ],
                        [
                            0.74,
                            6.293e-07
                        ],
                        [
                            0.73,
                            3.667e-07
                        ],
                        [
                            0.72,
                            1.063e-07
                        ],
                        [
                            0.71,
                            -1.209e-07
                        ],
                        [
                            0.7,
                            -2.9e-07
                        ],
                        [
                            0.69,
                            -3.912e-07
                        ],
                        [
                            0.68,
                            -4.298e-07
                        ],
                        [
                            0.67,
                            -4.222e-07
                        ],
                        [
                            0.66,
                            -3.875e-07
                        ],
                        [
                            0.65,
                            -3.409e-07
                        ],
                        [
                            0.64,
                            -2.931e-07
                        ],
                        [
                            0.63,
                            -2.496e-07
                        ],
                        [
                            0.62,
                            -2.127e-07
                        ],
                        [
                            0.61,
                            -1.826e-07
                        ],
                        [
                            0.6,
                            -1.582e-07
                        ],
                        [
                            0.59,
                            -1.389e-07
                        ],
                        [
                            0.58,
                            -1.235e-07
                        ],
                        [
                            0.57,
                            -1.113e-07
                        ],
                        [
                            0.56,
                            -1.016e-07
                        ],
                        [
                            0.55,
                            -9.373e-08
                        ],
                        [
                            0.54,
                            -8.749e-08
                        ],
                        [
                            0.53,
                            -8.235e-08
                        ],
                        [
                            0.52,
                            -7.818e-08
                        ],
                        [
                            0.51,
                            -7.475e-08
                        ],
                        [
                            0.5,
                            -7.186e-08
                        ],
                        [
                            0.49,
                            -6.937e-08
                        ],
                        [
                            0.48,
                            -6.703e-08
                        ],
                        [
                            0.47,
                            -6.494e-08
                        ],
                        [
                            0.46,
                            -6.3e-08
                        ],
                        [
                            0.45,
                            -6.119e-08
                        ],
                        [
                            0.44,
                            -5.946e-08
                        ],
                        [
                            0.43,
                            -5.79e-08
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.75,
                        2.021e-06
                    ]
                ],
                "reverse": [
                    [
                        0.68,
                        -4.298e-07
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            -1.75e-09,
                            -7.178e-10,
                            3.149e-10,
                            1.407e-09,
                            2.528e-09,
                            3.679e-09,
                            4.741e-09,
                            5.98e-09,
                            7.337e-09,
                            8.93e-09,
                            1.108e-08,
                            1.413e-08,
                            1.858e-08,
                            2.491e-08,
                            3.401e-08,
                            4.718e-08,
                            6.622e-08,
                            9.351e-08,
                            1.326e-07,
                            1.877e-07,
                            2.651e-07,
                            3.715e-07,
                            5.135e-07,
                            6.959e-07,
                            9.158e-07,
                            1.164e-06,
                            1.417e-06,
                            1.648e-06,
                            1.832e-06,
                            1.952e-06,
                            2.011e-06,
                            2.021e-06,
                            1.999e-06,
                            1.962e-06,
                            1.92e-06,
                            1.88e-06,
                            1.845e-06,
                            1.816e-06,
                            1.792e-06,
                            1.774e-06,
                            1.759e-06,
                            1.747e-06,
                            1.737e-06,
                            1.73e-06,
                            1.724e-06,
                            1.718e-06,
                            1.714e-06,
                            1.71e-06,
                            1.707e-06,
                            1.704e-06,
                            1.702e-06,
                            1.699e-06,
                            1.697e-06,
                            1.695e-06,
                            1.694e-06,
                            1.68e-06,
                            1.675e-06,
                            1.672e-06,
                            1.669e-06,
                            1.666e-06,
                            1.664e-06,
                            1.661e-06,
                            1.659e-06,
                            1.656e-06,
                            1.654e-06,
                            1.65e-06,
                            1.645e-06,
                            1.639e-06,
                            1.632e-06,
                            1.62e-06,
                            1.605e-06,
                            1.583e-06,
                            1.551e-06,
                            1.506e-06,
                            1.443e-06,
                            1.354e-06,
                            1.233e-06,
                            1.072e-06,
                            8.696e-07,
                            6.293e-07,
                            3.667e-07,
                            1.063e-07,
                            -1.209e-07,
                            -2.9e-07,
                            -3.912e-07,
                            -4.298e-07,
                            -4.222e-07,
                            -3.875e-07,
                            -3.409e-07,
                            -2.931e-07,
                            -2.496e-07,
                            -2.127e-07,
                            -1.826e-07,
                            -1.582e-07,
                            -1.389e-07,
                            -1.235e-07,
                            -1.113e-07,
                            -1.016e-07,
                            -9.373e-08,
                            -8.749e-08,
                            -8.235e-08,
                            -7.818e-08,
                            -7.475e-08,
                            -7.186e-08,
                            -6.937e-08,
                            -6.703e-08,
                            -6.494e-08,
                            -6.3e-08,
                            -6.119e-08,
                            -5.946e-08,
                            -5.79e-08
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "quasi-reversible"
                ],
                "e_half": [
                    0.715
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            -1.75e-09
                        ],
                        [
                            0.45,
                            -7.178e-10
                        ],
                        [
                            0.46,
                            3.149e-10
                        ],
                        [
                            0.47,
                            1.407e-09
                        ],
                        [
                            0.48,
                            2.528e-09
                        ],
                        [
                            0.49,
                            3.679e-09
                        ],
                        [
                            0.5,
                            4.741e-09
                        ],
                        [
                            0.51,
                            5.98e-09
                        ],
                        [
                            0.52,
                            7.337e-09
                        ],
                        [
                            0.53,
                            8.93e-09
                        ],
                        [
                            0.54,
                            1.108e-08
                        ],
                        [
                            0.55,
                            1.413e-08
                        ],
                        [
                            0.56,
                            1.858e-08
                        ],
                        [
                            0.57,
                            2.491e-08
                        ],
                        [
                            0.58,
                            3.401e-08
                        ],
                        [
                            0.59,
                            4.718e-08
                        ],
                        [
                            0.6,
                            6.622e-08
                        ],
                        [
                            0.61,
                            9.351e-08
                        ],
                        [
                            0.62,
                            1.326e-07
                        ],
                        [
                            0.63,
                            1.877e-07
                        ],
                        [
                            0.64,
                            2.651e-07
                        ],
                        [
                            0.65,
                            3.715e-07
                        ],
                        [
                            0.66,
                            5.135e-07
                        ],
                        [
                            0.67,
                            6.959e-07
                        ],
                        [
                            0.68,
                            9.158e-07
                        ],
                        [
                            0.69,
                            1.164e-06
                        ],
                        [
                            0.7,
                            1.417e-06
                        ],
                        [
                            0.71,
                            1.648e-06
                        ],
                        [
                            0.72,
                            1.832e-06
                        ],
                        [
                            0.73,
                            1.952e-06
                        ],
                        [
                            0.74,
                            2.011e-06
                        ],
                        [
                            0.75,
                            2.021e-06
                        ],
                        [
                            0.76,
                            1.999e-06
                        ],
                        [
                            0.77,
                            1.962e-06
                        ],
                        [
                            0.78,
                            1.92e-06
                        ],
                        [
                            0.79,
                            1.88e-06
                        ],
                        [
                            0.8,
                            1.845e-06
                        ],
                        [
                            0.81,
                            1.816e-06
                        ],
                        [
                            0.82,
                            1.792e-06
                        ],
                        [
                            0.83,
                            1.774e-06
                        ],
                        [
                            0.84,
                            1.759e-06
                        ],
                        [
                            0.85,
                            1.747e-06
                        ],
                        [
                            0.86,
                            1.737e-06
                        ],
                        [
                            0.87,
                            1.73e-06
                        ],
                        [
                            0.88,
                            1.724e-06
                        ],
                        [
                            0.89,
                            1.718e-06
                        ],
                        [
                            0.9,
                            1.714e-06
                        ],
                        [
                            0.91,
                            1.71e-06
                        ],
                        [
                            0.92,
                            1.707e-06
                        ],
                        [
                            0.93,
                            1.704e-06
                        ],
                        [
                            0.94,
                            1.702e-06
                        ],
                        [
                            0.95,
                            1.699e-06
                        ],
                        [
                            0.96,
                            1.697e-06
                        ],
                        [
                            0.97,
                            1.695e-06
                        ],
                        [
                            0.98,
                            1.694e-06
                        ],
                        [
                            0.98,
                            1.68e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.675e-06
                        ],
                        [
                            0.96,
                            1.672e-06
                        ],
                        [
                            0.95,
                            1.669e-06
                        ],
                        [
                            0.94,
                            1.666e-06
                        ],
                        [
                            0.93,
                            1.664e-06
                        ],
                        [
                            0.92,
                            1.661e-06
                        ],
                        [
                            0.91,
                            1.659e-06
                        ],
                        [
                            0.9,
                            1.656e-06
                        ],
                        [
                            0.89,
                            1.654e-06
                        ],
                        [
                            0.88,
                            1.65e-06
                        ],
                        [
                            0.87,
                            1.645e-06
                        ],
                        [
                            0.86,
                            1.639e-06
                        ],
                        [
                            0.85,
                            1.632e-06
                        ],
                        [
                            0.84,
                            1.62e-06
                        ],
                        [
                            0.83,
                            1.605e-06
                        ],
                        [
                            0.82,
                            1.583e-06
                        ],
                        [
                            0.81,
                            1.551e-06
                        ],
                        [
                            0.8,
                            1.506e-06
                        ],
                        [
                            0.79,
                            1.443e-06
                        ],
                        [
                            0.78,
                            1.354e-06
                        ],
                        [
                            0.77,
                            1.233e-06
                        ],
                        [
                            0.76,
                            1.072e-06
                        ],
                        [
                            0.75,
                            8.696e-07
                        ],
                        [
                            0.74,
                            6.293e-07
                        ],
                        [
                            0.73,
                            3.667e-07
                        ],
                        [
                            0.72,
                            1.063e-07
                        ],
                        [
                            0.71,
                            -1.209e-07
                        ],
                        [
                            0.7,
                            -2.9e-07
                        ],
                        [
                            0.69,
                            -3.912e-07
                        ],
                        [
                            0.68,
                            -4.298e-07
                        ],
                        [
                            0.67,
                            -4.222e-07
                        ],
                        [
                            0.66,
                            -3.875e-07
                        ],
                        [
                            0.65,
                            -3.409e-07
                        ],
                        [
                            0.64,
                            -2.931e-07
                        ],
                        [
                            0.63,
                            -2.496e-07
                        ],
                        [
                            0.62,
                            -2.127e-07
                        ],
                        [
                            0.61,
                            -1.826e-07
                        ],
                        [
                            0.6,
                            -1.582e-07
                        ],
                        [
                            0.59,
                            -1.389e-07
                        ],
                        [
                            0.58,
                            -1.235e-07
                        ],
                        [
                            0.57,
                            -1.113e-07
                        ],
                        [
                            0.56,
                            -1.016e-07
                        ],
                        [
                            0.55,
                            -9.373e-08
                        ],
                        [
                            0.54,
                            -8.749e-08
                        ],
                        [
                            0.53,
                            -8.235e-08
                        ],
                        [
                            0.52,
                            -7.818e-08
                        ],
                        [
                            0.51,
                            -7.475e-08
                        ],
                        [
                            0.5,
                            -7.186e-08
                        ],
                        [
                            0.49,
                            -6.937e-08
                        ],
                        [
                            0.48,
                            -6.703e-08
                        ],
                        [
                            0.47,
                            -6.494e-08
                        ],
                        [
                            0.46,
                            -6.3e-08
                        ],
                        [
                            0.45,
                            -6.119e-08
                        ],
                        [
                            0.44,
                            -5.946e-08
                        ],
                        [
                            0.43,
                            -5.79e-08
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 2.021e-06,
                "current_anodic": -4.298e-07,
                "diffusion": 1.451e-09
            }
        },
        {
            "_id": "040802dc-92a1-4de9-bcd6-b130b723984d",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:55.639210",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv07_12_22_47.bin",
                "header": "CV cycle01_cv07_12_22_47",
                "note": "",
                "date_recorded": "2024-01-31T12:23:17",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.05,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            -4.586e-09
                        ],
                        [
                            0.45,
                            -3.633e-09
                        ],
                        [
                            0.46,
                            -2.641e-09
                        ],
                        [
                            0.47,
                            -1.698e-09
                        ],
                        [
                            0.48,
                            -7.251e-10
                        ],
                        [
                            0.49,
                            3.655e-10
                        ],
                        [
                            0.5,
                            1.436e-09
                        ],
                        [
                            0.51,
                            2.625e-09
                        ],
                        [
                            0.52,
                            3.824e-09
                        ],
                        [
                            0.53,
                            5.337e-09
                        ],
                        [
                            0.54,
                            7.243e-09
                        ],
                        [
                            0.55,
                            9.867e-09
                        ],
                        [
                            0.56,
                            1.34e-08
                        ],
                        [
                            0.57,
                            1.844e-08
                        ],
                        [
                            0.58,
                            2.572e-08
                        ],
                        [
                            0.59,
                            3.622e-08
                        ],
                        [
                            0.6,
                            5.122e-08
                        ],
                        [
                            0.61,
                            7.28e-08
                        ],
                        [
                            0.62,
                            1.036e-07
                        ],
                        [
                            0.63,
                            1.472e-07
                        ],
                        [
                            0.64,
                            2.084e-07
                        ],
                        [
                            0.65,
                            2.929e-07
                        ],
                        [
                            0.66,
                            4.065e-07
                        ],
                        [
                            0.67,
                            5.538e-07
                        ],
                        [
                            0.68,
                            7.346e-07
                        ],
                        [
                            0.69,
                            9.426e-07
                        ],
                        [
                            0.7,
                            1.162e-06
                        ],
                        [
                            0.71,
                            1.371e-06
                        ],
                        [
                            0.72,
                            1.549e-06
                        ],
                        [
                            0.73,
                            1.679e-06
                        ],
                        [
                            0.74,
                            1.762e-06
                        ],
                        [
                            0.75,
                            1.803e-06
                        ],
                        [
                            0.76,
                            1.815e-06
                        ],
                        [
                            0.77,
                            1.809e-06
                        ],
                        [
                            0.78,
                            1.794e-06
                        ],
                        [
                            0.79,
                            1.777e-06
                        ],
                        [
                            0.8,
                            1.761e-06
                        ],
                        [
                            0.81,
                            1.746e-06
                        ],
                        [
                            0.82,
                            1.733e-06
                        ],
                        [
                            0.83,
                            1.723e-06
                        ],
                        [
                            0.84,
                            1.714e-06
                        ],
                        [
                            0.85,
                            1.707e-06
                        ],
                        [
                            0.86,
                            1.701e-06
                        ],
                        [
                            0.87,
                            1.696e-06
                        ],
                        [
                            0.88,
                            1.692e-06
                        ],
                        [
                            0.89,
                            1.688e-06
                        ],
                        [
                            0.9,
                            1.685e-06
                        ],
                        [
                            0.91,
                            1.683e-06
                        ],
                        [
                            0.92,
                            1.68e-06
                        ],
                        [
                            0.93,
                            1.678e-06
                        ],
                        [
                            0.94,
                            1.676e-06
                        ],
                        [
                            0.95,
                            1.674e-06
                        ],
                        [
                            0.96,
                            1.672e-06
                        ],
                        [
                            0.97,
                            1.671e-06
                        ],
                        [
                            0.98,
                            1.669e-06
                        ],
                        [
                            0.98,
                            1.66e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.655e-06
                        ],
                        [
                            0.96,
                            1.652e-06
                        ],
                        [
                            0.95,
                            1.65e-06
                        ],
                        [
                            0.94,
                            1.647e-06
                        ],
                        [
                            0.93,
                            1.645e-06
                        ],
                        [
                            0.92,
                            1.643e-06
                        ],
                        [
                            0.91,
                            1.64e-06
                        ],
                        [
                            0.9,
                            1.638e-06
                        ],
                        [
                            0.89,
                            1.635e-06
                        ],
                        [
                            0.88,
                            1.632e-06
                        ],
                        [
                            0.87,
                            1.628e-06
                        ],
                        [
                            0.86,
                            1.623e-06
                        ],
                        [
                            0.85,
                            1.616e-06
                        ],
                        [
                            0.84,
                            1.607e-06
                        ],
                        [
                            0.83,
                            1.593e-06
                        ],
                        [
                            0.82,
                            1.574e-06
                        ],
                        [
                            0.81,
                            1.546e-06
                        ],
                        [
                            0.8,
                            1.507e-06
                        ],
                        [
                            0.79,
                            1.452e-06
                        ],
                        [
                            0.78,
                            1.375e-06
                        ],
                        [
                            0.77,
                            1.27e-06
                        ],
                        [
                            0.76,
                            1.131e-06
                        ],
                        [
                            0.75,
                            9.577e-07
                        ],
                        [
                            0.74,
                            7.52e-07
                        ],
                        [
                            0.73,
                            5.275e-07
                        ],
                        [
                            0.72,
                            3.044e-07
                        ],
                        [
                            0.71,
                            1.067e-07
                        ],
                        [
                            0.7,
                            -4.606e-08
                        ],
                        [
                            0.69,
                            -1.472e-07
                        ],
                        [
                            0.68,
                            -2.006e-07
                        ],
                        [
                            0.67,
                            -2.179e-07
                        ],
                        [
                            0.66,
                            -2.133e-07
                        ],
                        [
                            0.65,
                            -1.968e-07
                        ],
                        [
                            0.64,
                            -1.767e-07
                        ],
                        [
                            0.63,
                            -1.568e-07
                        ],
                        [
                            0.62,
                            -1.39e-07
                        ],
                        [
                            0.61,
                            -1.239e-07
                        ],
                        [
                            0.6,
                            -1.112e-07
                        ],
                        [
                            0.59,
                            -1.009e-07
                        ],
                        [
                            0.58,
                            -9.231e-08
                        ],
                        [
                            0.57,
                            -8.524e-08
                        ],
                        [
                            0.56,
                            -7.943e-08
                        ],
                        [
                            0.55,
                            -7.456e-08
                        ],
                        [
                            0.54,
                            -7.053e-08
                        ],
                        [
                            0.53,
                            -6.71e-08
                        ],
                        [
                            0.52,
                            -6.417e-08
                        ],
                        [
                            0.51,
                            -6.176e-08
                        ],
                        [
                            0.5,
                            -5.956e-08
                        ],
                        [
                            0.49,
                            -5.76e-08
                        ],
                        [
                            0.48,
                            -5.587e-08
                        ],
                        [
                            0.47,
                            -5.422e-08
                        ],
                        [
                            0.46,
                            -5.26e-08
                        ],
                        [
                            0.45,
                            -5.11e-08
                        ],
                        [
                            0.44,
                            -4.965e-08
                        ],
                        [
                            0.43,
                            -4.835e-08
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.76,
                        1.815e-06
                    ]
                ],
                "reverse": [
                    [
                        0.67,
                        -2.179e-07
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            -4.586e-09,
                            -3.633e-09,
                            -2.641e-09,
                            -1.698e-09,
                            -7.251e-10,
                            3.655e-10,
                            1.436e-09,
                            2.625e-09,
                            3.824e-09,
                            5.337e-09,
                            7.243e-09,
                            9.867e-09,
                            1.34e-08,
                            1.844e-08,
                            2.572e-08,
                            3.622e-08,
                            5.122e-08,
                            7.28e-08,
                            1.036e-07,
                            1.472e-07,
                            2.084e-07,
                            2.929e-07,
                            4.065e-07,
                            5.538e-07,
                            7.346e-07,
                            9.426e-07,
                            1.162e-06,
                            1.371e-06,
                            1.549e-06,
                            1.679e-06,
                            1.762e-06,
                            1.803e-06,
                            1.815e-06,
                            1.809e-06,
                            1.794e-06,
                            1.777e-06,
                            1.761e-06,
                            1.746e-06,
                            1.733e-06,
                            1.723e-06,
                            1.714e-06,
                            1.707e-06,
                            1.701e-06,
                            1.696e-06,
                            1.692e-06,
                            1.688e-06,
                            1.685e-06,
                            1.683e-06,
                            1.68e-06,
                            1.678e-06,
                            1.676e-06,
                            1.674e-06,
                            1.672e-06,
                            1.671e-06,
                            1.669e-06,
                            1.66e-06,
                            1.655e-06,
                            1.652e-06,
                            1.65e-06,
                            1.647e-06,
                            1.645e-06,
                            1.643e-06,
                            1.64e-06,
                            1.638e-06,
                            1.635e-06,
                            1.632e-06,
                            1.628e-06,
                            1.623e-06,
                            1.616e-06,
                            1.607e-06,
                            1.593e-06,
                            1.574e-06,
                            1.546e-06,
                            1.507e-06,
                            1.452e-06,
                            1.375e-06,
                            1.27e-06,
                            1.131e-06,
                            9.577e-07,
                            7.52e-07,
                            5.275e-07,
                            3.044e-07,
                            1.067e-07,
                            -4.606e-08,
                            -1.472e-07,
                            -2.006e-07,
                            -2.179e-07,
                            -2.133e-07,
                            -1.968e-07,
                            -1.767e-07,
                            -1.568e-07,
                            -1.39e-07,
                            -1.239e-07,
                            -1.112e-07,
                            -1.009e-07,
                            -9.231e-08,
                            -8.524e-08,
                            -7.943e-08,
                            -7.456e-08,
                            -7.053e-08,
                            -6.71e-08,
                            -6.417e-08,
                            -6.176e-08,
                            -5.956e-08,
                            -5.76e-08,
                            -5.587e-08,
                            -5.422e-08,
                            -5.26e-08,
                            -5.11e-08,
                            -4.965e-08,
                            -4.835e-08
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "quasi-reversible"
                ],
                "e_half": [
                    0.715
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            -4.586e-09
                        ],
                        [
                            0.45,
                            -3.633e-09
                        ],
                        [
                            0.46,
                            -2.641e-09
                        ],
                        [
                            0.47,
                            -1.698e-09
                        ],
                        [
                            0.48,
                            -7.251e-10
                        ],
                        [
                            0.49,
                            3.655e-10
                        ],
                        [
                            0.5,
                            1.436e-09
                        ],
                        [
                            0.51,
                            2.625e-09
                        ],
                        [
                            0.52,
                            3.824e-09
                        ],
                        [
                            0.53,
                            5.337e-09
                        ],
                        [
                            0.54,
                            7.243e-09
                        ],
                        [
                            0.55,
                            9.867e-09
                        ],
                        [
                            0.56,
                            1.34e-08
                        ],
                        [
                            0.57,
                            1.844e-08
                        ],
                        [
                            0.58,
                            2.572e-08
                        ],
                        [
                            0.59,
                            3.622e-08
                        ],
                        [
                            0.6,
                            5.122e-08
                        ],
                        [
                            0.61,
                            7.28e-08
                        ],
                        [
                            0.62,
                            1.036e-07
                        ],
                        [
                            0.63,
                            1.472e-07
                        ],
                        [
                            0.64,
                            2.084e-07
                        ],
                        [
                            0.65,
                            2.929e-07
                        ],
                        [
                            0.66,
                            4.065e-07
                        ],
                        [
                            0.67,
                            5.538e-07
                        ],
                        [
                            0.68,
                            7.346e-07
                        ],
                        [
                            0.69,
                            9.426e-07
                        ],
                        [
                            0.7,
                            1.162e-06
                        ],
                        [
                            0.71,
                            1.371e-06
                        ],
                        [
                            0.72,
                            1.549e-06
                        ],
                        [
                            0.73,
                            1.679e-06
                        ],
                        [
                            0.74,
                            1.762e-06
                        ],
                        [
                            0.75,
                            1.803e-06
                        ],
                        [
                            0.76,
                            1.815e-06
                        ],
                        [
                            0.77,
                            1.809e-06
                        ],
                        [
                            0.78,
                            1.794e-06
                        ],
                        [
                            0.79,
                            1.777e-06
                        ],
                        [
                            0.8,
                            1.761e-06
                        ],
                        [
                            0.81,
                            1.746e-06
                        ],
                        [
                            0.82,
                            1.733e-06
                        ],
                        [
                            0.83,
                            1.723e-06
                        ],
                        [
                            0.84,
                            1.714e-06
                        ],
                        [
                            0.85,
                            1.707e-06
                        ],
                        [
                            0.86,
                            1.701e-06
                        ],
                        [
                            0.87,
                            1.696e-06
                        ],
                        [
                            0.88,
                            1.692e-06
                        ],
                        [
                            0.89,
                            1.688e-06
                        ],
                        [
                            0.9,
                            1.685e-06
                        ],
                        [
                            0.91,
                            1.683e-06
                        ],
                        [
                            0.92,
                            1.68e-06
                        ],
                        [
                            0.93,
                            1.678e-06
                        ],
                        [
                            0.94,
                            1.676e-06
                        ],
                        [
                            0.95,
                            1.674e-06
                        ],
                        [
                            0.96,
                            1.672e-06
                        ],
                        [
                            0.97,
                            1.671e-06
                        ],
                        [
                            0.98,
                            1.669e-06
                        ],
                        [
                            0.98,
                            1.66e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.655e-06
                        ],
                        [
                            0.96,
                            1.652e-06
                        ],
                        [
                            0.95,
                            1.65e-06
                        ],
                        [
                            0.94,
                            1.647e-06
                        ],
                        [
                            0.93,
                            1.645e-06
                        ],
                        [
                            0.92,
                            1.643e-06
                        ],
                        [
                            0.91,
                            1.64e-06
                        ],
                        [
                            0.9,
                            1.638e-06
                        ],
                        [
                            0.89,
                            1.635e-06
                        ],
                        [
                            0.88,
                            1.632e-06
                        ],
                        [
                            0.87,
                            1.628e-06
                        ],
                        [
                            0.86,
                            1.623e-06
                        ],
                        [
                            0.85,
                            1.616e-06
                        ],
                        [
                            0.84,
                            1.607e-06
                        ],
                        [
                            0.83,
                            1.593e-06
                        ],
                        [
                            0.82,
                            1.574e-06
                        ],
                        [
                            0.81,
                            1.546e-06
                        ],
                        [
                            0.8,
                            1.507e-06
                        ],
                        [
                            0.79,
                            1.452e-06
                        ],
                        [
                            0.78,
                            1.375e-06
                        ],
                        [
                            0.77,
                            1.27e-06
                        ],
                        [
                            0.76,
                            1.131e-06
                        ],
                        [
                            0.75,
                            9.577e-07
                        ],
                        [
                            0.74,
                            7.52e-07
                        ],
                        [
                            0.73,
                            5.275e-07
                        ],
                        [
                            0.72,
                            3.044e-07
                        ],
                        [
                            0.71,
                            1.067e-07
                        ],
                        [
                            0.7,
                            -4.606e-08
                        ],
                        [
                            0.69,
                            -1.472e-07
                        ],
                        [
                            0.68,
                            -2.006e-07
                        ],
                        [
                            0.67,
                            -2.179e-07
                        ],
                        [
                            0.66,
                            -2.133e-07
                        ],
                        [
                            0.65,
                            -1.968e-07
                        ],
                        [
                            0.64,
                            -1.767e-07
                        ],
                        [
                            0.63,
                            -1.568e-07
                        ],
                        [
                            0.62,
                            -1.39e-07
                        ],
                        [
                            0.61,
                            -1.239e-07
                        ],
                        [
                            0.6,
                            -1.112e-07
                        ],
                        [
                            0.59,
                            -1.009e-07
                        ],
                        [
                            0.58,
                            -9.231e-08
                        ],
                        [
                            0.57,
                            -8.524e-08
                        ],
                        [
                            0.56,
                            -7.943e-08
                        ],
                        [
                            0.55,
                            -7.456e-08
                        ],
                        [
                            0.54,
                            -7.053e-08
                        ],
                        [
                            0.53,
                            -6.71e-08
                        ],
                        [
                            0.52,
                            -6.417e-08
                        ],
                        [
                            0.51,
                            -6.176e-08
                        ],
                        [
                            0.5,
                            -5.956e-08
                        ],
                        [
                            0.49,
                            -5.76e-08
                        ],
                        [
                            0.48,
                            -5.587e-08
                        ],
                        [
                            0.47,
                            -5.422e-08
                        ],
                        [
                            0.46,
                            -5.26e-08
                        ],
                        [
                            0.45,
                            -5.11e-08
                        ],
                        [
                            0.44,
                            -4.965e-08
                        ],
                        [
                            0.43,
                            -4.835e-08
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 1.815e-06,
                "current_anodic": -2.179e-07,
                "diffusion": 1.451e-09
            }
        },
        {
            "_id": "7d081550-0c53-4bae-8f37-c45e60291240",
            "mol_id": "05MYHH",
            "submission_info": {
                "processing_id": "None",
                "source": "d3tales_robot",
                "author": "d3tales_robot",
                "author_email": "d3tales@gmail.com",
                "upload_time": "2024-01-31T12:36:57.590192",
                "file_type": "txt",
                "data_category": "experimentation",
                "data_type": "cv"
            },
            "data": {
                "file_name": "c/users/lab/d3talesrobotics/data/basiccvtest_meept/20240131/exp01_05myhh\\cycle01_cv08_12_23_22.bin",
                "header": "CV cycle01_cv08_12_23_22",
                "note": "",
                "date_recorded": "2024-01-31T12:24:14",
                "conditions": {
                    "data_source": "cv",
                    "scan_rate": {
                        "value": 0.025,
                        "unit": "V/s"
                    },
                    "num_scans": 2,
                    "initial_potential": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "high_e": {
                        "value": 0.985,
                        "unit": "V"
                    },
                    "low_e": {
                        "value": 0.425,
                        "unit": "V"
                    },
                    "comp_r": {
                        "value": 430.4,
                        "unit": "ohm"
                    },
                    "working_electrode": "screen-printed micro-electrode (working)",
                    "counter_electrode": "screen-printed micro-electrode (counter)",
                    "reference_electrode": "Ag/Ag+ wire electrode",
                    "solvent": [
                        {
                            "name": "Acetonitrile",
                            "purity": ""
                        }
                    ],
                    "instrument": "robotics_potentiostat_A_01",
                    "working_electrode_surface_area": {
                        "value": 0.031415926535897934
                    },
                    "redox_mol_concentration": {
                        "value": 0.01999245136525313,
                        "unit": "molar"
                    },
                    "temperature": {
                        "value": 1.0,
                        "unit": "dimensionless"
                    },
                    "experiment_run_id": "43251109-d74a-4816-9d32-b24af7b7aaac"
                },
                "segment": 2,
                "sample_interval": {
                    "value": 0.01,
                    "unit": "V"
                },
                "quiet_time": {
                    "value": 2.0,
                    "unit": "sec"
                },
                "sensitivity": {
                    "value": 1e-05,
                    "unit": "A/V"
                },
                "peak_potential": {
                    "value": 0.425,
                    "unit": "V"
                },
                "scan_data": [
                    [
                        [
                            0.44,
                            -7.259e-09
                        ],
                        [
                            0.45,
                            -6.276e-09
                        ],
                        [
                            0.46,
                            -5.372e-09
                        ],
                        [
                            0.47,
                            -4.41e-09
                        ],
                        [
                            0.48,
                            -3.447e-09
                        ],
                        [
                            0.49,
                            -2.425e-09
                        ],
                        [
                            0.5,
                            -1.393e-09
                        ],
                        [
                            0.51,
                            -2.437e-10
                        ],
                        [
                            0.52,
                            9.845e-10
                        ],
                        [
                            0.53,
                            2.448e-09
                        ],
                        [
                            0.54,
                            4.148e-09
                        ],
                        [
                            0.55,
                            6.379e-09
                        ],
                        [
                            0.56,
                            9.267e-09
                        ],
                        [
                            0.57,
                            1.328e-08
                        ],
                        [
                            0.58,
                            1.881e-08
                        ],
                        [
                            0.59,
                            2.671e-08
                        ],
                        [
                            0.6,
                            3.795e-08
                        ],
                        [
                            0.61,
                            5.407e-08
                        ],
                        [
                            0.62,
                            7.707e-08
                        ],
                        [
                            0.63,
                            1.098e-07
                        ],
                        [
                            0.64,
                            1.558e-07
                        ],
                        [
                            0.65,
                            2.194e-07
                        ],
                        [
                            0.66,
                            3.057e-07
                        ],
                        [
                            0.67,
                            4.187e-07
                        ],
                        [
                            0.68,
                            5.599e-07
                        ],
                        [
                            0.69,
                            7.26e-07
                        ],
                        [
                            0.7,
                            9.082e-07
                        ],
                        [
                            0.71,
                            1.091e-06
                        ],
                        [
                            0.72,
                            1.258e-06
                        ],
                        [
                            0.73,
                            1.396e-06
                        ],
                        [
                            0.74,
                            1.501e-06
                        ],
                        [
                            0.75,
                            1.574e-06
                        ],
                        [
                            0.76,
                            1.62e-06
                        ],
                        [
                            0.77,
                            1.648e-06
                        ],
                        [
                            0.78,
                            1.663e-06
                        ],
                        [
                            0.79,
                            1.67e-06
                        ],
                        [
                            0.8,
                            1.672e-06
                        ],
                        [
                            0.81,
                            1.671e-06
                        ],
                        [
                            0.82,
                            1.67e-06
                        ],
                        [
                            0.83,
                            1.667e-06
                        ],
                        [
                            0.84,
                            1.664e-06
                        ],
                        [
                            0.85,
                            1.662e-06
                        ],
                        [
                            0.86,
                            1.659e-06
                        ],
                        [
                            0.87,
                            1.656e-06
                        ],
                        [
                            0.88,
                            1.653e-06
                        ],
                        [
                            0.89,
                            1.651e-06
                        ],
                        [
                            0.9,
                            1.649e-06
                        ],
                        [
                            0.91,
                            1.647e-06
                        ],
                        [
                            0.92,
                            1.645e-06
                        ],
                        [
                            0.93,
                            1.643e-06
                        ],
                        [
                            0.94,
                            1.642e-06
                        ],
                        [
                            0.95,
                            1.64e-06
                        ],
                        [
                            0.96,
                            1.639e-06
                        ],
                        [
                            0.97,
                            1.638e-06
                        ],
                        [
                            0.98,
                            1.637e-06
                        ],
                        [
                            0.98,
                            1.631e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.628e-06
                        ],
                        [
                            0.96,
                            1.626e-06
                        ],
                        [
                            0.95,
                            1.624e-06
                        ],
                        [
                            0.94,
                            1.622e-06
                        ],
                        [
                            0.93,
                            1.62e-06
                        ],
                        [
                            0.92,
                            1.618e-06
                        ],
                        [
                            0.91,
                            1.616e-06
                        ],
                        [
                            0.9,
                            1.615e-06
                        ],
                        [
                            0.89,
                            1.612e-06
                        ],
                        [
                            0.88,
                            1.61e-06
                        ],
                        [
                            0.87,
                            1.607e-06
                        ],
                        [
                            0.86,
                            1.603e-06
                        ],
                        [
                            0.85,
                            1.598e-06
                        ],
                        [
                            0.84,
                            1.592e-06
                        ],
                        [
                            0.83,
                            1.583e-06
                        ],
                        [
                            0.82,
                            1.57e-06
                        ],
                        [
                            0.81,
                            1.551e-06
                        ],
                        [
                            0.8,
                            1.524e-06
                        ],
                        [
                            0.79,
                            1.484e-06
                        ],
                        [
                            0.78,
                            1.428e-06
                        ],
                        [
                            0.77,
                            1.352e-06
                        ],
                        [
                            0.76,
                            1.251e-06
                        ],
                        [
                            0.75,
                            1.121e-06
                        ],
                        [
                            0.74,
                            9.637e-07
                        ],
                        [
                            0.73,
                            7.846e-07
                        ],
                        [
                            0.72,
                            5.966e-07
                        ],
                        [
                            0.71,
                            4.165e-07
                        ],
                        [
                            0.7,
                            2.603e-07
                        ],
                        [
                            0.69,
                            1.374e-07
                        ],
                        [
                            0.68,
                            5.015e-08
                        ],
                        [
                            0.67,
                            -6.945e-09
                        ],
                        [
                            0.66,
                            -4.123e-08
                        ],
                        [
                            0.65,
                            -5.962e-08
                        ],
                        [
                            0.64,
                            -6.839e-08
                        ],
                        [
                            0.63,
                            -7.124e-08
                        ],
                        [
                            0.62,
                            -7.096e-08
                        ],
                        [
                            0.61,
                            -6.909e-08
                        ],
                        [
                            0.6,
                            -6.612e-08
                        ],
                        [
                            0.59,
                            -6.289e-08
                        ],
                        [
                            0.58,
                            -5.962e-08
                        ],
                        [
                            0.57,
                            -5.649e-08
                        ],
                        [
                            0.56,
                            -5.355e-08
                        ],
                        [
                            0.55,
                            -5.077e-08
                        ],
                        [
                            0.54,
                            -4.821e-08
                        ],
                        [
                            0.53,
                            -4.594e-08
                        ],
                        [
                            0.52,
                            -4.385e-08
                        ],
                        [
                            0.51,
                            -4.204e-08
                        ],
                        [
                            0.5,
                            -4.035e-08
                        ],
                        [
                            0.49,
                            -3.882e-08
                        ],
                        [
                            0.48,
                            -3.732e-08
                        ],
                        [
                            0.47,
                            -3.595e-08
                        ],
                        [
                            0.46,
                            -3.466e-08
                        ],
                        [
                            0.45,
                            -3.341e-08
                        ],
                        [
                            0.44,
                            -3.228e-08
                        ],
                        [
                            0.43,
                            -3.119e-08
                        ]
                    ]
                ],
                "forward": [
                    [
                        0.8,
                        1.672e-06
                    ]
                ],
                "reverse": [
                    [
                        0.63,
                        -7.124e-08
                    ]
                ],
                "plot_data": [
                    {
                        "x": [
                            0.44,
                            0.45,
                            0.46,
                            0.47,
                            0.48,
                            0.49,
                            0.5,
                            0.51,
                            0.52,
                            0.53,
                            0.54,
                            0.55,
                            0.56,
                            0.57,
                            0.58,
                            0.59,
                            0.6,
                            0.61,
                            0.62,
                            0.63,
                            0.64,
                            0.65,
                            0.66,
                            0.67,
                            0.68,
                            0.69,
                            0.7,
                            0.71,
                            0.72,
                            0.73,
                            0.74,
                            0.75,
                            0.76,
                            0.77,
                            0.78,
                            0.79,
                            0.8,
                            0.81,
                            0.82,
                            0.83,
                            0.84,
                            0.85,
                            0.86,
                            0.87,
                            0.88,
                            0.89,
                            0.9,
                            0.91,
                            0.92,
                            0.93,
                            0.94,
                            0.95,
                            0.96,
                            0.97,
                            0.98,
                            0.98,
                            0.97,
                            0.96,
                            0.95,
                            0.94,
                            0.93,
                            0.92,
                            0.91,
                            0.9,
                            0.89,
                            0.88,
                            0.87,
                            0.86,
                            0.85,
                            0.84,
                            0.83,
                            0.82,
                            0.81,
                            0.8,
                            0.79,
                            0.78,
                            0.77,
                            0.76,
                            0.75,
                            0.74,
                            0.73,
                            0.72,
                            0.71,
                            0.7,
                            0.69,
                            0.68,
                            0.67,
                            0.66,
                            0.65,
                            0.64,
                            0.63,
                            0.62,
                            0.61,
                            0.6,
                            0.59,
                            0.58,
                            0.57,
                            0.56,
                            0.55,
                            0.54,
                            0.53,
                            0.52,
                            0.51,
                            0.5,
                            0.49,
                            0.48,
                            0.47,
                            0.46,
                            0.45,
                            0.44,
                            0.43
                        ],
                        "y": [
                            -7.259e-09,
                            -6.276e-09,
                            -5.372e-09,
                            -4.41e-09,
                            -3.447e-09,
                            -2.425e-09,
                            -1.393e-09,
                            -2.437e-10,
                            9.845e-10,
                            2.448e-09,
                            4.148e-09,
                            6.379e-09,
                            9.267e-09,
                            1.328e-08,
                            1.881e-08,
                            2.671e-08,
                            3.795e-08,
                            5.407e-08,
                            7.707e-08,
                            1.098e-07,
                            1.558e-07,
                            2.194e-07,
                            3.057e-07,
                            4.187e-07,
                            5.599e-07,
                            7.26e-07,
                            9.082e-07,
                            1.091e-06,
                            1.258e-06,
                            1.396e-06,
                            1.501e-06,
                            1.574e-06,
                            1.62e-06,
                            1.648e-06,
                            1.663e-06,
                            1.67e-06,
                            1.672e-06,
                            1.671e-06,
                            1.67e-06,
                            1.667e-06,
                            1.664e-06,
                            1.662e-06,
                            1.659e-06,
                            1.656e-06,
                            1.653e-06,
                            1.651e-06,
                            1.649e-06,
                            1.647e-06,
                            1.645e-06,
                            1.643e-06,
                            1.642e-06,
                            1.64e-06,
                            1.639e-06,
                            1.638e-06,
                            1.637e-06,
                            1.631e-06,
                            1.628e-06,
                            1.626e-06,
                            1.624e-06,
                            1.622e-06,
                            1.62e-06,
                            1.618e-06,
                            1.616e-06,
                            1.615e-06,
                            1.612e-06,
                            1.61e-06,
                            1.607e-06,
                            1.603e-06,
                            1.598e-06,
                            1.592e-06,
                            1.583e-06,
                            1.57e-06,
                            1.551e-06,
                            1.524e-06,
                            1.484e-06,
                            1.428e-06,
                            1.352e-06,
                            1.251e-06,
                            1.121e-06,
                            9.637e-07,
                            7.846e-07,
                            5.966e-07,
                            4.165e-07,
                            2.603e-07,
                            1.374e-07,
                            5.015e-08,
                            -6.945e-09,
                            -4.123e-08,
                            -5.962e-08,
                            -6.839e-08,
                            -7.124e-08,
                            -7.096e-08,
                            -6.909e-08,
                            -6.612e-08,
                            -6.289e-08,
                            -5.962e-08,
                            -5.649e-08,
                            -5.355e-08,
                            -5.077e-08,
                            -4.821e-08,
                            -4.594e-08,
                            -4.385e-08,
                            -4.204e-08,
                            -4.035e-08,
                            -3.882e-08,
                            -3.732e-08,
                            -3.595e-08,
                            -3.466e-08,
                            -3.341e-08,
                            -3.228e-08,
                            -3.119e-08
                        ],
                        "modeplot_data": "lines",
                        "name": "cv",
                        "line": {
                            "color": "#003396",
                            "width": 3
                        }
                    }
                ],
                "reversibility": [
                    "quasi-reversible"
                ],
                "e_half": [
                    0.715
                ],
                "peak_splittings": [],
                "middle_sweep": [
                    [
                        [
                            0.44,
                            -7.259e-09
                        ],
                        [
                            0.45,
                            -6.276e-09
                        ],
                        [
                            0.46,
                            -5.372e-09
                        ],
                        [
                            0.47,
                            -4.41e-09
                        ],
                        [
                            0.48,
                            -3.447e-09
                        ],
                        [
                            0.49,
                            -2.425e-09
                        ],
                        [
                            0.5,
                            -1.393e-09
                        ],
                        [
                            0.51,
                            -2.437e-10
                        ],
                        [
                            0.52,
                            9.845e-10
                        ],
                        [
                            0.53,
                            2.448e-09
                        ],
                        [
                            0.54,
                            4.148e-09
                        ],
                        [
                            0.55,
                            6.379e-09
                        ],
                        [
                            0.56,
                            9.267e-09
                        ],
                        [
                            0.57,
                            1.328e-08
                        ],
                        [
                            0.58,
                            1.881e-08
                        ],
                        [
                            0.59,
                            2.671e-08
                        ],
                        [
                            0.6,
                            3.795e-08
                        ],
                        [
                            0.61,
                            5.407e-08
                        ],
                        [
                            0.62,
                            7.707e-08
                        ],
                        [
                            0.63,
                            1.098e-07
                        ],
                        [
                            0.64,
                            1.558e-07
                        ],
                        [
                            0.65,
                            2.194e-07
                        ],
                        [
                            0.66,
                            3.057e-07
                        ],
                        [
                            0.67,
                            4.187e-07
                        ],
                        [
                            0.68,
                            5.599e-07
                        ],
                        [
                            0.69,
                            7.26e-07
                        ],
                        [
                            0.7,
                            9.082e-07
                        ],
                        [
                            0.71,
                            1.091e-06
                        ],
                        [
                            0.72,
                            1.258e-06
                        ],
                        [
                            0.73,
                            1.396e-06
                        ],
                        [
                            0.74,
                            1.501e-06
                        ],
                        [
                            0.75,
                            1.574e-06
                        ],
                        [
                            0.76,
                            1.62e-06
                        ],
                        [
                            0.77,
                            1.648e-06
                        ],
                        [
                            0.78,
                            1.663e-06
                        ],
                        [
                            0.79,
                            1.67e-06
                        ],
                        [
                            0.8,
                            1.672e-06
                        ],
                        [
                            0.81,
                            1.671e-06
                        ],
                        [
                            0.82,
                            1.67e-06
                        ],
                        [
                            0.83,
                            1.667e-06
                        ],
                        [
                            0.84,
                            1.664e-06
                        ],
                        [
                            0.85,
                            1.662e-06
                        ],
                        [
                            0.86,
                            1.659e-06
                        ],
                        [
                            0.87,
                            1.656e-06
                        ],
                        [
                            0.88,
                            1.653e-06
                        ],
                        [
                            0.89,
                            1.651e-06
                        ],
                        [
                            0.9,
                            1.649e-06
                        ],
                        [
                            0.91,
                            1.647e-06
                        ],
                        [
                            0.92,
                            1.645e-06
                        ],
                        [
                            0.93,
                            1.643e-06
                        ],
                        [
                            0.94,
                            1.642e-06
                        ],
                        [
                            0.95,
                            1.64e-06
                        ],
                        [
                            0.96,
                            1.639e-06
                        ],
                        [
                            0.97,
                            1.638e-06
                        ],
                        [
                            0.98,
                            1.637e-06
                        ],
                        [
                            0.98,
                            1.631e-06
                        ]
                    ],
                    [
                        [
                            0.97,
                            1.628e-06
                        ],
                        [
                            0.96,
                            1.626e-06
                        ],
                        [
                            0.95,
                            1.624e-06
                        ],
                        [
                            0.94,
                            1.622e-06
                        ],
                        [
                            0.93,
                            1.62e-06
                        ],
                        [
                            0.92,
                            1.618e-06
                        ],
                        [
                            0.91,
                            1.616e-06
                        ],
                        [
                            0.9,
                            1.615e-06
                        ],
                        [
                            0.89,
                            1.612e-06
                        ],
                        [
                            0.88,
                            1.61e-06
                        ],
                        [
                            0.87,
                            1.607e-06
                        ],
                        [
                            0.86,
                            1.603e-06
                        ],
                        [
                            0.85,
                            1.598e-06
                        ],
                        [
                            0.84,
                            1.592e-06
                        ],
                        [
                            0.83,
                            1.583e-06
                        ],
                        [
                            0.82,
                            1.57e-06
                        ],
                        [
                            0.81,
                            1.551e-06
                        ],
                        [
                            0.8,
                            1.524e-06
                        ],
                        [
                            0.79,
                            1.484e-06
                        ],
                        [
                            0.78,
                            1.428e-06
                        ],
                        [
                            0.77,
                            1.352e-06
                        ],
                        [
                            0.76,
                            1.251e-06
                        ],
                        [
                            0.75,
                            1.121e-06
                        ],
                        [
                            0.74,
                            9.637e-07
                        ],
                        [
                            0.73,
                            7.846e-07
                        ],
                        [
                            0.72,
                            5.966e-07
                        ],
                        [
                            0.71,
                            4.165e-07
                        ],
                        [
                            0.7,
                            2.603e-07
                        ],
                        [
                            0.69,
                            1.374e-07
                        ],
                        [
                            0.68,
                            5.015e-08
                        ],
                        [
                            0.67,
                            -6.945e-09
                        ],
                        [
                            0.66,
                            -4.123e-08
                        ],
                        [
                            0.65,
                            -5.962e-08
                        ],
                        [
                            0.64,
                            -6.839e-08
                        ],
                        [
                            0.63,
                            -7.124e-08
                        ],
                        [
                            0.62,
                            -7.096e-08
                        ],
                        [
                            0.61,
                            -6.909e-08
                        ],
                        [
                            0.6,
                            -6.612e-08
                        ],
                        [
                            0.59,
                            -6.289e-08
                        ],
                        [
                            0.58,
                            -5.962e-08
                        ],
                        [
                            0.57,
                            -5.649e-08
                        ],
                        [
                            0.56,
                            -5.355e-08
                        ],
                        [
                            0.55,
                            -5.077e-08
                        ],
                        [
                            0.54,
                            -4.821e-08
                        ],
                        [
                            0.53,
                            -4.594e-08
                        ],
                        [
                            0.52,
                            -4.385e-08
                        ],
                        [
                            0.51,
                            -4.204e-08
                        ],
                        [
                            0.5,
                            -4.035e-08
                        ],
                        [
                            0.49,
                            -3.882e-08
                        ],
                        [
                            0.48,
                            -3.732e-08
                        ],
                        [
                            0.47,
                            -3.595e-08
                        ],
                        [
                            0.46,
                            -3.466e-08
                        ],
                        [
                            0.45,
                            -3.341e-08
                        ],
                        [
                            0.44,
                            -3.228e-08
                        ],
                        [
                            0.43,
                            -3.119e-08
                        ]
                    ]
                ],
                "n": 1,
                "current_cathodic": 1.672e-06,
                "current_anodic": -7.124e-08,
                "diffusion": 1.451e-09
            }
        }
    ]

    # g2c = Gaus2FrontCharacterization(
    #     _id=mol_id,
    #     calculation_type=calculation_type,
    #     conditions=conditions,
    #     charge=charge,
    #     insert=False,
    # )
    # all_data = g2c.get_all_data()

    print(CV2Front(backend_data=ex_data, run_anodic=False, insert=False).meta_dict)
