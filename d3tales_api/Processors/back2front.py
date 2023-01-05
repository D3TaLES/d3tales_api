import copy
import hashlib
from d3tales_api.database_info import db_info
from rdkit.Chem.AllChem import ComputeMolVolume
from d3tales_api.Calculators.calculators import *
from d3tales_api.D3database.d3database import FrontDB
from d3tales_api.D3database.d3database import DBconnector
from ocelot.routines.conformerparser import pmgmol_to_rdmol


class Gaus2FrontCharacterization:
    """
    Update frontend db with backend db data for a particular set of data.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, _id, calculation_type, conditions, charge, data=None, insert=True, all_props=False, rmsd=True):
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
        self.all_calcs = ['opt_groundState', 'freq_groundState', 'solv_energy_gsgs', 'tddft_groundState', 'opt_cation1',
                          'freq_cation1', 'solv_energy_c1c1', 'tddft_cation1', 'opt_cation2', 'freq_cation2',
                          'solv_energy_c2c2',
                          'tddft_cation2', 'opt_anion1', 'freq_anion1', 'solv_energy_a1a1', 'tddft_anion1',
                          'opt_anion2',
                          'freq_anion2', 'solv_energy_a2a2', 'tddft_anion2', 'energy_gsc1', 'energy_gsa1',
                          'energy_c1gs', 'energy_a1gs']

        for calc_type in self.all_calcs:
            if "solv_" in calc_type:
                setattr(self, calc_type, self.find_data(calc_type=calc_type, solvent=self.solvent))
            else:
                setattr(self, calc_type, self.find_data(calc_type=calc_type))

        if insert:
            # generate data dictionary
            self.character_dict = {}
            properties = [
                self.hole_reorganization_energy, self.electron_reorganization_energy,
                self.relaxation_groundState_cation1,
                self.relaxation_cation1_groundState, self.relaxation_groundState_anion1,
                self.relaxation_anion1_groundState,
                self.vertical_ionization_energy, self.vertical_electron_affinity, self.adiabatic_ionization_energy,
                self.adiabatic_ionization_energy_2, self.adiabatic_electron_affinity,
                self.adiabatic_electron_affinity_2,
                self.oxidation_potential, self.reduction_potential
            ]
            if rmsd:
                properties.extend([self.rmsd_groundState_cation1, self.rmsd_cation1_cation2, self.rmsd_groundState_anion1, self.rmsd_anion1_anion2])
            for prop in properties:
                try:
                    self.character_dict.update(prop())
                except ValueError:
                    pass
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
    def from_data(cls, processing_data):
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
        return cls(_id=_id, calculation_type=calculation_type, conditions=conditions, charge=charge, data=data)

    def return_descriptor_dict(self, value, unit="", hashes=None, name="", order=1, condition_addition=None):
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

    def vertical_electron_affinity(self):
        """Descriptor dict for the """
        h_ids, c_data = self.get_data(["opt_groundState", "energy_gsa1"])
        connector = {"opt_energy": "opt_groundState.scf_total_energy.value",
                     "energy": "energy_gsa1.scf_total_energy.value"}
        energy = RelaxationCalc(connector=connector).calculate(c_data)
        return self.return_descriptor_dict(energy, unit='eV', hashes=h_ids, name="vertical_electron_affinity")

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
        all_calcs = ['opt_groundState', 'freq_groundState', 'solv_energy_gsgs', 'tddft_groundState', 'opt_cation1',
                     'freq_cation1', 'solv_energy_c1c1', 'tddft_cation1', 'opt_cation2', 'freq_cation2',
                     'solv_energy_c2c2', 'tddft_cation2', 'opt_anion1', 'freq_anion1', 'solv_energy_a1a1',
                     'tddft_anion1', 'opt_anion2', 'freq_anion2', 'solv_energy_a2a2', 'tddft_anion2', 'energy_gsc1',
                     'energy_gsa1', 'energy_c1gs', 'energy_a1gs']
        h_ids, c_data = self.get_data(all_calcs, disregard_missing=True)
        c_data.update({"h_ids": h_ids})
        return c_data


if __name__ == "__main__":
    back_dbc = DBconnector(db_info.get("backend"))
    back_coll = back_dbc.get_collection("computation")
    _id = "28ef46bbcbb793fc0a8a4656e01aba3e"  # 502dc467c4780db94cc0c324a12c2b6b
    data = back_coll.find_one({"_id": _id})
    mol_id = data["mol_id"]
    calculation_type = data["calculation_type"]
    conditions = data["data"]["conditions"]
    charge = data["data"]["charge"]

    g2c = Gaus2FrontCharacterization(
        _id=mol_id,
        calculation_type=calculation_type,
        conditions=conditions,
        charge=charge,
        insert=False,
    )
    all_data = g2c.get_all_data()
