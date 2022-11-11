import os
import json
import copy
import hashlib
import numpy as np
from rdkit.Chem import rdMolAlign
from pymatgen.core.sites import Site
from pymatgen.core.structure import Molecule
from d3tales_api.database_info import db_info
from rdkit.Chem.AllChem import ComputeMolVolume
from ocelot.routines.conformerparser import pmgmol_to_rdmol
from d3tales_api.D3database.db_connector import DBconnector
from d3tales_api.D3database.d3database import FrontDB, ParamsDB


class Gaus2FrontCharacterization:
    """
    Update frontend db with backend db data for a particular set of data.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, _id, calculation_type, conditions, charge, data=None, insert=True, all_props=False):
        # connect to databases
        self.front_dbc = DBconnector(db_info("frontend"))
        self.front_coll = self.front_dbc.get_collection("base")
        self.back_dbc = DBconnector(db_info("backend"))
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
                          'freq_cation1', 'solv_energy_c1c1', 'tddft_cation1', 'opt_cation2', 'freq_cation2', 'solv_energy_c2c2',
                          'tddft_cation2', 'opt_anion1', 'freq_anion1', 'solv_energy_a1a1', 'tddft_anion1', 'opt_anion2',
                          'freq_anion2', 'solv_energy_a2a2', 'tddft_anion2', 'energy_gsc1', 'energy_gsa1', 'energy_c1gs', 'energy_a1gs']

        for calc_type in self.all_calcs:
            if "solv_" in calc_type:
                setattr(self, calc_type, self.find_data(calc_type=calc_type, solvent=self.solvent))
            else:
                setattr(self, calc_type, self.find_data(calc_type=calc_type))

        if insert:
            # generate data dictionary
            self.character_dict = {}
            properties = [
                self.hole_reorganization_energy, self.electron_reorganization_energy, self.relaxation_groundState_cation1,
                self.relaxation_cation1_groundState, self.relaxation_groundState_anion1, self.relaxation_anion1_groundState,
                self.vertical_ionization_energy, self.vertical_electron_affinity, self.adiabatic_ionization_energy,
                self.adiabatic_ionization_energy_2, self.adiabatic_electron_affinity, self.adiabatic_electron_affinity_2,
                self.rmsd_groundState_cation1, self.rmsd_cation1_cation2, self.rmsd_groundState_anion1,
                self.rmsd_anion1_anion2, self.oxidation_potential, self.reduction_potential
            ]
            for prop in properties:
                try:
                    self.character_dict.update(prop())
                except ValueError:
                    pass
            FrontDB(schema_layer="species_characterization", instance=self.species_descriptors, _id=self.id)
            FrontDB(schema_layer="mol_characterization", instance=self.character_dict, _id=self.id)

    @property
    def species_descriptors(self):
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
                        self.data[opt_descriptor].update({"source_hash_ids": [self._hash], "conditions": self.conditions})
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
        _id = processing_data.get("mol_id")
        calculation_type = processing_data.get("calculation_type")
        data = processing_data.get("data", {})
        conditions = data.get("conditions")
        charge = data.get("charge")
        return cls(_id=_id, calculation_type=calculation_type, conditions=conditions, charge=charge, data=data)

    def solvation_energy(self, species, idx=0):
        species_abbrev_dict = {
            "anion2": "a2",
            "anion1": "a1",
            "groundState": "gs",
            "cation1": "c1",
            "cation2": "c2",
        }
        species_abbrev = species_abbrev_dict.get(species)
        gas_phase = "opt_{}".format(species)
        solv = "solv_energy_{}{}".format(species_abbrev, species_abbrev)
        h_ids, unit = self.get_hash_ids([gas_phase, solv])
        solv_eng = getattr(self, solv)[idx]["scf_total_energy"]["value"] - getattr(self,gas_phase)[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(solv_eng, unit=unit, hashes=h_ids, name="solvation_energy")

    def return_descriptor_dict(self, value, unit="", hashes=None, name="", order=1, condition_addition=None):
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

    def hole_reorganization_energy(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "energy_gsc1", "opt_cation1", "energy_c1gs"])
        lambda1 = self.energy_gsc1[idx]["scf_total_energy"]["value"] - self.opt_cation1[idx]["scf_total_energy"]["value"]
        lambda2 = self.energy_c1gs[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(lambda1 + lambda2, unit=unit, hashes=h_ids, name="hole_reorganization_energy")

    def electron_reorganization_energy(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "energy_gsa1", "opt_anion1", "energy_a1gs"])
        lambda1 = self.energy_gsa1[idx]["scf_total_energy"]["value"] - self.opt_anion1[idx]["scf_total_energy"]["value"]
        lambda2 = self.energy_a1gs[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(lambda1 + lambda2, unit=unit, hashes=h_ids, name="hole_reorganization_energy")

    def relaxation_groundState_cation1(self, idx=0):
        h_ids, unit = self.get_hash_ids(["energy_gsc1", "opt_cation1"])
        lambda1 = self.energy_gsc1[idx]["scf_total_energy"]["value"] - self.opt_cation1[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(lambda1, unit=unit, hashes=h_ids, name="relaxation_groundState_cation1")

    def relaxation_cation1_groundState(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "energy_c1gs"])
        lambda2 = self.energy_c1gs[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(lambda2, unit=unit, hashes=h_ids, name="relaxation_cation1_groundState")

    def relaxation_groundState_anion1(self, idx=0):
        h_ids, unit = self.get_hash_ids(["energy_gsa1", "opt_anion1"])
        lambda1 = self.energy_gsa1[idx]["scf_total_energy"]["value"] - self.opt_anion1[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(lambda1, unit=unit, hashes=h_ids, name="relaxation_groundState_anion1")

    def relaxation_anion1_groundState(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "energy_a1gs"])
        lambda2 = self.energy_a1gs[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(lambda2, unit=unit, hashes=h_ids, name="relaxation_anion1_groundState")

    def vertical_ionization_energy(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "energy_gsc1"])
        vip = self.energy_gsc1[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(vip, unit=unit, hashes=h_ids, name="vertical_ionization_energy")

    def vertical_electron_affinity(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "energy_gsa1"])
        vea = self.energy_gsa1[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(vea, unit=unit, hashes=h_ids, name="vertical_electron_affinity")

    def adiabatic_ionization_energy(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "opt_cation1"])
        aip = self.opt_cation1[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(aip, unit=unit, hashes=h_ids, name="adiabatic_ionization_energy")

    def adiabatic_ionization_energy_2(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_cation1", "opt_cation2"])
        aip = self.opt_cation2[idx]["scf_total_energy"]["value"] - self.opt_cation1[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(aip, unit=unit, hashes=h_ids, name="adiabatic_ionization_energy_2")

    def adiabatic_electron_affinity(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "opt_anion1"])
        aea = self.opt_anion1[idx]["scf_total_energy"]["value"] - self.opt_groundState[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(aea, unit=unit, hashes=h_ids, name="adiabatic_electron_affinity")

    def adiabatic_electron_affinity_2(self, idx=0):
        h_ids, unit = self.get_hash_ids(["opt_anion1", "opt_anion2"])
        aea = self.opt_anion2[idx]["scf_total_energy"]["value"] - self.opt_anion1[idx]["scf_total_energy"]["value"]
        return self.return_descriptor_dict(aea, unit=unit, hashes=h_ids, name="adiabatic_electron_affinity_2")

    def rmsd_groundState_cation1(self, idx=0):
        return self.rmsd("opt_groundState", "opt_cation1", name="rmsd_groundState_cation1", idx=idx)

    def rmsd_cation1_cation2(self, idx=0):
        return self.rmsd("opt_cation1", "opt_cation2", name="rmsd_groundState_cation1", idx=idx)

    def rmsd_groundState_anion1(self, idx=0):
        return self.rmsd("opt_groundState", "opt_anion1", name="rmsd_groundState_anion1", idx=idx)

    def rmsd_anion1_anion2(self, idx=0):
        return self.rmsd("opt_anion1", "opt_anion2", name="rmsd_groundState_anion1", idx=idx)

    def oxidation_potential(self, idx=0, electrode="standard_hydrogen_electrode", num_electrons=1):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "freq_groundState", "solv_energy_gsgs", "opt_cation1", "freq_cation1", "solv_energy_c1c1"])
        delta_g_solv_ox = self.delta_g_solv(init_eng=self.opt_groundState[idx]["scf_total_energy"]["value"],
                                            init_corr=self.freq_groundState[idx]["gibbs_correction"]["value"],
                                            init_eng_solv=self.solv_energy_gsgs[idx]["scf_total_energy"]["value"],
                                            fin_eng=self.opt_cation1[idx]["scf_total_energy"]["value"],
                                            fin_corr=self.freq_cation1[idx]["gibbs_correction"]["value"],
                                            fin_eng_solv=self.solv_energy_c1c1[idx]["scf_total_energy"]["value"])
        return self.return_descriptor_dict(
            -(-delta_g_solv_ox / num_electrons + self.get_electrode_potential(electrode)), unit=unit, hashes=h_ids,
            name="oxidation_potential", condition_addition={"reference_electrode": electrode})

    def reduction_potential(self, idx=0, electrode="standard_hydrogen_electrode", num_electrons=1):
        h_ids, unit = self.get_hash_ids(["opt_groundState", "freq_groundState", "solv_energy_gsgs", "opt_anion1", "freq_anion1", "solv_energy_a1a1"])
        delta_g_solv_red = self.delta_g_solv(init_eng=self.opt_groundState[idx]["scf_total_energy"]["value"],
                                            init_corr=self.freq_groundState[idx]["gibbs_correction"]["value"],
                                            init_eng_solv=self.solv_energy_gsgs[idx]["scf_total_energy"]["value"],
                                            fin_eng=self.opt_anion1[idx]["scf_total_energy"]["value"],
                                            fin_corr=self.freq_anion1[idx]["gibbs_correction"]["value"],
                                            fin_eng_solv=self.solv_energy_a1a1[idx]["scf_total_energy"]["value"])
        return self.return_descriptor_dict(
            -delta_g_solv_red / num_electrons - self.get_electrode_potential(electrode), unit=unit, hashes=h_ids,
            name="reduction_potential", condition_addition={"reference_electrode": electrode})

    def find_data(self, calc_type, solvent=None):
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

    def rmsd(self, data1, data2, name="", idx=0):
        h_ids, unit = self.get_hash_ids([data1, data2])
        geom1 = pmgmol_to_rdmol(Molecule.from_sites([Site.from_dict(sd) for sd in getattr(self, data1)[idx]["geometry"]]))[0]
        geom2 = pmgmol_to_rdmol(Molecule.from_sites([Site.from_dict(sd) for sd in getattr(self, data2)[idx]["geometry"]]))[0]
        try:
            rmsd = rdMolAlign.GetBestRMS(geom1, geom2)
        except:
            raise ValueError("Error finding RMSE")
        return {name: [{
            "source_hash_ids": h_ids,
            "conditions": copy.deepcopy(self.conditions),
            "order": 1, "value": round(rmsd, 6), "unit": "A"}]}

    def get_hash_ids(self, calculations, check_unit=True):
        if self.calculation_type not in calculations and not self.all_props:
            # print("Calculation does not use this calculation type")
            raise ValueError("Calculation does not use this calculation type")
        calcs_data = [getattr(self, c) for c in calculations]
        jobs, h_ids = [c[0] for c in calcs_data], [c[1] for c in calcs_data]
        if None in jobs:
            # print("At least one calculation needed for this property is missing.")
            raise ValueError("At least one calculation needed for this property is missing.")
        if check_unit:
            unit = self.check_units(jobs)
            if not unit:
                # print("Error. Energies must have the same units to be compared.")
                raise TypeError("Error. Energies must have the same units to be compared.")
        else:
            unit = ""
        return h_ids, unit

    @staticmethod
    def hash_id(_id, calc_type, gaus_conditions):
        hash_dict = {
            "_id": _id,
            "calculation_type": calc_type,
            "conditions": gaus_conditions,
        }
        dhash = hashlib.md5()
        encoded = json.dumps(hash_dict, sort_keys=True).encode()
        dhash.update(encoded)
        return dhash.hexdigest()

    @staticmethod
    def get_electrode_potential(electrode):
        params_db = ParamsDB(collection_name="electrode", schema_directory="materials")
        electrode_data = params_db.coll.find_one({"_id": electrode})
        abs_potential = electrode_data.get("absolute_potential")
        if abs_potential:
            return abs_potential.get("value")

    @staticmethod
    def check_units(calculations, prop="scf_total_energy"):
        """
        Check if all units for a specific property in a list of jobs are the same
        :param calculations: list of jobs
        :param prop: the property to check for
        :return: bool
        """
        units = []
        for calculation in calculations:
            try:
                units.append(calculation[prop]["unit"])
            except KeyError:
                pass

        if all(x == units[0] for x in units):
            return units[0]

    @staticmethod
    def delta_g_solv(init_eng, init_corr, init_eng_solv, fin_eng, fin_corr, fin_eng_solv):
        """
        Calculate the change in free energy moving from the solvated reduced state to solvated the oxidized state
        """
        g_gas_init = init_eng + init_corr
        g_gas_fin = fin_eng + fin_corr
        delta_g_init_solv = init_eng_solv - init_eng  # entropy correction cancels out because the species are the same
        delta_g_fin_solv = fin_eng_solv - fin_eng  # entropy correction cancels out because the species are the same
        return g_gas_fin - g_gas_init + delta_g_fin_solv - delta_g_init_solv


if __name__ == "__main__":
    back_dbc = DBconnector(db_info("backend"))
    back_coll = back_dbc.get_collection("computation")
    _id = "502dc467c4780db94cc0c324a12c2b6b"
    data = back_coll.find_one({"_id": _id})
    mol_id = data["mol_id"]
    calculation_type = data["calculation_type"]
    conditions = data["data"]["conditions"]
    charge = data["data"]["charge"]

    Gaus2FrontCharacterization(
        _id=mol_id,
        calculation_type=calculation_type,
        conditions=conditions,
        charge=charge,
    )
