import os
import sys
import pint
import uuid
import json
import hashlib
from d3tales_api.database_info import db_info
from d3tales_api.Processors.parser_cv import *
from d3tales_api.Processors.parser_dft import *
from d3tales_api.Processors.parser_uvvis import *
from d3tales_api.D3database.db_connector import DBconnector


class ProcessDFT:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, _id=None, submission_info=None, metadata=None, parsing_class=ProcessGausLog):
        submission_info = submission_info or {}
        self.id = _id
        self.parsing_class = parsing_class
        self.submission_info = json.loads(json.dumps(submission_info, ))

        self.DFTData = parsing_class(metadata=metadata)

    @property
    def data_dict(self):
        if self.DFTData.calculation_type == 'wtuning':
            if not hasattr(self.DFTData, 'omega'):
                raise TypeError("{} parsing class is not a WTuning parsing class.".format(self.parsing_class))
            all_data_dict = {
                "_id": self.hash_id,
                "mol_id": self.id,
                "submission_info": self.submission_info,
                "calculation_type": self.DFTData.calculation_type,
                "data": {
                    "conditions": self.conditions,
                    "omega": self.DFTData.omega
                }
            }
            json_data = json.dumps(all_data_dict)
            return json.loads(json_data)
        data_dict = {
            "conditions": self.conditions,
            "charge": self.DFTData.charge,
            "spin_multiplicity": self.DFTData.spin_multiplicity,
            "is_groundState": self.is_groundState,
            "number_of_electrons": sum(self.DFTData.electrons)
        }
        if 'freq' in self.DFTData.calculation_type:
            data_dict.update({
                "gibbs_correction": {
                    "value": self.DFTData.gibbs_correction * 27.2114,  # convert to eV
                    "unit": "eV"
                },
                "frequency_dict": self.DFTData.frequency_dicts
            })
        elif 'opt' in self.DFTData.calculation_type or 'energy' in self.DFTData.calculation_type:
            data_dict.update({
                "scf_total_energy": {
                    "value": self.DFTData.final_energy,
                    "unit": "eV"
                },
                "scf_dipole_moment": {
                    "value": self.DFTData.dipole_moments[-1],
                    "unit": "Debye"
                },
                "homo": {
                    "value": self.DFTData.homo,
                    "unit": "eV"
                },
                "lumo": {
                    "value": self.DFTData.lumo,
                    "unit": "eV"
                },
            })
            if 'opt' in self.DFTData.calculation_type:
                if (int(self.DFTData.spin_multiplicity) % 2) == 0:
                    rss_dict = self.DFTData.get_radical_stability_score(spin_type="mulliken")
                    if rss_dict:
                        data_dict.update(rss_dict)
                data_dict.update({"geometry": self.DFTData.final_structure})
        elif 'tddft' in self.DFTData.calculation_type:
            data_dict.update({
                "excitations": self.DFTData.tddft_excitations,
                "singlet_plotting": self.DFTData.tddft_spectra_data,
                "scf_dipole_moment": {
                    "value": self.DFTData.dipole_moments[-1],
                    "unit": "Debye"
                },
            })
        all_data_dict = {
            "_id": self.hash_id,
            "mol_id": self.id,
            "submission_info": self.submission_info,
            "calculation_type": self.DFTData.calculation_type,
            "runtime": self.DFTData.runtime,
            "data": data_dict
        }
        json_data = json.dumps(all_data_dict)
        return json.loads(json_data)

    @property
    def hash_id(self):
        hash_dict = {
            "_id": self.id,
            "calculation_type": self.DFTData.calculation_type,
            "conditions": self.conditions,
        }
        dhash = hashlib.md5()
        encoded = json.dumps(hash_dict, sort_keys=True).encode()
        dhash.update(encoded)
        return dhash.hexdigest()

    @property
    def conditions(self):
        data_dict = {
            "data_source": 'dft',
            "code_name": self.DFTData.code_name,
            "code_version": self.DFTData.code_version,
            "functional": self.DFTData.functional,
            "basis_set": self.DFTData.basis_set,
        }
        if self.DFTData.tuning_parameter:
            data_dict['tuning_parameter'] = self.DFTData.tuning_parameter
        if getattr(self.DFTData, "solvent", None):
            data_dict['solvent'] = {
                'name': self.DFTData.solvent,
                'model': 'implicit_solvent',
                'dielectric_constant': self.DFTData.dielectric_constant
            }
        return data_dict

    @property
    def mol_info(self):
        base_coll = DBconnector(db_info.get("frontend")).get_collection("base")
        document = base_coll.find_one({"_id": self.id})
        if document:
            return document['mol_info']
        else:
            raise IOError("No molecule with id {} exists in the frontend database. Create an instance in the frontend "
                          "database first.".format(self.id))

    @property
    def is_groundState(self):
        mol_info = self.mol_info
        if 'groundState_charge' in mol_info.keys():
            gs_charge = mol_info['groundState_charge']
            return True if gs_charge == self.DFTData.charge else False
        else:
            raise IOError("The molecule does not have a specified groundState_charge.")


class ProcessCV:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, _id=None, submission_info=None, metadata=None, parsing_class=ParseChiCV):
        self.id = _id
        self.hash_id = str(uuid.uuid4())
        self.data_path = filepath
        self.submission_info = submission_info or {}
        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.experiment_run_id = metadata.get("experiment_run_id", '')
        self.working_electrode = metadata.get("electrode_working", '')
        self.counter_electrode = metadata.get("electrode_counter", '')
        self.reference_electrode = metadata.get("electrode_reference", '')
        self.temperature = metadata.get("temperature", '')
        self.redox_mol_concentration = metadata.get("redox_mol_concentration", '')
        self.working_electrode_surface_area = metadata.get("working_electrode_surface_area", '')
        self.solvents = metadata.get("solvent") if isinstance(metadata.get("solvent"), list) else [
            metadata.get("solvent")] if metadata.get("solvent") else []
        self.electrolytes = metadata.get("electrolyte") if isinstance(metadata.get("electrolyte"), list) else [
            metadata.get("electrolyte")] if metadata.get("electrolyte") else []
        self.ionic_liquids = metadata.get("ionic_liquid") if isinstance(metadata.get("ionic_liquid"), list) else [
            metadata.get("ionic_liquid")] if metadata.get("ionic_liquid") else []

        self.CVData = parsing_class(filepath)

    @property
    def data_dict(self):
        cv_data = self.CVData.cv_data
        cv_data.update(self.CVData.calculate_prop("peaks"))
        cv_data.update(dict(conditions=self.cv_conditions,
                            plot_data=self.CVData.calculate_plotting("plot_data"),
                            reversibility=self.CVData.calculate_prop("reversibility", return_type=list),
                            e_half=self.CVData.calculate_prop("e_half", return_type=list),
                            peak_splittings=self.CVData.calculate_prop("peak_splittings", return_type=list),
                            middle_sweep=self.CVData.calculate_prop("middle_sweep", return_type=list),
                            ))
        all_data_dict = {
            "_id": self.hash_id,
            "mol_id": self.id,
            "submission_info": self.submission_info,
            "data": cv_data
        }
        print("CV conditions: ", self.cv_conditions)
        json_data = json.dumps(all_data_dict)
        return json.loads(json_data)

    @property
    def mol_info(self):
        base_coll = DBconnector(db_info.get("frontend")).get_collection("base")
        document = base_coll.find({"_id": self.id})
        if document:
            return document['mol_info']
        else:
            raise IOError(
                "No molecule with id {} exists in the frontend database. Create an instance in the frontend database "
                "first.".format(id))

    @property
    def cv_conditions(self):
        conditions_data = self.CVData.conditions_data
        conditions_data.update({
            "data_source": 'cv',
            "working_electrode": self.working_electrode,
            "counter_electrode": self.counter_electrode,
            "reference_electrode": self.reference_electrode,
            "solvent": self.reagent_format(self.solvents),
            "electrolyte": self.reagent_format(self.electrolytes),
            "ionic_liquid": self.reagent_format(self.ionic_liquids),
            "instrument": self.instrument,
            "working_electrode_surface_area": self.data_format(self.working_electrode_surface_area),
            "redox_mol_concentration": self.data_format(self.redox_mol_concentration),
            "temperature": self.data_format(self.temperature),
            "experiment_run_id": self.experiment_run_id
        })
        return {k: v for k, v in conditions_data.items() if v}

    @staticmethod
    def reagent_format(reagents, purity=None):
        def format_i(r):
            if isinstance(r, dict):
                r_dict = r
            else:
                r_dict = {"name": str(r)}
                if purity:
                    r_dict["purity"] = purity
            return r_dict
        return [format_i(i) for i in reagents] if isinstance(reagents, list) else format_i(reagents)

    @staticmethod
    def data_format(data):
        if isinstance(data, dict):
            return {"value": float(data.get("value")), "unit": data.get("unit")}
        elif not isinstance(data, (str, float, int, dict)):
            value, unit = getattr(data, "value"), getattr(data, "unit")
            return {"value": value, "unit": unit}
        elif isinstance(data, float) or isinstance(data, float) or str(data).replace('.', '', 1).replace('-', '', 1).isdigit():
            return {"value": float(data)}
        else:
            ureg = pint.UnitRegistry()
            pint_m = ureg(str(data))
            return {"value": float(getattr(pint_m, 'magnitude')), "unit": str(getattr(pint_m, 'units'))}


class ProcessUvVis:
    """
    Class to process UV-Vis data files.
    Copyright 2021, University of Kentucky
    Args:
        filepath (str) : filepath to data file
        mol_id (str) : identifier for the molecule this data belongs to
        metadata (dict) : dictionary containing any metadata for this molecule, e.g., {"solvent": "acetonitrile"}
        parsing_class (class) : a UV-Vis parsing class with which to parse the data
    """

    def __init__(self, filepath, mol_id, metadata=None, parsing_class=ParseExcel):
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())

        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.solvent = metadata.get("solvent", '')

        self.UvVisData = parsing_class(filepath)

    @property
    def no_sql_data(self):
        """
        Returns UV-Vis information in a dictionary that matches the No-SQL schema
        """
        all_data_dict = {
            "date_recorded": self.UvVisData.date_recorded,
            "solvent": self.solvent,
            "instrument": self.instrument,
            "integration_time": self.UvVisData.integration_time,
            "absorbance_data": self.UvVisData.absorbance_data,
        }
        json_data = json.dumps(all_data_dict, default=str)
        return json.loads(json_data)

    @property
    def sql_absorbance_data(self):
        """
        Returns UV-Vis information in a dictionary that matches the SQL AbsorbanceData Table schema
        """
        data = self.UvVisData.absorbance_data
        return [{"uvvis_id": self.uuid,
                 "mol_id": self.mol_id,
                 "wavelength": wavelength,
                 "absorbance": absorbance}
                for wavelength, absorbance in zip(data["wavelength"], data["absorbance"])]

    @property
    def sql_data(self):
        """
        Returns UV-Vis information in a dictionary that matches the SQL UVVisData Table schema
        """
        data_dict = {
            "uvvis_id": self.uuid,
            "mol_id": self.mol_id,
            "date_recorded": self.UvVisData.date_recorded,
            "solvent": self.solvent,
            "instrument": self.instrument,
            "integration_time": self.UvVisData.integration_time,
            "absorbance_data": self.UvVisData.absorbance_data,
        }
        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)


if __name__ == "__main__":
    # data = ProcessDFT(sys.argv[1], parsing_class=ProcessGausLog).data_dict
    data = ProcessCV(sys.argv[1], parsing_class=ParseChiCV).data_dict