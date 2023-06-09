import os
import sys
import uuid
import hashlib
import warnings
from datetime import datetime
from elsapy.elsdoc import AbsDoc
from elsapy.elsclient import ElsClient
from articledownloader.articledownloader import ArticleDownloader

from d3tales_api.database_info import db_info
from d3tales_api.Processors.parser_cv import *
from d3tales_api.Processors.parser_dft import *
from d3tales_api.Processors.parser_uvvis import *
from d3tales_api.D3database.d3database import DBconnector
from d3tales_api.D3database.schema2class import Schema2Class

try:
    from chemdataextractor2 import Document
    from chemdataextractor2.doc.text import *
except ImportError:
    print("WARNING. ChemDataExtractor2 not installed! Install ChemDataExtractor if you plan on performing NLP parsing.")


class ProcessDFT:
    """
    Class to process molecular DFT output files.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, _id: str = None, filepath: str = None, submission_info: dict = None, metadata: dict = None,
                 parsing_class=ProcessGausLog):
        """
        :param _id: str, molecule ID
        :param submission_info: dict, submission information
        :param metadata: dict, metadata (`mol_file` should be a key to the filepath of the file to be processed)
        :param parsing_class: class, class to use for file parsing (EX: ProcessGausLog, ProcessCCLIB, or ProcessPsi4Log)
        """
        submission_info = submission_info or {}
        self.data_path = filepath or metadata.get("mol_file")
        self.id = _id
        self.parsing_class = parsing_class
        self.submission_info = json.loads(json.dumps(submission_info, ))

        self.DFTData = parsing_class(metadata=metadata)

    @property
    def data_dict(self):
        """
        Dictionary of processed data (in accordance with D3TalES backend schema)
        """
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
            "number_of_electrons": sum(self.DFTData.electrons)
        }
        try:
            data_dict.update({"is_groundState": self.is_groundState})
        except ConnectionError:
            print("Warning. Could not connect to the database, so no 'is_groundState' property was specified. DB_INFO_FILE may not be defined.")
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
                "homo_1": {
                    "value": self.DFTData.homo_1,
                    "unit": "eV"
                },
                "lumo_1": {
                    "value": self.DFTData.lumo_1,
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
        """
        Hash ID
        """
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
        """
        Dictionary of conditions (in accordance with D3TaLES backend schema)
        """
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
        """
        Dictionary containing basic molecule information using the ID from the D3TalES database
        """
        base_coll = DBconnector(db_info.get("frontend")).get_collection("base")
        document = base_coll.find_one({"_id": self.id})
        if document:
            return document['mol_info']
        else:
            raise IOError("No molecule with id {} exists in the frontend database. Create an instance in the frontend "
                          "database first.".format(self.id))

    @property
    def is_groundState(self):
        """
        True if current species is the ground state, else False
        """
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

    def __init__(self, filepath, _id: str = None, submission_info: dict = None, metadata: dict = None,
                 parsing_class=ParseChiCV, **kwargs):
        """
        :param filepath: str, filepath to
        :param _id: str, molecule ID
        :param submission_info: dict, submission info
        :param metadata: dict, metadata
        :param parsing_class: class, class to use for file parsing (EX: ParseChiCV)
        """
        self.id = _id
        self.hash_id = str(uuid.uuid4())
        self.data_path = filepath or metadata.get("mol_file")
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

        self.CVData = parsing_class(filepath, **kwargs)

    @property
    def data_dict(self):
        """
        Dictionary of processed data (in accordance with D3TalES backend schema)
        """
        cv_data = self.CVData.cv_data
        cv_data.update(self.CVData.calculate_prop("peaks"))
        cv_data.update(dict(conditions=self.cv_conditions,
                            plot_data=self.CVData.calculate_plotting("plot_data").get("abs_plot"),
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
        json_data = json.dumps(all_data_dict)
        return json.loads(json_data)

    @property
    def mol_info(self):
        """
        Dictionary containing basic molecule information using the ID from the D3TalES database
        """
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
        """
        Dictionary of conditions (in accordance with D3TaLES backend schema)
        """
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
        """
        Convert reagent data to D3TaLES schema format

        :param reagents: reagent data or list of reagent data
        :param purity: default purity value
        :return: formatted reagent data
        """
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
        """
        Format measurement data

        :param data: measurement data
        :return: formatted measurement data, EX: {"value": 0.5, "unit": "V"}
        """
        if isinstance(data, dict):
            return {"value": float(data.get("value")), "unit": data.get("unit")}
        elif not isinstance(data, (str, float, int, dict)):
            value, unit = getattr(data, "value"), getattr(data, "unit")
            return {"value": value, "unit": unit}
        elif isinstance(data, float) or isinstance(data, float) or str(data).replace('.', '', 1).replace('-', '',
                                                                                                         1).isdigit():
            return {"value": float(data)}
        else:
            ureg = pint.UnitRegistry()
            pint_m = ureg(str(data))
            return {"value": float(getattr(pint_m, 'magnitude')), "unit": str(getattr(pint_m, 'units'))}


class ProcessUvVis:
    """
    Class to process UV-Vis data files.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, mol_id, metadata=None, parsing_class=ParseExcel):
        """
        :param filepath: str, filepath to data file
        :param mol_id: str, molecule ID
        :param metadata: dict, dictionary containing any metadata for this molecule, e.g., {"solvent": "acetonitrile"}
        :param parsing_class: class, a UV-Vis parsing class with which to parse the data
        """
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())

        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.solvent = metadata.get("solvent", '')

        self.UvVisData = parsing_class(filepath)

    @property
    def no_sql_data(self):
        """
        UV-Vis information in a dictionary that matches the No-SQL schema
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
        UV-Vis information in a dictionary that matches the SQL AbsorbanceData Table schema
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
        UV-Vis information in a dictionary that matches the SQL UVVisData Table schema
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

class ProcessNlp:
    """
    Class to process NLP data for backend database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, doi=None, instance=None, els_apikey="3a5b9e26201041e25eee70953c65782c", article_download=False, download_dir="temp/"):
        """
        :param doi: str, DOI or scopus id
        :param instance: dict, NLP data. Instance data will override web-scrapped data.
        :param els_apikey: str, API access key for api.elsevier.com
        :param article_download: bool, download article and insert article main text if True
        :param els_apikey: str, directory path into which to download article
        """
        self.doi = (doi or instance.get("doi", instance.get("_id"))).strip("https://doi.org/")
        if not self.doi:
            raise ValueError("ProcessNlp requires a DOI. Either include doi as an argument or include doi as a key in the instance data.")
        self.els_apikey = els_apikey
        self.instance = instance or {}
        self.publisher, self.main_text = "", ""
        self.article_download = article_download
        self.download_path = "{}/{}.html".format(download_dir, self.doi.replace("/", "_"))
        self.s2c = Schema2Class(database="backend", schema_name="nlp")

        # Get basic article info
        self.basic_info = self.get_basic_info()
        for key in self.basic_info:
            setattr(self, key, self.basic_info[key])

        # Download full article
        if self.article_download:
            os.makedirs(download_dir, exist_ok=True)
            self.download_article()
            self.main_text = self.get_main_text()

    @property
    def data_dict(self):
        data_dict  = self.basic_info
        data_dict.update(self.instance)
        data_dict.update({"_id": self.doi})
        if self.article_download and self.main_text:
            data_dict.update({"main_text": self.main_text, "pdf_location": os.path.abspath(self.download_path)})
        json_data = json.dumps(data_dict)
        return json.loads(json_data)

    def get_main_text(self):
        try:
            doc = Document.from_file(self.download_path)
            elements = doc.elements

            paragraphs, refs, i = [], False, 0
            while i < len(elements) and not refs:
                elem = elements[i]
                if elem.text in ["Abstract"]:
                    i += 1
                    while i < len(elements):
                        elem = elements[i]
                        try:
                            if elem.text in ["References"]:
                                refs = True
                                break
                        except:
                            pass

                        if isinstance(elem, Paragraph):
                            paragraphs.append(elem.text)
                        i += 1
                    break
                i += 1

            data_obj = self.s2c.Data()
            data_obj.text = paragraphs
            return data_obj.as_dict()

        except:
            print('Failed to extract main text from {}'.format(self.doi))

    def get_basic_info(self):
        """
        Dictionary of processed data (in accordance with D3TalES backend schema)
        """
        nlp_obj = self.s2c.Nlp()
        author_obj = self.s2c.Author()

        client = ElsClient(self.els_apikey)
        ab = AbsDoc(uri="https://api.elsevier.com/content/abstract/doi/{}".format(self.doi))
        ab.read(client)

        if ab.read(client):
            print("Gathering data from the web...")
            # put the data into nlp raw data schema
            nlp_obj.publisher = ab.data["coredata"].get('dc:publisher', "")
            nlp_obj.journal = ab.data["coredata"].get('prism:publicationName', "")
            nlp_obj.publish_date = ab.data["coredata"].get('prism:coverDate', "")
            nlp_obj.authors = [author_obj.from_json(json.dumps(self.make_author_entry(author=author)))
                               for author in ab.data["authors"]["author"]]
            nlp_obj._id = ab.data["coredata"].get('prism:doi', "")
            nlp_obj.date_accessed = datetime.now().isoformat()
            nlp_obj.abstract = ab.data["coredata"].get('dc:description', "")
            nlp_obj.title = ab.title
        else:
            raise ConnectionError("Not able to retrieve abstract for {}.".format(self.doi))

        return nlp_obj.as_dict()

    def make_author_entry(self, author):
        return {
            "auid": author.get("@auid", ""),
            "indexed_name": author.get("ce:indexed-name", ""),
            "surname": author.get('ce:surname', ""),
            "given_name": author.get("ce:given-name", ""),
            "affiliations": self.get_affiliations(author.get("affiliation", ""))
        }

    def download_article(self):
        ad = ArticleDownloader(els_api_key=self.els_apikey)

        if self.publisher == "American Chemical Society":
            mode = "acs"
        elif self.publisher == "Royal Society of Chemistry":
            mode = "rsc"
        elif self.publisher == 'Elsevier B.V.':
            mode = "elsevier"
        elif self.publisher == 'Nature Research':
            mode = 'nature'
        elif self.publisher == 'IOP Publishing Ltd':
            mode = "ecs"
        elif self.publisher == 'Wiley-VCH Verlag':
            mode = "wiley"
        elif self.publisher == 'Korean Electrochemical Society':
            mode = "ecs"
        elif self.publisher == 'Institute of Physics Publishing':
            mode = "ecs"
        else:
            mode = "elsevier"

        with open(self.download_path, "wb") as f:
            ad.get_html_from_doi(doi=self.doi, writefile=f, mode=mode)

    @staticmethod
    def get_affiliations(affiliation):
        if isinstance(affiliation, list):
            affiliations = [x.get("@id", "") for x in affiliation]
        elif isinstance(affiliation, dict):
            affiliations = [affiliation.get("@id", "")]
        else:
            affiliations = []
        return affiliations



if __name__ == "__main__":
    # data = ProcessDFT(sys.argv[1], parsing_class=ProcessGausLog).data_dict
    # data = ProcessCV(sys.argv[1], parsing_class=ParseChiCV).data_dict
    data = ProcessNlp(sys.argv[1], article_download=False).data_dict

    print(data)
