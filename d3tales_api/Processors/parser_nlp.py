import os
import sys
from datetime import datetime
from articledownloader.articledownloader import ArticleDownloader

from d3tales_api.Processors.parser_dft import *
from d3tales_api.Processors.parser_uvvis import *
from d3tales_api.Calculators.utils import float_from_str
from d3tales_api.D3database.schema2class import Schema2Class

try:
    from chemdataextractor2 import Document
    from chemdataextractor2.doc.text import *
except ImportError:
    raise ImportError("WARNING. ChemDataExtractor2 not installed! Install ChemDataExtractor if you plan on "
                      "performing NLP parsing.")
try:
    from elsapy.elsdoc import AbsDoc, Paragraph
    from elsapy.elsclient import ElsClient
except ImportError:
    warnings.warn("WARNING. Elsevier's elsapy has not imported correctly! If you plan on performing NLP parsing, "
                  "please resolve this issue.")


class ProcessNlp:
    """
    Class to process NLP data for backend database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, doi=None, instance=None, els_apikey="3a5b9e26201041e25eee70953c65782c", article_download=False,
                 download_dir="temp/"):
        """
        :param doi: str, DOI or scopus id
        :param instance: dict, NLP data. Instance data will override web-scrapped data.
        :param els_apikey: str, API access key for api.elsevier.com
        :param article_download: bool, download article and insert article main text if True
        :param els_apikey: str, directory path into which to download article
        """
        instance = instance or {}
        self.doi = (doi or instance.get("doi", instance.get("_id")) or "").strip("https://doi.org/")
        if not self.doi:
            raise ValueError("ProcessNlp requires a DOI. Either include doi as an argument or include doi as a "
                             "key in the instance data.")
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
        data_dict = self.basic_info
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

    @classmethod
    def from_parsed_json(cls, json_path, **kwargs):
        """
        Generate data class from JSON file

        :param json_path: str, path to JSON file
        :return: data class
        """
        with open(json_path) as fn:
            processing_data = json.load(fn)
        return cls(instance=processing_data, **kwargs)

    @classmethod
    def from_dataframe(cls, nlp_df, nlp_model, doi=None, base_instance=None, date_generated=datetime.now(), **kwargs):
        """
        Translate pandas DataFrame from NLP extraction to JSON data conforming to the D3TaLES backend NLP
        schema: https://github.com/D3TaLES/schema/blob/main/schema_backend/nlp.schema.json. Note that the DataFrame
        should contain the all the columns listed in the expected_columns variable, and no other columns.

        :param nlp_df: pandas DataFrame, extracted NLP data
        :param nlp_model: str, name of NLP model used for data extraction
        :param doi: str, DOI associated with extracted data
        :param base_instance: dict, base instance containing additional info for backend data insertion
        :param date_generated: datetime obj, date and time the data was generated
        :return: dict, JSON formatted data conforming to the D3TaLES backend   NLP schema
        """
        # Check for DOI
        doi_str = (doi or base_instance.get("doi", base_instance.get("_id"))).strip("https://doi.org/")
        if not doi_str:
            raise ValueError(
                "ProcessNlp.from_dataframe requires a DOI. Either include doi as an argument or include doi as a key in the base_instance data.")

        # Check DF columns
        expected_columns = ["molecule", "property", "value", "unit", "line_number", "parent_sentence", "notes"]
        unexpected_columns = [c for c in nlp_df.columns if c not in expected_columns]
        if unexpected_columns:
            raise SyntaxError(
                "The column titles for the NLP DataFrame contain unexpected columns:  " + ", ".join(unexpected_columns))
        missing_columns = [c for c in expected_columns if c not in nlp_df.columns]
        if missing_columns:
            raise SyntaxError(
                "The column titles for the NLP DataFrame do not contain columns: " + ", ".join(missing_columns))

        # Check values in "properties" column
        nlp_df.property = nlp_df.property.str.replace(' ', '_')
        nlp_df.property = nlp_df.property.str.lower()
        expected_properties = ["oxidation_potential", "reduction_potential", "solubility", "stability", "conductivity",
                               "diffusion_coefficient", "charge_transfer_rate"]
        unexpected_properties = [p for p in set(nlp_df.property.tolist()) if p not in expected_properties]
        if unexpected_properties:
            raise ValueError("There were unexpected properties in the NLP DataFrame 'properties' column: ",
                             ", ".join(unexpected_properties))

        # Generate output JSON
        extracted_molecules = []
        date_generated = date_generated.isoformat() if isinstance(date_generated, datetime) else date_generated
        mol_names = set(nlp_df.molecule.tolist())
        for mol in mol_names:
            extracted_properties = {}
            mol_df = nlp_df[nlp_df.molecule == mol]
            for prop in set(mol_df.property.tolist()):
                props = []
                prop_df = mol_df[mol_df.property == prop]
                for i, row in prop_df.iterrows():
                    row.fillna("", inplace=True)
                    prop_data = {
                        "conditions": {
                            "data_source": "nlp",
                            "nlp_model": nlp_model,
                            "date_generated": date_generated,
                            "doi": doi_str,
                        },
                        "value": float_from_str(row.value),
                        "unit": row.unit
                    }
                    prop_data.update({p: row.get(p) for p in ["line_number", "parent_sentence", "notes"] if row.get(p)})
                    props.append(prop_data)
                extracted_properties[prop] = props
            extracted_molecules.append({"molecule_name": mol, "extracted_properties": extracted_properties})
        instance = base_instance or {}
        instance.update({"_id": doi_str, "extracted_molecules": extracted_molecules})
        return cls(instance=instance, **kwargs)

    @classmethod
    def from_html(cls, html_txt, nlp_model, doi=None, base_instance=None, date_generated=datetime.now(), **kwargs):
        nlp_dfs = pd.read_html(html_txt)
        if len(nlp_dfs) > 1:
            raise ValueError(
                "ProcessNlp.from_html found {} tables in the HTML string. There should be only one.".format(
                    len(nlp_dfs)))
        if len(nlp_dfs) < 1:
            raise ValueError("ProcessNlp.from_html no tables in the HTML string. There should be  one.")
        return cls.from_dataframe(nlp_dfs[0], nlp_model, doi=doi, base_instance=base_instance, date_generated=date_generated, **kwargs)

    @classmethod
    def from_prop_list(cls, prop_list, nlp_model, doi=None, base_instance=None, date_generated=datetime.now(), **kwargs):
        if len(prop_list) < 1:
            raise ValueError("ProcessNlp.from_prop_list no propreties in prop_list. There should be at least one.")
        nlp_df = pd.DataFrame(prop_list)
        return cls.from_dataframe(nlp_df, nlp_model, doi=doi, base_instance=base_instance, date_generated=date_generated, **kwargs)


if __name__ == "__main__":
    data = ProcessNlp(sys.argv[1], article_download=False).data_dict

