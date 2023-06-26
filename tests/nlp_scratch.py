import pandas as pd
from datetime import datetime

def nlp_pd2json(nlp_df, doi, nlp_model, date_generated=datetime.now()):
    """
    Translate pandas DataFrame from NLP extraction to JSON data conforming to the D3TaLES backend NLP
    schema: https://github.com/D3TaLES/schema/blob/main/schema_backend/nlp.schema.json. Note that the DataFrame
    should contain the all the columns listed in the expected_columns variable, and no other columns.

    :param nlp_df: pandas DataFrame, extracted NLP data
    :param doi: str, DOI associated with extracted data
    :param nlp_model: str, name of NLP model used for data extraction
    :param date_generated: datetime obj, date and time the data was generated
    :return: dict, JSON formatted data conforming to the D3TaLES backend   NLP schema
    """
    # Check DF columns
    expected_columns = ["molecule", "property", "value", "unit", "line_number", "parent_sentence", "notes"]
    unexpected_columns = [c for c in nlp_df.columns if c not in expected_columns]
    if unexpected_columns:
        raise SyntaxError("The column titles for the NLP DataFrame do not match the expected column titles. Columns {}"
                          "were unexpected.".format(", ".join(unexpected_columns)))
    if nlp_df.columns != expected_columns:
        raise SyntaxError("The column titles for the NLP DataFrame do not contain columns: ".format(
            ", ".join([c for c in expected_columns if c not in nlp_df.columns])))

    # Check values in "properties" column
    expected_properties = ["oxidation_potential", "reduction_potential", "solubility", "stability", "conductivity",
                           "diffusion_coefficient", "charge_transfer_rate"]
    unexpected_properties = [p for p in set(nlp_df.properties.tolist()) if p not in expected_properties]
    if unexpected_properties:
        raise ValueError("There were unexpected properties in the NLP DataFrame 'properties' column: ",
                         ", ".join(unexpected_properties))

    # Generate output JSON
    extracted_molecules = []
    mol_names = set(nlp_df.molecule.tolist())
    for mol in mol_names:
        extracted_properties = {}
        mol_df = nlp_df[nlp_df.molecule == mol]
        for prop in set(mol_df.property.tolist()):
            props = []
            prop_df = mol_df[mol_df.property == prop]
            for i, row in prop_df.iterrows():
                props.append({
                    "conditions": {
                        "data_source": "nlp",
                        "nlp_model": nlp_model,
                        "date_generated": date_generated,
                        "doi": doi,
                        "line_number": row.line_number,
                        "parent_sentence": row.parent_sentence,
                        "notes": row.notes
                    },
                    "value": row.value,
                    "unit": row.unit
                })
            extracted_properties[prop] = props
        extracted_molecules.append({"molecule_name": mol, "extracted_properties": extracted_properties})

    return {"extracted_molecules": extracted_molecules}



# NLP D3TALES API CODE
from d3tales_api.Processors.d3tales_parser import ProcessNlp

DOI = "10.1039/c7cs00569e"
NLP_MODLE = ""

# Generate initial data, download article
process_obj = ProcessNlp(DOI, article_download=True, download_dir="temp/")
nlp_data = process_obj.data_dict
article_path = process_obj.download_path

# PERFORM NLP
nlp_df = pd.DataFrame()  # NLP OUTPUT!
generation_data = ""  # Some datetime object

# Convert NLP data to backend database entry
backend_data = ProcessNlp.from_dataframe(nlp_df, NLP_MODLE, base_instance=nlp_data, date_generated=generation_data).data_dict
