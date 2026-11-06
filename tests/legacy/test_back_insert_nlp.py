import os
import pandas as pd
from d3tales_api.D3database.d3database import BackDB
from d3tales_api.Processors.parser_nlp import ProcessNlp

FROM_TEST_FILE = False
FROM_DATAFRAME= False
FROM_HTML= True

DOI = "10.1016/j.chempr.2019.04.021"
NLP_MODEL = "test_model"

if FROM_TEST_FILE:
    # Get nlp database data from test file
    test_file = os.path.join("raw_data", "nlp_backend_test.json")
    instance = ProcessNlp.from_json(test_file, article_download=True).data_dict
elif FROM_DATAFRAME:
    # Get nlp database data from test file
    nlp_df = pd.read_csv(os.path.join("raw_data", "molecule_properties.csv"))
    instance = ProcessNlp.from_dataframe(nlp_df, NLP_MODEL, doi=DOI).data_dict
elif FROM_HTML:
    # Get nlp database data from test file of HTML
    html_txt = os.path.join("raw_data", "molecule_properties.html")
    instance = ProcessNlp.from_html(html_txt, NLP_MODEL, doi=DOI).data_dict
else:
    # Get (scrape) nlp database data from DOI
    instance = ProcessNlp(DOI, article_download=True).data_dict

# print(instance)
back_id = BackDB(collection_name='nlp', instance=instance, last_updated=True).id
print(back_id)
