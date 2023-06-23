import os
import pandas as pd
from d3tales_api.D3database.d3database import BackDB
from d3tales_api.Processors.d3tales_parser import ProcessNlp

FROM_TEST_FILE = True
FROM_DATAFRAME= False

if FROM_TEST_FILE:
    # Get nlp database data from test file
    test_file = os.path.join("raw_data", "nlp_backend_test.json")
    instance = ProcessNlp.from_json(test_file, article_download=True).data_dict
elif FROM_DATAFRAME:
    # Get nlp database data from test file
    df = pd.from_csv(os.path.join("raw_data", ""))
    # instance = ProcessNlp.from_dataframe(df, article_download=True).data_dict
else:
    # Get (scrape) nlp database data from DOI
    doi = "10.1039/c7cs00569e"
    instance = ProcessNlp(doi, article_download=True).data_dict

back_id = BackDB(collection_name='nlp', instance=instance, last_updated=True).id
print(back_id)
