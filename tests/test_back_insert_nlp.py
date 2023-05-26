import os, json
from d3tales_api.D3database.d3database import BackDB
from d3tales_api.Processors.d3tales_parser import ProcessNlp

FROM_TEST_FILE = True

if FROM_TEST_FILE:
    # Get nlp database data from test file
    test_file = os.path.join("raw_data", "nlp_backend_test.json")
    with open(test_file) as fn:
        file_data = json.load(fn)
    instance = ProcessNlp(instance=file_data, article_download=True).data_dict
else:
    # Get (scrape) nlp database data from DOI
    doi = "10.1039/c7cs00569e"
    instance = ProcessNlp(doi, article_download=True).data_dict

back_id = BackDB(collection_name='nlp', instance=instance, last_updated=True).id
print(back_id)
