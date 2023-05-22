import os, json
from d3tales_api.D3database.d3database import BackDB

test_file = os.path.join("raw_data", "nlp_backend_test.json")
with open(test_file) as fn:
    instance = json.load(fn)

back_id = BackDB(collection_name='nlp', instance=instance).id
print(back_id)
