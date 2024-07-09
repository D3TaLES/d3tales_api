from d3tales_api.D3database.back2front import DOI2Front
import os

DOI = '10.1016/j.chempr.2019.04.021'
TEST_FILE = False

if TEST_FILE:
    test_file = os.path.join("raw_data", "nlp_backend_test.json")
    print(DOI2Front.from_json(test_file, insert=False).extracted_mol_data)
else:
    print(DOI2Front(doi=DOI, insert=True).mol_ids)
