from d3tales_api.Processors.back2front import DOI2Front
import os

IDS = ['0b5640fb61e165c8117cff02c59e9af5']

if __name__ == "__main__":
    test_file = os.path.join("raw_data", "nlp_backend_test.json")
    print(DOI2Front.from_json(test_file, insert=False).extracted_mol_data)