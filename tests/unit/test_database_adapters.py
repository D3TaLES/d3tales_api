import pytest

from d3tales_api.D3database.d3database import D3Database


class FakeCollection:
    def __init__(self):
        self.updates = []
        self.documents = {}

    def update_one(self, query, update, upsert=False):
        self.updates.append((query, update, upsert))
        return None

    def count_documents(self, query):
        return 0

    def find_one(self, query, projection=None):
        return self.documents.get(query.get("_id"))

    def find(self, query, projection=None):
        return iter([])


@pytest.mark.unit
def test_dot2dict_expands_dotted_keys():
    result = D3Database.dot2dict({"mol_info.smiles": "C", "plain": 1})
    assert result["mol_info"]["smiles"] == "C"
    assert result["plain"] == 1


@pytest.mark.unit
def test_field_checker_returns_empty_for_missing_field():
    db = D3Database.__new__(D3Database)
    db.coll = FakeCollection()
    assert db.field_checker("missing", "mol_info.smiles") == {}


@pytest.mark.unit
def test_path_insert_scalar_and_list():
    db = D3Database.__new__(D3Database)
    db.coll = FakeCollection()
    db.path_insert("M1", "plain", 42)
    db.path_insert("M1", "tags", ["a", "b"])
    assert db.coll.updates[0][1] == {"$set": {"plain": 42}}
    assert db.coll.updates[-1][1]["$addToSet"]["tags"]["$each"] == ["a", "b"]
