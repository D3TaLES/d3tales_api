import pytest

from d3tales_api.D3database.d3database import D3Database, ParamsDB, QueryDB


class FakeCollection:
    def __init__(self, records=None):
        self.records = records or []
        self.updates = []

    def update_one(self, query, update, upsert=False):
        self.updates.append((query, update, upsert))

    def count_documents(self, query):
        count = 0
        for record in self.records:
            if record.get("_id") != query.get("_id"):
                continue
            exists_check = query.get("tags", {}).get("$exists")
            if exists_check is None or ("tags" in record) == exists_check:
                count += 1
        return count

    def find_one(self, query, projection=None):
        for record in self.records:
            if all(record.get(key) == value for key, value in query.items()):
                return record
        return None

    def find(self, query, projection=None):
        return FakeCursor(self.records)


class FakeCursor(list):
    def limit(self, count):
        return FakeCursor(self[:count])


@pytest.mark.unit
def test_dot2dict_converts_dotted_paths_and_preserves_nested_values():
    assert D3Database.dot2dict({
        "experiment_data.experiment_ids": ["one"],
        "mol_info": {"smiles": "C"},
    }) == {
        "experiment_data": {"experiment_ids": ["one"]},
        "mol_info": {"smiles": "C"},
    }
    assert D3Database.dot2dict(None) == {}
    assert D3Database.dot2dict(["unchanged"]) == ["unchanged"]


@pytest.mark.unit
def test_insert_nested_fields_and_append_unique_lists():
    db = D3Database.__new__(D3Database)
    db.coll = FakeCollection()
    db.instance = {"details": {"name": "sample"}, "tags": ["a", "b"], "count": 2}
    db.public = True
    db.database = "test"
    db.verbose = 0

    db.insert("item-1", nested=True, override_lists=False)

    assert db.coll.updates == [
        ({"_id": "item-1"}, {"$set": {"public": True}}, True),
        ({"_id": "item-1"}, {"$set": {"details.name": "sample"}}, True),
        ({"_id": "item-1"}, {"$set": {"tags": []}}, True),
        ({"_id": "item-1"}, {"$addToSet": {"tags": {"$each": ["a", "b"]}}}, True),
        ({"_id": "item-1"}, {"$set": {"count": 2}}, True),
    ]


@pytest.mark.unit
def test_path_insert_and_array_checker_cover_existing_and_new_arrays():
    db = D3Database.__new__(D3Database)
    db.coll = FakeCollection(records=[{"_id": "existing", "tags": []}])

    db.path_insert("new", "tags", ["x"])
    db.path_insert("existing", "tags", ["x"])
    db.path_insert("existing", "label", "value")

    assert db.coll.updates == [
        ({"_id": "new"}, {"$set": {"tags": []}}, True),
        ({"_id": "new"}, {"$addToSet": {"tags": {"$each": ["x"]}}}, True),
        ({"_id": "existing"}, {"$addToSet": {"tags": {"$each": ["x"]}}}, True),
        ({"_id": "existing"}, {"$set": {"label": "value"}}, True),
    ]


@pytest.mark.unit
def test_make_query_json_projection_and_single_record():
    db = D3Database.__new__(D3Database)
    db.coll = FakeCollection([{"_id": "a", "value": 1}, {"_id": "b", "value": 2}])

    assert db.make_query(output="json", projection={"value": 1}) == [
        {"_id": "a", "value": 1},
        {"_id": "b", "value": 2},
    ]
    assert db.make_query(output="records", multi=False) == [{"_id": "a", "value": 1}]


@pytest.mark.unit
def test_params_check_id_creates_suffix_and_update_requires_existing():
    params = ParamsDB.__new__(ParamsDB)
    params.coll = FakeCollection([{"_id": "common"}, {"_id": "common_1"}])
    params.update = False

    with pytest.warns(UserWarning, match="new entry"):
        assert params.check_id("common") == "common_2"
    assert params.check_id("available") == "available"

    params.update = True
    assert params.check_id("common") == "common"
    with pytest.raises(IOError, match="no parameter entry"):
        params.check_id("missing")


@pytest.mark.unit
def test_querydb_returns_json_and_passthrough_cursor(monkeypatch):
    collection = FakeCollection([{"_id": "a", "value": 1}])
    monkeypatch.setattr("d3tales_api.D3database.d3database.DBconnector.get_collection",
                        lambda self, name: collection)

    query = QueryDB(collection_name="records", output="json", db_file={})
    assert query.make_query() == '[{"_id": "a", "value": 1}]'

    query.output = "cursor"
    assert list(query.make_query()) == [{"_id": "a", "value": 1}]
