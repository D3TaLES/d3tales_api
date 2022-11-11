import os
import json
import warnings
import jsonschema
import pandas as pd
from nanoid import generate
from dotty_dict import dotty
from monty.json import jsanitize
from d3tales_api.D3database.db_connector import DBconnector
from d3tales_api.D3database.schema2class import Schema2Class
from d3tales_api.D3database.info_from_smiles import GenerateMolInfo

db_file = os.environ.get('DB_INFO_FILE') or os.getenv('DB_INFO_FILE')
if not db_file:
    raise NameError("Environment variable DB_INFO_FILE not defined.")
with open(db_file, "r") as f:
    db_info = json.load(f)


class D3Database:
    # Copyright 2021, University of Kentucky
    def __init__(self, database=None, collection_name=None, instance=None, schema_layer="", schema_directory=None,
                 public=None, schema_db=None, default_id=None, validate_schema=True):
        self.collection_name = collection_name
        self.instance = {schema_layer: self.dot2dict(instance)} if schema_layer else self.dot2dict(instance)
        self.database = database
        self.public = public
        self.dbc = DBconnector(db_info(self.database))
        self.coll = self.dbc.get_collection(self.collection_name)
        schema_db = schema_db or database

        # validate
        if instance and validate_schema:
            self.instance['_id'] = self.instance.get("_id") or default_id
            self.s2c = Schema2Class(schema_name=collection_name, database=schema_db, schema_directory=schema_directory)
            jsonschema.validate(schema=self.s2c.schema, instance=self.instance)

    def insert(self, _id, nested=False, update_public=True, instance=None):
        """
        Upserts a dictionary into a frontend collection
        :param _id: str, _id for insertion
        :param nested: bool,
        :param update_public: bool,
        :param instance: dict,

        :return: None
        """
        if not instance:
            instance = jsanitize(self.instance, allow_bson=True)

        # Update public status
        if isinstance(self.public, bool) and update_public:
            self.coll.update_one({"_id": _id}, {"$set": {"public": self.public}}, upsert=True)

        for path, v in instance.items():
            if nested and isinstance(v, dict):
                for nest_k, nest_v in v.items():
                    new_path = ".".join(path.split(".") + [nest_k])
                    self.insert(_id, nested=True, update_public=False, instance={new_path: nest_v})

            if isinstance(v, list):
                self.array_checker(path, _id)
                self.coll.update_one({"_id": _id}, {"$addToSet": {path: {"$each": v}}}, upsert=True)
            else:
                self.coll.update_one({"_id": _id}, {"$set": {path: v}}, upsert=True)

        print("{} inserted into the {} database.".format(_id, self.database))

    def path_insert(self, _id, path, value):
        if isinstance(value, list):
            self.array_checker(path, _id)
            self.coll.update_one({"_id": _id}, {"$addToSet": {path: {"$each": value}}}, upsert=True)
        else:
            self.coll.update_one({"_id": _id}, {"$set": {path: value}}, upsert=True)

    @staticmethod
    def dot2dict(dot_dict):
        if isinstance(dot_dict, dict):
            final_dict = {}
            for k, v in dot_dict.items():
                if isinstance(k, str) and len(k.split(".")) > 1:
                    dot = dotty()
                    dot[k] = v
                    dot_dict = dot.to_dict()
                    final_dict.update(dot_dict)
                else:
                    final_dict[k] = v
            return final_dict
        return dot_dict

    def array_checker(self, field_path, _id):
        """
        Create empty array in filed if this field does not already exists. Prepare for set insertion (avoid duplicates)
        """
        print(_id, field_path)
        if not self.coll.count_documents({"_id": _id, field_path: {"$exists": True}}):
            self.coll.update_one({"_id": _id}, {"$set": {field_path: []}}, upsert=True)

    def checker(self, entry, field):
        """

        :param entry: value
        :param field: field name
        :return: {} or the match
        """
        result = self.coll.find_one({field: entry})
        return result if result else {}

    def make_query(self, query: dict = None, projection: dict = None, output="pandas", multi=True, limit=200):
        """
        Make MongoDB database query
        :param query: dict,
        :param projection: dict,
        :param output: str, output type
        :param multi: bool,
        :param limit: int, limit to the number of responses the query will return
        :return: 1) A dataframe if output="pandas"
                 2) A list if multi=False and a pymongo dursor if multi=True; output != "pandas
        """
        query = query or {}
        projection = projection or {}

        if multi:
            if projection:
                cursor = self.coll.find(query, projection).limit(limit)
            else:
                cursor = self.coll.find(query).limit(limit)
        else:
            if projection:
                cursor = [self.coll.find_one(query, projection).limit(limit)]
            else:
                cursor = [self.coll.find_one(query).limit(limit)]

        if output == "pandas":
            return pd.DataFrame.from_records(cursor)
        elif output == "json":
            json_data = json.loads(json.dumps(list(cursor), ))
            return json_data
        else:
            return list(cursor)


class FrontDB(D3Database):
    # Copyright 2021, University of Kentucky
    def __init__(self, schema_layer=None, instance=None, _id=None, smiles=None, group=None, collection_name="base",
                 generate_mol=False, public=None):
        super().__init__("frontend", collection_name=collection_name, instance=instance, schema_layer=schema_layer,
                         public=public, default_id=_id, validate_schema=True if _id else False)
        self.new_molecule = generate_mol

        # Identifier properties: must have either id or smiles and group
        self.smiles = smiles
        self.group = group

        if instance:
            self.id = _id if _id else self.generate_id()
            self.insert(self.id, nested=True)

        if self.new_molecule and schema_layer != "mol_info":
            self.id = _id if _id else self.generate_id()
            self.public = public
            instance = GenerateMolInfo(self.smiles, origin_group=self.group, database=self.database).mol_info_dict
            FrontDB(schema_layer='mol_info', instance=instance, _id=self.id)

    def check_if_in_db(self):
        if self.smiles:
            if self.coll.count_documents({"mol_info.smiles": self.smiles}) == 0:
                return False
            else:
                _id = self.coll.find_one({"mol_info.smiles": self.smiles})["_id"]
            return _id
        else:
            raise TypeError("If no smiles is given, molecule smiles must be given")

    def generate_id(self):
        group_dict = {
            "": '00',
            "unknown": '00',
            "Eaton": '01',
            "Robotics": '11',
            "Ganapathysubramanian": '02',
            "Jenkins": '03',
            "Mason": '04',
            "Odom": '05',
            "Odom_Hussein": '05',
            "Odom_Aman": '05',
            "Risko": '06',
            "Sarkar": '07',
            "Shaw": '08',
            "Teague": '09',
            "Csd": '80',
            "Zinc": '81',
            "Risko_Benchmark": '90',
            "Risko_Diheds": '91'
        }
        id_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        if self.smiles and self.group:
            group_num = group_dict.get(self.group.title(), '10')
            if self.coll.count_documents({"mol_info.smiles": self.smiles}) == 0:
                self.new_molecule = True
                rand = generate(id_chars, size=4)
                _id = "".join([group_num, rand])
                while self.coll.count_documents({"_id": _id}) > 0:
                    rand = generate(id_chars, size=6)
                    _id = "".join([group_num, rand])
            else:
                _id = self.coll.find_one({"mol_info.smiles": self.smiles})["_id"]
            return _id
        else:
            raise TypeError("If no id is given, molecule smiles and origin group must be given")


class BackDB(D3Database):
    # Copyright 2021, University of Kentucky
    def __init__(self, collection_name=None, instance=None, create_hash=False):
        super().__init__("backend", collection_name, instance)

        if instance:
            # Identifier properties: must have id and frontend instance
            self.id = self.instance.get("_id")
            if self.id is None and not create_hash:
                raise IOError("Hash ID is required for backend database insertion. Create an instance in the "
                              "frontend database first.")
            self.insert(self.id)

    def get_geometry(self, hash_id):
        array_data = self.coll.find_one({"_id": hash_id}, {"characterization.geometry": 1}
                                        )["characterization"]['geometry']
        for i in array_data:
            if i["source_hash_ids"] == hash_id:
                return i["sites"]


class ParamsDB(D3Database):
    # Copyright 2021, University of Kentucky
    def __init__(self, database='frontend', collection_name=None, instance=None, schema_directory=None, _id=None,
                 update=False):
        super().__init__(database, collection_name, instance, schema_directory=schema_directory, schema_db='frontend')
        self.update = update
        if instance:
            # Identifier properties: must have id and frontend instance
            try:
                id_to_check = _id or self.instance["_id"]
                self.id = self.check_id(id_to_check)
                self.instance["_id"] = self.id
            except KeyError:
                raise IOError("Parameter ID (common name) is required for parameter database insertion.")
            self.insert(self.id)

    def check_id(self, id_to_check):
        if self.coll.count_documents({"_id": id_to_check}) == 0:
            if self.update:
                raise IOError(
                    "There exists no parameter entry with the id {}. Insert a new entry or change update to False.".format(
                        id_to_check))
            return id_to_check
        else:
            if not self.update:
                id_addon = 1
                _id = id_to_check
                while self.coll.count_documents({"_id": _id}) > 0:
                    _id = "_".join([id_to_check, str(id_addon)])
                    id_addon += 1
                warnings.warn(
                    "A parameter already exists with the id {}. A new entry was created with the id {}".format(
                        id_to_check, _id))
                return _id
            return id_to_check


class Expflow(D3Database):
    # Copyright 2021, University of Kentucky
    def __init__(self, collection_name='experimentation', instance=None, create_hash=False):
        super().__init__("expflow", collection_name, instance, schema_db='backend')
        if instance:
            # Identifier properties: must have id and frontend instance
            self.id = self.instance.get("_id")
            if self.id is None and not create_hash:
                raise IOError("Hash ID is required for database insertion. Be sure to create an instance in the "
                              "molecule database first.")
            self.insert(self.id)


class QueryDB:
    """
    class to query the database. Use mongoDB style query and projections.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, query: dict = None, projection: dict = None, collection_name=None, db_file=None, output="pandas",
                 multi=True):
        self.dbc = DBconnector(db_file)
        self.query = query or {}
        self.collection_name = collection_name
        self.output = output
        self.projection = projection or {}
        self.multi = multi

    def make_query(self):
        """

        :return: 1) A dataframe if output="pandas"
                 2) A list if multi=False and a pymongo dursor if multi=True; output != "pandas

        """
        coll = self.dbc.get_collection(self.collection_name)

        if self.multi:
            if self.projection:
                cursor = coll.find(self.query, self.projection)
            else:
                cursor = coll.find(self.query)
        else:
            if self.projection:
                cursor = [coll.find_one(self.query, self.projection)]
            else:
                cursor = [coll.find_one(self.query)]

        if self.output == "pandas":
            return pd.DataFrame.from_records(cursor)
        elif self.output == "json":
            return json.dumps(list(cursor))
        else:
            return cursor


if __name__ == "__main__":
    FrontDB(instance={"experiment_data.experiment_ids": ["test"]}, _id="06TNKR")
