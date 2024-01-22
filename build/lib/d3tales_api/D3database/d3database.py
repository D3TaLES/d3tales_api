import json
import warnings
import jsonschema
import pandas as pd
from nanoid import generate
from dotty_dict import dotty
from datetime import datetime
from monty.json import jsanitize
from d3tales_api.database_info import db_info, source_groups
from d3tales_api.D3database.schema2class import Schema2Class
from d3tales_api.D3database.info_from_smiles import GenerateMolInfo

from pymongo import MongoClient


class DBconnector:
    """
    Class to retrieve a collection from a database and insert new entry.
    Requires a db_infos.json file with credentials
    Copyright 2021, University of Kentucky
    """

    def __init__(self, db_information: dict):
        """
        :param db_information: dictionary of database info
        """

        self.host = db_information.get("host", )
        self.password = db_information.get("admin_password", )
        self.user = db_information.get("admin_user", )
        self.port = int(db_information.get("port", 0))
        self.database = db_information.get("database", )
        self.collection = db_information.get("collection", )

    def get_database(self, **kwargs):
        """
        Returns a database object

        :return: a database object
        """
        try:
            if "@" in self.host:
                conn = MongoClient(self.host)
            else:
                conn = MongoClient(host=self.host, port=self.port,
                                   username=self.user,
                                   password=self.password, **kwargs)

            db = conn[self.database]
        except:
            raise ConnectionError

        return db

    def get_collection(self, coll_name=None):
        """
        Returns a collection from the database

        :param coll_name: name of collection
        :return: db.collection
        """
        db = self.get_database()
        if not coll_name:
            coll_name = self.collection
            if not coll_name:
                raise IOError("No collection specified")

        return db[coll_name]


class D3Database:
    """
    Base class for connecting to a database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, database=None, collection_name=None, instance=None, schema_layer="", schema_directory=None,
                 public=None, schema_db=None, default_id=None, validate_schema=True):
        """
        :param database: str, name of database (should be a key in the DB_INFO_FILE)
        :param collection_name: str, name of collection
        :param instance: dict, instance to insert or validate
        :param schema_layer: str, schema layer
        :param schema_directory: str, schema directory
        :param public: bool, instance is marked as public if True
        :param schema_db: str, database containing schema
        :param default_id: str, default instance ID
        :param validate_schema: bool, validate schema if True
        """
        self.collection_name = collection_name
        self.instance = {schema_layer: self.dot2dict(instance)} if schema_layer else self.dot2dict(instance)
        self.database = database
        self.public = public
        self.dbc = DBconnector(db_info.get(self.database))
        self.coll = self.dbc.get_collection(self.collection_name)
        schema_db = schema_db or database

        # validate
        if instance and validate_schema:
            self.instance['_id'] = self.instance.get("_id") or default_id
            self.s2c = Schema2Class(schema_name=collection_name, database=schema_db, schema_directory=schema_directory)
            jsonschema.validate(schema=self.s2c.schema, instance=self.instance)

    def insert(self, _id, nested=False, update_public=True, instance=None, override_lists=False):
        """
        Upserts a dictionary into a collection

        :param _id: str, _id for insertion
        :param nested: bool, insert nested attributes if True
        :param update_public: bool, update the public status if true
        :param instance: dict, instance to be inserted
        :param override_lists: bool, override existing lists in insertion if True
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

            elif isinstance(v, list) and not override_lists:
                self.array_checker(path, _id)
                self.coll.update_one({"_id": _id}, {"$addToSet": {path: {"$each": v}}}, upsert=True)
            else:
                self.coll.update_one({"_id": _id}, {"$set": {path: v}}, upsert=True)

        print("{} {}... inserted into the {} database.".format(_id, str(instance)[:15], self.database))

    def path_insert(self, _id, path, value):
        """
        Insert a piece of data to a specific path, updating array if the path leads to an array

        :param _id: str, instance ID
        :param path: str, path to insertion
        :param value: value to insert
        """
        if isinstance(value, list):
            self.array_checker(path, _id)
            self.coll.update_one({"_id": _id}, {"$addToSet": {path: {"$each": value}}}, upsert=True)
        else:
            self.coll.update_one({"_id": _id}, {"$set": {path: value}}, upsert=True)

    def array_checker(self, field_path, _id):
        """
        Create empty array in filed if this field does not already exists. Prepare for set insertion (avoid duplicates)

        :param field_path: str, path to check
        :param _id: str, instance ID
        """
        if not self.coll.count_documents({"_id": _id, field_path: {"$exists": True}}):
            self.coll.update_one({"_id": _id}, {"$set": {field_path: []}}, upsert=True)

    def field_checker(self, entry, field):
        """
        Check if field exists and return result or {}

        :param entry: value
        :param field: field name
        :return: {} or the match
        """
        result = self.coll.find_one({field: entry})
        return result if result else {}

    def make_query(self, query: dict = None, projection: dict = None, output="pandas", multi=True, limit=200):
        """
        Make MongoDB database query

        :param query: dict, query
        :param projection: dict, projection
        :param output: str, output type
        :param multi: bool, return multiple query responses if True
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
                cursor = [self.coll.find_one(query, projection)]
            else:
                cursor = [self.coll.find_one(query)]

        if output == "pandas":
            return pd.DataFrame.from_records(cursor)
        elif output == "json":
            json_data = json.loads(json.dumps(list(cursor), ))
            return json_data
        else:
            return list(cursor)

    @staticmethod
    def dot2dict(dot_dict):
        """
        Convert a dotted dictionary to a nested dictionary

        :param dot_dict: dotted dictionary
        :return: nested final dictionary
        """
        dot_dict = dot_dict or {}
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


class FrontDB(D3Database):
    """
    Class for accessing the frontend D3TaLES database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, schema_layer=None, instance=None, _id=None, smiles=None, group=None, collection_name="base",
                 generate_mol=False, public=None):
        """
        :param schema_layer: str, schema layer
        :param instance: dict, instance to insert or validate
        :param _id: str, default instance ID
        :param collection_name: str, name of collection
        :param public: bool, instance is marked as public if True
        :param smiles: str, SMILES string for a new molecule
        :param group: str, origin group for instance
        :param generate_mol: bool, generate new molecule instance if True
        """
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
        """
        Check if molecule with object SMILES string exists in the database

        :return: molecule ID if exists, else False
        """
        if self.smiles:
            if self.coll.count_documents({"mol_info.smiles": self.smiles}) == 0:
                return False
            else:
                _id = self.coll.find_one({"mol_info.smiles": self.smiles})["_id"]
            return _id
        else:
            raise TypeError("If no smiles is given, molecule smiles must be given")

    def generate_id(self):
        """
        Generate a new molecule ID with origin group

        :return: new molecule ID
        """
        id_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        if self.smiles and self.group:
            group_num = source_groups.get(self.group.title(), '10')
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
    """
    Class for accessing the backend D3TaLES database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, collection_name=None, instance=None, create_hash=False, last_updated=False):
        """
        :param collection_name: str, name of collection
        :param instance: dict, instance to insert or validate
        :param create_hash: bool, generate new instance hash id if True
        :param last_updated: bool, insert attribute last_updated to instance with the current date and time if True
        """
        super().__init__("backend", collection_name, instance)

        if instance:
            # Identifier properties: must have id and frontend instance
            self.id = self.instance.get("_id")
            if self.id is None and not create_hash:
                raise IOError("Hash ID is required for backend database insertion. Create an instance in the "
                              "frontend database first.")
            if last_updated:
                self.instance.update({"last_updated": datetime.now().isoformat()})
            self.insert(self.id)

    def get_geometry(self, hash_id):
        """
        Get molecule geometry form hash ID

        :param hash_id: ID
        :return: molecule geometry sites
        """
        array_data = self.coll.find_one({"_id": hash_id}, {"characterization.geometry": 1}
                                        )["characterization"]['geometry']
        for i in array_data:
            if i["source_hash_ids"] == hash_id:
                return i["sites"]


class ParamsDB(D3Database):
    """
    Class for accessing the parameters D3TaLES database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, database='frontend', collection_name=None, instance=None, schema_directory=None, _id=None,
                 update=False):
        """
        :param database: str, name of database (should be a key in the DB_INFO_FILE)
        :param collection_name: str, name of collection
        :param instance: dict, instance to insert or validate
        :param schema_directory: str, schema directory
        :param _id: str, parameter ID
        :param update: bool, create new parameter ID if False
        """
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
        """
        Check if parameter ID already exists, create new one if it already exists and object is not set to update
        Copyright 2021, University of Kentucky

        :param id_to_check: parameter ID to check
        :return: parameter ID (checked ID or new ID)
        """
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
    """
    Class for accessing the ExpFlow experimental database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, collection_name='experimentation', instance=None, create_hash=False):
        """
        :param collection_name: str, name of collection
        :param instance: dict, instance to insert or validate
        :param create_hash: bool, generate new instance hash id if True
        """
        super().__init__("expflow", collection_name, instance, schema_db='backend')
        if instance:
            # Identifier properties: must have id and frontend instance
            self.id = self.instance.get("_id")
            if self.id is None and not create_hash:
                raise IOError("Hash ID is required for database insertion. Be sure to create an instance in the "
                              "molecule database first.")
            self.insert(self.id)


class RobotStatusDB(D3Database):
    """
    Class for accessing the Robot Status database
    Copyright 2021, University of Kentucky
    """

    def __init__(self, apparatus_type: str, _id: str = None, instance: dict = None, override_lists: bool = True,
                 wflow_name: str = None):
        """
        Initiate class
        :param apparatus_type: str, name of apparatus
        :param _id: str, _id
        :param instance: dict, instance to insert or validate
        :param override_lists: bool,
        :param wflow_name: str, name of active workflow; checks if database instance has appropriate wflow_name if set
        """
        super().__init__("robotics", 'status_' + apparatus_type, instance, schema_db='robot')
        self.id = _id or self.instance.get("_id")
        self.wflow_name = wflow_name
        if instance:
            instance["_id"] = self.id
            if not self.id:
                raise IOError("ID is required for {} status database insertion.".format(apparatus_type))
            self.insert(self.id, override_lists=override_lists)
        if wflow_name and self.id:
            self.check_wflow_name()

    def check_wflow_name(self):
        current_wflow = self.coll.find_one({"_id": self.id}).get("current_wflow_name")
        if current_wflow == self.wflow_name:
            return True
        raise NameError("Argument wflow_name ({}) does not match instance current_wflow_name {}.".format(self.wflow_name, current_wflow))

    def get_prop(self, prop: str):
        """
        Get database prop for instance with _id
        :param prop: str, property name
        :return: prop
        """
        return self.coll.find_one({"_id": self.id}).get(prop)

    def update_status(self, new_status, status_name="location"):
        """
        Update status for a vial location or station vial
        :param new_status: str, new vial location or new vial in station
        :param status_name: str, name for status properties
        """
        current_status = self.get_prop("current_" + status_name)
        history = self.get_prop(status_name + "_history")
        history.append(current_status)
        self.insert(self.id, override_lists=True, instance={
            "current_" + status_name: new_status,
            status_name + "_history": history
        })


class QueryDB:
    """
    Class to query the database. Use mongoDB style query and projections.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, query: dict = None, projection: dict = None, collection_name=None, db_file=None, output="pandas",
                 multi=True):
        """
        :param query: dict, query
        :param projection: dict, projection
        :param collection_name: str, name of collection
        :param db_file: str, database info file
        :param output: str, output type
        """
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
