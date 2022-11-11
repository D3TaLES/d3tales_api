from pymongo import MongoClient


class DBconnector:
    """
    Class to retrieve a collection from a database and insert new entry.
    Requires a db_infos.json file with credentials
    Copyright 2021, University of Kentucky
    """

    def __init__(self, db_info=None):

        self.host = db_info.get("host", )
        self.password = db_info.get("admin_password", )
        self.user = db_info.get("admin_user", )
        self.port = int(db_info.get("port", 0))
        self.database = db_info.get("database", )
        self.collection = db_info.get("collection", )

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
            raise Exception

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
