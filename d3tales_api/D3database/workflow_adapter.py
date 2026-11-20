"""Database and website adapters used by Workflows.

Database-facing operations are kept under ``D3database`` and loaded lazily by
workflow helpers. Local workflows do not require credentials, network access,
or a database connection.
"""
import os
import warnings

# Copyright 2021-2026, University of Kentucky

D3TALES_URL = os.environ.get("D3TALES_URL", "https://d3tales.as.uky.edu")
FRONT_COLL = "base"


# ---------------------------------------------------------------------------
# REST API client (moved from D3database/restapi.py)
# ---------------------------------------------------------------------------
class RESTAPI(object):
    def __init__(self, method=None, url=D3TALES_URL, endpoint=None,
                 login_endpoint='login', username='', password='',
                 upload_file=None, params=None, expected_endpoint=None, return_json=False):
        """
        Upload a file to through d3tales.as.uky.edu file upload feature.
        :param method: str, html method (such as post or get)
        :param url: str, base url
        :param endpoint: str, post or get endpoint url (not containing base url)
        :param login_endpoint: str, login url (not containing base url)
        :param username: str, user username
        :param password: str, user password
        :param upload_file: str, path to file to be uploaded
        :param params: dict, form parameters for post
        :param return_json: bool, get or post method returns json if true

        Copyright 2021, University of Kentucky
        """
        self.method = method
        self.successful = False
        self.endpoint = "{}/{}/".format(url, endpoint).replace("//", "/").replace(':/', '://')
        print("Endpoint: ", self.endpoint)
        self.login_endpoint = "{}/{}/".format(url, login_endpoint).replace("//", "/").replace(':/', '://')
        if expected_endpoint:
            self.expected_endpoint = "{}/{}/".format(url, expected_endpoint).replace("//", "/").replace(':/', '://')
        default_username = username or os.environ.get('UPLOAD_USER') or os.getenv('UPLOAD_USER')
        default_password = password or os.environ.get('UPLOAD_PASS') or os.getenv('UPLOAD_PASS')
        self.user_data = dict(username=default_username,
                              password=default_password) if default_username and default_password else None
        self.client = self.get_client()
        params.update(dict(csrfmiddlewaretoken=self.csrftoken, next='/')) if params else {}
        self.params = params or {}
        self.upload_file = upload_file
        self.return_json = return_json

        if self.method in ["get", "GET", "Get"]:
            self.response = self.get_process()
        elif self.method in ["POST", "post", "Post"]:
            self.response = self.post_process()

        if expected_endpoint:
            if self.response.request.url == self.expected_endpoint:
                self.successful = True
            else:
                warnings.warn("The {} response url for {} to {} did not match the expected response url".format(
                    self.upload_file, self.endpoint, self.method))

    @property
    def cookies(self):
        return self.client.get(self.endpoint).cookies  # sets cookie

    @property
    def csrftoken(self):
        # Retrieve the CSRF token for data post
        return self.cookies['csrftoken'] if 'csrftoken' in self.cookies else self.cookies.get(['csrf'], )

    def get_client(self):
        import requests
        with requests.Session() as client:
            if self.login_endpoint and self.user_data:
                # Login
                client.get(self.login_endpoint)  # sets cookie
                csrftoken = client.cookies.get('csrftoken') or client.cookies.get('csrf')
                self.user_data.update(dict(csrfmiddlewaretoken=csrftoken, next='/'))
                # Submit login form
                client.post(self.login_endpoint, data=self.user_data, headers=dict(Referer=self.login_endpoint))
            return client

    def post_process(self):
        # Submit data form
        file_data = dict(file=open(self.upload_file, 'rb')) if self.upload_file else None
        req = self.client.post(self.endpoint, data=self.params, files=file_data,
                               headers=dict(Referer=self.endpoint), cookies=self.cookies)
        return req.json() if self.return_json else req

    def get_process(self):
        if self.params:
            req = self.client.get(self.endpoint, data=self.params, headers=dict(Referer=self.endpoint),
                                  cookies=self.cookies)
        else:
            req = self.client.get(self.endpoint, headers=dict(Referer=self.endpoint))
        return req.json() if self.return_json else req


# ---------------------------------------------------------------------------
# MongoDB connector (moved from D3database/d3database.py)
# ---------------------------------------------------------------------------
class DBconnector:
    """
    Class to retrieve a collection from a database and insert new entry.
    Requires a db_infos.json file with credentials
    Copyright 2021, University of Kentucky
    """

    def __init__(self, db_information: dict):
        self.db_information = db_information
        self.host = db_information.get("host", )
        self.password = db_information.get("admin_password", )
        self.user = db_information.get("admin_user", )
        self.port = int(db_information.get("port", 0))
        self.database = db_information.get("database", )
        self.collection = db_information.get("collection", )

    def get_database(self, **kwargs):
        from pymongo import MongoClient
        try:
            if "@" in self.host:
                conn = MongoClient(self.host)
            else:
                conn = MongoClient(host=self.host, port=self.port, username=self.user, password=self.password,
                                   **kwargs)
            db = conn[self.database]
        except Exception:
            raise ConnectionError("Error connecting to the MongoDB database with the following credentials: ",
                                  self.db_information)
        return db

    def get_collection(self, coll_name=None):
        db = self.get_database()
        if not coll_name:
            coll_name = self.collection
            if not coll_name:
                raise IOError("No collection specified")
        return db[coll_name]


# ---------------------------------------------------------------------------
# REST API helpers (read)
# ---------------------------------------------------------------------------
def _rest_get(endpoint):
    return RESTAPI(method='get', endpoint=endpoint, url=D3TALES_URL, return_json=True).response


def _mol_info_field(identifier, field):
    response = _rest_get("restapi/molecules/_id={}/mol_info.{}=1".format(identifier, field))
    if response:
        return response[0].get("mol_info", {}).get(field)
    return None


def d3t_get_smiles(identifier):
    """Return the SMILES stored in the D3TaLES frontend DB for `identifier` ('' if none)."""
    return _mol_info_field(identifier, 'smiles') or ''


def d3t_get_init_structure(identifier):
    """Return the initial structure (pymatgen Molecule) stored in the frontend DB. Raises LookupError if missing."""
    from pymatgen.core.structure import Molecule
    xyz = _mol_info_field(identifier, 'init_structure')
    if not xyz:
        raise LookupError("No init_structure found in D3TaLES for {}".format(identifier))
    return Molecule.from_str(xyz, 'xyz')


def d3t_get_ground_state(identifier, prop='charge'):
    """Return ground state charge/spin (int) from the frontend DB, or None if not found."""
    value = _mol_info_field(identifier, 'groundState_{}'.format(prop))
    return int(value) if value is not None else None


def d3t_get_omega(identifier):
    """Return the tuned omega value from the frontend DB, or None if not found."""
    response = _rest_get("restapi/molecules/_id={}/mol_characterization.omega=1".format(identifier))
    try:
        return response[0]["mol_characterization"]['omega']['value']
    except Exception:
        return None


def d3t_get_raw_computation(_hash, fields=None):
    """Return the `data` dict of a backend computation entry, or None if the entry does not exist."""
    endpoint = "restapi/rawdata/computation/_id={}".format(_hash)
    if fields:
        endpoint += "/{}=1".format(fields)
    response = _rest_get(endpoint)
    if response:
        return response[0].get("data", {})
    return None


# ---------------------------------------------------------------------------
# REST API / DB helpers (write)
# ---------------------------------------------------------------------------
def d3t_submit_calculation(zip_path, identifier, calc_type):
    """Upload a zipped calculation to the D3TaLES website for processing."""
    submission = RESTAPI(method='post', endpoint='tools/upload/computation-gaussian',
                         url=D3TALES_URL, expected_endpoint="tools/user_uploads",
                         upload_file=zip_path,
                         params=dict(molecule_id=identifier, calculation_type=calc_type))
    if not submission.successful:
        raise Exception("Calculation files not successfully submitted with endpoint {}. Response endpoint "
                        "was {}, not {}. \nSubmission Params: {} \n Zip path: {}".format(
                            submission.endpoint, submission.response.request.url, submission.expected_endpoint,
                            submission.params, zip_path))
    return submission


def _nested_set(coll, _id, path, value):
    """Same behavior as D3Database.insert(_id, nested=True) with override_lists=True."""
    if isinstance(value, dict):
        for k, v in value.items():
            _nested_set(coll, _id, "{}.{}".format(path, k), v)
    else:
        coll.update_one({"_id": _id}, {"$set": {path: value}}, upsert=True)


def _generate_mol_id(coll, smiles, group, source_groups):
    """Same behavior as FrontDB.generate_id."""
    from nanoid import generate
    id_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    existing = coll.find_one({"mol_info.smiles": smiles}, {"_id": 1})
    if existing:
        return existing["_id"]
    group_num = source_groups.get(group.title(), '10')
    _id = group_num + generate(id_chars, size=4)
    while coll.count_documents({"_id": _id}) > 0:
        _id = group_num + generate(id_chars, size=6)
    return _id


def d3t_register_molecule(smiles, names=None, origin_group='Risko', public=False):
    """
    Generate molecule info and insert it into the D3TaLES frontend DB.
    Equivalent to FrontDB(schema_layer='mol_info', instance=instance, smiles=smiles, group=group, public=public),
    which does not run schema validation when no _id is given.
    :return: (identifier, groundState_charge, groundState_spin)
    """
    from monty.json import jsanitize
    from d3tales_api.database_info import db_info, source_groups
    from d3tales_api.Processors.info_from_smiles import GenerateMolInfo  # Processors: next module to decouple
    instance = GenerateMolInfo(smiles, names=names or [], origin_group=origin_group).mol_info_dict

    coll = DBconnector(db_info().get("frontend")).get_collection(FRONT_COLL)
    _id = _generate_mol_id(coll, smiles, origin_group, source_groups())
    if isinstance(public, bool):
        coll.update_one({"_id": _id}, {"$set": {"public": public}}, upsert=True)
    _nested_set(coll, _id, "mol_info", jsanitize(instance, allow_bson=True))
    return _id, instance.get('groundState_charge'), instance.get('groundState_spin')

