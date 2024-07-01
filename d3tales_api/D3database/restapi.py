import os
import re
import warnings
import requests
import functools
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from pathlib import Path
from d3tales_api.D3database.d3database import FrontDB

BASE_DIR = Path(__file__).resolve().parent.parent
USERNAME = ''  # for pulling data
PASSWORD = ''


class RESTAPI(object):
    def __init__(self, method=None, url="https://d3tales.as.uky.edu", endpoint=None,
                 login_endpoint='login', username=USERNAME, password=PASSWORD,
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
        with requests.Session() as client:
            if self.login_endpoint and self.user_data:
                # Login
                client.get(self.login_endpoint)  # sets cookie
                csrftoken = client.cookies.get('csrftoken') or client.cookies.get('csrf')
                self.user_data.update(dict(csrfmiddlewaretoken=csrftoken, next='/'))
                # Submit login form
                req = client.post(self.login_endpoint, data=self.user_data, headers=dict(Referer=self.login_endpoint))
            return client

    def post_process(self):
        # Submit data form
        file_data = dict(file=open(self.upload_file, 'rb')) if self.upload_file else None
        req = self.client.post(self.endpoint, data=self.params, files=file_data,
                               headers=dict(Referer=self.endpoint), cookies=self.cookies)
        return_data = req.json() if self.return_json else req
        return return_data

    def get_process(self):
        if self.params:
            req = self.client.get(self.endpoint, data=self.params, headers=dict(Referer=self.endpoint),
                                  cookies=self.cookies)
        else:
            req = self.client.get(self.endpoint, headers=dict(Referer=self.endpoint))
        return_data = req.json() if self.return_json else req
        return return_data


class D3talesData:
    def __init__(self, username=USERNAME, password=PASSWORD):
        """
        This class pulls data from the D3Tales database and outputs plots or Pandas dataframes

        :param username: D3TaLES website username (must have REST API permissions)
        :param password: D3TaLES website password (must have REST API permissions)
        """
        self.username = username
        self.password = password

    def rgetkeys(self, _dict, keys, **kwargs):
        """
        Functions for getting property data

        :param _dict:
        :param keys:
        :return:
        """
        def _getkey(_dict, key):
            _dict = _dict or {}
            if isinstance(_dict, dict):
                return _dict.get(key, **kwargs)
            if isinstance(_dict, list) and key.isdigit():
                return _dict[int(key)]

        return functools.reduce(_getkey, [_dict] + keys.split('.'))
    def get_prop_data_raw(self, query, database='molecules', limit=0):
        """
        Get property data from D3TaLES database based on RESTAPI query

        :param query: str, D3TaLES REST API query
        :param database: str, name of database to query
        :param limit: limit query items to return

        :return: pandas DataFrame with query data
        """
        # Gather property data from REST API
        rest_query = re.split(r"\.0\.", query)[0].strip('.')
        response = RESTAPI(method='get',
                           endpoint="restapi/{}/{}==true/{}=1/limit={}".format(database, rest_query, rest_query, limit),
                           username=self.username, password=self.password,
                           url="https://d3tales.as.uky.edu", login_endpoint='login', return_json=True).response

        # Clean data
        data_df = pd.DataFrame(response)
        data_df.set_index('_id', inplace=True)
        return data_df

    def get_prop_data(self, query, max_cutoff=None, min_cutoff=None, database='molecules', limit=0):
        """
        Get property data from D3TaLES database based on RESTAPI query

        :param query: str, D3TaLES REST API query
        :param max_cutoff: float, maximum value to return for specified property
        :param min_cutoff: float, minimum value to return for specified property
        :param database: str, name of database to query
        :param limit: limit query items to return

        :return: pandas DataFrame with query data
        """
        split_query = re.split(r"\.0\.", query)
        clean_keys = '0.' + split_query[-1] if len(split_query) > 1 else None
        rest_query = split_query[0].strip('.')
        prop_category = rest_query.split('.')[0]
        prop_name = rest_query.split('.')[-1]
        column_name = rest_query.split('.')[1] + "_" + prop_name if prop_category == "species_characterization" else prop_name

        data_df = self.get_prop_data_raw(query=query, database=database, limit=limit)
        if clean_keys:
            data_df[column_name] = data_df[prop_category].apply(lambda x: self.rgetkeys(x.get(prop_name), clean_keys))

        data_df = data_df[[column_name]]
        data_df.dropna(inplace=True)
        data_df = data_df[pd.to_numeric(data_df[column_name], errors='coerce').notna()]

        # Remove outliers
        if pd.api.types.is_float_dtype(data_df[column_name]):
            data_df = data_df[(np.abs(stats.zscore(data_df)) < 3).all(axis=1)]  # drop outliers
            if min_cutoff:
                data_df = data_df[(data_df > min_cutoff).all(axis=1)]
            if max_cutoff:
                data_df = data_df[(data_df < max_cutoff).all(axis=1)]
        return data_df

    def get_master_df(self, master_fn='d3tales_props.csv'):
        """
        Get all major properties from D3TaLES database.

        :param master_fn: str, filepath to CSV file in which to save data
        uses the D3database FrontDB module by default, which required the DB_INFO_FILE to be defined to work.

        :return: pandas DataFrame with all data
        """
        props = [
            "mol_info.smiles",
            "mol_info.source_group",
            "mol_info.groundState_charge",
            "mol_info.number_of_atoms",
            "mol_info.molecular_weight",
            "mol_characterization.reorganization_energy.0.value",
            "mol_characterization.vertical_ionization_energy.0.value",
            "mol_characterization.vertical_ionization_energy_2.0.value",
            "mol_characterization.vertical_electron_affinity.0.value",
            "mol_characterization.redox_potential.0.value",
            "mol_characterization.rmsd_groundState_cation1.0.value",
            "mol_characterization.rmsd_cation1_cation2.0.value",
            "mol_characterization.omega.0.value",
            "species_characterization.ground_state.globular_volume.0.value",
            "species_characterization.ground_state.homo_lumo_gap.0.value",
            "species_characterization.ground_state.dipole_moment.0.value",
            "species_characterization.ground_state.solvation_energy.0.value",
            "species_characterization.cation1.globular_volume.0.value",
            "species_characterization.cation1.homo_lumo_gap.0.value",
            "species_characterization.cation1.dipole_moment.0.value",
            "species_characterization.cation1.solvation_energy.0.value",
            "species_characterization.cation2.globular_volume.0.value",
            "species_characterization.cation2.homo_lumo_gap.0.value",
            "species_characterization.cation2.dipole_moment.0.value",
            "species_characterization.cation2.solvation_energy.0.value",
        ]

        master = pd.DataFrame()
        for p in props:
            print("Getting Prop: ", p, "...")
            df = self.get_prop_data(p)
            master = pd.concat([master, df], axis=1)

        master.to_csv(master_fn)
        return master

    def hist_1d(self, query, **kwargs):
        """
        Plot histogram data from D3TaLES database based on RESTAPI query

        :param query: str, D3TaLES REST API query

        :return: seaborn histogram plot
        """
        df = self.get_prop_data(query, **kwargs)
        sns.histplot(data=df, x=df.columns[0])
        return df

    def hist_2d(self, query1, query2, db1='molecules', db2='molecules', **kwargs):
        """
        Plot histogram data from D3TaLES database based on RESTAPI query

        :param query1: str, D3TaLES REST API query for x axis
        :param query2: str, D3TaLES REST API query for y axis
        :param db1: str, name of database for query1
        :param db2: str, name of database for query2

        :return: seaborn 2D histogram plot
        """
        df1 = self.get_prop_data(query1, database=db1, **kwargs)
        df2 = self.get_prop_data(query2, database=db2, **kwargs)

        final_df = df1.join(df2, rsuffix='_2', lsuffix='_1')
        sns.histplot(data=final_df, x=final_df.columns[0], y=final_df.columns[1])
        return df1, df2


if __name__ == "__main__":
    # r = RESTAPI(method='post', endpoint='tools/upload/computation-gaussian', expected_endpoint="tools/user_uploads",
    #             url="https://d3tales.as.uky.edu", login_endpoint='login',
    #             upload_file="/mnt/research/~scratch~/gau_files.zip",
    #             params=dict(molecule_id='05DIRJ', calculation_type='opt_groundState'))
    # p = RESTAPI(method='post', endpoint='tools/upload/e505a51f58ccc8550a772eadf59eeb18',
    #             expected_endpoint="tools/user_uploads",
    #             url="https://d3tales.as.uky.edu", login_endpoint='login',
    #             params=dict(approved='on'))

    print(D3talesData().get_prop_data('mol_characterization.oxidation_potential.0.value', limit=2))
    # D3talesData().hist_1d('mol_characterization.oxidation_potential.0.value', min_cutoff=-10, max_cutoff=10)
    # master_df = D3talesData().get_master_df()
