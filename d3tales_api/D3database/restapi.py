import os
import warnings
import requests
from pathlib import Path

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


if __name__ == "__main__":
    # r = RESTAPI(method='post', endpoint='tools/upload/computation-gaussian', expected_endpoint="tools/user_uploads",
    #             url="https://d3tales.as.uky.edu", login_endpoint='login',
    #             upload_file="/mnt/research/~scratch~/gau_files.zip",
    #             params=dict(molecule_id='05DIRJ', calculation_type='opt_groundState'))
    p = RESTAPI(method='post', endpoint='tools/upload/e505a51f58ccc8550a772eadf59eeb18',
                expected_endpoint="tools/user_uploads",
                url="https://d3tales.as.uky.edu", login_endpoint='login',
                params=dict(approved='on'))
