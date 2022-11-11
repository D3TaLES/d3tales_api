from __future__ import unicode_literals
import os
import time
import pymongo
import zipfile
import traceback
from atomate.utils.utils import env_chk
from d3tales_api.D3database.d3database import *
from d3tales_api.D3database.restapi import RESTAPI
from d3tales_api.Processors.d3tales_parser import *
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from d3tales_api.Processors.back2front import Gaus2FrontCharacterization

# Copyright 2021, University of Kentucky
TESTING = False

@explicit_serialize
class FileTransferTask(FiretaskBase):
    """
    A Firetask to Transfer files. Note that
    Required params:
        - mode: (str) - move, mv, copy, cv, copy2, copytree, copyfile, rtransfer
        - files: ([str]) or ([(str, str)]) - list of source files, or dictionary containing
                'src' and 'dest' keys
    Optional params:
        - dest: (str) destination directory, if not specified within files parameter (else optional)
        - server: (str) server host for remote transfer
        - user: (str) user to authenticate with on remote server
        - key_filename: (str) optional SSH key location for remote transfer
        - max_retry: (int) number of times to retry failed transfers; defaults to `0` (no retries)
        - retry_delay: (int) number of seconds to wait between retries; defaults to `10`
        - ignore_errors (bool): Optional. Whether to ignore errors. Defaults to False.
    """
    _fw_name = 'FileTransferTask'

    def run_task(self, fw_spec):
        ignore_errors = self.get('ignore_errors', False)
        max_retry = self.get('max_retry', 0)
        retry_delay = self.get('retry_delay', 10)
        skip_transfer = self.get('skip_transfer', False)
        method = self.get("method", "get")
        remote_server = self.get("remote_server", )
        key_file = env_chk(self.get("proc_vm_key"), fw_spec) or fw_spec["_fw_env"]["processing_vm_key"]
        remote_server["key_filename"] = key_file

        if TESTING:
            skip_transfer=True

        if skip_transfer:
            return None

        import paramiko
        ssh = paramiko.SSHClient()
        ssh.load_system_host_keys()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(**remote_server)
        sftp = ssh.open_sftp()

        for f in self["files"]:
            try:
                dest = self['dest']
                if method == "put":
                    sftp.put(f, dest)
                else:
                    sftp.get(f, dest)
                sftp.remove(f)

            except Exception:
                traceback.print_exc()
                if max_retry:

                    # we want to avoid hammering either the local or remote machine
                    time.sleep(retry_delay)
                    self['max_retry'] -= 1
                    self.run_task(fw_spec)

                elif not ignore_errors:
                    raise ValueError(
                        "There was an error performing operation {} from {} "
                        "to {}".format("rtansfer", self["files"], self["dest"]))

        sftp.close()
        ssh.close()


@explicit_serialize
class ProcessFile(FiretaskBase):
    _fw_name = "ProcessFile"

    def run_task(self, fw_spec):
        filepath = self["filepath"]
        mol_id = self["mol_id"]
        submission_data = self["submission_data"]
        metadata = self["metadata"]
        data_type = submission_data["data_type"]
        automatic_approval = self.get("automatic_approval", False)

        if TESTING:
            print("TESTING")
            os.chdir('/home/ubuntu/fireworks/zip_files/')

        with zipfile.ZipFile(filepath, "r") as target_zip:
            names = target_zip.namelist()
            submission_data["all_files_in_zip"] = names
            if data_type == 'gaussian':
                calculation_type = metadata.get('calculation_type', )
                if calculation_type == 'wtuning':
                    output_file = target_zip.extract(member="output.log")
                    mol_file = target_zip.extract(member="wtuning.com")
                    metadata["wtuning_output"] = output_file
                else:
                    logfile = [f for f in names if ".log" in f][0]
                    mol_file = target_zip.extract(member=logfile)
                metadata["mol_file"] = mol_file
                data = ProcessDFT(_id=mol_id, submission_info=submission_data, metadata=metadata).data_dict
                if automatic_approval:
                    Gaus2FrontCharacterization.from_data(data)
                    pass
            if data_type == 'cv':
                file = [f for f in names if ".txt" in f or "csv" in f][0]
                datafile = target_zip.extract(member=file)
                data = ProcessCV(datafile, _id=mol_id, submission_info=submission_data,
                                 metadata=metadata).data_dict
        target_zip.close()

        return FWAction(update_spec={"processed_data": data,
                                     "submission_data": submission_data,
                                     "hash_id": data["_id"],
                                     "mol_id": mol_id,
                                     "automatic_approval": automatic_approval
                                     })


@explicit_serialize
class SendToStorage(FiretaskBase):
    _fw_name = "SendToStorage"

    def run_task(self, fw_spec):
        max_retry = self.get('max_retry', 0)
        retry_delay = self.get('retry_delay', 10)
        ignore_errors = self.get('ignore_errors', False)

        filepath = self["filepath"]
        remote_server = self.get("remote_server", )
        hash_id = self.get("hash_id") or fw_spec.get("hash_id")
        storage_base_loc = env_chk(self.get("storage_base_loc"), fw_spec) or fw_spec["_fw_env"]["storage_base_loc"]
        key_file = env_chk(self.get("storage_key"), fw_spec) or fw_spec["_fw_env"].get("storage_key")
        if key_file:
            remote_server["key_filename"] = key_file

        submission_data = fw_spec["submission_data"]
        data_category = submission_data["data_category"]
        data_type = submission_data["data_type"]
        destination_path = "{}/d3tales/{}/{}/{}".format(storage_base_loc, fw_spec.get("mol_id"), data_category, data_type)
        destination = os.path.join(destination_path, filepath.split("/")[-1]).replace("//", "/")

        def mkdir_p(sftp, remote_directory):
            dir_path = str()
            for dir_folder in remote_directory.split("/"):
                if dir_folder == "":
                    continue
                dir_path += r"/{0}".format(dir_folder)
                try:
                    sftp.listdir(dir_path)
                except IOError:
                    sftp.mkdir(dir_path)

        import paramiko
        ssh = paramiko.SSHClient()
        ssh.load_system_host_keys()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        print('rs', remote_server)
        ssh.connect(**remote_server)
        sftp = ssh.open_sftp()

        try:
            mkdir_p(sftp, destination_path)
            sftp.put(filepath, destination)

        except Exception:
            traceback.print_exc()
            if max_retry:

                # we want to avoid hammering either the local or remote machine
                time.sleep(retry_delay)
                self['max_retry'] -= 1
                self.run_task(fw_spec)

            elif not ignore_errors:
                raise ValueError(
                    "There was an error performing operation {} from {} "
                    "to {}".format("rtansfer", self["files"], self["dest"]))

        # Update submission DB stored_location field
        dbc = DBconnector(db_info("backend"))
        coll = dbc.get_collection(data_category)
        identifier_dict = {"_id": hash_id}
        coll.update_one(identifier_dict, {"$set": {'submission_info.stored_location': destination}}, upsert=True)

        # Update processed data fireworks spec
        processed_data = fw_spec["processed_data"]
        processed_data["submission_info"].update({"stored_location": destination})

        if not TESTING:
            try:
                sftp.remove(filepath)
                pass
            except FileNotFoundError:
                pass
        sftp.close()
        ssh.close()

        return FWAction(update_spec={"processed_data": processed_data})


@explicit_serialize
class UpdateBackendDB(FiretaskBase):
    _fw_name = "UpdateBackendDB"

    def run_task(self, fw_spec):
        instance = fw_spec["processed_data"]
        submission_data = fw_spec["submission_data"]
        data_category = submission_data["data_category"]
        BackDB(collection_name=data_category, instance=instance)


@explicit_serialize
class UpdateUserDB(FiretaskBase):
    _fw_name = "UpdateUserDB"

    def run_task(self, fw_spec):
        processing_id = str(self.fw_id)
        hash_id = fw_spec["hash_id"]
        automatic_approval = fw_spec.get("automatic_approval") or self.get("automatic_approval", False)
        processed_data = json.loads(json.dumps(fw_spec["processed_data"]["data"]))

        dbc = DBconnector(db_info("frontend"))
        coll = dbc.get_collection('submission')
        try:
            _id = coll.find_one({"processing_id": processing_id})["_id"]
        except TypeError:
            return FWAction(update_spec={"already_updated": True})
        coll.update_one({"_id": _id}, {"$set": {"processing_id": hash_id}}, upsert=True)
        coll.update_one({"_id": _id}, {"$set": {"processed_data": processed_data}}, upsert=True)
        if automatic_approval:
            pass
            coll.update_one({"_id": _id}, {"$set": {"approved": True}}, upsert=True)


@explicit_serialize
class UpdateApprovalStatus(FiretaskBase):
    _fw_name = "UpdateApprovalStatus"

    def run_task(self, fw_spec):
        if self.get("already_updated", ):
            return None
        hash_id = self["hash_id"]
        approved_status = self['approved_status']
        data_category = self["data_category"]

        dbc = DBconnector(db_info("backend"))
        coll = dbc.get_collection(data_category)
        identifier_dict = {"_id": hash_id}
        coll.update_one(identifier_dict, {"$set": {'submission_info.approved': approved_status}}, upsert=True)


@explicit_serialize
class UpdateFrontendDB(FiretaskBase):
    _fw_name = "UpdateFrontendDB"

    def run_task(self, fw_spec):
        hash_id = self["hash_id"]
        data_category = self["data_category"]

        dbc = DBconnector(db_info("backend"))
        coll = dbc.get_collection(data_category)
        backend_data = coll.find_one({"_id": hash_id})

        data_type = backend_data.get('submission_info', {}).get("data_type", )
        data = backend_data.get("data", )
        mol_id = backend_data.get("mol_id", )
        calc_type = backend_data.get("calculation_type", )

        if data_type == 'gaussian':
            conditions = data.get("conditions", )
            charge = data.get("charge", )
            Gaus2FrontCharacterization(mol_id, calc_type, conditions, charge)


@explicit_serialize
class ApproveJobs(FiretaskBase):
    _fw_name = "ApproveJobs"

    def run_task(self, fw_spec):
        data_category = self["data_category"] or "computation"

        frontend_db = FrontDB(collection_name="submission")
        unapproved_cursor = frontend_db.coll.find({'$and': [
            {"approved": False},
            {"data_category": data_category},
            {"processed_data": {'$ne': {}}}
        ]
        })
        frontend_db = FrontDB(collection_name='submission')
        for mol in unapproved_cursor:
            h_id = mol.get('processing_id')
            while frontend_db.coll.count_documents({'processing_id': h_id}) > 1:
                frontend_db.coll.find_one_and_delete({'processing_id': h_id}, sort=[('_id', pymongo.ASCENDING)])
                print("---one submission document deleted for processing id {}".format(h_id))
            try:
                RESTAPI(method='post', endpoint='tools/upload/' + h_id, expected_endpoint="tools/user_uploads",
                        url="https://d3tales.as.uky.edu", login_endpoint='login', params=dict(approved='on'))
                print("Hash id {} approved".format(h_id))
            except Exception:
                print("ERROR approving hash id {}".format(h_id))

