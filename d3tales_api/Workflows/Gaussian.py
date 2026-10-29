import abc
import subprocess
import multiprocessing
import uuid
from six import add_metaclass
from atomate.utils.utils import get_logger, env_chk
from fireworks import FiretaskBase, explicit_serialize, FWAction
from d3tales_api.Workflows.utils import *
from d3tales_api.Workflows.wtuning import WtuningJob
from d3tales_api.Workflows.d3tales_utils import d3t_submit_calculation

logger = get_logger(__name__)
cpus = [multiprocessing.cpu_count() if multiprocessing.cpu_count() < 16 else 16]
nprocs = str(cpus[0])


# Copyright 2021, University of Kentucky


@add_metaclass(abc.ABCMeta)
class GaussianBase(FiretaskBase):
    _fw_name = "GaussianBase"

    def setup_files(self, fw_spec, calc_type='opt'):
        # get parameters
        self.data_store = get_data_store(self, fw_spec)
        self.gaussian_cmd = env_chk(self.get("g16_cmd"), fw_spec)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.full_name = name_tag + self['name']
        self.calc_name = self['name']
        self.paramset = self["paramset"]
        self.identifier = fw_spec.get("identifier", ) or self.get("identifier", uuid.uuid4().__str__())
        self.smiles = fw_spec.get("smiles", ) or self.get("smiles")
        self.runfile_log = env_chk(self.get('runfile_log'), fw_spec)
        self.skip_freq = fw_spec.get("skip_freq", ) or self.get("skip_freq") or False
        self.skip_freq = True if self.calc_name == "opt_mol" else self.skip_freq
        self.check_if_already_run = spec_or_task(self, fw_spec, "check_if_already_run", True)
        self.submit = spec_or_task(self, fw_spec, "submit", True)
        self.restricted = spec_or_task(self, fw_spec, "restricted", True)
        self.run_nto = fw_spec.get("run_nto") or self.get("run_nto") or False

        self.gaussian_file_name = fw_spec.get("gaussian_file_name") or self.get("gaussian_file_name") or "gaussian"

        # set prefix for the file names
        if self.get("prefix"):
            prefix = name_tag + self["prefix"]
        elif self.get("subtype", ):
            prefix = name_tag + calc_type + "_" + self["subtype"]
        else:
            prefix = name_tag + calc_type + "_mol"

        # get calc directory
        if self.get("path", ):
            self.calc_dir = "{}/{}/{}".format(env_chk(self.get('path'), fw_spec), self.identifier,
                                              self.gaussian_file_name)
        else:
            self.calc_dir = self.get("calc_dir", os.getcwd())

        # folder for local JSON records (local data store only)
        self.store_dir = get_store_dir(self.identifier, get_local_output_dir(self, fw_spec), self.calc_dir)

        # get the type of gaussian optimization calculation
        if self.get("type", ):
            self.type = self["type"]
            if self.type == "coupling":
                self.working_dir = '{}/{}/{}/{}/{}'.format(self.calc_dir, self.type, calc_type, self["label"],
                                                           self["subtype"])
            elif self.get("label"):
                self.working_dir = '{}/{}/{}/{}'.format(self.calc_dir, self.type, calc_type, self["label"])
            elif self.type.startswith("solv_"):
                self.working_dir = '{}/{}/{}/{}'.format(self.calc_dir, self.type, calc_type, self["subtype"])
            else:
                self.working_dir = '{}/{}/{}'.format(self.calc_dir, self.type, calc_type)
        else:
            self.type = 'mol_{}'.format(calc_type)
            self.working_dir = "{}/{}".format(self.calc_dir, self.type)

        # Set default files
        self.file_com = "{}/{}.com".format(self.working_dir, prefix)
        self.file_chk = "{}/{}.chk".format(self.working_dir, prefix)
        self.file_fchk = "{}/{}.fchk".format(self.working_dir, prefix)
        self.file_log = "{}/{}.log".format(self.working_dir, prefix)
        self.freq_com = "{}/freq_{}.com".format(self.working_dir, prefix)
        self.freq_chk = "{}/freq_{}.chk".format(self.working_dir, prefix)
        self.freq_fchk = "{}/freq_{}.fchk".format(self.working_dir, prefix)
        self.freq_log = "{}/freq_{}.log".format(self.working_dir, prefix)

        # create folders
        os.makedirs(self.working_dir, exist_ok=True)
        os.chdir(self.working_dir)

    def setup_calc(self, fw_spec, calc_type='opt'):
        self.setup_files(fw_spec, calc_type=calc_type)
        name_tag = fw_spec.get("name_tag", ) or self.get("name_tag") or ""
        self.solvent = fw_spec.get("solvent", ) or self.get("solvent", )
        self.gs_charge = fw_spec.get("gs_charge") or self.get("gs_charge") or get_groundState(
            self.identifier, self.smiles, data_store=self.data_store) or 0
        self.gs_spin = fw_spec.get("gs_spin") or self.get("gs_spin") or get_groundState(
            self.identifier, self.smiles, prop='spin', data_store=self.data_store) or 1
        use_iop = fw_spec.get("use_iop", True) if self.get("use_iop", True) else self.get("use_iop")
        run_from_com = fw_spec.get("run_from_com") or self.get("run_from_com") or False

        # get iop string (omega value) for the calculation
        self.iop_str = fw_spec.get("iop_str", ) or self.get("iop_str", )
        if not self.iop_str and use_iop:
            tuned_w = get_omega(self.identifier, data_store=self.data_store, store_dir=self.store_dir)
            if tuned_w:
                self.iop_str = omega_to_iop(tuned_w)
        self.result_iop_str = self.iop_str  # fixed copy used for local record hashes

        # get starting geometry for calculation
        geometry = fw_spec.get("geometry", ) or self.get("geometry", )
        geometry_sites = fw_spec.get("{}_geom".format(geometry), )
        if run_from_com:
            self.mol = GaussianInput.from_file(self.file_com).molecule
        else:
            if geometry_sites:
                self.mol = Molecule.from_sites([Site.from_dict(sd) for sd in geometry_sites])
            else:
                geometry_hash = fw_spec.get("{}_hash".format(geometry),
                                            orig_hash_id(self.identifier, self.calc_name, self.paramset.functional,
                                                         self.paramset.basis_set, tuning_parameter=self.iop_str,
                                                         solvent=self.solvent))
                self.mol = get_db_geom(geometry_hash, data_store=self.data_store, store_dir=self.store_dir) \
                    or start_from_smiles(self.identifier, self.smiles, data_store=self.data_store)

            # End job if the total number of atoms is greater than 200
            num_atoms = len(self.mol.sites)
            if num_atoms > 200:
                return FWAction(defuse_children=True,
                                update_spec={"defuser_reason": "The molecule has {} atoms".format(num_atoms)})

        # Update route parameters with omega value and pcm solvent
        radical_electrons = (self.gs_spin - 1 - self.paramset.charge) % 2
        self.paramset.multiplicity = radical_electrons + 1  # calculate spin multiplicity with Hand's rule
        self.paramset.charge += self.gs_charge
        # Update functional for unrestricted calculation
        self.functional = self.paramset.functional if self.restricted else "R{}".format(self.paramset.functional)
        print("Charge: ", self.paramset.charge, "\t Multiplicity: ", self.paramset.multiplicity)
        self.paramset.link0_parameters.update({"%chk": self.file_chk, "%mem": "48GB", "%nprocshared": nprocs})
        if self.iop_str and use_iop:
            self.paramset.route_parameters.update({"iop(3/107={}, 3/108={})".format(self.iop_str, self.iop_str): ""})
        if self.solvent:
            self.paramset.route_parameters.update({"SCRF": "(PCM,Solvent={})".format(self.solvent)})

    def local_conditions(self):
        return orig_conditions(self.paramset.functional, self.paramset.basis_set,
                               tuning_parameter=self.result_iop_str, solvent=self.solvent)

    def local_key(self, calc_name=None):
        # key for local (database-free) records; same hash scheme as the D3TaLES backend
        return orig_hash_id(self.identifier, calc_name or self.calc_name, self.paramset.functional,
                            self.paramset.basis_set, tuning_parameter=self.result_iop_str, solvent=self.solvent)

    def existing_data(self, _hash=None, calc_name=None):
        # check if this job has already run and been stored in the data store
        if self.data_store == "local":
            return get_existing_data(self.local_key(calc_name), data_store=self.data_store, store_dir=self.store_dir)
        _hash = _hash or orig_hash_id(self.identifier, self.calc_name, self.paramset.functional,
                                      self.paramset.basis_set,
                                      tuning_parameter=self.iop_str, solvent=self.solvent)
        return get_existing_data(_hash, data_store=self.data_store)

    def post_job(self, upload_files=None, delete_files=None, calc_type=None, name_tag=''):
        # Store data (D3TaLES: upload through website; local: write JSON record)
        upload_files = upload_files or [self.file_log, self.file_fchk]
        delete_files = []  # delete_files or [self.file_chk, self.file_com, self.file_fchk, self.file_log]
        calc_type = calc_type or self.calc_name
        _hash = get_hash_id(self.identifier, self.file_log, self.calc_name) if os.path.isfile(self.file_log) else None
        if self.data_store == "d3tales":
            upload_names = [f.split('/')[-1] for f in upload_files]
            zip_path = zip_files(upload_names, zip_name='{}_{}.zip'.format(self.identifier, name_tag + self.full_name))
            if self.submit:
                d3t_submit_calculation(zip_path, self.identifier, calc_type)
                print("File {}_{}.zip successfully posted!".format(self.identifier, name_tag + self.full_name))
        else:
            from d3tales_api.Processors.parser_dft import ProcessGausLog
            logfile = next((f for f in upload_files if f.endswith(".log")), self.file_log)
            record = ProcessGausLog(_id=self.identifier, metadata={
                "mol_file": logfile, "calculation_type": calc_type,
                "mol_info": {"groundState_charge": self.gs_charge}}).data_dict
            print("Record saved to", save_local_result(self.store_dir, record, key=self.local_key(calc_type)))
        # Write runfile to runfile_log so the runfile can be deleted after calculation
        with open(self.runfile_log, 'a') as fn:
            fn.write("{}, {}\n".format(_hash, self.working_dir))
        # TODO write script that checks process status, approves, and deletes run dir
        # Remove excess files
        os.chdir(self.working_dir)
        if delete_files:
            for file in delete_files:
                os.system("rm -fr {}".format(file))


@explicit_serialize
class RunGaussianEnergy(GaussianBase):

    def run_task(self, fw_spec):
        setup_obj = self.setup_calc(fw_spec, calc_type='energy')
        if isinstance(setup_obj, FWAction):
            return setup_obj

        # generate the input for gaussian energy run
        gauss_inp = generate_gaussian_input(paramset=self.paramset, mol=self.mol)
        gauss_inp.write_file(self.file_com, cart_coords=True)
        if self.check_if_already_run:
            existing_data = self.existing_data()
            if existing_data:
                return FWAction(
                    update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier,
                                 "gs_charge": self.gs_charge, "gs_spin": self.gs_spin,
                                 "{}_eng".format(self.full_name): existing_data.get("scf_total_energy", {}).get(
                                     "value"),
                                 "{}_hash".format(self.full_name): get_hash_id(self.identifier, self.file_com,
                                                                               self.calc_name),
                                 "iop_str": self.iop_str})

        # run gaussian16
        print('SUBMITTING GAUSSIAN  ENERGY JOB {} for {}'.format(self.full_name, self.identifier))
        subprocess.call(self.gaussian_cmd + " " + self.file_com, shell=True)

        # check for normal termination of gaussian job
        gout = GaussianOutput(self.file_log)
        final_energy = gout.final_energy
        if not gout.properly_terminated:
            raise RuntimeError(self.file_com + " not terminated. calc_dir " + self.working_dir)
        else:
            return_code = subprocess.call("formchk " + self.file_chk, shell=True)

            if return_code != 0:
                raise RuntimeError("formatted checkpoint file could not be created")

            # Clean up files and transfer to website processing
            self.post_job(delete_files=[self.file_log, self.file_chk, self.file_fchk, self.file_com])
        return FWAction(
            update_spec={"gaussrun_dir": self.working_dir, "identifier": self.identifier, "iop_str": self.iop_str,
                         "gs_charge": self.gs_charge, "gs_spin": self.gs_spin,
                         "{}_eng".format(self.full_name): final_energy})


@explicit_serialize
class RunGaussianOpt(GaussianBase):

    def run_task(self, fw_spec):
        setup_obj = self.setup_calc(fw_spec, calc_type='opt')
        if isinstance(setup_obj, FWAction):
            return setup_obj

        # run gaussian optimize until normal termination is achieved. max runs is 5
        runs = 0
        converged = False
        final_structure, final_energy, gibbs_correction = None, None, None
        freq_name = "freq_{}".format(self.full_name.split('_')[-1])

        while not converged and runs < 5:
            runs += 1

            # write input files for gaussian optimization calculation
            gauss_inp = generate_gaussian_input(paramset=self.paramset, mol=self.mol)
            gauss_inp.write_file(self.file_com, cart_coords=True)
            opt_hash = get_hash_id(self.identifier, self.file_com, self.calc_name)
            freq_hash = get_hash_id(self.identifier, self.file_com, freq_name)
            if self.check_if_already_run:
                existing_data = self.existing_data()
                if existing_data:
                    return FWAction(
                        update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier,
                                     "gs_charge": self.gs_charge, "gs_spin": self.gs_spin,
                                     "{}_hash".format(self.full_name): opt_hash,
                                     "{}_hash".format(freq_name): freq_hash,
                                     "{}_geom".format(self.full_name): existing_data.get("geometry"),
                                     "{}_eng".format(self.full_name): existing_data.get("scf_total_energy", {}).get(
                                         "value"),
                                     "{}_gibb".format(self.full_name): self.existing_data(
                                         _hash=freq_hash, calc_name=freq_name).get(
                                         "gibbs_correction", {}).get("value"),
                                     "iop_str": self.iop_str})

            # run gaussian optimization
            if not os.path.isfile(self.file_log) or not os.path.isfile(self.file_chk) or runs != 1:
                print('SUBMITTING GAUSSIAN OPT JOB {} for {}'.format(self.full_name, self.identifier))
                subprocess.call(self.gaussian_cmd + " " + self.file_com, shell=True)
            gout = GaussianOutput(self.file_log)

            # check for normal termination
            if not gout.properly_terminated:
                try:
                    self.mol = gout.final_structure
                    continue
                except IndexError:  # if previous job never even finished 1 structure
                    continue
            return_code = subprocess.call("formchk " + self.file_chk, shell=True)
            if return_code != 0:
                raise RuntimeError("formatted checkpoint file for opt calculation could not be created")

            if self.skip_freq:
                final_structure = gout.final_structure.as_dict()['sites']
                final_energy = gout.final_energy
                converged = True
                break

            # generate input for frequency calculation
            ginp = gout.to_input()

            # update route parameters
            if ginp.route_parameters.get("iop(3107", ):
                self.iop_str = ginp.route_parameters["iop(3107"].split(",")[0]
                ginp.route_parameters = {"iop(3/107={}, 3/108={})".format(self.iop_str, self.iop_str): "", "freq": "",
                                         "Geom": "AllCheck", "Guess": "TCheck", "SCRF": "Check"}
            else:
                ginp.route_parameters = {"freq": "", "Geom": "AllCheck", "Guess": "TCheck", "SCRF": "Check"}
            ginp.link0_parameters.update({"%chk": self.freq_chk, "%oldchk": self.file_chk})
            ginp.dieze_tag = '#P'
            ginp.write_file(self.freq_com, cart_coords=True)

            # run gaussian frequency
            if not os.path.isfile(self.freq_log) or not os.path.isfile(self.freq_chk) or runs != 1:
                print('SUBMITTING GAUSSIAN FREQUENCY JOB {} for {}'.format(self.full_name, self.identifier))
                subprocess.call("{} {} > {}".format(self.gaussian_cmd, self.freq_com, self.freq_log), shell=True)
            fgout = GaussianOutput(self.freq_log)

            if not fgout.properly_terminated:
                raise RuntimeError("Frequency calculation failed Normal termination")

            elif fgout.frequencies[0][0]["frequency"] < 0:
                self.mol.perturb(0.1)  # perturb structure to find structure with no negative frequencies
                continue
            else:
                converged = True
                return_code = subprocess.call("formchk " + self.freq_chk, shell=True)
                if return_code != 0:
                    raise RuntimeError("formatted checkpoint file for frequency calculation could not be created")

            final_structure = gout.final_structure.as_dict()['sites']
            final_energy = gout.final_energy
            gibbs_correction = gout.corrections.get("Gibbs Free Energy")

        if not converged:
            raise RuntimeError("Structure not converged. calc_dir " + self.working_dir)

        # Clean up files and transfer to website processing
        if not self.skip_freq:
            self.post_job(upload_files=[self.freq_log, self.freq_fchk], calc_type=freq_name,
                          delete_files=[self.freq_log, self.freq_chk, self.freq_fchk, self.freq_com], name_tag="freq_")
        self.post_job()

        return FWAction(
            update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier, "gs_charge": self.gs_charge,
                         "gs_spin": self.gs_spin,
                         "{}_hash".format(self.full_name): opt_hash,
                         "{}_hash".format(freq_name): freq_hash,
                         "{}_geom".format(self.full_name): final_structure,
                         "{}_eng".format(self.full_name): final_energy,
                         "{}_gibb".format(self.full_name): gibbs_correction,
                         "iop_str": self.get('iop_str')})


@explicit_serialize
class RunWtuning(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameters
        data_store = get_data_store(self, fw_spec)
        path = env_chk(self.get('path'), fw_spec)
        runfile_log = env_chk(self.get('runfile_log'), fw_spec)
        paramset = self["paramset"]
        identifier = fw_spec.get("identifier", ) or self.get("identifier")
        smiles = fw_spec.get("smiles", ) or self.get("smiles")
        gs_charge = fw_spec.get("gs_charge") or self.get("gs_charge") or get_groundState(
            identifier, smiles, data_store=data_store)
        gs_spin = fw_spec.get("gs_spin") or self.get("gs_spin") or get_groundState(
            identifier, smiles, prop='spin', data_store=data_store)
        gaussian_file_name = fw_spec.get("gaussian_file_name") or self.get("gaussian_file_name") or "gaussian"
        calc_dir = "{}/{}/{}".format(path, identifier or smiles, gaussian_file_name)
        store_dir = get_store_dir(identifier or smiles, get_local_output_dir(self, fw_spec), calc_dir)
        check_if_already_run = spec_or_task(self, fw_spec, "check_if_already_run", True)
        submit = spec_or_task(self, fw_spec, "submit", True)
        restricted = spec_or_task(self, fw_spec, "restricted", True)

        radical_electrons = (gs_spin - 1 - paramset.charge) % 2
        paramset.multiplicity = radical_electrons + 1  # calculate spin multiplicity with Hand's rule
        paramset.charge += gs_charge
        geometry = fw_spec.get("geometry", ) or self.get("geometry", )
        geometry_sites = fw_spec.get("{}_geom".format(geometry), )
        if geometry_sites:
            mol = Molecule.from_sites([Site.from_dict(sd) for sd in geometry_sites])
        else:
            geometry_hash = fw_spec.get("{}_hash".format(geometry), )
            mol = get_db_geom(geometry_hash, data_store=data_store, store_dir=store_dir) \
                or start_from_smiles(identifier, smiles, data_store=data_store)

        # create folder
        working_dir = "{}/{}".format(calc_dir, "wtuning")
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)
        gauss_inp = generate_gaussian_input(paramset=paramset, mol=mol)
        gauss_inp.write_file('wtuning.com', cart_coords=True)

        # check if job has already run
        if check_if_already_run:
            tuned_w = get_omega(identifier, data_store=data_store, store_dir=store_dir)
            if tuned_w:
                iop_str = omega_to_iop(tuned_w)
                return FWAction(
                    update_spec={"iop_str": iop_str, "gaussrun_dir": calc_dir, "identifier": identifier,
                                 "gs_charge": gs_charge, "gs_spin": gs_spin})

        # generate the input for wtuning in gaussian and run tuning
        functional = paramset.functional if restricted else "R{}".format(paramset.functional)
        wt_mol = WtuningJob(func=functional, basis=paramset.basis_set, name="wtuning", nproc=nprocs, mem=48,
                            n_charge=paramset.charge, n_spin=paramset.multiplicity, wdir='./', scheme='Jh', wbmin=0.05,
                            wbmax=0.5)
        wt_mol.mol = mol
        wt_mol.wtuning_cycle(max_cycles=0)

        # get the tuned w value from the run
        os.chdir(working_dir)
        with open("output.log", "r") as f:
            output = f.read()
        tuned_w = float(output.strip().split()[-3])
        iop_str = omega_to_iop(tuned_w)

        # Store data (D3TaLES: upload through website; local: write JSON record)
        if data_store == "d3tales":
            file_list = ['wtuning.com', 'output.log']
            zip_name = zip_files(file_list, zip_name='{}/{}_wtuning.zip'.format(working_dir, identifier))
            if submit:
                d3t_submit_calculation(zip_name, identifier, self['name'])
        else:
            conditions = orig_conditions(paramset.functional, paramset.basis_set)  # no tuning_parameter, as in back2front
            record = local_record(WTUNING_KEY, identifier, self['name'],
                                  {"omega": tuned_w, "conditions": conditions,
                                   "output_file": os.path.join(working_dir, "output.log")})
            print("Record saved to", save_local_result(store_dir, record))
        data_hash = get_hash_id(identifier, 'wtuning.com', self['name'], output_file='output.log')
        with open(runfile_log, 'a') as fn:
            fn.write("{}, {}\n".format(data_hash, working_dir))
        os.system("rm -rf *.chk *.fchk tuning_wtuning_ocycle*")
        return FWAction(
            update_spec={"iop_str": iop_str, "gaussrun_dir": calc_dir, "identifier": identifier, "gs_charge": gs_charge,
                         "gs_spin": gs_spin})


@explicit_serialize
class RunGaussianTDDFT(GaussianBase):

    # Does not yet have solvent capability
    def run_task(self, fw_spec):
        setup_obj = self.setup_calc(fw_spec, calc_type='energy')
        if isinstance(setup_obj, FWAction):
            return setup_obj

        # generate the input for gaussian energy run
        gauss_inp = generate_gaussian_input(paramset=self.paramset, mol=self.mol)
        gauss_inp.write_file(self.file_com, cart_coords=True)
        if self.check_if_already_run:
            existing_data = self.existing_data()
            if existing_data:
                return FWAction(
                    update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier, "iop_str": self.iop_str,
                                 "gs_charge": self.gs_charge, "gs_spin": self.gs_spin})

        # run gaussian16
        print('SUBMITTING GAUSSIAN TDDFT JOB {} for {}'.format(self.full_name, self.identifier))
        subprocess.call(self.gaussian_cmd + " " + self.file_com, shell=True)

        # check for normal termination of gaussian job
        gout = GaussianOutput(self.file_log)
        if not gout.properly_terminated:
            raise RuntimeError(self.file_com + " not terminated. calc_dir " + self.working_dir)
        else:
            return_code = subprocess.call("formchk " + self.file_chk, shell=True)

            if return_code != 0:
                raise RuntimeError("TDDFT calculation formatted checkpoint file could not be created")

            if self.run_nto:
                nto_com, nto_chk = write_nto(self.file_log, 1, tddft_chk=self.file_chk)
                subprocess.call(self.gaussian_cmd + " " + nto_com, shell=True)
                nto_gout = GaussianOutput(nto_com.replace(".com", ".log"))
                if not nto_gout.properly_terminated:
                    raise RuntimeError(nto_com + " not terminated. calc_dir " + self.working_dir)
                return_code = subprocess.call("formchk " + nto_chk, shell=True)
                if return_code != 0:
                    raise RuntimeError("NTO analysis formatted checkpoint file could not be created")

            self.post_job()
        return FWAction(
            update_spec={"gaussrun_dir": self.calc_dir, "identifier": self.identifier, "iop_str": self.iop_str,
                         "gs_charge": self.gs_charge, "gs_spin": self.gs_spin})
