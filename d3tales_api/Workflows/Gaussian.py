import abc
import os
import subprocess
import multiprocessing
import re
import uuid
from copy import deepcopy
from six import add_metaclass
from atomate.utils.utils import get_logger
from fireworks import FiretaskBase, explicit_serialize, FWAction
from d3tales_api.Workflows.utils import (
    ElectronicStructureTaskMixin,
    GaussianInput,
    GaussianOutput,
    Molecule,
    Site,
    WTUNING_KEY,
    generate_gaussian_input,
    get_db_geom,
    get_groundState,
    get_hash_id,
    get_local_output_dir,
    get_omega,
    get_store_dir,
    geometry_key_candidates,
    local_record,
    omega_to_iop,
    orig_conditions,
    orig_hash_id,
    save_local_result,
    resolve_charge_multiplicity,
    spec_or_task,
    start_from_smiles,
    forcefield_conformer_for_task,
    write_nto,
    zip_files
)
from d3tales_api.Workflows.wtuning import WtuningJob
from d3tales_api.D3database.workflow_adapter import (
    d3t_submit_calculation,
    d3t_direct_store_calculation,
)
from d3tales_api.Workflows.runtime import WorkerRuntime

logger = get_logger(__name__)
cpus = [multiprocessing.cpu_count() if multiprocessing.cpu_count() < 16 else 16]
nprocs = str(cpus[0])


# Copyright 2021, University of Kentucky


@add_metaclass(abc.ABCMeta)
class GaussianBase(ElectronicStructureTaskMixin, FiretaskBase):
    engine_name = "Gaussian"
    _fw_name = "GaussianBase"

    @staticmethod
    def _link0_value(parameters, name):
        normalized = name.lower().lstrip("%")
        for key, value in parameters.items():
            if str(key).lower().lstrip("%") == normalized:
                return value
        return None

    @staticmethod
    def _expand_resource_indices(value, name):
        if not isinstance(value, str) or not value.strip():
            raise ValueError("{} must be a non-empty comma-separated list of indices or ranges".format(name))
        indices = []
        for part in value.split(","):
            match = re.fullmatch(r"\s*(\d+)(?:-(\d+))?\s*", part)
            if not match:
                raise ValueError("{} contains an invalid index or range: {!r}".format(name, part.strip()))
            start = int(match.group(1))
            end = int(match.group(2) or start)
            if end < start:
                raise ValueError("{} ranges must be ascending".format(name))
            indices.extend(range(start, end + 1))
        if len(indices) != len(set(indices)):
            raise ValueError("{} cannot contain duplicate indices".format(name))
        return indices

    def configure_link0_resources(self, fw_spec):
        """Set Gaussian CPU/GPU resource directives from FireWorker configuration."""
        if not hasattr(self, "_base_link0_parameters"):
            self._base_link0_parameters = dict(self.paramset.link0_parameters)
        parameters = dict(self._base_link0_parameters)
        runtime = getattr(self, "worker_runtime", None)
        if runtime is None:
            runtime = WorkerRuntime.from_task(self, fw_spec, software="gaussian")
        memory = runtime.engine_options.get("memory")
        nprocshared = runtime.engine_options.get("nprocshared")
        cpu_list = runtime.engine_options.get("cpu_list")
        gpu_cpu = runtime.engine_options.get("gpu_cpu")

        if memory is None:
            memory = self._link0_value(parameters, "mem") or "48GB"
        if nprocshared is None:
            nprocshared = self._link0_value(parameters, "nprocshared") or nprocs
        if not isinstance(memory, str) or not memory.strip() or "\n" in memory:
            raise ValueError("memory must be a non-empty single-line Gaussian memory value")

        if gpu_cpu is not None:
            if cpu_list is None:
                raise ValueError("cpu_list must be set when gpu_cpu is configured")
            allocated_cpus = self._expand_resource_indices(cpu_list, "cpu_list")
            if not isinstance(gpu_cpu, str) or gpu_cpu.count("=") != 1:
                raise ValueError("gpu_cpu must use Gaussian's gpu-list=control-cpus format")
            gpu_ids, control_cpus = gpu_cpu.split("=")
            gpu_ids = self._expand_resource_indices(gpu_ids, "gpu_cpu GPU list")
            control_cpus = self._expand_resource_indices(control_cpus, "gpu_cpu control CPU list")
            if len(gpu_ids) != len(control_cpus):
                raise ValueError("gpu_cpu must specify exactly one control CPU per GPU")
            if not set(control_cpus).issubset(allocated_cpus):
                raise ValueError("cpu_list must include every GPU control CPU")

        resource_names = {"mem", "nprocshared", "cpu", "gpucpu"}
        parameters = {
            key: value for key, value in parameters.items()
            if str(key).lower().lstrip("%") not in resource_names
        }
        parameters["%mem"] = memory.strip()
        if gpu_cpu is not None:
            parameters["%CPU"] = cpu_list.strip()
            parameters["%GPUCPU"] = gpu_cpu.strip()
        elif cpu_list is not None:
            self._expand_resource_indices(cpu_list, "cpu_list")
            parameters["%CPU"] = cpu_list.strip()
        else:
            try:
                nprocshared = int(nprocshared)
            except (TypeError, ValueError):
                raise ValueError("nprocshared must be a positive integer")
            if nprocshared < 1:
                raise ValueError("nprocshared must be a positive integer")
            parameters["%nprocshared"] = str(nprocshared)
        self.paramset.link0_parameters = parameters

    def setup_files(self, fw_spec, calc_type='opt'):
        # get parameters
        self.worker_runtime = WorkerRuntime.from_task(
            self, fw_spec, software="gaussian", require_executables=True
        )
        self.data_store = self.worker_runtime.data_store
        self.gaussian_cmd = self.worker_runtime.engine_options.get("g16_cmd")
        name_tag = spec_or_task(self, fw_spec, "name_tag", "")
        self.full_name = name_tag + self['name']
        self.calc_name = self['name']
        self.paramset = deepcopy(self["paramset"])
        self.identifier = spec_or_task(self, fw_spec, "identifier", uuid.uuid4().__str__())
        self.smiles = spec_or_task(self, fw_spec, "smiles")
        self.runfile_log = self.worker_runtime.runfile_log
        self.skip_freq = spec_or_task(self, fw_spec, "skip_freq", False)
        self.skip_freq = True if self.calc_name == "opt_mol" else self.skip_freq
        self.check_if_already_run = spec_or_task(self, fw_spec, "check_if_already_run", True)
        self.restricted = spec_or_task(self, fw_spec, "restricted", True)
        self.run_nto = spec_or_task(self, fw_spec, "run_nto", False)

        self.gaussian_file_name = spec_or_task(self, fw_spec, "gaussian_file_name", "gaussian")

        # set prefix for the file names
        if self.get("prefix"):
            prefix = name_tag + self["prefix"]
        elif self.get("subtype", ):
            prefix = name_tag + calc_type + "_" + self["subtype"]
        else:
            prefix = name_tag + calc_type + "_mol"

        # get calc directory
        if self.worker_runtime.path:
            self.calc_dir = "{}/{}/{}".format(self.worker_runtime.path, self.identifier,
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

    def setup_calc(self, fw_spec, calc_type='opt', geom_from_smiles=False):
        self.setup_files(fw_spec, calc_type=calc_type)
        name_tag = spec_or_task(self, fw_spec, "name_tag", "")
        self.solvent = spec_or_task(self, fw_spec, "solvent")
        self.gs_charge = spec_or_task(self, fw_spec, "gs_charge")
        if self.gs_charge is None:
            self.gs_charge = get_groundState(
                self.identifier, self.smiles, data_store=self.data_store
            )
        self.gs_spin = spec_or_task(self, fw_spec, "gs_spin")
        if self.gs_spin is None:
            self.gs_spin = get_groundState(
                self.identifier, self.smiles, prop='spin', data_store=self.data_store
            )
        use_iop = spec_or_task(self, fw_spec, "use_iop", True)
        run_from_com = spec_or_task(self, fw_spec, "run_from_com", False)

        # get iop string (omega value) for the calculation
        self.iop_str = spec_or_task(self, fw_spec, "iop_str")
        if not self.iop_str and use_iop:
            tuned_w = get_omega(
                self.identifier, data_store=self.data_store, store_dir=self.store_dir, engine="Gaussian"
            )
            if tuned_w:
                self.iop_str = omega_to_iop(tuned_w)
        self.result_iop_str = self.iop_str  # fixed copy used for local record hashes

        # get starting geometry for calculation
        geometry = spec_or_task(self, fw_spec, "geometry")
        geometry_key = next(
            (key for key in geometry_key_candidates(self, fw_spec, geometry)
             if fw_spec.get("{}_geom".format(key)) is not None),
            geometry,
        )
        geometry_sites = fw_spec.get("{}_geom".format(geometry_key))
        forcefield_mol = forcefield_conformer_for_task(self, fw_spec)
        if forcefield_mol is not None:
            self.mol = forcefield_mol
        elif run_from_com:
            self.mol = GaussianInput.from_file(self.file_com).molecule
        else:
            if geometry_sites:
                self.mol = Molecule.from_sites([Site.from_dict(sd) for sd in geometry_sites])
            else:
                geometry_hash = fw_spec.get("{}_hash".format(geometry_key),
                                            orig_hash_id(self.identifier, geometry, self.paramset.functional,
                                                         self.paramset.basis_set, tuning_parameter=self.iop_str,
                                                         solvent=self.solvent))
                hash_geom = get_db_geom(geometry_hash, data_store=self.data_store, store_dir=self.store_dir)
                if hash_geom:
                    self.mol = hash_geom
                elif geom_from_smiles:
                    self.mol = start_from_smiles(self.identifier, self.smiles, data_store=self.data_store)
                else:
                    raise Exception(f"No geometry found for {self.identifier} and geom_from_smiles is set to False.")

            # End job if the total number of atoms is greater than 200
            num_atoms = len(self.mol.sites)
            if num_atoms > 200:
                return FWAction(defuse_children=True,
                                update_spec={"defuser_reason": "The molecule has {} atoms".format(num_atoms)})

        # Update route parameters with omega value and pcm solvent
        self.paramset.charge, self.paramset.multiplicity = resolve_charge_multiplicity(
            self, self.paramset, fw_spec, self.gs_charge, self.gs_spin
        )
        self.functional = self.paramset.functional if self.restricted else "U{}".format(self.paramset.functional)
        self.paramset.functional = self.functional
        print("Charge: ", self.paramset.charge, "\t Multiplicity: ", self.paramset.multiplicity)
        self.configure_link0_resources(fw_spec)
        self.paramset.link0_parameters["%chk"] = self.file_chk
        if self.iop_str and use_iop:
            self.paramset.route_parameters.update({"iop(3/107={}, 3/108={})".format(self.iop_str, self.iop_str): ""})
        if self.solvent:
            self.paramset.route_parameters.update({"SCRF": "(PCM,Solvent={})".format(self.solvent)})




    def post_job(self, upload_files=None, delete_files=None, calc_type=None, name_tag=''):
        # Store data (D3TaLES: upload through website; local: write JSON record)
        upload_files = upload_files or [self.file_log, self.file_fchk]
        delete_files = delete_files or [self.file_chk, self.file_com]
        calc_type = calc_type or self.calc_name
        _hash = get_hash_id(self.identifier, self.file_log, self.calc_name) if os.path.isfile(self.file_log) else None
        if self.data_store == "d3tales_web":
            upload_names = [f.split('/')[-1] for f in upload_files]
            zip_path = zip_files(upload_names, zip_name='{}_{}.zip'.format(self.identifier, name_tag + self.full_name))
            d3t_submit_calculation(zip_path, self.identifier, calc_type)
            print("File {}_{}.zip successfully posted!".format(self.identifier, name_tag + self.full_name))
        elif self.data_store == "d3tales_direct":
            from d3tales_api.Processors.parser_dft import ProcessGausLog
            logfile = next((f for f in upload_files if f.endswith(".log")), self.file_log)
            record = ProcessGausLog(_id=self.identifier, metadata={
                "mol_file": logfile, "calculation_type": calc_type,
                "mol_info": {"groundState_charge": self.gs_charge}}).data_dict
            archive = os.path.join(self.store_dir, "{}_{}.zip".format(self.identifier, record["_id"]))
            d3t_direct_store_calculation(record, files=upload_files, archive_path=archive)
        else:
            from d3tales_api.Processors.parser_dft import ProcessGausLog
            logfile = next((f for f in upload_files if f.endswith(".log")), self.file_log)
            record = ProcessGausLog(_id=self.identifier, metadata={
                "mol_file": logfile, "calculation_type": calc_type,
                "mol_info": {"groundState_charge": self.gs_charge}}).data_dict
            print("Record saved to", save_local_result(self.store_dir, record, key=self.local_key(calc_type)))
        # Write runfile to runfile_log so the runfile can be deleted after calculation
        if self.runfile_log:
            with open(self.runfile_log, 'a') as fn:
                fn.write("{}, {}\n".format(_hash, self.working_dir))
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
                return FWAction(update_spec=self.common_update_spec(
                    gaussrun_dir=self.calc_dir,
                    **{
                        "{}_eng".format(self.full_name): existing_data.get("scf_total_energy", {}).get("value"),
                        "{}_hash".format(self.full_name): get_hash_id(
                            self.identifier, self.file_com, self.calc_name
                        ),
                    }
                ))

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
        return FWAction(update_spec=self.common_update_spec(
            gaussrun_dir=self.working_dir,
            **{"{}_eng".format(self.full_name): final_energy}
        ))


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
                    return FWAction(update_spec=self.common_update_spec(
                        gaussrun_dir=self.calc_dir,
                        **{
                            "{}_hash".format(self.full_name): opt_hash,
                            "{}_hash".format(freq_name): freq_hash,
                            "{}_geom".format(self.full_name): existing_data.get("geometry"),
                            "{}_eng".format(self.full_name): existing_data.get("scf_total_energy", {}).get("value"),
                            "{}_gibb".format(self.full_name): self.existing_data(
                                _hash=freq_hash, calc_name=freq_name
                            ).get("gibbs_correction", {}).get("value"),
                        }
                    ))

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

        return FWAction(update_spec=self.common_update_spec(
            gaussrun_dir=self.calc_dir,
            **{
                "{}_hash".format(self.full_name): opt_hash,
                "{}_hash".format(freq_name): freq_hash,
                "{}_geom".format(self.full_name): final_structure,
                "{}_eng".format(self.full_name): final_energy,
                "{}_gibb".format(self.full_name): gibbs_correction,
            }
        ))


@explicit_serialize
class RunWtuning(FiretaskBase):

    def run_task(self, fw_spec):
        # get parameters
        worker_runtime = WorkerRuntime.from_task(self, fw_spec, software="gaussian")
        data_store = worker_runtime.data_store
        path = worker_runtime.path
        runfile_log = worker_runtime.runfile_log
        gaussian_cmd = worker_runtime.engine_options.get("g16_cmd")
        if not path:
            raise ValueError("A worker path is required for Gaussian w-tuning")
        paramset = deepcopy(self["paramset"])
        identifier = spec_or_task(self, fw_spec, "identifier")
        smiles = spec_or_task(self, fw_spec, "smiles")
        gs_charge = spec_or_task(self, fw_spec, "gs_charge")
        if gs_charge is None:
            gs_charge = get_groundState(identifier, smiles, data_store=data_store)
        gs_spin = spec_or_task(self, fw_spec, "gs_spin")
        if gs_spin is None:
            gs_spin = get_groundState(identifier, smiles, prop='spin', data_store=data_store)
        gaussian_file_name = spec_or_task(self, fw_spec, "gaussian_file_name", "gaussian")
        calc_dir = "{}/{}/{}".format(path, identifier or smiles, gaussian_file_name)
        store_dir = get_store_dir(identifier or smiles, get_local_output_dir(self, fw_spec), calc_dir)
        check_if_already_run = spec_or_task(self, fw_spec, "check_if_already_run", True)
        restricted = spec_or_task(self, fw_spec, "restricted", True)

        paramset.charge, paramset.multiplicity = resolve_charge_multiplicity(
            self, paramset, fw_spec, gs_charge, gs_spin
        )
        geometry = spec_or_task(self, fw_spec, "geometry")
        geometry_key = next(
            (key for key in geometry_key_candidates(self, fw_spec, geometry)
             if fw_spec.get("{}_geom".format(key)) is not None),
            geometry,
        )
        geometry_sites = fw_spec.get("{}_geom".format(geometry_key))
        if geometry_sites:
            mol = Molecule.from_sites([Site.from_dict(sd) for sd in geometry_sites])
        else:
            geometry_hash = fw_spec.get("{}_hash".format(geometry_key), )
            mol = get_db_geom(geometry_hash, data_store=data_store, store_dir=store_dir) \
                or (start_from_smiles(identifier, smiles, data_store=data_store)
                    if self.get("geom_from_smiles", False) else None)
            if mol is None:
                raise LookupError("No geometry found for requested geometry {!r}".format(geometry))

        # FireWorks only propagates update_spec to direct children. Forward the molecule inputs so the
        # optimizations downstream of w-tuning see the same geometry/SMILES as they would with wtune=False.
        forward_spec = {"data_store": data_store}
        if smiles:
            forward_spec["smiles"] = smiles
        if geometry and geometry_sites:
            forward_spec["{}_geom".format(geometry_key)] = geometry_sites
        if geometry_key and fw_spec.get("{}_hash".format(geometry_key)) is not None:
            forward_spec["{}_hash".format(geometry_key)] = fw_spec["{}_hash".format(geometry_key)]

        # create folder
        working_dir = "{}/{}".format(calc_dir, "wtuning")
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)
        gauss_inp = generate_gaussian_input(paramset=paramset, mol=mol)
        gauss_inp.write_file('wtuning.com', cart_coords=True)

        # check if job has already run
        if check_if_already_run:
            tuned_w = get_omega(identifier, data_store=data_store, store_dir=store_dir, engine="Gaussian")
            if tuned_w:
                iop_str = omega_to_iop(tuned_w)
                return FWAction(
                    update_spec={"iop_str": iop_str, "gaussrun_dir": calc_dir, "identifier": identifier,
                                 "gs_charge": gs_charge, "gs_spin": gs_spin, **forward_spec})

        # generate the input for wtuning in gaussian and run tuning
        if not gaussian_cmd:
            raise KeyError("FireWorker environment variable 'g16_cmd' is required to run Gaussian w-tuning")
        functional = paramset.functional if restricted else "U{}".format(paramset.functional)
        wt_mol = WtuningJob(
            func=functional,
            basis=paramset.basis_set,
            name="wtuning",
            nproc=nprocs,
            mem=48,
            n_charge=paramset.charge,
            n_spin=paramset.multiplicity,
            wdir='./',
            scheme='Jh',
            wbmin=0.05,
            wbmax=0.5,
            gauss_cmd=gaussian_cmd,
        )
        wt_mol.mol = mol
        wt_mol.wtuning_cycle(max_cycles=0)

        # get the tuned w value from the run
        os.chdir(working_dir)
        with open("output.log", "r") as f:
            output = f.read()
        tuned_w = float(output.strip().split()[-3])
        iop_str = omega_to_iop(tuned_w)

        # Store data through the website, directly in MongoDB, or in a local JSON record.
        if data_store == "d3tales_web":
            file_list = ['wtuning.com', 'output.log']
            zip_name = zip_files(file_list, zip_name='{}/{}_wtuning.zip'.format(working_dir, identifier))
            d3t_submit_calculation(zip_name, identifier, self['name'])
        elif data_store == "d3tales_direct":
            conditions = orig_conditions(paramset.functional, paramset.basis_set)
            record = local_record(
                orig_hash_id(identifier, self['name'], paramset.functional, paramset.basis_set),
                identifier,
                self['name'],
                {"engine": "Gaussian", "omega": tuned_w, "conditions": conditions,
                 "output_file": os.path.join(working_dir, "output.log")},
            )
            d3t_direct_store_calculation(
                record,
                files=[os.path.join(working_dir, name) for name in ("wtuning.com", "output.log")],
                archive_path=os.path.join(store_dir, "{}_{}.zip".format(identifier, record["_id"])),
            )
        else:
            conditions = orig_conditions(paramset.functional, paramset.basis_set)  # no tuning_parameter, as in back2front
            record = local_record(WTUNING_KEY, identifier, self['name'],
                                  {"omega": tuned_w, "conditions": conditions,
                                   "output_file": os.path.join(working_dir, "output.log")})
            print("Record saved to", save_local_result(store_dir, record))
        data_hash = get_hash_id(identifier, 'wtuning.com', self['name'], output_file='output.log')
        if runfile_log:
            with open(runfile_log, 'a') as fn:
                fn.write("{}, {}\n".format(data_hash, working_dir))
        os.system("rm -rf *.chk *.fchk tuning_wtuning_ocycle*")
        return FWAction(
            update_spec={"iop_str": iop_str, "gaussrun_dir": calc_dir, "identifier": identifier, "gs_charge": gs_charge,
                         "gs_spin": gs_spin, **forward_spec})


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
                return FWAction(update_spec=self.common_update_spec(
                    gaussrun_dir=self.calc_dir
                ))

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
        return FWAction(update_spec=self.common_update_spec(
            gaussrun_dir=self.calc_dir
        ))
