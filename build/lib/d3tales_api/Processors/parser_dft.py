import re
from d3tales_api.Calculators.plotters import *
from d3tales_api.Calculators.calculators import *
from pymatgen.io.gaussian import GaussianOutput, GaussianInput


class ProcessCCLIB:
    """
    Class to process most logfiles using https://cclib.github.io/index.html.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, metadata=None):
        """
        :param metadata: dict, metadata (`mol_file` should be a key to the filepath of the file to be processed)
        """
        self.log_path = metadata.get('mol_file', )
        self.wtuning_output = metadata.get('wtuning_output', '')
        self.calculation_type = metadata.get('calculation_type', )

    def cclib_parse(self):
        """
        Use CCLIB to parse data file
        """
        if not self.wtuning_output:
            self.cmol = cclib.io.ccopen(self.log_path).parse()
            # self.functional = self.cmol.functional  # TODO basis set
            # self.basis_set = self.cmol.basis_set  # TODO basis set
            self.charge = self.cmol.charge
            # self.spin_multiplicity = self.cmol.atomspin  # TODO spin

            # Skip calculation attributes if this is a Gaussian Input file
            file_extension = self.log_path.split('.')[-1]
            # self.electrons = getattr(self.cmol, "electrons", None)  # TODO electrons
            self.final_energy = self.cmol.scfenergies[-1]  # already in eV
            # self.final_structure = self.get_sites_dict(self.cmol) # TODO fix
            self.homo = self.cmol.moenergies[0][self.cmol.homos[0]]
            self.lumo = self.cmol.moenergies[0][int(self.cmol.homos[0]) + 1]
            self.frequency_dicts = self.get_freq_dicts(self.cmol)
            # self.gibbs_correction = self.cmol.corrections.get("Gibbs Free Energy") # TODO gibbs_correction

    # TODO runtime
    # TODO TDDFT
    # TODO dipole
    # TODO solvents stuff
    # TODO omega
    @staticmethod
    def get_sites_dict(cmol):
        """
        Get a pymatgen sites dictionary list from a cclib molecule object
        :return: list of atom site dicts
        """
        coords = cmol.atomcoords[-1]
        atomic_nums = cmol.atomnos
        atoms = [periodictable[a] for a in atomic_nums]
        return [{'name': atoms[i], 'species': [{'element': atoms[i], 'occu': 1}], 'xyz': c, 'properties': {}} for i, c
                in enumerate(coords)]

    @staticmethod
    def get_freq_dicts(cmol):
        """
        Get a pymatgen frequency dictionary list from a cclib molecule object
        :return: list of frequency dicts
        """
        try:
            freqs = cmol.vibfreqs
            masses = cmol.vibrmasses
            vibfconsts = cmol.vibfconsts
            vibirs = cmol.vibirs
            vibsyms = cmol.vibsyms
            return [{'frequency': f, 'r_mass': masses[i], 'f_constant': vibfconsts[i], 'IR_intensity': vibirs[i],
                     'symmetry': vibsyms[i]} for i, f in enumerate(freqs)]
        except:
            return []

    def get_radical_stability_score(self, spin_type="mulliken"):
        """
        Get radical stability score
        :param spin_type: str, type of spin to use
        :return: dictionary of radical stability score and its components
        """
        c_data = {"log_file": self.log_path, "spin_type": spin_type}
        connector = {"log_file": "log_file", "spin_type": "spin_type"}
        return dict(radical_buried_vol={"value": RadBuriedVolCalc(connector=connector).calculate(c_data), "unit": "A^3"},
                    radical_spin={"value": RadicalSpinCalc(connector=connector).calculate(c_data)},
                    radical_stability_score={"value": RSSCalc(connector=connector).calculate(c_data)})


class ProcessPsi4Log(ProcessCCLIB):
    """
    Class to process Psi4 logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, metadata=None):
        super().__init__(metadata=metadata)

        # TODO finish


class ProcessGausLog(ProcessCCLIB):
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, metadata=None):
        """
        :param metadata: dict, metadata (`mol_file` should be a key to the filepath of the file to be processed)
        """
        super().__init__(metadata=metadata)

        self.log_path = metadata.get('mol_file', )
        self.wtuning_output = metadata.get('wtuning_output', '')
        self.calculation_type = metadata.get('calculation_type', )
        self.code_name = metadata.get('code_name', 'Gaussian')
        self.code_version = metadata.get('code_version', '16')

        self.pmgmol = self.get_pmgmol(self.log_path)
        self.functional = self.pmgmol.functional
        self.basis_set = self.pmgmol.basis_set
        self.charge = self.pmgmol.charge
        self.spin_multiplicity = self.pmgmol.spin_multiplicity

        # Skip calculation attributes if this is a Gaussian Input file
        file_extension = self.log_path.split('.')[-1]
        if file_extension != 'gjf' and file_extension != 'com':
            self.electrons = getattr(self.pmgmol, "electrons", None)
            self.final_energy = self.pmgmol.final_energy * 27.2114  # convert to eV
            self.final_structure = self.pmgmol.final_structure.as_dict()['sites']
            self.homo_1 = self.homo_lumo[0]
            self.homo = self.homo_lumo[1]
            self.lumo = self.homo_lumo[2]
            self.lumo_1 = self.homo_lumo[3]
            self.frequency_dicts = [f for f in self.pmgmol.frequencies[0]] if getattr(self.pmgmol, "frequencies",
                                                                                      None) else {}
            self.gibbs_correction = self.pmgmol.corrections.get("Gibbs Free Energy")

            self.get_solvent_info()
            self.parse_tddft()

    @staticmethod
    def get_pmgmol(file):
        """
        Get Pymatgen molecule object from file

        :param file: str, file path
        :return: Pymatgen molecule object
        """
        file_extension = file.split('.')[-1]
        if file_extension == 'gjf' or file_extension == 'com':
            mol = GaussianInput.from_file(file)
        elif file_extension == 'log':
            mol = GaussianOutput(file)
        else:
            mol = Molecule.from_file(file)

        return mol

    @property
    def runtime(self):
        """Runtime in core hours from a logfile"""
        time_patt = re.compile(r"\d+\.d+|\d+")
        time_data = []
        with open(self.log_path, "r") as f:
            line = f.readline()
            while line != "":
                if re.match("Job cpu time", line.strip()):
                    time_data.extend(time_patt.findall(line))
                line = f.readline()
        if time_data:
            time_data = [float(time) for time in time_data]
            runtime = (time_data[0] * 86400 + time_data[1] * 3600 + time_data[2] * 60 + time_data[3]) / 3600
        else:
            runtime = 0
        return round(runtime, 3)

    @staticmethod
    def get_tddft_excitations(log_path):
        """
        Read excitation energies after a TD-DFT calculation.

        :param log_path: str, filepath to log file
        :return: A list of tuple for each transition such as [(energy (eV), lambda (nm), oscillatory strength), ... ]
        """
        float_patt = re.compile(r"\s*([+-]?\d+\.\d+)")
        state_patt = re.compile(r"[a-zA-Z]*let")
        transitions = {"Singlet": [],
                       "Doublet": [],
                       "Triplet": []}

        # read in file
        with open(log_path, "r") as f:
            line = f.readline()
            td = False
            while line != "":
                if re.search(r"^\sExcitation energies and oscillator strengths:", line):
                    td = True

                if td:
                    if re.search(r"^\sExcited State\s*\d", line):
                        val = [float(v) for v in float_patt.findall(line)]
                        try:
                            state = state_patt.findall(line)[0]
                        except Exception:
                            state_val = val.pop(0)
                            if round(state_val) == 1:
                                state = "Singlet"
                            elif round(state_val) == 2:
                                state = "Doublet"
                            elif round(state_val) == 3:
                                state = "Triplet"
                            else:
                                raise ValueError(
                                    "Calculation has a spin state greater than triplet -- spin {}".format(state_val))
                        transitions[state].append(tuple(val))
                line = f.readline()
        return transitions

    @property
    def dipole_moments(self):
        """A list of tuples, each containing teh (X, Y, Z) coordinate for a dipole moment such as [(X1, Y1, Z1), (X2, Y2, Z2)... ]"""

        coord_patt = re.compile(r"[A-Z]+=\s+[+-]?[0-9]*\.[0-9]+")
        dipole_moments = []

        # read in file
        with open(self.log_path, "r") as f:
            line = f.readline()
            td = False
            while line != "":
                if re.search(r"^\sDipole moment \(field-independent basis, Debye\):", line):
                    td = True
                    line = f.readline()
                if td:
                    coord_dict = {}
                    coordinates = [c for c in coord_patt.findall(line)]
                    for coord in coordinates:
                        dimension = re.compile(r"[A-Z]").findall(coord)[0]
                        value = float(re.compile(r"[+-]?[0-9]*\.[0-9]+").findall(coord)[0])
                        coord_dict[dimension] = value
                    coord_list = [coord_dict["X"], coord_dict["Y"], coord_dict["Z"]]
                    dipole_moments.append(tuple(coord_list))
                    td = False
                line = f.readline()
        return dipole_moments

    @property
    def pcm_info(self):
        """Solvent and dielectric constant from Gaussian logfile:  [solvent_name, dielectric_constant] - list containing solvent name and dielectric constant]"""
        name_patt = re.compile(r":(.*?),\s*Eps")
        dc_patt = re.compile(r"Eps=\s*(.*?)\s*Eps")
        with open(self.log_path, "r") as f:
            line = f.readline()
            pcm = False
            while line != "":
                if re.search(r"^\sPolarizable Continuum Model", line):
                    pcm = True
                if pcm:
                    if re.search(r"^\sSolvent", line):
                        solvent_name = name_patt.findall(line)[0].strip()
                        dielectric_constant = float(dc_patt.findall(line)[0].strip())
                        return [solvent_name, dielectric_constant]
                line = f.readline()

    @property
    def homo_lumo(self):
        """HOMO and LUMO energies from a Gaussian molecule: [homo, lumo] - list containing homo then lumo in eV"""
        num_electrons = self.pmgmol.electrons[0]
        eigens = list(self.pmgmol.eigenvalues.values())[0]
        homo_1 = eigens[num_electrons - 2] * 27.2114  # convert to eV
        homo = eigens[num_electrons - 1] * 27.2114  # convert to eV
        lumo = eigens[num_electrons] * 27.2114  # convert to eV
        lumo_1 = eigens[num_electrons + 1] * 27.2114  # convert to eV

        return [homo_1, homo, lumo, lumo_1]

    @property
    def omega(self):
        """Omega value from w-tuning output file (generated with the OCELOT w-tuning module)"""
        try:
            with open(self.wtuning_output, 'r') as fn:
                w_data = fn.readlines()[-2].split()[1]
            return float("0.{}".format(w_data.split('.')[1]))
        except FileNotFoundError:
            return None

    @property
    def tuning_parameter(self):
        """Omega value from Pymatgen molecule object"""
        route_params = self.pmgmol.route_parameters
        route_params = {k.lower(): v for k, v in route_params.items()}
        try:
            omega_string = re.sub("[^0-9]", "", route_params['3108'])
            return float('0.{}'.format(omega_string[1:5]))
        except KeyError:
            try:
                omega_string = re.sub("[^0-9]", "", route_params['iop(3107'].split(',')[0])
                return float('0.{}'.format(omega_string[1:5]))
            except KeyError:
                return None

    def get_solvent_info(self):
        """
        Get solvent information
        """
        try:
            if self.pmgmol.is_pcm:
                self.solvent = self.pcm_info[0]
                self.dielectric_constant = self.pcm_info[1]
        except AttributeError:
            return None

    def parse_tddft(self, sigma=0.10, step=0.01):
        """
        Parse TDDFT data

        :param sigma:
        :param step:
        """
        # read the log files
        self.tddft_excitations = self.get_tddft_excitations(self.log_path)

        # get descriptors
        transitions = self.tddft_excitations["Singlet"] + self.tddft_excitations["Doublet"] + self.tddft_excitations["Triplet"]

        # clean out negative absorptions
        negative_absorptions = [transitions.pop(index) for index, val in enumerate(transitions) if val[0] < 0]
        if negative_absorptions:
            print("WARNING: {} calculation contained a negative excitation energy".format(self.calculation_type))

        # Get plot data
        if transitions:
            connector = {"transitions": "transitions", "sigma": "sigma", "step": "step"}
            c_data = {"transitions": transitions, "sigma": sigma, "step": step}
            self.tddft_spectra_data = DFTSpecPlotter(connector=connector).plot_data(c_data)
