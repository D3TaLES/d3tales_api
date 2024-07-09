import uuid
import json
import math
import warnings
import dateutil.parser
from scipy.ndimage import gaussian_filter1d

from d3tales_api.Processors.utils import *
from d3tales_api.Calculators.plotters import CVPlotter
from d3tales_api.Calculators.utils import unit_conversion
from d3tales_api.Calculators.calculators import CVDescriptorCalculator, CAResistanceCalculator

MIN_SCAN_LEN = 5


class ProcessPotBase:
    """
    Base class for processing potentiostat files
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, _id: str = None, submission_info: dict = None, metadata: dict = None):
        """
        :param filepath: str, filepath to
        :param _id: str, molecule ID
        :param submission_info: dict, submission info
        :param metadata: dict, metadata
        """
        self.id = _id
        self.data_source = ""
        self.filepath = filepath or metadata.get("mol_file")
        self.hash_id = str(uuid.uuid4())
        self.submission_info = submission_info or {}

        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.experiment_run_id = metadata.get("experiment_run_id", '')
        self.working_electrode = metadata.get("electrode_working", '')
        self.counter_electrode = metadata.get("electrode_counter", '')
        self.reference_electrode = metadata.get("electrode_reverence", '')

        self.e_rev = unit_conversion(metadata.get("e_rev"), default_unit="V", return_dict=True)
        self.temperature = unit_conversion(metadata.get("temperature"), default_unit="K", return_dict=True)
        self.redox_mol_concentration = unit_conversion(metadata.get("redox_mol_concentration"),
                                                       default_unit="M", return_dict=True)
        self.electrolyte_concentration = unit_conversion(metadata.get("electrolyte_concentration"),
                                                         default_unit="M", return_dict=True)
        self.working_electrode_radius = unit_conversion(metadata.get("working_electrode_radius"),
                                                        default_unit="cm", return_dict=True)
        self.working_electrode_surface_area = unit_conversion(metadata.get("working_electrode_surface_area"),
                                                              default_unit="cm^2", return_dict=True)
        if self.working_electrode_radius and not self.working_electrode_surface_area:
            self.working_electrode_surface_area = dict(value=math.pi * (self.working_electrode_radius["value"] ** 2),
                                                       unit="cm^2")
        self.solvents = metadata.get("solvent") if isinstance(metadata.get("solvent"), list) else [
            metadata.get("solvent")] if metadata.get("solvent") else []
        self.electrolytes = metadata.get("electrolyte") if isinstance(metadata.get("electrolyte"), list) else [
            metadata.get("electrolyte")] if metadata.get("electrolyte") else []
        self.ionic_liquids = metadata.get("ionic_liquid") if isinstance(metadata.get("ionic_liquid"), list) else [
            metadata.get("ionic_liquid")] if metadata.get("ionic_liquid") else []

    @property
    def data_dict(self):
        """
        Dictionary of processed data (in accordance with D3TalES backend schema)
        """
        all_data_dict = {
            "_id": self.hash_id,
            "mol_id": self.id,
            "submission_info": self.submission_info,
            "data": self.parsed_data
        }
        json_data = json.dumps(all_data_dict)
        return json.loads(json_data)

    @property
    def parsed_data(self):
        """
        Dictionary of parsed data
        """
        p_data = {
            "file_name": getattr(self, 'file_name', None),
            "header": getattr(self, 'header', None),
            "note": getattr(self, 'note', None),
            "date_recorded": getattr(self, 'date_recorded', None),
            "conditions": self.conditions,
            "segment": getattr(self, 'segment', None),
            "sample_interval": getattr(self, 'sample_interval', None),
            "quiet_time": getattr(self, 'quiet_time', None),
            "sensitivity": getattr(self, 'sensitivity', None),
            "scan_data": getattr(self, 'scan_data', None),
        }
        return {k: v for k, v in p_data.items() if v is not None}

    @property
    def conditions(self):
        """
        Dictionary of conditions (in accordance with D3TaLES backend schema)
        """
        conditions_data = {
            "data_source": self.data_source,
            "working_electrode": self.working_electrode,
            "counter_electrode": self.counter_electrode,
            "reference_electrode": self.reference_electrode,
            "solvent": self.reagent_format(self.solvents),
            "electrolyte": self.reagent_format(self.electrolytes),
            "ionic_liquid": self.reagent_format(self.ionic_liquids),
            "instrument": self.instrument,
            "working_electrode_radius": self.working_electrode_radius,
            "working_electrode_surface_area": self.working_electrode_surface_area,
            "redox_mol_concentration": self.redox_mol_concentration,
            "electrolyte_concentration": self.electrolyte_concentration,
            "temperature": self.temperature,
            "experiment_run_id": self.experiment_run_id,

            "scan_rate": getattr(self, 'scan_rate', None),
            "num_scans": getattr(self, 'num_scans', None),
            "initial_potential": getattr(self, 'init_e', None),
            "high_e": getattr(self, 'high_e', None),
            "low_e": getattr(self, 'low_e', None),
            "comp_r": getattr(self, 'comp_R', None),
            "step": getattr(self, 'step', None),
            "pulse_width": getattr(self, 'pulse_width', None),
            "high_f": getattr(self, 'high_f', None),
            "low_f": getattr(self, 'low_f', None),
            "amp": getattr(self, 'amp', None),
        }
        return {k: v for k, v in conditions_data.items() if v is not None}

    @staticmethod
    def reagent_format(reagents, purity=None):
        """
        Convert reagent data to D3TaLES schema format

        :param reagents: reagent data or list of reagent data
        :param purity: default purity value
        :return: formatted reagent data
        """

        def format_i(r):
            if isinstance(r, dict):
                r_dict = r
            else:
                r_dict = {"name": str(r)}
                if purity:
                    r_dict["purity"] = purity
            return r_dict

        return [format_i(i) for i in reagents] if isinstance(reagents, list) else format_i(reagents)


class ParseChiMixin:
    """
    Extract data from raw Chi experiment files
    """

    def parse(self, data_header="Potential", delimiter=","):
        """
        :param filepath: str, filepath to experiment data file
        :param data_header: str,
        :param delimiter: str, 
        """

        with open(self.filepath, "r") as f:
            self.lines = f.readlines()

        line_count = 0
        for line_count, line in enumerate(self.lines):
            if line_count == 0:
                self.date_recorded = dateutil.parser.parse(line.strip()).isoformat()

            if line.startswith("File"):
                self.file_name = self.extract_value_unit(line, value_break=":")
            elif line.startswith("Instrument Mode"):
                self.instrument = self.extract_value_unit(line, value_break=":")
            elif line.startswith("Header"):
                self.header = self.extract_value_unit(line, value_break=":")
            elif line.startswith("Note"):
                self.note = self.extract_value_unit(line, value_break=":")
            elif line.startswith("Init E"):
                self.init_e = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("High E"):
                self.high_e = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Low E"):
                self.low_e = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Init P/N"):
                self.init_p_n = self.extract_value_unit(line)
            elif line.startswith("Scan Rate"):
                self.scan_rate = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Segment ="):
                self.segment = self.extract_value_unit(line, value_type='int')
            elif line.startswith("Sample Interval"):
                self.sample_interval = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Quiet Time"):
                self.quiet_time = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Sensitivity"):
                self.sensitivity = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Comp R"):
                self.comp_R = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("High Frequency"):
                self.high_f = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Low Frequency"):
                self.low_f = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Amplitude"):
                self.amp = self.extract_value_unit(line, value_type='float', return_dict=True)
            elif line.startswith("Step"):
                self.step = self.extract_value_unit(line, value_type='float', return_dict=False)
            elif line.startswith("Pulse Width"):
                self.pulse_width = self.extract_value_unit(line, value_type='float', return_dict=True)

            if line.startswith(data_header):
                self.x_unit = line.split(",")[0].split("/")[1].strip()
                self.y_unit = line.split(",")[1].split("/")[1].strip()
                break

        self.all_data = np.loadtxt(self.filepath, delimiter=delimiter, skiprows=line_count + 1)

    @staticmethod
    def extract_value_unit(line, value_break="=", value_type='str', return_dict=False):
        """
        Extract data value and unit from a Chi CV data line

        :param line: str, data line
        :param value_break: str, type of character which breaks the line
        :param value_type: str, value datatype
        :param return_dict: bool, return data dict if True
        :return: value or data dict
        """
        if value_type == 'float':
            value = float("".join(line.split(value_break)[1:]).strip())
        elif value_type == 'int':
            value = int("".join(line.split(value_break)[1:]).strip())
        else:
            value = ("".join(line.split(value_break)[1:]).strip())
        unit = line[line.find("(") + 1:line.find(")")]
        if return_dict:
            return {"value": value, "unit": unit}
        else:
            return value


class ProcessChiCV(ParseChiMixin, ProcessPotBase):
    """
    Extract data from raw Chi CV experiment files
    """

    def __init__(self, filepath, _id: str = None, submission_info: dict = None, metadata: dict = None,
                 min_scan_len=MIN_SCAN_LEN, micro_electrodes=False, **kwargs):
        """
        :param filepath: str, filepath to experiment data file
        """
        super().__init__(filepath, _id=_id, submission_info=submission_info, metadata=metadata)
        self.data_source = "cv"
        self.micro_electrodes = micro_electrodes
        self.parse(data_header="Potential", delimiter=",")

        self.potential = self.all_data[:, 0]
        self.current = self.all_data[:, 1]

        if not getattr(self, "sample_interval", None):
            self.sample_interval = {"value": np.average(np.absolute(np.diff(self.potential))), "unit": self.x_unit}
        if not getattr(self, "high_e", None):
            self.high_e = {"value": max(self.potential), "unit": self.x_unit}
        if not getattr(self, "low_e", None):
            self.low_e = {"value": min(self.potential), "unit": self.x_unit}
        if not getattr(self, "low_e", None):
            self.init_e = {"value": self.potential[0], "unit": self.x_unit}

        self.scan_data = [s for s in break_scan_data(self.all_data) if len(s) > min_scan_len]
        self.num_scans = len(self.scan_data)
        self.plot_data = CVPlotter().plot_data(self.__dict__).get("abs_plot")

        if micro_electrodes:
            self.i_ss = {"value": max(self.current) - min(self.current), "unit": self.y_unit}
            self.e_half = {"value": self.get_micro_e_half(self.scan_data[0]), "unit": self.x_unit}
        else:
            self.reversibility = CVDescriptorCalculator().reversibility(self.__dict__)
            self.e_half = {"value": CVDescriptorCalculator().e_half(self.__dict__), "unit": self.x_unit}
            self.peak_splittings = {"value": CVDescriptorCalculator().peak_splittings(self.__dict__), "unit": self.x_unit}
            self.middle_sweep = CVDescriptorCalculator().middle_sweep(self.__dict__)

    @property
    def parsed_data(self):
        p_data = super().parsed_data
        p_data.update({
            "peak_potential": getattr(self, 'high_e', None),
            "e_half": getattr(self, 'e_half', None),
            "plot_data": getattr(self, 'plot_data', None),
        })
        if self.micro_electrodes:
            p_data.update({
                "i_ss": getattr(self, 'i_ss', None),
                "e_rev": getattr(self, 'e_rev', None),
            })
        else:
            p_data.update({
                "reversibility": getattr(self, 'reversibility', None),
                "peak_splittings": getattr(self, 'peak_splittings', None),
                "middle_sweep": getattr(self, 'middle_sweep', None),
            })
        return p_data

    @staticmethod
    def inflection_points(y_values, sigma=10):
        # smooth current data
        smooth_current = np.array(gaussian_filter1d(y_values, sigma))
        # compute second derivative
        current_d2 = np.gradient(np.gradient(smooth_current))
        # Find the inflection point
        return np.where(np.diff(np.sign(current_d2)))[0]

    def get_micro_e_half(self, sweep, sigma=None):
        potential = [i[0] for i in sweep]
        current = [i[1] for i in sweep]
        test_sigmas = [sigma] if sigma else [5, 10, 20, 50, 100, 150, 200]
        for sigma in test_sigmas:
            inflection_points = self.inflection_points(current, sigma=sigma)
            if len(inflection_points) == 1:
                print(f"Inflection point {potential[inflection_points[0]]} found with sigma {sigma} ")
                return potential[inflection_points[0]]
            print(f"{[potential[i] for i in inflection_points]} inflection points found with {sigma} sigma.")
        warnings.warn(f"Error calculating E1/2 with microelectrode data. Sigma {test_sigmas} tested. ")


class ProcessChiCA(ParseChiMixin, ProcessPotBase):
    """
    Extract data from raw Chi CV experiment files
    """

    def __init__(self, filepath, _id: str = None, submission_info: dict = None, metadata: dict = None, **kwargs):
        """
        :param filepath: str, filepath to experiment data file
        """
        super().__init__(filepath, _id=_id, submission_info=submission_info, metadata=metadata)
        self.data_source = "ca"
        self.cell_constant = metadata.get("cell_constant", 1)
        self.parse(data_header="Time", delimiter=",")

        self.t = self.all_data[:, 0]
        self.i = self.all_data[:, 1]
        if self.x_unit != "s":
            pass  # TODO convert to s and A if not already there.
        if self.y_unit != "A":
            pass
        self.f_slp = self.get_data_calcs(calc_name="Slp", header="Forward:")
        self.f_int = self.get_data_calcs(calc_name="Int", header="Forward:")
        self.f_cor = self.get_data_calcs(calc_name="Cor", header="Forward:")
        self.r_slp = self.get_data_calcs(calc_name="Slp", header="Reverse:")
        self.r_int = self.get_data_calcs(calc_name="Int", header="Reverse:")
        self.r_cor = self.get_data_calcs(calc_name="Cor", header="Reverse:")

        self.measured_resistance = self.get_resistance()
        self.measured_conductance = 1 / self.measured_resistance

    @property
    def parsed_data(self):
        p_data = super().parsed_data
        p_data.update({
            "f_slp": self.f_slp,
            "f_int": self.f_int,
            "f_cor": self.f_cor,
            "r_slp": self.r_slp,
            "r_int": self.r_int,
            "r_cor": self.r_cor,
            "measured_resistance": {"value": self.measured_resistance, "unit": "Ohm"},
            "measured_conductance": {"value": self.measured_conductance, "unit": "S"},
            "conductivity": self.conductivity,
            "time": self.t.tolist(),
            "current": self.i.tolist(),
        })
        return p_data

    @property
    def conductivity(self):
        """
        Get conductivity based on cell constant
        """
        c_measured = unit_conversion(self.measured_conductance, default_unit="mS")
        cell_constant = unit_conversion(self.cell_constant, default_unit="cm^-1")
        return {"value": c_measured * cell_constant, "unit": "mS/cm"}

    def get_resistance(self, **kwargs):
        connector = {
            "i_s": "i",
            "t_s": "t",
            "pulse_width": "pulse_width",
            "steps": "step",
            "low_e": "low_e",
        }
        r_calculator = CAResistanceCalculator(connector=connector)
        return r_calculator.calculate(self.__dict__, **kwargs)

    def get_data_calcs(self, calc_name="Slp", header="Forward"):
        h_lines = [i for i, l in enumerate(self.lines) if l.startswith(header)]
        if h_lines:
            calcs_lines = self.lines[h_lines[0] + 1:h_lines[0] + 4]
            calc_line = [l for l in calcs_lines if calc_name in l]
            if calc_line:
                calc = calc_line[0].split()[2]
                return float(calc)


class ProcessChiESI(ParseChiMixinProcessPotBase):
    """
    Extract data from raw Chi ESI experiment files
    """

    def __init__(self, filepath, _id: str = None, submission_info: dict = None, metadata: dict = None, **kwargs):
        """
        :param filepath: str, filepath to experiment data file
        """
        super().__init__(filepath, _id=_id, submission_info=submission_info, metadata=metadata)
        self.data_source = "esi"
        self.parse(data_header="Freq", delimiter=",")
        self.scan_data = self.all_data

        self.freq = self.all_data[:, 0]
        self.z_real = self.all_data[:, 1]
        self.z_img = self.all_data[:, 2]
        self.z = self.all_data[:, 3]
        self.phase = self.all_data[:, 4]

        self.resistance = min(zip(self.z_real, self.z_img), key=lambda x: abs(x[1]))[0]

        if not getattr(self, "high_f", None):
            self.high_f = {"value": max(self.freq), "unit": 'Hz'}
        if not getattr(self, "low_f", None):
            self.low_f = {"value": min(self.freq), "unit": 'Hz'}

    @property
    def parsed_data(self):
        p_data = super().parsed_data
        p_data.update({"resistance": self.resistance})
        return p_data
