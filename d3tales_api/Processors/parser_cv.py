import numpy as np
import dateutil.parser
from scipy.ndimage import gaussian_filter1d
from d3tales_api.Calculators.plotters import CVPlotter
from d3tales_api.Calculators.calculators import CVDescriptorCalculator, CAResistanceCalculator

MIN_SCAN_LEN = 5


# noinspection PyTypeChecker
def break_list(lst, return_brks=False):
    if len(lst) < 2:
        return lst
    direction = lst[1] - lst[0] / abs(lst[1] - lst[0])
    brks = []
    for i in range(1, len(lst)):
        previous_diff = lst[i] - lst[i - 1]
        if not previous_diff:
            continue
        i_direction = previous_diff / abs(previous_diff)
        if i_direction != direction:
            brks.append(i)
            direction = i_direction
    if return_brks:
        return brks
    return [lst[x:y] for x, y in zip([0] + brks, brks + [None])]


# noinspection PyTypeChecker
def break_scan_data(scan_data):
    scan_data = np.array(scan_data).tolist()
    volts = [v for v, i in scan_data]
    brks = break_list(volts, return_brks=True)
    return [scan_data[x:y] for x, y in zip([0] + brks, brks + [None])]


class ParseChiBase:
    """
    Extract data from raw Chi experiment files
    """

    def __init__(self, file_path, data_header="Potential", delimiter=","):
        """
        :param file_path: str, filepath to experiment data file
        """
        self.file_path = file_path

        with open(self.file_path, "r") as f:
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
                break

        self.all_data = np.loadtxt(self.file_path, delimiter=delimiter, skiprows=line_count+1)

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


class ParseChiCV(ParseChiBase):
    """
    Extract data from raw Chi CV experiment files
    """

    def __init__(self, file_path, data_header="Potential", delimiter=",", min_scan_len=MIN_SCAN_LEN):
        """
        :param file_path: str, filepath to experiment data file
        """
        super().__init__(file_path, data_header, delimiter)

        potentials = self.all_data[:, 0]

        if not getattr(self, "sample_interval", None):
            self.sample_interval = {"value": np.average(np.absolute(np.diff(potentials))), "unit": ''}
        if not getattr(self, "high_e", None):
            self.high_e = {"value": max(potentials), "unit": ''}
        if not getattr(self, "low_e", None):
            self.low_e = {"value": min(potentials), "unit": ''}
        if not getattr(self, "low_e", None):
            self.init_e = {"value": potentials[0], "unit": ''}

        scan_data = [s for s in break_scan_data(self.all_data) if len(s) > min_scan_len]
        self.num_scans = len(scan_data)
        self.scan_data = scan_data
        self.peak_potential = max(potentials)

        self.conditions_data = {
            "data_source": 'cv',
            "scan_rate": getattr(self, 'scan_rate', {}),
            "num_scans": getattr(self, 'num_scans', None),
            "initial_potential": getattr(self, 'init_e', {}),
            "high_e": getattr(self, 'high_e', {}),
            "low_e": getattr(self, 'low_e', {}),
            "comp_r": getattr(self, 'comp_R', {}),
        }
        self.parsed_data = {
            "file_name": getattr(self, 'file_name'),
            "header": getattr(self, 'header'),
            "note": getattr(self, 'note'),
            "date_recorded": getattr(self, 'date_recorded'),
            "conditions": self.conditions_data,
            "segment": getattr(self, 'segment', None),
            "sample_interval": getattr(self, 'sample_interval', {}),
            "quiet_time": getattr(self, 'quiet_time', {}),
            "sensitivity": getattr(self, 'sensitivity', {}),
            "peak_potential": getattr(self, 'init_e', {}),
            "scan_data": getattr(self, 'scan_data', {}),
        }

        self.low_e_value = self.low_e.get("value")
        self.sample_interval_value = self.sample_interval.get("value")

    def calculate_prop(self, prop_name, return_type=dict):
        """
        Calculate a given property using the D3TaLES CV Calculators

        :param prop_name: str, property name
        :param return_type: str, datatype to return
        :return: property or empty datatype
        """
        connector = {
            "intv": "sample_interval_value",
            "low_e": "low_e_value",
            "scan_data": "scan_data"
        }
        try:
            func = getattr(CVDescriptorCalculator(connector=connector), prop_name)
            return func(self)
        except Exception as e:
            print("CVDescriptorCalculator does not have function ", prop_name)
            print("\t" + str(e))
            return return_type()

    def calculate_plotting(self, prop_name, return_type=dict):
        """
        Calculate a given plotting property using the D3TaLES Plotters

        :param prop_name: str, property name
        :param return_type: str, datatype to return
        :return: plotting property or empty datatype
        """
        connector = {"scan_data": "scan_data"}
        try:
            func = getattr(CVPlotter(connector=connector), prop_name)
            return func(self)
        except Exception:
            print("CVPlotter does not have function ", prop_name)
            return return_type()


class ParseChiCVMicro(ParseChiBase):
    """
    Extract data from raw Chi CV experiment files
    """

    def __init__(self, file_path, data_header="Potential", delimiter=",", min_scan_len=MIN_SCAN_LEN):
        """
        :param file_path: str, filepath to experiment data file
        """
        super().__init__(file_path, data_header, delimiter)

        self.potential = self.all_data[:, 0]
        self.current = self.all_data[:, 1]

        if not getattr(self, "sample_interval", None):
            self.sample_interval = {"value": np.average(np.absolute(np.diff(self.potential))), "unit": ''}
        if not getattr(self, "high_e", None):
            self.high_e = {"value": max(self.potential), "unit": ''}
        if not getattr(self, "low_e", None):
            self.low_e = {"value": min(self.potential), "unit": ''}
        if not getattr(self, "low_e", None):
            self.init_e = {"value": self.potential[0], "unit": ''}

        scan_data = [s for s in break_scan_data(self.all_data) if len(s) > min_scan_len]
        self.num_scans = len(scan_data)
        self.scan_data = scan_data
        self.peak_potential = max(self.potential)
        self.i_ss = max(self.current) - min(self.current)
        self.e_half = self.get_e_half(scan_data[0])

        self.conditions_data = {
            "data_source": 'cv',
            "scan_rate": getattr(self, 'scan_rate', {}),
            "num_scans": getattr(self, 'num_scans', None),
            "initial_potential": getattr(self, 'init_e', {}),
            "high_e": getattr(self, 'high_e', {}),
            "low_e": getattr(self, 'low_e', {}),
            "comp_r": getattr(self, 'comp_R', {}),
        }
        self.parsed_data = {
            "file_name": getattr(self, 'file_name'),
            "header": getattr(self, 'header'),
            "note": getattr(self, 'note'),
            "date_recorded": getattr(self, 'date_recorded'),
            "conditions": self.conditions_data,
            "segment": getattr(self, 'segment', None),
            "sample_interval": getattr(self, 'sample_interval', {}),
            "quiet_time": getattr(self, 'quiet_time', {}),
            "sensitivity": getattr(self, 'sensitivity', {}),
            "peak_potential": self.peak_potential,
            "i_ss": self.i_ss,
            "e_half": [self.e_half],
            "scan_data": self.scan_data,
        }

        self.low_e_value = self.low_e.get("value")
        self.sample_interval_value = self.sample_interval.get("value")

    @staticmethod
    def get_e_half(sweep):
        potential = [i[0] for i in sweep]
        current = [i[1] for i in sweep]
        # smooth current data
        smooth_current = gaussian_filter1d(current, 100)
        # compute second derivative
        current_d2 = np.gradient(np.gradient(smooth_current))
        # Find the inflection point
        inflection_points = np.where(np.diff(np.sign(current_d2)))[0]
        if len(inflection_points) != 1:
            raise Exception(f"Error calculating E1/2 with microelectrode. {inflection_points} inflection points found.")
        return potential[inflection_points[0]]

    def calculate_prop(self, prop_name, return_type=dict):
        """
        Calculate a given property using the D3TaLES CV Calculators

        :param prop_name: str, property name
        :param return_type: str, datatype to return
        :return: property or empty datatype
        """
        connector = {
            "intv": "sample_interval_value",
            "low_e": "low_e_value",
            "scan_data": "scan_data"
        }
        try:
            func = getattr(CVDescriptorCalculator(connector=connector), prop_name)
            return func(self)
        except Exception as e:
            print("CVDescriptorCalculator does not have function ", prop_name)
            print("\t" + str(e))
            return return_type()

    def calculate_plotting(self, prop_name, return_type=dict):
        """
        Calculate a given plotting property using the D3TaLES Plotters

        :param prop_name: str, property name
        :param return_type: str, datatype to return
        :return: plotting property or empty datatype
        """
        connector = {"scan_data": "scan_data"}
        try:
            func = getattr(CVPlotter(connector=connector), prop_name)
            return func(self)
        except Exception:
            print("CVPlotter does not have function ", prop_name)
            return return_type()

class ParseChiCA(ParseChiBase):
    """
    Extract data from raw Chi CV experiment files
    """

    def __init__(self, file_path, data_header="Time", delimiter=","):
        """
        :param file_path: str, filepath to experiment data file
        """
        super().__init__(file_path, data_header, delimiter)

        self.t = self.all_data[:, 0]
        self.i = self.all_data[:, 1]
        self.f_slp = self.get_data_calcs(calc_name="Slp", header="Forward:")
        self.f_int = self.get_data_calcs(calc_name="Int", header="Forward:")
        self.f_cor = self.get_data_calcs(calc_name="Cor", header="Forward:")
        self.r_slp = self.get_data_calcs(calc_name="Slp", header="Reverse:")
        self.r_int = self.get_data_calcs(calc_name="Int", header="Reverse:")
        self.r_cor = self.get_data_calcs(calc_name="Cor", header="Reverse:")

        self.conditions_data = {
            "data_source": 'ca',
            "high_e": getattr(self, 'high_e', {}),
            "low_e": getattr(self, 'low_e', {}),
            "step": getattr(self, 'step', {}),
            "pulse_width": getattr(self, 'pulse_width', {}),
        }

        self.parsed_data = {
            "file_name": getattr(self, 'file_name'),
            "header": getattr(self, 'header'),
            "note": getattr(self, 'note'),
            "date_recorded": getattr(self, 'date_recorded'),
            "conditions": self.conditions_data,
            "quiet_time": getattr(self, 'quiet_time', {}),
            "f_slp": self.f_slp,
            "f_int": self.f_int,
            "f_cor": self.f_cor,
            "r_slp": self.r_slp,
            "r_int": self.r_int,
            "r_cor": self.r_cor,
            "resistance": self.resistance,
            "time": self.t.tolist(),
            "current": self.i.tolist(),
        }

    @property
    def resistance(self, **kwargs):
        connector = {
            "i_s": "i",
            "t_s": "t",
            "pulse_width": "pulse_width.value",
            "steps": "step",
            "low_e": "low_e.value",
        }
        r_calculator = CAResistanceCalculator(connector=connector)
        return r_calculator.calculate(self.__dict__, **kwargs)

    def get_data_calcs(self, calc_name="Slp", header="Forward"):
        h_lines = [i for i, l in enumerate(self.lines) if l.startswith(header)]
        if h_lines:
            calcs_lines = self.lines[h_lines[0]+1:h_lines[0]+4]
            calc_line = [l for l in calcs_lines if calc_name in l]
            if calc_line:
                calc = calc_line[0].split()[2]
                return float(calc)


class ParseChiESI(ParseChiBase):
    """
    Extract data from raw Chi CV experiment files
    """

    def __init__(self, file_path, data_header="Freq", delimiter=","):
        """
        :param file_path: str, filepath to experiment data file
        """
        super().__init__(file_path, data_header, delimiter)

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

        self.conditions_data = {
            "data_source": 'esi',
            "high_f": getattr(self, 'high_f', {}),
            "low_f": getattr(self, 'low_f', {}),
            "amp": getattr(self, 'amp', {}),
        }
        self.parsed_data = {
            "file_name": getattr(self, 'file_name'),
            "header": getattr(self, 'header'),
            "note": getattr(self, 'note'),
            "date_recorded": getattr(self, 'date_recorded'),
            "conditions": self.conditions_data,
            "quiet_time": getattr(self, 'quiet_time', {}),
            "all_data": self.all_data,
            "resistance": self.resistance,
        }

