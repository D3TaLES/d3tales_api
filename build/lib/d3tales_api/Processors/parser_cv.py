import numpy as np
import dateutil.parser
from d3tales_api.Calculators.plotters import CVPlotter
from d3tales_api.Calculators.calculators import CVDescriptorCalculator


class ParseChiCV:
    """
    Extract data from raw Chi CV experiment files
    """
    def __init__(self, file_path):
        """
        :param file_path: str, filepath to experiment data file
        """
        self.file_path = file_path

        with open(self.file_path, "r") as f:
            line = f.readline()
            self.date_recorded = dateutil.parser.parse(line.strip()).isoformat()

            segment = 1
            line = f.readline()
            while not line.startswith("Potential"):
                if line.startswith("File"):
                    self.file_name = self.extract_value_unit(line, value_break=":")
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

                line = f.readline()

            line = f.readline()
            while not line.strip().split():
                line = f.readline()
            all_data, potentials = [], []
            while line.strip().split():
                all_data.append([float(x) for x in line.strip().split(",")])
                potentials.append([float(x) for x in line.strip().split(",")][0])
                line = f.readline()

            if not getattr(self, "sample_interval", None):
                self.sample_interval = {"value": np.average(np.absolute(np.diff(potentials))), "unit": ''}
            if not getattr(self, "high_e", None):
                self.high_e = {"value": max(potentials), "unit": ''}
            if not getattr(self, "low_e", None):
                self.low_e = {"value": min(potentials), "unit": ''}
            if not getattr(self, "low_e", None):
                self.init_e = {"value": potentials[0], "unit": ''}

            self.data_points_per_scan = int(abs(self.high_e["value"] - self.low_e["value"]) / self.sample_interval["value"])
            scan_data = [np.array(all_data[i:i + self.data_points_per_scan]).tolist() for i in range(0, len(all_data), self.data_points_per_scan)]

        self.num_scans = len(scan_data)
        self.scan_data = scan_data
        self.peak_potential = max(potentials)

        self.cv_data = {
            "file_name": getattr(self, 'file_name'),
            "header": getattr(self, 'header'),
            "note": getattr(self, 'note'),
            "date_recorded": getattr(self, 'date_recorded'),
            "data_points_per_scan": getattr(self, 'data_points_per_scan', None),
            "segment": getattr(self, 'segment', None),
            "sample_interval": getattr(self, 'sample_interval', {}),
            "quiet_time": getattr(self, 'quiet_time', {}),
            "sensitivity": getattr(self, 'sensitivity', {}),
            "comp_r": getattr(self, 'comp_R', {}),
            "peak_potential": getattr(self, 'init_e', {}),
            "scan_data": getattr(self, 'scan_data', {}),
        }

        self.conditions_data = {
            "scan_rate": getattr(self, 'scan_rate', {}),
            "num_scans": getattr(self, 'num_scans', None),
            "initial_potential": getattr(self, 'init_e', {}),
            "high_e": getattr(self, 'high_e', {}),
            "low_e": getattr(self, 'low_e', {}),
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
        except Exception:
            print("CVDescriptorCalculator does not have function ", prop_name)
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

