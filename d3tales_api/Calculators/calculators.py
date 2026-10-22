from d3tales_api.Calculators.utils import *
from d3tales_api.Calculators.utils import periodictable
import cclib
import abc
import json
import math
import warnings
import numpy as np
import pandas as pd
from rdkit import Chem
import dbstep.Dbstep as db
from rdkit.Chem import rdMolAlign
from interruptingcow import timeout
from scipy.signal import find_peaks
from pymatgen.core.sites import Site
from pymatgen.core.structure import Molecule
from rdkit.Chem.Descriptors import ExactMolWt
from d3tales_api.Calculators.ocelot_transform import pmgmol_to_rdmol
from scipy.stats import linregress

class D3Calculator(abc.ABC):
    """
    D3Calculators base class

    :param connector: dictionary describing connections between calculator variables (keys) and location in `data`
    variable. EX: {"smiles": "mol_info.smiles", "volume": "experimental_data.volume"}. If no connector is provided,
    the connections default to the connections found in `default_connector.json`.
    """

    def __init__(self, connector=None):
        if connector:
            self.key_pairs = connector
        else:
            self.key_pairs = {k: k for k in self.all_connections}

    def make_connections(self, obj):
        d = {}
        for key, connection in self.key_pairs.items():
            try:
                reg_val = rgetattr(obj, connection)
            except:
                try:
                    reg_val = rgetkeys(obj, connection)
                except:
                    reg_val = None
            value = reg_val
            if isinstance(reg_val, str):
                if key == reg_val:
                    value = None
            if value is not None:
                d.update({key: value})
        return d

    def calculate(self, data):
        pass

    def description(self):
        print(self.__class__.__name__)

    @property
    def all_connections(self):
        return {
            "A": "electrode area (default = cm^2)",
            "C": "concentration of the solution (default = mol/cm^3)",
            "D": "Diffusion constant (default = cm^2/s)",
            "T": "Temperature (default = 293 K)",
            "X": "peak shift (default = V)",
            "e": "E1/2 (default = V)",
            "e_ref": "formal potential (default = V)",
            "electrode": "electrode name as str or potential as float (default = standard_hydrogen_electrode)",
            "energy": "energy another geometry (default = eV)",
            "energy_final": "energy final (default = eV)",
            "energy_initial": "energy initial (default = eV)",
            "fin_corr": "final entropy correction (default = eV)",
            "fin_eng": "final energy (default = eV)",
            "fin_eng_solv": "final energy of solvation (default = eV)",
            "geom_final": "geometry final (default = A)",
            "geom_initial": "geometry initial (default = A)",
            "gs_energy": "ground state energy at ion geometry (default = eV)",
            "gs_opt": "ground state optimized energy (default = eV)",
            "i_s": "list, current points (s)",
            "init_corr": "initial entropy correction (default = eV)",
            "init_eng": "initial energy (default = eV)",
            "init_eng_solv": "initial energy of solvation (default = eV)",
            "ion_energy": "ion energy at ground state geometry (default = eV)",
            "ion_opt": "ion optimized energy (default = eV)",
            "log_file": "calculation output file. Must be readable with CCLIB",
            "low_e": "float, lowest voltage (V)",
            "middle_scan": "middle of scan data",
            "n": "number of electrons, default 1",
            "num_electrons": "number of electrons (default = 1)",
            "opt_energy": "optimized energy (default = eV)",
            "pulse_width": "float, time width of a pulse (s)",
            "redox_density": "density for redox-active molecule",
            "sample_interval": "sample interval value (V)",
            "scan_data": "optional, if e not provided, scan data will be used to find e (default = None)",
            "smiles": "SMILES string",
            "solv_density": "density for solvent",
            "spin_type": "type of CCLIB spin to extract (default = Mulliken)",
            "steps": "int, number of potentiostat sweeps",
            "t_s": "list, time points (s)",
            "v": "scan rate (default = V/s)",
            "volume": "Volume",
            "weight": "actual weight",
        }


# ----------------------- Cyclic Voltammetry Calculators ------------------------------
class ConcentrationCalculator(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Diffusion constant using Randles-Scidwick equation

        Connection Points:
            :smiles: SMILES string
            :volume: Volume
            :weight: actual weight
            :solv_density: density for solvent
            :redox_density: density for redox-active molecule

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: concentrations (eV)
        """
        self.data = data

        conns = self.make_connections(data)
        solv_density = conns.get("solv_density")
        redox_density = conns.get("redox_density")
        mol_weight = ExactMolWt(Chem.MolFromSmiles(conns["smiles"]))
        concentration = ((unit_conversion(conns["weight"], default_unit='g', density=redox_density) / unit_conversion(
            mol_weight, default_unit='g/mol')) /
                         unit_conversion(conns["volume"], default_unit='L', density=solv_density))
        return float(np.format_float_scientific(concentration, precision=precision))


class CVDescriptorCalculator(D3Calculator):

    def peaks(self, data: dict, width: float = 1, middle_sweep=True):
        """
        Gather CV peaks

        Connection Points:
            :scan_data: scanned data from CV file

        :param data: data for calculation (rows of voltage, current)
        :param width: required width of peaks to identify
        :param middle_sweep: ensure the middle sweep in analyzed if True
        :type data: dict
        :type width: float
        :type middle_sweep: bool

        :return: dictionary containing list of forward peaks and list of reverse peaks
        """

        if middle_sweep:
            scan_data = self.middle_sweep(data)
        else:
            self.data = data
            conns = self.make_connections(data)
            scan_data = conns["scan_data"]

        peak_dict = {}
        for data_list in scan_data:
            data = np.array(data_list)
            if data[0, 0] < data[-1, 0]:
                try:
                    peaks_data = find_peaks(data[:, 1], width=width)
                    f_peaks = peak_dict.get("forward", []) + self.prominent_peaks(peaks_data, data)
                    peak_dict.update({"forward": f_peaks})
                except ValueError:
                    pass
            else:
                try:
                    peaks_data = find_peaks(-data[:, 1], width=width)
                    r_peaks = peak_dict.get("reverse", []) + self.prominent_peaks(peaks_data, data)
                    peak_dict.update({"reverse": r_peaks})
                except ValueError:
                    pass
        return peak_dict

    def peaks_for_analysis(self, data: dict, get_max=True, cut_extras=False, **kwargs):
        """
        Get peaks for analysis

        :param data: data for calculation
        :param cut_extras: cut extra peaks off the end of the appropriate array if True and there are
            not the same number of forward and reverse peaks
        :param kwargs:
        :return: forward_peaks, reverse_peaks
        """
        peaks = self.peaks(data, **kwargs)
        forward_peaks = [item[0] for item in peaks.get('forward', [])]
        reverse_peaks = [item[0] for item in peaks.get('reverse', [])]

        # Check if there are the same number of forward and reverse peaks
        if len(forward_peaks) != len(reverse_peaks):
            num_peaks = min([len(forward_peaks), len(reverse_peaks)])
            if get_max:
                f_idx = np.argsort([item[1] for item in peaks.get('forward', [])])[-num_peaks:]
                r_idx = np.argsort([item[1] for item in peaks.get('reverse', [])])[-num_peaks:]
                forward_peaks = [forward_peaks[i] for i in f_idx]
                reverse_peaks = [reverse_peaks[i] for i in r_idx]
            elif cut_extras:
                forward_peaks, reverse_peaks = forward_peaks[:num_peaks], reverse_peaks[:num_peaks]
            else:
                raise ValueError("Error. There are {} forward peaks and {} reverse peaks. Either examine data or set "
                                 "get_max or cut_extras kwarg to True.".format(len(forward_peaks), len(reverse_peaks)))

        return sorted(forward_peaks), sorted(reverse_peaks)

    def reversibility(self, data: dict, rev_upperbound: float = 63, quasi_rev_upperbound: float = 200, **kwargs):
        """
        Categorization of CV reversibility

        Connection Points:
            :scan_data: scanned data from CV file (potentials in V)
            :sample_interval: sample interval value (V)
            :low_e: lowest energy value (V)

        :param data: data for calculation
        :param rev_upperbound : upperbound for reversibility (mV)
        :param quasi_rev_upperbound : upperbound for quasi reversibility (mV)
        :type data: dict
        :type rev_upperbound: float
        :type quasi_rev_upperbound: float

        :return: list of reversibility categorizations for peaks
        """

        forward_peaks, reverse_peaks = self.peaks_for_analysis(data, **kwargs)
        peaks_reversibility = []
        for f_peak, r_peak in zip(forward_peaks, reverse_peaks):
            delta_e = abs(f_peak - r_peak)
            if delta_e <= rev_upperbound / 1000:
                peaks_reversibility.append('reversible')
            elif rev_upperbound / 1000 <= delta_e <= quasi_rev_upperbound / 1000:
                peaks_reversibility.append('quasi-reversible')
            else:
                peaks_reversibility.append('irreversible')
        return peaks_reversibility

    def e_half(self, data: dict, **kwargs):
        """
        Get CV E 1/2

        Connection Points:
            :scan_data: scanned data from CV file (potentials in V)
            :sample_interval: sample interval value (V)
            :low_e: lowest energy value (V)

        :param data: data for calculation
        :type data: dict

        :return: list of E 1/2 for peaks
        """

        forward_peaks, reverse_peaks = self.peaks_for_analysis(data, **kwargs)
        e_halfs = []
        for f_peak, r_peak in zip(forward_peaks, reverse_peaks):
            e_halfs.append(round((f_peak + r_peak) / 2, 3))
        return e_halfs

    def peak_splittings(self, data: dict, **kwargs):
        """
        CV peak splitting

        Connection Points:
            :scan_data: scanned data from CV file (potentials in V)
            :sample_interval: sample interval value (V)
            :low_e: lowest energy value (V)

        :param data: data for calculation
        :type data: dict

        :return: list of peak splittings for peaks
        """

        forward_peaks, reverse_peaks = self.peaks_for_analysis(data, **kwargs)
        splittings = []
        for f_peak, r_peak in zip(forward_peaks, reverse_peaks):
            splittings.append(round(abs(f_peak - r_peak), 3))
        return splittings

    def middle_sweep(self, data):
        """
        CV middle sweep

        Connection Points:
            :scan_data: scanned data from CV file (potentials in V)

        :param data: data for calculation
        :type data: dict

        :return: middle sweep from the CV
        """

        self.data = data
        conns = self.make_connections(data)

        middle_idx = len(conns["scan_data"]) / 2
        return conns["scan_data"][int(middle_idx - 1):int(middle_idx + 1)]

    @staticmethod
    def prominent_peaks(peaks_data: list, orig_data: dict, cutoff: float = 0.0999):
        """
        Get prominent peaks from the CV data
        :param peaks_data: output data from scipy peak `find_peaks` function
        :param orig_data: original peak data
        :param cutoff: percentage as decimal for peaks to disregard
        :type peaks_data: list
        :type orig_data: dict
        :type cutoff: float
        :return:
        """
        prominences = peaks_data[1]["prominences"]
        peaks_list = list(np.where(prominences / max(prominences) > cutoff)[0])
        values = []
        for idx in peaks_list:
            values.append(list(orig_data[[peaks_data[0][idx]], :][0]))
        return values

    def peak_currents(self, data: dict, cathodic_anodic: str = 'cathodic', percent_for_baseline: float = 0.2):
        """
        Gather CV peaks

        Connection Points:
            :scan_data: scanned data from CV file

        :param data: data for calculation (rows of voltage, current)
        :type data: dict
        :param cathodic_anodic:
        :type cathodic_anodic: str
        :param percent_for_baseline:
        :type percent_for_baseline: float

        :return: dictionary containing list of forward peaks and list of reverse peaks
        """

        forward_data, reverse_data = self.middle_sweep(data)
        forward_sweep = pd.DataFrame(forward_data, columns=['potential', 'current'])
        reverse_sweep = pd.DataFrame(reverse_data, columns=['potential', 'current'])

        full_cv_data = pd.concat([forward_sweep, reverse_sweep], ignore_index=True)

        if cathodic_anodic == "cathodic":
            # Identify baseline region and fit regression for cathodic peak determination
            forward_baseline_region = forward_sweep.iloc[:int(len(forward_sweep) * percent_for_baseline)]
            forward_slope, forward_intercept, _, _, _ = linregress(forward_baseline_region['potential'],
                                                                   forward_baseline_region['current'])

            # Extend forward baseline for plotting
            forward_baseline_potentials = np.linspace(forward_sweep['potential'].min(),
                                                      forward_sweep['potential'].max(), 100)

            # Find cathodic peak (maximum current for IUPAC convention)
            cathodic_peak_idx = full_cv_data['current'].idxmax()
            cathodic_peak_potential = full_cv_data['potential'][cathodic_peak_idx]
            cathodic_peak_current = full_cv_data['current'][cathodic_peak_idx]

            # Calculate baseline current at cathodic peak potential
            cathodic_baseline_current = forward_slope * cathodic_peak_potential + forward_intercept
            cathodic_peak_relative_current = cathodic_peak_current - cathodic_baseline_current
            return cathodic_peak_relative_current
        if cathodic_anodic == "anodic":
            # Identify baseline region and fit regression for anodic peak determination
            reverse_baseline_region = reverse_sweep.iloc[:int(len(reverse_sweep) * percent_for_baseline)]
            reverse_slope, reverse_intercept, _, _, _ = linregress(reverse_baseline_region['potential'],
                                                                   reverse_baseline_region['current'])

            # Extend reverse baseline for plotting
            reverse_baseline_potentials = np.linspace(reverse_sweep['potential'].min(),
                                                      reverse_sweep['potential'].max(), 100)

            # Find anodic peak (minimum current for IUPAC convention)
            anodic_peak_idx = full_cv_data['current'].idxmin()
            anodic_peak_potential = full_cv_data['potential'][anodic_peak_idx]
            anodic_peak_current = full_cv_data['current'][anodic_peak_idx]

            # Calculate baseline current at anodic peak potential
            anodic_baseline_current = reverse_slope * anodic_peak_potential + reverse_intercept
            anodic_peak_relative_current = anodic_peak_current - anodic_baseline_current
            return anodic_peak_relative_current
        else:
            raise ValueError(f"Function peak_currents requires argument cathodic_anodic to be either 'cathodic' or "
                             f"'anodic'. Instead, cathodic_anodic={cathodic_anodic}.")


class CVDiffusionCalculator(D3Calculator):

    def calculate(self, data: list, precision: int = 3, sci_notation: bool = False, cathodic_anodic: str = 'cathodic'):
        """
        Diffusion constant using Randles-Scidwick equation

        Connection Points:
            :i_p: peak current (default = A)
            :A: electrode area (default = cm^2)
            :v: scan rate (default = V/s)
            :n: number of electrons, default 1
            :C: concentration of the solution (default = mol/cm^3)
            :scan_data: if i_p is not provided and middle_scan not provided, scan data will be used to find i_p

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :param sci_notation: return in scientific notation if True
        :param cathodic_anodic:
        :type data: list
        :type precision: int
        :type sci_notation: bool
        :type sci_notation: str

        :return: average diffusion constant for single redox event (cm^2/s)
        """
        self.data = data
        self.n = data.__len__()

        diffusion_constants = np.zeros(self.n)
        i_ps = np.zeros(self.n)
        vs = np.zeros(self.n)
        for idx, obj in enumerate(self.data):
            conns = self.make_connections(obj)
            descriptor_cal = CVDescriptorCalculator(connector=self.key_pairs)
            if conns.get("i_p"):
                i_p_raw = conns["i_p"]
            else:
                i_p_raw = descriptor_cal.peak_currents(obj, cathodic_anodic=cathodic_anodic)
            i_p = unit_conversion(i_p_raw, default_unit='A')
            A = unit_conversion(conns["A"], default_unit='cm^2')
            v = unit_conversion(conns["v"], default_unit='V/s')
            C = unit_conversion(conns["C"], default_unit='mol/cm^3')
            n = conns.get("n", 1)
            i_ps[idx] = i_p
            vs[idx] = pow(v, 1 / 2)
            diffusion_constants[idx] = (
                pow(i_p / (2.692e5 * pow(n, (3 / 2)) * A * C * pow(v, 1 / 2)), 2))
        slope = np.polyfit(vs, i_ps, 1)[0]
        diffusion_fitted = pow(slope / (2.692e5 * pow(n, (3 / 2)) * A * C), 2)

        results = [np.format_float_scientific(diffusion_constants.mean(), precision=precision),
                   np.format_float_scientific(diffusion_fitted, precision=precision)]
        print("DIFFUSION: ", results)
        return results if sci_notation else [float(x) for x in results]


class CVDiffusionCalculatorMicro(D3Calculator):

    def calculate(self, data: dict, precision: int = 3, sci_notation: bool = False):
        """
        Diffusion constant calculation for CV gathered with ultramicroelectrodes


        Connection Points:
            :i_ss: steady state current (default = A)
            :n: number of electrons, default 1
            :C: concentration of the solution (default = mol/cm^3)
            :r: radius of ultramicroelectrode (default = cm)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :param sci_notation: return in scientific notation if True
        :type data: list
        :type precision: int
        :type sci_notation: bool

        :return: average diffusion constant for single redox event (cm^2/s)
        """
        self.data = data
        self.n = data.__len__()
        conns = self.make_connections(data)

        i_ss = unit_conversion(conns["i_ss"], default_unit='A')
        n = conns["n"]
        C = unit_conversion(conns["C"], default_unit='mol/cm^3')
        r = unit_conversion(conns["r"], default_unit='cm')
        F = 96485.3321  # Faraday constant, s A / mol
        D = i_ss / (4 * n * F * C * r)

        print(conns)
        print(i_ss, 1000*i_ss/(2*r*math.pi))
        print(D)
        return np.format_float_scientific(D) if sci_notation else float(D)


class CVChargeTransferCalculator(D3Calculator):

    def calculate(self, data: list, precision: int = 3, sci_notation: bool = False):
        """
        Charge transfer rate calculation

        Connection Points:
            :T: Temperature (default = 293 K)
            :D: Diffusion constant (default = cm^2/s)
            :v: scan rate (default = V/s)
            :n: number of electrons
            :X: peak shift (default = V)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :param sci_notation: return in scientific notation if True
        :type data: list
        :type precision: int
        :type sci_notation: bool

        :return: array[["scan_rate", "k_0_i"]] for ith CV (cm/s)
        """
        self.data = data
        self.n = data.__len__()

        psis = []
        long_terms = []
        for idx, obj in enumerate(self.data):
            conns = self.make_connections(obj)
            X = unit_conversion(conns["X"], default_unit='V') * 1000  # convert to mV
            D = unit_conversion(conns["D"], default_unit='cm^2/s')
            v = unit_conversion(conns["v"], default_unit='V/s')
            T = unit_conversion(conns["T"], default_unit='K') or 298
            psis.append((-0.6288 + 0.0021 * X) / (1 - 0.017 * X))
            long_terms.append(pow((math.pi * D * conns.get("n", 1) * 96485.34 * v) / (8.3145 * T), - 1 / 2))
        chargetranfer_rate = np.polyfit(long_terms, psis, 1)[0]

        result = np.format_float_scientific(chargetranfer_rate, precision=precision)
        return result if sci_notation else round(float(result), precision)


class CVChargeTransferCalculatorMicro(D3Calculator):

    def calculate(self, data: dict, precision: int = 3, sci_notation: bool = False):
        """
        Charge transfer rate calculation for CVs gathered with ultramicroelectrodes

        Connection Points:
            :e_half: measured half-wave potential of steady-state sigmoid (default = V)
            :e_ref: formal potential of redox couple in reversible conditions (default = V)
            :n: number of electrons
            :T: Temperature (default = 293 K)
            :D: Diffusion constant (default = cm^2/s)
            :r: radius of ultra-microelectrode (default = cm)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :param sci_notation: return in scientific notation if True
        :type data: list
        :type precision: int
        :type sci_notation: bool

        :return: charge transfer rate (cm/s)
        """
        self.data = data
        self.n = data.__len__()
        conns = self.make_connections(data)

        e_half = unit_conversion(conns["e_half"], default_unit='V')
        e_ref = unit_conversion(conns["e_ref"], default_unit='V')
        n = conns["n"]
        T = unit_conversion(conns["T"], default_unit='K')
        D = unit_conversion(conns["D"], default_unit='cm^2/s')
        r = unit_conversion(conns["r"], default_unit='cm')
        F = 96485.3321  # Faraday constant, s A / mol
        R = 8.314  # molar gas constant, J / K / mol
        k = 2 * D / (r * (1 - (math.exp(n * F * (e_half - e_ref) / (R * T)))))
        return np.format_float_scientific(k) if sci_notation else float(k)


class AvgEHalfCalculator(D3Calculator):

    def calculate(self, data: list, precision: int = 3, sci_notation: bool = False):
        """
        Average E half calculation

        Connection Points
            :e: E1/2 (default = V)
            :scan_data: optional, if e not provided, scan data will be used to find e (default = None)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :param sci_notation: return in scientific notation if True
        :type data: dict
        :type precision: int
        :type sci_notation: bool

        :return: average E1/2 (in units V)
        """
        self.data = data
        self.n = data.__len__()

        e_halfs = []
        for idx, obj in enumerate(self.data):
            conns = self.make_connections(obj)
            descriptor_cal = CVDescriptorCalculator(connector=self.key_pairs)
            e_raw = conns.get("e", descriptor_cal.e_half(obj)[0])
            E = unit_conversion(e_raw, default_unit='V')
            e_halfs.append(E)
        avg_e_half = np.average(e_halfs)

        result = np.format_float_scientific(avg_e_half, precision=precision)
        return result if sci_notation else float(result)


class DirtyElectrodeDetector(D3Calculator):

    def calculate(self, data: dict, max_current_range: float = 0.00001):
        """
        Detect dirty electrode

        Connection Points
            :scan_data: optional, if e not provided, scan data will be used to find e (default = None)

        :param data: data for calculation
        :type data: dict
        :param max_current_range: maximum rango of current points allowed for a clean electrode
        :type max_current_range: float

        :return: average E1/2 (in units V)
        """
        self.data = data
        self.n = data.__len__()

        descriptor_cal = CVDescriptorCalculator(connector=self.key_pairs)
        peaks_dict = descriptor_cal.peaks(data)
        forward_peak = max([p[1] for p in peaks_dict.get("forward", [])])
        reverse_peak = min([p[1] for p in peaks_dict.get("reverse", [])])
        current_range = forward_peak - reverse_peak
        print("Dirty Electrode Detection CURRENT RANGE: ", current_range)

        return False if current_range < max_current_range else True


class CAResistanceCalculator(D3Calculator):

    def calculate(self, data: dict, offset_factor: float = 5, return_error: bool = False):
        """
        Calculator for calculating resistance after a CA experiment.

        Connection Points
            :i_s: list, current points (A)
            :t_s: list, time points (s)
            :pulse_width: float, time width of a pulse (s)
            :steps: int, number of potentiostat sweeps
            :low_e: float, lowest voltage (V)

        :param data: data for calculation
        :type data: dict
        :param offset_factor: factor by which to consider the voltage offset
        :type offset_factor: int
        :param return_error: return error and resistance if True
        :type return_error: bool

        :return: array [resistance, resistance error] if return_error is True, else resistance in Ohm
        """
        self.data = data
        self.n = data.__len__()

        conns = self.make_connections(data)
        pulse_width = unit_conversion(conns["pulse_width"], default_unit='s')
        low_e = unit_conversion(conns["low_e"], default_unit='V')

        dt = conns["t_s"][-1] / len(conns["t_s"])  # time interval per time measurement
        n = float(pulse_width / dt)  # number of time measurements in a pulse
        n_pulses = math.floor(conns["steps"] / 2) - 1  # number of pulses
        offset = n / offset_factor  # a slight offset to measure voltage after switching

        max_i = [max([abs(conns["i_s"][int(2 * (i + 1) * n - offset + j)]) for j in range(int(n))]) for i in
                 range(n_pulses)]

        last_R = abs(float(low_e) / np.mean(max_i))
        last_dR = last_R * np.std(max_i) / np.mean(max_i)

        return [last_R, last_dR] if return_error else last_R


# ------------------------- Molecular DFT Calculators --------------------------------


class EnergyDiffCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Difference between two energies calculation

        Connection Points:
            :energy_final: energy final (default = eV)
            :energy_initial: energy initial (default = eV)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: solvation energy in units eV
        """
        conns = self.make_connections(data)
        energy = (unit_conversion(conns["energy_final"], default_unit='eV') - unit_conversion(conns["energy_initial"],
                                                                                              default_unit='eV'))
        return float(np.format_float_scientific(energy, precision=precision))


class ReorganizationCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Reorganization energy calculation

        Connection Points:
            :gs_opt: ground state optimized energy (default = eV)
            :ion_opt: ion optimized energy (default = eV)
            :gs_energy: ground state energy at ion geometry (default = eV)
            :ion_energy: ion energy at ground state geometry (default = eV)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: reorganization energy (in units eV)
        """
        conns = self.make_connections(data)
        lambda1 = (unit_conversion(conns["ion_energy"], default_unit='eV') - unit_conversion(conns["ion_opt"],
                                                                                             default_unit='eV'))
        lambda2 = (unit_conversion(conns["gs_energy"], default_unit='eV') - unit_conversion(conns["gs_opt"],
                                                                                            default_unit='eV'))
        return float(np.format_float_scientific(lambda1 + lambda2, precision=precision))


class RelaxationCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Relaxation energy calculation from high energy geometry to optimized geometry

        Connection Points:
            :opt_energy: optimized energy (default = eV)
            :energy: energy another geometry (default = eV)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: relaxation energy (in units eV)
        """
        conns = self.make_connections(data)
        lambdax = unit_conversion(conns["energy"], default_unit='eV') - unit_conversion(conns["opt_energy"],
                                                                                        default_unit='eV')
        return float(np.format_float_scientific(lambdax, precision=precision))


class RMSDCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Root mean squared error calculator

        Connection Points:
            :geom_final: geometry final (default = A)
            :geom_initial: geometry initial (default = A)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: RMSD (in units A)
        """
        conns = self.make_connections(data)
        try:
            geom1 = pmgmol_to_rdmol(Molecule.from_sites([Site.from_dict(sd) for sd in conns["geom_initial"]]))[0]
            geom2 = pmgmol_to_rdmol(Molecule.from_sites([Site.from_dict(sd) for sd in conns["geom_final"]]))[0]
            print("Finding best RMS...this may take a few minutes...")
            with timeout(120, exception=RuntimeError):
                rmsd = rdMolAlign.GetBestRMS(geom1, geom2)
        except:
            raise ValueError("Error finding RMSD")
        return float(np.format_float_scientific(rmsd, precision=precision))


class DeltaGSolvCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Change in Gibbs energy for solvation calculator

        Connection Points:
            :init_eng: initial energy (default = eV)
            :init_corr: initial entropy correction (default = eV)
            :init_eng_solv: initial energy of solvation (default = eV)
            :fin_eng: final energy (default = eV)
            :fin_corr: final entropy correction (default = eV)
            :fin_eng_solv: final energy of solvation (default = eV)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: delta G solv  (in units A)
        """
        conns = self.make_connections(data)

        g_gas_init = unit_conversion(conns["init_eng"], default_unit='eV') + unit_conversion(conns["init_corr"],
                                                                                             default_unit='eV')
        g_gas_fin = unit_conversion(conns["fin_eng"], default_unit='eV') + unit_conversion(conns["fin_corr"],
                                                                                           default_unit='eV')

        # entropy correction cancels out because the species are the same
        delta_g_init_solv = unit_conversion(conns["init_eng_solv"], default_unit='eV') - unit_conversion(
            conns["init_eng"], default_unit='eV')
        delta_g_fin_solv = unit_conversion(conns["fin_eng_solv"], default_unit='eV') - unit_conversion(conns["fin_eng"],
                                                                                                       default_unit='eV')

        delta_g = g_gas_fin - g_gas_init + delta_g_fin_solv - delta_g_init_solv

        return float(np.format_float_scientific(delta_g, precision=precision))


class RedoxPotentialCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Redox potential calculator

        Connection Points:
            :init_eng: initial energy (default = eV)
            :init_corr: initial entropy correction (default = eV)
            :init_eng_solv: initial energy of solvation (default = eV)
            :fin_eng: final energy (default = eV)
            :fin_corr: final entropy correction (default = eV)
            :fin_eng_solv: final energy of solvation (default = eV)

            :num_electrons: number of electrons (default = 1)
            :electrode: electrode name as str or potential as float (default = standard_hydrogen_electrode)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: redox potential (in units A)
        """
        conns = self.make_connections(data)
        delta_g = DeltaGSolvCalc(connector=self.key_pairs).calculate(data)

        std_potential = get_electrode_potential(conns["electrode"]) if conns.get("electrode") else 4.42
        potential = -delta_g / conns.get("num_electrons", 1) + std_potential

        return float(np.format_float_scientific(potential, precision=precision))


class RadBuriedVolCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Radical buried volume for the atom with the atom with the highest portion
        of spin. Uses DBSTEP.

        Connection Points:
            :log_file: calculation output file. Must be readable with CCLIB
            :spin_type: type of CCLIB spin to extract (default = Mulliken)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: radical buried volume (in units A^3)
        """
        conns = self.make_connections(data)

        cmol = cclib.io.ccopen(conns["log_file"]).parse()
        spins = cmol.atomspins[conns.get("spin_type", "mulliken")]
        all_data = pd.DataFrame({"atoms": cmol.atomnos, "spin_density": spins})
        all_data["atom_sym"] = all_data.apply(lambda x: periodictable[int(x.atoms)], axis=1)
        all_data["atom_idx"] = all_data.index
        all_data["atom_idx"] = all_data.apply(lambda x: x.atom_idx + 1, axis=1)
        all_data = all_data[all_data.atom_sym != "H"]
        all_data["fractional_spin"] = all_data["spin_density"].abs() / all_data["spin_density"].abs().sum()

        self.max_cdf = all_data.loc[all_data["fractional_spin"].idxmax()]
        rad_bur_vol = float(db.dbstep(conns["log_file"], atom1=self.max_cdf["atom_idx"], volume=True).bur_vol)

        return float(np.format_float_scientific(rad_bur_vol, precision=precision))


class RadicalSpinCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Radical spin density. Uses DBSTEP.

        Connection Points:
            :log_file: calculation output file. Must be readable with CCLIB
            :spin_type: type of CCLIB spin to extract (default = Mulliken)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: radical stability score
        """
        self.make_connections(data)
        rad_bur_obj = RadBuriedVolCalc(connector=self.key_pairs)
        rad_bur_obj.calculate(data)
        max_cdf = rad_bur_obj.max_cdf

        return float(np.format_float_scientific(max_cdf["fractional_spin"], precision=precision))


class RSSCalc(D3Calculator):

    def calculate(self, data: dict, precision: int = 3):
        """
        Radical stability score. Uses DBSTEP.

        Connection Points:
            :log_file: calculation output file. Must be readable with CCLIB
            :spin_type: type of CCLIB spin to extract (default = Mulliken)

        :param data: data for calculation
        :param precision: number of significant figures (in scientific notation)
        :type data: dict
        :type precision: int

        :return: radical stability score
        """
        self.make_connections(data)
        rad_bur_obj = RadBuriedVolCalc(connector=self.key_pairs)
        rad_bur_vol = rad_bur_obj.calculate(data)
        max_cdf = rad_bur_obj.max_cdf
        rss = rad_bur_vol + 50 * (1 - max_cdf["fractional_spin"])

        return float(np.format_float_scientific(rss, precision=precision))
