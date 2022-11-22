from d3tales_api.Calculators.utils import *
from d3tales_api.Calculators.calculators_dft import *

class ConcentrationCalculator(D3Calculator):

    def calculate(self, data, precision=3):
        """
        Diffusion constatnt using Randles-Scidwick equation

        :connection points
        "smiles": "SMILES string",
        "volume": "Volume",
        "weight": "actual weight"

        "solv_density": "density for solvent"
        "redox_density": "density for redox-active molecule"

        :return: concentrations eV
        """
        self.data = data

        conns = self.make_connections(data)
        solv_density = conns.get("solv_density")
        redox_density = conns.get("redox_density")
        mol_weight = ExactMolWt(Chem.MolFromSmiles(conns["smiles"]))
        concentration = ((unit_conversion(conns["weight"], default_unit='g', density=redox_density) / unit_conversion(mol_weight, default_unit='g/mol')) /
                         unit_conversion(conns["volume"], default_unit='L', density=solv_density))
        return float(np.format_float_scientific(concentration, precision=precision))


class CVDescriptorCalculator(D3Calculator):

    def peaks(self, data, width=10):
        """
        CV peaks

        :connection points
        "scan_data" : scanned data from CV file

        :return: dictionary containing list of forward peaks and list of reverse peaks
        """

        self.data = data
        conns = self.make_connections(data)
        scan_dict = {}
        for data_list in conns["scan_data"]:
            data = np.array(data_list)
            if data[0, 0] < data[1, 0]:
                peaks_data = find_peaks(data[:, 1], width=width)
                scan_dict.update(
                    {"forward": self.prominent_peaks(peaks_data, data)})
            else:
                peaks_data = find_peaks(-data[:, 1], width=width)
                scan_dict.update(
                    {"reverse": self.prominent_peaks(peaks_data, data)})
        return scan_dict

    def reversibility(self, data, rev_upperbound=63, quasi_rev_upperbound=200, **kwargs):
        """
        Categorization of CV reversibility

        :rev_upperbound : upperbound for reversibility (mV)
        :quasi_rev_upperbound : upperbound for quasi reversibility (mV)
        connection points
        "scan_data" : scanned data from CV file (potentials in V)
        "sample_interval" : sample interval value (V)
        "low_e" : lowest energy value (V)

        :return: list of reversibility categorizations for peaks
        """

        self.data = data

        peaks = self.peaks(data, **kwargs)
        forward_peaks = sorted([item[0] for item in peaks.get('forward', [])])
        reverse_peaks = sorted([item[0] for item in peaks.get('reverse', [[]])], reverse=True)
        # If there are not the same number of forward and reverse peaks, default to irreversible.
        if len(forward_peaks) != len(reverse_peaks):
            return ['irreversible']
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

    def e_half(self, data, **kwargs):
        """
        CV E 1/2

        :connection points
        "scan_data" : scanned data from CV file (potentials in V)
        "sample_interval" : sample interval value (V)
        "low_e" : lowest energy value (V)

        :return: list of E 1/2 for peaks
        """

        self.data = data
        conns = self.make_connections(data)

        peaks = self.peaks(data, **kwargs)
        forward_peaks = sorted([item[0] for item in peaks.get('forward', [])])
        reverse_peaks = sorted([item[0] for item in peaks.get('reverse', [[]])], reverse=True)
        e_halfs = []
        for f_peak, r_peak in zip(forward_peaks, reverse_peaks):
            e_halfs.append(round((f_peak + r_peak) / 2, 3))
        return e_halfs

    def peak_splittings(self, data, **kwargs):
        """
        CV peak splitting

        :connection points
        "scan_data" : scanned data from CV file (potentials in V)
        "sample_interval" : sample interval value (V)
        "low_e" : lowest energy value (V)

        :return: list of peak splittings for peaks
        """

        self.data = data

        peaks = self.peaks(data, **kwargs)
        forward_peaks = sorted([item[0] for item in peaks.get('forward', [])])
        reverse_peaks = sorted([item[0] for item in peaks.get('reverse', [[]])], reverse=True)
        splittings = []
        for f_peak, r_peak in zip(forward_peaks, reverse_peaks):
            splittings.append(round(abs(f_peak - r_peak), 3))
        return splittings

    def middle_sweep(self, data):
        """
        CV middle sweep

        :connection points
        "scan_data" : scanned data from CV file (potentials in V)

        :return: middle sweep from the CV
        """

        self.data = data
        conns = self.make_connections(data)

        middle_idx = len(conns["scan_data"]) / 2
        return conns["scan_data"][int(middle_idx - 1):int(middle_idx + 1)]

    @staticmethod
    def prominent_peaks(peaks_data, orig_data, cutoff=0.0999):
        prominences = peaks_data[1]["prominences"]
        peaks_list = list(np.where(prominences / max(prominences) > cutoff)[0])
        values = []
        for idx in peaks_list:
            values.append(list(orig_data[[peaks_data[0][idx]], :][0]))
        return values


class CVDiffusionCalculator(D3Calculator):

    def calculate(self, data, precision=3, sci_notation=False):
        """
        Diffusion constant using Randles-Scidwick equation

        :connection points
        "i_p" : peak current (default = A)
        "A" : electrode area (default = cm^2)
        "v" : scan rate (default = V/s)
        "n" : number of electrons
        "C" : concentration of the solution (default = mol/cm^3)

        :return: average diffusion constant for single redox event, cm^2/s
        """
        self.data = data
        self.n = data.__len__()
        diffusion_constants = np.zeros(self.n)

        i_ps = np.zeros(self.n)
        vs = np.zeros(self.n)
        for idx, obj in enumerate(self.data):
            conns = self.make_connections(obj)
            print(conns)
            i_p = unit_conversion(conns["i_p"], default_unit='A')
            A = unit_conversion(conns["A"], default_unit='cm^2')
            v = unit_conversion(conns["v"], default_unit='V/s')
            C = unit_conversion(conns["C"], default_unit='mol/cm^3')
            i_ps[idx] = i_p
            vs[idx] = pow(v, 1/2)
            diffusion_constants[idx] = (pow(i_p / (2.692e5 * pow(conns["n"], 3 / 2) * A * C * pow(v, 1 / 2)), 2))
        slope = np.polyfit(vs, i_ps, 1)[0]
        # vs = vs[:, np.newaxis]
        # slope = np.linalg.lstsq(vs, i_ps, rcond=None)[0][0]
        diffusion_fitted = pow(slope / (2.692e5 * pow(conns["n"], 3 / 2) * A * C), 2)

        results = [np.format_float_scientific(diffusion_constants.mean(), precision=precision),
                   np.format_float_scientific(diffusion_fitted, precision=precision)]
        return results if sci_notation else [float(x) for x in results]


class CVChargeTransferCalculator(D3Calculator):

    def calculate(self, data, precision=3, sci_notation=False):
        """
        Charge tranfer rate calculation

        :connection points
        "T" : Temperature (default = 293 K)
        "D" : Diffusion constant (default = cm^2/s)
        "v" : scan rate (default = V/s)
        "n" : number of electrons
        "X" : peak shift (default = V)

        :return: array[["scan_rate", "k_0_i"]] for ith CV, cm/s
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
            long_terms.append(pow((math.pi * D * conns["n"] * 96485.34 * v) / (8.3145 * T), - 1 / 2))
        chargetranfer_rate = np.polyfit(long_terms, psis, 1)[0]

        result = np.format_float_scientific(chargetranfer_rate, precision=precision)
        return result if sci_notation else round(float(result), precision)


class AvgEHalfCalculator(D3Calculator):

    def calculate(self, data, precision=3, sci_notation=False):
        """
        Charge tranfer rate calculation

        :connection points
        "e" : E1/2 (default = V)

        :return: average E1/2 in units V
        """
        self.data = data
        self.n = data.__len__()

        e_halfs = []
        for idx, obj in enumerate(self.data):
            conns = self.make_connections(obj)
            E = unit_conversion(conns["e"], default_unit='V')
            e_halfs.append(E)
        avg_e_half = np.average(e_halfs)

        result = np.format_float_scientific(avg_e_half, precision=precision)
        return result if sci_notation else float(result)
