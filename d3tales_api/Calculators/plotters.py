from scipy.stats import norm
import scipy.constants as cst
import matplotlib.pyplot as plt
from d3tales_api.Calculators.utils import *


class CVPlotter(D3Plotter):

    def plot_data(self, data):
        """
        CV plot data for plotly

        :connection points
        "scan_data" : scanned data from CV file

        :return: plot data for plotly
        """

        self.data = data
        conns = self.make_connections(data)

        x = []
        y = []
        for scan in conns["scan_data"]:
            x.extend([i[0] for i in scan])
            y.extend([i[1] for i in scan])

        plotting_data = [{
            "x": x,
            "y": y,
            "modeplot_data": 'lines',
            "name": "cv",
            "line": {
                "color": "#003396",
                "width": 3
            }
        }]
        return plotting_data

    def live_plot(self, data):
        self.data = data
        conns = self.make_connections(data)
        scan_data = conns["scan_data"]
        for data in scan_data:
            plt.plot(data[:, 0], data[:, 1], color="red")

        plt.show()


class DFTSpecPlotter(D3Plotter):

    def plot_data(self, data):
        """
        Spectrum plot data for plotly

        :connection points
        "transitions" : A list of tuple for each transition such as
                        [(energy (eV), lambda (nm), oscillatory strength), ... ]
        "sigma" : (default = 0.10)
        "step" : (default = 0.01)

        :return: plot data for plotly
        """

        self.data = data
        conns = self.make_connections(data)
        transitions = conns["transitions"]
        sigma = conns["sigma"]
        step = conns["step"]

        negative_absorptions = [transitions.pop(index) for index, val in enumerate(transitions) if val[0] < 0]
        if negative_absorptions:
            print("WARNING: {} calculation contained a negative excitation energy".format(self.calculation_type))

        minval = min([val[0] for val in transitions]) - 5.0 * sigma
        maxval = max([val[0] for val in transitions]) + 5.0 * sigma
        npts = int((maxval - minval) / step) + 1

        eneval = np.linspace(minval, maxval, npts)  # in eV
        lambdaval = [cst.h * cst.c / (val * cst.e) * 1.e9
                     for val in eneval]  # in nm

        # sum of gaussian functions
        spectre = np.zeros(npts)
        for trans in transitions:
            spectre += trans[2] * norm.pdf(eneval, trans[0], sigma)
        spectre /= spectre.max()
        spectra_data = {"energies": list(eneval), "lambda": lambdaval, "xas": list(spectre)}

        # for plotting in frontend
        abs_plot = [{"x": [round(x, 0) for x in spectra_data["lambda"]],
                     "y": spectra_data["xas"],
                     "mode": 'lines',
                     "name": "absorption",
                     "line": {
                         "color": "#003396",
                         "width": 3
                     }}]
        spectra_data.update({"abs_plot": abs_plot})
        return spectra_data

    def live_plot(self, data):
        spectra_data = self.calculate(data)
        plt.plot(spectra_data["lambda"], spectra_data["xas"], color="red")

        plt.show()

