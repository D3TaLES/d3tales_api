import numpy as np
from scipy.stats import norm
import scipy.constants as cst
import matplotlib.pyplot as plt
from d3tales_api.Calculators.calculators import D3Calculator


class D3Plotter(D3Calculator):
    """    
    D3Plotters base class, based on D3Calculators base class
    """

    def plot_data(self, data):
        pass


class CVPlotter(D3Plotter):

    def plot_data(self, data):
        """
        CV plot data for plotly

        Connection Points:
            :scan_data: scanned data from CV file

        :param data: data for calculation
        :type data: dict

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
        return {"abs_plot": plotting_data, "x": x, "y": y}

    def live_plot(self, data, fig_path=None):
        """
        Live Matplotlib plot for data

        Connection Points:
            :scan_data: scanned data from CV file

        :param data: data for calculation
        :type data: dict
        :param fig_path: path to which to save the figure
        :type fig_path: str

        :return: shows matplotlib plot
        """
        plt_data = self.plot_data(data)
        plt.scatter(plt_data["x"], plt_data["y"], color="red")

        if fig_path:
            plt.savefig(fig_path, dpi=300)
            plt.close()
        else:
            plt.show()


class DFTSpecPlotter(D3Plotter):

    def plot_data(self, data):
        """
        Spectrum plot data for plotly

        Connection Points:
            :transitions: A list of tuple for each transition such as
                            [(energy (eV), lambda (nm), oscillatory strength), ... ]
            :sigma: (default = 0.10)
            :step: (default = 0.01)

        :param data: data for calculation
        :type data: dict

        :return: plot data for plotly
        """

        self.data = data
        conns = self.make_connections(data)
        transitions = conns["transitions"]
        sigma = conns["sigma"]
        step = conns["step"]

        negative_absorptions = [transitions.pop(index) for index, val in enumerate(transitions) if val[0] < 0]
        if negative_absorptions:
            print("WARNING: Calculation contains a negative excitation energy")

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

    def live_plot(self, data, fig_path=None):
        """
        Live Matplotlib plot for data

        Connection Points:
            :transitions: A list of tuple for each transition such as
                            [(energy (eV), lambda (nm), oscillatory strength), ... ]
            :sigma: (default = 0.10)
            :step: (default = 0.01)

        :param data: data for calculation
        :type data: dict
        :param fig_path: path to which to save the figure
        :type fig_path: str


        :return: shows matplotlib plot
        """
        spectra_data = self.calculate(data)
        plt.plot(spectra_data["lambda"], spectra_data["xas"], color="red")

        if fig_path:
            plt.savefig(fig_path, dpi=300)
        else:
            plt.show()

