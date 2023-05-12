import numpy as np
from scipy.stats import norm
import scipy.constants as cst
import matplotlib.pyplot as plt
from collections import OrderedDict
from d3tales_api.Calculators.calculators import *


class D3Plotter(D3Calculator):
    """    
    D3Plotters base class, based on D3Calculators base class
    """

    def plot_data(self, data):
        pass


class CVPlotter(D3Plotter):

    def plot_data(self, data, self_standard=False):
        """
        CV plot data for plotly

        Connection Points:
            :scan_data: scanned data from CV file

        :param data: data for calculation
        :type data: dict
        :param self_standard: establish self standard (e_half=0V) if True
        :type self_standard: bool

        :return: plot data for plotly
        """

        self.data = data
        conns = self.make_connections(data)

        x = []
        y = []
        for scan in conns["scan_data"]:
            x.extend([i[0] for i in scan])
            y.extend([i[1] for i in scan])

        if self_standard:
            print("Performing self-standard adjustment...")
            e_half = CVDescriptorCalculator(self.key_pairs).e_half(data)[0]
            x = [i-e_half for i in x]

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

    def live_plot(self, data, fig_path=None, self_standard=False, **plt_kwargs):
        """
        Live Matplotlib plot for data

        Connection Points:
            :scan_data: scanned data from CV file

        :param data: data for calculation
        :type data: dict
        :param fig_path: path to which to save the figure
        :type fig_path: str
        :param self_standard: establish self standard (e_half=0V) if True
        :type self_standard: bool

        :return: shows matplotlib plot
        """
        plt_data = self.plot_data(data, self_standard=self_standard)
        plt.scatter(plt_data["x"], plt_data["y"], color="red", s=10)
        plt.gca().update(dict(**plt_kwargs))
        plt.tight_layout()

        if fig_path:
            plt.savefig(fig_path, dpi=300)
            plt.close()
        else:
            plt.show()

    def live_plot_multi(self, data, fig_path=None, sort=True, self_standard=False, **plt_kwargs):
        """
        Live Matplotlib plot for data

        Connection Points:
            :scan_data: scanned data from CV file
            :variable_prop: property that varies between CVs

        :param data: data for calculation
        :type data: list
        :param fig_path: path to which to save the figure
        :type fig_path: str
        :param sort: sort by variable_prop if True
        :type sort: bool
        :param self_standard: establish self standard (e_half=0V) if True
        :type self_standard: bool

        :return: shows matplotlib plot
        """
        # Get data_dict for sorting
        data_dict = {}
        for d in data:
            plt_data = self.plot_data(d, self_standard=self_standard)
            conns = self.make_connections(d)
            var_prop = conns["variable_prop"]
            data_dict[var_prop] = plt_data

        if sort:
            data_dict = dict(OrderedDict(sorted(data_dict.items())))

        [plt.scatter(d["x"], d["y"], label=p, s=10) for p, d in data_dict.items()]
        legend_title = plt_kwargs.pop("legend_title") if "legend_title" in plt_kwargs.keys() else "Legend"
        plt.gca().update(dict(**plt_kwargs))
        plt.tight_layout()
        plt.legend(title=legend_title)

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

    def live_plot(self, data, fig_path=None, **plt_kwargs):
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
        plt.gca().update(dict(**plt_kwargs))
        plt.tight_layout()

        if fig_path:
            plt.savefig(fig_path, dpi=300)
        else:
            plt.show()

