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
