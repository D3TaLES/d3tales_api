import abc
import json
import functools
import matplotlib.pyplot as plt


def rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)

    return functools.reduce(_getattr, [obj] + attr.split('.'))


class D3Plotter(abc.ABC):
    def __init__(self, connector=None):
        self.connector(key_pairs=connector)

    def connector(self, key_pairs=None):
        if key_pairs:
            self.key_pairs = key_pairs
            return self.key_pairs
        else:
            with open("connector.json") as f:
                connectors = json.load(f)
            self.key_pairs = connectors[self.__class__.__name__]

    def make_connections(self, obj):
        d = {}
        for key, connection in self.key_pairs.items():
            try:
                d.update({key: rgetattr(obj, connection)})
            except:
                d.update({key: obj[connection]})
        return d

    def plot_data(self, data):
        pass

    def description(self):
        print(self.__class__.__name__)


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
