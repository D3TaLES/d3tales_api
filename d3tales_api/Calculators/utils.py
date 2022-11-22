import abc
import math
import json
import pint
import functools
import numpy as np
from rdkit import Chem
from scipy.signal import find_peaks
from rdkit.Chem.Descriptors import ExactMolWt
from d3tales_api.D3database.d3database import ParamsDB


def rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)

    return functools.reduce(_getattr, [obj] + attr.split('.'))


def rgetkeys(_dict, keys, *args):
    def _getkey(_dict, key):
        _dict = _dict or {}
        if isinstance(_dict, dict):
            return _dict.get(key, *args)
        if isinstance(_dict, list) and key.isdigit():
            return _dict[int(key)]

    return functools.reduce(_getkey, [_dict] + keys.split('.'))


def unit_conversion(measurement, default_unit, density=None):
    if not measurement:
        return None
    # Set context in case conversion include mass-->volume or volume-->mass
    ureg = pint.UnitRegistry()
    c = pint.Context('mol_density')
    if density:
        c.add_transformation('[mass]', '[volume]', lambda ureg_c, x: x / ureg_c(density))
        c.add_transformation('[volume]', '[mass]', lambda ureg_c, x: x * ureg_c(density))
    ureg.add_context(c)
    # Get measurement value and unit
    if not isinstance(measurement, (str, float, int, dict)):
        value, unit = getattr(measurement, "value"), getattr(measurement, "unit")
    else:
        value = measurement.get("value") if isinstance(measurement, dict) else measurement
        unit = ""
        if isinstance(value, float) or str(value).replace('.', '', 1).replace('-', '', 1).isdigit():
            unit = measurement.get("unit", default_unit) if isinstance(measurement, dict) else default_unit
    # Convert measurement to default unit
    pint_unit = ureg("{}{}".format(value, unit))
    return pint_unit.to(default_unit, 'mol_density').magnitude

def get_electrode_potential(electrode):
    if str(electrode).isdigit():
        return fload(electrode)
    params_db = ParamsDB(collection_name="electrode", schema_directory="materials")
    electrode_data = params_db.coll.find_one({"_id": electrode})
    abs_potential = electrode_data.get("absolute_potential")
    if abs_potential:
        return abs_potential.get("value")
    else:
        raise ValueError(f"Electrode {electrode} note found in the D3TaLES prameters database")

class D3Calculator(abc.ABC):
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
                d.update({key: rgetkeys(obj, connection)})
        return d

    def calculate(self, data):
        pass

    def description(self):
        print(self.__class__.__name__)


class D3Plotter(D3Calculator):

    def plot_data(self, data):
        pass
