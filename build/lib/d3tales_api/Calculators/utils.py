import pint
import functools
from monty.serialization import loadfn
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


def dict2obj(d: dict, master_obj: object = None):
    """
    Convert a dictionary to a base object if `mater_obj` is specified, else generic object
    Adapted from: https://www.geeksforgeeks.org/convert-nested-python-dictionary-to-object/

    :param d: dictionary to be converted
    :param master_obj: base object to which to convert the dictionary
    :type d: dict
    :type master_obj: object
    :return: object from dictionary
    """
    # checking whether object d is an instance of class list
    if isinstance(d, list):
        d = [dict2obj(x) for x in d]

    # if d is not an instance of dict then directly object is returned
    if not isinstance(d, dict):
        return d

    # declaring a class
    class C(object):
        pass

    # constructor of the class passed to obj
    obj = master_obj if master_obj else C()

    for k in d:
        setattr(obj, k, dict2obj(d[k]))
        # obj.__dict__[k] = dict2obj(d[k])

    return obj


def json2obj(json_file, **kwargs):
    """
    Apply `dict2obj` to contents of a JSON file
    :param json_file: path to JSON file
    :type json_file: str
    :return: object from dictionary
    """
    d = loadfn(json_file)
    return dict2obj(d, **kwargs)


def unit_conversion(measurement, default_unit: str, density=None):
    """
    Convert a measurement into a default unit using pint. 
    
    :param measurement: Measurements can be pint object, int or float(in which case it will be assumed to already be in the default unit), string of magnitude and unit, or a measurement dictionary (EX: {"value": 0.5, "unit": "eV"}
    :param default_unit: default unit / unit to be converted to
    :param density: molecular density (in case needed for conversion)  
    :type default_unit: str
    :type density: str
    :return: float magnitude for the converted measurement 
    """
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
    unit = default_unit if unit == "dimensionless" else unit
    pint_unit = ureg("{}{}".format(value, unit))
    return pint_unit.to(default_unit, 'mol_density').magnitude


def get_electrode_potential(electrode):
    """
    Get electrode potential by searching the D3TaLES electrode parameters database
    :param electrode: name of an electrode or the electrode potential
    :return:
    """
    if str(electrode).replace(".", "").replace("-", "").isdigit():
        return float(electrode)
    params_db = ParamsDB(collection_name="electrode", schema_directory="materials")
    electrode_data = params_db.coll.find_one({"_id": electrode})
    abs_potential = electrode_data.get("absolute_potential")
    if abs_potential:
        return abs_potential.get("value")
    else:
        raise ValueError(f"Electrode {electrode} not found in the D3TaLES parameters database")


def get_periodic_table():
    """
    List elements in the periodic table
    :return: List of element abbreviations for elements in the periodic table 
    """
    return ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
            "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
            "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
            "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
            "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
            "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
            "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub", "Uut",
            "Uuq", "Uup", "Uuh", "Uus", "Uuo", ]


periodictable = get_periodic_table()
