# https://www.geeksforgeeks.org/convert-nested-python-dictionary-to-object/
from monty.serialization import loadfn


def dict2obj(d: dict, master_obj: object = None):
    """
    Convert a dictionary to a base object if `mater_obj` is specified, else generic object

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
