from monty.json import MSONable
from monty.serialization import loadfn


def sett(obj, cls):
    if not obj:
        return cls.__class__()
    elif isinstance(obj, dict):
        return cls.__class__.from_dict(obj)
    elif isinstance(obj, cls.__class__):
        return obj
    else:
        raise TypeError(type(obj))


class GaussianParameters(MSONable):
    def __init__(self, route_parameters=None, charge=None, multiplicity=None, basis_set=None,
                 functional=None, input_parameters=None, link0_parameters=None):
        self.route_parameters = route_parameters or {}
        self.input_parameters = input_parameters or {}
        self.link0_parameters = link0_parameters or {}
        self.charge = charge or 0
        self.multiplicity = multiplicity or 1
        self.basis_set = basis_set
        self.functional = functional


class GausParamSet(MSONable):
    def __init__(self, name=None, **kwargs):
        self.name = name or ""
        for calc, params in kwargs.items():
            exec("self.{} = sett(params, GaussianParameters())".format(calc))

    @staticmethod
    def from_json(jsonfile):
        d = loadfn(jsonfile)
        return GausParamSet.from_dict(d)


class VaspParameters(MSONable):
    def __init__(self, incar=None, kpoints=None):
        self.incar = incar or {}
        self.kpoints = kpoints or {}


class VaspParamSet(MSONable):
    def __init__(self, name=None, **kwargs):
        self.name = name or ""
        for calc, params in kwargs.items():
            exec("self.{} = sett(params, VaspParameters())".format(calc))

    @staticmethod
    def from_json(jsonfile):
        d = loadfn(jsonfile)
        return VaspParamSet.from_dict(d)
