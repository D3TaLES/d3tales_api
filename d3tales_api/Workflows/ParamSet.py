import warnings

from monty.json import MSONable
from monty.serialization import loadfn


def sett(obj, cls):
    if obj is None or obj == {}:
        return cls.__class__()
    elif isinstance(obj, dict):
        return cls.__class__.from_dict(obj)
    elif isinstance(obj, cls.__class__):
        return obj
    else:
        raise TypeError(type(obj))


def _set_calculation_parameters(instance, parameters, parameter_type):
    for calculation, values in parameters.items():
        if (not isinstance(calculation, str) or not calculation.isidentifier()
                or hasattr(instance, calculation)):
            raise ValueError("Invalid or reserved calculation name: {!r}".format(calculation))
        setattr(instance, calculation, sett(values, parameter_type()))


class GaussianParameters(MSONable):
    def __init__(self, route_parameters=None, charge=None, multiplicity=None, basis_set=None,
                 functional=None, input_parameters=None, link0_parameters=None):
        self.route_parameters = {} if route_parameters is None else route_parameters
        self.input_parameters = {} if input_parameters is None else input_parameters
        self.link0_parameters = {} if link0_parameters is None else link0_parameters
        self.charge = 0 if charge is None else charge
        self.multiplicity = 1 if multiplicity is None else multiplicity
        self.basis_set = basis_set
        self.functional = functional


class GausParamSet(MSONable):
    def __init__(self, name=None, **kwargs):
        self.name = name or ""
        _set_calculation_parameters(self, kwargs, GaussianParameters)

    def as_dict(self):
        d = {"@module": type(self).__module__, "@class": type(self).__name__, "name": self.name}
        d.update({k: v.as_dict() for k, v in vars(self).items() if isinstance(v, GaussianParameters)})
        return d

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
        warnings.warn(
            "VaspParamSet is not supported by the current workflow registry and will be removed "
            "in a future release.",
            DeprecationWarning,
            stacklevel=2,
        )
        self.name = name or ""
        _set_calculation_parameters(self, kwargs, VaspParameters)

    @staticmethod
    def from_json(jsonfile):
        d = loadfn(jsonfile)
        return VaspParamSet.from_dict(d)
