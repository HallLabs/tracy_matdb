from .vasp import AsyncVasp as Vasp
from .aflux import AsyncAflow as Aflow
from .qe import AsyncQe as Qe
from .tracy import Tracy_QE as TracyQE
from matdb import msg

try:
    from .quip import SyncQuip as Quip
except:
    msg.info("Could not import the Quip calculator.")


def get_calc_class(name):
    """Gets the class definition objects for the calculator that has the
    specified name.
    """
    globs = globals()
    try:
        target = globs[name]
    except KeyError:
        msg.err("Cannot import calculator {}. ".format(name) +
                "Does not exist at package level.")

    return target

def get_calculator_module(calcargs):
    """Returns the module corresponding to the calculator mentioned in
    `calcargs`.

    Args:
        calcargs (dict): the "calculator" dictionary that is part of the
          arguments for a db group.
    """
    from matdb import calculators
    from inspect import getmodule
    cls = getattr(calculators, calcargs["name"])
    mod = None

    try:
        mod = getmodule(cls)
    except: #pragma: no cover
        pass

    return mod
