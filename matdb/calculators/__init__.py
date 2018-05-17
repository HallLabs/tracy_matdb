from .vasp import AsyncVasp as Vasp
from .aflux import AsyncAflow as Aflow
from .qe import AsyncQe as Qe
from matdb.msg import info
try:
    from .quip import SyncQuip as Quip
except:
    info("Could not import the Quip calculator.")

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
