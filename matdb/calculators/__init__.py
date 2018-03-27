from .vasp import AsyncVasp as Vasp
from .quip import SyncQuip as Quip
from .aflux import AsyncAflow as Aflow

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
