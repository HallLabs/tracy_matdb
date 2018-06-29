from .vasp import AsyncVasp as Vasp
from .aflux import AsyncAflow as Aflow
from .qe import AsyncQe as Qe
from matdb import msg
from matdb.utility import chdir

try:
    from .quip import SyncQuip as Quip
except:
    msg.info("Could not import the Quip calculator.")

def build_calc(name, relpath, *args, **kwargs):
    """Builds a calculator instance using sensible defaults for *interatomic potentials*
    that do *not* require a temporary directory to dump files.

    .. note:: The `Vasp`, `Aflow` and `Qe` calculators are not buildable using this
      function.

    .. warning:: The calculator produced by this function *cannot* be hashed or serialized
      as part of a full `matdb`; it is intended for _local_ use only.

    Args:
        name (str): name of the calculator in this package; one of ['Quip'].
        relpath (str): path to the directory in which to instantiate the calculator.

    Notes:
        atoms (matdb.Atoms): default atoms object for the calculator. An
          empty `Atoms` object is created. This shouldn't impact calculations since the
          calculator does not require a folder to run.
        workdir (str): path to a working directory. Calculators that don't need this
          folder are the only ones supported by this function, so we default to local
          directory.
        contr_dir (str): path to the `matdb` controller's root directory. That folder is
          used only when the calculator is hashed or serialized, for local use it doesn't
          matter so we default to the local directory.
        ran_seed (int): random seed used to initialized the calculator.

    Raises:

    ValueError: if the `name` is not a folder-independent interatomic potential.
    """
    globs = globals()
    try:
        target = globs[name]
    except KeyError:
        msg.err("Cannot import calculator {}. ".format(name) +
                "Does not exist at package level.")

    from matdb.atoms import Atoms
        atoms = Atoms()
    if relpath is not None:
        with chdir(relpath):
            result = target(atoms, '.', '.', 0, *args, **kwargs)
    else:
        result = target(atoms, '.', '.', 0, *args, **kwargs)
    return result

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
