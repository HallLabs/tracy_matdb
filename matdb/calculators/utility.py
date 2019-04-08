from six import string_types
from os import path
from hashlib import sha1
from operator import itemgetter

from matdb import msg
from matdb.utility import chdir, recursive_getattr
import matdb.calculators
    
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

    Note:
        atoms (matdb.atoms.Atoms): 
            default atoms object for the calculator. An empty `Atoms` object is created. 
            This shouldn't impact calculations since the calculator does not require a folder to run.

        workdir (str): 
            path to a working directory. Calculators that don't need this folder are the only ones 
            supported by this function, so we default to local directory.

        contr_dir (str): 
            path to the `matdb` controller's root directory. That folder is used only when 
            the calculator is hashed or serialized, for local use it doesn't matter so we default to the local directory.

        ran_seed (int): 
            random seed used to initialized the calculator.

    Raises:
        ValueError: if the `name` is not a folder-independent interatomic potential.
    """
    #This import is purposefully here to avoid recursive import loops.
    from matdb.atoms import Atoms
    target = matdb.calculators.get_calc_class(name)
    atoms = Atoms()
    if relpath is not None:
        with chdir(relpath):
            result = target(atoms, '.', '.', 0, *args, **kwargs)
    else:
        result = target(atoms, '.', '.', 0, *args, **kwargs)
    return result


def get_calculator_hashes(key, value, breadcrumb, result):
    """Looks for calculator absolute paths recursively in configuration dicts.
    
    .. warning:: The `result` parameter will be mutated.
    
    Args:
        key (str): key that this entry was under in the configuration dict.
        value: value object may be any supported type.
        breadcrumb (str): recursive path description for this key-value pair.
        result (dict): keys are breadcrumb hashes; values are absolute paths.
    """
    #Update the breadcrumb based on the context.
    if isinstance(key, int):
        breadcrumb += "[{0:d}]".format(key)
    elif isinstance(key, string_types):
        breadcrumb += "{0}{1}".format('.' if breadcrumb != "" else "", key)

    #We only recurse for list and dict types; other scalar types can't be
    #dictionary configs.
    if isinstance(value, list):
        for i, o in enumerate(value):
            get_calculator_hashes(i, o, breadcrumb, result)
    elif isinstance(value, dict):
        for k, v in sorted(value.items(), key=itemgetter(0)):
            get_calculator_hashes(k, v, breadcrumb, result)
        
    if (key == "calculator" and isinstance(value, dict)):
        if "name" in value:
            calcname = value["name"].lower()
            calccls = matdb.calculators.get_calc_class(calcname)            
        else:
            #Get the global calculator name/class instance.
            calcname = next(iter(result.keys())).lower()
            calccls = matdb.calculators.get_calc_class(calcname)
            
        for pathattr in calccls.pathattrs:
            attrval = recursive_getattr(value, pathattr)
            if attrval is not None:
                abspath = path.abspath(attrval)
                icrumb = "{0}.{1}".format(breadcrumb, pathattr)
                phash = str(sha1(icrumb.encode("ASCII")).hexdigest())
                #Add a new dictionary for the calculator if it is the first time
                if calcname not in result:
                    result[calcname] = {}
                result[calcname][phash] = abspath

paths = {}
"""dict: keys are hashed `matdb` names from multiple `matdb.yml` files. Values
are dictionaries populated by a call to :func:`get_calculator_hashes`.
"""

def set_paths(configyml):
    """Sets the absolute paths to any calculator-specific file resources based
    on the specified *live* configuration settings derived from a `matdb.yml`
    file.

    Args:
        configyml (dict): contents of the `matdb.yml` file to set paths for.
    """
    global paths
    name = configyml["title"].strip().replace(' ', '_')
    namehash = str(sha1(name.encode("ASCII")).hexdigest())
    calcpaths = {}
    get_calculator_hashes("matdb", configyml, "", calcpaths)
    paths[namehash] = calcpaths
