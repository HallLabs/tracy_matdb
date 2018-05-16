"""Utility functions for interacting with file system, shell, etc.
"""
from os import path
from six import string_types
from matdb import msg
import six
import numpy as np
import h5py

import sys
from contextlib import contextmanager
@contextmanager
def redirect_stdout(new_target):
    old_target, sys.stdout = sys.stdout, new_target # replace sys.stdout
    try:
        yield new_target # run some code with the replaced stdout
    finally:
        sys.stdout = old_target # restore to the previous value

@contextmanager
def chdir(target):
    """Context manager for executing some code within a different
    directory after which the current working directory will be set
    back to what it was before.

    Args:
        target (str): path to the directory to change into.
    """
    from os import getcwd, chdir
    current = getcwd()
    try:
        chdir(target)
        yield target
    finally:
        chdir(current)

def import_fqdn(fqdn):
    """Returns the object from the specified FQDN. Any exceptions raised will
    bubble up.

    Args:
        fqdn (str): '.'-separated list of `package.module.callable` to
          import. The callable will *not* be called.

    Returns:
        tuple: `(module, callable)`, where `module` is the module object that
        `callable` resides in.
    """
    from importlib import import_module
    parts = fqdn.split('.')
    call = parts[-1]
    module = '.'.join(parts[0:-1])
    module = import_module(module)
    call = getattr(module, call)
    return (module, call)

def execute(args, folder, wait=True, nlines=100, venv=None,
            printerr=True, env_vars=None, errignore=None, **kwargs):
    """Executes the specified tuple that should include the command as
    first item and additional arguments afterward. See the
    documentation for :class:`subprocess.Popen` for details.

    Args:
        args (list): of `str`; first item should be command to
          execute; additional arguments following.
        folder (str): directory to switch into before executing the
          command.
        wait (bool): when True, block the current thread until
          execution completes; otherwise, returns immediately.
        nlines (int): by default, `stdout` and `stderr` are redirected to
          :data:`subprocess.PIPE`. This is the maximum number of lines that will
          be returned for large outputs (so that memory doesn't get overwhelmed
          by large outputs).
        venv (str): when not `None`, the name of a virtualenv to
          activate before running the command.
        printerr (bool): when True, if `stderr` is not empty, print
          the lines automatically.
        env_vars (dict): of environment variables to set before calling the
          execution. The variables will be set back after execution.
        errignore (str): if errors are produced that include this pattern, then
          they will *not* be printed to `stdout`.
        kwargs (dict): additional arguments that are passed directly
          to the :class:`subprocess.Popen` constructor.

    Returns:
        dict: with keys ['process', 'stdout', 'stderr'], where 'process' is the
        instance of the subprocess that was created; 'stdout' and 'stderr' are
        only included if they were set to :data:`subprocess.PIPE`.

    .. note:: If the output from 'stdout' and 'stderr' are too large, only the
      first 100 lines will be returned. Use parameter `nlines` to control output
      size.
    """
    from subprocess import Popen, PIPE
    if "stdout" not in kwargs:
        kwargs["stdout"] = PIPE
    if "stderr" not in kwargs:
        kwargs["stderr"] = PIPE
    kwargs["cwd"] = folder

    if venv is not None: # pragma: no cover No guarantee that virtual
                         # envs exist on testing machine.
        if isinstance(venv, string_types): 
            vargs = ["virtualenvwrapper_derive_workon_home"]
            vres = execute(vargs, path.abspath("."))
            prefix = path.join(vres["output"][0].strip(), venv, "bin")
        elif venv == True:
            import sys
            prefix = path.dirname(sys.executable)
        args[0] = path.join(prefix, args[0])

    from os import environ
    if env_vars is not None:
        oldvars = {}
        for name, val in env_vars.items():
            oldvars[name] = environ[name] if name in environ else None
            environ[name] = val
        
    msg.std("Executing `{}` in {}.".format(' '.join(args), folder), 2)
    pexec = Popen(' '.join(args), shell=True, executable="/bin/bash", **kwargs)
    
    if wait:
        from os import waitpid
        waitpid(pexec.pid, 0)

    if env_vars is not None:
        #Set the environment variables back to what they used to be.
        for name, val in oldvars.items():
            if val is None:
                del environ[name]
            else:
                environ[name] = val
        
    #Redirect the output and errors so that we don't pollute stdout.
    output = None
    if kwargs["stdout"] is PIPE:
        output = []
        for line in pexec.stdout:
            output.append(line)
            if len(output) >= nlines:
                break
        pexec.stdout.close()

    error = None
    if kwargs["stderr"] is PIPE:
        error = []
        for line in pexec.stderr:
            if errignore is None or errignore not in line:
                error.append(line)
            if len(error) >= nlines:
                break
        pexec.stderr.close()
        if printerr and len(error) > 0:
            msg.err(''.join(error))

    return {
        "process": pexec,
        "output": output,
        "error": error
    }    

def cat(files, target):
    """Combines the specified list of files into a single file.

    Args:
        files (list): of `str` file paths to combine.
        target (str): name/path of the output file that will include all of the
          combined files.
    """
    with open(target, 'w') as outfile:
        for fname in files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def symlink(target, source):
    """Creates a symbolic link from `source` to `target`.
    """
    from os import path, symlink, remove
    from matdb import msg
    if path.isfile(target) or path.islink(target):
        remove(target)
    elif path.isdir(target):
        msg.warn("Cannot auto-delete directory '{}' for symlinking.".format(target))
        return
    
    symlink(source, target)

def linecount(filename):
    """Counts the number of lines in file.

    Args:
        filename (str): full path to the file to count lines for.
    """
    if not path.isfile(filename):
        return 0
    
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def dict_update(a, b):
    """Inserts keys and values from `b` into `a` if they don't already
    exist. This function performs this updated recursively for nested `dict`.

    .. warning:: This function modifies the variable `a`.
    """
    for k, v in b.items():
        if k not in a:
            a[k] = v
        elif isinstance(a[k], dict):
            dict_update(a[k], v)

def safe_update(obj, kv):
    """Updates the key-value pairs on `obj` only if they are not
    `None`.

    Args:
        obj: attributes will be set on this object according to the
          keys in `kv`.
        kv (dict): keys are target attribute names; values are desired
          values.
    """
    for k, v in kv.items():
        if hasattr(obj, k) and getattr(obj, k) is not None:
            continue
        setattr(obj, k, v)

def obj_update(obj, k, v, copy=True):
    """Updates a particular value in a set of nested objects.

    .. note:: If `copy=True`, a *copy* of the object will be
      manipulated, though it will only use `copy()`, not any deep
      copying.

    Args:
        obj: list or list of dict to update some sub-key on.
        k (str): a '.'-separated list of keys to traverse down to get
          to the value to set.
        v: value to set at the final key.
        copy (bool): when True, return a copy of the dictionary.
    """
    from copy import copy as ocopy
    chain = list(reversed(k.split('.')))
    result = ocopy(obj) if copy else obj
    target = obj

    if isinstance(target, list):
        firstkey = chain[-1]
        target = next(d for d in target if firstkey in d)
    
    while len(chain) > 1:
        key = chain.pop()
        target = (target[key] if isinstance(target, dict)
                  else getattr(target, chain.pop()))

    if isinstance(target, dict):
        target[chain[0]] = v
    else:
        setattr(target, chain[0], v)
        
    return result
        
def _get_reporoot():
    """Returns the absolute path to the repo root directory on the current
    system.
    """
    import matdb
    medpath = path.abspath(matdb.__file__)
    return path.dirname(path.dirname(medpath))

def relpath(s):
    """Returns the *repository-relative* path for the string `s`.

    Args:
        s (str): repository-relative path, see the examples.

    Examples:

        Suppose I have a repository at `/usr/opt/repo`, then:

        >>> relpath("./tests") == "/usr/opt/repo/tests"
        True
        >>> relpath("../other/docs") == "/usr/opt/other/docs"
        True
    """
    with chdir(reporoot):
        result = path.abspath(s)
    return result

def copyonce(src, dst):
    """Copies the specified file to the target *only* if it doesn't
    already exist.
    """
    if not path.isfile(dst):
        from shutil import copyfile
        copyfile(src, dst)

def compare_tree(folder, model):
    """Compares to directory trees to determine if they have the same
    contents. This does *not* compare the contents of each file, but rather just
    that the file exists.

    Args:
        folder (str): path to the root directory that has folders/files in the
          first level of `model`.
        model (dict): keys are either folder names or `__files__`. For folder
          names, :func:`compare_tree` is called recursively on the next level
          down. For files, the existing of each file is checked.
    """
    if "__files__" in model:
        for fname in model["__files__"]:
            assert path.isfile(path.join(folder, fname))

    for foldname, tree in model.items():
        if foldname == "__files__":# pragma: no cover
            continue

        target = path.join(folder, foldname)
        assert path.isdir(target)
        compare_tree(target, tree)

def which(program):
    """Tests whether the specified program is anywhere in the environment
    PATH so that it probably exists."""
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

def touch(fpath):
    """Mimics the `touch` command in the unix to create an empty file with the
    given path.
    """
    import os
    with open(fpath, 'a'):
        os.utime(fpath, None)

def pgrid(options, ignore=None):
    """Creates a parameter grid over the specified options.

    .. note:: This function treats keys that end in `*` specially. When the key
      ends in `*`, the values in the specified list contribute to the cartesian
      product.

    Args:
        options (dict): key-value pairs to iterate across.
        ignore (list): of `str` keys to ignore in the options dict.

    Returns:
        tuple: `(grid, keys)` where grid is a list of tuples with a value for
        each key in `options` and keys is a listed of the keys in the order they
        appear in the tuples.
    """
    from itertools import product
    from collections import OrderedDict
    from operator import itemgetter

    #Sort the dict so that we get the same ordering of parameters in the
    #product.
    params = OrderedDict(sorted(options.items(), key=itemgetter(0)))
    values = []
    keys = []
    for k, v in params.items():
        if ignore is None or k in ignore:
            continue
        
        if k[-1] != '*':
            values.append([v])
        else:
            values.append(v)
        keys.append(k.strip('*'))
        
    grid = list(product(*values))
    return (grid, keys)

import pytz
from datetime import datetime
from six import string_types

epoch = datetime(1970,1,1, tzinfo=pytz.utc)
"""datetime.datetime: 1/1/1970 for encoding UNIX timestamps.
"""

def datetime_handler(x):
    """Prepares a :class:`datetime.datetime` for JSON serialization by turning
    it into the ISO format for time-zone sensitive data.

    Args:
        x (datetime.datetime): to be stringified.

    Returns:
        str: ISO format; `None` if `x` is not a date.
    """
    if isinstance(x, datetime):
        return x.isoformat()

def parse_date(v):
    """Parses the date from the specified single value or list of values.

    Args:
        v (str): string representation of the :class:`datetime` returned by
          :func:`datetime_handler`.
    """
    from dateutil import parser
    if isinstance(v, (list, tuple)):
        return [parse_date(vi) for vi in v]
    elif isinstance(v, string_types):
        return parser.parse(v)
    else:
        raise ValueError("Not a valid datetime string.")

def load_datetime(pairs):
    """Deserialize a JSON string with dates encoded by
    :func:`datetime_handler`.
    """
    d = {}
    for k, v in pairs:
        if isinstance(v, (list, tuple, string_types)):
            try:
                d[k] = parse_date(v)
            except ValueError:# pragma: no cover
                d[k] = v
        else:
            d[k] = v             
    return d

def getattrs(obj, chain):
    """Recursively gets attributes in a chain from object to object until the
    train terminates.

    Args:
        obj: to get attributes from.
        chain (str): of '.'-separated attribute names.

    Examples:

        Get the `energy` attribute of the `atoms` attribute of a database
          (`obj`).

        >>> getattrs(obj, "atoms.energy")
    """
    o = obj
    for attr in chain.split('.'):
        if isinstance(o, dict):
            o = o[attr]
        else:
            o = getattr(o, attr)
    return o
        
reporoot = _get_reporoot()
"""The absolute path to the repo root on the local machine.
"""

def slicer(obj, args):
    """Slices the object along the path defined in args.
    
    Args:
        obj (iterable): an object to be sliced or divided.
        args (iterable): the locations that the slices should be at.
    """
    from itertools import islice
    if not isinstance(args,(list,tuple)):
        msg.err("The slicer args must be a list or a tuple.")
        return
    result = []
    if len(args)%2==0:
        for a in range(0,len(args),2):
            for b, c in islice(((args[a],args[a+1]),),0,None,2):
                result.extend(obj[slice(b,c)])
    else:
        msg.err("Could not parse slice {} without start and stop values.".format(args))
    return result

def _py_execute(module, function, params):
    """Executes the specified function by importing it and then using `eval` on
    the exist parameter string given.

    Args:
        module (str): FQDN of the module in which `function` resides.
        function (str): name of the callable in `module` to execute.
        params (str): exact parameter string (including parentheses) to pass to
          the function call with `eval`.
    """
    from importlib import import_module
    module = import_module(module)
    call = getattr(module, function)

    xstr = "call{}".format(params)
    return eval(xstr)

def special_values(vs, seed=None):
    """Converts the specified special value string into its python
    representation. We allow the following "special" directives for
    parameter values:

    1. "linspace(start, stop, length)" calls the :func:`numpy.linspace` with the
       given parameters to produce the list of values for the grid search.
    2. "logspace(start, stop, length)" calls :func:`numpy.logspace` to produce
       the parameter list.
    3. "range(start, stop, step)" creates a range of numbers using the built-in
       python `range` iterator. It will have type `list`.
    4. "random:{id}(\*\*params)" draws samples from one of the distributions on
       :class:`numpy.random.RandomState`. `id` should be one of the instance
       methods for that distribution and `\*\*params` should be the exact python
       code you would pass to the method call.
    5. "distr:{id}(\*\*params)" draws samples from a discrete/continuous
       distribution in `scipy.stats`. `id` should be one of the distribution
       names in that module (which has a `rvs` method); `\*\*params` should be the
       exact python code would pass to the object constructor.

    Args:
        vs (str): special value string.

    """
    if vs is None or not isinstance(vs, string_types):
        return vs
    
    sdict = {
        "linspace": ("numpy", "linspace"),
        "logspace": ("numpy", "logspace"),
        "range": ("numpy", "arange"),
        "random:" : ("numpy.random", None),
        "distr:" : ("scipy.stats", None),
        "[": ("slicer",None)
    }

    #We allow for |nogs| to be appended if the person is just specifying weights
    #for each input value (for example).
    v = vs.replace("|nogs|", "")
    
    for k, f in sdict.items():
        if k in v:
            module, function = f
            rest = v[len(k):]
            if function is not None:
                result = _py_execute(module, function, rest)
            else:
                #We still need to determine the name of the callable.
                if k in ["random:","distr:"]:
                    first = rest.index('(')
                    caller = rest[:first]
                if k == "random:":
                    from numpy.random import RandomState
                    rs = RandomState(seed)
                    d = getattr(rs, caller)
                    result = eval("d{}".format(rest[first:]))
                elif k == "distr:":
                    result = _py_execute(module, caller, rest[first:])
                elif k== "[":
                    temp = eval(v)
                    result = slicer(range(1,max(temp)),temp)
            break
    else:
        result = v    
        
    return result

import collections

def special_functions(sf,values):
    """Converts the specified function value string into its python
    representation and evaluates it for each item in the values list. 
    We allow the following "special" directives for parameter values:

    1. "linalg:{id}" 'id' is an operation to be performed on a matrix 
        from the numpy.linalg package.
    2. "math:{id}" 'id' is an operation from the math package.
    3. "numpy:{id}" 'id' is any operation from numpy.

    Args:
        sf (str): special function string.
        values (list): a list that the function is to be applied to.

    .. note:: the value returned by the special function must be an integer or a float.
    """
    import numpy as np
    import math
    if sf is None or not isinstance(sf, (string_types,dict)):
        raise ValueError("The special function must be a string.")
    
    sdict = {
        "linalg": np.linalg,
        "math": math,
        "numpy": np,
    }

    result = []
    reshape = None
    if isinstance(sf, dict):
        modname, func = sf["func"].split(':')
        reshape = sf.get("reshape")
    else:
        modname, func = sf.split(':')

    if reshape is not None:
        arg = np.array(values).reshape(reshape)
    else:
        arg = values
        
    call = getattr(sdict[modname], func)
    return call(arg)
    
def is_number(s):
    """Determines if the given string is a number.
    
    Args:
        s (str): the string to checke.
    """
    
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try: # pragma: no cover
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False
    
def is_nested(d):
    """Determines if a dictoinary is nested, i.e. contains another dictionary.

    Args:
        d (dict): dictionary to test.
    """
    for k, v in d.items():
        if k[-1] == '*':
            return True
        elif isinstance(v, dict) and is_nested(v):
            return True

    return False

def get_suffix(d, k, index, values):
    """Returns the suffix for the specified key in the dictionary that
    is creating a parameter grid.

    Args:
        d (dict): the dictionary being turned into a grid.
        k (str): the key in the dictionary.
        index (int): the index for the value (gets used as the default suffix).
        values (str, list, float): the value for the parameter.
    """
    from matdb.utility import special_functions
    nk = k[0:-1]
    suffix = "{0}_suffix".format(nk)
    ssuff = suffix + '*'
    
    if suffix in d and (isinstance(d[suffix], dict) or ':' in d[suffix]):
        keyval = special_functions(d[suffix], values)
    elif suffix in d and isinstance(suffix, six.string_types):
        keyval = d[suffix].format(values)
    elif ssuff in d:
        keyval = d[ssuff][index]
    else:
        keyval = index+1
    
    if isinstance(keyval, float):
        return "{0}-{1:.2f}".format(nk[:3], keyval)
    else:
        return "{0}-{1}".format(nk[:3], keyval)

def get_grid(d, suffices=None):
    """Recursively generates a grid of parameters from the dictionary of parameters
    that has duplicates or wildcars in it. 
    
    Args:
       d (dict): the dictionary to be turned into a grid.
       suffices (list): an optional list of suffices.
    
    Returns:
       A dictionary of (key: value) where the key is the suffix string for
       the parameters and the value are the exact parameters for each
       entry in the grid.
    """
    dcopy = d.copy()
    stack = [(dcopy, None)]
    result = {}
    
    if suffices is None:
        suffices = {k: v for k, v in d.items() if "suffix" in k[-8:]}
        for k in suffices:
            del dcopy[k]
    else:
        for k,v in d.items():
            if "suffix" in k[-8:]:
                suffices[k] = v
        for k in suffices:
            if k in dcopy:
                del dcopy[k]            
        
    while len(stack) > 0:
        oned, nsuffix = stack.pop()
        for k, v in sorted(oned.items()):
            if k[-1] == '*':
                nk = k[0:-1]                    
                for ival, value in enumerate(v):
                    suffix = get_suffix(suffices, k, ival, value)
                    dc = oned.copy()
                    del dc[k]
                    dc[nk] = value
                    compsuffix = suffix if nsuffix is None else '-'.join(map(str, (nsuffix, suffix)))
                    stack.append((dc, compsuffix))
                break
            elif isinstance(v, dict) and is_nested(v):
                blowup = get_grid(v, suffices)
                for rsuffix, entry in sorted(blowup.items()):
                    dc = oned.copy()
                    dc[k] = entry
                    compsuffix = rsuffix if nsuffix is None else '-'.join(map(str, (nsuffix, rsuffix)))
                    stack.append((dc, compsuffix))
                break
        else:
            result[nsuffix] = oned
            
    return result

class ParameterGrid(collections.MutableSet):
    """An ordered list of the paramater combinations for the database.
    Values are the suffixes of the combinations of parameters as
    tuples: e.g. (8, "dog", 1.2) for "dim", "animal", "temperature"
    ({"animal*": ["dog", "cat", "cow"], "dim*": [[],[],[]],
    "temperature": 1.2})
    Args:
        params (dict): the paramaters needed to build the database.
    
    Attributes:
        values (dict): keys are the suffix tuple and the values are the 
          actual values needed by the database.
        keys (list): the `str` names of the different parameters in
          the database.
    """

    def __init__(self, params):
        for k in ["root","parent","atoms"]:
            if k in params:
                params.pop(k)
        grid = get_grid(params)
        #add these items to the set.
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        self.values = {}
        self.params = params
        for i, v in grid.items():
            if i is not None:
                self.add(i,v)
            
    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def __getitem__(self, key):
        return self.values[key]
                        
    def add(self, key, value):
        """Adds key to the set if it is not already in the set.

        Args:
            key (tuple): Anything that could be added to the set.
            value (tuple): The actual values that the suffix's 
              correspond to.
        """
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]
            self.values[key] = value
        else:
            msg.warn("The key {} already exists in the set, ignoring addition.".format(key))

    def discard(self, key):
        """Removes the key from the set.

        Args:
            key (tuple): An element of the set.
        """        
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev
            self.values.pop(key,None)

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            if curr[0] is not None:
                yield curr[0]
            curr = curr[2]

    def pop(self, key):
        """Removes an element from the set.

        Args:
            key (tuple): An element of the set.
        """
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other): 
        if isinstance(other, ParameterGrid):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other) 

def is_uuid4(uuid_string):
    """Determines of the string passed in is a valid uuid4 string.
    """
    from uuid import UUID
    
    try:
        val = UUID(uuid_string, version=4)
    except:
        return False

    # If the uuid_string is a valid hex code, 
    # but an invalid uuid4,
    # the UUID.__init__ will convert it to a 
    # valid uuid4. This is bad for validation purposes.

    return val.hex == uuid_string.replace('-','')    

def dbcat(files, output, sources=None, docat=True, **params):
    """Constructs a new database file using a set of existing database files.

    .. note:: This function is important because it enforces reproducibility. It
      assigns a version number to the new database file and stores a config file
      with the specific details of how it was created.

    Args:
        files (list): of `str` paths to files to combine to create the new
          file.
        sources (list): of `str` sources as a reference for a created file.
        output (str): path to the file to write the combined files to.
        docat (bool): when True, perform the concatenation; otherwise, just
          create the config file.
        params (dict): key-value pairs that characterize *how* the database was
          created using the source files.
    """
    from uuid import uuid4
    from datetime import datetime
    from matdb import __version__
    from matdb.database.utility import dbconfig
    import json
    
    confpath = output + ".json"
    config = {
        "version": str(uuid4()),
        "sources": [],
        "timestamp": datetime.utcnow(),
        "matdb": __version__
    }
    if sources is None:
        for dbpath in files:
            _dbconf = dbconfig(dbpath)
            config["sources"].append((dbpath, _dbconf.get("version")))
    else:
        config["sources"] = sources

    config.update(params)

    try:
        with open(confpath, 'w') as f:
            json.dump(config, f, default=datetime_handler)

        if docat:
            if len(files) > 1 or (len(files) == 1 and files[0] != output):
                cat(files, output)
    except:
        from os import remove
        remove(confpath)                

def convert_dict_to_str(dct):
    """Recursively converts a dictionary to a string.

    Args:
        dct (dict): the dictionary to be converted.

    Returns:
        A string of the dictionaries values.
    """

    dict_str = ''
    for key, val in sorted(dct.items()):
        if isinstance(val,dict):
            dict_str += "'%s':'%s';"%(key,convert_dict_to_str(val))
        else:
            dict_str += "'%s':'%s';"%(key,val)

    return dict_str

def check_deps():
    """Checks that the needed dependencies have been installed and that
    they were all installed by pip.

    Returns:
        A dictionary of the dependencies and their version numbers.
    """

    req_pckgs = required_packages()    
    versions = {}
    instld_pckgs = [l.strip() for l in execute(["pip freeze"], ".", venv=True)["output"]]

    for pkg in instld_pckgs:
        if "==" in pkg:
            name, version = pkg.strip().split("==")
        else:
            name, version = pkg, None
        if name in req_pckgs:
            count = version.count(".")
            if version is not None and (count in [1,2,3] and ("/" not in version
                                                              and ":" not in version
                                                              and  "-" not in version
                                                              and ".com" not in version)):
                versions[name] = version
                req_pckgs.remove(name)
            else: #pragma: no cover There won't be locally installed
                  #packages on the test machines.
                msg.err("Cannot run `matdb` with locally installed package {} "
                        "as it would not be reproducable.".format(pkg))
    if len(req_pckgs) >= 1:
        for pkg in req_pckgs:
            msg.err("Could not find required package {} in environment.".format(pkg))

    return versions

def required_packages():
    """Returns the list of required packages for matdb. These have to be
    hard coded before each commit.
    """

    return ["ase", "beautifulsoup4", "certifi",
            "chardet", "cycler", "h5py", "html5lib", "idna", "matplotlib", "mpld3",
            "numpy", "phenum", "phonopy", "pyparsing", "python-dateutil", "pytz",
            "PyYAML", "requests", "six", "subprocess32", "termcolor",
            "tqdm", "urllib3", "webencodings", "lazy-import", "seekpath"]

