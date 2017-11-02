"""Utility functions for interacting with file system, shell, etc.
"""
from os import path
from six import string_types
from matdb import msg

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
        
def execute(args, folder, wait=True, nlines=100, venv=None,
            printerr=True, env_vars=None, **kwargs):
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

    if venv is not None:
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
    if path.isfile(target) or path.islink(target):# pragma: no cover
        #This will never fire for normal unit testing because we used
        #new temporary directories each time.
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
        if hasattr(obj, k) and getattr(k, v) is not None:
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
    """Prepares a :class:`datetime.datetime` for JSON serialization by encoding
    it in UNIX timestamp relative to zero-offset UTC.
    Args:
        x (datetime.datetime): to be turned into a `float`.
    Returns:
        float: number of seconds since 1/1/1970, UTC zero offset.
    """
    if isinstance(x, datetime):
        delta = x.astimezone(pytz.utc) - epoch
        return "dt:{0:.7f}".format(delta.total_seconds())

def parse_date(v):
    """Parses the date from the specified single value or list of values.
    Args:
        v (str): string representation of the :class:`datetime` returned by
          :func:`datetime_handler`.
    """
    if isinstance(v, (list, tuple)):
        return [parse_date(vi) for vi in v]
    elif isinstance(v, string_types) and v[0:3] == "dt:":
        ts = float(v[3:])
        return datetime.fromtimestamp(ts, pytz.utc)
    else:# pragma: no cover
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

def dbconfig(dbfile):
    """Returns the database configuration `dict` of the specified database file.

    Args:
        dbfile (str): path to the database file to get config information for.
    """
    import json
    confpath = dbfile + ".json"
    with open(confpath) as f:
        config = json.load(f, object_pairs_hook=load_datetime)

    return config

def dbcat(files, output, **params):
    """Constructs a new database file using a set of existing database files.

    .. note:: This function is important because it enforces reproducibility. It
      assigns a version number to the new database file and stores a config file
      with the specific details of how it was created.

    Args:
        files (list): of `str` paths to files to combine to create the new
          file. 
        output (str): path to the file to write the combined files to.
        params (dict): key-value pairs that characterize *how* the database was
          created using the source files.
    """
    from uuid import uuid4
    from datetime import datetime
    from matdb import __version__
    import json
    confpath = output + ".json"
    config = {
        "version": uuid4(),
        "sources": [],
        "timestamp": datetime.now(),
        "matdb": __version__
    }
    for dbpath in files:
        _dbconf = dbconfig(dbpath)
        config["sources"].append((dbpath, _dbconf["version"]))

    config.update(params)

    try:
        with open(confpath, 'w') as f:
            json.dump(config, f, default=datetime_handler)
        cat(files)
    except:
        from os import remove
        remove(confpath)

reporoot = _get_reporoot()
"""The absolute path to the repo root on the local machine.
"""
