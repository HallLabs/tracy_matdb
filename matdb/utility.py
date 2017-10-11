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
            printerr=True, **kwargs):
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

    msg.std("Executing `{}` in {}.".format(' '.join(args), folder), 2)
    pexec = Popen(' '.join(args), shell=True, executable="/bin/bash", **kwargs)
    
    if wait:
        from os import waitpid
        waitpid(pexec.pid, 0)

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
        
reporoot = _get_reporoot()
"""The absolute path to the repo root on the local machine.
"""
