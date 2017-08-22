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
    
def _get_reporoot():
    """Returns the absolute path to the repo root directory on the current
    system.
    """
    from os import path
    import matdb
    medpath = path.abspath(matdb.__file__)
    return path.dirname(path.dirname(medpath))

reporoot = _get_reporoot()
"""The absolute path to the repo root on the local machine.
"""
