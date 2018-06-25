 #!/usr/bin/python
import argparse
from os import path, remove, listdir
from shutil import move
import sys

import yaml

from matdb import msg
from matdb import msg, base
from matdb.utility import execute
from matdb.database import Controller

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    # from matdb import msg
    script = "MATDB Mover and Duplicator"
    explain = ("During regular research prototyping, it is often useful to "
               "retry a fit with different parameters (for example while "
               "changing the code). This script allows an existing trainer "
               "to be renamed or duplicated.")
    contents = [(("Rename a soap fit."),
                 "matdb_move.py system -p soap.mb --to soap-2.mb",
                 "Use `matdb_find` to test your patterns function.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-t": {"help": ("When specified, search the trainer context."),
           "action": "store_true"},
    "--to": {"help": "New name for the fit.", "required": True},
    "-p": {"help": ("Specify the search pattern(s)"), "nargs": '+',
           "required": True},
    "--dupe": {"help": "Duplicate the fit to the new location.",
               "action": "store_true"}
    }
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    # import argparse
    # import sys
    # from matdb import base
    pdescr = "MATDB Mover and Duplicator"
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)

    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args

def _get_targets(root, dbspec):
    """Returns a list of trainer config definitions from the yml files.

    Args:
        root (str): path to the root folder of the `matdb`.
        dbspec (dict): *raw* specification loaded from the YML (i.e., without
          recursive replacement of file directives).

    Returns:

    tuple: `(fits, targets, fitroot)`, where `fits` is a list of the *raw* config
    directives for fitters (which may include file directives) and `targets` is
    a `dict` with keys being trainer names and values being config
    dictionaries. `fitroot` is the path to the context directory where file
    directive config specifications are saved.
    """
    fits = dbspec["fitting"].get("fits", [])
    fitroot = None
    targets = {}

    if "context" in dbspec:
        context = dbspec["context"].get("fitting")
        fitroot = path.join(root, context)
        for k in fits:
            if isinstance(k, dict):
                name = k.get("name")
                targets[name] = k
            else:
                with open(path.join(fitroot, "{}.yml".format(k[1:]))) as f:
                    config = yaml.load(f)
                name = config.get("name")
                targets[name] = config
    else:
        for k in fits:
            name = k.get("name")
            targets[name] = k

    return fits, targets, fitroot

def _move_trainer(trainer, dbspec=None, dupe=False, to=None, **kwargs):
    """Moves the specified trainer.

    Args:
        trainer (matdb.fitting.basic.Trainer): trainer whose directory structure
          should be moved.
        dbspec (str): path to the spec yml file for the entire `matdb`.
        dupe (bool): when True, duplicate the trainer so that the existing
          trainer remains at its current location and a copy is placed under the
          new name.
        to (str): name of the new trainer; all references, including those in
          configuration files, will be changed and folders will be moved.
    """
    #We need to:
    # 1) Rename the trainer in the yml settings that defined the trainer.
    # 2) If duplicating, copy the yml settings file, add its reference to the
    #    system yml file.
    # 3) Move and/or duplicate the folder structure for the trainer.
    origspec = dbspec
    dbspec = path.abspath("{}.yml".format(dbspec))
    with open(dbspec) as f:
        spec = yaml.load(f)
    root = path.dirname(dbspec)

    fits, targets, fitroot = _get_targets(root, spec)
    if to in fits:
        raise ValueError("The specified fit name is already taken! "
                         "Choose a different new name for the fit.")
    mvname = trainer.parent.name
    try:
        settings = targets[mvname]
    except KeyError:
        msg.err("Cannot find trainer specification for {}.".format(mvname))

    #Handle duplication if that was specified. Otherwise, just move the config
    #specification and rename the trainer within it.
    if dupe:
        dsettings = settings.copy()
        dsettings["name"] = to
        if fitroot is not None:
            with open(path.join(fitroot, "{}.yml".format(to)), 'w') as f:
                yaml.dump(dsettings, f)
            spec["fitting"]["fits"].append(":{}".format(to))
        else:
            spec["fitting"]["fits"].append(dsettings)
    else:
        settings["name"] = to
        if fitroot is not None:
            with open(path.join(fitroot, "{}.yml".format(to)), 'w') as f:
                yaml.dump(settings, f)
            remove(path.join(fitroot, "{}.yml".format(mvname)))
            spec["fitting"]["fits"].append(":{}".format(to))
        else:
            spec["fitting"]["fits"].append(settings)

    #Next, move the actual directory. We have to be careful because the first
    #level of subdirectories also have the name in them, but they also may have
    #seed directives for the special `*` and `^` notation. If duplication is
    #switched on, we *still* move the folder because a new one with the same
    #settings will be created to replace the original one.
    src = path.join(root, mvname)
    dst = path.join(root, to)
    move(src, dst)

    for objname in listdir(dst):
        _src = path.join(dst, objname)
        if path.isdir(_src):
            if mvname == objname[0:len(mvname)]:
                newname = to + objname[len(mvname):]
                _dst = path.join(dst, newname)
                move(_src, _dst)

    #Rewrite the `matdb.yml` file to include the new fit
    #specifications.
    if not dupe:
        spec["fitting"]["fits"].remove(":{}".format(mvname))
    with open(dbspec, 'w') as f:
        yaml.dump(spec, f)

def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    if args is None:
        return

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    # from matdb.database import Controller
    cdb = Controller(args["dbspec"])

    if args["t"]:
        matches = []
        for pattern in args["p"]:
            for entry in cdb.trainers.find(pattern):
                matches.append(entry)
        if len(matches) != 1:
            msg.err("Can only move a single trainer at a time."
                    "Found {} trainers.".format(len(matches)))
            exit(-1)

        _move_trainer(matches[0], **args)

if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
