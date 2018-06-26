#!/usr/bin/python
import argparse
from os import path, remove
from shutil import move
import sys

import h5py
import numpy as np

from matdb import msg, base
from matdb.atoms import Atoms
from matdb.io import load_dict_from_h5, save_dict_to_h5
from matdb.database import Controller

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    # from matdb import msg
    script = "MATDB Precomp Atoms Fixer"
    explain = ("In a previous version of `matdb`, the pre-comp-atoms.h5 "
               "file would be removed after cleanup. But, if the cleanup "
               "produced any unnoticed errors, we may lose the pre-comp "
               "file and it is hard to reproduce. This script remedies that.")
    contents = [(("Reconstitute the missing pre-comp file."),
                 "matdb_move.py system -p distort-1/Distortion/*",
                 "Use `matdb_find` to test your patterns function.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-p": {"help": ("Specify the search pattern(s)"), "nargs": '+',
           "required": True},
    "--purge": {"help": ("If true, the existing `atoms.h5` file is removed "
                         "after the pre_comp file is created."),
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
    pdescr = "MATDB Pre-Comp Fixer"
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)

    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args

def _fix_precomp(db, purge=False):
    """Reconstitutes the `pre_comp_atoms.h5` file by copying the existing
    `atoms.h5`, and applying some modifications to the atoms object.

    Args:
        db (matdb.database.Group): a group to create pre_comp files for.
        purge (bool): when True, remove the `atoms.h5` files after the pre-comp
          file is created.
    """
    for aid, apath in db.configs.items():
        target = path.join(apath, "atoms.h5")
        with h5py.File(target,"r") as hf:
            data = load_dict_from_h5(hf)
        data["pbc"] = np.array([1,1,1], dtype=bool)

        newtarg = path.join(apath, "pre_comp_atoms.h5")
        with h5py.File(newtarg, "w") as hf:
            save_dict_to_h5(hf, data, '/')

        if path.isfile(newtarg) and purge:
            remove(target)

def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    if args is None:
        return

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    # from matdb.database import Controller
    cdb = Controller(args["dbspec"])

    matches = []
    for pattern in args["p"]:
        for entry in cdb.find(pattern):
            matches.append(entry)

    for db in matches:
        _fix_precomp(db)

if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
