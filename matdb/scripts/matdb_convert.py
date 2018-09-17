#!/usr/bin/python
from os import path
from matdb import msg
import argparse
import sys

from matdb import base
from matdb.database import Controller
from matdb.atoms import AtomsList

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    script = "MATDB Database Converter"
    explain = ("Because matdb supports many different interatomic potential "
               "trainers, and since each possibly has its own custom format "
               "for the configuration databases, conversion is inevitable. "
               "Matdb uses HDF5 for compactness, but can convert between "
               "formats using this script.")
    contents = [(("Convert the neb-1 database to XYZ."), 
                 "matdb_convert.py system -p neb-3/* --format xyz -o neb-3.xyz",
                 "Use `matdb_find` to test your patterns function.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "--format": {"help": "New format to generate.", "required": True,
                 "choices": ["xyz", "cfg"]},
    "-p": {"help": ("Specify the search pattern(s)"), "nargs": '+',
           "required": True},
    "-o": {"help": "Specify the name of the output file to convert to.",
           "required": True},
    "--overwrite": {"help": ("When specified, overwrite the output file "
                             "if it already exists."), "action": "store_true"}
    }
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    pdescr = "MATDB Database Converter"
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args
        
def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    if args is None:
        return

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.

    cdb = Controller(args["dbspec"])

    matches = []
    configs = AtomsList()
    for pattern in args["p"]:
        for entry in cdb.find(pattern):
            for iatoms in entry.iconfigs:
                configs.append(iatoms)

    if args["format"] == "xyz":
        from matdb.conversion import to_xyz
        target = path.abspath(path.expanduser(args["o"]))
        to_xyz(configs, target, args["overwrite"])
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
