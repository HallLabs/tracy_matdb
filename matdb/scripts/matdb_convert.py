#Copyright (C) 2019  HALL LABS
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#If you have any questions contact: wmorgan@tracy.com
"""Because matdb supports many different interatomic potential
   trainers, and since each possibly has its own custom format
   for the configuration databases, conversion is inevitable.
   Matdb uses HDF5 for compactness, but can convert between 
   formats using this script.
"""

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

_script_options = {
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

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    pdescr = "MATDB Database Converter"

    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in _script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args
        
def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    print("matdb  Copyright (C) 2019  HALL LABS")
    print("This program comes with ABSOLUTELY NO WARRANTY.")
    print("This is free software, and you are welcome to redistribute it under "
          "certain conditions.")
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
