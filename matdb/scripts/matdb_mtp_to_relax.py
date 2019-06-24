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
"""Constructs the to-relax.cfg file for the mtp potential
   to learn on for the desired system.
"""

#!/usr/bin/python
def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
    script = "MATDB mtp to-relax.cfg file construction"
    explain = ("Constructs the to-relax.cfg file for the mtp potential "
               "to learn on for the desired system.")
    contents = [(("Construct the to-relax.cfg file."), 
                 "matdb_mtp_to_relax.py",
                 "This constructs the to-relax.cfg file.")]
    required = ("'to_relax.json' file with construction settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

_script_options = {}
"""
dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse
    import sys
    from matdb import base
    pdescr = "MATDB MTP to-relax.cfg constructor"
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
    
    import numpy as np
    import json
    from matdb import msg
    from os import path
    from matdb.fitting.mtp import create_to_relax

    print("matdb  Copyright (C) 2019  HALL LABS")
    print("This program comes with ABSOLUTELY NO WARRANTY.")
    print("This is free software, and you are welcome to redistribute it under "
          "certain conditions.")
    if path.isfile('to_relax.json'):
        with open('to_relax.json', "r") as f:
            settings = json.load(f)
            create_to_relax(settings)

    else:
        msg.err("Could not find 'to_relax.json' file needed for computation.")
        return
        
    if "status" in args:
        cdb.trainers.status()
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
