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
