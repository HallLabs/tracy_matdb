#!/usr/bin/python

"""This is another stub for unit testing that mimics the behavior of `module
load/unload` at the HPC center. It accepts *extra* arguments that define the
location of the model `vasprun.xml` output that should be copied to the
execution directory created by the unit tests.
"""
def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
    script = "MATDB Module Loader/Unloader Testing Stub"
    explain = ("Since unit tests run in an unconfigured environment, we need to"
               "generate stubs to handle any commands that aren't normally "
               "there. This includes `vasp` and `module`. This script mimcs the"
               "behavior of `module` (but doesn't actually do anything). The "
               "`vasp.py` stub mimics `vasp` by copying output." )
    contents = []
    required = ("")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

_script_options = {
    "cmd": {"help": "Which command to mimic; may not do anything.",
            "nargs": "+"},
    }
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse
    import sys
    from matdb import base
    pdescr = "MATDB `module` STUB"
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
    if args is None:
        return

    with open(".matdb.module", 'w') as f:
        f.write(' '.join(args["cmd"]))
    
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
