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
    contents = [(("Setup the `.matdb.json` file for the `vasprun.xml` in the "
                  "current directory; place it in `/tests/Pd-2`."), 
                 "module.py load arbitrary/module --copy vasprun.xml --xdir /tests/Pd-2",
                 "This simply copies the full path to `vasprun.xml` in the "
                 "current folder to a `.matdb.json` file in the execution dir.")]
    required = ("")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "cmd": {"help": "Which command to mimic; may not do anything.",
            "choices": ["list", "load", "unload"],
            "nargs": "+"},
    "--copy": {"help": ("Specify the name(s) of a list of files in the current "
                        "directory that should be written to a `.matdb.json` "
                        "file in the *execution* directory."),
               "nargs": "*"},
    "--xdir": {"help": ("Path to the execution directory where the "
                        "`.matdb.json` file will be created.")}
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

    if not args["xdir"]:
        return
    
    from os import path
    matdb = {}
    current = path.abspath('.')
    for f in args["copy"]:
        matdb[target] = path.join(current, f)

    import json
    target = path.abspath(path.expanduser(args["xdir"]))
    with open(path.join(target, ".matdb.json"), 'w') as f:
        json.dump(matdb, f)
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
