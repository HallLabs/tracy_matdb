"""This is another stub for unit testing that mimics the behavior of
the `mlp` executable. 
"""

from os import path

from matdb.utility import touch, _get_reporoot, copyonce

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
    script = "MATDB Module mlp Testing Stub"
    explain = ("Since unit tests run in an unconfigured environment, we need to"
               "generate stubs to handle any commands that aren't normally "
               "there. This includes `vasp` and `module`. This script mimcs the"
               "behavior of `mlp` (but just creates file place holders). The "
               "`vasp.py` stub mimics `vasp` by copying output." )
    contents = []
    required = ("")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
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

    cmd = args["cmd"]
    print(cmd)
    keyword = cmd[0].strip()
    print(keyword)
    if keyword == "calc-grade":
        touch("state.mvs")
    elif keyword == "select-add":
        template_root = path.join(_get_reporoot(), "tests", "fitting", "files")
        touch("new_training.cfg")
        # src = path.join(template_root,"mtp_new_training.cfg")
        # dest = "new_training.cfg"
        # copyonce(src, dest)
    elif keyword == "convert-cfg":
        template_root = path.join(_get_reporoot(), "tests", "fitting", "files")
        for i in range(1,11):
            touch("POSCAR{0}".format(i))
            # src = path.join(template_root,"mtp_convert_POSCAR{0}".format(i))
            # dest = "POSCAR{0}".format(i)
            # copyonce(src, dest)
    
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
