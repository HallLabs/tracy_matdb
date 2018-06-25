#!/usr/bin/python
import argparse
from os import path
import sys

from glob import glob
from tqdm import tqdm

from matdb import msg, base
from matdb.atoms import Atoms
from matdb.transforms import _get_supers
from matdb.utility import chdir

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    # from matdb import msg
    script = "MATDB Supercell Generator"
    explain = ("When producing hessians with the Hessian group, a supercell "
               "needs to be selected. This script streamlines the selection "
               "process for multiple seeds and sizes.")
    contents = [(("Select supercells for all the seeds and sizes 32 and 64."),
                 "matdb_supercell.py * --sizes 32 64",
                 "")]
    required = ("Seed files in the `seed` directory.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "seeds": {"help": ("File patterns for choosing seeds."),
              "nargs": '+'},
    "--sizes": {"help": ("Target cell sizes to find for."),
                "nargs": '+', "type": int, "required": True},
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
    pdescr = "MATDB Supercell Selector"
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

    targets = {}
    with chdir("seed"):
        for pattern in args["seeds"]:
            #Handle the default file type, which is vasp.
            if ':' in pattern:
                fmt, pat = pattern.split(':')
            else:
                fmt, pat = "vasp", pattern
            for filename in glob(pat):
                targets[filename] = Atoms(filename, format=fmt)

    result = {}
    for filename, at in tqdm(list(targets.items())):
        result[filename] = _get_supers(at, args["sizes"])

    items = [
        ("Filename", 20, "cokay"),
        ("Supercell", 40, "cstds"),
        ("Req.", 6, "cinfo"),
        ("Act.", 6, "cgens"),
        ("rmin", 8, "cerrs"),
        ("pg", 6, "cwarn")
    ]

    msg.blank(2)
    heading = '|'.join(["{{0: ^{0}}}".format(size).format(name)
                        for name, size, color in items])
    msg.arb(heading, [msg.cenum[i[2]] for i in items], '|')
    msg.std(''.join('-' for i in range(len(heading)+1)))
    for filename, hs in result.items():
        for size, hnf in hs.items():
            names = (filename, hnf.hnf.flatten().tolist(), size,
                     hnf.size, hnf.rmin, hnf.pg)
            text = '|'.join(["{{0: <{0}}}".format(item[1]).format(name)
                                for name, item in zip(names, items)])
            msg.arb(text, [msg.cenum[i[2]] for i in items], '|')
        msg.blank(2)

    return result

if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
