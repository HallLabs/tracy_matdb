"""Depending on the `matdb` context, it is often useful to watch
   the output that is being piped to the batch output file. This
   script is a simple handler around `tail -f` that keeps track
   of file paths for output files.

Copyright (C) 2019  HALL LABS

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

If you have any questions contact: wmorgan@tracy.com
"""

#!/usr/bin/python
from os import path
from matdb import msg
import argparse
import sys

from matdb import base

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    script = "MATDB Watcher for Groups/Trainers and General Contexts"
    explain = ("Depending on the `matdb` context, it is often useful to watch "
               "the output that is being piped to the batch output file. This "
               "script is a simple handler around `tail -f` that keeps track "
               "of file paths for output files.")
    contents = [(("Watch a GAP fit's output stream."), 
                 "matdb_watch.py system.yaml -t -p soap.mb",
                 "Use `matdb_find` to test your patterns function uses the "
                 "format dbname/group/seed")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

_script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-d": {"help": ("When specified, search the database context."),
           "action": "store_true"},
    "-t": {"help": ("When specified, search the fitter/trainer context."),
           "action": "store_true"},
    "-p": {"help": ("Specify the search pattern(s)"), "nargs": '+',
           "required": True}
    }
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    pdescr = "MATDB Context Execution Watcher"
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in _script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args

def _generic_find(controller, heading, patterns):
    """Performs a generic find operation on the specified controller and formats
    the output in color.

    Args:
        controller: an instance of :class:`matdb.database.Controller` or
          :class:`matdb.fitting.Controller`. The specified controller's `find`
          method is used for the lookup.
        heading (str): title to print before the table of discovered values.
        patterns (list): list of `str` patterns to search for.
    """
    msg.info(heading)
    msg.info("--------------------------")
    msg.blank()
    for pattern in patterns:
        for entry in controller.find(pattern):
            if hasattr(entry, "uuid"):
                eid = entry.uuid
            elif hasattr(entry, "fqn"):
                eid = entry.fqn
            else:
                eid = entry.name
            text = "{} | {} ".format(eid, entry.root)
            msg.arb(text, [msg.cenum["cwarn"],
                           msg.cenum["cstds"]], '|')

def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    print("matdb  Copyright (C) 2019  HALL LABS")
    print("This program comes with ABSOLUTELY NO WARRANTY."
    print("This is free software, and you are welcome to redistribute it under "
          "certain conditions.")
    if args is None:
        return

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    from matdb.database import Controller
    cdb = Controller(args["dbspec"])

    if args["d"]:
        _generic_find(cdb, "Database Context Instances", args["p"])
    if args["t"]:
        _generic_find(cdb.trainers, "Fitter Context Instances", args["p"])
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
