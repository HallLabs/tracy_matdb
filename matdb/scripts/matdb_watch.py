 #!/usr/bin/python
from os import path
from matdb import msg

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
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

script_options = {
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
    import argparse
    import sys
    from matdb import base
    pdescr = "MATDB Context Execution Watcher"
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
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
        patterns (list): of `str` patterns to search for.
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
