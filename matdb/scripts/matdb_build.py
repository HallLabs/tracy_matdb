 #!/usr/bin/python
def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
    script = "MATDB Automated Materials Database Constructor"
    explain = ("In order to apply machine learning to produce potentials "
               "we first need a collection of atomic configurations from "
               "which the energies, forces and virials can be learned. "
               "Ideally, the creation process for such configurations should "
               "be scientific and reproducible. This script uses the "
               "`matdb.database package to produce just such a set of "
               "configurations.")
    contents = [(("Setup as many DFT folders as possible, taking dependencies "
                  "between databases into account."), 
                 "matdb.py system.yaml -s",
                 "This creates a directory for each configuration in the YAML "
                 "configuration file 'system.yaml'. Within in configuration's "
                 "directory, a folder for each database type will be ready, "
                 "depending on how many calculations have already been "
                 "performed.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-s": {"action": "store_true",
           "help": "Run the setup method for each database."},
    "-x": {"action": "store_true",
           "help": ("Submit the job array file for each database that has "
                    "folders ready to run. If --recover is specified, then "
                    "the recovery jobfile is submitted instead.")},
    "-e": {"action": "store_true",
           "help": ("Process the databases that have completed execution "
                    "so that results can be extracted.")},
    "--status": {"action": "store_true",
                "help": ("Determines status of the databases "
                         "based on presence of completed VASP "
                         "directories. Sanity check before `-x`.")},
    "--rerun": {"type": int, "default": 0,
                "help": ("Re-run the specified option, even if it has already "
                         "been done before. Higher values re-run at a deeper level.")},
    "--dfilter": {"nargs": "+",
                  "help": ("Specify a list of patterns to match against _database_ "
                           "names that should be *included*.")},
    "--clean": {"default": "default",
                "type": str,
                "help": ("Specify the cleanup level for the database, 'aggressive', "
                         "'default', or 'light'.")},
    "--busy": {"action": "store_true",
               "help": ("Display a list of configurations that haven't "
                        "finished running in DFT yet.")},
    "--recover": {"action": "store_true",
                  "help": ("Creates a jobfile for those configs that didn't "
                           "finish computing so that they can be re-run.")},
    "--dryrun": {"action": "store_true",
                  "help": ("For execution, doesn't actually submit any jobs, "
                           "just prints what would be done.")},
    "--asis": {"action": "store_true",
               "help": ("When specified, extract results even if the"
                        " calculations didn't converge completely.")}
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
    pdescr = "MATDB Database Constructor"
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
    from matdb.database import Controller
    cdb = Controller(args["dbspec"])
    if args["s"]:
        cdb.setup(args["rerun"], args["dfilter"])
    if args["x"]:
        cdb.execute(args["recover"], args["dfilter"], dryrun=args["dryrun"])
    if args["e"]:
        cdb.extract(args["dfilter"], cleanup=args["clean"], asis=args["asis"])

    if args["recover"] and not args["x"]:
        cdb.recover(args["rerun"], args["dfilter"])
        
    if args["status"]:
        cdb.status(args["busy"], args["dfilter"])
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
