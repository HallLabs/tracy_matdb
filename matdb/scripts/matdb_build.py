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
           "help": "Run the setup method for each database."}
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
    from matdb.database.controller import Controller
    cdb = Controller(args["dbspec"])
    if args["s"]:
        cdb.setup()
    
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
