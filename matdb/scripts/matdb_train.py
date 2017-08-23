 #!/usr/bin/python
def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
    script = "MATDB Automated Materials Database Constructor"
    explain = ("Once a database of configurations has been built, the next "
               "step is to train some interatomic potentials. This script "
               "wraps the `teach_sparse` functionality in QUIP and connects "
               "the database of configurations to it according to best "
               "practices. Note that 2+3+xxx training is an iterative process "
               "and so this script needs to be run multiple times.")
    contents = [(("Start training a new potential for the first time."), 
                 "matdb.py system.yaml -t",
                 "This generates the XYZ database of configurations, and sets "
                 "up the jobfile for execution."),
                (("Continue training the next GAP in the potential after "
                  "execution of the previous example is completed."), 
                 "matdb.py system.yaml -c",
                 "This generates a new jobfile for execution and updates the "
                 "`delta` parameter based on fitting errors. This should be "
                 "called once for each additional GAP in the overall potential.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-t": {"action": "store_true",
           "help": ("Create the training XYZ files and the initial job script"
                    "for training the first GAP in the potential.")},
    "-x": {"action": "store_true",
           "help": ("Submit the job array file for the next GAP.")},
    "-v": {"action": "store_true",
           "help": "Validate the potential's energy and force predictions."},
    "--status": {"action": "store_true",
                "help": ("Determines status of the GAP fitting routines and "
                         "job files. Sanity check before `-x`.")},
    "--xyz": {"action": "store_true",
              "help": ("Recreate the XYZ training and validation files using "
                       "the latest settings in `matdb.yml`")},
    "--recalc": {"type": int, "default": 0,
                 "help": ("Specify the stack depth for recalculation; this "
                          "value is decreased as the stack is traversed. In "
                          "any method with `recalc > 0`, the operation is done "
                          "over from scratch.")}
}
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
    pdescr = "MATDB Potential Fitter"
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

    import numpy as np
    from matdb import msg

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    from matdb.database.controller import Controller
    cdb = Controller(args["dbspec"])
    if args["xyz"]:
        cdb.split(cdb.trainer.split, recalc=args["recalc"])
    if args["t"]:
        cdb.trainer.write_jobfile()
    if args["x"]:
        cdb.trainer.execute()
    if args["v"]:
        vdict = cdb.trainer.validate()
        e_err = np.std(vdict["e_dft"]-vdict["e_gap"])
        f_err = np.std(vdict["f_dft"].flatten()-vdict["f_gap"].flatten())
        msg.info("Energy RMS: {0:.4f}".format(e_err))
        msg.info("Force RMS: {0:.4f}".format(f_err))
        
    if args["status"]:
        cdb.trainer.status()
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
