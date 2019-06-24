"""Once a database of configurations has been built, the next
   step is to train some interatomic potentials. This script
   wraps the `teach_sparse` functionality in QUIP and connects
   the database of configurations to it according to best
   practices. Note that 2+3+xxx training is an iterative process
   and so this script needs to be run multiple times.

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
                (("Continue training the next step in the sequence after "
                  "execution of the previous example is completed."), 
                 "matdb.py system.yaml -c",
                 "This generates a new jobfile for execution and updates the "
                 "parameters based on fitting errors. This should be "
                 "called once for each additional step in the overall "
                 "sequence of potential training steps.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

_script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-t": {"action": "store_true",
           "help": ("Create the database splits and the initial job script"
                    "for training the first steps in each sequence.")},
    "-x": {"action": "store_true",
           "help": ("Submit the job array files for the next fits.")},
    "-v": {"action": "store_true",
           "help": ("Validate the potential's energy, force and virial"
                    " predictions.")},
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
                          "over from scratch.")},
    "--data": {"help": "Specify which data file to validate."},
    "--tfilter": {"nargs": "+",
                 "help": ("Specify a list of patterns to match against _fit_ "
                          "names that should be *included*.")},
    "--sfilter": {"nargs": "+",
                  "help": ("Specify a list of patterns to match against _step_ "
                           "names that should be *included*.")},
    "--energy": {"action": "store_true",
                 "help": "With --validate, calculate the energies."},
    "--force": {"action": "store_true",
                 "help": "With --validate, calculate the forces."},
    "--virial": {"action": "store_true",
                 "help": "With --validate, calculate the virial."},
    "--dryrun": {"action": "store_true",
                 "help": ("With --execute, only do a dry-run to see what"
                          " would happen.")},
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
    for arg, options in _script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args

def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    print("matdb  Copyright (C) 2019  HALL LABS")
    print("This program comes with ABSOLUTELY NO WARRANTY."
    print("This is free software, and you are welcome to redistribute it under "
          "certain conditions.")
    if args is None:
        return

    import numpy as np
    from matdb import msg

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    from matdb.database import Controller
    cdb = Controller(args["dbspec"])
    if args["xyz"]:
        cdb.split(cdb.trainers.split, recalc=args["recalc"])
    if args["t"]:
        cdb.trainers.jobfiles()
    if args["x"]:
        cdb.trainers.execute(args["tfilter"], args["sfilter"], args["dryrun"])
    if args["v"]:
        vdict = cdb.trainers.validate(args["data"], args["tfilter"],
                                      args["sfilter"], args["energy"],
                                      args["force"], args["virial"])
        if "e_ref" in vdict:
            e_err = np.std(vdict["e_ref"]-vdict["e_pot"])
            msg.info("Energy RMS: {0:.4f}".format(e_err))
        if "f_ref" in vdict:
            f_err = np.std(vdict["f_ref"].flatten()-vdict["f_pot"].flatten())
            msg.info("Force RMS: {0:.4f}".format(f_err))
        if "v_ref" in vdict:
            v_err = np.std(vdict["v_ref"].flatten()-vdict["v_pot"].flatten())
            msg.info("Virial RMS: {0:.4f}".format(v_err))
        
    if args["status"]:
        cdb.trainers.status()
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
