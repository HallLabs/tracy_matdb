 #!/usr/bin/python
def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
    script = "MATDB Automated Materials Database Plotter"
    explain = ("While constructing the database, there are several times when "
               "plots are needed: phonon bands, prediction errors, etc. This "
               "script provides a unified interface for all such plotting.")
    contents = [(("Plot the phonon bands for the structure labeled as 'PdAg25' "
                  "in the `system.yaml` database specification."), 
                 "matdb_plot.py system.yaml -d PdAg25.phonon.dynmatrix --bands",
                 "If the `bands.yaml` file hasn't been created yet using "
                 "`phonopy`, then this will create that file. Note that "
                 "this does *not* plot any predictions; use -p to specify "
                 "the names of MTP/GAP potentials to calculate phonons with."),
                ("Plot the supercell convergence of phonon bands for the "
                 "structure labeled as 'PdAg25' in the `system.yaml` database"
                 " specification.", 
                 "matdb_plot.py system.yaml -d PdAg25.phonon*.dynmatrix --bands",
                 "Notice that * can be used as a pattern. For example, to plot "
                 "all configurations you could use *.phonon*.dynmatrix, "
                 "assuming that the database was called `phonon`.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-d": {"help": ("Specify the pattern of the databases to work with. "
                    "For example `Pd.phonon-2.dynmatrix` specifies the "
                    "dynamical matrix *step* for suffix `2` of the phonon"
                    " *database* of configuration `Pd`.")},
    "--bands": {"action": "store_true",
                "help": ("Plot the phonon bands for the structures given "
                         "by '-d'")},
    "--dim": {"nargs": "+", "default": 2,
              "help": ("Specify the supercell dimensions for the phonon "
                       "calculations.")},
    "--pots": {"nargs": "+",
               "help": ("Specify the names of GAP/MTP potential parameter "
                        "files to compute properties for.")},
    "--title": {"help": ("Override the default title for plotting; "
                         "use {} for formatting chemical formula."),
                "default": "{} Phonon Spectrum"},
    "--npts": {"type": int, "default": 100,
               "help": ("Specify the number of points to sample along the "
                        "special path in k-space.")},
    "--figsize": {"nargs": 2, "type": float, "default": (10, 8),
                  "help": "Specify the size of the figure in inches."},
    "--save": {"help": "Specify the name of a file to save the plot to." },
    "--nbands": {"help": "Number of bands to plot.",
                 "default": 4}
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
    pdescr = "MATDB Database Plotting"
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

    #We allow any number of databases to be plotted together at the same
    #tame. The `s` argument gives a `fnmatch` pattern for database paths.
    #For example `Pd.phonon-2.dynmatrix` specifies the dynamical matrix *step*
    #for suffix `2` of the phonon *database* of configuration `Pd`.
    dbs = cdb.find(args["d"])
    
    if args["bands"]:
        from matdb.plotting.comparative import band_plot
        band_plot(phondbs, **args)
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
