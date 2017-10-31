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
                    " *database* of configuration `Pd`."),
           "nargs": "+"},
    "--bands": {"action": "store_true",
                "help": ("Plot the phonon bands for the structures given "
                         "by '-d'")},
    "--dim": {"nargs": "+", "default": 2,
              "help": ("Specify the supercell dimensions for the phonon "
                       "calculations.")},
    "--pots": {"nargs": "+",
               "help": ("Specify a list of patterns for trainer sequences "
                        "to compute/plot properties for.")},
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
                 "default": 4},
    "--generate": {"action": "store_true",
                   "help": ("Generate the sub-plots for HTML interactive"
                            " package.")},
    "--plots": {"default": "efvodt",
                "help": ("Specify the kinds of plots to generate as a string: "
                         "1. Energy Correlation Plot (`e`); "
                         "2. Force Correlation Plot (`f`); "
                         "3. Virial Correlation Plot (`v`); "
                         "4. Energy vs. Volume Plot (`o`); "
                         "5. Dimer Plots for 2-body Potentials (`d`); "
                         "6. Trimer Plots for 3-body Potentials (`t`); ")},
    "--base64": {"action": "store_true",
                 "help": ("When true, don't save plots, rather base64 "
                          "encode them for HTML.")},
    "--folder": {"help": ("Specify the folder to store the plots in.")},
    "--html": {"action": "store_true",
               "help": ("Generate the interactive HTML plotter for all "
                        "databases an potentials specified.")},
    "--gzip": {"action": "store_true",
               "help": ("For --html, gzip the resulting directory of plots "
                        "and the html page.")},
    "--splits": {"nargs": "+",
                "help": "Specify a list of global splits to use for plotting."},
    "--subset": {"choices": ["train", "holdout", "super"], "default": "holdout",
                 "help": "Specify which database subset to use."}
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
    from os import getcwd, path
    cdb = Controller(args["dbspec"])

    #We allow any number of databases to be plotted together at the same
    #tame. The `s` argument gives a `fnmatch` pattern for database paths.
    #For example `Pd.phonon-2.dynmatrix` specifies the dynamical matrix *step*
    #for suffix `2` of the phonon *database* of configuration `Pd`.
    dbs = []
    if args["d"]:
        for dbp in args["d"]:
            dbs.extend(cdb.find(dbp))

    pots = []
    if args["pots"]:
        for potp in args["pots"]:
            pots.extend(cdb.trainers.find(potp)) 
    if len(pots) == 0:
        #Try and get a potential from the current directory.
        dot = getcwd()
        step = path.dirname(dot)
        seq = path.dirname(step)
        probable = "{}.{}".format(seq, step)
        _pot = cdb.trainers.find(probable)
        if _pot:
            pots.append(_pot)

    if args["bands"]:
        from matdb.plotting.comparative import band_plot
        raise NotImplementedError("Still need to refactor since db/pot major refactor.")
        band_plot(dbs, **args)

    if args["generate"]:
        from matdb.plotting.potentials import generate
        from cPickle import dump
        
        if len(pots) > 1:
            raise ValueError("Generate only operates for a single trainer "
                             "at a time; don't specify so many patterns.")
        pot = pots[0]
        pdis = generate(args["plots"], pot.calculator,
                        pot.configs(args["subset"]),
                        args["folder"], args["base64"])

        pklname = "{}-{}-plotgen.pkl".format(args["subset"], args["plots"])
        target = path.join(pot.root, pklname)
        with open(target, 'w') as f:
            dump(pdis, f)
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
