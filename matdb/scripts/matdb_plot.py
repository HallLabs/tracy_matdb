 #!/usr/bin/python
import argparse
from cPickle import dump, load
from shutil import copyfile
from os import path, mkdir
import sys

import matplotlib

from matdb import msg, base
from matdb.atoms import AtomsList
from matdb.database import Controller
from matdb.plotting.potentials import generate
from matdb.plotting.plotter import PlotManager
from matdb.plotting.comparative import band_plot
from matdb.plotting.matd3 import html

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    # from matdb import msg
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
    "-x": {"help": ("Specify a parameter name to plot along `x` axis in "
                    "interactive plots.")},
    "-y": {"help": ("Specify a parameter name to plot along `y` axis in "
                    "interactive plots.")},
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
    "--nbands": {"help": "Number of bands to plot.", "type": int},
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
                 "help": "Specify which database subset to use."},
    "--generic": {"help": ("Specify the configuration file to use for a "
                           "generic plot of the data."),
                  "nargs": '+'},
    "--valkey": {"help": ("Specify the key for parameter and property names "
                          "that should be used to compare to the potential "
                          "in comparative plots.")}
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
    pdescr = "MATDB Database Plotting"
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)

    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args

def _generate_pkl(pot, dbs=None, **args):
    """Generates a pickle file for a single potential and its default databases.
    """
    # from matdb.plotting.potentials import generate
    # from matdb.atoms import AtomsList
    # from cPickle import dump
    outdir = path.join(args["folder"], pot.fqn)
    if not path.isdir(outdir):
        mkdir(outdir)

    if dbs is not None:
        configs = AtomsList()
        for db in dbs:
            configs.extend(list(db.iconfigs))

        pdis = generate(args["plots"], pot.calculator,
                        configs, outdir, args["base64"],
                        valkey=args["valkey"])
    else:
        pdis = generate(args["plots"], pot.calculator,
                        pot.configs(args["subset"]),
                        outdir, args["base64"])

    pklname = "{}-{}-plotgen.pkl".format(args["subset"], args["plots"])
    target = path.join(outdir, pklname)
    with open(target, 'w') as f:
        dump(pdis, f)

def _generate_html(potlist, **args):
    """Generates the interactive HTML page for the potentials and databases.
    """
    # from cPickle import load
    #We create a new folder for the HTML plot inside the overall system
    #directory.
    plotdir = path.join(potlist[0].controller.db.plotdir, args["save"])
    if not path.isdir(plotdir):
        mkdir(plotdir)

    data = {}
    pklname = "{}-{}-plotgen.pkl".format(args["subset"], args["plots"])
    for ipot, pot in enumerate(potlist):
        target = path.join(args["folder"], pot.fqn, pklname)
        if not path.isfile(target):
            _generate_pkl(pot, **args)
        with open(target) as f:
            pdis = load(f)

        xf, yf = "{%s}" % args["x"], "{%s}" % args["y"]
        x, y = xf.format(**pot.params), yf.format(**pot.params)
        #Actually, use the parameters dict for the fitting.
        data[pot.fqn] = {
            "location": (x, y),
            "index": ipot
        }
        data[pot.fqn].update(pdis)

        # from shutil import copyfile
        for pdi in pdis.values():
            newname = "{}__{}".format(pot.fqn, pdi.url)
            trg = path.join(plotdir, newname)
            copyfile(path.join(args["folder"], pot.fqn, pdi.url), trg)
            pdi.filename = newname

    # from matdb.plotting.matd3 import html
    subplot_kw = {
        "title": "Interactive Potential Explorer",
        "xlabel": args["x"],
        "ylabel": args["y"]
    }
    plot_kw = {
        "alpha": 0.3
    }

    html(data, plotdir, subplot_kw=subplot_kw, plot_kw=plot_kw)

def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    if args is None:
        return

    matplotlib.use('Agg')
    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    # from matdb.database import Controller
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

    if args["generic"]:

        objs = dbs + pots
        # from matdb.plotting.plotter import PlotManager
        manager = PlotManager(cdb)
        for cname in args["generic"]:
            manager.plot(objs, cname)

    if args["bands"]:
        # from matdb.plotting.comparative import band_plot
        band_plot(dbs, pots, **args)

    if args["generate"]:
        if len(pots) > 1:
            raise ValueError("Generate only operates for a single trainer "
                             "at a time; don't specify so many patterns.")
        pot = pots[0]
        _generate_pkl(pot, dbs, **args)

    if args["html"]:
        _generate_html(pots, **args)

if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
