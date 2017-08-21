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
                 "matdb_plot.py system.yaml -s PdAg25 --bands",
                 "If the `bands.yaml` file hasn't been created yet using "
                 "`phonopy`, then this will create that file. Note that "
                 "this does *not* plot any predictions; use -p to specify "
                 "the names of MTP/GAP potentials to calculate phonons with.")]
    required = ("'matdb.yaml' file with database settings.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "dbspec": {"help": "File containing the database specifications."},
    "-s": {"help": "Specify the name of the structure to work with."},
    "--bands": {"action": "store_true",
                "help": ("Plot the phonon bands for the structure given "
                         "by '-s'")},
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

def _band_plot(phondb, args):
    """Plots the phonon bands for the specified CLI args.

    Args:
        phondb (matdb.database.phonon.PhononDFT): `phonopy` calculation
          database instance that has DFT-accurate band information.
    """
    from matdb.phonons import bandplot
    from os import path
    DFT = phondb.bands
    names, kpath = phondb.kpath
    #matplotlib needs the $ signs for latex if we are using special characters.
    names = ["${}$".format(n) if '\\' in n else n for n in names]

    from matdb.phonons import calc as phon_calc
    from tqdm import tqdm
    bands = {"DFT": DFT}
    colors = ['k', 'b', 'g', 'r', 'c', 'm', 'y' ]
    style = {"DFT": {"color": colors[0], "lw": 2}}
    
    if args["pots"]:
        for poti, potpath in enumerate(tqdm(args["pots"])):
            potname = path.basename(potpath)
            potkey = potname[:-4]
            pottype = "GAP" if potname[0:2] == "gp" else "MTP"
            bands[potkey] = phon_calc(phondb.atoms, potpath, kpath,
                                      phondb.phonocache, supercell=args["dim"],
                                      Npts=args["npts"], potname=pottype)
            style[potkey] = {"color": colors[1+poti], "lw": 2}

    title = args["title"].format(phondb.atoms.get_chemical_formula())
    savefile = None
    if args["save"]:
        savefile = path.join(phondb.parent.plotdir, args["save"])
                             
    bandplot(bands, names, title=title, outfile=savefile,
             figsize=args["figsize"], style=style, nbands=args["nbands"])

def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    if args is None:
        return

    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    from matdb.database.controller import Controller
    cdb = Controller(args["dbspec"])
    target = cdb.collections[args["s"]]
    
    if args["bands"]:
        from matdb.database.phonon import PhononDFT
        phondb = target.databases[PhononDFT.name]
        phondb.calc_bands()
        _band_plot(phondb, args)        
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())
