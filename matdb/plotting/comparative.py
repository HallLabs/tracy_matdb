"""`matdb` generates many databases and potentials en route to the final
product. In order to adjust parameters it is useful to plot potentials and
convergence runs against each other.
"""
from os import path

from ase.build import make_supercell
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from matdb import msg
from matdb.atoms import Atoms
from matdb.calculators import build_calc
from matdb.kpoints import parsed_kpath
from matdb.phonons import bandplot
from matdb.phonons import from_yaml, _calc_bands, calc as phon_calc
from matdb.transforms import conform_supercell
from matdb.utility import chdir

def band_plot(dbs, fits=None, npts=100, title="{} Phonon Spectrum", save=None,
              figsize=(10, 8), nbands=None, delta=0.01, quick=True, **kwargs):
    """Plots the phonon bands for the specified CLI args.

    Args:
        dbs (list): of :class:`matdb.database.hessian.Hessian` `phonopy`
          calculation database instances that have DFT-accurate band
          information.
        fits (list): of :class:`~matdb.fitting.basic.Trainer` to calculate bands
          for.
        dim (list): of `int`; supercell dimensions for the phonon calculations.
        npts (int): number of points to sample along the special path in
          k-space.
        title (str): Override the default title for plotting; use `{}` for
          formatting chemical formula.
        save (str): name of a file to save the plot to; otherwise the plot is
          shown in a window.
        figsize (tuple): of `float`; the size of the figure in inches.
        delta (float): size of displacement for finite difference derivative.
        nbands (int): number of bands to plot.
        quick (bool): when True, use symmetry to speed up the Hessian
          calculation for the specified potentials.
        kwargs (dict): additional "dummy" arguments so that this method can be
          called with arguments to other functions.
    """
    if len(dbs) == 0:
        raise ValueError("Must specify at least one Hessian group "
                         "to plot phonons for.")

    db = dbs[0]
    if db.atoms is not None:
        ratoms = db.atoms
    else:
        #Use the parent group; this one must have been one of the sub-sequence
        #recursively generated groups.
        ratoms = db.parent.atoms
    title = title.format(ratoms.get_chemical_formula())

    nlines = len(dbs) + (0 if fits is None else len(fits))
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, nlines))

    bands, style = {}, {}
    names, kpath = parsed_kpath(ratoms)
    #matplotlib needs the $ signs for latex if we are using special
    #characters. We only get names out from the first configuration; all
    #the others have to use the same one.
    names = ["${}$".format(n) if '\\' in n or '_' in n
             else n for n in names]

    for dbi, db in enumerate(dbs):
        db.calc_bands()
        bands[db.key] = db.bands
        style[db.key] = {"color": colors[dbi], "lw": 2}

    #All of the phonon calculations use the same base atoms configuration. The
    #last `phondb` in the enumerated list is as good as any other.
    if fits is not None:
        for fiti, fit in enumerate(tqdm(fits)):
            gi = len(dbs) + fiti
            ai = ratoms.copy()
            ai.set_calculator(fit.calculator)
            H = phon_calc(ai, supercell=db.supercell, delta=delta, quick=quick)
            bands[fit.fqn] = _calc_bands(db.atoms, H, db.supercell)
            style[fit.fqn] = {"color": colors[gi], "lw": 2}

    savefile = None
    if save:
        savefile = path.join(db.database.parent.plotdir, save)

    bandplot(bands, names, title=title, outfile=savefile,
             figsize=figsize, style=style, nbands=nbands)

def band_raw(primitive, bandfiles=None, pots=None, supercell=None, npts=100,
             title="{} Phonon Spectrum", save=None, figsize=(10, 8), nbands=4,
             line_names=None, delta=0.01, quick=True, **kwargs):
    """Plots the phonon bands from raw `band.yaml` files.

    Args:
        primitive (str): path to the atoms file for the *primitive* to plot
          bands for. Use the ASE format string as a prefix,
          e.g. `vasp-xml:vasprun.xml` or `extxyz:atoms.xyz`. Default assumes
          `vasp:{}` if no format is specified.
        bandfiles (list): of `str` file paths to the plain `band.yaml` files.
        supercell (list): of `int`; supercell dimensions for the phonon
          calculations.
        npts (int): number of points to sample along the special path in
          k-space.
        title (str): Override the default title for plotting; use `{}` for
          formatting chemical formula.
        save (str): name of a file to save the plot to; otherwise the plot is
          shown in a window.
        figsize (tuple): of `float`; the size of the figure in inches.
        nbands (int): number of bands to plot.
        delta (float): size of displacement for finite difference derivative.
        quick (bool): when True, use symmetry to speed up the Hessian
          calculation for the specified potentials.
        kwargs (dict): additional "dummy" arguments so that this method can be
          called with arguments to other functions.
    """
    nlines = len(bandfiles) + (0 if pots is None else len(pots))
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, nlines))
    bands, style = {}, {}

    #Handle DSL format import on the file path for the primitive cell.
    if ':' not in primitive:
        atoms = Atoms(primitive, format="vasp")
    else:
        fmt, atpath = primitive.split(':')
        atoms = Atoms(atpath, format=fmt)

    names, kpath = parsed_kpath(atoms)

    #matplotlib needs the $ signs for latex if we are using special
    #characters. We only get names out from the first configuration; all
    #the others have to use the same one.
    names = ["${}$".format(n) if '\\' in n or '_' in n
             else n for n in names]

    for ifile, fpath in enumerate(bandfiles):
        if line_names is not None:
            key = line_names[ifile]
        else:
            key = "File {}".format(ifile)

        bands[key] = from_yaml(fpath)
        style[key] = {"color": colors[ifile], "lw": 2}

    if pots is not None:
        for fiti, fit in enumerate(tqdm(pots)):
            gi = len(bandfiles) + fiti
            atoms.set_calculator(fit)
            H = phon_calc(atoms, supercell=supercell, delta=delta, quick=quick)
            bands[line_names[gi]] = _calc_bands(atoms, H, supercell)
            style[line_names[gi]] = {"color": colors[gi], "lw": 2}

    title = title.format(atoms.get_chemical_formula())
    savefile = None
    if save:
        savefile = save

    bandplot(bands, names, title=title, outfile=savefile,
             figsize=figsize, style=style, nbands=nbands)
