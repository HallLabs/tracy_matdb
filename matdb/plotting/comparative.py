"""`matdb` generates many databases and potentials en route to the final
product. In order to adjust parameters it is useful to plot potentials and
convergence runs against each other.
"""
from os import path
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

from matdb.atoms import Atoms
from matdb.phonons import bandplot
from matdb.utility import chdir
from matdb.phonons import from_yaml, calc as phon_calc, _calc_bands
from matdb.kpoints import parsed_kpath
from matdb import msg
from matdb.transforms import conform_supercell

def band_plot(dbs, fits=None, npts=100, title="{} Phonon Spectrum", save=None,
              figsize=(10, 8), nbands=12, delta=0.01, **kwargs):
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
        kwargs (dict): additional "dummy" arguments so that this method can be
          called with arguments to other functions.
    """
    if len(dbs) == 0:
        raise ValueError("Must specify at least one Hessian group "
                         "to plot phonons for.")

    db = dbs[0]
    title = title.format(db.atoms.get_chemical_formula())
    nlines = len(dbs) + (0 if fits is None else len(fits))
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, nlines))
    
    bands, style = {}, {}
    names, kpath = parsed_kpath(db.atoms)
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
            #make_supercell already returns a copy of the atoms object.
            scell = conform_supercell(db.supercell)
            ai = db.atoms.make_supercell(scell)
            ai.set_calculator(fit.calculator)
            #Note that for these kinds of calculations, we don't want to use the
            #default cache directory for the potential; rather use a temporary
            #directory (default when None is specified).
            H = phon_calc(ai, None, delta)
            bands[fit.fqn] = _calc_bands(db.atoms, H, scell)
            style[fit.fqn] = {"color": colors[gi], "lw": 2}

    savefile = None
    if save:
        savefile = path.join(db.database.parent.plotdir, save)
                             
    bandplot(bands, names, title=title, outfile=savefile,
             figsize=figsize, style=style, nbands=nbands)

def band_raw(poscar, bandfiles, pots=None, supercell=None, npts=100,
             title="{} Phonon Spectrum", save=None, figsize=(10, 8), nbands=4,
             line_names=None, delta=0.01, **kwargs):
    """Plots the phonon bands from raw `band.yaml` files.

    Args:
        poscar (str): path to the POSCAR file for the *primitive* to plot bands for.
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
        kwargs (dict): additional "dummy" arguments so that this method can be
          called with arguments to other functions.
    """
    nlines = len(bandfiles) + 0 if pots is None else len(pots)
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, nlines))        
    bands, style = {}, {}
    atoms = Atoms(poscar, format="vasp")
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
        for fiti, potfile in enumerate(tqdm(pots)):
            gi = len(bandfiles) + fiti
            aprim = atoms.copy()
            ai = aprim.make_supercell(supercell)
            potdir, potname = path.split(potfile)
            with chdir(potdir):
                fit = quippy.Potential("IP GAP", param_filename=potname)
                cachedir = path.join(potdir, "cache")
                ai.set_calculator(fit)
                
            Hess = phon_calc(ai, cachedir, delta)
            bands[line_names[gi]] = _calc_bands(aprim, Hess, supercell)
            style[line_names[gi]] = {"color": colors[gi], "lw": 2}
        
    title = title.format(atoms.get_chemical_formula())
    savefile = None
    if save:
        savefile = save
                             
    bandplot(bands, names, title=title, outfile=savefile,
             figsize=figsize, style=style, nbands=nbands)
    
