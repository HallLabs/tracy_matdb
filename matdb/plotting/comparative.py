"""`matdb` generates many databases and potentials en route to the final
product. In order to adjust parameters it is useful to plot potentials and
convergence runs against each other.
"""
def band_plot(phondbs, fits=None, dim=2, npts=100, title="{} Phonon Spectrum",
              save=None, figsize=(10, 8), nbands=4, **kwargs):
    """Plots the phonon bands for the specified CLI args.

    Args:
        phondbs (list): of :class:`matdb.database.phonon.DynMatrix` `phonopy`
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
        nbands (int): number of bands to plot.
        kwargs (dict): additional "dummy" arguments so that this method can be
          called with arguments to other functions.
    """
    from matdb.phonons import bandplot
    from os import path
    from matdb.phonons import calc as phon_calc
    from tqdm import tqdm    

    #Make sure we have bands calculated for each of the databases passed in.
    for phondb in phondbs:
        phondb.calc_bands()

    colors = ['k', 'b', 'g', 'r', 'c', 'm', 'y' ]
    bands, style = {}, {}
    names, kpath = None, None
    
    for dbi, phondb in enumerate(phondbs):
        if names is None:
            names, kpath = phondb.kpath
            #matplotlib needs the $ signs for latex if we are using special
            #characters. We only get names out from the first configuration; all
            #the others have to use the same one.
            names = ["${}$".format(n) if '\\' in n else n for n in names]
        bands[phondb.parent.name] = phondb.bands
        style[phondb.parent.name] = {"color": colors[dbi], "lw": 2}

    #All of the phonon calculations use the same base atoms configuration. The
    #last `phondb` in the enumerated list is as good as any other.
    if fits is not None:
        for fiti, fit in enumerate(tqdm(fits)):
            bands[fit.fqn] = phon_calc(phondb.atoms, fit, kpath,
                                      phondb.phonocache, supercell=dim,
                                      Npts=npts, potname=fit.fqn)
            style[fit.fqn] = {"color": colors[len(phondbs)+fiti], "lw": 2}

    title = title.format(phondb.atoms.get_chemical_formula())
    savefile = None
    if save:
        savefile = path.join(phondb.parent.plotdir, save)
                             
    bandplot(bands, names, title=title, outfile=savefile,
             figsize=figsize, style=style, nbands=nbands)
