"""Methods for calculating the phonon spectra of materials.
"""
import numpy as np
from os import path

def _ordered_unique(items):
    """Returns the list of unique items in the list while still *preserving* the
    order in which they appeared.
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def _get_calculator(name, potfile):
    """Returns an ASE-compatible calculator to use for the phonon calculations.

    Args:
        name (str): name of the potential type; valid options are ["GAP",
          "MTP"].
        potfile (str): full path to the potential file that the potential should
          use. For GAP, this is a `.xml` file; for MTP, it is a `.ini`.
    """
    from os import getcwd, chdir
    from quippy.potential import Potential
    
    #GAP doesn't like relative paths when their are multiple binary files
    #representing the various inter-atomic potentials. We want to switch to the
    #phonon directory so that the phonon files are cached in the right place;
    #but first we need to change into the potential directory.
    curdir = getcwd()
    potdir = path.dirname(potfile)
    chdir(potdir)

    try:
        if name.lower() == "mtp":
            from lammpslib import LAMMPSlib
            otypes = _ordered_unique(AgPd.get_chemical_symbols())
            atom_types = {s: i+1 for i, s in enumerate(otypes)}
            header = ["units metal",
                      "dimension 3",
                      "boundary p p p",
                      "atom_style atomic",
                      "atom_modify map array"]
            cmds = ["pair_style MLIP {}".format(potfile),
                    "pair_coeff * *",
                    "mass 1 106.42", #Pd
                    "mass 2 107.87", #Ag
                    "neighbor 2.0 bin",
                    "neigh_modify delay 10 check yes"]
            pot = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types,
                            lammps_header=header, keep_alive=True)
        else:
            pot = Potential("IP GAP", param_filename=potfile)
    finally:
        #If there is an initialization error in either GAP or MTP, we don't want
        #the current directory to be different than it was before we came into
        #this method.
        chdir(curdir)

    return pot

def calc(atoms, fit, kpath, cachedir, supercell=4, delta=0.01, Npts=100,
         potname="GAP"):
    """Calculates the phonon band structure for the specified QUIP IP parameter
    file.

    .. warning:: If the delta given to this routine does *not* match the delta
      used in the phonon calculations (`phonopy` uses the
      `DISPLACEMENT_DISTANCE` tag to set this value). Since the `phonopy`
      default is 0.01 A, we use that same value here.
    
    Args:
        atoms (quippy.Atoms): atomic structure information.
        fit (matdb.fitting.basic.Trainer): potential to run the calculation for.
        kpath (list): of :class:`numpy.ndarray` specifying the points in k-space
          that should be visited.
        cachedir (str): path to the directory where phonon calculations are
          cached. The cache is sorted by system (using chemical formula), *and*
          by the name of the potential being used.
        supercell (int or list): number of times to duplicate the cell when
          forming the supercell; if an `int` is given, the cell is duplicated
          equally in all three dimensions.
        delta (float): displacement in Angstroms of each atom when computing the
          phonons. 
        Npts (int): number of points to sample across the k-path.
        potname (str): one of ["MTP", "GAP"] specifying whether to use the
          `lammpslib` and MTP :class:`ase.Calculator` to interface with LAMMPS
          instead of a `quippy` GAP potential.

    Returns:
        tuple: `(q, omega, path_k, Q)` where `q` is the 1D path vector for plotting;
        `omega` is the matrix of phonon modes; `path_k` is the matrix of actual
        k-points that were sampled (corresponding to `q`); `Q` is the list of
        special points (in `q`) corresponding to the given `kpath`.
    """
    from matdb.utility import reporoot
    from os import path, mkdir, getcwd, chdir
    #The phonon calculator caches the displacements and force sets for each
    #atomic displacement using pickle. This generates three files for each
    #atomic degree of freedom (one for each cartesian direction). We want to
    #save these in a special directory.
    pot = fit.calculator
    dirname = path.join(cachedir, fit.fqn)
    if not path.isdir(dirname):
        mkdir(dirname)
    curdir = getcwd()

    from ase.optimize.precon import Exp, PreconLBFGS
    from matdb.utility import redirect_stdout
    try:
        chdir(dirname)
            
        #Relax the cell before we try and get phonons out.
        atoms.set_calculator(pot)
        atoms.set_cutoff(pot.cutoff())
        minim = PreconLBFGS(atoms, precon=Exp(A=3.0))
        minim.run(fmax=1e-7)

        #Now run the phonon calculation
        from ase.phonons import Phonons
        if isinstance(supercell, int):
            supercell = [supercell]*3
        ph = Phonons(atoms, pot, supercell=tuple(supercell), delta=delta)

        from sys import stdout
        with open("phonons.log", 'w') as f:
            with redirect_stdout(f):
                ph.run()

                #Read forces and assemble the dynamical matrix
                ph.read(acoustic=True)

        #Conversion from eV to THz. We want to get out the phonon modes for the
        #full path. They are in units of eV, which we want to convert to THz for
        #comparison with the published values. We use the speed of light $c$:

        #$c = f \lambda$ is the simple relation for wave velocity vs. frequency
        #and wavelength; we combine that with the photon energy relation $E =
        #\frac{h c}{\lambda}$ to obtain $E = h f$.

        #Simplifying, we obtain: $f = \frac{E}{h}$
        #The value for $h$ in units of `eV.s` is $h = 0.004135668 \times
        #10^{-12}$ where the $10^{-12} gives us `eV/THz`.
        h = 0.004135668
        import ase.dft.kpoints as kpt
        path_kc, q, Q = kpt.get_bandpath(kpath, atoms.cell, Npts)
        omega = ph.band_structure(path_kc)/h #Now in units of THz
    finally:
        #Switch the directory back to what it was.
        chdir(curdir)

    result = {
        "q": q,
        "w": omega,
        "path": path_kc,
        "Q": Q
    }
    return result

def from_yaml(filepath):
    """Converts the `band.yaml` file generated by phonopy into a 
    format we can plot against the IPs
    
    Args:
        filepath (str): path to the `band.yaml` file.
            
    Returns:
        tuple: (q, w) where `q` is the distance along the special path 
        and `w` is the array of band values at that distance.
    """
    import yaml
    with open(filepath, 'r') as stream:
        dbands = yaml.load(stream)
    q = np.array([e["distance"] for e in dbands["phonon"]])
    w = np.array([[k["frequency"] for k in e["band"]]
                     for e in dbands["phonon"]])
    path = np.array([e["q-position"] for e in dbands["phonon"]])
    Q = np.array([e["distance"] for e in dbands["phonon"]
                  if "label" in e])
    
    result = {
        "q": q,
        "w": w,
        "path": path,
        "Q": Q
    }
    return result

def bandplot(phonons, names, nbands=None, style=None, figsize=(8, 6),
             ptype=None, title="Phonon Spectrum", outfile=None):
    """Plots the phonon spectrum along the given path.

    Args:
        phonons (dict): keys are arbitrary ids, though they should
          correspond to `style` and `ptype`; values are `dict`
          instances returned by either :func:`calc` or
          :func:`from_yaml`. 
        names (list): of `str`; specifies the labels for each of the points in `Q`.
        nbands (int): number of phonon bands (dispersion curves) to
          plot. If not specified, all of them are used.
        style (dict): keys are same as for `phonons`; values are the line/scatter
          styles to use for the plotting, given as a keyword argument dictionary that
          is passed directly to the plotting method.
        ptype (dict): keys are same as for `q` an `w`; values are one of 
          ["plot", "scatter"], specifying whether to produce a line plot or scatter
          plot. If not specified, then scatter plots are used by
          default.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    plt.figure(1, figsize, frameon=False)
    plt.axes([.1, .07, .67, .85])

    first = list(phonons.keys())[0]
    maxq = np.max(phonons[first]["q"])
    minq = np.min(phonons[first]["q"])

    Q = np.unique(phonons[first]["Q"])
    plt.xlim(minq-.1, maxq+.1)
    plt.xticks(Q, names, fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=22)
    plt.grid('on')
    
    for label, pdict in phonons.items():
        omega = pdict["w"]
        qs = pdict["q"].copy()
        scaling = (maxq-minq)/(np.max(qs)-np.min(qs))
        qs = scaling*(qs-np.min(qs))
        
        if len(omega.shape) == 1:
            omega.shape = (len(omega), 1)
        if nbands is None:
            nbands = omega.shape[1]

        for n in range(min((nbands, omega.shape[1]))):
            omega_n = omega[:,n]
            
            if style is not None and label in style:
                kwargs = style[label]
            else:
                #The default dictionary is for the line plots.
                kwargs = {"lw": 2}
            
            #We only want the label to show up once in the legend.
            if ptype is not None and label in ptype:
                plotter = getattr(plt, ptype[label])
            else:
                plotter = plt.plot
                
            if n == 0:
                plotter(qs, omega_n, label=label, **kwargs)
            else:        
                plotter(qs, omega_n, **kwargs)

    #Just use the first dict's special points.
    plt.title(title, fontsize=24)
    plt.legend(loc=8)
    if outfile is not None:
        plt.savefig(outfile)
    else:
        plt.show()
