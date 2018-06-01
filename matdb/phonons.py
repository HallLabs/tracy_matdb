"""Methods for calculating the phonon spectra of materials.
"""
import numpy as np
from os import path, mkdir
from matdb.base import testmode
from matdb.utility import chdir, redirect_stdout
from matdb import msg

def _ordered_unique(items):
    """Returns the list of unique items in the list while still *preserving* the
    order in which they appeared.
    """
    seen = set()
    seen_add = seen.add
    return [x for x in items if not (x in seen or seen_add(x))]

def roll(hessian):
    """Rolls the specified hessian into the `phonopy` format.
    
    Args:
        hessian (numpy.ndarray): of shape `n_atoms * 3`.
    """
    n = hessian.shape[0]/3
    result = np.zeros((n, n, 3, 3), dtype='double')
    
    for i in range(n):
        for j in range(n):
            result[i, j] = hessian[i*3:(i+1)*3, j*3:(j+1)*3]

    return result

def _calc_bands(atoms, hessian, supercell=(1, 1, 1), outfile=None, grid=None):
    """Calculates the band structure for the given Hessian matrix.

    Args:
        atoms (matdb.Atoms): atoms object corresponding to the *primitive*
          cell. The specified supercell matrix should  result in a number
          of atoms that matches the dimensionality of the Hessian.
        supercell (tuple): of `int` supercell matrix components; can have either
          3 or 9 components.
        hessian (numpy.ndarray): with shape `(natoms*3, natoms*3)`.
        grid (list): of `int` specifying the number of divisions in k-space
          along each reciprocal unit vector.
        outfile (str): path to the output `band.yaml` file that should be
          created by this function.

    Returns:
        If `outfile` is None, then this method returns a dictionary that has the
        same format as :func:`from_yaml`.
    """
    from tempfile import mkdtemp
    from phonopy.file_IO import write_FORCE_CONSTANTS
    from matdb.kpoints import parsed_kpath
    from matdb.utility import execute
    from shutil import rmtree
    
    #Create a temporary directory in which to work.
    target = mkdtemp()
    bandfile = path.join(target, "band.yaml")
    
    if grid is None:
        grid = [13, 13, 13]
    if isinstance(supercell, np.ndarray):
        supercell = supercell.flatten()
    
    #First, roll up the Hessian and write it as a FORCE_CONSTANTS file.
    with chdir(target):
        HR = roll(hessian)
        write_FORCE_CONSTANTS(HR)
        atoms.write("POSCAR", format="vasp")

        #We need to create the band.conf file and write the special
        #paths in k-space at which the phonons should be calculated.
        atom_types = _ordered_unique(atoms.get_chemical_symbols())
        settings = [
            ("FORCE_CONSTANTS", "READ"),
            ("ATOM_NAME", ' '.join(atom_types)),
            ("DIM", ' '.join(map(str, supercell))),
            ("MP", ' '.join(map(str, grid)))
        ]

        labels, bands = parsed_kpath(atoms)
        bandfmt = "{0:.3f} {1:.3f} {2:.3f}"
        sband = []        
        for Q in bands:
            sband.append(bandfmt.format(*Q))

        settings.append(("BAND", "  ".join(sband)))
        settings.append(("BAND_LABELS", ' '.join(labels)))

        with open("band.conf", 'w') as f:
            for k, v in settings:
                f.write("{} = {}\n".format(k, v))

    sargs = ["phonopy", "band.conf"]
    xres = execute(sargs, target, venv=True)

    if not path.isfile(bandfile): #pragma: no cover
        msg.err("could not calculate phonon bands; see errors.")
        msg.std(''.join(xres["output"]))

    result = None
    if outfile is not None:
        #Move the band.yaml file to the new target location.
        from shutil import move
        move(bandfile, outfile)
    else:
        result = from_yaml(bandfile)

    #Remove the temporary directory that we created and return the result.
    rmtree(target)
    return result            

def calc(atoms, cachedir=None, delta=0.01):
    """Calculates the Hessian for a given atoms object (which *must* have an
    attached calculator).

    .. note:: We choose to use the Hessian as the fundamental quantity in
      vibrational analysis in `matdb`.

    .. note:: `atoms` will be relaxed before calculating the Hessian.

    Args:
        atoms (quippy.Atoms): atomic structure of the *supercell*.
        cachedir (str): path to the directory where phonon calculations are
          cached. If not specified, a temporary directory will be used.
        supercell (tuple): number of times to duplicate the cell when
          forming the supercell.
        delta (float): displacement in Angstroms of each atom when computing the
          phonons. 

    Returns:
        numpy.ndarray: Hessian matrix that has dimension `(natoms*3, natoms*3)`,
        where `natoms` is the number of atoms in the *supercell*.
    """
    from ase.optimize.precon import Exp, PreconLBFGS
    from ase.vibrations import Vibrations
        
    #The phonon calculator caches the displacements and force sets for each
    #atomic displacement using pickle. This generates three files for each
    #atomic degree of freedom (one for each cartesian direction). We want to
    #save these in a special directory.
    if cachedir is None:
        from tempfile import mkdtemp
        cachedir = mkdtemp()
    else:
        cachedir = path.abspath(path.expanduser(cachedir))
    if not path.isdir(cachedir):
        mkdir(cachedir)

    result = None
    precon = Exp(A=3)
    with chdir(cachedir):           
        #Relax the cell before we calculate the Hessian; this gets the forces
        #close to zero before we make harmonic approximation.
        try:
            with redirect_stdout("phonons.log"):
                minim = PreconLBFGS(atoms, precon=precon, use_armijo=True,
                                    logfile="phonons.log")
                minim.run(fmax=1e-5)
        except:
            #The potential is unstable probably. Issue a warning.
            msg.warn("Couldn't optimize the atoms object. Potential may be unstable.")

        vib = Vibrations(atoms, delta=delta)
        with open("phonons.log", 'a') as f:
            with redirect_stdout(f):
                vib.run()

        #Read forces and assemble the Hessian matrix.
        vib.summary(log="vibsummary.log")
        result = vib.H

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
          instances returned by :func:`from_yaml`. 
        names (list): of `str`; specifies the labels for each of the points in `Q`.
        nbands (int): number of phonon bands (dispersion curves) to
          plot. If not specified, all of them are used.
        style (dict): keys are same as for `phonons`; values are the line/scatter
          styles to use for the plotting, given as a keyword argument dictionary that
          is passed directly to the plotting method.
        ptype (dict): keys are same as for `phonons` an `style`; values are one of 
          ["plot", "scatter"], specifying whether to produce a line plot or scatter
          plot. If not specified, then scatter plots are used by
          default.
        figsize (tuple): specifying the width and height (in inches) of the plot
          that will be produced.
        title (str): plot title override; used as-is.
        outfile (str): if the plot should be saved, the path to the output
          file. Output format is detected from file extension.
    """
    import matplotlib.pyplot as plt
    if testmode:
        import matplotlib
        matplotlib.use('Agg')

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
        qs = pdict["q"]
        
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
