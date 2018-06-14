"""Methods for calculating the phonon spectra of materials.
"""
import numpy as np
from hashlib import sha1
from os import path, mkdir, remove
from glob import glob
from tempfile import mkdtemp
from phonopy import Phonopy
from phonopy.file_IO import write_FORCE_CONSTANTS
from shutil import rmtree
from ase.optimize.precon import Exp, PreconLBFGS
from ase.optimize import FIRE

from matdb.base import testmode
from matdb.utility import chdir, redirect_stdout, convert_dict_to_str, execute
from matdb import msg
from matdb.database.hessian import matdb_to_phonopy, unroll_fc, phonopy_to_matdb
from matdb.transforms import conform_supercell
from matdb.atoms import Atoms
from matdb.kpoints import parsed_kpath
    
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

def _calc_quick(atoms, supercell=(1, 1, 1), delta=0.01):
    """Calculates the Hessian for a given atoms object just like :func:`calc`,
    *but*, it uses symmetry to speed up the calculation. Depending on the
    calculator being used, it is possible that the symmetrized result may be
    different from the full result with all displacements, done manually by
    :func:`calc`.

    Args:
        atoms (matdb.Atoms): atomic structure of the *primitive*.
        supercell (list): or `tuple` or :class:`numpy.ndarray` specifying the
          integer supercell matrix.
        delta (float): displacement in Angstroms of each atom when computing the
          phonons. 

    Returns:
        numpy.ndarray: Hessian matrix that has dimension `(natoms*3, natoms*3)`,
        where `natoms` is the number of atoms in the *supercell*.
    """
    #We need to make sure we are at the zero of the potential before
    ratoms = atoms.copy()
    try:
        with open("phonons.log", 'w') as f:
            with redirect_stdout(f):
                precon = None
                if ratoms.n > 100:
                    precon = Exp(A=3)
                #minim = PreconLBFGS(ratoms, precon=precon)
                minim = FIRE(ratoms)
                minim.run(fmax=1e-4)
    except:
        #The potential is unstable probably. Issue a warning.
        msg.warn("Couldn't optimize the atoms object. Potential may be unstable.")
    
    primitive = matdb_to_phonopy(ratoms)
    phonon = Phonopy(primitive, conform_supercell(supercell))
    phonon.generate_displacements(distance=delta)
    supercells = phonon.get_supercells_with_displacements()
    pot = atoms.get_calculator()
    assert pot is not None
    
    forces = []
    for scell in supercells:
        matoms = phonopy_to_matdb(scell)
        #Call a manual reset of the calculator so that we explicitly recalculate
        #the forces for the current atoms object.
        pot.reset()
        matoms.set_calculator(pot)
        forces.append(matoms.get_forces())

    phonon.set_forces(forces)
    phonon.produce_force_constants()
    return unroll_fc(phonon._force_constants)

def calc(primitive, cachedir=None, supercell=(1, 1, 1), delta=0.01, quick=True):
    """Calculates the Hessian for a given atoms object (which *must* have an
    attached calculator).

    .. note:: We choose to use the Hessian as the fundamental quantity in
      vibrational analysis in `matdb`.

    .. note:: `atoms` will be relaxed before calculating the Hessian.

    Args:
        primitive (matdb.Atoms): atomic structure of the *primitive*.
        cachedir (str): path to the directory where phonon calculations are
          cached. If not specified, a temporary directory will be used.
        supercell (tuple): number of times to duplicate the cell when
          forming the supercell.
        delta (float): displacement in Angstroms of each atom when computing the
          phonons. 
        quick (bool): when True, use symmetry to speed up the Hessian
          calculation. See :func:`_calc_quick`.

    Returns:
        numpy.ndarray: Hessian matrix that has dimension `(natoms*3, natoms*3)`,
        where `natoms` is the number of atoms in the *supercell*.
    """
    if quick:
        return _calc_quick(primitive, supercell, delta)
    else:
        atoms = primitive.make_supercell(supercell)
        atoms.set_calculator(primitive.get_calculator())

    from ase.vibrations import Vibrations
        
    #The phonon calculator caches the displacements and force sets for each
    #atomic displacement using pickle. This generates three files for each
    #atomic degree of freedom (one for each cartesian direction). We want to
    #save these in a special directory.
    tempcache = False
    if cachedir is None:
        cachedir = mkdtemp()
        tempcache = True
    else:
        cachedir = path.abspath(path.expanduser(cachedir))
    if not path.isdir(cachedir):
        mkdir(cachedir)

    result = None
    precon = Exp(A=3)
    aphash = None
        
    #Calculate a hash of the calculator and atoms object that we are calculating
    #for. If the potential doesn't have a `to_dict` method, then we ignore the
    #hashing.
    if not tempcache and hasattr(atoms, "to_dict") and hasattr(atoms._calc, "to_dict"):
        atoms_pot = {"atoms": atoms.to_dict(), "pot": atoms._calc.to_dict()}
        #This UUID will probably be different, even if the positions and species
        #are identical.
        del atoms_pot["atoms"]["uuid"]
        hash_str = convert_dict_to_str(atoms_pot)
        aphash = str(sha1(hash_str).hexdigest())

    if not tempcache:
        #Check whether we should clobber the cache or not.
        extras = ["vibsummary.log", "vib.log", "phonons.log"]
        
        with chdir(cachedir):
            hash_match = False
            if path.isfile("atomspot.hash"):
                with open("atomspot.hash") as f:
                    xhash = f.read()
                hash_match = xhash == aphash

            hascache = False
            if not hash_match:
                for vibfile in glob("vib.*.pckl"):
                    remove(vibfile)
                    hascache = True

                for xfile in extras:
                    if path.isfile(xfile):
                        remove(xfile)
                        hascache = True

            if hascache:
                msg.warn("Using hard-coded cache directory. We were unable to "
                         "verify that the atoms-potential combination matches "
                         "the one for which existing cache files exist. So, we "
                         "clobbered the existing files to get the science "
                         "right. You can fix this by using `matdb.Atoms` "
                         "and `matdb.calculators.*Calculator` objects.")
            
    with chdir(cachedir):           
        #Relax the cell before we calculate the Hessian; this gets the forces
        #close to zero before we make harmonic approximation.
        try:
            with open("phonons.log") as f:
                with redirect_stdout(f):
                    minim = PreconLBFGS(atoms, precon=precon, use_armijo=True)
                    minim.run(fmax=1e-5)
        except:
            #The potential is unstable probably. Issue a warning.
            msg.warn("Couldn't optimize the atoms object. Potential may be unstable.")

        vib = Vibrations(atoms, delta=delta)
        with open("vib.log", 'a') as f:
            with redirect_stdout(f):
                vib.run()

        vib.summary(log="vibsummary.log")
        result = vib.H

        #Cache the hash of the atoms object and potential that we were using so
        #that we can check next time whether we should clobber the cache or not.
        if aphash is not None and not tempcache:
            with open(path.join(cachedir, "atomspot.hash"), 'w') as f:
                f.write(aphash)
        
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
