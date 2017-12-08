"""Functions for plotting interatomic potential performance vs. correct (DFT)
answers within a particular :class:`quippy.AtomsList`.
"""
import numpy as np
import quippy
from matdb.plotting.matd3 import PointDetailImage as PDI

def generate(plots, pot, atoms, folder=None, base64=False, index=0, ndimer=50,
             ntrimer=75):
    """Generates a set of plots for the given potential and atoms list. The
    following are possible options; for each one, a shorthand character is
    specified:

    1. Energy Correlation Plot (`e`)
    2. Force Correlation Plot (`f`)
    3. Virial Correlation Plot (`v`)
    4. Energy vs. Volume Plot (`o`)
    5. Dimer Plots for 2-body Potentials (`d`)
    6. Trimer Plots for 3-body Potentials (`t`)

    You can specify which plots to include in the set by adding the character
    codes to the `plots` argument. For example, to plot dimer, trimer and energy
    vs. volume plots, use `plots="dte"`.

    .. note:: The order that the `plots` are specified in is the order that they
      will shop up in the interactive plot.

    Args:
        plots (str): of character codes as described above.
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        seed_index (int): integer index to use for the first plot. It is
          incremented automatically for each plot added.
        ndimer (int): number of samples to take along the dimer trajectory.
        ntrimer (int): number of samples to take along the trimer angle
          trajectory. 
    """
    #We need to determine the unique list of elements in the atoms list. Lets
    elements = set()
    for a in atoms:
        elements |= set(a.get_chemical_symbols())
    
    plotmap = {
        'e': lambda i: energy(pot, atoms, folder, base64, i),
        'f': lambda i: force(pot, atoms, folder, base64, i),
        'v': lambda i: virial(pot, atoms, folder, base64, i),
        'o': lambda i: EvsV(pot, atoms, folder, base64, i),
        'd': lambda i: dimer(pot, atoms, elements, folder, base64, i, ndimer),
        'z': lambda i: dimer(pot, atoms, elements, folder, base64, i, ndimer, zoom=True),
        't': lambda i: trimer(pot, atoms, elements, folder, base64, i, ndimer)
    }

    result = {}
    for char in plots:
        image = plotmap[char](index)
        if char in ['d', 't', 'z']:
            for k, img in image.items():
                result[k] = img
        else:
            result[image.imgtype] = image

    return result

def trimer(pot, atoms, elements, folder=None, base64=False, index=None,
           nsamples=50):
    """Plots the potential behavior as the angle between three atoms changes in a
    trimer.

    .. note:: This produces a *dict* of all possible 3-body interactions between
      distinct element types. 

    Args:
        elements (list): of chemical symbols in the system.
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent
          collection.
        nsamples (int): number of samples to take along the trajectory.

    Returns:

        dict: keys are concatenated element names (e.g. `AgPd`) and values are
        :class:`PointDetailImage` instances.
    """
    #First, determine the set of unique element pairs in the list.
    from itertools import product
    possible = list(product(elements, repeat=3))
    target = list(set(map(lambda l: tuple(sorted(l)), possible)))

    subplot_kw = {
        "xlabel": "Angle (Rad)",
        "ylabel": "IP Energy (eV)"
    }
    
    result = {}
    for elements in target:
        #Look up the lattice parameters for each of the elements. Use Vegard's
        #law to get a decent domain to plot the energy over.
        from matdb.data import vegard
        uelements = list(set(elements))
        concdict = [elements.count(e)/3. for e in uelements]
        rvegard = vegard(uelements, concdict)

        trimer = quippy.Atoms(n=3, lattice=np.eye(3)*100)
        trimer.set_cutoff(pot.cutoff())
        trimer.pos = 0.
        trimer.set_chemical_symbols(elements)
        #Set the position of the second atom to be at equilibrium with respect
        #to the Vegard's law calculation.
        trimer.pos[1,2] = rvegard

        #Now, vary the angle of the third atom with respect to the original two
        #and see how the angle changes.
        theta = np.linspace(np.pi/6, np.pi, nsamples)
        energy = []
        for t in theta:
            x = rvegard*np.cos(t)
            y = rvegard*np.sin(t)
            trimer.pos[1,3] = x
            trimer.pos[2,3] = y
            pot.calc(trimer,energy=True)
            energy.append(trimer.energy)

        elemstr = trimer.get_chemical_formula()
        img = PDI(theta, np.array(energy), "plot", subplot_kw=subplot_kw,
                  index=index, base64=base64, name=elemstr, imgtype=elemstr,
                  folder=folder)
        result[elemstr] = img

    return result        

def _dimer_range(elements):
    """Calculates a reasonable range over which to vary distance in a dimer.

    Args:
        elements (list): of chemical symbols in the system.

    Returns:
        tuple: of `(rmin, rvegard, rmax)` where `rmin` and `rmax` are the minimum
        and maximum lattice parameters for the *pure* elements and `rvegard` is the
        linearly interpolated lattice parameter for the *alloy* using Vegard's law.
    """
    from matdb.data import vegard, latpars
    #For a dimer, the concentrations are always equal.
    concs = [0.5, 0.5]
    alloy = vegard(elements, concs)
    lats = [latpars[e] for e in elements]
    rmin, rmax = min(lats), max(lats)
    return (rmin, alloy, rmax)

def dimer(pot, atoms, elements, folder=None, base64=False, index=None,
          nsamples=50, zoom=False):
    """Plots the potential behavior as the distance between two atoms in a
    dimer.

    .. note:: This produces a *dict* of all possible 2-body interaction between
      distinct element types. 

    Args:
        elements (list): of chemical symbols in the system.
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent
          collection.
        nsamples (int): number of samples to take along the trajectory.
        zoom (bool): when True, show a zoomed in view of the potential.

    Returns:

        dict: keys are concatenated element names (e.g. `AgPd`) and values are
        :class:`PointDetailImage` instances.
    """
    #First, determine the set of unique element pairs in the list.
    from itertools import product
    possible = list(product(elements, repeat=2))
    target = list(set(map(lambda l: tuple(sorted(l)), possible)))

    subplot_kw = {
        "xlabel": "Distance (Ang)",
        "ylabel": "IP Energy (eV)"
    }
    
    result = {}
    for elements in target:
        #Look up the lattice parameters for each of the elements. Use Vegard's
        #law to get a decent domain to plot the energy over.
        rmin, rvegard, rmax = _dimer_range(elements)

        dimer = quippy.Atoms(n=2, lattice=np.eye(3)*100)
        dimer.set_cutoff(pot.cutoff())
        dimer.pos = 0.
        dimer.set_chemical_symbols(elements)

        if zoom:
            rs = np.linspace(0.7*rmin, pot.cutoff(), nsamples)
        else:
            rs = np.linspace(0.4*rmin, pot.cutoff(), nsamples)
            
        energy = []
        for r in rs:
            dimer.pos[1,2] = r
            pot.calc(dimer, energy=True)
            energy.append(dimer.energy)

        elemstr = ''.join(elements) + ('z' if zoom else "")
        img = PDI(rs, np.array(energy), "plot", subplot_kw=subplot_kw,
                  index=index, base64=base64, name=elemstr, imgtype=elemstr,
                  folder=folder)
        result[elemstr] = img

    return result        

def _get_xy(pot, atoms, prop, peratom=False, **calcargs):
    """Gets the x and y values for the specified potential and property.
    
    Args:
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        prop (str): name of the property on each atoms object.
        peratom (bool): when True, plot per atom quantities.
        calcargs (dict): key-value pairs passed to
          :meth:`~quippy.Potential.calc`.

    Returns:
        tuple: of `(dft, pot)` values.
    """
    from tqdm import tqdm
    if peratom:
        dft = [getattr(a, prop)/float(a.n) for a in atoms]
    else:
        dft = getattr(atoms, prop)
    ipprop = next(calcargs.iterkeys())
    ip = []
    
    for a in tqdm(atoms):
        a.set_cutoff(pot.cutoff())
        a.calc_connect()
        pot.calc(a, **calcargs)
        if peratom:
            ip.append(getattr(a, ipprop)/float(a.n))
        else:
            ip.append(getattr(a, ipprop))
            
    return (dft, ip)

def energy(pot, atoms, folder=None, base64=False, index=None):
    """Produces an energy correlation plot for the specified potential.

    Args:
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    dft, ip = _get_xy(pot, atoms, "dft_energy", energy=True, peratom=True)
    dft, ip = np.array(dft), np.array(ip)
    
    #Setup the keyword arguments for the axes labels, etc.
    subplot_kw = {
        "xlabel": "DFT Energy (eV)",
        "ylabel": "IP Energy (eV)"
    }
    if hasattr(atoms[0],"config_type"):
        ctypes = np.array(atoms.config_type)
    else:
        ctypes = None
    title = "RMSE {0:.4f}".format(np.std(dft-ip))
    return PDI(dft, ip, subplot_kw=subplot_kw, index=index, base64=base64,
               name="Energy", imgtype="energy", folder=folder, withdiag=True,
               title=title, partition=ctypes)

def force(pot, atoms, folder=None, base64=False, index=None):
    """Produces a force correlation plot for the specified potential.

    Args:
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    dft, ip = _get_xy(pot, atoms, "dft_force", force=True)
    dft = np.hstack([np.array(f).flatten() for f in dft]).flatten()
    ip = np.hstack([np.array(f).flatten() for f in ip]).flatten()

    if hasattr(atoms[0],"config_type"):
        aconfigs = atoms.config_type
        N = np.array(atoms.n)*3
        ctypes = np.array([ac for ac, n in zip(aconfigs, N) for i in range(n)])
        assert len(ctypes) == len(dft)
    else:
        ctypes = None
    
    #Setup the keyword arguments for the axes labels, etc.
    subplot_kw = {
        "xlabel": "DFT Forces (eV/A)",
        "ylabel": "IP Forces (eV/A)"
    }
    title = "RMSE {0:.4f}".format(np.std(dft-ip))
    return PDI(dft.flatten(), ip.flatten(), subplot_kw=subplot_kw, index=index,
               base64=base64, name="Force", imgtype="force",
               folder=folder, withdiag=True, title=title, partition=ctypes)

def virial(pot, atoms, folder=None, base64=False, index=None):
    """Produces a virial correlation plot for the specified potential.

    Args:
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    dft, ip = _get_xy(pot, atoms, "dft_virial", virial=True, peratom=True)
    dft = np.hstack([np.array(v).flatten() for v in dft]).flatten()
    ip = np.hstack([np.array(v).flatten() for v in ip]).flatten()

    if hasattr(atoms[0],"config_type"):
        aconfigs = atoms.config_type
        ctypes = np.array([ac for ac in aconfigs for i in range(9)])
        assert len(ctypes) == len(dft)
    else:
        ctypes = None
    
    #Setup the keyword arguments for the axes labels, etc.
    subplot_kw = {
        "xlabel": "DFT Virial (eV/A^3)",
        "ylabel": "IP Virial (eV/A^3)"
    }
    title = "RMSE {0:.4f}".format(np.std(dft-ip))
    return PDI(dft.flatten(), ip.flatten(), subplot_kw=subplot_kw, index=index,
               base64=base64, name="Virial", imgtype="virial",
               folder=folder, withdiag=True, title=title, partition=ctypes)

def EvsV(pot, atoms, folder=None, base64=False, index=None):
    """Produces an energy vs. volume plot for the energies predicted by the
    potential.

    Args:
        pot (quippy.Potential): IP to calculate properties with.
        atoms (quippy.AtomsList): list of atoms to calculate correlations with
          respect to.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    dft, ip = _get_xy(pot, atoms, "dft_energy", energy=True)
    vol = np.array([a.get_volume()/a.n for a in atoms])
    
    subplot_kw = {
        "xlabel": "Volume/Atom (A^3)",
        "ylabel": "IP Energy (eV)"
    }

    return PDI(np.array(dft), np.array(vol), subplot_kw=subplot_kw, index=index,
               base64=base64, name="EvsV", imgtype="evsv", folder=folder)
