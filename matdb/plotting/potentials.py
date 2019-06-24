#Copyright (C) 2019  HALL LABS
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#If you have any questions contact: wmorgan@tracy.com
"""Functions for plotting interatomic potential performance vs. correct (REF)
answers within a particular :class:`matdb.atoms.AtomsList`.
"""
from tqdm import tqdm
import numpy as np

from matdb.atoms import Atoms
from matdb.plotting.matd3 import PointDetailImage as PDI

def generate(plots, pot, atoms, folder=None, base64=False, index=0, ndimer=50,
             ntrimer=75, valkey="ref"):
    """Generates a set of plots for a given potential and atoms list. The
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

    .. note:: The specified order of the `plots` will be the same order they will show up in the interactive plot.

    Args:
        plots (str): string of character codes as described above.
        fit (matdb.calculators.basic.AsyncCalculator): trainer that has an IP to calculate
          properties with.
        atoms (matdb.atoms.AtomsList): list of atoms from which to calculate correlations.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        seed_index (int): integer index to use for the first plot. It is
          incremented automatically for each plot added.
        ndimer (int): number of samples to take along the dimer trajectory.
        ntrimer (int): number of samples to take along the trimer angle
          trajectory.
        valkey (str): prefix on parameter and property names that should be used
          for comparison.
    """
    #We need to determine the unique list of elements in the atoms list. Lets
    elements = set()
    for a in atoms:
        elements |= set(a.get_chemical_symbols())

    plotmap = {
        'e': lambda i: energy(pot, atoms, folder, base64, i, valkey=valkey),
        'f': lambda i: force(pot, atoms, folder, base64, i, valkey=valkey),
        'v': lambda i: virial(pot, atoms, folder, base64, i, valkey=valkey),
        'o': lambda i: EvsV(pot, atoms, folder, base64, i, valkey=valkey),
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
        elements (list): list of chemical symbols in the system.
        pot (matdb.calculators.basic.AsyncCalculator): IP to calculate properties with.
        atoms (matdb.atoms.AtomsList): list of atoms from which to calculate correlations.
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

        trimer = Atoms(positions=np.zeros((3, 3)), cell=np.eye(3)*100)
        trimer.positions = 0.
        trimer.set_chemical_symbols(elements)
        trimer.set_calculator(pot)
        #Set the position of the second atom to be at equilibrium with respect
        #to the Vegard's law calculation.
        trimer.positions[1,2] = rvegard

        #Now, vary the angle of the third atom with respect to the original two
        #and see how the angle changes.
        theta = np.linspace(np.pi/6, np.pi, nsamples)
        energy = []
        for t in theta:
            x = rvegard*np.cos(t)
            y = rvegard*np.sin(t)
            trimer.positions[1,2] = x
            trimer.positions[1,2] = y
            energy.append(trimer.get_potential_energy())

        elemstr = trimer.get_chemical_formula()
        img = PDI(theta, np.array(energy), "plot", subplot_kw=subplot_kw,
                  index=index, base64=base64, name=elemstr, imgtype=elemstr,
                  folder=folder)
        result[elemstr] = img

    return result

def _dimer_range(elements):
    """Calculates a reasonable range over which to vary distance in a dimer.

    Args:
        elements (list): list of chemical symbols in the system.

    Returns:
        tuple: array of `(rmin, rvegard, rmax)` where `rmin` and `rmax` are the minimum
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

    .. note:: This produces a **dict** of all possible 2-body interaction between
      distinct element types.

    Args:
        elements (list): list of chemical symbols in the system.
        pot (matdb.calculators.basic.AsyncCalculator): IP to calculate properties with.
        atoms (matdb.atoms.AtomsList): list of atoms from which to calculate correlations.
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

        dimer = Atoms(positions=np.zeros((2, 3)), cell=np.eye(3)*100)
        dimer.positions = 0.
        dimer.set_calculator(pot)
        dimer.set_chemical_symbols(elements)

        if zoom:
            rs = np.linspace(0.7*rmin, pot.cutoff(), nsamples)
        else:
            rs = np.linspace(0.4*rmin, pot.cutoff(), nsamples)

        energy = []
        for r in rs:
            dimer.positions[1,2] = r
            energy.append(dimer.get_potential_energy())

        elemstr = ''.join(elements) + ('z' if zoom else "")
        img = PDI(rs, np.array(energy), "plot", subplot_kw=subplot_kw,
                  index=index, base64=base64, name=elemstr, imgtype=elemstr,
                  folder=folder)
        result[elemstr] = img

    return result

def _get_xy(pot, atoms, prop, peratom=False, energy=False, force=False,
            virial=False):
    """Gets the x and y values for the specified potential and property.

    Args:
        pot (matdb.calculators.basic.AsyncCalculator): IP to calculate properties with.
        atoms (matdb.atoms.AtomsList): list of atoms to calculate correlations with
          respect to.
        prop (str): name of the property on each atoms object.
        peratom (bool): when True, plot per atom quantities.
        energy (bool): when True, calculate the energy.
        force (bool): when True, calculate the forces on each atom.
        virial (bool): when True, calculate the virial tensor..

    Returns:
        tuple: tuple of `(ref, pot)` values.
    """
    ipprop = None
    if energy:
        ipprop = lambda a: a.get_potential_energy()
    if force:
        ipprop = lambda a: a.get_forces()
    if virial:
        ipprop = lambda a: a.get_stress(False)*a.get_volume()

    ref, ip = [], []
    import pudb
    for i, a in tqdm(enumerate(atoms)):
        if hasattr(a, prop):
            a.set_calculator(pot)
            if peratom:
                ref.append(getattr(a, prop)/float(a.n))
                ip.append(ipprop(a)/float(a.n))
            else:
                ref.append(getattr(a, prop))
                ip.append(ipprop(a))
    return (ref, ip)

def energy(pot, atoms, folder=None, base64=False, index=None, valkey="ref"):
    """Produces an energy correlation plot for the specified potential.

    Args:
        pot (matdb.calculators.basic.AsyncCalculator): IP to calculate properties with.
        atoms (matdb.atoms.AtomsList): list of atoms from which to calculate correlations.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    ref, ip = _get_xy(pot, atoms, "{}_energy".format(valkey), energy=True,
                      peratom=True)
    ref, ip = np.array(ref), np.array(ip)

    #Setup the keyword arguments for the axes labels, etc.
    subplot_kw = {
        "xlabel": "REF Energy (eV)",
        "ylabel": "IP Energy (eV)"
    }
    if hasattr(atoms[0],"config_type"):
        ctypes = np.array(atoms.config_type)
    else:
        ctypes = None
    title = "RMSE {0:.4f}".format(np.std(ref-ip))

    return PDI(ref, ip, subplot_kw=subplot_kw, index=index, base64=base64,
               name="Energy", imgtype="energy", folder=folder, withdiag=True,
               title=title, partition=ctypes)

def force(pot, atoms, folder=None, base64=False, index=None, valkey="ref"):
    """Produces a force correlation plot for the specified potential.

    Args:
        pot (matdb.calculators.basic.AsyncCalculator): IP to calculate properties with.
        atoms (matdb.atoms.AtomsList): list of atoms from which to calculate correlations.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    ref, ip = _get_xy(pot, atoms, "{}_force".format(valkey), force=True)
    ref = np.concatenate([np.array(f).flatten() for f in ref]).flatten()
    ip = np.concatenate([np.array(f).flatten() for f in ip]).flatten()

    if hasattr(atoms[0],"config_type"):
        aconfigs = atoms.config_type
        N = np.array(atoms.n)*3
        ctypes = np.array([ac for ac, n in zip(aconfigs, N) for i in range(n)])
        assert len(ctypes) == len(ref)
    else:
        ctypes = None

    #Setup the keyword arguments for the axes labels, etc.
    subplot_kw = {
        "xlabel": "REF Forces (eV/A)",
        "ylabel": "IP Forces (eV/A)"
    }
    title = "RMSE {0:.4f}".format(np.std(ref-ip))
    return PDI(ref.flatten(), ip.flatten(), subplot_kw=subplot_kw, index=index,
               base64=base64, name="Force", imgtype="force",
               folder=folder, withdiag=True, title=title, partition=ctypes)

def virial(pot, atoms, folder=None, base64=False, index=None, valkey="ref"):
    """Produces a virial correlation plot for the specified potential.

    Args:
        pot (matdb.calculators.basic.AsyncCalculator): IP to calculate properties with.
        atoms (matdb.atoms.AtomsList): list of atoms from which to calculate correlations.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    ref, ip = _get_xy(pot, atoms, "{}_virial".format(valkey), virial=True, peratom=True)
    ref = np.hstack([np.array(v).flatten() for v in ref]).flatten()
    ip = np.hstack([np.array(v).flatten() for v in ip]).flatten()

    if hasattr(atoms[0],"config_type"):
        aconfigs = atoms.config_type
        ctypes = np.array([ac for ac in aconfigs for i in range(9)])
        assert len(ctypes) == len(ref)
    else:
        ctypes = None

    #Setup the keyword arguments for the axes labels, etc.
    subplot_kw = {
        "xlabel": "REF Virial (eV/A^3)",
        "ylabel": "IP Virial (eV/A^3)"
    }
    title = "RMSE {0:.4f}".format(np.std(ref-ip))
    return PDI(ref.flatten(), ip.flatten(), subplot_kw=subplot_kw, index=index,
               base64=base64, name="Virial", imgtype="virial",
               folder=folder, withdiag=True, title=title, partition=ctypes)

def EvsV(pot, atoms, folder=None, base64=False, index=None, valkey="ref"):
    """Produces an energy vs. volume plot for the energies predicted by the
    potential.

    Args:
        pot (matdb.calculators.basic.AsyncCalculator): IP to calculate properties with.
        atoms (matdb.atoms.AtomsList): list of atoms from which to calculate correlations.
        folder (str): path to the folder where the saved image should be stored.
        base64 (bool): when True, use a base-64 encoded `src` for the output of
          this image; otherwise, the image is saved to `folder` and only its
          name is returned.
        index (int): integer index of this plot in the parent collection.
    """
    ref, ip = _get_xy(pot, atoms, "{}_energy".format(valkey), energy=True)
    vol = np.array([a.get_volume()/a.n for a in atoms])

    subplot_kw = {
        "xlabel": "Volume/Atom (A^3)",
        "ylabel": "IP Energy (eV)"
    }

    return PDI(np.array(ref), np.array(vol), subplot_kw=subplot_kw, index=index,
               base64=base64, name="EvsV", imgtype="evsv", folder=folder)
