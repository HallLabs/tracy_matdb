"""Database of configurations that is created by displacing a given
unit cell along phonon modes using the eigenvectors.
"""
from .basic import Database
from matdb import msg
from os import path
import numpy as np

class PhononBase(Database):
    """Sets up the displacement calculations needed to construct the dynamical
    matrix. The dynamical matrix is required by :class:`PhononDatabase` to
    create the individual modulations.

    Args:
        atoms (quippy.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        incar (dict): specify additional settings for the INCAR file (i.e.,
          differing from, or in addition to those in the global set).
        kpoints (dict): specify additional settings for the PRECALC file (i.e.,
          differing from, or in addition to those in the global set).
        phonons (dict): specifying additional settings for `phonopy`
          configuration files (i.e., differing from, or in addition to those in
          the global set).

    .. note:: Additional attributes are also exposed by the super class
      :class:`Database`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
        supercell (list): of `int`; number of cells in each direction for
          generating the supercell.
        phonodir (str): directory in which all `phonopy` executions take place.
        grid (list): of `int`; number of splits in each reciprocal
          lattice vector according to Monkhorst-Pack scheme.
    """
    name = "phonbase"
    
    def __init__(self, atoms=None, root=None, parent=None,
                 kpoints={}, incar={}, phonons={}, execution={}):
        super(PhononBase, self).__init__(atoms, incar, kpoints, execution,
                                         path.join(root, self.name),
                                         parent, "W")
        self.supercell = list(phonons.get("dim", [2, 2, 2]))
        self.grid = list(phonons.get("mp", [20, 20, 20]))
        self.phonodir = path.join(self.root, "phonopy")
        
        from os import mkdir
        if not path.isdir(self.phonodir):
            mkdir(self.phonodir)

        self._update_incar()

    def _update_incar(self):
        """Adds the usual settings for the INCAR file when performing
        frozen-phonon calculations. They are only added if they weren't already
        specified in the config file.
        """
        usuals = {
            "encut": 500,
            "ibrion": -1,
            "ediff": '1.0e-08',
            "ialgo": 38,
            "ismear": 0,
            "lreal": False
        }
        for k, v in usuals.items():
            if k not in self.incar:
                self.incar[k] = v
            
    def ready(self):
        """Returns True if all the phonon calculations have been completed, the
        force sets have been created, and the DOS has been calculated.
        """
        #If the DOS has been calculated, then all the other steps must have
        #completed correctly.
        dosfile = path.join(self.phonodir, "mesh.yaml")
        return path.isfile(dosfile)

    def calc_DOS(self, recalc=False):
        """Calculates the *total* density of states.

        Args:
            recalc (bool): when True, recalculate the DOS, even if the
              file already exists.
        """
        dosfile = path.join(self.phonodir, "mesh.yaml")
        if not recalc and path.isfile(dosfile):
            return

        #Make sure we have calculated the force sets already.
        self.calc_forcesets(recalc)
        settings = {
            "ATOM_NAME": ' '.join(self.parent.species),
            "DIM": ' '.join(map(str, self.supercell)),
            "MP": ' '.join(map(str, self.grid))
        }
        with open(path.join(self.phonodir, "dos.conf"), 'w') as f:
            for k, v in settings.items():
                f.write("{} = {}\n".format(k, v))

        from matdb.utility import execute
        sargs = ["phonopy", "-p", "dos.conf"]
        xres = execute(sargs, self.phonodir, venv=True)
        if not path.isfile(dosfile):
            msg.err("could not calculate the DOS; see errors.")

        if len(xres["error"]) > 0:
            msg.std(''.join(xres["error"]))
            
    def calc_forcesets(self, recalc=False):
        """Extracts the force sets from the displacement calculations.

        Args:
            recalc (bool): when True, recalculate the force sets, even if the
              file already exists.
        """
        fsets = path.join(self.phonodir, "FORCE_SETS")
        if not recalc and path.isfile(fsets):
            return
        
        #First, make sure we have `vasprun.xml` files in each of the
        #directories.
        vaspruns = []
        for i, folder in self.configs.items():
            vasprun = path.join(folder, "vasprun.xml")
            if not path.isfile(vasprun):
                msg.err("vasprun.xml does not exist for {}.".format(folder))
            else:
                vaspruns.append(vasprun)

        if len(vaspruns) == len(self.configs):
            from matdb.utility import execute
            sargs = ["phonopy", "-f"] + vaspruns
            xres = execute(sargs, self.phonodir, venv=True)

        if not path.isfile(fsets):
            msg.err("Couldn't create the FORCE_SETS:")
        if len(xres["error"]) > 0:
            msg.std(''.join(xres["error"]))
            
    def setup(self):
        """Displaces the seed configuration preparatory to calculating the force
        sets for phonon spectra.
        """
        if super(PhononBase, self).setup():
            return
        
        from ase.io import write
        from matdb.utility import execute
        write(path.join(self.phonodir, "POSCAR"), self.atoms, "vasp")        
        scell = ' '.join(map(str, self.supercell))
        sargs = ["phonopy", "-d", '--dim="{}"'.format(scell)]
        pres = execute(sargs, self.phonodir, venv=True)

        from os import getcwd, chdir, mkdir, rename
        from glob import glob
        current = getcwd()
        chdir(self.phonodir)

        try:
            from quippy.atoms import Atoms
            for dposcar in glob("POSCAR-*"):
                dind = int(dposcar.split('-')[1])
                datoms = Atoms(dposcar, format="POSCAR")
                self.create(datoms)
        finally:
            chdir(current)

        # Last of all, create the job file to execute the job array.
        self.jobfile()

    def cleanup(self, recalc=False):
        """Runs post-DFT execution routines to calculate the force-sets and the
        density of states.

        Args:
            recalc (bool): when True, redo any calculations that use the DFT
              outputs to find other quantities.

        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.
        """
        self.calc_forcesets(recalc)
        self.calc_DOS(recalc)
        return self.ready()

def sample_dos(meshfile, sampling="uniform", nfreqs=100):
    """Samples the DOS to extract frequencies from which modulations can be
    generated.
    
    Args:
        meshfile (str): path to the `mesh.yaml` file that the frequencies and
          q-vectors can be extracted from.
        sampling (str): one of ['uniform', 'sample', 'top'], where the method
          dictates how frequencies are selected from the DOS for the seed
          configuration's phonon spectrum.
        nfreqs (int): number of frequencies to return by sampling the DOS.

    - *uniform*: frequencies are selected uniformly from the list of *unique*
       frequencies in the DOS.
    - *sample*: frequencies are chosen randomly *and* weighted by the q-point
       weight in the BZ and the number of times the frequency shows up.
    - *top*: the top N *unique* frequencies are selected.

    .. note:: Because each atomic degree of freedom produces 3 phonon bands, a
      cell with 4 unique atoms will produce 12 different frequencies. When the
      DOS is sampled, we consider each of these frequencies as independent in
      the overall BZ. Thus, if a 20x20x20 q-point grid is used for sampling in
      the `mesh.yaml` file, then we would have 8,000 * 12 = 96,000 frequencies
      to choose from, each with a corresponding q-vector.

    Returns:
        numpy.ndarray: where each row is a q-vector in reciprocal space
        corresponding to a frequency that was selected using the specified
        sampling method.
    """
    import yaml
    with open(meshfile, 'r') as stream:
        dmesh = yaml.load(stream)
        
    freqs = []
    lookup = {}
    for ph in dmesh["phonon"]:
        #Depending on the type of sampling we do, we are interested in
        #either the mere existing of a frequency or how often it actually
        #appears in the BZ.
        fs = [b["frequency"] for b in ph["band"]]
        q = ph["q-position"]
        w = ph["weight"]
        
        if sampling in ["uniform", "top"]:
            freqs.extend(fs)
        elif sampling == "sample":
            freqs.extend(fs*w)

        for f in fs:
            rf = np.round(f, 6)
            if rf not in lookup:
                lookup[rf] = q

    #Now we can do the actual sampling.
    if sampling == "uniform":
        dfreqs = np.unique(freqs)
        sample = np.random.choice(dfreqs, size=nfreqs)
    elif sampling == "top":
        sfreqs = np.sort(np.unique(freqs))
        sample = sfreqs[-nfreqs:]
    elif sampling == "sample":
        sample = np.random.choice(freqs, size=nfreqs)

    return [lookup[f] for f in np.round(sample, 6)]

def update_phonons(basic):
    """Updates the `basic` phonon settings using the usual defaults. The update
    only happens if the user didn't already specify a value in the config file.

    Args:
        basic (dict): user-specified phonon settings that should be updated to
          include defaults.
    """
    usuals = {
        "mesh": [13, 13, 13],
        "dim": [2, 2, 2]
    }
    for k, v in usuals.items():
        if k not in basic:
            basic[k] = v

def modulate_atoms(db):
    """Generates modulated configurations using the dynamical matrix of the
    :class:`PhononBase` instance.

    Args:
        db (Database): database with parameters needed to module the atoms.
    """
    #Generating the modulation file. We need to sample the DOS in order to
    #compute that correctly.
    dosfile = path.join(db.base.phonodir, "mesh.yaml")
    qvecs = sample_dos(dosfile, sampling=db.sampling, nfreqs=db.nconfigs)
    conffile = path.join(db.base.phonodir, db.confname)

    modstr = [' '.join(map(str, db.phonons["dim"]))]
    for iq, qvec in enumerate(qvecs):
        mstr = "{0:.7f} {1:.7f} {2:.7f} {3:d} {4:.7f} {5:.7f}"
        if hasattr(db, "_amplitude"):
            A = np.random.normal(1, 0.25)*db._amplitude[iq]
        else:
            A = np.random.normal(1, 0.25)*db.amplitude
        phi = np.random.uniform(0, 180)
        args = qvec + [bandi, A, phi]
        modstr.append(mstr.format(*args))

    phondict = db.phonons.copy()
    phondict["atom_name"] = ' '.join(db.base.parent.species)
    phondict["modulation"] = ', '.join(modstr)
    with open(conffile, 'w') as f:
        for k, v in phondict.items():
            f.write("{} = {}\n".format(k.upper(), v))

    from matdb.utility import execute
    sargs = ["phonopy", db.confname]
    xres = execute(sargs, db.base.phonodir, venv=True)
            
class PhononCalibration(Database):
    """Represents a set of modulated sub-configurations of differing amplitude,
    used to determine the maximum modulation amplitude where the force is still
    in the linear regime.

    Args:
        atoms (quippy.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        incar (dict): specify additional settings for the INCAR file (i.e.,
          differing from, or in addition to those in the global set).
        kpoints (dict): specify additional settings for the PRECALC file (i.e.,
          differing from, or in addition to those in the global set).
        phonons (dict): specifying additional settings for `phonopy`
          configuration files (i.e., differing from, or in addition to those in
          the global set).
        nconfigs (int): the number of different *amplitudes* to try out in
          calibrating.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Database`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
        phonons (dict): specifying additional settings for `phonopy`
          configuration files (i.e., differing from, or in addition to those in
          the global set).
        base (PhononBase): reference database that computed the dynamical matrix
          for the seed configuration.
        confname (str): name of the phonopy configuration file used for the
          modulations in this database.
        amplitudes (dict): of logarithmically-spaced amplitudes to
          modulate the seed configuration with. Keys are the `cid` keys in
          :attr:`configs`.
        outfile (str): path to the `calibration.dat` file that contains the list
          of amplitudes and displacements from the calibration run.
    """
    name = "phoncalib"
    confname = "calibrate.conf"

    def __init__(self, atoms=None, root=None, parent=None,
                 kpoints={}, incar={}, phonons={}, execution={}, nconfigs=10):
        super(PhononCalibration, self).__init__(atoms, incar, kpoints, execution,
                                                path.join(root, self.name),
                                                parent, "C", nconfigs)
        self.base = self.parent.databases[PhononBase.name]
        self.phonons = phonons
        
        update_phonons(self.phonons)

        #Calculate which amplitudes to use for the calibration based on the
        #number of desired configurations (calibration points).
        self._amplitudes = np.logspace(1, 500, nconfigs)
        self.amplitudes = {}
        self.outfile = path.join(self.root, "calibration.dat")
        
    def ready(self):
        """Determines if this database is finished calculating by testing the
        existence of the xyz database file in the root folder.
        """
        return path.isfile(self.outfile)

    def cleanup(self):
        """Extracts the calibration information from the configurations to
        determine the maiximum allowable amplitude to maintain linear force
        regime.

        Returns:
           bool: True if the amplitude calibration is ready.
        """
        success = self.xyz()
        if not success:
            return False

        #Read in the XYZ file and extract the forces on each atom in each
        #configuration.
        import quippy
        forces = {}
        success = None
        for cid, folder in self.configs.items():
            #Find the mean, *absolute* force in each of the directions. There
            #will only be one atom in the atoms list.
            try:
                al = quippy.AtomsList(path.join(folder, "output.xyz"))
                forces[cid] = np.mean(np.abs(np.array(al[0].force)), axis=1)
            except:
                success = False
                break
        else:
            success = True

        #The success local will only be True if we were able to extract a force
        #matrix for each of the configurations generated by this database.
        if success:
            fmt = "{0:.7f}  {1:.7f}  {2:.7f}  {3:.7f}"
            with open(self.outfile, 'w') as f:
                for cid in self.configs:
                    A, F = self.amplitudes[cid], forces[cid]
                    f.write(fmt.format(A, *F))
    
    def setup(self):
        """Displaces the seed configuration with varying amplitudes so that the
        resulting forces can be calibrated sensibly.
        """
        if super(PhononCalibration, self).setup():
            return

        #We can't module atoms unless the phonon base is ready.
        if not self.base.ready():
            return
        
        modulate_atoms(self)
        
        from os import getcwd, chdir, remove
        from glob import glob
        from quippy.atoms import Atoms
        current = getcwd()
        chdir(self.base.phonodir)

        try:
            for mi, mposcar in enumerate(sorted(list(glob("MPOSCAR-*")))):
                cid = int(mposcar.split('-')[1])
                matoms = Atoms(mposcar, format="POSCAR")
                self.create(matoms, cid)
                self.amplitudes[cid] = self._amplitudes[mi]
                #Remove the MPOSCAR file so that the directory isn't polluted.
                remove(mposcar)
        finally:
            chdir(current)
            
        # Last of all, create the job file to execute the job array.
        self.jobfile()

    def infer_amplitude(self):
        """Tries to infer the maximum amplitude that can be used for modulation
        such that the average forces experienced by all atoms in the seed
        configuration's supercell are still in the linear regime.
        """
        #Our approach is to interpolate linearly starting with the two closest
        #points and then moving away one point at a time until the error starts
        #to diverge between the linear interpolation and the actual points.
        if not self.ready():
            return None
        else:
            return 1.
            
class PhononDatabase(Database):
    """Represents a set of displaced configurations where atoms are
    moved, within a supercell, according to phonon eigenmodes.

    Args:
        atoms (quippy.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        incar (dict): specify additional settings for the INCAR file (i.e.,
          differing from, or in addition to those in the global set).
        kpoints (dict): specify additional settings for the PRECALC file (i.e.,
          differing from, or in addition to those in the global set).
        phonons (dict): specifying additional settings for `phonopy`
          configuration files (i.e., differing from, or in addition to those in
          the global set).
        calibrate (bool): when True, the maximum amplitude possible will be
          selected ensuring that the force is still in the linear regime.
        amplitude (float): amplitude of displacement :math:`A` for eigenmode
          modulation. 
        sampling (str): on of ['uniform', 'sample', 'top'], where the method
          dictates how frequencies are selected from the DOS for the seed
          configuration's phonon spectrum.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Database`.

    Attributes:
        sampling (str): on of ['uniform', 'sample', 'top'], where the method
          dictates how frequencies are selected from the DOS for the seed
          configuration's phonon spectrum.
        supercell (list): of `int`; number of cells in each direction for
          generating the supercell.
        phonodir (str): directory in which all `phonopy` executions take place.
        calibrate (bool): when True, the maximum amplitude possible will be
          selected ensuring that the force is still in the linear regime.
        amplitude (float): amplitude of displacement :math:`A` for eigenmode
          modulation. 
        calibrator (PhononCalibration): instance used to calculate force
          vs. amplitude for the seed configuration.
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
        confname (str): name of the phonopy configuration file used for the
          modulations in this database.
    """
    name = "phonons"
    confname = "modulate.conf"

    def __init__(self, atoms=None, root=None, parent=None,
                 kpoints={}, incar={}, phonons={}, execution={}, nconfigs=100,
                 calibrate=True, amplitude=1., sampling="uniform"):
        super(PhononDatabase, self).__init__(atoms, incar, kpoints, execution,
                                             path.join(root, self.name),
                                             parent, "M", nconfigs)
        self.sampling = sampling
        self.calibrate = calibrate

        #Setup a calibrator if automatic calibration was selected.
        if calibrate:
            self.calibrator = PhononCalibration(atoms, root, parent, kpoints,
                                                incar, phonons, execution, calibrate)
            self.parent.databases[PhononCalibration.name] = self.calibrator
            self.amplitude = self.calibrator.infer_amplitude()
        else:
            self.calibrator = None
            self.amplitude = amplitude
            
        self.base = self.parent.databases[PhononBase.name]
        self.phonons = phonons        
        update_phonons(self.phonons)

    def ready(self):
        """Determines if this database is finished calculating by testing the
        existence of the xyz database file in the root folder.
        """
        return (path.isfile(path.join(self.root, "output.xyz")) and
                len(self.configs) == self._nsuccess)

    def cleanup(self):
        """Generates the XYZ database file for all the sub-configs in this
        phonon database.

        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.
        """
        return self.xyz()
    
    def setup(self):
        """Displaces the seed configuration preparatory to calculating the force
        sets for phonon spectra.
        """
        if super(PhononDatabase, self).setup():
            return

        #We can't module atoms unless the phonon base is ready.
        if not self.base.ready():
            return

        #We can't module in calibrate mode unless the calibrator is also ready.
        if self.amplitude is None:
            return
        
        modulate_atoms(self)

        from os import getcwd, chdir, remove
        from glob import glob
        from quippy.atoms import Atoms
        current = getcwd()
        chdir(self.base.phonodir)

        try:
            for mi, mposcar in enumerate(sorted(list(glob("MPOSCAR-*")))):
                cid = int(mposcar.split('-')[1])
                matoms = Atoms(mposcar, format="POSCAR")
                self.create(matoms, cid)
                #Remove the MPOSCAR file so that the directory isn't polluted.
                remove(mposcar)
        finally:
            chdir(current)
            
        # Last of all, create the job file to execute the job array.
        self.jobfile()
