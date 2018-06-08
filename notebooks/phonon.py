"""Group of configurations that is created by displacing a given
unit cell along phonon modes using the eigenvectors.
"""
from matdb import msg
from os import path, remove
import numpy as np
from operator import itemgetter
from .basic import Group

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

        for bandi, f in enumerate(fs):
            rf = np.round(f, 6)
            if rf not in lookup:
                lookup[rf] = (bandi, q)

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

def update_phonons(basic, base):
    """Updates the `basic` phonon settings using the usual defaults. The update
    only happens if the user didn't already specify a value in the config file.
    Args:
        basic (dict): user-specified phonon settings that should be updated to
          include defaults.
        base (DynMatrix): dynamical matrix step to take default values from.
    """
    usuals = {
        "mesh": base.dosmesh,
        "dim": base.supercell
    }
    for k, v in usuals.items():
        if k not in basic:
            basic[k] = v

def modulate_atoms(db):
    """Generates modulated configurations using the dynamical matrix of the
    :class:`DynMatrix` instance.
    Args:
        db (Group): database with parameters needed to module the
            atoms, (Calibration or Modulation databases).
    """
    #Generating the modulation file. We need to sample the DOS in order to
    #compute that correctly.
    dosfile = path.join(db.base.phonodir, "mesh.yaml")
    qvecs = sample_dos(dosfile, sampling=db.sampling, nfreqs=db.nconfigs)
    conffile = path.join(db.base.phonodir, db.confname)

    modstr = [' '.join(map(str, db.base.supercell))]
    for iq, (bandi, qvec) in enumerate(qvecs):
        mstr = "{0:.7f} {1:.7f} {2:.7f} {3:d} {4:.7f} {5:.7f}"
        if hasattr(db, "_amplitudes"):
            A = db._amplitudes[iq]
        else:
            A = np.random.normal(1, 0.25)*db.amplitude
            
        phi = np.random.uniform(0, 180)
        args = qvec + [bandi, A, phi]
        modstr.append(mstr.format(*args))

    phondict = db.phonons.copy()
    phondict["dim"] = db.base.supercell
    phondict["atom_name"] = ' '.join(db.base.parent.species)
    phondict["modulation"] = ', '.join(modstr)
    with open(conffile, 'w') as f:
        for k, v in phondict.items():
            if isinstance(v, (list, tuple, set)):
                value = ' '.join(map(str, v))
                f.write("{} = {}\n".format(k.upper(), value))
            else:
                f.write("{} = {}\n".format(k.upper(), v))

    from matdb.utility import execute
    sargs = ["phonopy", db.confname]
    xres = execute(sargs, db.base.phonodir, venv=True)

    #Make sure that phonopy actually produced files; otherwise show the output
    #(phonopy doesn't write to stderr, only stdout).
    testmod = path.join(db.base.phonodir, "MPOSCAR-001")
    if not path.isfile(testmod):#pragma: no cover
        msg.err(''.join(xres["output"]))
            
class Calibration(Group):
    """Represents a set of modulated sub-configurations of differing amplitude,
    used to determine the maximum modulation amplitude where the force is still
    in the linear regime .
    Args:
        atoms (matdb.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.Sequence): parent sequence
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
      :class:`Group`.
    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
        phonons (dict): specifying additional settings for `phonopy`
          configuration files (i.e., differing from, or in addition to those in
          the global set).
        base (DynMatrix): reference database that computed the dynamical matrix
          for the seed configuration.
        confname (str): name of the phonopy configuration file used for the
          modulations in this database.
        amplitudes (dict): of logarithmically-spaced amplitudes to
          modulate the seed configuration with. Keys are the `cid` keys in
          :attr:`configs`.
        outfile (str): path to the `calibration.dat` file that contains the list
          of amplitudes and displacements from the calibration run.
        sampling (str): specifies how the DOS is sampled for this database to
          produce modulated sub-configurations.
    """
    confname = "calibrate.conf"
    sampling = "uniform"
    
    def __init__(self, atoms=None, root=None, parent=None,
                 kpoints={}, incar={}, phonons={}, execution={},
                 nconfigs=10, name="calib", dynmat="dynmatrix",
                 config_type="calib"):
        self.name = name
        super(Calibration, self).__init__(atoms, incar, kpoints, execution,
                                                path.join(root, self.name),
                                                parent, "C", nconfigs, config_type)
        self.base = self.parent.steps[dynmat]
        self.phonons = phonons
        update_phonons(self.phonons, self.base)

        #Calculate which amplitudes to use for the calibration based on the
        #number of desired configurations (calibration points).
        self._amplitudes = np.logspace(0, 1.7, nconfigs)
        self.outfile = path.join(self.root, "calibration.dat")

        self.amplitudes = {}
        if len(self.configs) == self.nconfigs:
            for cid in self.configs:
                self.amplitudes[cid] = self._amplitudes[cid-1]
        
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
        if not super(Calibration, self).cleanup():
            msg.warn("cannot cleanup calibration; not all configs ready.")
            return False

        success = self.xyz()
        if not success:
            msg.warn("could not extract the calibration XYZ configurations.")
            return False
        else:
            imsg = "Extracted calibration configs from {0:d} folders."
            msg.okay(imsg.format(len(self.configs)))

        #Read in the XYZ file and extract the forces on each atom in each
        #configuration.
        from matdb.atoms import AtomsList
        forces = {}
        failed = 0
        for cid, folder in self.configs.items():
            #Find the mean, *absolute* force in each of the directions. There
            #will only be one atom in the atoms list. If the calculation didn't
            #finish, then we exclude it. This happens for some of the
            #calibration runs if the atoms are too close together.
            try:
                al = AtomsList(path.join(folder, "output.xyz"))
                forces[cid] = np.mean(np.abs(np.array(al[0].dft_force)), axis=1)
            except:
                failed += 1
                pass

        if failed > 0:
            msg.warn("couldn't extract forces for {0:d} configs.".format(failed))

        if len(forces) > 0:
            fmt = "{0:.7f}  {1:.7f}  {2:.7f}  {3:.7f}\n"
            with open(self.outfile, 'w') as f:
                for cid in forces:
                    A, F = self.amplitudes[cid], forces[cid]
                    f.write(fmt.format(A, *F))
        else:
            msg.warn("no forces available to write {}.".format(self.outfile))

        return len(forces) > 3
    
    def setup(self, rerun=False):
        """Displaces the seed configuration with varying amplitudes so that the
        resulting forces can be calibrated sensibly.
        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """
        if super(Calibration, self).setup():
            return

        #We can't module atoms unless the phonon base is ready.
        if not self.base.ready():
            return

        #Don't reproduce the folders if they already exist or we have already
        #computed.
        if self.ready():
            return
        
        modulate_atoms(self)
        
        from os import getcwd, chdir, remove
        from glob import glob
        from matdb.atoms import Atoms
        current = getcwd()
        chdir(self.base.phonodir)

        try:
            for mi, mposcar in enumerate(sorted(list(glob("MPOSCAR-*")))):
                if mposcar == "MPOSCAR-orig":
                    remove(mposcar)
                    continue
                
                cid = int(mposcar.split('-')[1])
                matoms = Atoms(mposcar, format="vasp")
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
            raise NotImplementedError("Automatic calibration not configured.")
            
class Modulation(Group):
    """Represents a set of displaced configurations where atoms are
    moved, within a supercell, according to phonon eigenmodes.
    Args:
        atoms (matdb.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.Sequence): parent collection
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
      :class:`Group`.
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
        calibrator (Calibration): instance used to calculate force
          vs. amplitude for the seed configuration.
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
        confname (str): name of the phonopy configuration file used for the
          modulations in this database.
    """
    confname = "modulate.conf"

    def __init__(self, atoms=None, root=None, parent=None,
                 phonons={}, nconfigs=100, calibrate=True,
                 amplitude=None, sampling="uniform",
                 name="modulations", dynmat="dynmatrix",
                 calibrator="calib", config_type="ph"):
        self.name = name
        super(Modulation, self).__init__(atoms, incar, kpoints, execution,
                                             path.join(root, self.name),
                                             parent, "M", nconfigs, config_type)
        self.sampling = sampling
        self.calibrate = calibrate
        
        #Setup a calibrator if automatic calibration was selected.
        if calibrate and amplitude is None:
            self.calibrator = Calibration(atoms, root, parent, kpoints, incar,
                                          phonons, execution, calibrate,
                                          name=calibrator, dynmat=dynmat)
            self.parent.steps[calibrator] = self.calibrator
            self.amplitude = self.calibrator.infer_amplitude()
            calibrated = "*calibrated* "
        else:
            self.calibrator = None
            self.amplitude = amplitude
            calibrated = ""

        imsg = "Using {} as modulation {}amplitude for {}."
        msg.info(imsg.format(self.amplitude, calibrated, self.parent.name), 2)

        self.base = self.parent.steps[dynmat]
        self.phonons = phonons
        update_phonons(self.phonons, self.base)

    def ready(self):
        """Determines if this database is finished calculating by testing the
        existence of the xyz database file in the root folder.
        """
        target = path.join(self.root, "output.xyz")
        result = False
        if path.isfile(target):
            from matdb.utility import linecount
            #We add +2 for the parameter line and the number of atoms.
            #This doesn't work as advertised (it's off by a factor). debug after refactor.
            lpconfig = self.base.atoms.n * np.linalg.det(np.array(self.base.supercell).reshape(3,3)) + 2
            nlines = linecount(target)
            nconfigs = nlines/lpconfig
            result = nconfigs  == len(self.configs)
            if not result:
                wmsg = ("Number of structures in `output.xyz` does not match "
                        "number of requested configs. Found {0} configs in"
                        " {1} lines.")
                msg.warn(wmsg.format(nconfigs, nlines))
        return result

    def cleanup(self):
        """Generates the XYZ database file for all the sub-configs in this
        phonon database.
        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.
        """
        if not super(Modulation, self).cleanup():
            return
        
        return self.xyz(combine=True)
    
    def setup(self, rerun=False):
        """Displaces the seed configuration preparatory to calculating the force
        sets for phonon spectra.
        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """
        folders_ok = super(Modulation, self).setup()
        if folders_ok and not rerun:
            return

        #We can't modulate atoms unless the phonon base is ready.
        if not self.base.ready():
            return

        #We can't module in calibrate mode unless the calibrator is also ready.
        if self.amplitude is None:
            return

        if not folders_ok:
            modulate_atoms(self)

            from os import getcwd, chdir, remove
            from glob import glob
            from matdb.atoms import Atoms
            current = getcwd()
            chdir(self.base.phonodir)

            from tqdm import tqdm
            try:
                for mi, mposcar in tqdm(enumerate(sorted(list(glob("MPOSCAR-*"))))):
                    if mposcar == "MPOSCAR-orig":
                        continue
                    
                    cid = int(mposcar.split('-')[1])
                    matoms = Atoms(mposcar, format="vasp")
                    self.create(matoms, cid)
                    #Remove the MPOSCAR file so that the directory isn't polluted.
                    remove(mposcar)
            finally:
                chdir(current)
            
        # Last of all, create the job file to execute the job array.
        self.jobfile()
