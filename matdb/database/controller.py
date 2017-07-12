"""Exposes classes and functions for interacting with the database
folders via a simple configuration file.
"""
from os import path
from matdb import msg
import numpy as np

class DatabaseCollection(object):
    """Represents a collection of databases (all inheriting from
    :class:`Database`) that are all related be the atomic species that they
    model.

    .. note:: See the list of attributes below which correspond to the sections
    in the YAML database specification file.

    Args:
        name (str): name of the configuration that this database collection is
          operating for.
        poscar (str): name of the POSCAR file in `root` to extract atomic
          configuration information from.
        root (str): root directory in which all other database collections for
          the configurations in the same specification will be stored.
        parent (Controller): instance controlling multiple configurations.
        database (dict): with keys and values describing the kinds of
          sub-configuration databases to setup.

    Attributes:
        title (str): title for the alloy system that these databases work with.
        species (list): of `str` element names in the alloy system.
        potcars (dict): keys are lower-case element names; values are *suffixes*
          for the various pseudo-potentials that ship with VASP.
        root (str): root directory in which all other database directories will
          be stored. Defaults to the current directory.
        plotdir (str): path to the directory to store plots in for this
          database. Defaults to the parent controller's plot directory.
        incar (dict): keys are valid settings in the INCAR file that should be
          used as defaults for *all* database calculations (with their
          corresponding values).
        kpoints (dict): keys are valid settings for the Mueller k-point PRECALC
          file that should be used as default for *all* database calculations.
        databases (dict): keys are database types (e.g. `liquid`, `phonon`,
          etc.); values are the corresponding class instances.
        order (list): of string database names (attribute `name` of each
          database class) in the order in which they should be processed.
        parent (Controller): instance controlling multiple configurations.
        configdbs (list): of `str` name attributes from databases that can
          actually be used for potential training. Several support databases
          (like PhononBase and PhononCalibration) don't provide configurations
          for training and testing.
    """
    order = ["phonbase", "phoncalib", "phonons"]
    configdbs = ["phonons"]
    
    def __init__(self, name, poscar, root, parent, database):
        from quippy.atoms import Atoms
        self.name = name
        self.atoms = Atoms(path.join(root, poscar), format="POSCAR")
        self.root = path.join(root, name)
        
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        parrefs = ["species", "potcars", "incar", "kpoints", "execution",
                   "plotdir"]
        for ref in parrefs:
            setattr(self, ref, getattr(parent, ref))
        self.parent = parent

        self._trainfile = None
        """str: name of the file that stores the XYZ *training* database.
        """
        self._holdoutfile = None
        """str: name of the file that stores the XYZ *validation* database.
        """        

        #See if we have a training and holdout file specified already.
        from glob import glob
        from matdb.utility import chdir
        with chdir(self.root):
            for xyzname in glob("*.xyz"):
                if "train" in xyzname:
                    self._trainfile = xyzname
                if "hold" in xyzname:
                    self._holdoutfile = xyzname
        
        from importlib import import_module
        self._dbsettings = database
        """dict: with keys and values describing the kinds of sub-configuration
        databases to setup.
        """

        self.databases = {}
        for dbspec in database:
            modname, clsname = dbspec["type"].split('.')
            fqdn = "matdb.database.{}".format(modname)
            module = import_module(fqdn)
            if not hasattr(module, clsname):
                #We haven't implemented this database type yet, just skip the
                #initialization for now.
                continue
            
            cls = getattr(module, clsname)

            #Make a copy of the original dictionary so that we don't mess up the
            #pointers; then add in the keyword arguments that are missing.
            cpspec = dbspec.copy()
            del cpspec["type"]
            cpspec["atoms"] = self.atoms
            cpspec["root"] = self.root
            cpspec["parent"] = self

            #Handle the special cases where settings are specified uniquely for
            #each of the configurations separately.
            for k in list(cpspec.keys()):
                if isinstance(cpspec[k], dict) and self.name in cpspec[k]:
                    cpspec[k] = cpspec[k][self.name]
            
            instance = cls(**cpspec)
            self.databases[instance.name] = instance

    def status(self):
        """Prints a status message for each of the databases relative
        to VASP execution status.
        """
        for dbname in self.order:
            if dbname not in self.databases:
                continue

            db = self.databases[dbname]
            imsg = "{}:{} => {}".format(self.name, dbname, db.status())
            msg.info(imsg)
        msg.blank(level=1)
            
    def execute(self):
        """Submits job array files for any of the databases that are ready to
        execute, but which haven't been submitted yet.
        """
        ready = True
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            if ready:
                ready = (self.databases[dbname].ready() or
                         self.databases[dbname].execute())
            else:
                imsg = ("Database {}:{} is not ready to execute yet, or is "
                        "already executing. Done.")
                msg.info(imsg.format(self.name, dbname))
                break
        msg.blank()

    @property
    def train_file(self):
        """Returns the full path to the XYZ database file that can be
        used for training.
        """
        return path.join(self.root, self._trainfile)

    @property
    def holdout_file(self):
        """Returns the full path to the XYZ database file that can be
        used to validate the potential fit.
        """
        return path.join(self.root, self._holdoutfile)

    def split(self, train_perc, trainfile="train.xyz", holdfile="holdout.xyz",
              recalc=False):
        """Splits the total available data in all databases into a
        training and holdout set.

        Args:
            train_perc (float): percentage of the data to use for training.
            trainfile (str): name of the training XYZ file to create.
            holdfile (str): name of the holdout/validation XYZ file to create.
            recalc (bool): when True, re-split the data and overwrite any
              existing *.xyz files.
        """
        if self._trainfile is None:
            self._trainfile = trainfile
        if self._holdoutfile is None:
            self._holdoutfile = holdfile

        if (path.isfile(self.train_file) and path.isfile(self.holdout_file)
            and not recalc):
            return
        
        #Compile a list of all the sub-configurations we can include in the
        #training.
        Ntot = 0
        subconfs = {}
        for dbname in self.configdbs:
            for fi, f in enumerate(self.databases[dbname].configs.values()):
                subconfs[fi + Ntot] = f
            Ntot += len(subconfs)

        Ntrain = int(np.ceil(Ntot*train_perc))
        ids = np.arange(len(subconfs))
        tids = np.random.choice(ids, size=Ntrain)
        hids = np.setdiff1d(ids, tids)

        from matdb.io import vasp_to_xyz
        from tqdm import tqdm
        trainfiles = []
        for tid in tqdm(tids):
            if vasp_to_xyz(subconfs[tid]):
                trainfiles.append(path.join(subconfs[tid], "output.xyz"))
                
        holdfiles = []
        for hid in tqdm(hids):
            if vasp_to_xyz(subconfs[hid]):
                holdfiles.append(path.join(subconfs[hid], "output.xyz"))

        from matdb.utility import cat
        cat(trainfiles, self.train_file)
        cat(holdfiles, self.holdout_file)
        
    def cleanup(self):
        """Runs the cleanup methods of each database in the collection, in the
        correct order.
        """
        ready = True
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            
            if ready:
                ready = self.databases[dbname].cleanup()
                
            if not ready:
                imsg = "Database {}:{} is not ready yet. Done."
                msg.info(imsg.format(self.name, dbname))
                break
        msg.blank()
            
    def setup(self):
        """Sets up the database collection by generating the POTCAR file and
        initializing any databases that haven't already been initialized.

        .. note:: The db setup functions are setup to only execute once, and then
           only if their dependencies have completed their calculations. This
           method can, therefore, be safely called repeatedly between different
           terminal sessions.
        """
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            db = self.databases[dbname]
            msg.info("Setting up database {}:{}".format(self.name, dbname))
            db.setup()
        msg.blank()
        
class Controller(object):
    """Implements methods for tying a configuration dictionary (in
    YAML format) to instances of various databases.

    Args:
        config (str): path to the configuration YAML file that
          specifies all information for constructing the set of databases.

    Attributes:
        specs (dict): the parsed settings from the YAML configuration file.
        collections (dict): keys are configuration names listed in attribute
          `configs` of the YAML file. Values are the :class:`DatabaseCollection`
          instances.
        plotdir (str): path to the directory to store plots in for all
          databases.
    """
    def __init__(self, config):
        import yaml
        with open(config, 'r') as stream:
            self.specs = yaml.load(stream)

        #We allow the user to specify paths relative the matdb repo; this is
        #mainly to allow unit testing to run smoothly.
        from matdb.utility import chdir, reporoot
        with chdir(reporoot):
            self.root = path.abspath(path.expanduser(self.specs["root"]))
            
        self.plotdir = path.join(self.root, "plots")
        self.title = self.specs["title"]
        self.collections = {}
        self.species = [s for s in self.specs["species"]]
        self.potcars = self.specs["potcars"]
        self.incar = self.specs.get("incar", {})
        self.kpoints = self.specs.get("kpoints", {})
        self.execution = self.specs["execution"]

        for cspec in self.specs["configs"]:
            name, poscar = cspec["name"], cspec["poscar"]
            instance = DatabaseCollection(name, poscar, self.root, self,
                                          self.specs.get("databases", []))
            self.collections[name] = instance

        from os import mkdir
        if not path.isdir(self.plotdir):
            mkdir(self.plotdir)
            
        #Extract the POTCAR that all the databases are going to use. TODO: pure
        #elements can't use this POTCAR, so we have to copy the single POTCAR
        #directly for those databases; update the create().
        self.POTCAR()

        #If the controller is going to train any potentials, we also need to 
        self.trainer = None
        if "training" in self.specs:
            from matdb.fitting.gap import GAPTrainer
            gpdict = self.specs["training"].copy()
            gpdict["root"] = self.root
            gpdict["db"] = self
            self.trainer = GAPTrainer(**gpdict)

        self._trainfile = None
        """str: name of the file that stores the XYZ *training* database.
        """
        self._holdoutfile = None
        """str: name of the file that stores the XYZ *validation* database.
        """
            
    def setup(self):
        """Sets up each of configuration's databases.
        """
        for name, coll in self.collections.items():
            coll.setup()

    def cleanup(self):
        """Runs cleanup on each of the configuration's databases.
        """
        for name, coll in self.collections.items():
            coll.cleanup()

    def execute(self):
        """Submits job array scripts for each database collection.
        """
        for name, coll in self.collections.items():
            coll.execute()

    def status(self):
        """Prints status messages for each of the configuration's databases.
        """
        for name, coll in self.collections.items():
            coll.status()   
            
    def POTCAR(self):
        """Creates the POTCAR file using the pseudopotential and version
        specified in the file.
        """
        target = path.join(self.root, "POTCAR")
        if not path.isfile(target):
            # Make sure that the POTCAR version and pseudopotential type match
            # up so that we don't get nasty surprises.
            potsrc = path.join(self.potcars["directory"],
                               "pot{}".format(self.potcars["pseudo"]))
            potcars = []
            
            for element in self.species:
                lel = element.lower()
                pp = element
                if lel in self.potcars:
                    pp += self.potcars[lel]

                ppot = path.join(potsrc, pp, "POTCAR")
                with open(ppot) as f:
                    first = f.readline()
                    if "version" in self.potcars:
                        assert self.potcars["version"] in first

                potcars.append(ppot)

            if len(potcars) == len(self.species):
                from matdb.utility import cat
                cat(potcars, target)
            else:
                msg.err("Couldn't create POTCAR for system.")

    @property
    def train_file(self):
        """Returns the full path to the XYZ database file that can be
        used for training.
        """
        return path.join(self.root, self._trainfile)

    @property
    def holdout_file(self):
        """Returns the full path to the XYZ database file that can be
        used to validate the potential fit.
        """
        return path.join(self.root, self._holdoutfile)

    def split(self, train_perc, trainfile="train.xyz", holdfile="holdout.xyz",
              recalc=False):
        """Splits the total available data in all databases into a
        training and holdout set.

        Args:
            train_perc (float): percentage of the data to use for training.
            trainfile (str): name of the training XYZ file to create.
            holdfile (str): name of the holdout/validation XYZ file to create.
            recalc (bool): when True, re-split the data and overwrite any
              existing *.xyz files.
        """
        if self._trainfile is None:
            self._trainfile = trainfile
        if self._holdoutfile is None:
            self._holdoutfile = holdfile

        if (path.isfile(self.train_file) and path.isfile(self.holdout_file)
            and not recalc):
            return
        
        trainfiles = []
        holdfiles = []
        for name, db in self.collections.items():
            db.split(train_perc, trainfile, holdfile, recalc)
            trainfiles.append(db.train_file)
            holdfiles.append(db.holdout_file)

        from matdb.utility import cat
        cat(trainfiles, self.train_file)
        cat(holdfiles, self.holdout_file)
