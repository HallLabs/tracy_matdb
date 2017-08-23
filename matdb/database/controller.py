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
          (like PhononDFT and PhononCalibration) don't provide configurations
          for training and testing.
        byorder (dict): keys are database type qualified names; values
          are a list of actual database names of that type.
        clsorder (list): of database type qualified names and the
          order in which they should be processed.
    """
    clsorder = ["phonon.PhononDFT", "phonon.PhononCalibration",
                "phonon.PhononDatabase", "md.DynamicsDatabase",
                "liquid.LiquidDatabase"]
    configdbs = ["phonon.PhononDatabase", "liquid.LiquidDatabase"]
    
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
        self._superfile = None
        """str: name of the file that stores the XYZ *super*-validation
        database.
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
                if "super" in xyzname:
                    self._superfile = xyzname
        
        from importlib import import_module
        self._dbsettings = database
        """dict: with keys and values describing the kinds of sub-configuration
        databases to setup.
        """

        self.databases = {}
        self.bytype = {k: [] for k in self.clsorder}        
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
            self.bytype[dbspec["type"]].append(instance.name)

        #Finally, set the order that the databases should be processed
        #in by iterating over names from each database type.
        self.order = []                
        for dbtype in self.clsorder:
            self.order.extend(self.bytype[dbtype])

    def recover(self):
        """Runs recovery on this database to determine which configs failed and
        then create a jobfile to requeue them for compute.
        """
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            db = self.databases[dbname]
            db.recover()
            
    def status(self, busy=False):
        """Prints a status message for each of the databases relative
        to VASP execution status.

        Args:
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
        """
        from matdb.msg import verbosity
        for dbname in self.order:
            if dbname not in self.databases:
                continue

            db = self.databases[dbname]
            if not busy:
                imsg = "{}:{} => {}".format(self.name, dbname, db.status(verbosity<2))
                msg.info(imsg)
            else:
                detail = db.status(False)
                running = [k for k, v in detail["done"].items() if not v]
                for config in running:
                    msg.std(config.replace(self.root, ""))
                
        msg.blank(level=1)
            
    def execute(self, recovery=False):
        """Submits job array files for any of the databases that are ready to
        execute, but which haven't been submitted yet.

        Args:
            recovery (bool): when True, submit the script for running recovery
              jobs.
        """
        ready = True
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            if ready:
                ready = (self.databases[dbname].ready() or
                         self.databases[dbname].execute(recovery=recovery))
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

    @property
    def super_file(self):
        """Returns the full path to the XYZ database file that can be
        used to *super* validate the potential fit.
        """
        return path.join(self.root, self._superfile)
    
    def split(self, train_perc, trainfile="train.xyz", holdfile="holdout.xyz",
              recalc=0, superfile="super.xyz"):
        """Splits the total available data in all databases into a
        training and holdout set.

        Args:
            train_perc (float): percentage of the data to use for training.
            trainfile (str): name of the training XYZ file to create.
            holdfile (str): name of the holdout/validation XYZ file to create.
            recalc (int): when non-zero, re-split the data and overwrite any
              existing *.xyz files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level XYZ files, increase this value.
        """
        if self._trainfile is None:
            self._trainfile = trainfile
        if self._holdoutfile is None:
            self._holdoutfile = holdfile
        if self._superfile is None:
            self._superfile = superfile

        if (path.isfile(self.train_file) and path.isfile(self.holdout_file)
            and path.isfile(self.super_file) and recalc <= 0):
            return
        
        #Compile a list of all the sub-configurations we can include in the
        #training.
        from cPickle import dump, load
        idfile = path.join(self.root, "ids.pkl")
        if path.isfile(idfile):
            with open(idfile, 'rb') as f:
                data = load(f)
            subconfs = data["subconfs"]
            ids = data["ids"]
            Ntrain = data["Ntrain"]
            Nhold = data["Nhold"]
            Ntot = data["Ntot"]
            Nsuper = data["Nsuper"]
        else:
            Ntot = 0
            subconfs = {}
            for dbtype in self.configdbs:
                relevant = [self.databases[n].configs.values()
                            for n in self.bytype[dbtype]]
                for dbconfigs in relevant:
                    for fi, f in enumerate(dbconfigs):
                        subconfs[fi + Ntot] = f
                    Ntot += len(subconfs)

            Ntrain = int(np.ceil(Ntot*train_perc))
            ids = np.arange(len(subconfs))
            Nhold = int(np.ceil((Ntot-Ntrain)*train_perc))
            Nsuper = Ntot-Ntrain-Nhold
            np.random.shuffle(ids)

            #We need to save these ids so that we don't mess up the statistics on
            #the training and validation sets.
            data = {
                "subconfs": subconfs,
                "ids": ids,
                "Ntrain": Ntrain,
                "Nhold": Nhold,
                "Ntot": Ntot,
                "Nsuper": Nsuper
            }
            with open(idfile, 'wb') as f:
                dump(data, f)
        
        tids = ids[0:Ntrain]
        hids = ids[Ntrain:-Nsuper]
        sids = ids[-Nsuper:]

        def subset(subconfs, idlist):
            from matdb.io import vasp_to_xyz
            from tqdm import tqdm
            files = []
            for aid in tqdm(idlist):
                if vasp_to_xyz(subconfs[aid], recalc=recalc-1):
                    files.append(path.join(subconfs[aid], "output.xyz"))

            return files

        trainfiles = subset(subconfs, tids)
        holdfiles = subset(subconfs, hids)
        superfiles = subset(subconfs, sids)

        from matdb.utility import cat
        cat(trainfiles, self.train_file)
        cat(holdfiles, self.holdout_file)
        cat(superfiles, self.super_file)
        
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
            
    def setup(self, rerun=False):
        """Sets up the database collection by generating the POTCAR file and
        initializing any databases that haven't already been initialized.

        .. note:: The db setup functions are setup to only execute once, and then
           only if their dependencies have completed their calculations. This
           method can, therefore, be safely called repeatedly between different
           terminal sessions.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """
        for dbname in self.order:
            if dbname not in self.databases:
                continue
            db = self.databases[dbname]
            msg.info("Setting up database {}:{}".format(self.name, dbname))
            db.setup(rerun)
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
        self._superfile = None
        """str: name of the file that stores the XYZ *super*-validation
        database.
        """        

            
    def setup(self, rerun=False):
        """Sets up each of configuration's databases.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.
        """
        for name, coll in self.collections.items():
            coll.setup(rerun)

    def cleanup(self):
        """Runs cleanup on each of the configuration's databases.
        """
        for name, coll in self.collections.items():
            coll.cleanup()

    def execute(self, cfilter=None, recovery=False):
        """Submits job array scripts for each database collection.

        Args:
            cfilter (list): of `str` patterns to match against configuration
              names. This limits which configs status is retrieved for.
            recovery (bool): when True, submit the script for running recovery
              jobs.
        """
        from fnmatch import fnmatch
        for name, coll in self.collections.items():
            if cfilter is None or any(fnmatch(name, p) for p in cfilter):
                coll.execute(recovery)

    def recover(self, cfilter=None):
        """Runs recovery on this database to determine which configs failed and
        then create a jobfile to requeue them for compute.

        Args:
            cfilter (list): of `str` patterns to match against configuration
              names. This limits which configs status is retrieved for.
        """
        from fnmatch import fnmatch
        for name, coll in self.collections.items():
            if cfilter is None or any(fnmatch(name, p) for p in cfilter):
                coll.recover() 
                
    def status(self, cfilter=None, busy=False):
        """Prints status messages for each of the configuration's
        databases.

        Args:
            cfilter (list): of `str` patterns to match against configuration
              names. This limits which configs status is retrieved for.
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
        """
        from fnmatch import fnmatch
        for name, coll in self.collections.items():
            if cfilter is None or any(fnmatch(name, p) for p in cfilter):
                coll.status(busy) 
            
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

    @property
    def super_file(self):
        """Returns the full path to the XYZ database file that can be
        used to *super* validate the potential fit.
        """
        return path.join(self.root, self._superfile)
    
    def split(self, train_perc, trainfile="train.xyz", holdfile="holdout.xyz",
              recalc=0, superfile="super.xyz"):
        """Splits the total available data in all databases into a
        training and holdout set.

        Args:
            train_perc (float): percentage of the data to use for training.
            trainfile (str): name of the training XYZ file to create.
            holdfile (str): name of the holdout/validation XYZ file to
              create.
            recalc (int): when non-zero, re-split the data and overwrite any
              existing *.xyz files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level XYZ files, increase this value.
        """
        if self._trainfile is None:
            self._trainfile = trainfile
        if self._holdoutfile is None:
            self._holdoutfile = holdfile
        if self._superfile is None:
            self._superfile = superfile

        if (path.isfile(self.train_file) and path.isfile(self.holdout_file)
            and path.isfile(self._superfile) and recalc <= 0):
            return
        
        trainfiles = []
        holdfiles = []
        superfiles = []
        from fnmatch import fnmatch
        for name, db in self.collections.items():
            if self.trainer.configs is not None:
                #Skip all those databases that are not explicitly included in
                #the training spec.
                if not any([fnmatch(name, p) for p in self.trainer.configs]):
                    continue
            
            db.split(train_perc, trainfile, holdfile, recalc-1)
            trainfiles.append(db.train_file)
            holdfiles.append(db.holdout_file)
            superfiles.append(db.super_file)

        from matdb.utility import cat
        cat(trainfiles, self.train_file)
        cat(holdfiles, self.holdout_file)
        cat(superfiles, self.super_file)
