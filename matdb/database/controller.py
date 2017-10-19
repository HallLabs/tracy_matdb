"""Exposes classes and functions for interacting with the database
folders via a simple configuration file.
"""
from os import path
from matdb import msg
import numpy as np
import six

class DatabaseSequence(object):
    """Represents a sequence of databases (all inheriting from
    :class:`Database`) that are all related be the atomic
    configuration that they model.

    .. note:: See the list of attributes below which correspond to the sections
      in the YAML database specification file.

    Args:
        name (str): name of the configuration that this database sequence is
          operating for.
        repeater (SequenceRepeater): repeater for a set of sequences that all
          share the same root `atoms` object.
        root (str): root directory in which all other database sequences for
          the configurations in the same specification will be stored.
        parent (Controller): instance controlling multiple configurations.
        steps (list): of `dict` describing the kinds of sub-configuration
          database steps to setup.

    Attributes:
        title (str): title for the alloy system that these databases work with.
        config (str): name of the configuration that this database sequence is
          running for.
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
        steps (OrderedDict): keys are step names (e.g. `dft`, `calibration`,
          etc.); values are the corresponding class instances.
        parent (Controller): instance controlling multiple configurations.
    """
    def __init__(self, name, repeater, root, parent, steps):
        self.name = name
        self.config = name.split('.')[0]
        self.atoms = repeater.atoms.copy()
        self.root = path.join(root, name)
        self.repeater = repeater
        
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        parrefs = ["species", "potcars", "incar", "kpoints", "execution",
                   "plotdir"]
        for ref in parrefs:
            setattr(self, ref, getattr(parent, ref))
        self.parent = parent

        self._trainfile = None
        """str: name of the file that stores the sequence's XYZ *training*
        database.
        """
        self._holdoutfile = None
        """str: name of the file that stores the sequence's XYZ *validation*
        database.
        """        
        self._superfile = None
        """str: name of the file that stores the sequence's XYZ
        *super*-validation database.
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
        self._settings = steps
        """dict: with keys and values describing the kinds of step databases to setup.
        """

        from collections import OrderedDict
        self.steps = OrderedDict()
        for dbspec in steps:
            if isinstance(dbspec, six.string_types):
                #This is a reference to an existing database instance that was
                #defined previously.
                instance = self.parent[dbspec]
                self.steps[instance.name] = instance
                continue
                
            modname, clsname = dbspec["type"].split('.')
            fqdn = "matdb.database.{}".format(modname)
            module = import_module(fqdn)
            if not hasattr(module, clsname):# pragma: no cover
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
                if isinstance(cpspec[k], dict) and self.config in cpspec[k]:
                    cpspec[k] = cpspec[k][self.config]
            
            instance = cls(**cpspec)
            self.steps[instance.name] = instance

    @property
    def isteps(self):
        """Returns a generator over steps in this sequence. The generator yields
        the next step *only* if the previous one is already finished (i.e., the
        `ready()` method returns True.
        """
        previous = None
        for dbname, db in self.steps.items():
            if previous is None or previous[1].ready():
                previous = (dbname, db)
                yield previous
            else:
                raise StopIteration()
            
    def recover(self, rerun=False):
        """Runs recovery on this database to determine which configs failed and
        then create a jobfile to requeue them for compute.

        Args:
            rerun (bool): when True, recreate the jobfile even if it
              already exists. 
        """
        for dbname, db in self.steps.items():
            db.recover(rerun)
            
    def status(self, busy=False):
        """Prints a status message for each of the databases relative
        to VASP execution status.

        Args:
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
        """
        from matdb.msg import verbosity
        for dbname, db in self.isteps:
            if not busy:
                imsg = "{}:{} => {}".format(self.name, dbname, db.status(verbosity<2))
                msg.info(imsg)
            else:
                detail = db.status(False)
                running = [k for k, v in detail["done"].items() if not v]
                for config in running:
                    msg.std(config.replace(self.root, ""))
                
        msg.blank(level=1)
            
    def execute(self, recovery=False, env_vars=None):
        """Submits job array files for any of the databases that are ready to
        execute, but which haven't been submitted yet.

        Args:
            recovery (bool): when True, submit the script for running recovery
              jobs.
            env_vars (dict): of environment variables to set before calling the
              execution. The variables will be set back after execution.
        """
        ready = True
        for dbname, db in self.isteps:
            if ready:
                ready = (db.ready() or db.execute(recovery=recovery,
                                                  env_vars=env_vars))
            if not ready:
                imsg = ("Database {}.{} is not ready to execute yet, or is "
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
        from matdb.utility import safe_update
        safe_update(self, {"_trainfile": trainfile, "_holdoutfile": holdfile,
                           "_superfile": superfile})

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
            fi = 0
            for dbname, db in self.isteps:
                if not db.has_data:
                    continue
                    
                for dbconfigs in db.configs.values():
                    for f in dbconfigs:
                        subconfs[fi] = f
                        fi += 1
                        
            Ntot = len(subconfs)
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
        for dbname, db in self.isteps:
            if not db.cleanup():
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
        for dbname, db in self.isteps:
            msg.info("Setting up database {}:{}".format(self.name, dbname))
            db.setup(rerun)
        msg.blank()

class SequenceRepeater(object):
    """Repeats sequences of database steps across a parameter grid.

    Args:
        name (str): name of the configuration that this database sequence is
          operating for.
        poscar (str): name of the POSCAR file in `root` to extract atomic
          configuration information from.
        root (str): root directory in which all other database sequences for
          the configurations in the same specification will be stored.
        parent (Controller): instance controlling multiple configurations.
        steps (list): of `dict` describing the kinds of sub-configuration
          database steps to setup.

    Attributes:
        atoms (quippy.atoms.Atoms): a single atomic configuration from
          which many others may be derived using MD, phonon
          displacements, etc.
        poscar (str): path to the POSCAR file for the seed configuration that
          all sequences in this repeater will use.
        kfile (str): path to the kpath.json cached file for this repeater.        
    """
    def __init__(self, name, poscar, root, parent, steps, niterations=None):
        from collections import OrderedDict
        from copy import copy
        from quippy.atoms import Atoms
        
        self.name = name
        self.sequences = OrderedDict()
        self.poscar = path.join(root, poscar)
        self.atoms = Atoms(self.poscar, format="POSCAR")
        self.kfile = path.join(parent.kpathdir, "{0}.json".format(self.name))

        self._kpath = None
        """tuple: result of querying the materialscloud.org special path
        service. First term is a list of special point labels; second is the
        list of points corresponding to those labels.
        """
        
        if niterations is not None:
            from matdb.utility import obj_update
            
            for i, repeater in enumerate(niterations):
                nname = None
                isteps = copy(steps)
                for k, v in repeater.items():
                    if k == "suffix":
                        nname = self.name + v
                        continue
                    
                    obj_update(isteps, k, v, False)
                    
                if nname is None:
                    nname = self.name + "-{0:d}".format(i)

                iobj = DatabaseSequence(nname, self, root, parent, isteps)
                self.sequences[nname] = iobj
        else:
            single = DatabaseSequence(name, self, root, parent, steps)
            self.sequences[name] = single

    @property
    def kpath(self):
        """Returns the materialscloud.org special path in k-space for the seed
        configuration of this database.

        Returns:
            tuple: result of querying the materialscloud.org special path
            service. First term is a list of special point labels; second is the
            list of points corresponding to those labels.
        """
        if self._kpath is None:
            import json
            #We use some caching here so that we don't have to keep querying the
            #server and waiting for an identical response.
            if path.isfile(self.kfile):
                with open(self.kfile) as f:
                    kdict = json.load(f)
            else:
                from .phonon import _parsed_kpath
                labels, band = _parsed_kpath(self.poscar)
                kdict = {"labels": labels, "band": band}
                with open(self.kfile, 'w') as f:
                    json.dump(kdict, f)
            
            self._kpath = (kdict["labels"], kdict["band"])
            
        return self._kpath
            
    def recover(self, rerun=False):
        """Runs recovery on each step in the sequence to determine which configs failed
        and then create a jobfile to requeue them for compute.

        Args:
            rerun (bool): when True, recreate the jobfile even if it
              already exists.
        """
        for db in self.sequences.values():
            db.recover(rerun)

    def status(self, busy=False):
        """Prints a status message for each step in the sequence relative
        to VASP execution status.

        Args:
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
        """
        for db in self.sequences.values():
            db.status(busy)

    def execute(self, recovery=False, env_vars=None):
        """Submits job array files for any of steps in the sequence that are ready to
        execute, but which haven't been submitted yet.

        Args:
            recovery (bool): when True, submit the script for running recovery
              jobs.
            env_vars (dict): of environment variables to set before calling the
              execution. The variables will be set back after execution.
        """
        for db in self.sequences.values():
            db.execute(recovery, env_vars=env_vars)

    def split(self, train_perc, trainfile="train.xyz", holdfile="holdout.xyz",
              recalc=0, superfile="super.xyz"):
        """Splits the total available data in each step's databases into a training and
        holdout set.

        Args:
            train_perc (float): percentage of the data to use for training.
            trainfile (str): name of the training XYZ file to create.
            holdfile (str): name of the holdout/validation XYZ file to create.
            recalc (int): when non-zero, re-split the data and overwrite any
              existing *.xyz files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level XYZ files, increase this value.
        """
        for db in self.sequences.values():
            nrecalc = recalc - (0 if len(self.steps) == 0 else 1)
            db.split(train_perc, trainfile, holdfile, superfile, nrecalc)

    def cleanup(self):
        """Runs the cleanup methods of each step's databases in the sequence.
        """
        for db in self.sequences.values():
            db.cleanup()

    def setup(self, rerun=False):
        """Sets up the each step's database in the sequence in order. If an
        earlier step is not ready yet, the next item in the step won't setup.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """
        for db in self.sequences.values():
            db.setup(rerun)
            
class Controller(object):
    """Implements methods for tying a configuration dictionary (in
    YAML format) to instances of various databases.

    Args:
        config (str): path to the configuration YAML file that
          specifies all information for constructing the set of databases.
        tmpdir (str): path to a temporary directory to use for the
          database. This is for unit testing purposes.

    Attributes:
        specs (dict): the parsed settings from the YAML configuration file.
        collections (dict): keys are configuration names listed in attribute
          `configs` of the YAML file. Values are the :class:`DatabaseCollection`
          instances.
        plotdir (str): path to the directory to store plots in for all
          databases.
        kpathdir (str): path to the directory where cached k-paths are stored
          for root configurations.
    """
    def __init__(self, config, tmpdir=None):
        import yaml
        with open(config, 'r') as stream:
            self.specs = yaml.load(stream)

        #We allow the user to specify paths relative the matdb repo.
        from matdb.utility import relpath
        self.root = relpath(path.expanduser(self.specs["root"]))
        if tmpdir is not None:
            self.root = tmpdir
            
        self.plotdir = path.join(self.root, "plots")
        self.kpathdir = path.join(self.root, "kpaths")
        self.title = self.specs["title"]
        self.collections = {}
        self.species = [s for s in self.specs["species"]]
        self.potcars = self.specs["potcars"]
        self.incar = self.specs.get("incar", {})
        self.kpoints = self.specs.get("kpoints", {})
        self.execution = self.specs["execution"]

        for cspec in self.specs["configs"]:
            name, poscar = cspec["name"], cspec["poscar"]
            self.collections[name] = {}
            
            # We need to split out the databases by user-given name to create
            # the sequences.
            for dbspec in self.specs.get("databases", []):
                dbname = '.'.join((name, dbspec["name"]))
                steps = dbspec["steps"]
                seq = SequenceRepeater(dbname, poscar, self.root, self,
                                       steps, dbspec.get("niterations"))
                self.collections[name][dbspec["name"]] = seq

        from os import mkdir
        if not path.isdir(self.plotdir):
            mkdir(self.plotdir)
        if not path.isdir(self.kpathdir):
            mkdir(self.kpathdir)
            
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

    def ifiltered(self, cfilter=None, dfilter=None):
        """Returns a *filtered* generator over sequences in the config collections of
        this object.

        Args:
            cfilter (list): of `str` patterns to match against *configuration*
              names. This limits which configs are returned.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        from fnmatch import fnmatch
        for name, coll in self.collections.items():
            if cfilter is None or any(fnmatch(name, c) for c in cfilter):
                for dbn, seq in coll.items():
                    if dfilter is None or any(fnmatch(dbn, d) for d in dfilter):
                        yield (dbn, seq)

    def find(self, pattern):
        """Finds a list of database steps that match the given pattern. The
        pattern is formed using `config.dbname-suffix.step`. `*` can be used as
        a wildcard for any portion of the '.'-separated path.

        .. note:: Actually, an :func:`~fnmatch.fnmatch` pattern can be used.

        Args:
            pattern (str): fnmatch pattern that follows the convention of the DB
              key.

        Examples:
        
            Get all the dynamical matrix databases for the `Pd`
            configuration. The example assumes that the database name is
            `phonon` and that it includes a dynamical matrix step.

            >>> Pd = Controller("Pd.yml")
            >>> Pd.find("Pd.phonon*.dynmatrix")
        """
        from fnmatch import fnmatch
        config, parent, db = pattern.split('.')
        colls = [v for k, v in self.collections.items() if fnmatch(k, config)]

        #For databases without repeaters, there is no suffix.
        if '-' in parent:
            dbname, suffix = parent.split('-')
        else:
            dbname, suffix = parent, None

        steps = []
        for coll in colls:
            dbs = [dbi for dbn, dbi in coll.items() if fnmatch(dbn, dbname)]
            for dbi in dbs:
                seqs = [seqi for seqn, seqi in dbi.sequences.items()
                        if fnmatch(seqn, '.'.join((config, parent)))]
                for seq in seqs:
                    steps.extend([si for sn, si in seq.steps.items()
                                  if fnmatch(sn, db)])
        return steps

    def steps(self):
        """Compiles a list of all steps in this set of databases.
        """
        result = []
        for config, coll in self.collections.items():
            for repeater in coll.values():
                for parent, seq in repeater.sequences.items():
                    for step in seq.steps:
                        result.append("{0}.{1}".format(parent, step))

        return sorted(result)        
    
    def sequences(self):
        """Compiles a list of all sequences in this set of databases.
        """
        result = []
        for config, coll in self.collections.items():
            for repeater in coll.values():
                result.extend(repeater.sequences.keys())

        return sorted(result)        
    
    def __getitem__(self, key):
        """Returns the database object associated with the given key. This is
        necessary because of the hierarchy of objects needed to implement
        sequence repitition.
        """
        config, parent, db = key.split('.')
        coll = self.collections[config]
        if '-' in parent:
            dbname, suffix = parent.split('-')
        else:
            dbname = parent
        seq = coll[dbname].sequences['.'.join((config, parent))]
        return seq.steps[db]
                        
    def setup(self, rerun=False, cfilter=None, dfilter=None):
        """Sets up each of configuration's databases.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.
            cfilter (list): of `str` patterns to match against *configuration*
              names. This limits which configs are returned.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        for dbname, seq in self.ifiltered(cfilter, dfilter):
            seq.setup(rerun)

    def cleanup(self, cfilter=None, dfilter=None):
        """Runs cleanup on each of the configuration's databases.

        Args:
            cfilter (list): of `str` patterns to match against *configuration*
              names. This limits which configs are returned.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        for dbname, seq in self.ifiltered(cfilter, dfilter):
            seq.cleanup()

    def execute(self, recovery=False, cfilter=None, dfilter=None, env_vars=None):
        """Submits job array scripts for each database collection.

        Args:
            recovery (bool): when True, submit the script for running recovery
              jobs.
            cfilter (list): of `str` patterns to match against configuration
              names. This limits which configs status is retrieved for.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
            env_vars (dict): of environment variables to set before calling the
              execution. The variables will be set back after execution.
        """
        for dbname, seq in self.ifiltered(cfilter, dfilter):
            seq.execute(recovery, env_vars=env_vars)

    def recover(self, rerun=False, cfilter=None, dfilter=None):
        """Runs recovery on this database to determine which configs failed and
        then create a jobfile to requeue them for compute.

        Args:
            rerun (bool): when True, recreate the jobfile even if it
              already exists. 
            cfilter (list): of `str` patterns to match against configuration
              names. This limits which configs status is retrieved
              for.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        for dbname, seq in self.ifiltered(cfilter, dfilter):
            seq.recover(rerun) 
                
    def status(self, busy=False, cfilter=None, dfilter=None):
        """Prints status messages for each of the configuration's
        databases.

        Args:
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
            cfilter (list): of `str` patterns to match against configuration
              names. This limits which configs status is retrieved for.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        for dbname, seq in self.ifiltered(cfilter, dfilter):
            seq.status(busy)
            
    def POTCAR(self):
        """Creates the POTCAR file using the pseudopotential and version
        specified in the file.
        """
        from matdb.utility import relpath
        target = path.join(self.root, "POTCAR")
        if not path.isfile(target):
            # Make sure that the POTCAR version and pseudopotential type match
            # up so that we don't get nasty surprises.
            potsrc = path.join(relpath(self.potcars["directory"]),
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
                        if element in self.potcars["version"]:
                            version = self.potcars["version"][element]
                        else:
                            version = self.potcars["version"]
                        assert version in first

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
        from matdb.utility import safe_update
        safe_update(self, {"_trainfile": trainfile, "_holdoutfile": holdfile,
                           "_superfile": superfile})

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
