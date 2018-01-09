"""Exposes classes and functions for interacting with the database
folders via a simple configuration file.
"""
from os import path
from matdb import msg
import numpy as np
import six
import collections
from glob import glob
from uuid import uuid4

def parse_path(root,seeds,rseed=None):
    """Finds the full path to the seed files for this system.
    Args:
        root (str): the root directory for the databes.
        seeds (str or list of str): the seeds for the database that 
            need to be parsed and have the root folder found.
        rseed (optional): the seed for the random generator if used.
    
    Returns:
        seed_files (list): a list of the seed files for the database.
    """
    from matdb.utility import special_values
    from itertools import product
    
    seed_files = []
    for seed in seeds:
        # if there is a '/' in the seed then this is a path to the
        # seed file configuration, otherwise we assume the user put
        # the seed file in the `seed` directory in root.
        if len(seed.split("/")) > 1:
            this_seeds = []
            seed_path = root
            for segment in seed.split("/"):
                if "*" in segment:
                    if len(this_seeds) >=1:
                        this_level = []
                        for t_path in res:
                            this_level.extend(glob(path.join(seed_path,t_path,segment)))
                    else:
                        this_level = glob(path.join(seed_path,segment))
                else:
                    this_level = [segment]
                if len(this_seeds) >= 1:
                    this_seeds.extend([path.join(*i) for i in product(this_seeds,this_level)])
                else:
                    this_seeds.extend(this_level)                    

        else:
            seed_path = path.join(root,"seed")
            to_parse = seed
            if "*" in to_parse:
                this_seeds = glob(path.join(seed_path,to_parse))
            else:
                this_seeds = [to_parse]

        for ts in this_seeds:
            t_seed = path.join(seed_path,ts)
            if path.isfile(t_seed):
                seed_files.append(t_seed)
            else:
                msg.err("The seed file {} could not be found.".format(t_seed))

        return seed_files
    
def is_nested(d):
    """Determines if a dictoinary is nested, i.e. contains another dictionary.

    Args:
        d (dict): dictionary to test.
    """
    for k, v in d.items():
        if k[-1] == '*':
            return True
        elif isinstance(v, dict) and is_nested(v):
            return True

    return False

def get_suffix(d, k, index, values):
    """Returns the suffix for the specified key in the dictionary that
    is creating a parameter grid.

    Args:
        d (dict): the dictionary being turned into a grid.
        k (str): the key in the dictionary.
        index (int): the index for the value (gets used as the default suffix).
        values (str, list, float): the value for the parameter.
    """
    from matdb.utility import special_functions
    nk = k[0:-1]
    suffix = "{0}_suffix".format(nk)
    ssuff = suffix + '*'
    
    if suffix in d and (isinstance(d[suffix], dict) or ':' in d[suffix]):
        keyval = special_functions(d[suffix], values)
    elif suffix in d and isinstance(suffix, six.string_types):
        keyval = d[suffix].format(values)
    elif ssuff in d:
        keyval = d[ssuff][index]
    else:
        keyval = index
    
    if isinstance(keyval, float):
        return "{0}-{1:.2f}".format(nk[:3], keyval)
    else:
        return "{0}-{1}".format(nk[:3], keyval)

def get_grid(d, suffices=None):
    """Recursively generates a grid of parameters from the dictionary of parameters
    that has duplicates or wildcars in it. 
    
    Args:
       d (dict): the dictionary to be turned into a grid.
       suffices (list): an optional list of suffices.
    
    Returns:
       A dictionary of (key: value) where the key is the suffix string for
       the parameters and the value are the exact parameters for each
       entry in the grid.
    """
    dcopy = d.copy()
    stack = [(dcopy, None)]
    result = {}
    
    if suffices is None:
        suffices = {k: v for k, v in d.items() if "suffix" in k[-8:]}
        for k in suffices:
            del dcopy[k]
    else:
        for k,v in d.items():
            if "suffix" in k[-8:]:
                suffices[k] = v
        for k in suffices:
            if k in dcopy:
                del dcopy[k]            
        
    while len(stack) > 0:
        oned, nsuffix = stack.pop()
        for k, v in sorted(oned.items()):
            if k[-1] == '*':
                nk = k[0:-1]                    
                for ival, value in enumerate(v):
                    suffix = get_suffix(suffices, k, ival+1, value)
                    dc = oned.copy()
                    del dc[k]
                    dc[nk] = value
                    compsuffix = suffix if nsuffix is None else '-'.join(map(str, (nsuffix, suffix)))
                    stack.append((dc, compsuffix))
                break
            elif isinstance(v, dict) and is_nested(v):
                blowup = get_grid(v, suffices)
                for rsuffix, entry in sorted(blowup.items()):
                    dc = oned.copy()
                    dc[k] = entry
                    compsuffix = rsuffix if nsuffix is None else '-'.join(map(str, (nsuffix, rsuffix)))
                    stack.append((dc, compsuffix))
                break
        else:
            result[nsuffix] = oned
            
    return result

class ParameterGrid(collections.MutableSet):
    """An ordered list of the paramater combinations for the database. 
    Values are the suffixes of the combinations of parameters as tuples:
    e.g. (8, "dog", 1.2) for "dim", "animal", "temperature"
    ({"animal*": ["dog", "cat", "cow"], "dim*": [[],[],[]], "temperature": 1.2})
    Args:
        params (dict): the paramaters needed to build the database.
    
    Attributes:
        values (dict): keys are the suffix tuple and the values are the 
            actual values needed by the database.
        keys (list): the `str` names of the different parameters in the database.
    """
    def __init__(self, params):
        for k in ["root","parent","atoms"]:
            if k in params:
                params.pop(k)
        grid = get_grid(params)
        #add these items to the set.
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        self.values = {}
        self.params = params
        for i, v in grid.items():
            self.add(i,v)
            
    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def __getitem__(self, key):
        return self.values[key]
                        
    def add(self, key, value):
        """Adds key to the set if it is not already in the set.
        Args:
            key (tuple): Anything that could be added to the set.
            value (tuple): The actual values that the suffix's 
                correspond to.
        """
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]
            self.values[key] = value
        else:
            msg.warn("The key {} already exists in the set, ignoring addition.".format(key))

    def discard(self, key):
        """Removes the key from the set.
        Args:
            key (tuple): An element of the set.
        """        
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev
            self.values.pop(key,None)

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            if curr[0] is not None:
                yield curr[0]
            curr = curr[2]

    def pop(self, key):
        """Removes an element from the set.
        Args:
            key (tuple): An element of the set.
        """
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, ParameterGrid):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

class Database(object):
    """Represents a Database of groups (all inheriting from :class:`Group`) that 
    are all related be the atomic configuration that they model.
    .. note:: See the list of attributes below which correspond to the sections
      in the YAML database specification file.
    Args:
        name (str): name of the configuration that this database sequence is
          operating for.
        root (str): root directory in which all other database sequences for
          the configurations in the same specification will be stored.
        parent (Controller): instance controlling multiple configurations.
        steps (list): of `dict` describing the kinds of sub-configuration
          database steps to setup.
        splits (dict): keys are split names; values are `float` *training*
          percentages to use.
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
    def __init__(self, name, root, parent, steps, splits):
        self.name = name
        self.config = name.split('.')[0]
        self.root = root
        self.splits = {} if splits is None else splits
        
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        parrefs = ["species", "execution", "plotdir", "calculator"]
        for ref in parrefs:
            setattr(self, ref, getattr(parent, ref))
        self.parent = parent
        
        from importlib import import_module
        self._settings = steps
        """dict: with keys and values describing the kinds of step databases to setup.
        """

        from collections import OrderedDict
        from os import mkdir
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
            
            cpspec["pgrid"] = ParameterGrid(cpspec.copy())
            if len(cpspec["pgrid"]) ==0:
                cpspec["parameters"] = cpspec["pgrid"].params
            for k in list(cpspec.keys()):
                if "suffix" in k:
                    del cpspec[k]
                elif "*" == k[-1]:
                    cpspec[k[:-1]] = None
                    del cpspec[k]
            
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
                imsg = ("Group {}.{} is not ready to execute yet, or is "
                        "already executing. Done.")
                msg.info(imsg.format(self.name, dbname))
                break
        msg.blank()

    def train_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used for training.
        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{0}-train.xyz".format(split))

    def holdout_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used to validate the potential fit.
        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{0}-holdout.xyz".format(split))

    def super_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used to *super* validate the potential fit.
        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{0}-super.xyz".format(split))

    def split(self, recalc=0):
        """Splits the database multiple times, one for each `split` setting in
        the database specification.
        """
        for name in self.splits:
            self._split(name, recalc)
    
    def _split(self, name, recalc=0):
        """Splits the total available data in all databases into a training and holdout
        set.
        Args:
            name (str): name of the split to perform.
            recalc (int): when non-zero, re-split the data and overwrite any
              existing files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level XYZ files, increase this value.
        """
        train_file = self.train_file(name)
        holdout_file = self.holdout_file(name)
        super_file = self.super_file(name)
        if (path.isfile(train_file) and path.isfile(holdout_file)
            and path.isfile(super_file) and recalc <= 0):
            return

        train_perc = self.splits[name]
        
        #Compile a list of all the sub-configurations we can include in the
        #training.
        from cPickle import dump, load
        idfile = path.join(self.root, "{0}-ids.pkl".format(name))
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
                if len(db.rset) == 0 or not db.trainable:
                    continue
                    
                for configpath in db.rset.values():
                    subconfs[fi] = configpath
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
                "uuid": uuid4(),
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

        def subset(subconfs, idlist, recalc):
            from tqdm import tqdm
            files = []
            for aid in tqdm(idlist):
                target = path.join(subconfs[aid], "output.xyz")
                if path.isfile(target):
                    files.append(target)

            return files

        trainfiles = subset(subconfs, tids, recalc)
        holdfiles = subset(subconfs, hids, recalc)
        superfiles = subset(subconfs, sids, recalc)

        #We don't use dbcat for the concatenation because there are hundreds of
        #files and the pkl file above keeps track of all the splits and ids,
        #etc. Instead, we dbcat after the fact to create the configuration files
        #and version numbers.        
        from matdb.utility import cat, dbcat
        cat(trainfiles, train_file)
        cat(holdfiles, holdout_file)
        cat(superfiles, super_file)

        sources = self.steps.keys()
        dbcat([], train_file, sources=sources, N=Ntrain)
        dbcat([], holdout_file, sources=sources, N=Nhold)
        dbcat([], super_file, sources=sources, N=Nsuper)
        
    def cleanup(self):
        """Runs the cleanup methods of each database in the collection, in the
        correct order.
        """
        for dbname, db in self.isteps:
            if not db.cleanup():
                imsg = "Group {}:{} is not ready yet. Done."
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
        
class Controller(object):
    """Implements methods for tying a configuration dictionary (in
    YAML format) to instances of various databases.
    Args:
        config (str): name of the YML file (without the .yml) that
          specifies all information for constructing the set of databases.
        tmpdir (str): path to a temporary directory to use for the
          database. This is for unit testing purposes.

    Attributes:
        specs (dict): the parsed settings from the YAML configuration file.
        collections (dict): keys are configuration names listed in attribute
          `configs` of the YAML file. Values are the :class:`DatabaseCollection`
          instances.
        legacy (dict): keys are legacy database names; values are
          :class:`~matdb.database.legacy.LegacyDatabase`.
        plotdir (str): path to the directory to store plots in for all
          databases.
        venv (str): name of a virtual environment to activate for plotting
          potentials after fitting.
    """
    def __init__(self, config, tmpdir=None):
        from matdb.io import read
        self.config = path.expanduser(path.abspath(config))
        if path.isabs(config):
            root, config = path.split(config)
        else:
            root, config = path.dirname(self.config), config
        self.specs = read(root, config)

        #We allow the user to specify paths relative the matdb repo.
        from matdb.utility import relpath
        self.root = relpath(path.expanduser(self.specs["root"]))
        if tmpdir is not None:
            self.root = tmpdir
            
        self.plotdir = path.join(self.root, "plots")
        self.title = self.specs["title"]
        self.legacy = {}
        self.collections = {}
        self.species = [s for s in self.specs["species"]]
        self.execution = self.specs.get("execution", {})
        self.calculator = self.specs.get("calculator", {})
        self.potcars = self.specs["potcars"]
        if "Vasp" in self.calculator["name"]:
            from os import environ
            environ["VASP_PP_PATH"] = relpath(path.expanduser(self.potcars["directory"]))
            from matdb import calculators
            from ase import Atoms, Atom
            calcargs = self.calculator.copy()
            calc = getattr(calculators, calcargs["name"])
            del calcargs["name"]
            potargs = self.potcars.copy()
            del potargs["directory"]
            calcargs.update(potargs)
            elems = sorted(self.species)
            this_atom = Atoms([Atom(a,[0,0,i]) for i,a
                               in enumerate(elems)],cell=[1,1,len(elems)+1])
            calc = calc(this_atom,self.root,**calcargs)
            calc.write_potcar(directory=self.root)

        self.venv = self.specs.get("venv")
        self.random_seed = self.specs.get("random seed")

        # We need to split out the databases by user-given name to create
        # the sequences.
        from matdb.database.legacy import LegacyDatabase
        for dbspec in self.specs.get("databases", []):
            if dbspec.get("legacy", False):
                cpspec = dbspec.copy()
                cpspec["root"] = self.root
                cpspec["controller"] = self
                cpspec["splits"] = self.specs.get("splits")
                #We allow the user to specify the folder relative to repository
                #root; this is mainly for unit tests.
                cpspec["folder"] = relpath(cpspec["folder"])
                del cpspec["legacy"]
                self.legacy[cpspec["name"]] = LegacyDatabase(**cpspec)
            else:
                dbname = dbspec["name"]
                if dbname not in self.collections:
                    self.collections[dbname] = {}
                steps = dbspec["steps"]
                db = Database(dbname, self.root, self,
                              steps, self.specs.get("splits"))
                self.collections[dbname][dbspec["name"]] = db

        from os import mkdir
        if not path.isdir(self.plotdir):
            mkdir(self.plotdir)
            
        #Extract the POTCAR that all the databases are going to use. TODO: pure
        #elements can't use this POTCAR, so we have to copy the single POTCAR
        #directly for those databases; update the create().
        # if "potcars" in self.specs:
        #     self.potcars = self.specs["potcars"]
        #     self.POTCAR()

        #If the controller is going to train any potentials, we also need to 
        self.trainers = None
        if "fitting" in self.specs:
            from matdb.fitting.controller import TController
            tdict = self.specs["fitting"].copy()
            tdict["root"] = self.root
            tdict["db"] = self
            self.trainers = TController(**tdict)

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

    def relpaths(self, pattern):
        """Finds the relative paths for the seed configurations within the databases that 
        match to the pattern.
        Args:
            pattern (str): the pattern to match.
        """
        
        return parse_path(self.root,pattern,rseed=self.random_seed)
    
    def find(self, pattern):
        """Finds a list of :class:`matdb.database.basic.Group` that match the given
        pattern. The pattern is formed using `group.dbname[[.seed].params]`. `*`
        can be used as a wildcard for any portion of the '.'-separated path.
        .. note:: Actually, an :func:`~fnmatch.fnmatch` pattern can be used.
        Args: pattern (str): fnmatch pattern that follows the convention of the
        DB key.  Examples:
        
            Get all the dynamical matrix databases for the `Pd`
            configuration. The example assumes that the database name is
            `phonon` and that it includes a dynamical matrix step.
            >>> Pd = Controller("Pd.yml")
            >>> Pd.find("DynMatrix.phonon.Pd.*")
            Get all the database sequences for liquids across all configurations in
            the database.
            >>> CdWO4 = Controller("CdWO4.yml")
            >>> CdWO4.find("*.liquid*")
        """
        if pattern == '*':
            return self.find('*.*')
        
        from fnmatch import fnmatch
        if pattern.count('.') == 2:
            parent, db, config = pattern.split('.')
        elif pattern.count('.') == 1:
            parent, config = pattern.split('.')
            db = None
        else:
            #We must be searching legacy databases; match the pattern against
            #those.
            return [li for ln, li in self.legacy.items() if fnmatch(ln, pattern)]
        
        colls = [v for k, v in self.collections.items() if fnmatch(k, config)]

        #For databases without repeaters, there is no suffix.
        if '-' in parent:
            dbname, suffix = parent.split('-')
        else:
            dbname, suffix = parent, None

        result = []
        for coll in colls:
            dbs = [dbi for dbn, dbi in coll.items() if fnmatch(dbn, dbname)]
            for dbi in dbs:
                seqs = [seqi for seqn, seqi in dbi.items()
                        if fnmatch(seqn, '.'.join((config, parent)))]

                if db is not None:
                    for seq in seqs:
                        result.extend([si for sn, si in seq.steps.items()
                                       if fnmatch(sn, db)])
                else:
                    result.extend(seqs)

        if config == '*':
            #Add all the possible legacy databases.
            result.extend([li for ln, li in self.legacy.items()
                           if fnmatch(ln, config)])
                    
        return result

    def steps(self):
        """Compiles a list of all steps in this set of databases.
        """
        result = []
        for config, coll in self.collection.items():
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
        if key.count('.') == 2:
            config, parent, db = key.split('.')
        else:
            config, parent = key.split('.')
            db = None
            
        coll = self.collections[config]
        if '-' in parent:
            dbname, suffix = parent.split('-')
        else:
            dbname = parent
        seq = coll[dbname].sequences['.'.join((config, parent))]

        if db is not None:
            return seq.steps[db]
        else:
            return seq
                        
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
            else: #pragma: no cover
                msg.err("Couldn't create POTCAR for system.")

    def split(self, recalc=0, cfilter=None, dfilter=None):
        """Splits the total available data in all databases into a training and holdout
        set.
        Args:
            recalc (int): when non-zero, re-split the data and overwrite any
              existing *.xyz files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level XYZ files, increase this value.
            cfilter (list): of `str` patterns to match against *configuration*
              names. This limits which configs are returned.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        for dbname, seq in self.ifiltered(cfilter, dfilter):
            seq.split(recalc)
