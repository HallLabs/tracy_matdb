"""Implements a basic trainer that train material models.
"""
import abc
from matdb import atoms
from os import path
from matdb import msg

class Trainer(object):
    """Represents an algorithm that can use a :class:`DatabaseSequence` (or set of
    sequences) to produce a multi-component potential.

    .. note:: This object provides the following methods that must be overriden:
    
      1. :meth:`command` must setup any dependencies in the folder that the
         training executable requires.
      2. :meth:`status` should return either a dictionary or a printed status
         message describing how things are going with the fitting.
      3. :meth:`get_calculator` should return an :ase:`Calculator` that can
         compute energies, forces and virial tensors. **IMPORTANT**: this method
         should be able to produce a calculator using only whatever computation
         is done by :meth:`command` and :meth:`execute`.

      The basic workflow is:

      1. Call :meth:`jobfile` to construct a jobfile for HPC resources. This
         method calls :meth:`command` as part of its templating procedures.
      2. Call :meth:`execute`, which is a generic method that executes whatever
         command is in the `jobfile`.
      3. Use :attr:`calculator` to calculate energies etc. This method
         internally calls :meth:`get_calculator` only once and then caches the
         object for the lifetime of the Trainer.

    .. note:: For the `dbfilter` parameter, you can use one of the standard
      comparison `operator` ['>', '<', '==', '!=', '>=', '<=']; `value` should be the RHS
      of the operator. Finally, `dbs` specifies a list of databases to apply the
      filter to.

    Args:
        name (str): name of this model; used as the output folder name.
        controller (matdb.fitting.controller.TController): fitting controller
          whose data will be used to train the potentials.
        dbs (list): of `str` database patterns from the `db` that should be
          included in the training and validation.
        execution (dict): settings needed to configure the jobfile for running
          the fit.
        parent (TrainingSequence): training sequence that this trainer belongs
          to.
        dbfilter (dict): keys are attributes on individual :class:`Atoms`
          objects; values are dictionaries with keys `operator`, `value` and
          `dbs` as described above.

    Attributes:
        fqn (str): fully-qualified name of the trainer.
        cust_splits (dict): keys are database sequence names; values are custom
          splits that should be used for that particular database sequence.
    """
    def __init__(self, controller=None, dbs=None, execution={}, split=None,
                 root=None, parent=None, dbfilter=None):
        self.controller = controller
        self.parent = parent
        self.fqn = "{}.{}".format(self.parent.name, self.name)
        self.execution = {} if execution is None else execution.copy()
        self.split = split
        self.params = {}
        self._dbs = ['*.*'] if dbs is None else dbs
        if not isinstance(self._dbs, (list, set, tuple)):
            self._dbs = [self._dbs]

        #Configure the fitting directory for this particular set of
        #potentials. This way, we can get separate directories if the parameters
        #change.
        from os import mkdir
        self.root = path.join(root, self.name)
        if not path.isdir(self.root):
            mkdir(self.root)

        #Find all the database sequences that match the patterns supplied to us.
        self.dbs = []
        self.cust_splits = {}
        self._find_dbs()

        #Invert the dbfilter so that it is keyed by database; this is for
        #optimization purposes.
        self.dbfilter = {} if dbfilter is None else dbfilter
        self._dbfilters = {}
        """dict: keys are database names; values are also dictionaries but with
        keys being :class:`Atoms` attribute names and values functions that
        return a bool based on a reference comparison value.
        """
        
        #We allow plotting to use split percentages; calculate the average split
        #percentage.
        self._update_split_params()
            
        self._calculator = None
        """ase.Calculator: potential for the fitted file in this Trainer.
        """
        self._trainfile = path.join(self.root, "train.xyz")
        """str: path to the XYZ training file that will be passed to the training
        command.
        """
        self._holdfile = path.join(self.root, "holdout.xyz")
        """str: path to the XYZ validation file that will be passed to the training
        command.
        """
        self._superfile = path.join(self.root, "super.xyz")
        """str: path to the XYZ super validation file that will be passed to the
        training command.
        """

        self._jobfile = path.join(self.root, "jobfile.sh")
        """str: path to the jobfile for the current trainer.
        """

    def _find_dbs(self):
        """Finds all dbs that match the patterns supplied for training.
        """
        for dbpat in self._dbs:
            if '/' in dbpat:
                pat, split = dbpat.split('/')
            else:
                pat, split = dbpat, None

            matches = self.controller.db.find(pat)
            self.dbs.extend(matches)
            if split is not None:
                for match in matches:
                    self.cust_splits[match.name] = split        
        
    def _invert_filters(self):
        """Inverts any available db filters for this trainer for optimization
        purposes.
        """
        from fnmatch import fnmatch
        from six import string_types
        import operator
        from matdb.utility import getattrs
        
        alldbs = [db.name for db in self.dbs]        
        for attr, spec in self.dbfilter.items():
            if "dbs" in spec:
                dbs = [db for db in alldbs
                       for dbpat in spec["dbs"]
                       if fnmatch(db, dbpat)]
            else:
                dbs = alldbs

            vals = {}
            for k, v in spec.items():
                if k not in ["operator", "dbs"]:
                    vals[k] = v
                if isinstance(v, string_types):
                    if len(v) > 0 and v[0] == '|':
                        otype, oname, chain = v[1:].split('|')
                        if otype == "db":
                            obj = self.controller.db.find(oname)
                        elif otype == "ip":
                            if oname == "self":
                                obj = [self]
                            else:
                                obj = self.controller.find(oname)

                        if len(obj) == 0:
                            emsg = "Cannot find object {} with type {}."                            
                            raise ValueError(emsg.format(oname, otype))
                        vals[k] = getattrs(obj[0], chain)

            estr = ("{__v}" + spec["operator"])
            for db in dbs:
                if db not in self._dbfilters:
                    self._dbfilters[db] = {}

                opf = lambda v: eval(estr.format(__v=v, **vals))
                self._dbfilters[db][attr] = (opf, vals)
        
    def _update_split_params(self):
        """Updates the parameters set for the trainer to include the splits.
        """
        if len(self.dbs) > 0:
            _splitavg = []
            for db in self.dbs:
                if db in self.cust_splits:
                    splt = self.cust_splits[db]
                    _splitavg.append(1 if splt == '*' else splt)
                else:
                    _splitavg.append(db.database.splits[self.split])
            self.params["split"] = sum(_splitavg)/len(self.dbs)
        
    def compile(self):
        """Compiles the cumulative database for this particular fit.
        """
        if path.isfile(self._trainfile):
            #No need to recompile the file. Since the splits are defined
            #globally and cached, we don't have to worry about inconsistencies.
            return
        
        #Setup the train.xyz file for the set of databases specified in the
        #source patterns.
        self._invert_filters()
        self.configs("train", False)
        
    @abc.abstractmethod
    def get_calculator(self):
        """Constructs an :ase:`Calculator` instance using the fitted potential
        in this :class:`Trainer` sub-class.
        """
        pass

    @abc.abstractmethod
    def command(self):
        """Returns the command that is needed to train the potentials specified by this
        object.

        .. note:: This method *must* also configure the directory that the
          command will run in so that it has the relevant files.
        """
        pass

    def extras(self):
        """Returns a list of extra XYZ training files that should be included in
        the database training. This is used by some trainers that create
        additional configurations as part of the training.

        .. note:: This method returns an empty list; override it if necessary.
        """
        return []
    
    def ready(self):
        """Returns True if the training step is complete. This method usually
        just tests whether the :attr:`calculator` is a valid object. However,
        for training steps that don't produce calculators, this method can be
        overridden.
        """
        return self.calculator is not None
    
    @abc.abstractmethod
    def status(self, printed=True):
        """Returns or prints the current status of the training.

        Args:
            printed (bool): when True, print the status to screen; otherwise,
              return a dictionary with the relevant quantities.
        """
        pass
    
    @property
    def calculator(self):
        """Returns an instance of :class:`ase.Calculator` using the latest
        fitted GAP potential in this trainer.
        """
        if self._calculator is None:
            from matdb.utility import chdir
            with chdir(self.root):
                self._calculator = self.get_calculator()
        return self._calculator
            
    @property
    def validation(self):
        """Returns a :class:`atoms.AtomsList` of configurations that can be
        used for potential validation.
        """
        return self.configs("holdout")

    def _filter_dbs(self, seqname, dbfiles):
        """Filters each of the database files specified so that they conform to
        any specified filters.

        Args:
            seqname (str): name of the sequence that the database files are
              from.
            dbfiles (list): of `str` paths to database files to filter.

        Returns:
            list: of `str` paths to include in the database from this sequence.
        """
        from matdb.utility import dbcat
        if len(self.dbfilter) > 0 and seqname in self._dbfilters:
            filtered = []
            #The filters values have a function and a list of the actual values
            #used in the formula replacement. Extract the parameters; we can't
            #serialize the actual functions.
            filters = self._dbfilters[seqname].items()
            params = {k: v[1] for k, v in filters}
            
            for dbfile in dbfiles:
                dbname = path.basename(path.dirname(dbfile))
                filtdb = path.join(self.root, "__{}.xyz".format(dbname))
                if path.isfile(filtdb):
                    continue
                
                al = atoms.AtomsList(dbfile)
                nl = atoms.AtomsList()
                for a in al:
                    #The 0 index here gets the function out; see comment above
                    #about the filters dictionary.
                    if not any(opf[0](getattr(a, attr)) for attr, opf in filters):
                        nl.append(a)

                if len(nl) != len(al):
                    nl.write(filtdb)
                    dN, N = (len(al)-len(nl), len(nl))
                    dbcat([dbfile], filtdb, filters=params, dN=dN, N=N)
                    filtered.append(filtdb)
                else:
                    filtered.append(nfile)
        else:
            filtered = dbfiles

        return filtered
    
    def configs(self, kind, asatoms=True):
        """Loads a list of configurations of the specified kind.

        Args:
            kind (str): on of ['train', 'holdout', 'super'].
            asatoms (bool): when True, return a :class:`atoms.AtomsList`
              object; otherwise just compile the file.

        Returns:
            atoms.AtomsList: for the specified configuration class.
        """
        fmap = {
            "train": lambda seq, splt: seq.train_file(splt),
            "holdout": lambda seq, splt: seq.holdout_file(splt),
            "super": lambda seq, splt: seq.super_file(splt)
        }
        smap = {
            "train": self._trainfile,
            "holdout": self._holdfile,
            "super": self._superfile
        }
        cfile = smap[kind]

        if not path.isfile(cfile):
            from matdb.utility import dbcat
            from fnmatch import fnmatch
            cfiles = []
            for seq in self.dbs:
                #We need to split to get training data. If the split has already
                #been done as part of a different training run, then it won't be
                #done a second time.
                msg.info("Compiling database {} for {}.".format(seq.name, self.fqn))
                seq.split()
                if seq.name in self.cust_splits:
                    splt = self.cust_splits[seq.name]
                else:
                    splt = self.split

                #We grab a list of all the files that match the particular split
                #pattern. Then we apply any filters to individual atoms objects
                #within each of the databases.
                if splt == '*':
                    nfiles = [f(seq) for f in fmap.values()]
                else:
                    nfiles = [fmap[kind](seq, splt)]

                filtered = self._filter_dbs(seq.name, nfiles)
                cfiles.extend(filtered)

            #First, save the configurations to a single file.
            dbcat(cfiles, cfile)

        if asatoms:
            return atoms.AtomsList(cfile)
    
    def validate(self, configs=None, energy=True, force=True, virial=True):
        """Validates the calculator in this training object against the `holdout.xyz`
        configurations in the source databases.

        Args:
            configs (atoms.AtomsList): list of configurations to validate
              against. If not provided, built-in holdout set will be used.
            energy (bool): when True, validate the energies of each
              configuration.
            forces (bool): when True, validate the force *components* of each
              configuration.
            virial (bool): when True, validate the virial *components* of each
              configuration.         

        Returns:
            dict: with keys ["*_dft", "*_pot"] where the values are
            reference-calculated and potential-calculated predictions for the
            energies, force components or virial components.
        """
        from tqdm import tqdm
        for a in tqdm(self.validation if configs is None else configs):
            a.set_cutoff(self.calculator.cutoff())
            a.calc_connect()
            self.calculator.calc(a, energy=energy, force=force, virial=virial)

        result = {}
        if energy:
            result["e_ref"] = np.array(al.dft_energy)
            result["e_pot"] = np.array(al.energy)
        if forces:
            result["f_ref"] = np.array(al.dft_force.flatten())
            result["f_pot"] = np.array(al.force.flatten())
        if virial:
            result["v_ref"] = np.array(al.dft_virial.flatten())
            result["v_pot"] = np.array(al.virial.flatten())

        return result
    
    def execute(self, dryrun=False):
        """Submits the job script for the currently configured potential training.

        Args:
            dryrun (bool): when True, simulate the submission without actually
              submitting.

        Returns:
            bool: True if the submission generated a job id (considered
            successful).
        """
        if not path.isfile(self._jobfile):
            return False

        if not path.isfile(self._trainfile):
            from matdb.msg import std
            std("train.xyz missing in {}; can't execute.".format(self.root))
            return False
        
        # We must have what we need to execute. Compile the command and submit.
        from matdb.utility import execute
        cargs = ["sbatch", self._jobfile]
        if dryrun:
            from matdb.msg import okay
            okay("Executed {} in {}".format(' '.join(cargs), self.root))
            return True
        else:
            xres = execute(cargs, self.root)

        if len(xres["output"]) > 0 and "Submitted" in xres["output"][0]:
            from matdb.msg import okay
            okay("{}: {}".format(self.root, xres["output"][0].strip()))
            return True
        else:
            return False
            
    def jobfile(self):
        """Creates the job file for training the potential. This includes creating
        folders, the training XYZ file and setting up any other dependencies
        required by the trainer.

        .. note:: This method also runs :meth:`command`, which configures the
          directory that the executable will run in. You are responsible for
          sub-classing that method correctly.
        """
        # We use the global execution parameters and then any updates
        # locally. We need to add the execution directory (including prefix) and
        # the number of jobs in the array.
        settings = self.execution.copy()
        settings["execution_path"] = self.root
        settings["exec_path"] = self.command()

        #Add the settings for virtual env and plot function so that we can
        #generate the plots in parallel.
        if self.controller.plotting is not None:
            settings.update(self.controller.plotting)
            settings["fqn"] = self.fqn
            settings["matdbyml"] = self.controller.db.config
        if self.controller.db.venv is not None:
            settings["venv"] = self.controller.db.venv
        
        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template(settings["template"])
        with open(self._jobfile, 'w') as f:
            f.write(template.render(**settings))
