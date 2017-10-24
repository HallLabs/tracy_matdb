"""Implements a basic trainer that train material models.
"""
import abc
from os import path
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

    Args:
        name (str): name of this model; used as the output folder name.
        controller (matdb.fitting.controller.Controller): fitting controller
          whose data will be used to train the potentials.
        dbs (list): of `str` database patterns from the `db` that should be
          included in the training and validation.
        execution (dict): settings needed to configure the jobfile for running
          the fit.

    Attributes:
    """
    def __init__(self, controller=None, dbs=None, execution={}, split=None,
                 root=None):
        self.controller = controller
        self.execution = {} if execution is None else execution.copy()
        self.split = split
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
        for dbpat in self._dbs:
            self.dbs.extend(self.controller.db.find(dbpat))
            
        self._calculator = None
        """ase.Calculator: potential for the fitted file in this Trainer.
        """
        self._trainfile = path.join(self.root, "train.xyz")
        """str: path to the XYZ training file that will be passed to the training
        command.
        """
        self._jobfile = path.join(self.root, "jobfile.sh")
        """str: path to the jobfile for the current trainer.
        """

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
        """Returns a :class:`quippy.AtomsList` of configurations that can be
        used for potential validation.
        """
        result = quippy.AtomsList()
        for seq in self.dbs:
            #We need to split to get validation data. If the split has already
            #been done as part of a different training run, then it won't be
            #done a second time.
            seq.split(self.split)
            result.extend(quippy.AtomsList(seq.holdout_file))

        return result
    
    def validate(self, energy=True, force=True, virial=True):
        """Validates the calculator in this training object against the `holdout.xyz`
        configurations in the source databases.

        Args:
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
        for a in tqdm(self.validation):
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
        if path.isfile(self._jobfile):
            return

        #Setup the train.xyz file for the set of databases specified in the
        #source patterns.
        from matdb.utility import cat
        tfiles = []
        for seq in self.dbs:
            #We need to split to get training data. If the split has already
            #been done as part of a different training run, then it won't be
            #done a second time.
            seq.split(self.split)
            tfiles.append(seq.train_file)
        tfiles.extend(self.extras)
        cat(tfiles, self._trainfile)
        
        # We use the global execution parameters and then any updates
        # locally. We need to add the execution directory (including prefix) and
        # the number of jobs in the array.
        settings = self.execution.copy()
        settings["execution_path"] = self.root
        settings["exec_path"] = self.command()

        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template(settings["template"])
        with open(target, 'w') as f:
            f.write(template.render(**settings))
