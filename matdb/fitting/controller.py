"""Since some potential fitting requires multi-step processes, we need to
implement controller objects that keep the various training steps ordered
correctly, and provide easy access to relevant data from previous steps. This
module supplies objects (similar in spirit to :mod:`matdb.database.__init__`
for producting sequences of training objects that can be repeated across
parameter grids.
"""
from collections import OrderedDict
from os import path
import six
from matdb import msg

class TrainingSequence(object):
    """Represents a sequence of training steps (each sub-classing
    :class:`~matdb.training.basic.Trainer`) where each subsequent step depends
    on the results of the previous step.

    Args:
        name (str): name of the training sequence so it can be referenced by FQN
          notation.
        repeater (TSequenceRepeater): repeater for a set of training sequences
          that all share the same root training and validation databases.
        root (str): root directory in which all other training sequences for
          the configurations in the same specification will be stored.
        controller (TController): instance controlling multiple trainers.
        steps (list): of `dict` describing the kinds of potential training steps
          to setup.
        kwargs (dict): key-value pairs to pass to the individual trainer
          objects (only used if they don't override the values).
    """
    def __init__(self, name, repeater, root, controller, steps, **kwargs):
        self.name = name
        self.root = path.join(root, name)
        self.repeater = repeater
        self.controller = controller
        
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        from importlib import import_module
        self._settings = steps
        """dict: with keys and values describing the kinds of training steps to setup.
        """

        from collections import OrderedDict
        self.steps = OrderedDict()
        for tspec in steps:
            if isinstance(tspec, six.string_types):
                #This is a reference to an existing database instance that was
                #defined previously.
                instance = self.controller[tspec]
                self.steps[instance.name] = instance
                continue
                
            modname, clsname = tspec["type"].split('.')
            fqdn = "matdb.fitting.{}".format(modname)
            module = import_module(fqdn)
            if not hasattr(module, clsname):# pragma: no cover
                #We haven't implemented this database type yet, just skip the
                #initialization for now.
                msg.warn("Cannot find trainer of type {}.".format(tspec["type"]))
                continue
            
            cls = getattr(module, clsname)

            #Make a copy of the original dictionary so that we don't mess up the
            #pointers; then add in the keyword arguments that are missing.
            cpspec = tspec.copy()
            del cpspec["type"]
            cpspec["root"] = self.root
            cpspec["parent"] = self
            cpspec["controller"] = self.controller

            #Add in the default values passed in from the parent instances, but
            #only update them if they weren't specified.
            for k, v in kwargs.items():
                if k not in cpspec:
                    cpspec[k] = v

            instance = cls(**cpspec)
            self.steps[instance.name] = instance

    @property
    def isteps(self):
        """Returns a generator over steps in this sequence. The generator yields
        the next step *only* if the previous one is already finished (i.e., the
        `ready()` method returns True.
        """
        previous = None
        for tname, trainer in self.steps.items():
            if previous is None or previous[1].ready():
                previous = (tname, trainer)
                yield previous
            else:
                raise StopIteration()

    def jobfiles(self):
        """Runs the jobfile creation methods of the steps in this sequence, in
        the correct order and only as each step completes.
        """
        for fitname, fit in self.isteps:
            if fit.ready():
                continue        

            fit.jobfile()
            if path.isfile(path.join(fit._jobfile)):
                fqn = "{}.{}".format(self.name, fitname)
                msg.okay("Completed jobfile creation for {}.".format(fqn))

    def status(self, printed=True):
        """Determines status information for each of the training steps within
        this sequence.

        Args:
            printed (bool): when True, print the status to screen; otherwise,
              return a dictionary with the relevant quantities.
        """
        for fitname, fit in self.isteps:
            fit.status(printed)

    def execute(self, dryrun=False):
        """Submits the job script for the each training step in the sequence.

        Args:
            dryrun (bool): when True, simulate the submission without actually
              submitting.
        """
        for fitname, fit in self.isteps:
            fit.execute(dryrun)
            
class TSequenceRepeater(object):
    """Represents a group of training sequences that share the same underlying
    database of configurations for training and validation, but that vary only
    by standard parameter grids.
    """
    def __init__(self, name=None, root=None, controller=None, steps=None,
                 niterations=None, **kwargs):
        self.name = name
        self.sequences = OrderedDict()
        self.root = path.join(root, name)
        
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)
        
        if niterations is not None:
            from matdb.utility import obj_update, pgrid
            from copy import copy
            
            for i, sequence in enumerate(niterations):
                suffix = sequence.get("suffix")
                grid, keys = pgrid(sequence, ["suffix"])
                
                for ival, vals in enumerate(grid):
                    isteps = copy(steps)
                    for k, oval in zip(keys, vals):
                        obj_update(isteps, k, oval, False)

                    nsuffix = "" if len(vals) == 1 else "-{0:d}".format(ival)
                    if suffix is not None:
                        nn = self.name + "-{0}{1}".format(suffix, nsuffix)
                    else:
                        nn = self.name + "-{0:d}{1}".format(i, nsuffix)

                    iobj = TrainingSequence(nn, self, self.root,
                                            controller, isteps, **kwargs)
                    self.sequences[nn] = iobj
        else:
            single = TrainingSequence(name, self, self.root, controller, steps,
                                      **kwargs)
            self.sequences[name] = single

class TController(object):
    """Controller for handling sequences of training steps.

    Args:
        db (matdb.database.Controller): controller for the available data sets
          to use.
        root (str): path to the root directory for the system.
        fits (dict): key-value pairs defining the list of fitters to try out.
        e0 (list): of `float` values indicating the *isolated* atom energy in
          `eV`.
        plotting (dict): keyword arguments for plotting to generate for each
          trainer after it is finished fitting. This uses
          :func:`~matdb.plotting.potentials.generate`. 
        kwargs (dict): additional key-value pairs to pass down to the individual
          trainer objects.
    """
    def __init__(self, db=None, root=None, fits=None, e0=None, plotting=None,
                 **kwargs):
        from matdb.utility import dict_update
        self.db = db
        self.e0 = e0
        self.root = root
        self.plotting = plotting
        self.fits = OrderedDict()

        for fspec in fits:
            cpspec = fspec.copy()
            cpspec["root"] = self.root
            cpspec["controller"] = self

            dict_update(cpspec, kwargs)
            seq = TSequenceRepeater(**cpspec)
            self.fits[fspec["name"]] = seq

    def ifiltered(self, tfilter=None, sfilter=None):
        """Returns a *filtered* generator over sequences of trainers in the :attr:`fits`
        collection of this object.

        Args:
            tfilter (list): of `str` patterns to match against *fit* names.
            tfilter (list): of `str` patterns to match against *step* names.
        """
        from fnmatch import fnmatch
        for name, fit in self.fits.items():
            if tfilter is None or any(fnmatch(name, t) for t in tfilter):
                for fitn, seq in fit.sequences.items():
                    if sfilter is None or any(fnmatch(fitn,s) for s in sfilter):
                        yield (fitn, seq)

    def steps(self):
        """Compiles a list of all steps in this set of trainers.
        """
        result = []
        for name, fits in self.fits.items():
            for repeater in fits.values():
                for parent, seq in repeater.sequences.items():
                    for step in seq.steps:
                        result.append("{0}.{1}".format(parent, step))

        return sorted(result)        
    
    def sequences(self):
        """Compiles a list of all sequences in this set of trainers.
        """
        result = []
        for name, fits in self.fits.items():
            for repeater in fits.values():
                result.extend(repeater.sequences.keys())

        return sorted(result) 

    def validate(self, datafile=None, tfilter=None, sfilter=None, energy=True,
                 force=True, virial=True):
        """Validates all potentials/fitters in this controller against the
        specified data file. If not specified, then the built-in hold out sets
        are used.

        Args:
            datafile (str): path to the data file to read the atoms list from.
            tfilter (list): of `str` patterns to match against *fit* names.
            tfilter (list): of `str` patterns to match against *step* names.        
            energy (bool): when True, validate the energies of each
              configuration.
            forces (bool): when True, validate the force *components* of each
              configuration.
            virial (bool): when True, validate the virial *components* of each
              configuration.
        """
        if datafile is not None:
            configs = quippy.AtomsList(datafile)
        else:
            configs = None

        for trainer in self.ifiltered(tfilter, sfilter):
            trainer.validate(configs, energy, force, virial)
    
    def jobfiles(self, tfilter=None, sfilter=None):
        """Creates the jobfiles for every trainer in the system.

        Args:
            tfilter (list): of `str` patterns to match against *fit* names.
            tfilter (list): of `str` patterns to match against *step* names.
        """
        for fitname, fit in self.ifiltered(tfilter, sfilter):
            fit.jobfiles()

    def execute(self, tfilter=None, sfilter=None, dryrun=False):
        """Executes the jobfiles for every trainer in the system.

        Args:
            tfilter (list): of `str` patterns to match against *fit* names.
            tfilter (list): of `str` patterns to match against *step* names.
        """
        for fitname, fit in self.ifiltered(tfilter, sfilter):
            fit.execute(dryrun)

    def status(self, tfilter=None, sfilter=None, printed=True):
        """Prints status information for training steps in the controller.

        Args:
            tfilter (list): of `str` patterns to match against *fit* names.
            tfilter (list): of `str` patterns to match against *step* names.
            printed (bool): when True, print the status to screen; otherwise,
              return a dictionary with the relevant quantities.
        """
        for fitname, fit in self.ifiltered(tfilter, sfilter):
            fit.status(printed)
            
    def find(self, pattern):
        """Finds a list of trainers that match the specified pattern. Trainer FQNs
        follow the pattern `trainer-suffix.step` where the suffix may be
        optional and the step names are defined by the `name` of each
        :class:`Trainer` sub-class. Wildcard `*` or any other
        :func:`~fnmatch.fnmatch` patterns may be used.
        """
        if pattern.count('.') == 2:
            repname, trainer, step = pattern.split('.')
        else:
            trainer, step = pattern.split('.')
            if '-' in trainer:
                repname = trainer.split('-')[0]
            else:
                repname = trainer

        from fnmatch import fnmatch
        result = []
        reps = [irep for krep, irep in self.fits.items()
                if fnmatch(krep, repname)]
        for rep in reps:
            seqs = [iseq for kseq, iseq in rep.sequences.items()
                    if fnmatch(kseq, trainer)]
            for seq in seqs:
                result.extend([istep for kstep, istep in seq.steps.items()
                               if fnmatch(kstep, step)])

        return result
            
    def __getitem__(self, key):
        """Gets the trainer or sequence with the specified key. The keys follow
        the pattern `trainer-suffix.step` where the suffix may be optional and
        the step names are defined by the `name` of each :class:`Trainer`
        sub-class.
        """
        trainer, step = key.split('.')
        if '-' in trainer:
            repname = trainer.split('-')[0]
        else:
            repname = trainer

        result = None
        if repname in self.fits:
            rep = self.fits[repname]
            if trainer in rep.sequences:
                seq = rep.sequences[trainer]
                result = seq.steps.get(step)

        return result
