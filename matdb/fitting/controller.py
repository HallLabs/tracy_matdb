"""Since some potential fitting requires multi-step processes, we need to
implement controller objects that keep the various training steps ordered
correctly, and provide easy access to relevant data from previous steps. This
module supplies objects (similar in spirit to :mod:`matdb.database.controller`
for producting sequences of training objects that can be repeated across
parameter grids.
"""
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
        parent (TController): instance controlling multiple trainers.
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
                continue
            
            cls = getattr(module, clsname)

            #Make a copy of the original dictionary so that we don't mess up the
            #pointers; then add in the keyword arguments that are missing.
            cpspec = tspec.copy()
            del cpspec["type"]
            cpspec["root"] = self.root
            cpspec["parent"] = self

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
        
        if niterations is not None:
            from matdb.utility import obj_update
            
            for i, sequence in enumerate(niterations):
                nname = None
                isteps = copy(steps)
                for k, v in sequence.items():
                    if k == "suffix":
                        nname = self.name + v
                        continue
                    
                    obj_update(isteps, k, v, False)
                    
                if nname is None:
                    nname = self.name + "-{0:d}".format(i)

                iobj = TrainingSequence(nname, self, self.root, parent, isteps,
                                        **kwargs)
                self.sequences[nname] = iobj
        else:
            single = TrainingSequence(name, self, self.root, parent, steps,
                                      **kwargs)
            self.sequences[name] = single

class TController(object):
    """Controller for handling sequences of training steps.

    Args:
        db (matdb.database.controller.Controller): controller for the available
          data sets to use.
        root (str): path to the root directory for the system.
        fits (dict): key-value pairs defining the list of fitters to try out.
        e0 (list): of `float` values indicating the *isolated* atom energy in
          `eV`.
        kwargs (dict): additional key-value pairs to pass down to the individual
          trainer objects.
    """
    def __init__(self, db=None, root=None, fits=None, e0=None, **kwargs):
        self.db = db
        self.e0 = e0
        self.fits = OrderedDict()

        for fspec in fits:
            cpspec = fspec.copy()
            cpspec["root"] = self.root
            cpspec["controller"] = self
            seq = TSequenceRepeater(**cpspec)
            self.fits[fspec["name"]] = seq

    def find(self, pattern):
        """Finds a list of trainers that match the specified pattern. Trainer FQNs
        follow the pattern `trainer-suffix.step` where the suffix may be
        optional and the step names are defined by the `name` of each
        :class:`Trainer` sub-class. Wildcard `*` or any other
        :func:`~fnmatch.fnmatch` patterns may be used.
        """
        trainer, step = key.split('.')
        if '-' in trainer:
            repname, suffix = trainer.split('-')
        else:
            repname, suffix = trainer, None

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
            repname, suffix = trainer.split('-')
        else:
            repname, suffix = trainer, None

        result = None
        if repname in self.fits:
            rep = self.fits[repname]
            if trainer in rep.sequences:
                seq = rep.sequences[trainer]
                result = seq.steps.get(step)

        return result
