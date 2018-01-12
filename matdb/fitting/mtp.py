"""Implements classes and methods for performing a GAP fit over a
database defined in the YAML specification file.
"""
from os import path, rename
from matdb import msg
from collections import OrderedDict
import numpy as np
import quippy
from .basic import Trainer
from matdb.utility import cat
from glob import glob

def dict_to_str(settings, spacer=""):
    """Converts the specified dictionary of QUIP-compatible settings into a
    single string.

    Args:
        settings (dict): key-value pairs that are recognized settings by QUIP
          for descriptors, teach_sparse, eval, etc.
    """
    result = []
    for k, v in settings.items():
        if isinstance(v, (list, set, tuple)):
            value = '{' + ' '.join(map(str, v)) + '}'
        else:
            value = v

        result.append("{}={}".format(k, value))

    joiner = " \\\n  {}".format(spacer)
    return joiner.join(result)

class MTP(Trainer):
    """Implements a simple wrapper around the MTP training functionality for
    creating MTP potentials.

    Args:
        controller (matdb.fitting.controller.Controller): fitting controller
          provides access to previous fitting steps and training/validation data.
        split (str): name of the split specification to use for training.
        execution (dict): settings needed to configure the jobfile for running
          the GAP fit.
        dbs (list): of `str` patterns from the database that should be included
          in the training and validation.

    Attributes:
        name (str): name of the folder in which all the fitting for this trainer
          takes place; defaults to `{n}b`.
        root (str): root directory that this trainer operates in.
        params (dict): key-value pairs that are parameters for the model
          fitting.
    """
    def __init__(self, controller=None, dbs=None, execution=None,
                 split=None, root=None, parent=None, dbfilter=None, **mtpargs):
        self.name = "{}b".format(nb) if isinstance(nb, int) else "mtp"
        super(MTP, self).__init__(controller, dbs, execution, split, root,
                                  parent, dbfilter)
        self.controller = controller

        self.mtp_file = "{}.xml".format(self.name)
        
        #Configure the fitting directory for this particular potential.
        from os import mkdir
        if not path.isdir(self.root):
            mkdir(self.root)

    def get_calculator(self):
        """Returns an instance of :class:`ase.Calculator` using the latest
        fitted GAP potential in this trainer.
        """
        from quippy.potential import Potential
        if path.isfile(self.mtp_file):
            return Potential("IP MTP", param_filename=self.mtp_file)

    def ready(self):
        return path.isfile(path.join(self.root, self.mtp_file))

    def _make_train_cfg(self):
        """Creates the 'train.cfg' file needed to train the potential.
        """
        pass

    def _make_pot_initial(self):
        """Creates the initial 'pot.mtp' file.
        """
        pass
    
    def _make_relax_ini(self):
        """Creates the 'relax.ini' file for relaxing the structures.
        """
        pass

    def _make_to_relax_cfg(self):
        """Creates the list of files to relax to check the mtp against.
        """
        pass
        
    def command(self):
        """Returns the command that is needed to train the GAP
        potentials specified by this object.

        .. note:: This method also configures the directory that the command
          will run in so that it has the relevant files.
        """

        self._make_relax_ini()
        self._make_pot_initial()
        self._make_train_cfg(self.configs)

        self._make_to_relax_cfg()

        #command to train the potential
        template = "mpirun -n {} mlp train pot.mtp train.cfg > training.txt"

        rename("Trained_mtp_","pot.mtp")
        
        # command to relax structures
        template = "mpirun -n {} mlp relax relax.ini --cfg-filename=to-relax.cfg"

        cat(glob("selected.cfg_*"),"selected.cfg")

        # command to select next training set.
        template = "mlp select-add pot.mtp traic.cfg selected.cfg diff.cfg -selection-limit={}"
        
        return template.format(**fields)

    def status(self, printed=True):
        """Returns or prints the current status of the MTP training.

        Args:
            printed (bool): when True, print the status to screen; otherwise,
              return a dictionary with the relevant quantities.
        """
        # Our interest is in knowing which GAP model is the latest (if any) and
        # whether the job script has been created for the next one in the
        # sequence.
        result = {
            "trained": path.isfile(self.mtp_file),
            "file": self.mtp_file,
            "jobfile": path.isfile(self._jobfile)
        }
        
        if printed:
            fqn = "{}.{}".format(self.parent.name, self.name)
            msg.info("{} => Model ready: {}".format(fqn, result["trained"]))
            x = "exists" if result["jobfile"] else "does not exist"
            msg.info("{} => Next jobfile '{}' {}".format(fqn, self._jobfile, x))
            msg.blank()
        else:
            return result
