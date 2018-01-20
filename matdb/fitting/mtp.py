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
        self.ncores = execution["ntasks"]
        self.root = root 

        if mtpargs["relax"] is not None:
            self._set_relax_ini(mtpargs["relax"])
        else:
            self._set_relax_ini({})

        self.species = mtpargs["species"]        
        
        self.mtp_file = "{}.xml".format(self.name)
        
        #Configure the fitting directory for this particular potential.
        from os import mkdir
        if not path.isdir(self.root):
            mkdir(self.root)

    def _set_relax_ini(self,relaxargs):
        """Sets the arguments for the relax.ini file.

        Args:
            relaxargs (dict): A dictionary containing keys for the weights.
        """
        relax_args = {}

        if "calc-efs" in relaxargs:
            relax_args["calc_efs"] = relaxargs["calc-efs"]
        else:
            relax_args["calc_efs"] = "TRUE"

        if "active-learn" in relaxargs:
            relax_args["active_learn"] = relaxargs["active-learn"]
        else:
            relax_args["active_learn"] = "TRUE"

        if "fit" in relaxargs:
            relax_args["fit_setting"] = relaxargs["fit"]
        else:
            relaxargs["fit_setting"] = "FALSE"
        
        if "site-weight" in relaxargs:
            relax_args["site_weight"] = relaxargs["site-weight"]
        else:
            relax_args["site_weight"] = "0.0"

        if "energy-weight" in relaxargs:
            relax_args["energy_weight"] = relaxargs["energy-weight"]
        else:
            relax_args["energy_weight"] = "1.0"

        if "force-weight" in relaxargs:
            relax_args["force_weight"] = relaxargs["force-weight"]
        else:
            relax_args["force_weight"] = "0.001"

        if "stress-weight" in relaxargs:
            relax_args["stress_weight"] = relaxargs["stress-weight"]
        else:
            relax_args["stress_weight"] = "0.0001"

        if "threshold" in relaxargs:
            relax_args["extrap_threshold"] = relaxargs["threshold"]
        else:
            relax_args["extrap_threshold"] = "2.0"

        if "threshold-break" in relaxargs:
            relax_args["threshold_break"] = relaxargs["threshold-break"]
        else:
            relax_args["threshold_break"] = "10.0"

        self.relax = relax_args
        
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

        target = path.join(self.root, "pot.mtp")
        
        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template("pot.mtp")

        with open(target,'w') as f:
            f.write(template.render("n_species"=len(self.species)))
    
    def _make_relax_ini(self):
        """Creates the 'relax.ini' file for relaxing the structures.
        """

        target = path.join(self.root, "relax.ini")
        
        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template("relax.ini")

        with open(target,'w') as f:
            f.write(template.render(**self.relax))

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

        if not path.isfile(path.join(self.root,"relax.ini")):
            self._make_relax_ini()
        if not path.isfile(path.join(self.root,"pot.mtp")):
            self._make_pot_initial()
        
        self._make_train_cfg(self.configs)

        self._make_to_relax_cfg()

        if not path.isfile(path.join(self.root,"status.txt")):
            status = "train"
            iter_count = 1
        else:
            from os import remove
            with open(path.join(self.root,"status.txt"),"r") as f:
                for line in f:
                    old_status = line.strip().split()

            if old_status[0] == "select":
                status = "train"
            else:
                status = old_status[0]
            iter_count = int(old_status[1]) + 1
                
                    
            
        # IF at start
        #command to train the potential
        if status == "train":
            template = "mpirun -n {} mlp train pot.mtp train.cfg > training.txt"
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("relax {1}".format(iter_count))

        # If pot has been trained
        rename("Trained_mtp_","pot.mtp")

        # command to relax structures
        template = "mpirun -n {} mlp relax relax.ini --cfg-filename=to-relax.cfg"

        # if relaxation is done
        cat(glob("selected.cfg_*"),"selected.cfg")

        # command to select next training set.
        template = "mlp select-add pot.mtp traic.cfg selected.cfg diff.cfg -selection-limit={}"

        # if active set has been selected add it to the `Active` group.
        
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
