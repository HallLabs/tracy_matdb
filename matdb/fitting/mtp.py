"""Implements classes and methods for performing a GAP fit over a
database defined in the YAML specification file.
"""
from os import path, rename
from matdb import msg
from collections import OrderedDict
import numpy as np
from .basic import Trainer
from matdb.utility import cat
from glob import glob
import os
from tqdm import tqdm

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

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
        self.name = "mtp"
        super(MTP, self).__init__(controller, dbs, execution, split, root,
                                  parent, dbfilter)
        self.controller = controller
        self.ncores = execution["ntasks"]
        if "mtp" not in self.root:
            self.root = path.join(root,"mtp")
        else:
            self.root = root

        if "relax" in mtpargs and mtpargs["relax"] is not None:
            self._set_relax_ini(mtpargs["relax"])
        else:
            self._set_relax_ini({})

        self.selection_limit = mtpargs["selection-limit"]
        self.species = controller.db.species
        
        self.mtp_file = "pot.mtp"
        if path.isfile(path.join(self.root,"status.txt")):
            with open(path.join(self.root,"status.txt"),"r") as f:
                for line in f:
                    old_status = line.strip().split()
                self.iter_status = old_status[0]
        else:
            self.iter_status = None

        # we need access to the active learning set
        from matdb.database.active import Active
        from matdb.database import Database
        db_root = self.controller.db.root
        steps = [{"type":"active.Active"}]
        dbargs = {"root":db_root,
                  "parent":Database("active", db_root, self.controller.db, steps, {}),
                  "calculator":self.controller.db.calculator}
        self.active = Active(**dbargs)
        self._trainfile = path.join(self.root, "train.cfg")
        
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
        
    # def get_calculator(self):
    #     """Returns an instance of :class:`ase.Calculator` using the latest
    #     fitted GAP potential in this trainer.
    #     """
    #     from quippy.potential import Potential
    #     if path.isfile(self.mtp_file):
    #         return Potential("IP MTP", param_filename=self.mtp_file)

    def ready(self):
        """Determines if the potential is ready for use.
        """

        pot_file = (path.isfile(path.join(self.root, self.mtp_file)))
        stat = (self.iter_status == "done")
        next_iter = (len(self.active.last_iteration) == 0 if
                     self.active.last_iteration is not None else False)
        return all([pot_file,stat,next_iter])

    def _make_train_cfg(self,iteration):
        """Creates the 'train.cfg' file needed to train the potential from the
        databeses used.

        Args:
            iteration (int): the number of iterations of MTP has been 
                through.
        """
        from matdb.utility import cat
        if iteration == 1:
            for db in self.dbs:
                pbar = tqdm(total=len(db.fitting_configs))
                for config in db.fitting_configs:
                    self._create_train_cfg(config)
                    pbar.update(1)
        else:
            pbar = tqdm(total=len(self.active.last_iteration.values()))
            for config in self.active.last_iteration.values():
                self._create_train_cfg(path.join(self.root,"train.cfg"),config)
                pbar.update(1)

    def _create_train_cfg(self,target):
        """Creates a 'train.cfg' file for the calculation stored at the target
        directory.

        Args:
            target (str): the path to the directory in which a calculation 
                was performed.
        """
        from os import rename
        from matdb.utility import cat
        if path.isfile(path.join(target,"OUTCAR")):
            mapping = self._get_mapping(target)
            os.system("mlp convert-cfg {0}/OUTCAR {1}/diff.cfg --input-format=vasp-outcar >> outcar.txt".format(target,self.root))
            if path.isfile(path.join(self.root,"diff.cfg")):
                rename(path.join(self.root,"diff.cfg"),
                       path.join(self.root,"diff_orig.cfg"))
                with open(path.join(self.root,"diff_orig.cfg"),"r") as f_in:
                    with open(path.join(self.root,"diff.cfg"),"w+") as f_out:
                        for i, line_in in enumerate(f_in):
                            if i==2:
                                n_atoms = int(line_in.strip())
                                f_out.write(line_in)
                            elif i >= 8 and i < 8+n_atoms:
                                temp_line = line_in.strip().split()
                                temp_line[1] = mapping[int(temp_line[1])]
                                temp_line = "            {0}    {1}       {2}      {3}      {4}     {5}    {6}    {7}".format(*temp_line)
                                f_out.write(temp_line)
                            else:
                                f_out.write(line_in)
                if path.isfile(path.join(self.root,"train.cfg")):
                    cat([path.join(self.root,"train.cfg"), path.join(self.root,"diff.cfg")],
                        path.join(self.root,"temp.cfg"))
                    rename(path.join(self.root,"temp.cfg"), path.join(self.root,"train.cfg"))
                else:
                    rename(path.join(self.root,"diff.cfg"), path.join(self.root,"train.cfg"))
            else:
                msg.err("There was an error making the config file for folder "
                        "{}".format(path.join(self.root,target)))
        else:
            msg.err("The folder {} didn't run.".format(path.join(self.root,target)))

    def _get_mapping(self,target):
        """Finds the species mappings for the atomic numbers found in the
        trani.cfg file so that is will be correct.

        Args:
            target (str): the path to the directory in which a calculation 
                was performed.

        """
        from os import rename

        if not path.isfile(path.join(target,"POSCAR")):
            msg.err("Setup failed for {0}.".format(target))
        else:
            # We need to grab the first line of the POSCAR to
            # determine the species present.
            with open(path.join(target,"POSCAR"),"r") as f:
                specs = f.readline().strip().split()

            mapping = {}
            j = 0
            for i, s  in enumerate(specs):
                mapping[i] = self.species.index(s)

        return mapping

    def _make_pot_initial(self):

        """Creates the initial 'pot.mtp' file.
        """

        target = path.join(self.root, "pot.mtp")
        
        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template("pot.mtp")

        with open(target,'w') as f:
            f.write(template.render(n_species=str(len(self.species))))
    
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
        from matdb.utility import _get_reporoot
        from os import path, remove
        from phenum.makeStr import _make_structures
        from itertools import combinations

        target = path.join(self.root,"to-relax.cfg")

        if path.isfile(target):
            remove(target)
        
        args = {"config":"t",
                "species":self.species,
                "structures": "all",
                "debug": False,
                "example": False,
                "displace":0.0,
                "mink":True,
                "outfile":target,
                "verbose": None,
                "rattle":0.0,
                "mapping": None,
        }

        if len(self.species) >4:
            msg.err("The MTP relaxation isn't setup to create a to-relax.cfg "
                    "file for systems with more than 4 elements.")
        
        msg.info("Setting up to-relax.cfg file.")
        for crystal in ["bcc","fcc","sc","hcp"]:
            for size in range(2,len(self.species)+1):
                # if the size we're currently on is smaller than the
                # system in question then we need to loop over the
                # different species mappings possible to correctly form
                # all the edges/faces of the phase diagram.
                if size != len(self.species):
                    infile = path.join(_get_reporoot(),"matdb","templates",
                                       "struct_enum.out_{0}_{1}_sub".format(size,crystal))
                    args["input"] = infile
                    for edge in combinations(range(len(self.species)),size):
                        args["mapping"] = {i:j for i, j in enumerate(edge)}
                        _make_structures(args)
                else:
                    infile = path.join(_get_reporoot(),"matdb","templates",
                                       "struct_enum.out_{0}_{1}".format(size,crystal))
                    args["input"] = infile
                    _make_structures(args)
                    
        msg.info("to-relax.cfg file completed.")                   
        
    def command(self):
        """Returns the command that is needed to train the GAP
        potentials specified by this object.

        .. note:: This method also configures the directory that the command
          will run in so that it has the relevant files.
        """
        
        if not path.isfile(path.join(self.root,"status.txt")):
            self.iter_status = "train"
            iter_count = 1
        else:
            from os import remove
            with open(path.join(self.root,"status.txt"),"r") as f:
                for line in f:
                    old_status = line.strip().split()

            if old_status[0] == "done":
                self.iter_status = "train"
                iter_count = int(old_status[1]) + 1
            else:
                self.iter_status = old_status[0]
                iter_count = int(old_status[1])
                
        #if we're at the start of a training iteration use the command to train the potential
        if self.iter_status == "train":
            self._make_train_cfg(iter_count)

            #remove the selected.cfg and rexaed.cfg from the last iteration if they exist
            if path.isfile(path.join(self.root,"selected.cfg")):
                from os import remove
                remove(path.join(self.root,"selected.cfg"))
            if path.isfile(path.join(self.root,"relaxed.cfg")):
                from os import remove
                remove(path.join(self.root,"relaxed.cfg"))        
        
            if not path.isfile(path.join(self.root,"relax.ini")):
                self._make_relax_ini()
            if not path.isfile(path.join(self.root,"pot.mtp")):
                self._make_pot_initial()
            template = "mpirun -n {} mlp train pot.mtp train.cfg > training.txt".format(self.ncores)
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("relax {0}".format(iter_count))

        if self.iter_status == "relax":
            # if the unrelaxed.cfg file exists we need to move it to
            # replace the existing 'to-relax.cfg' otherwise we need to
            # create the 'to-relax.cfg' file.
            os.system("mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg")
            if path.isfile(path.join(self.root,"unrelaxed.cfg")):
                from os import rename
                rename(path.join(self.root,"unrelaxed.cfg"),path.join(self.root,"to-relax.cfg"))
            else:
                self._make_to_relax_cfg()

            # If pot has been trained
            rename("Trained_mtp_","pot.mtp")

            # command to relax structures
            template = "mpirun -n {} mlp relax relax.ini --cfg-filename=to-relax.cfg".format(self.ncores)
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("select {0}".format(iter_count))

        # if relaxation is done
        if self.iter_status == "select":
            cat(glob("selected.cfg_*"),"selected.cfg")

            # command to select next training set.
            template = "mlp select-add pot.mtp traic.cfg selected.cfg diff.cfg -selection-limit={}".format(self.selection_limit)

        # if selection is done
        if self.iter_status == "add":
            from glob import glob
            from quippy.atoms import Atoms
            os.system("mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar")
            os.system("cp selected.cfg selected.cfg_iter_{}".format(iter_count))
            if path.isfile("relaxed.cfg"):
                os.system("cp relaxed.cfg relaxed.cfg_iter_{}".format(iter_count))

            new_confgs = []
            new_POSCARS = glob("POSCAR*")
            for POSCAR in new_POSCARS:
                new_configs.append(Atoms(POSCAR,format="POSCAR"))

            self.active.add_configs(new_configs)
            self.active.setup()
            
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("done {0}".format(iter_count))
            
        return template

    def status(self, printed=True):
        """Returns or prints the current status of the MTP training.

        Args:
            printed (bool): when True, print the status to screen; otherwise,
              return a dictionary with the relevant quantities.
        """
        # Our interest is in knowing which MTP model is the latest (if any) and
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
