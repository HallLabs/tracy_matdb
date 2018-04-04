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

        self.selection_limit = mtpargs["selection-limit"]
        self.species = mtpargs["species"]        
        
        self.mtp_file = "pot.mtp"
        if path.isfile(path.join(self.root,"status.txt")):
            with open(path.join(self.root,"status.txt"),"r") as f:
                for line in f:
                    old_status = line.strip().split()
                self.iter_status = old_status[0]
        else:
            self.stat_iter = None

        # we need access to the active learning set
        from matdb.database.active import Active
        self.active = Active()
        
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
        """Determines if the potential is ready for use.
        """

        pot_file = (path.isfile(path.join(self.root, self.mtp_file)))
        stat = (self.iter_state == "done")
        next_iter = (len(self.last_iteration) == 0)
        return all(pot_file,stat,next_iter)

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
                for config in db.configs.values():
                    self._create_train_cfg(path.join(self.root,"train.cfg"),config)
        else:
            for config in self.active.last_iter.values():
                self._create_train_cfg(path.join(self.root,"train.cfg"),config)

    def _create_train_cfg(self,target):
        """Creates a 'train.cfg' file for the calculation stored at the target
        directory.

        Args:
            target (str): the path to the directory in which a calculation 
                was performed.
        """
        from os import system, rename
        from matdb.utility import cat
        if path.isfile(path.join(target,"OUTCAR")):
            self._parse_POSCAR(target)
            system("mlp convert-cfg {0}/OUTCAR {1}/diff.cfg --input-format=vasp-outcar >> outcar.txt".format(target,self.root))
            if path.isfile(path.join(self.root,"diff.cfg")):
                if path.isfile(path.join(self.root,"train.cfg")):
                    cat([path.join(self.root,"train.cfg"),path.join(self.root,"diff.cfg")],
                        path.join(self.root,"temp.cfg"))
                    rename(path.join(self.root,"temp.cfg"),path.join(self.root,"train.cfg"))
                else:
                    rename(path.join(self.root,"diff.cfg"),path.join(self.root,"train.cfg"))
            else:
                msg.err("There was an error making the config file for folder "
                        "{}".format(path.join(self.root,target)))
        else:
            msg.err("The folder {} didn't run.".format(path.join(self.root,target)))

    def _parse_POSCAR(self,target):
        """Changes the POTCAR and CONTCAR to have the correct concentration
        string with the zeros intact.

        Args:
            target (str): the path to the directory in which a calculation 
                was performed.
        """
        from os import rename

        if not path.isfile(path.join(target,"CONTCAR")):
            msg.err("Calculations for {0} directory didn't finish.".format(target))
        else:
            rename(path.join(target,"CONTCAR"),path.join(target,"CONTCAR_ase"))
            rename(path.join(target,"POSCAR"),path.join(target,"POSCAR_ase"))

            # We need to grad the first line of the POSCAR to
            # determine the species present and the concentration
            # string of the POSCAR.
            with open(path.join(target,"POSCAR_ase"),"r") as f:
                specs = f.readline().strip().split()
                temp = f.readline()
                for i in range(3):
                    temp = f.readline()
                concs = f.readline().strip().split()
                if not RepresentsInt(concs[0]):
                    concs = f.readline().strip().split()
                if not RepresentsInt(concs[0]) or len(specs) != len(concs):
                    msg.err("Could not Parse concentration from POSCAR in  "
                            "{0}".format(target))

            if self.species != specs:
                new_concs = []
                j = 0
                for s in self.species:
                    # if s isn't the same as the species in specs
                    # then we've found one of the missing species
                    if s != specs[j]:
                        new_concs.append(0)
                    else:
                        new_concs.append(concs[j])
                        j += 1
            else:
                new_concs = concs

            with open(path.join(target,"POSCAR"),"w+") as new_f:
                with open(path.join(target,"POSCAR_ase"),"w+") as old_f:
                    for i, old_line in enumerate(old_f):
                        if i== 0:
                            new_f.write("{}\n".format(" ".join(self.species)))
                        elif i in [5,6]:
                            old_concs = line.strip().split()
                            if RepresentsInt(old_concs[0]):
                                new_f.write("  {}\n".format("   ".join(new_concs)))
                            else:
                                new_f.write(old_line)
                        else:
                            new_f.write(old_line)

            with open(path.join(target,"CONTCAR"),"w+") as new_f:
                with open(path.join(target,"CONTCAR_ase"),"w+") as old_f:
                    for i, old_line in enumerate(old_f):
                        if i== 0:
                            new_f.write("{}\n".format(" ".join(self.species)))
                        elif i in [5,6]:
                            old_concs = line.strip().split()
                            if RepresentsInt(old_concs[0]):
                                new_f.write("  {}\n".format("   ".join(new_concs)))
                            else:
                                new_f.write(old_line)
                        else:
                            new_f.write(old_line)

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
                    for edge in combinations(range(len(species)),size):
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
            if os.path.isfile(path.join(self.root,"selected.cfg")):
                from os import remove
                remove(path.join(self.root,"selected.cfg"))
            if os.path.isfile(path.join(self.root,"relaxed.cfg")):
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
            system("mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg")
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
            from os import system
            from glob import glob
            from quippy.atoms import Atoms
            system("mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar")
            system("cp selected.cfg selected.cfg_iter_{}".format(iter_count))
            if path.isfile("relaxed.cfg"):
                system("cp relaxed.cfg relaxed.cfg_iter_{}".format(iter_count))

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
