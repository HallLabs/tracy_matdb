"""Implements classes and methods for performing a GAP fit over a
database defined in the YAML specification file.
"""
from os import path, rename, remove, mkdir
import numpy as np
from glob import glob
from tqdm import tqdm
from itertools import product, combinations
import tarfile

from jinja2 import Environment, PackageLoader
from phenum.makeStr import _make_structures

from matdb import msg
from matdb.utility import cat, chdir, _get_reporoot
from matdb.fitting.basic import Trainer
from matdb.exceptions import MlpError, LogicError
from matdb.atoms import Atoms
from matdb.database.active import Active
from matdb.database import Database
from matdb.io import atoms_to_cfg

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

        self.use_unrelaxed = mtpargs["use_unrelaxed"] if "use_unrelaxed" in mtpargs else False

        self.selection_limit = mtpargs["selection-limit"] if "selection-limit" in mtpargs else 300
        self.species = controller.db.species
        self.crystals_to_relax = mtpargs["crystals_to_relax"] if "crystals_to_relax" in mtpargs else ["sc", "fcc", "bcc", "hcp", "prototypes"]

        self.mtp_file = "pot.mtp"
        if path.isfile(path.join(self.root,"status.txt")):
            with open(path.join(self.root,"status.txt"),"r") as f:
                for line in f:
                    old_status = line.strip().split()
                self.iter_status = old_status[0]
        else:
            self.iter_status = None

        # we need access to the active learning set
        db_root = self.controller.db.root
        steps = [{"type":"active.Active"}]
        dbargs = {"root":db_root,
                  "parent":Database("active", db_root, self.controller.db, steps, {},
                                    self.ran_seed),
                  "calculator":self.controller.db.calculator}
        self.active = Active(**dbargs)
        self._trainfile = path.join(self.root, "train.cfg")
        
        #Configure the fitting directory for this particular potential.
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
            relax_args["fit_setting"] = "FALSE"
        
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
        if iteration == 1:
            for db in self.dbs:
                pbar = tqdm(total=len(db.fitting_configs))
                for atm in db.config_atoms:
                    self._create_train_cfg(atm, path.join(self.root, "train.cfg"))
                    pbar.update(1)
        else:
            if self.active.last_iteration is None or len(self.active.last_iteration) < 1:
                if path.isfile(self.active.iter_file):
                    self.active._load_last_iter()
                else:
                    msg.err("File {0} containing most recently added "
                            "structures is missing.".format(self.active.iter_file))
            msg.info("Extracting from {0} folders".format(len(self.active.last_iteration)))
            self.active.extract()
            pbar = tqdm(total=len(self.active.last_iteration.values()))
            for atm in db.config_atoms:
                self._create_train_cfg(atm, path.join(self.root, "train.cfg"))
                pbar.update(1)

    def _create_train_cfg(self, atm, target):
        """Creates a 'train.cfg' file for the calculation stored at the target
        directory.

        Args:
            atm (matdb.atoms.Atoms): an atoms object to write to the cfg file.
            target (str): the path to the desierd "train.cfg" file.            
        """
        temp_cfg = path.join(self.root, "temp.cfg")
        atoms_to_cfg(atm, temp_cfg)

        if not path.isfile(temp_cfg):
            raise IOError("Failed to create cfg file for atmso object stored "
                          "at: {0}".format(atm.calc.folder))
        
        if path.isfile(target):
            cat([temp.cfg, target], path.join(self.root, "temp2.cfg"))
            rename(path.join(self.root, "temp2.cfg"), target)
        else:
            rename(path.join(self.root, "temp.cfg"), target)
            
    def _make_pot_initial(self):

        """Creates the initial 'pot.mtp' file.
        """

        target = path.join(self.root, "pot.mtp")
        
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template("pot.mtp")

        with open(target,'w') as f:
            f.write(template.render(n_species=str(len(self.species))))
    
    def _make_relax_ini(self):
        """Creates the 'relax.ini' file for relaxing the structures.
        """

        target = path.join(self.root, "relax.ini")
        
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template("relax.ini")

        with open(target,'w') as f:
            f.write(template.render(**self.relax))

    def _make_to_relax_cfg(self, testing=False):
        """Creates the list of files to relax to check the mtp against.

        Args:
            testing (bool): True if unit tests are being run.
        """

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
        prot_map = {0: "uniqueUnaries", 1: "uniqueBinaries", 2: "uniqueTernaries"}
        for crystal in self.crystals_to_relax:
            if crystal == "prototypes":

                # The prototypes are saved into the file prototypes.tar.gz, if
                # this is the first time prototypes has been run we need to unpack it.
                template_root = path.join(_get_reporoot(), "matdb", "templates")
                if not path.isdir(path.join(template_root, "uniqueUnaries")):
                    with chdir(template_root):
                        tarf = "prototypes.tar.gz"
                        tar = tarfile.open(tarf, "r:gz")
                        tar.extractall()
                        tar.close()

                for size in range(len(self.species)):
                    dry_count = 0
                    cand_path = path.join(template_root, prot_map(size))
                    structures = glob("{0}/*".format(cand_path))
                    perms = [list(i) for i in permutations(self.species, r=size+1)]
                    for fpath, perm in product(structures, perms):
                        type_map = {}
                        for i, s in enumerate(self.species):
                            if s in perm:
                                type_map[perm.index(s)] = i
                        self._prot_to_cfg(fpath, perm, target, type_map)
                        
                    #For unit testing
                    if dry_run:
                        if dry_count == 10:
                            break
                        dry_count += 1
                        
                    
            else:
                if dry_run: 
                    infile = path.join(_get_reporoot(),"matdb","templates",
                                       "test_enum.out_{0}".format(crystal))
                else: #pragma: no cover
                    infile = path.join(_get_reporoot(),"matdb","templates",
                                       "struct_enum.out_{0}_{1}".format(len(self.species),crystal))
                args["input"] = infile
                _make_structures(args)
                    
        msg.info("to-relax.cfg file completed.")                   
        

    def _prot_to_cfg(self, source, species, relax_file, type_map):
        """Corrects the POSCAR so that it has the correct lattice parameter
        and title string for the system to be read into ASE.
        
        Args:
            source (str): the path to the prototype POSCAR.
            species (list): a list of the species to be used for this POSCAR.
            relax_file (str): the full path to the to-relax.cfg file being
              written.
            type_map (dict): the type mapping to apply to the cfg file.
        """

        f_lines = []
        lat_vecs = []
        # read in the original POSCAR
        with open(source, "r") as f:
            for i, line in enumerate(f):
                f_lines.append(line)
                if i in [2, 3, 4]:
                    lat_vecs.append([float(j) for j in line.strip().split()])
                if i == 5:
                    concs = [int(j) for j in line.strip().split()]

        # fix the title.
        f_lines[0] = "{0} : {1}".format(" ".join(species), f_lines[0])

        # fix the lattice parameter
        lat_param, ttl = get_lattice_parameter(species, concs, lat_vecs, sum(concs), " ")
        f_lines[1] = "{} \n".format(lat_param)

        target = path.join(self.root, "PROT")
        with open(target, "w+") as f:
            for line in f_lines:
                f.write(line)

        atm = Atoms(target)
        
        atoms_to_cfg(atm, path.join(self.root, "prot.cfg"))
        cat([relax_file, path.join(self.root, "prot.cfg")])
        config_id = "{0}_{1}".format("".join(self.species), source.split("/")[-1])
        remove(path.join(self.root, "prot.cfg"), config_id = config_id, type_map = type_map)
        remove(target)
                
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
                remove(path.join(self.root,"selected.cfg"))
            if path.isfile(path.join(self.root,"relaxed.cfg")):
                remove(path.join(self.root,"relaxed.cfg"))        
        
            if not path.isfile(path.join(self.root,"relax.ini")):
                self._make_relax_ini()
            if not path.isfile(path.join(self.root,"pot.mtp")):
                self._make_pot_initial()
            template = "mpirun -n {} mlp train pot.mtp train.cfg > training.txt".format(self.ncores)
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("relax {0}".format(iter_count))

        if self.iter_status == "relax":
            # If pot has been trained
            rename(path.join(self.root,"Trained.mtp_"), path.join(self.root,"pot.mtp"))
            
            # if the unrelaxed.cfg file exists we need to move it to
            # replace the existing 'to-relax.cfg' otherwise we need to
            # create the 'to-relax.cfg' file.
            with chdir(self.root):
                os.system("mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg")

            if not path.isfile(path.join(self.root, "state.mvs")):
                raise MlpError("mlp failed to produce the 'state.mvs` file with command "
                               "'mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg'")
            
            if path.isfile(path.join(self.root,"unrelaxed.cfg")) and self.use_unrelaxed:
                rename(path.join(self.root,"unrelaxed.cfg"),path.join(self.root,"to-relax.cfg"))
            elif not path.isfile(path.join(self.root,"to-relax.cfg")):
                self._make_to_relax_cfg()

            # command to relax structures
            template = "mpirun -n {} mlp relax relax.ini --cfg-filename=to-relax.cfg".format(self.ncores)
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("select {0}".format(iter_count))

        # if relaxation is done
        if self.iter_status == "select":
            cat(glob(path.join(self.root,"selected.cfg_*")), path.join(self.root,"selected.cfg"))

            # command to select next training set.
            with chdir(self.root):
                os.system("mlp select-add pot.mtp train.cfg selected.cfg diff.cfg --selection-limit={}".format(self.selection_limit))

                # Now add the selected atoms to the Active database.
                os.system("mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar")
                os.system("cp selected.cfg selected.cfg_iter_{}".format(iter_count))
                if path.isfile("relaxed.cfg"):
                    os.system("cp relaxed.cfg relaxed.cfg_iter_{}".format(iter_count))

                new_configs = []
                new_POSCARS = glob("POSCAR*")
                # Here We need to re-write the POSCARs so that they
                # have the correct species.
                for POSCAR in new_POSCARS:
                    # We need to put the correct species into the
                    # title of the POSCAR so that ASE can read it
                    pos_file = []
                    with open(POSCAR, 'r') as f:
                        for line in f:
                            pos_file.append(line)
                    pos_file[0] = "{0} : {1}".format(" ".join(self.species), pos_file[0])
                    with open(POSCAR, 'w+') as f:
                        for line in pos_file:
                            f.write(line)
                            
                    new_configs.append(Atoms(POSCAR,format="vasp"))
                    remove(POSCAR)

                self.active.add_configs(new_configs, iter_count)
                self.active.setup()
                if len(self.active.configs) != self.active.nconfigs: #pragma: no cover
                    raise LogicError("The active database failed to setup the calculations "
                                     "for iteration {0}".format(iter_count))
                self.active.execute()
            
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("done {0}".format(iter_count))
            template = ''

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
