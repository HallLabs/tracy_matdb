"""Implements classes and methods for performing a GAP fit over a
database defined in the YAML specification file.
"""
from os import path, rename, remove, mkdir
import numpy as np
from glob import glob
from tqdm import tqdm
from itertools import product, combinations
import tarfile
import json

from jinja2 import Environment, PackageLoader
from phenum.phenumStr import _make_structures

from matdb import msg
from matdb.utility import cat, chdir, _get_reporoot, execute
from matdb.fitting.basic import Trainer
from matdb.exceptions import MlpError, LogicError
from matdb.atoms import Atoms
from matdb.database.active import Active
from matdb.database import Database
from matdb.io import atoms_to_cfg

def create_to_relax(setup_args):
    """Creates the to-relax.cfg file from the passed in args dictionary.
    
    Args:
        setup_args (dict): A dictionary containing the arguments needed to 
          construct the potenital.
    """

    args = setup_args["args"]
    species = setup_args["species"]
    crystals = setup_args["crystals"]
    min_atoms = setup_args["min_atoms"]
    max_atoms = setup_args["max_atoms"]
    root = setup_args["root"]

    for crystal in crystals:
        prot_map = {0: "uniqueUnaries", 1: "uniqueBinaries", 2: "uniqueTernaries"}
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

            for size in range(len(species)):
                cand_path = path.join(template_root, prot_map(size))
                structures = glob("{0}/*".format(cand_path))
                perms = [list(i) for i in permutations(species, r=size+1)]
                for fpath, perm in product(structures, perms):
                    type_map = {}
                    for i, s in enumerate(species):
                        if s in perm:
                            type_map[perm.index(s)] = i
                    _prot_to_cfg(fpath, perm, target, type_map, root, min_atoms, max_atoms)
                    
        else:
            # eventually we'll want to replace this with an actual
            # enumeration over the options the user specifies.
            infile = path.join(_get_reporoot(),"matdb","templates",
                               "struct_enum.out_{0}_{1}".format(len(species),crystal))
            if min_atoms != 1 or max_atoms is not None:
                with open(infile, "r") as f:
                    min_num = None
                    max_num = None
                    past_start = False
                    for line in f:
                        if line.split()[0] == "start":
                            past_start = True
                            continue
                        if past_start:
                            data = line.strip().split()
                            lab = data[-2]
                            if len(lab) == min_atoms and min_num is None:
                                min_num = int(data[0])
                            if len(lab) > max_atoms and max_num is None:
                                max_num = int(data[0])-1
                                break
                args["structures"] = [min_num, max_num]

            args["input"] = infile
            _make_structures(args)

def _prot_to_cfg(source, species, relax_file, type_map, root, min_atoms, max_atoms):
    """Corrects the POSCAR so that it has the correct lattice parameter
    and title string for the system to be read into ASE.
    
    Args:
        source (str): the path to the prototype POSCAR.
        species (list): a list of the species to be used for this POSCAR.
        relax_file (str): the full path to the to-relax.cfg file being
          written.
        type_map (dict): the type mapping to apply to the cfg file.
        root (str): the root directory where the cfg needs to be stored.
        min_atoms (int): the smallest number of atoms wanted.
        max_atoms (int): the largest number of atoms wanted.
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

    # If this prototype has more or less atoms in the cell than the
    # user wants then
    if (sum(concs) < min_atoms) or (max_atoms is not None and sum(concs) > max_atoms):
        return
    # fix the title.
    f_lines[0] = "{0} : {1}".format(" ".join(species), f_lines[0])

    # fix the lattice parameter
    lat_param, ttl = get_lattice_parameter(species, concs, lat_vecs, sum(concs), " ")
    f_lines[1] = "{} \n".format(lat_param)

    target = path.join(root, "PROT")
    with open(target, "w+") as f:
        for line in f_lines:
            f.write(line)

    atm = Atoms(target)
        
    atoms_to_cfg(atm, path.join(root, "prot.cfg"))
    cat([relax_file, path.join(root, "prot.cfg")])
    config_id = "{0}_{1}".format("".join(species), source.split("/")[-1])
    remove(path.join(root, "prot.cfg"), config_id = config_id, type_map = type_map)
    remove(target)

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

        self._set_attributes(self, mtpargs)

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
        db_root = self.controller.db.root
        steps = [{"type":"active.Active"}]
        # here self.controller is a TController istance and
        # self.controller.db is the database Controller instance.
        dbargs = {"root":db_root,
                  "parent":Database("active", db_root, self.controller.db, steps, {},
                                    self.ran_seed),
                  "calculator":self.controller.db.calculator}
        self.active = Active(**dbargs)
        self._trainfile = path.join(self.root, "train.cfg")
        
        #Configure the fitting directory for this particular potential.
        if not path.isdir(self.root):
            mkdir(self.root)

    def _set_attributes(self, mtpargs):
        """Sets the attributes of the mtp object from the input dictionary.

        Args:
            mtpargs (dict): the input dictionary of arguments.
        """

        if "relax_ini" in mtpargs and mtpargs["relax_ini"] is not None:
            self._set_relax_ini(mtpargs["relax_ini"])
        else:
            self._set_relax_ini({})

        if "relax" in mtpargs and mtpargs["relax"] is not None:
            self.relax_args = mtpargs["relax"]
        else:
            self.relax_args = {}

        self.ran_seed = mtpargs["ran_seed"] if "ran_seed" in mtpargs else 0
            
        if "train" in mtpargs and mtpargs["train"] is not None:
            self.train_args = mtpargs["train"]
        else:
            self.train_args = {}

        if "smallest_relax_cell" in mtpargs:
            self.relax_min_atoms = mtpargs["smallest_relax_cell"]
        else:
            self.relax_min_atoms = 1

        if "largest_relax_cell" in mtpargs:
            self.relax_max_atoms = mtpargs["larges_relax_cell"]
        else:
            self.relax_max_atoms = None            

        if "calc-grade" in mtpargs and mtpargs["calc-grade"] is not None:
            self.grade_args = mtpargs["calc-grade"]
        else:
            self.grade_args = {}

        if "select" in mtpargs and mtpargs["select"] is not None:
            self.select_args = mtpargs["select"]
        else:
            self.select_args = {}

        if "to-relax" in mtpargs and mtpargs["to-relax"] is not None:
            self.to_relax_args = mtpargs["to-relax"]
        else:
            self.to_relax_args = {}

        if "use_mpi" in mtpargs and mtpargs["use_mpi"] is not None:
            self.use_mpi = mtpargs["use_mpi"]
        else:
            self.use_mpi = True

        self.use_unrelaxed = mtpargs["use_unrelaxed"] if "use_unrelaxed" in mtpargs else False
        
        self.crystals_to_relax = mtpargs["crystals_to_relax"] if "crystals_to_relax" in mtpargs else ["sc", "fcc", "bcc", "hcp", "prototypes"]
        

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

        self.relax_ini = relax_args

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
            pbar = tqdm(total=len(self.active.last_iteration))
            for atm in db.last_config_atoms:
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
            rename(temp_cfg, target)

        if path.isfile(temp_cfg):
            remove(temp_cfg)
        if path.isfile(path.join(self.root, "temp2.cfg")):
            remove(path.join(self.root, "temp2.cfg"))
            
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
            f.write(template.render(**self.relax_ini))

    def _setup_to_relax_cfg(self, testing=False):
        """Creates the list of files to relax to check the mtp against.

        Args:
            testing (bool): True if unit tests are being run.
        """

        target = path.join(self.root, "to-relax.cfg")

        if path.isfile(target):
            remove(target)
            
        setup_args = {}
        
        if len(self.species) >4:
            msg.err("The MTP relaxation isn't setup to create a to-relax.cfg "
                    "file for systems with more than 4 elements.")
            return
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
        setup_args["phenum_args"] = args
        setup_args["crystals"] = self.crystals_to_relax
        setup_args["species"] = self.species
        setup_args["min_atoms"] = self.relax_min_atoms
        setup_args["max_atoms"] = self.relax_max_atoms
        setup_args["root"] = self.root

        with open(path.join(self.root, "to_relax.json"), "w+"):
            json.dump(setup_args)

    def _train_template(self):
        """Creates the train command template.
        """
        # mlp train potential.mtp train_set.cfg [options]:
        #   trains potential.mtp on the training set from train_set.cfg
        #   Options include:
        #     --energy-weight=<double>: weight of energies in the fitting. Default=1
        #     --force-weight=<double>: weight of forces in the fitting. Default=0.01
        #     --stress-weight=<double>: weight of stresses in the fitting. Default=0.001
        #     --scale-by-force=<double>: Default=0. If >0 then configurations near equilibrium
        #                                (with roughtly force < <double>) get more weight. 
        #     --valid-cfgs=<string>: filename with configuration to validate
        #     --max-iter=<int>: maximal number of iterations. Default=1000
        #     --curr-pot-name=<string>: filename for potential on current iteration.
        #     --trained-pot-name=<string>: filename for trained potential. Default=Trained.mtp_
        #     --bfgs-conv-tol=<double>: stopping criterion for optimization. Default=1e-8
        #     --weighting=<string>: how to weight configuration wtih different sizes
        #         relative to each other. Default=vibrations. Other=molecules, structures.
        #     --init-params=<string>: how to initialize parameters if a potential was not
        #         pre-fitted. Default is random. Other is same - this is when interaction
        #         of all species is the same (more accurate fit, but longer optimization)
        #     --skip-preinit: skip the 75 iterations done when params are not given

        if self.use_mpi:
            template = ("mpirun -n {} mlp train pot.mtp "
                        "train.cfg".format(self.ncores))
        else:
            template = "train pot.mtp train.cfg"
            
        for k, v in self.train_args.items:
            if k == "curr-pot-name" or k == "trained-pot-name":
                msg.warn("Renaming of the potential file is not enabled.")
                continue
            if k == "valid-cfgs":
                msg.warn("Validating configurations is not enabled.")
                continue

            template = template + " --{0}={1}".format(k,v)

        return template + " > training.txt"

    def _calc_grade_template(self):
        """Creates the template for the calc-grade command.
        """
        # mlp calc-grade pot.mtp train.cfg in.cfg out.cfg:
        # actively selects from train.cfg, generates state.mvs file from train.cfg, and
        # calculates maxvol grades of configurations located in in.cfg
        # and writes them to out.cfg
        #   Options:
        #   --init-threshold=<num>: set the initial threshold to 1+num, default=1e-5
        #   --select-threshold=<num>: set the select threshold to num, default=1.1
        #   --swap-threshold=<num>: set the swap threshold to num, default=1.0000001
        #   --energy-weight=<num>: set the weight for energy equation, default=1
        #   --force-weight=<num>: set the weight for force equations, default=0
        #   --stress-weight=<num>: set the weight for stress equations, default=0
        #   --nbh-weight=<num>: set the weight for site energy equations, default=0
        #   --mvs-filename =<filename>: name of mvs file        
        template = "mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg"
        for k, v in self.grade_args:
            if k == "mvs-filename":
                msg.warn("Renaming the mvs state file is not enabled.")
            template = template + " --{0}={1}".format(k,v)

        return template
    
    def _relax_template(self):
        """Creates the template for the relax command.
        """
        # mlp relax settings-file [options]:
        # settings file should contain settings for relaxation and for mlip regime.
        #  Options can be given in any order. Options include:
        # --pressure=<num>: external pressure (in GPa)
        # --iteration_limit=<num>: maximum number of iterations
        # --min-dist=<num>: terminate relaxation if atoms come closer than <num>
        # --force-tolerance=<num>: relaxes until forces are less than <num>(eV/Angstr.)
        #       Zero <num> disables atom relaxation (fixes atom fractional coordinates)
        # --stress-tolerance=<num>: relaxes until stresses are smaller than <num>(GPa)
        #       Zero <num> disables lattice relaxation
        # --max-step=<num>: Maximal allowed displacement of atoms and lattice vectors
        #       (in Angstroms)
        # --min-step=<num>: Minimal displacement of atoms and lattice vectors (Angstr.)
        #       If all actual displacements are smaller then the relaxation stops.
        # --bfgs-wolfe_c1
        # --bfgs-wolfe_c2
        # --cfg-filename=<str>: Read initial configurations from <str>
        # --save-relaxed=<str>: Save the relaxed configurations to <str>
        # --save-unrelaxed=<str>: If relaxation failed, save the configuration to <str>
        # --log=<str>: Write relaxation log to <str>

        if self.use_mpi:
            template = ("mpirun -n {0} mlp relax relax.ini "
                        "--cfg-filename=to-relax.cfg "
                        "--save-relaxed={1} --log=relax_{2}"
                        "--save-unrelaxed={3}".format(self.ncores,
                                                      "relaxed.cfg",
                                                      "log.txt", "candidate.cfg"))
        else:
            template = ("mlp relax relax.ini "
                        "--cfg-filename=to-relax.cfg "
                        "--save-relaxed={1} --log=relax_{2}"
                        "--save-unrelaxed={3}".format(self.ncores,
                                                      "relaxed.cfg",
                                                      "log.txt", "candidate.cfg"))

        for k, v in self.relax_args.items():
            if k in ["log", "save-unrelaxed", "save-relaxed", "cfg-filename"]:
                msg.warn("Changing the {0} file name is not supported.".format(k))
                continue
            if k in ["bfgs-wolfe_c1", "bfgs-wolfe_c2"]:
                template = template + " --{0}".format(k)
            else:                         
                template = template + " --{0}={1}".format(k,v)

        return template

    def _select_template(self):
        """Creates the select command template.
        """
        # mlp select-add pot.mtp train.cfg new.cfg diff.cfg:
        # actively selects configurations from new.cfg and save those
        # that need to be added to train.cfg to diff.cfg
        #   Options:
        #   --init-threshold=<num>: set the initial threshold to num, default=1e-5
        #   --select-threshold=<num>: set the select threshold to num, default=1.1
        #   --swap-threshold=<num>: set the swap threshold to num, default=1.0000001
        #   --energy-weight=<num>: set the weight for energy equation, default=1
        #   --force-weight=<num>: set the weight for force equations, default=0
        #   --stress-weight=<num>: set the weight for stress equations, default=0
        #   --nbh-weight=<num>: set the weight for site energy equations, default=0
        #   --mvs-filename=<filename>: name of mvs file
        #   --selected-filename=<filename>: file with selected configurations
        #   --selection-limit=<num>: swap limit for multiple selection, default=0 (disabled)
        #   --weighting=<string>: way of weighting the functional for better fitting of
        # properties. Default=vibrations. Others=molecules, structures.
        
        template = "mlp select-add pot.mtp train.cfg candidate.cfg new_training.cfg"

        for k, v in self.select_args.items():
            if k in ["mvs-filename", "selected-filename"]:
                msg.warn("Changing the {0} file name is not enabled.")
                continue
            template = template + " --{0}={1}".format(k,v)

        return template
    
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
            template = self._train_template
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("relax_setup {0}".format(iter_count))

        if self.iter_status == "relax_setup":
            # If pot has been trained
            rename(path.join(self.root,"Trained.mtp_"), path.join(self.root,"pot.mtp"))

            # Calculate the grad of the training configurations.
            calc_grade = self._calc_grade_template()
            execute(calc_grade.split(), self.root)
            remove(path.join(self.root, "temp1.cfg"))
            if not path.isfile(path.join(self.root, "state.mvs")):
                raise MlpError("mlp failed to produce the 'state.mvs` file with command "
                               "'mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg'")
            
            # if the unrelaxed.cfg file exists we need to move it to
            # replace the existing 'to-relax.cfg' otherwise we need to
            # create the 'to-relax.cfg' file.            
            if path.isfile(path.join(self.root,"unrelaxed.cfg")) and self.use_unrelaxed:
                rename(path.join(self.root,"unrelaxed.cfg"),path.join(self.root,"to-relax.cfg"))
                self.iter_status = "relax"
            elif not path.isfile(path.join(self.root,"to-relax.cfg")):
                self._setup_to_relax_cfg()
                template = "matdb_mtp_to_relax.py"
                with open(path.join(self.root,"status.txt"),"w+") as f:
                    f.write("relax {0}".format(iter_count))
            else:
                self.iter_status = "relax"

        if self.iter_status == "relax":
            # command to relax structures
            template = self._relax_template()
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("select {0}".format(iter_count))

        # if relaxation is done
        if self.iter_status == "select":
            cat(glob(path.join(self.root,"candidate.cfg_*")), path.join(self.root,"candidate.cfg"))
            # command to select next training set.
            select_template = self._select_template()
            execute(select_template.split(), self.root)

            with chdir(self.root):
                # Now add the selected atoms to the Active database.
                # input filename sould always be entered before output filename
                execute(["mlp", "convert-cfg", "new_training.cfg", "POSCAR",
                         "--output-format=vasp-poscar"], self.root)
                rename("new_training.cfg", "new_training.cfg_iter_{}".format(iter_count))
                if path.isfile(self.relaxed_file):
                    rename("relaxed.cfg", "relaxed.cfg_iter_{}".format(iter_count))

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
        
    def _update_split_params(self):
        """Updates the parameter set for the trainer to include the splits.
        """
        if self.split is None:
            self.params["split"] = 1
        else:
            if len(self.dbs) > 0:
                _splitavg = []
                for db in self.dbs:
                    if db in self.cust_splits:
                        splt = self.cust_splits[db]
                        _splitavg.append(1 if splt == '*' else splt)
                    else:
                        _splitavg.append(db.splits[self.split])
                self.params["split"] = sum(_splitavg)/len(self.dbs)
