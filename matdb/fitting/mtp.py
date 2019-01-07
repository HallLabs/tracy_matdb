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
from itertools import permutations

from jinja2 import Environment, PackageLoader
from phenum.phenumStr import _make_structures
from phenum.element_data import get_lattice_parameter

from matdb import msg
from matdb.utility import cat, chdir, _get_reporoot, execute, touch
from matdb.fitting.basic import Trainer
from matdb.exceptions import MlpError, LogicError
from matdb.atoms import Atoms
from matdb.database.active import Active
from matdb.database import Database
from matdb.io import atoms_to_cfg

def _is_float(var):
    """Determines if the input is a float.
    """

    if isinstance(var, float):
        return True
    else:
        try:
            float(var)
            return True
        except:
            return False

def create_to_relax(setup_args):
    """Creates the to-relax.cfg file from the passed in args dictionary.
    
    Args:
        setup_args (dict): A dictionary containing the arguments needed to 
          construct the potenital.
    """

    args = setup_args["phenum_args"]
    species = setup_args["species"]
    crystals = setup_args["crystals"]
    min_atoms = setup_args["min_atoms"]
    max_atoms = setup_args["max_atoms"]
    root = setup_args["root"]

    touch(args["outfile"])
    
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
                cand_path = path.join(template_root, prot_map[size])
                structures = glob("{0}/*".format(cand_path))
                perms = [list(i) for i in permutations(species, r=size+1)]
                for fpath, perm in product(structures, perms):
                    type_map = {}
                    for i, s in enumerate(species):
                        if s in perm:
                            type_map[perm.index(s)] = i
                    _prot_to_cfg(fpath, perm, args["outfile"], type_map,
                                 root, min_atoms, max_atoms)
                    
        else:
            # eventually we'll want to replace this with an actual
            # enumeration over the options the user specifies.
            infile = path.join(_get_reporoot(),"matdb","templates",
                               "struct_enum.out_{0}_{1}".format(len(species),crystal))
            expected_min_atoms = 1
            if crystal == "hcp":
                min_atoms = 2
                expected_min_atoms = 2
            if min_atoms != expected_min_atoms or max_atoms is not None:
                with open(infile, "r") as f:
                    min_num = None
                    max_num = None
                    last_num = 1
                    past_start = False
                    for line in f:
                        if line.split()[0] == "start":
                            past_start = True
                            continue
                        if past_start:
                            data = line.strip().split()
                            lab = data[-2]
                            last_num = int(data[0])
                            if len(lab) == min_atoms and min_num is None:
                                min_num = int(data[0])
                            if len(lab) > max_atoms and max_num is None:
                                max_num = int(data[0])
                                break
                # In case in the input file, there is no configuration that has the number of atoms
                # bigger than max_atoms, the max_num will be set to the number in the last line.  
                if max_num is None:
                    max_num = last_num
                args["structures"] = range(min_num, max_num)

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

    atm = Atoms(target, format="vasp")
        
    atoms_to_cfg(atm, path.join(root, "prot.cfg"))
    cat([relax_file, path.join(root, "prot.cfg")], path.join(root, "temp.cfg"))
    rename(path.join(root, "temp.cfg"), relax_file)
    config_id = "{0}_{1}".format("".join(species), source.split("/")[-1])
    remove(path.join(root, "prot.cfg"))
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
        self._root = root
        self._set_root()
        
        if path.isfile(path.join(self.root,"status.txt")):
            with open(path.join(self.root,"status.txt"),"r") as f:
                for line in f:
                    old_status = line.strip().split()
                self.iter_status = old_status[0]
                self.iter_count = int(old_status[1])
                self.cell_iter = int(old_status[2])
        else:
            self.iter_status = None
            self.iter_count = None
            self.cell_iter = 0
            
        self._set_local_attributes(mtpargs)
        self.species = controller.db.species


        self.mtp_file = "pot.mtp"

        # sets the resource usage for this run from the maximum
        # supplied in the yml file.
        self._set_resources()
        # we need access to the active learning set
        db_root = self.controller.db.root
        steps = [{"type":"active.Active"}]
        # here self.controller is a TController istance and
        # self.controller.db is the database Controller instance.
        dbargs = {"root":db_root,
                  "parent":Database("active", db_root, self.controller.db, steps, {},
                                    self.ran_seed),
                  "calculator": self.controller.db.calculator}
        self.active = Active(**dbargs)
        fix_static = getattr(self.active.calc, "set_static", None)
        if callable(fix_static):
            self.active.calcargs = fix_static(self.active.calcargs)
        self._trainfile = path.join(self.root, "train.cfg")

    def _set_resources(self):
        """determines the resource usage depending on which `mtp` step we're
        on.
        """
        if self.iter_status in ("select", "relax_setup", "add"):
            self.ncores = 1
            self.execution["ntasks"] = 1
            self.execution["mem_per_cpu"] = self.execution["total_mem"]
            
        elif self.iter_status in (None, "train", "relax"):
            self.ncores = self.execution["ntasks"]
            self.execution["mem_per_cpu"] = self._convert_mem()

    def _convert_mem(self):
        """Converts memory to the correct values including units."""

        init_mem, final_unit = self.execution["total_mem"][:-2], self.execution["total_mem"][-2:]
        final_mem = float(init_mem)/self.execution["ntasks"]
        while final_mem < 1:
            final_mem = final_mem*1000
            if final_unit.lower() == "gb":
                final_unit = "MB"
            elif final_uint.lower() == "mb":
                final_unit = "KB"
            else:
                msg.error("Unrecognized memory size {}".format(final_unit))

        final_mem = int(np.round(final_mem))

        return "{0}{1}".format(final_mem, final_unit)
        
    def _set_root(self):

        """Sets the root directory.
        """
        if "mtp" != self._root.split('/')[-1]:
            self.root = path.join(self._root,"mtp")
        else:
            self.root = self._root
        #Configure the fitting directory for this particular potential.
        if not path.isdir(self.root): #pragma: no cover
            # This never gets run. Kept as a fail safe.
            mkdir(self.root)        
        
    def _set_local_attributes(self, mtpargs):
        """Sets the attributes of the mtp object from the input dictionary.

        Args:
            mtpargs (dict): the input dictionary of arguments.
        """

        if "iteration_threshold" in mtpargs and mtpargs["iteration_threshold"] is not None:
            self.iter_threshold = mtpargs["iteration_threshold"]
        else:
            self.iter_threshold = 50

        if "next_cell_threshold" in mtpargs and mtpargs["next_cell_threshold"] is not None:
            self.next_cell_threshold = mtpargs["next_cell_threshold"]
        else:
            self.next_cell_threshold = 0
            
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
            if isinstance(mtpargs["largest_relax_cell"],list):
                self.cell_sizes = mtpargs["largest_relax_cell"]
                if len(self.cell_sizes) > self.cell_iter:
                    self.relax_max_atoms = self.cell_sizes[self.cell_iter]
                else:
                    msg.err("The MTP process has finished for the cell sizes {0}"
                            "specified in the YML file.".format(self.cell_sizes))
            else:
                self.cell_sizes = [mtpargs["largest_relax_cell"]]
                self.relax_max_atoms = mtpargs["largest_relax_cell"]
        else:
            self.relax_max_atoms = None
            self.cell_sizes = None

        if (self.relax_max_atoms is not None and
            (self.relax_max_atoms < self.relax_min_atoms or
             self.relax_max_atoms <0)) or self.relax_min_atoms < 0:
            raise ValueError("The max and min number of atoms to be relaxed must be "
                             "larger than 0 and the max must be larger than the min.")

        if "calc-grade" in mtpargs and mtpargs["calc-grade"] is not None:
            self.grade_args = mtpargs["calc-grade"]
        else:
            self.grade_args = {}

        if "select" in mtpargs and mtpargs["select"] is not None:
            self.select_args = mtpargs["select"]
        else:
            self.select_args = {}

        if "use_mpi" in mtpargs and mtpargs["use_mpi"] is not None:
            self.use_mpi = False if "false" in mtpargs["use_mpi"].lower() else True
        else:
            self.use_mpi = True

        if "use_unrelaxed" in mtpargs:
            self.use_unrelaxed = True if "true" in mtpargs["use_unrelaxed"].lower()  else False
        else:
            self.use_unrelaxed = False            
        
        self.crystals_to_relax = mtpargs["crystals_to_relax"] if "crystals_to_relax" in mtpargs else ["sc", "fcc", "bcc", "hcp", "prototypes"]
        

    def _set_relax_ini(self,relaxargs):
        """Sets the arguments for the relax.ini file.
        Args:
            relaxargs (dict): A dictionary containing keys for the weights.
        """
        relax_args = {}

        if "calc-efs" in relaxargs and relaxargs["calc-efs"] in ["TRUE", "FALSE"]:
            relax_args["calc_efs"] = relaxargs["calc-efs"]
        else:
            relax_args["calc_efs"] = "TRUE"

        if "efs-ignore" in relaxargs and relaxargs["efs-ignore"] in ["TRUE", "FALSE"]:
            relax_args["efs_ignore"] = relaxargs["efs-ignore"]
        else:
            relax_args["efs_ignore"] = "FALSE"

        if "active-learn" in relaxargs and relaxargs["active-learn"] in ["TRUE", "FALSE"]:
            relax_args["active_learn"] = relaxargs["active-learn"]
        else:
            relax_args["active_learn"] = "TRUE"

        if "fit" in relaxargs and relaxargs["fit"] in ["TRUE", "FALSE"]:
            relax_args["fit_setting"] = relaxargs["fit"]
        else:
            relax_args["fit_setting"] = "FALSE"
        
        if "site-weight" in relaxargs and _is_float(relaxargs["site-weight"]):
            relax_args["site_weight"] = relaxargs["site-weight"]
        else:
            relax_args["site_weight"] = "0.0"

        if "energy-weight" in relaxargs and _is_float(relaxargs["energy-weight"]):
            relax_args["energy_weight"] = relaxargs["energy-weight"]
        else:
            relax_args["energy_weight"] = "1.0"

        if "force-weight" in relaxargs  and _is_float(relaxargs["force-weight"]):
            relax_args["force_weight"] = relaxargs["force-weight"]
        else:
            relax_args["force_weight"] = "0.001"

        if "stress-weight" in relaxargs and _is_float(relaxargs["stress-weight"]):
            relax_args["stress_weight"] = relaxargs["stress-weight"]
        else:
            relax_args["stress_weight"] = "0.0001"

        if "threshold" in relaxargs and _is_float(relaxargs["threshold"]):
            relax_args["extrap_threshold"] = relaxargs["threshold"]
        else:
            relax_args["extrap_threshold"] = "2.0"

        if "threshold-break" in relaxargs and _is_float(relaxargs["threshold-break"]):
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
        from matdb.database.legacy import LegacyDatabase
        if iteration == 1:
            for db in self.dbs:
                if not isinstance(db, LegacyDatabase):
                    for step in db.steps.values():
                        pbar = tqdm(total=len(step.rset))
                        for atm in step.rset:
                            self._create_train_cfg(atm, path.join(self.root, "train.cfg"))
                            pbar.update(1)
                else:
                    pbar = tqdm(total=len(db.rset))
                    for atm in db.rset:
                        self._create_train_cfg(atm, path.join(self.root, "train.cfg"))
                        pbar.update(1)
                            
        else:
            if self.active.last_iteration is None or len(self.active.last_iteration) < 1:
                if path.isfile(self.active.iter_file):
                    self.active._load_last_iter()
                else:
                    raise IOError("File {0} containing most recently added "
                                  "structures is missing.".format(self.active.iter_file))
            msg.info("Extracting from {0} folders".format(len(self.active.last_iteration)))
            self.active.extract()
            pbar = tqdm(total=len(self.active.last_iteration))
            for atm in self.active.last_config_atoms.values():
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

        if not path.isfile(temp_cfg): #pragma: no cover
            raise IOError("Failed to create cfg file for atoms object stored "
                          "at: {0}".format(atm.calc.folder))
        
        if path.isfile(target):
            cat([temp_cfg, target], path.join(self.root, "temp2.cfg"))
            rename(path.join(self.root, "temp2.cfg"), target)
        else:
            rename(temp_cfg, target)

        if path.isfile(temp_cfg):
            remove(temp_cfg)
        if path.isfile(path.join(self.root, "temp2.cfg")):#pragma: no cover
            # This line should only get run if something goes terribly
            # wrong up above. I'm keeping it here purely as a fail
            # safe.
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

    def _setup_to_relax_cfg(self):
        """Creates the list of files to relax to check the mtp against.
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

        with open(path.join(self.root, "to_relax.json"), "w+") as f:
            json.dump(setup_args, f)

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
            template = ("mpirun --allow-run-as-root -n {} mlp train pot.mtp "
                        "train.cfg".format(self.ncores))
        else:
            template = "mlp train pot.mtp train.cfg"
            
        for k, v in self.train_args.items():
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
        for k, v in self.grade_args.items():
            if k == "mvs-filename":
                msg.warn("Renaming the mvs state file is not enabled.")
                continue
            template = template + " --{0}={1}".format(k,v)

        return template + " > training_calc_grade.txt"
    
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
            template = ("mpirun --allow-run-as-root -n {0} mlp relax relax.ini "
                        "--cfg-filename=to-relax.cfg "
                        "--save-relaxed={1} --log=relax_{2} "
                        "--save-unrelaxed={3}".format(self.ncores,
                                                      "relaxed.cfg",
                                                      "log.txt",
                                                      "unrelaxed.cfg"))
        else:
            template = ("mlp relax relax.ini "
                        "--cfg-filename=to-relax.cfg "
                        "--save-relaxed={1} --log=relax_{2} "
                        "--save-unrelaxed={3}".format(self.ncores,
                                                      "relaxed.cfg",
                                                      "log.txt",
                                                      "unrelaxed.cfg"))

        for k, v in self.relax_args.items():
            if k in ["log", "save-unrelaxed", "save-relaxed", "cfg-filename"]:
                msg.warn("Changing the {0} file name is not supported.".format(k))
                continue
            if k in ["bfgs-wolfe_c1", "bfgs-wolfe_c2"]:
                template = template + " --{0}".format(k)
            else:                         
                template = template + " --{0}={1}".format(k,v)

        return template + " > training_relax.txt"

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

        return template + " > training_select.txt"
    
    def command(self):
        """Returns the command that is needed to train the GAP
        potentials specified by this object.

        .. note:: This method also configures the directory that the command
          will run in so that it has the relevant files.
        """
        
        self._set_root()

        if not path.isfile(path.join(self.root,"status.txt")):
            self.iter_status = "train"
            self.iter_count = 1
            self.cell_iter = 0
        else:
            with open(path.join(self.root,"status.txt"),"r") as f:
                for line in f:
                    old_status = line.strip().split()

            if old_status[0] == "done":
                self.iter_status = "train"
                self.iter_count = int(old_status[1]) + 1
                    
                if int(old_status[3]) <= self.next_cell_threshold:
                    self.iter_count = 1
                    self.cell_iter = int(old_status[2]) + 1
                    if len(self.cell_sizes) < self.cell_iter:
                        self.relax_max_atoms = self.cell_sizes[self.cell_iter]
                    else:
                        msg.err("The MTP process has finished for the cell sizes "
                                "specified in the YML file.")
                    target1 = path.join(self.root,"unrelaxed.cfg")
                    target2 = path.join(self.root,"to-relax.cfg")
                    if path.isfile(target1):
                        remove(target1)
                    if path.isfile(target2):
                        remove(target2)                        

                if self.iter_count >= self.iter_threshold:
                    msg.err("The number of iterations for this size has exceeded "
                            "the allowed number of {0}. To increase this limit "
                            "use the iteration_threshold aption in the YML "
                            "file.".format(self.iter_threshold))
            else:
                self.iter_status = old_status[0]
                self.iter_count = int(old_status[1])
                self.cell_iter = int(old_status[2])
                
        #if we're at the start of a training iteration use the command to train the potential
        if self.iter_status == "train":
            self._make_train_cfg(self.iter_count)

            #remove the selected.cfg and rexaed.cfg from the last iteration if they exist
            if path.isfile(path.join(self.root,"new_training.cfg")):
                remove(path.join(self.root,"new_training.cfg"))
            if path.isfile(path.join(self.root,"relaxed.cfg")):
                remove(path.join(self.root,"relaxed.cfg"))        
        
            if not path.isfile(path.join(self.root,"relax.ini")):
                self._make_relax_ini()
            if not path.isfile(path.join(self.root,"pot.mtp")):
                self._make_pot_initial()
            template = self._train_template()
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("relax_setup {0} {1}".format(self.iter_count, self.cell_iter))

        if self.iter_status == "relax_setup":
            # If pot has been trained
            if path.isfile(path.join(self.root,"Trained.mtp_")):
                rename(path.join(self.root,"Trained.mtp_"), path.join(self.root,"pot.mtp"))

            # Calculate the grad of the training configurations.
            calc_grade = self._calc_grade_template()
            execute(calc_grade.split(), self.root)
            if path.isfile(path.join(self.root, "temp1.cfg")):
                remove(path.join(self.root, "temp1.cfg"))
            if not path.isfile(path.join(self.root, "state.mvs")): #pragma: no cover
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
                template = "matdb_mtp_to_relax.py > create-to-relax.txt"
                with open(path.join(self.root,"status.txt"),"w+") as f:
                    f.write("relax {0} {1}".format(self.iter_count, self.cell_iter))
            else:
                self.iter_status = "relax"

        if self.iter_status == "relax":
            # command to relax structures
            if path.isfile(path.join(self.root,"candidate.cfg")):
                remove(path.join(self.root,"candidate.cfg"))

            template = self._relax_template()
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("select {0} {1}".format(self.iter_count, self.cell_iter))

        # if relaxation is done
        if self.iter_status == "select":
            cand_files = glob(path.join(self.root,"candidate.cfg_*"))
            cat(cand_files, path.join(self.root,"candidate.cfg"))
            for cfile in cand_files:
                remove(cfile)

            # command to select next training set.
            template = self._select_template()
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("add {0} {1}".format(self.iter_count, self.cell_iter))

        if self.iter_status == "add":
            with chdir(self.root):
                # Now add the selected atoms to the Active database.
                # input filename sould always be entered before output filename
                execute(["mlp", "convert-cfg", "new_training.cfg", "POSCAR",
                         "--output-format=vasp-poscar"], self.root)
                if path.isfile("new_training.cfg"):
                    rename("new_training.cfg", "new_training.cfg_iter_{}".format(self.iter_count))
                if path.isfile("relaxed.cfg"):
                    rename("relaxed.cfg", "relaxed.cfg_iter_{}".format(self.iter_count))

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

            self.active.add_configs(new_configs, self.iter_count)
            self.active.setup()
            if len(self.active.configs) != self.active.nconfigs: #pragma: no cover
                raise LogicError("The active database failed to setup the calculations "
                                 "for iteration {0}".format(self.iter_count))
            self.active.execute()
            
            with open(path.join(self.root,"status.txt"),"w+") as f:
                f.write("done {0} {1} {2}".format(self.iter_count, self.cell_iter,
                                                  len(new_configs)))
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
        last_iter = self.active.last_iteration
        result = {
            "trained": self.ready(),
            "file": self.mtp_file,
            "jobfile": path.isfile(self._jobfile),
            "mtp step": self.iter_status,
            "last addition": len(last_iter) if last_iter is not None else 0
        }
        
        if printed:
            fqn = "{0}.{1}".format(self.parent.name, self.name)
            msg.info("{0} => Model ready: {1}".format(fqn, result["trained"]))
            x = "exists" if result["jobfile"] else "does not exist"
            msg.info("{0} => Next jobfile '{1}' {2}".format(fqn, self._jobfile, x))
            msg.info("{0} => Current MTP step {1} iteration {2}.".format(fqn, self.iter_status,
                                                                      self.iter_count))
            msg.info("{0} => {1} configurations added "
                     "to the training set.".format(fqn, result["last addition"]))
            msg.blank()
        else:
            return result
        
    def _update_split_params(self):
        """Updates the parameter set for the trainer to include the splits.
        """
        if self.split is None:
            self.params["split"] = 1
        else: #pragma: no cover (this feature is tested by the base
              #class. It also doesn't make sense to use splits within
              #the mtp context.
            super(MTP, self)._update_split_params()
