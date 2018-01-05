"""Abstract base class for creating and interacting with a database group of
configurations for machine learning materials.
"""
from os import path, mkdir
from matdb import msg
from .controller import ParameterGrid
import abc
import pickle
import datetime
from matdb import calculators
from contextlib import contextmanager
import ase.db
from uuid import uuid4
import numpy as np
import quippy
import six 

def atoms_to_json(atoms, folder):
    """Exports the specified atoms object, with its calculator's parameters, to
    a `atoms.json` file.
    Args:
        atoms (quippy.Atoms): configuration to write to JSON.
        folder (str): path to the folder to write the file in.
    """
    db = ase.db.connect(path.join(folder, "atoms.json"))
    #We need to extract out the additional parameters that a quippy.Atoms object
    #can have, but which are not supported by ASE.
    handled = ["id", "volume", "magmom", "age", "mass", "formula",
               "user", "charge", "unique_id", "fmax"]
    data = {}
    for param, value in atoms.params.items():
        if param not in handled:
            data[param] = value

    props = []
    for prop, value in atoms.properties.items():
        if prop not in ["pos", "species", "Z", "n_neighb", "map_shift"]:
            data[prop] = value
            props.append(prop)
    data["propnames"] = props

    db.write(atoms, data=data)

def atoms_from_json(folder):
    """Retrieves a :class:`quippy.Atoms` object from JSON in the specified
    folder.
    Args:
        folder (str): path to the folder that has the `atoms.json` file.
    """
    db = ase.db.connect(path.join(folder, "atoms.json"))
    row = db.get()
    atoms = quippy.Atoms()
    _atoms = row.toatoms()
    atoms.copy_from(_atoms)

    props = row.data.propnames
    for prop in props:
        atoms.properties[prop] = np.array(row.data[prop])

    for param in row.data:
        if param not in props:
            atoms.params[param] = row.data[param]

    calcargs = row.get('calculator_parameters', {})
    calculator = getattr(calculators, row.calculator.title())
    #We have to coerce all unicode strings into standard strings for
    #interoperability with the Fortran in QUIP.
    if row.calculator == "quip":
        caster = lambda a: str(a) if isinstance(a, unicode) else a
        calcargs["calcargs"] = list(map(caster, calcargs["calcargs"]))
    atoms.calc = calculator(atoms, folder, **calcargs)
            
    return atoms

class Group(object):
    """Represents a collection of material configurations (varying in
    structure and composition) from which a machine learning model can
    be created. Includes logic for generating the DFT directories that
    need to be run as well as extracting the relevant data from such
    computations.
    Args:
        execution (dict): key-value pairs of settings for the supercomputer job
          array batch file.
        root (str): full path to the root directory that this database will live
          in.
        parent (matdb.database.controller.Database): the database that this 
          group of calculations belong to.
        prefix (str): sub-sampled configurations will be stored using integer
          ids after this prefix; for example `S.1`, `S.2`, etc.
        nconfigs (int): number of displaced configurations to create.
        config_type (str): the type of configuration.
        calculator (dict): a dictionary containing the information for
          the calculator object.
       trainable (bool): True if the groups configs will be used for traning.
       pgrid (ParamaterGrid): The ParameterGrid for the database.
       seeds (list, str, quippy.Atoms): The location of the files that will be
          read into to make the atoms object or an atoms object.
    Attributes:
        atoms (quippy.atoms.Atoms): a single atomic configuration from
          which many others may be derived using MD, phonon
          displacements, etc.
        configs (dict): keys are integer identifiers for the particular
          configuration; values are paths (relative to the base atoms root
          directory) in which calculations are performed.
        root (str): full path to the root directory that this database will live
          in.
        database (matdb.database.controller.Database): parent database
          to which this group belongs.
        prefix (str): sub-sampled configurations will be stored using integer
          ids after this prefix; for example `S.1`, `S.2`, etc.
        nconfigs (int): number of displaced configurations to create.
    """
    seeded = False
    def __init__(self, root=None, parent=None, prefix='S', pgrid = None,
                 nconfigs=None, calculator=None, seeds=None, db_name=None,
                 config_type=None, parameters=None, execution={}, trainable=False):
        from collections import OrderedDict
        from quippy.atoms import Atoms
        from os import path

        self.index = {}
        self.root = root            
        self._read_index()
        self.parent = parent
        self.execution = execution
        self.params = None
        self.atoms = None
        
        self._trainable = trainable
        self.is_seed_listed = None
        if seeds is not None:
            self.is_seed_listed = isinstance(seeds, (list,six.string_types))

        self.seeds = None
        if self.is_seed_listed:
            self.seeds = OrderedDict()
            for atomspec in seeds:
                fmt, pattern = atomspec.split(':')
                for apath in self.parent.controller.relpaths(pattern):
                    self.seeds[path.basename(apath)] = Atoms(apath,format=fmt)
                    
        elif seeds is None and self.seeded:
            self.seeds = OrderedDict()
            for n_seeds, a in enumerate(self.prev.rset):
                seedname = "seed-{}".format(n_seeds)
                self.seeds[seedname] = a.copy()            
        
        self.sequence = OrderedDict()
        if self.seeds is not None:
            for seedname, at_seed in self.seeds.items():
                seed_root = path.join(root, seedname)                    
                if not path.isdir(seed_root):
                    mkdir(seed_root)                                        
                self.sequence[seedname] = Group(root=this_root,
                                                parent=parent,
                                                prefix=prefix,
                                                seeds=at_seed,
                                                db_name=db_name,
                                                calculator=calculator,
                                                pgrid=pgrid)
        else:
            if pgrid is not None and len(pgrid) >0:
                for params in pgrid:
                    this_root = path.join(root,pgrid.to_str(params))
                    if not path.isdir(this_root):
                        mkdir(this_root)
                    self.sequence[pgrid.to_str(params)] = Group(root=this_root,parent=parent,
                                                                prefix=prefix,
                                                                db_name=db_name,
                                                                calculator=calculator,
                                                                seeds = seeds,
                                                                parameters=pgrid[params])
            else:
                self.atoms = seeds
                self.params = parameters
                
        self.calc = None
        if calculator is not None:
            self.calc = getattr(calculators, calculator["name"])
        self.calcargs = calculator
        self.prefix = prefix
        self.nconfigs = nconfigs
        if self.params is not None and "nconfigs" in self.params:
            self.nconfigs = self.params["nconfigs"]
        self.config_type = config_type
        
        self._nsuccess = 0
        """int: number of configurations whose output files were successfully
        converted to XYZ format. Should be equal to :attr:`nconfigs` if the
        database is complete.
        """
        self._db_name = db_name
        
        #Try and load existing folders that match the prefix into the configs
        #list.
        from glob import glob
        from os import path
        from matdb.utility import chdir
        self.configs = {}
        self.config_atoms = {}
        self.uuid = uuid4()

        with chdir(self.root):
            for folder in glob("{}.*".format(prefix)):
                try:
                    cid = int(folder.split('.')[1])
                    self.configs[cid] = path.join(self.root, folder)
                    self.config_atoms[cid] = atoms_from_json(self.configs[cid])
                except:
                    #The folder name doesn't follow our convention.
                    pass

    @property
    def trainable(self):
        """Determines if the group configs should be used for training.
        """
        if self.dep is None:
            return True
        else:
            return self._trainable

    @abc.abstractproperty
    def fitting_configs(self):
        """Returns a :class:`quippy.AtomsList` for all configs in this group.
        """
        pass
    
    @abc.abstractproperty
    def rset(self):
        """Saves the rset for the group and all sequences of the group.
        """
        pass

    def save_pkl(self, obj, file_name, rel_path=None):
        """Saves the obj passed to the correct location on the path.
        
        Args:
            obj (dict): The dictionary to be written to file.
            file_name (str): the file name to be save too.
            rel_path (str): the relative path from self.root that the file will be saved to.
        """
        f_path = path.join(self.root,rel_path,file_name) if rel_path is not None else path.join(self.root,file_name)
        with open(f_path,"w+") as f:
            pickle.dum(obj,f)
            
    def save_index(self):
        """Writes the unique index for each of the configs to file along with
        the relative path to the atoms.json file        
        """
        with open(path.join(self.root,"index.json"),"w+") as f:
            json.dump(self.index,f)
    
    def _read_index(self):
        """Reads in the index from the index.json file if it exists.
        """
        if path.isfile(path.join(self.root,"index.json")):
            with open(path.join(self.root,"index.json"),"r") as f:
                self.index = json.load(f)
                
    @property  
    def prev(self):
        """Finds the previous group in the database.
        """
        keylist = self.parent.steps.keys()
        for i, v in enumerate(keylist):
            if v == self._db_name and i!=0:
                return self.parent.steps[keylist[i-1]]
            
    @property
    def dep(self):
        """Finds the next, or dependent, group in the databes.
        """
        keylist = self.parent.steps.keys()
        for i, v in enumerate(keylist):
            if v == self._db_name and i!=(len(keylist)-1):
                return self.parent.steps[keylist[i+1]]
        
    def is_executing(self):
        """Returns True if the database DFT calculations are in process of being
        executed.
        """
        if len(self.sequence) == 0:
            is_executing = False
            for atoms in self.config_atoms:
                is_executing = atoms.calculator.is_executing()
                if is_executing:
                    break                
        else:
            executing = [group.is_executing() for group in self.sequence.values()]
            is_executing = all(executing)
            
        return is_executing
            
    def execute(self, dryrun=False, recovery=False, env_vars=None):
        """Submits the job script for each of the folders in this
        database if they are ready to run.
        Args:
            dryrun (bool): when True, simulate the submission without
              actually submitting.
            recovery (bool): when True, submit the script for running recovery
              jobs.
            env_vars (dict): of environment variables to set before calling the
              execution. The variables will be set back after execution.
        Returns:
            bool: True if the submission generated a job id
            (considered successful).
        """

        if len(self.sequence) == 0:
            jobfile = "recovery.sh" if recovery else "jobfile.sh"
            if not path.isfile(path.join(self.root, jobfile)):
                return False

            if not recovery:
                if not all(a.calculator.can_execute()
                           for a in self.config_atoms.values()):
                    return False

                #We also need to check that we haven't already submitted this
                #job. Check to see if it is executing.
                if any(a.calculator.is_executing()
                       for a in self.config_atoms.values()):
                    return False

                #Make sure that the calculation isn't complete.
                if any(a.calculator.can_cleanup()
                       for a in self.config_atoms.values()):
                    return False                
        
            # We must have what we need to execute. Compile the command and
            # submit.
            from matdb.utility import execute
            cargs = ["sbatch", jobfile]
            if dryrun:
                from matdb.msg import okay
                okay("Executed {} in {}".format(' '.join(cargs), self.root))
                return True
            else:
                xres = execute(cargs, self.root, env_vars=env_vars)

            if len(xres["output"]) > 0 and "Submitted" in xres["output"][0]:
                from matdb.msg import okay
                okay("{}: {}".format(self.root, xres["output"][0].strip()))
                return True
            else:
                return False

        else:
            already_executed = []
            for group in self.sequence.values():
                already_executed.append(group.execute(dryrun=dryrun,
                                                      recovery=recover,
                                                      env_vars=env_vars))
            return all(already_executed)

    def recover(self, rerun=False):
        """Compiles a list of all DFT runs that didn't complete and compiles the
        `failures` file. Creates a jobfile to re-run the failed
        folders only.
        Args:
            rerun (bool): when True, recreate the jobfile even if it
              already exists. 
        """

        if len(self.sequence) == 0:
            detail = self.status(False)
            failed = [k for k, v in detail["done"].items() if not v]
            identity = "{0}|{1}".format(self._db_name, self.name)
            xpath = path.join(self.root, "failures")

            if len(failed) > 0:
                #Only write a failures file if we had failures.
                with open(xpath, 'w') as f:
                    f.write('\n'.join(failed))

                imsg = "{0}: queued {1:d} configs for recovery."
                msg.info(imsg.format(identity, len(failed)))
            else:
                msg.okay("{0}: no failures.".format(identity))
                
            #Only create a jobfile if there were actually failures
            if len(failed) > 0:
                self.jobfile(rerun, recovery=True)
            else:
                #Delete any existing recovery files from previous failures.
                from os import remove
                jobfile = path.join(self.root, "recovery.sh")
                if path.isfile(jobfile):
                    remove(jobfile)
                if path.isfile(xpath):
                    remove(xpath)
        else:
            for group in self.sequence.values():
                group.recover(rerun=rerun)
                    
    def jobfile(self, rerun=False, recovery=False):
        """Creates the job array file to run each of the sub-configurations in
        this database.
        Args:
            rerun (bool): when True, recreate the jobfile even if it
              already exists. 
            recovery (bool): when True, configure the jobfile to run
              recovery jobs for those that have previously failed. This uses a
              different template and execution path.
        """

        if len(self.sequence) == 0:
            if recovery:
                from matdb.utility import linecount
                target = path.join(self.root, "recovery.sh")
                xpath = path.join(self.root, "failures")
                asize = linecount(xpath)
            else:
                target = path.join(self.root, "jobfile.sh")
                xpath = path.join(self.root, "{}.".format(self.prefix))
                asize = len(self.configs)
            
            if path.isfile(target) and not rerun:
                return
        
            # We use the global execution parameters and then any updates
            # locally. We need to add the execution directory (including prefix) and
            # the number of jobs in the array.
            settings = self.parent.execution.copy()
            settings.update(self.execution.items())
            
            settings["execution_path"] = xpath
            settings["array_size"] = asize
            
            if "array_limit" in settings and asize < settings["array_limit"]:
                del settings["array_limit"]

            from jinja2 import Environment, PackageLoader
            env = Environment(loader=PackageLoader('matdb', 'templates'))
            if recovery:
                template = env.get_template(settings["template"].replace("array", "recovery"))
            else:
                template = env.get_template(settings["template"])
            
            with open(target, 'w') as f:
                f.write(template.render(**settings))
        else:
            for group in self.sequence.values():
                group.jobfile(rerun=rerun, recovery=recovery)
                    
    def create(self, atoms, cid=None, rewrite=False, sort=None):
        """Creates a folder within this group to calculate properties.
        Args:
            atoms (quippy.atoms.Atoms): atomic configuration to run.
            cid (int): integer configuration id; if not specified, defaults to
              the next available integer.
            rewrite (bool): when True, overwrite any existing files with the
              latest settings.
            sort (bool): when True, sort the atoms by type so that
              supercell writes work correctly.
        """

        if len(self.sequence)==0:
            uid = uuid4()
            if cid is None:
                cid = len(self.configs) + 1

            from os import path, mkdir
            target = path.join(self.root, "{}.{}".format(self.prefix, cid))
            if not path.isdir(target):
                mkdir(target)

            import pudb
            pudb.set_trace()
            calc = self.calc(atoms, target, **self.calcargs)
            calc.create()

            #Finally, store the configuration for this folder.
            self.configs[cid] = target
            self.config_atoms[cid] = atoms
        else:
            for group in self.sequence.values():
                group.create(atoms, cid=cid, rewrite=rewrite, sort=sort)

    def ready(self):
        """Determines if this database has been completely initialized *and* has
        all its computations' results ready.
        .. note:: This method should be overloaded by a sub-class.
        Raises:
            NotImplementedError: this method is intended to be overloaded by a
            sub-class.
        """
        raise NotImplementedError("Method `ready` must be overloaded by a "
                                  "sub-class.")        

    def is_setup(self):
        """Determines if all the necessary folders for sub-configurations of the seed
        atomic configuration exist.
        """        
        if len(self.sequence) == 0:
            #Test to see if we have already set the database up.
            confok = False
            if (len(self.configs) == self.nconfigs or
                len(self.configs) > 0 and self.nconfigs is None):
                imsg = "The {} database has already been setup."
                msg.info(imsg.format(self.name), 2)
                confok = True

            #Don't run setup if the program is currently executing.
            xok = False
            if self.is_executing():
                xok = True

            result = confok or xok
        else:
            already_setup = [group.is_setup() for group in self.sequence.values()]
            result = all(already_setup)

        return result
            
    def setup(self, db_setup, rerun =False):
        """Performs the setup of the database using the `db_setup` function
        passed in by the specific group instance.
        
        Args:
            db_setup (function): a function that will perform the setup for each
                group independently.
            rerun (bool): default value is False.
        """

        if self.prev is None or (self.prev is not None and self.prev.can_cleanup()):
            if len(self.sequence) == 0:
                ok = self.is_setup()
                if ok and not rerun:
                    return
                db_setup(rerun)
                with open(path.join(self.root,"compute.pkl"),"w+") as f:
                    pickle.dump({"params":self.params,"date":datetime.datetime.now(),
                                 "uuid":self.uuid},f)
            else:
                for group in self.sequence.values():
                    group.setup(db_setup,rerun=rerun)
        else:
            return False                    
            
    def status(self, print_msg=True):
        """Returns a status message for statistics of sub-configuration
        execution.
        Args:
            print_msg (bool): when True, return a text message with aggregate status
              information; otherwise, return a dict of the numbers involved
        """
        from numpy import count_nonzero as cnz
        from tqdm import tqdm
        ready = {}
        done = {}

        allatoms = self.config_atoms.values().copy()
        if len(self.sequence) > 0:
            for group in self.sequence.values():
                allatoms.update(group.config_atoms.values())
                
        for a in tqdm(allatoms):
            ready[f] = a.calculator.can_execute()
            done[f] = a.calculator.can_cleanup()
        
        rdata, ddata = cnz(ready.values()), cnz(done.values())
        N = len(self.configs)        
        is_busy = self.is_executing()

        rmsg = "ready to execute {}/{};".format(rdata, N)
        dmsg = "finished executing {}/{};".format(ddata, N)
        busy = " busy executing..." if is_busy else ""

        if print_msg:
            return "{} {}{}".format(rmsg, dmsg, busy)
        else:
            return {
                "ready": ready,
                "done": done,
                "stats": {
                    "ready": rdata,
                    "done": ddata,
                    "N": N
                },
                "busy": is_busy
            }
        
    def can_cleanup(self):
        """Runs post-execution routines to clean-up the calculations. This super class
        implementation only checks that each of the sub-config directories has
        the necessary files needed for cleanup.
        """
        if len(self.sequence) == 0:
            if (len(self.configs) != self.nconfigs and
                self.nconfigs is not None):
                #We need to have at least one folder for each config;
                #otherwise we aren't ready to go.
                return False
        
            cleanups = [a.calculator.can_cleanup() for a in self.config_atoms.values()]
        else: 
            cleanups = [group.can_cleanup() for group in self.sequence.values()]
        return all(cleanups)

    def tarball(self, filename="output.tar.gz"):
        """Creates a zipped tar archive that contains each of the specified
        files in sub-sampled configurations' output folders.
        
        Args:
            filename (str): name of the zipped archive to create.
        """
        if len(self.sequence) == 0:
            parts = []
            for fname in self.calc.tarball:
                parts.append("{}.*/{}".format(self.prefix, fname))

            targs = ["tar", "-cvzf", filename, ' '.join(parts)]
            from matdb.utility import execute
            execute(targs, self.root)
        else:
            for group in self.sequence.values():
                group.tarball(filename)

    def cleanup(self):
        """Creates a JSON file for each atoms object in the group.
        """
        if len(self.sequence) == 0 and self.can_cleanup():
            for cid, folder in self.configs.items():
                atoms = self.config_atoms[cid]
                if hasattr(atoms, calc):
                    atoms.calc.cleanup()
                atoms_to_json(atoms, folder)

        else:
            for group in self.sequence.values():
                group.cleanup()

        if self.can_cleanup():
            self.save_index()
