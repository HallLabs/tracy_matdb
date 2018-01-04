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
        atoms (optional, quippy.atoms.Atoms): a single atomic configuration from
          which many others may be derived using MD, phonon
          displacements, etc.

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
    def __init__(self, root=None, parent=None, prefix='S', atoms=None,
                 nconfigs=None, calculator=None, seed=None, db_name=None,
                 config_type=None, parameters=None, execution={}):
        from collections import OrderedDict
        from quippy.atoms import Atoms
        from os import path

        self.parent = parent
        self.execution = execution
        if parameters is not None:
            self.params = parameters
        else:
            self.params = None

        self.is_seed_explicit = None
        if seed is not None:
            self.is_seed_explicit = isinstance(seed, list)

        self.seeds = None
        if isinstance(seed, list):
            self.seeds = []
            for atomspec in seed:
                fmt, pattern = atomspec.split(':')
                for apath in self.parent.controller.relpaths(pattern):
                    self.seeds.append((apath, fmt))
                    
        elif seeds is not None:
            self.seeds = [(a.copy(), None) for a in self.prev.rset]
                                  
        self.sequence = OrderedDict()
        if self.seeds is not None:
            for n_seeds, (seedpath, seedfmt) in enumerate(self.seeds):
                if self.is_seed_explicit:
                    this_atoms = seedpath
                    seedname = "seed-{}".format(n_seeds)
                else:
                    this_atoms = Atoms(seedpath, format=seedfmt)
                    seedname = path.basename(this_seed)
                seed_root = path.join(root, seedname)
                    
                if not path.isdir(seed_root):
                    mkdir(seed_root)
                                        
                if atoms is not None and isinstance(atoms, ParameterGrid):
                    for params in atoms:
                        pkey = "{0}/{1}".format(seedname, atoms.to_str(params))
                        this_root = path.join(seed_root, pkey)
                        if not path.isdir(this_root):
                            mkdir(this_root)
                        self.sequence[pkey] = Group(root=this_root,
                                                    parent=parent,
                                                    prefix=prefix,
                                                    atoms=this_atoms,
                                                    db_name=db_name,
                                                    calculator=calculator,
                                                    parameters=atoms[params])
                        
                elif atoms is not None and isinstance(atoms,Atoms):
                    self.atoms = atoms
        else:
            if atoms is not None and isinstance(atoms, ParameterGrid):
                for params in atoms:
                    this_root = path.join(root,atoms.to_str(params))
                    if not path.isdir(this_root):
                        mkdir(this_root)
                    self.sequence[atoms.to_str(params)] = Group(root=this_root,parent=parent,
                                                                prefix=prefix,
                                                                db_name=db_name,
                                                                calculator=calculator,
                                                                parameters=atoms[params])
            elif atoms is not None and isinstance(atoms, Atoms):
                self.atoms = atoms
            else:
                self.atoms = None

        self.calc = getattr(calculators, calculator["name"])
        self.calcargs = calculator
        self.root = root            
        self.prefix = prefix
        self.nconfigs = nconfigs
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
        
        with chdir(self.root):
            for folder in glob("{}.*".format(prefix)):
                try:
                    cid = int(folder.split('.')[1])
                    self.configs[cid] = path.join(self.root, folder)
                    self.config_atoms[cid] = self.calc.from_folder(self.configs[cid]).atoms
                except:
                    #The folder name doesn't follow our convention.
                    pass

    @abc.abstractproperty
    def rset(self):
        pass
            
    def prev(self):
        """Finds the previous group in the database.
        """
        keylist = self.parent.steps.keys()
        for i, v in enumerate(keylist):
            if v == self._db_name and i!=0:
                return self.parent.steps[keylist[i-1]]

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
            executing = []
            for seq in self.sequence:
                executing.append(self.sequence[seq].is_executing())
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
            for seq in self.sequence:
                already_executed.append(self.sequence[seq].execute(dryrun=dryrun,
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
            for seq in self.sequence:
                self.sequence[seq].recover(rerun=rerun)
                    
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
            for seq in self.sequence:
                self.sequence[seq].jobfile(rerun=rerun, recovery=recovery)
                    
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

        if not bool(self.sequence):
            if cid is None:
                cid = len(self.configs) + 1

            from os import path, mkdir
            target = path.join(self.root, "{}.{}".format(self.prefix, cid))
            if not path.isdir(target):
                mkdir(target)

            calc = self.calc(atoms, target, **self.calcargs)
            calc.create()

            #Finally, store the configuration for this folder.
            self.configs[cid] = target
            self.config_atoms[cid] = atoms
        else:
            for seqkey, group in self.sequence.items():
                group.create(group.atoms, cid=cid, rewrite=rewrite, sort=sort)

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
            already_setup = []
            for seq in self.sequence:
                already_setup.append(self.sequence[seq].setup())
            result = all(already_setup)

        return result
            
    @contextmanager
    def setup(self):
        ok = self.is_setup()
        yield ok
        if not ok and self.is_setup():
            with open(path.join(self.root,"compute.pkl"),"w+") as f:
                pickle.dump({"params":self.params,"date":datetime.datetime.now(),
                             "uuid":None},f)
            
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
        if (len(self.configs) != self.nconfigs and
            self.nconfigs is not None):
            #We need to have at least one folder for each config;
            #otherwise we aren't ready to go.
            return False
        
        cleanups = [a.calculator.can_cleanup() for a in self.config_atoms.values()]
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
