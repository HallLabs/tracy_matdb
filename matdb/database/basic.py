"""Abstract base class for creating and interacting with a database of
configurations for machine learning materials.
"""
from os import path
from matdb import msg

def can_execute(folder):
    """Returns True if the specified folder is ready to execute VASP
    in.
    """
    if not path.isdir(folder):
        return False
    
    from os import stat
    required = ["INCAR", "POSCAR", "KPOINTS", "POTCAR"]
    present = {}
    for rfile in required:
        target = path.join(folder, rfile)
        sizeok = stat(target).st_size > 25
        present[rfile] = path.isfile(target) and sizeok

    if not all(present.values()):
        from matdb.msg import info
        for f, ok in present.items():
            if not ok:
                info("{} not present for VASP execution.".format(f), 2)
    return all(present.values())

def can_cleanup(folder):
    """Returns True if the specified VASP folder has completed
    executing and the results are available for use.
    """
    if not path.isdir(folder):
        return False
    
    #If we can extract a final total energy from the OUTCAR file, we
    #consider the calculation to be finished.
    outcar = path.join(folder, "OUTCAR")
    if not path.isfile(outcar):
        return False

    import mmap
    line = None
    with open(outcar, 'r') as f:
        # memory-map the file, size 0 means whole file
        m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)  
        i = m.rfind('free  energy')
        if i > 0:
            # seek to the location and get the rest of the line.
            m.seek(i)
            line = m.readline()

    if line is not None and "TOTEN" in line:
        return True
    else:
        #For some of the calibration tests, it is possible that an error was
        #raised because the ions are too close together. That still counts as
        #finishing, we just don't have output.
        line = None
        with open(outcar, 'r') as f:
            # memory-map the file, size 0 means whole file
            m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)  
            i = m.rfind('free  energy')
            if i > 0:
                # seek to the location and get the rest of the line.
                m.seek(i)
                line = m.readline()

        if line is not None and "Error" in line:
            return True
        else:
            return False
        
class Database(object):
    """Represents a collection of material configurations (varying in
    structure and composition) from which a machine learning model can
    be created. Includes logic for generating the DFT directories that
    need to be run as well as extracting the relevant data from such
    computations.

    Args:
        atoms (quippy.atoms.Atoms): a single atomic configuration from
          which many others may be derived using MD, phonon
          displacements, etc.
        incar (dict): key-value pairs of settings for INCAR when
          running this particular database.
        kpoints (dict): key-value pairs of settings for the Mueller
          KPOINTS selection for this database.
        execution (dict): key-value pairs of settings for the supercomputer job
          array batch file.
        root (str): full path to the root directory that this database will live
          in.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        prefix (str): sub-sampled configurations will be stored using integer
          ids after this prefix; for example `S.1`, `S.2`, etc.
        nconfigs (int): number of displaced configurations to create.

    Attributes:
        atoms (quippy.atoms.Atoms): a single atomic configuration from
          which many others may be derived using MD, phonon
          displacements, etc.
        incar (dict): key-value pairs of settings for INCAR when
          running this particular database.
        kpoints (dict): key-value pairs of settings for the Mueller
          KPOINTS selection for this database.
        execution (dict): key-value pairs of settings for the supercomputer job
          array batch file.
        configs (dict): keys are integer identifiers for the particular
          configuration; values are paths (relative to the base atoms root
          directory) in which calculations are performed.
        root (str): full path to the root directory that this database will live
          in.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        prefix (str): sub-sampled configurations will be stored using integer
          ids after this prefix; for example `S.1`, `S.2`, etc.
        nconfigs (int): number of displaced configurations to create.
    """
    def __init__(self, atoms, incar, kpoints, execution, root, parent,
                 prefix='S', nconfigs=100, config_type=None):
        self.atoms = atoms
        self.incar = incar.copy()
        self.kpoints = kpoints.copy()
        self.execution = execution.copy()
        self.root = root
        from os import path, mkdir
        if not path.isdir(self.root):
            mkdir(self.root)
            
        self.parent = parent
        self.prefix = prefix
        self.nconfigs = nconfigs
        self.config_type = config_type
        self._nsuccess = 0
        """int: number of configurations whose output files were successfully
        converted to XYZ format. Should be equal to :attr:`nconfigs` if the
        database is complete.
        """
        
        #Try and load existing folders that match the prefix into the configs
        #list.
        from glob import glob
        from os import path
        from matdb.utility import chdir
        self.configs = {}

        with chdir(self.root):
            for folder in glob("{}.*".format(prefix)):
                try:
                    cid = int(folder.split('.')[1])
                    self.configs[cid] = path.join(self.root, folder)
                except:
                    #The folder name doesn't follow our convention.
                    pass

    def is_executing(self):
        """Returns True if the database DFT calculations are in process of being
        executed.
        """
        outcars = False
        for f in self.configs.values():
            outcar = path.join(f, "OUTCAR")
            if path.isfile(outcar):
                outcars = outcars or True

        busy = False
        for f in self.configs.values():
            if not can_cleanup(f):
                busy = True
                break
            
        return outcars and busy
            
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
        jobfile = "recovery.sh" if recovery else "jobfile.sh"
        if not path.isfile(path.join(self.root, jobfile)):
            return False

        if not recovery:
            executors = {f: can_execute(f) for f in self.configs.values()}
            if not all(executors.values()):
                return False

            #We also need to check that we haven't already submitted this job. If
            #the OUTCAR file exists, then we don't want to resubmit.
            for f in self.configs.values():
                outcar = path.join(f, "OUTCAR")
                if path.isfile(outcar):
                    return False
        
        # We must have what we need to execute. Compile the command
        # and submit.
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
            
    def PRECALC(self, rewrite=False):
        """Creates the INCAR file using the default settings plus any extras
        that were passed in to this database.

        Args:
            rewrite (bool): when True, overwrite any existing PRECALC with the
              latest settings.
        """
        from os import path
        target = path.join(self.root, "PRECALC")

        if rewrite or not path.isfile(target):
            with open(target, 'w') as f:
                for k, v in self.parent.kpoints.items():
                    #Skip any keywords that are overridden by the
                    #local database settings.
                    if k not in self.kpoints:
                        f.write("{}={}\n".format(k.upper(), v))

                for k, v in self.kpoints.items():
                    f.write("{}={}\n".format(k.upper(), v))        
            
    def INCAR(self, rewrite=False):
        """Creates the INCAR file using the default settings plus any extras
        that were passed in to this database.

        Args:
            rewrite (bool): when True, overwrite any existing INCAR with the
              latest settings.
        """
        from os import path
        target = path.join(self.root, "INCAR")

        if rewrite or not path.isfile(target):
            with open(target, 'w') as f:
                for k, v in self.parent.incar.items():
                    #Skip any keywords that are overridden by the
                    #local database settings.
                    if k not in self.incar:
                        f.write("{}={}\n".format(k.upper(), v))

                for k, v in self.incar.items():
                    f.write("{}={}\n".format(k.upper(), v))

    def recover(self, rerun=False):
        """Compiles a list of all DFT runs that didn't complete and compiles the
        `failures` file. Creates a jobfile to re-run the failed
        folders only.

        Args:
            rerun (bool): when True, recreate the jobfile even if it
              already exists. 
        """
        detail = self.status(False)
        failed = [k for k, v in detail["done"].items() if not v]
        identity = "{0}|{1}".format(self.parent.name, self.name)
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
        from os import path
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
                    
    def create(self, atoms, cid=None, rewrite=False, sort=None):
        """Creates a folder within this database in which VASP may be run.

        Args:
            atoms (quippy.atoms.Atoms): atomic configuration to run.
            cid (int): integer configuration id; if not specified, defaults to
              the next available integer.
            rewrite (bool): when True, overwrite any existing files with the
              latest settings.
            sort (bool): when True, sort the atoms by type so that
              supercell writes work correctly.
        """
        if cid is None:
            cid = len(self.configs) + 1

        from os import path, mkdir
        target = path.join(self.root, "{}.{}".format(self.prefix, cid))
        if not path.isdir(target):
            mkdir(target)

        #Now, just generate the POSCAR file.
        from ase.io import write
        write(path.join(target, "POSCAR"), atoms, "vasp", sort=sort)

        #Create symbolic links to the INCAR and POTCAR files that we need. INCAR
        #is stored locally for each database type (in `self.root`) while the
        #POTCAR is for the entire system and lives two directories up.
        from matdb.utility import symlink, execute
        #Make sure that the INCAR and PRECALC for this database has been created
        #already.
        self.INCAR(rewrite)
        INCAR = path.join(target, "INCAR")
        symlink(INCAR, path.join(self.root, "INCAR"))
        POTCAR = path.join(target, "POTCAR")
        symlink(POTCAR, path.join(self.parent.parent.root, "POTCAR"))

        if "matdb" not in self.kpoints:
            self.PRECALC(rewrite)
            PRECALC = path.join(target, "PRECALC")
            symlink(PRECALC, path.join(self.root, "PRECALC"))
            execute(["getKPoints"], target)
        else:
            from matdb.kpoints import custom
            custom(self.root, self.kpoints["matdb"], self.atoms)
            KPOINTS = path.join(target, "KPOINTS")
            symlink(KPOINTS, path.join(self.root, "KPOINTS"))            
        
        #Finally, store the configuration for this folder.
        self.configs[cid] = target

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
        
    def setup(self):
        """Creates all the necessary folders for sub-configurations of the seed
        atomic configuration, preparatory to DFT computations.

        .. note:: This method should be overloaded by a sub-class, which also
          calls this method.
        """
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
            
        return confok or xok
            
    def status(self, print_msg=True):
        """Returns a status message for statistics of sub-configuration execution
        with VASP.

        Args:
            print_msg (bool): when True, return a text message with aggregate status
              information; otherwise, return a dict of the numbers involved
        """
        from numpy import count_nonzero as cnz
        from tqdm import tqdm
        ready = {}
        done = {}
        for f in tqdm(self.configs.values()):
            ready[f] = can_execute(f)
            done[f] = can_cleanup(f)
        
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
        
    def cleanup(self):
        """Runs post-DFT execution routines to clean-up the database. This super class
        implementation only checks that each of the sub-config directories has
        the necessary files needed for cleanup.
        """
        if (len(self.configs) != self.nconfigs and
            self.nconfigs is not None):
            #We need to have at least one folder for each config;
            #otherwise we aren't ready to go.
            return False
        
        cleanups = [can_cleanup(f) for f in self.configs.values()]
        return all(cleanups)
    
    def xyz(self, filename="output.xyz",
            properties=["species", "pos", "z", "dft_force"],
            parameters=["dft_energy", "dft_virial"],
            recalc=False, combine=False):
        """Creates an XYZ file for each of the sub-sampled configurations in this
        database.

        Args:
            filename (str): name of the XYZ file to create; this is created in
              each sub-sampled configurations directory.
            properties (list): of `str` *atom* property names (such as position,
              force, Z, etc.) to include in the XYZ file.
            parameters (list): of `str` *configuration* paramater names (such as
              energy, stress, etc.).
            recalc (bool): when True, re-create the XYZ files, even if they already
              exist. 
            combine (bool): when True, combine all the sub-configuration XYZ
              files into a single, giant XYZ file.

        Returns:
            bool: True if the number of xyz files created equals the number of
            configurations in the database, which means that the database is
            fully calculated in a usable way.
        """
        from matdb.io import vasp_to_xyz
        from tqdm import tqdm
        created = []
        for i, folder in tqdm(self.configs.items()):
            if vasp_to_xyz(folder, filename, recalc, properties, parameters):
                outpath = path.join(folder, filename)
                created.append(outpath)
                
                if config_type is not None:
                    #We need to load the XYZ file, add the config_type
                    #parameter and then save it again.
                    from quippy.atoms import Atoms
                    a = Atoms(outpath)
                    a.params["config_type"] = self.config_type
                    a.write(outpath)
                    
        self._nsuccess = len(created)
        
        #Finally, combine all of them together into a single
        if combine:
            from matdb.utility import cat
            dboutpath = path.join(self.root, filename)
            cat(created, dboutpath)
        
        return len(created) == len(self.configs)

    def tarball(self, filename="output.tar.gz", files=["OUTCAR"]):
        """Creates a zipped tar archive that contains each of the specified
        files in sub-sampled configurations' output folders.
        
        Args:
            filename (str): name of the zipped archive to create.
            files (list): of `str` files in each sub-sampled folder to include
              in the archive.
        """
        parts = []
        for fname in files:
            parts.append("{}.*/{}".format(self.prefix, fname))

        targs = ["tar", "-cvzf", filename, ' '.join(parts)]
        from matdb.utility import execute
        execute(targs, self.root)
