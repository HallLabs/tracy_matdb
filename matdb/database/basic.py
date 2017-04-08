"""Abstract base class for creating and interacting with a database of
configurations for machine learning materials.
"""
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
                 prefix='S', nconfigs=100):
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

        self._nsuccess = 0
        """int: number of configurations whose output files were successfully
        converted to XYZ format. Should be equal to :attr:`nconfigs` if the
        database is complete.
        """
        
        #Try and load existing folders that match the prefix into the configs
        #list.
        from glob import glob
        from os import path, getcwd, chdir
        self.configs = {}
        current = getcwd()
        chdir(self.root)
        
        try:
            for folder in glob("{}.*".format(prefix)):
                try:
                    cid = int(folder.split('.')[1])
                    self.configs[cid] = path.join(self.root, folder)
                except:
                    #The folder name doesn't follow our convention.
                    pass
        finally:
            chdir(current)
            
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

    def jobfile(self):
        """Creates the job array file to run each of the sub-configurations in
        this database.
        """
        from os import path
        target = path.join(self.root, "jobfile.sh")
        if path.isfile(target):
            return
        
        # We use the global execution parameters and then any updates
        # locally. We need to add the execution directory (including prefix) and
        # the number of jobs in the array.
        settings = self.parent.execution.copy()
        settings.update(self.execution.items())

        xpath = path.join(self.root, "{}.".format(self.prefix))
        asize = len(self.configs)
        settings["execution_path"] = xpath
        settings["array_size"] = asize

        if "array_limit" in settings and asize < settings["array_limit"]:
            del settings["array_limit"]
        
        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template(settings["template"])
        with open(target, 'w') as f:
            f.write(template.render(**settings))
                    
    def create(self, atoms, cid=None, rewrite=False):
        """Creates a folder within this database in which VASP may be run.

        Args:
            atoms (quippy.atoms.Atoms): atomic configuration to run.
            cid (int): integer configuration id; if not specified, defaults to
              the next available integer.
            rewrite (bool): when True, overwrite any existing files with the
              latest settings.
        """
        if cid is None:
            cid = len(self.configs) + 1

        from os import path, mkdir
        target = path.join(self.root, "{}.{}".format(self.prefix, cid))
        if not path.isdir(target):
            mkdir(target)

        #Now, just generate the POSCAR file.
        from ase.io import write
        write(path.join(target, "POSCAR"), atoms, "vasp")

        #Create symbolic links to the INCAR and POTCAR files that we need. INCAR
        #is stored locally for each database type (in `self.root`) while the
        #POTCAR is for the entire system and lives two directories up.
        from matdb.utility import symlink, execute
        #Make sure that the INCAR and PRECALC for this database has been created
        #already.
        self.INCAR(rewrite)
        self.PRECALC(rewrite)
        INCAR = path.join(target, "INCAR")
        POTCAR = path.join(target, "POTCAR")
        PRECALC = path.join(target, "PRECALC")
        symlink(INCAR, path.join(self.root, "INCAR"))
        symlink(PRECALC, path.join(self.root, "PRECALC"))
        symlink(POTCAR, path.join(self.parent.parent.root, "POTCAR"))

        execute(["getKPoints"], target)
        
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
        if len(self.configs) == self.nconfigs:
            msg.info("The phonon-base database has already been setup.", 2)
            return True
        else:
            return False

    def cleanup(self):
        """Runs post-DFT execution routines to clean-up the database.

        .. note:: This method should be overloaded by a sub-class.

        Raises:
            NotImplementedError: this method is intended to be overloaded by a
            sub-class.
        """
        raise NotImplementedError("Method `cleanup` must be overloaded by a "
                                  "sub-class.")
    
    def xyz(self, filename="output.xyz",
            properties=["species", "pos", "z", "force"],
            parameters=["energy", "stress"]):
        """Creates an XYZ file for all the sub-sampled configurations in this
        database.

        Args:
            filename (str): name of the XYZ file to create; this is created in
              each sub-sampled configurations directory.
            properties (list): of `str` *atom* property names (such as position,
              force, Z, etc.) to include in the XYZ file.
            parameters (list): of `str` *configuration* paramater names (such as
              energy, stress, etc.).

        Returns:
            bool: True if the number of xyz files created equals the number of
            configurations in the database, which means that the database is
            fully calculated in a usable way.
        """
        from os import path
        p = ','.join(properties)
        P = ','.join(parameters)
        sargs = ["convert.py", "-I", "OUTCAR", "-p", p, "-P", P, "-f", "xyz",
                 "OUTCAR", "-o", filename]

        from matdb.utility import execute
        created = []
        for i, folder in self.configs.items():
            relpath = path.join(self.root, folder)
            execute(sargs, relpath)
            outpath = path.join(relpath, filename)
            if path.isfile(outpath):
                created.append(outpath)            

        #Finally, combine all of them together into a single
        from matdb.utility import cat
        cat(created, path.join(self.root, filename))
        self._nsuccess = len(created)
        
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
