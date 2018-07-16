"""Implements a `matdb` compatible subclass of the
:class:`ase.calculators.vasp.Vasp` calculator.
.. note:: Because this calculator is intended to be run asynchronously as part
  of `matdb` framework, it does *not* include a method to actually execute the
  calculation. Although the ASE calculator provides an interface to do so,
  `matdb` uses templates to access HPC resources.
.. warning:: Because of the underlying implementation in ASE, you must use a separate
  instance of the :class:`AsyncVasp` for each :class:`ase.Atoms` object that you
  want to calculate for.
"""
import ase
from ase.calculators.vasp import Vasp
from os import path, stat, mkdir, remove, environ, rename
import mmap
from hashlib import sha1
import re

from matdb.calculators.basic import AsyncCalculator
from matdb import msg
from matdb.kpoints import custom as write_kpoints
from matdb.utility import chdir, execute, relpath, config_paths
from matdb.exceptions import VersionError, SpeciesError
from matdb.calculators import paths

def phonon_defaults(d, dfpt=False):
    """Adds the usual settings for the INCAR file when performing frozen-phonon
    calculations to `d`. They are only added if they weren't already specified
    in the config file.
    .. warning:: This method mutates `d`.
    Args:
        d (dict): "calculator" dictionary to updated arguments for.
        dfpt (bool): when True, perform a DFT perturbation theory calculation to
          extract the force constants; otherwise, do frozen phonons.
    """
    assert d["name"] == "Vasp"
    usuals = {
        "encut": 500,
        "ediff": '1.0e-08',
        "ialgo": 38,
        "ismear": 0,
        "lreal": False,
        "addgrid": True,
        "lwave": False,
        "lcharg": False
    }

    if dfpt:
        usuals["ibrion"] = 8
    else:
        usuals["ibrion"] = -1
        
    for k, v in usuals.items():
        if k not in d:
            d[k] = v

def extract_force_sets(configs, phonodir):
    """Extracts the force sets from a set of VASP frozen phonon calculations
    using `phonopy`.
    .. note:: This method uses `phonopy` to create the `FORCE_SETS` file; it
      does not actually return the force sets.
    Args:
        configs (dict): keys are config `id`; values are full paths to the
          folders where the configs were calculated.
        phonodir (str): full path to the `phonopy` directory for the set of
          calculations in `configs`.
    """
    #First, make sure we have `vasprun.xml` files in each of the
    #directories.
    vaspruns = []
    for i, folder in configs.items():
        vasprun = path.join(folder, "vasprun.xml")
        if not path.isfile(vasprun): #pragma: no cover
            msg.err("vasprun.xml does not exist for {}.".format(folder))
        else:
            vaspruns.append(vasprun)

    if len(vaspruns) == len(configs):
        sargs = ["phonopy", "-f"] + vaspruns
        xres = execute(sargs, phonodir, venv=True)
    else:
        xres = {"error":""}

    return xres

def extract_force_constants(configs, phonodir):
    """Extracts the force constants matrix from a single VASP DFPT calculation
    with support from `phonopy`.
    .. note:: This method uses `phonopy` to create the `FORCE_CONSTANTS` file;
      it does not actually return the force constants matrix.
    Args:
        configs (dict): keys are config `id`; values are full paths to the
          folders where the configs were calculated.
        phonodir (str): full path to the `phonopy` directory for the set of
          calculations in `configs`.
    """
    #There will only be a single config folder if we are running with
    #DFPT, since the displacements take place internally.
    vaspruns = []
    for i, folder in configs.items():
        vasprun = path.join(folder, "vasprun.xml")
        if not path.isfile(vasprun): #pragma: no cover
            msg.err("vasprun.xml does not exist for {}.".format(folder))
        else:
            vaspruns.append(vasprun)

    assert len(vaspruns) == 1

    sargs = ["phonopy", "--fc"] + vaspruns
    xres = execute(sargs, phonodir, venv=True)

    return xres
        
class AsyncVasp(Vasp, AsyncCalculator):
    """Represents a calculator that can compute material properties with VASP,
    but which can do so asynchronously.
    .. note:: The arguments and keywords for this object are identical to the
      :class:`~ase.calculators.vasp.Vasp` calculator that ships with ASE. We
      add some extra functions so that it plays nicely with `matdb`.
    Args:
        atoms (matdb.Atoms): configuration to calculate using VASP.
        folder (str): path to the directory where the calculation should take
          place.
        contr_dir (str): The absolute path of the controller's root directory.
        ran_seed (int or float): the random seed to be used for this calculator.
    
    Attributes:
        tarball (list): of `str` VASP output file names that should be included
          in an archive that represents the result of the calculation.
        folder (str): path to the directory where the calculation should take
          place.
    """
    key = "vasp"
    tarball = ["vasprun.xml"]
    pathattrs = ["potcars.directory"]

    def __init__(self, atoms, folder, contr_dir, ran_seed, *args, **kwargs):
        self.init_calc(kwargs)
        if contr_dir == '$control$':
            with open(path.join(self.contr_dir, "NAME"), "r") as f:
                name = f.readline().strip()
                self.contr_dir = config_paths[name]
        if path.isdir(contr_dir):
            self.contr_dir = contr_dir
        else:
            msg.err("{} is not a valid directory.".format(contr_dir))
                
        self.folder = path.abspath(path.expanduser(folder))
        self.kpoints = None

        # We will store all the unique potcars in a single
        # directory. If it doesn't exist then create it.
        if not path.isdir(path.join(contr_dir,"POTCARS")):
            mkdir(path.join(contr_dir,"POTCARS"))
            
        if "kpoints" in kwargs:
            self.kpoints = kwargs.pop("kpoints")            

        # remove the "potcars" section of the kwargs for latter use in
        # creation of the POTCAR file.
        self.potcars = kwargs.pop("potcars")
        if "exec_path" in kwargs:
            self.executable = kwargs.pop("exec_path")
        else:
            self.executable = None
            
        self.args = args
        self.kwargs = kwargs

        # if 'xc' was set on either the kwargs or the potcars
        # dictionaries but not the other then we need to copy it over.
        if "xc" in self.kwargs and "xc" not in self.potcars:
            self.potcars["xc"] = self.kwargs["xc"]
        elif "xc" in self.potcars and "xc" not in self.kwargs:
            self.kwargs["xc"] = self.potcars["xc"]
        elif "xc" not in self.potcars and "xc" not in self.kwargs:
            raise ValueError("'xc' must be provided as either a calculator keyword "
                             "or a potcar keyword.")

        self.ran_seed = ran_seed
        self.version = None
        super(AsyncVasp, self).__init__(**kwargs)
        if not path.isdir(self.folder):
            mkdir(self.folder)

        self.input_params["setups"] = self.potcars["setups"]
            
        self.atoms = atoms
        pot_args = self.potcars.copy()
        environ["VASP_PP_PATH"] = relpath(path.expanduser(pot_args["directory"]))
            
        self.initialize(atoms)
        # The POTCAR file is either stored in a file who's name is a
        # hashed string of the species and versions or needs to be
        # created at that location.
        this_potcar = ""
        for pot in self.ppp_list:
            for spec, ver in self.potcars["versions"].items():
                if spec in pot:
                    this_potcar += "{0}{1}".format(spec, ver)
        self.this_potcar = str(sha1(this_potcar).hexdigest())
        self._check_potcar()

    def write_input(self, atoms, directory='./'):
        """Overload of the ASE input writer that handles the k-points using our
        built-in routines.
        """
        from ase.io.vasp import write_vasp
        write_vasp(path.join(directory, 'POSCAR'),
                   self.atoms_sorted,
                   symbol_count=self.symbol_count)
        self.write_incar(atoms, directory=directory)
        self._write_potcar()
        if self.kpoints is not None:
            write_kpoints(directory, self.kpoints, self.atoms)
        elif "kspacing" not in self.kwargs:
            self.write_kpoints(directory=directory)
        self.write_sort_file(directory=directory)

    def _check_potcar(self):
        """Checks the directories needed to establish POTCAR files as symbolic
        links for computation.
        """
        if "versions" in self.potcars:
            versions = self.potcars["versions"]
        else:
            raise VersionError("POTCAR versions must be specified for each species in the "
                         "yml file.")

        if path.isfile(path.join(self.contr_dir, "POTCARS", self.this_potcar)):
            with open(path.join(self.contr_dir, "POTCARS", self.this_potcar), "r") as potfile:
                for line in potfile:
                    if self.potcars["xc"] in line and "=" not in line:
                        ver_line = line.strip().split()
                        spec = ver_line[1].split("_")[0]
                        ver = ver_line[2]
                        if spec in versions.keys():
                            if not versions[spec] == ver:
                                raise VersionError("The version {0} of the POTCAR for {1} "
                                             "does not match the requested version "
                                             "{2}".format(versions[spec], spec, ver))
                        else:
                            raise SpeciesError("The species found in the POTCAR {0} is "
                                            "not in the system being studied".format(line))

        else:
            for potfile in self.ppp_list:
                with open(potfile, "r") as f:
                    ver_line = f.readline().strip().split()
                spec = ver_line[1].split("_")[0]
                ver = ver_line[2]
                if spec in versions.keys():
                    if not versions[spec] == ver:
                        raise VersionError("The version {0} of the POTCAR for {1} "
                                           "does not match the requested version "
                                           "{2}".format(versions[spec], spec, ver))
                else:
                    raise SpeciesError("The species found in the POTCAR {0} is not in "
                                        "the system being studied".format(spec))
            
    def _write_potcar(self):
        """Makes a symbolic link between the main POTCAR file for the database
        and the folder VASP will execute in."""
        from matdb.utility import symlink
       
        POTCAR_DIR = path.join(self.contr_dir,"POTCARS")
        calc_args = self.kwargs.copy()
        calc_args["potcars"] = self.potcars
        
        # First we check to see if the POTCAR file already exists, if
        # it does then all we have to do is create the symbolic link.
        if not path.isfile(path.join(POTCAR_DIR, self.this_potcar)):
            calc = AsyncVasp(self.atoms, POTCAR_DIR, self.contr_dir, self.ran_seed, **calc_args)
            calc.write_potcar(directory=POTCAR_DIR)
            rename(path.join(POTCAR_DIR,"POTCAR"), path.join(POTCAR_DIR, self.this_potcar))
            
        symlink(path.join(self.folder, "POTCAR"), path.join(POTCAR_DIR, self.this_potcar))

    def can_execute(self, folder):
        """Returns True if the specified folder is ready to execute VASP
        in.
        """
        if not path.isdir(folder):
            return False

        sizeok = lambda x: stat(x).st_size > 25        
        required = ["INCAR", "POSCAR", "POTCAR"]
        if "kspacing" not in self.kwargs and "KSPACING" not in self.kwargs:
            required.append("KPOINTS")
            
        present = {}
        for rfile in required:
            target = path.join(folder, rfile)
            present[rfile] = path.isfile(target) and sizeok(target)

        if not all(present.values()):
            for f, ok in present.items():
                if not ok:
                    msg.info("{} not present for VASP execution.".format(f), 2)
        return all(present.values())

    def can_extract(self, folder):
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

        line = None
        with open(outcar, 'r') as f:
            # memory-map the file, size 0 means whole file
            m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)  
            # we look for 'free  energy' to verify that VASP wasn't
            # terminated during runtime for memory or time
            # restrictions
            i = m.rfind('free  energy')
            if i > 0:
                # seek to the location and get the rest of the line.
                m.seek(i)
                line = m.readline()

        if line is not None:
            return "TOTEN" in line or "Error" in line
        else:
            return False

    def is_executing(self, folder):
        """Returns True if the specified VASP folder is in process of executing.
        Args:
            folder (str): path to the folder in which the executable was run.
        """
        outcar = path.join(folder, "OUTCAR")
        outcars = path.isfile(outcar)
        busy = not self.can_extract(folder)            
        return outcars and busy

    def create(self, rewrite=False):
        """Creates all necessary input files for the VASP calculation.
        Args:
            rewrite (bool): when True, overwrite any existing files with the
              latest settings.
        """
        self.write_input(self.atoms, self.folder)

    @staticmethod
    def set_static(input_dict):
        """Converts the input dictionary to one for a static VASP calculation.
        """
        input_dict["nsw"] = 0
        return input_dict
        
    def extract(self, folder, cleanup="default"):
        """Extracts results from completed calculations and sets them on the
        :class:`ase.Atoms` object.
        Args:
            folder (str): path to the folder in which the executable was run.
            cleanup (str): the level of cleanup to perfor after extraction.
        """
        # Read output
        atoms_sorted = ase.io.read(path.join(folder,'CONTCAR'), format='vasp')

        if (self.int_params['ibrion'] is not None and
                self.int_params['nsw'] is not None):
            if self.int_params['ibrion'] > -1 and self.int_params['nsw'] > 0:
                # Update atomic positions and unit cell with the ones read
                # from CONTCAR.
                self.atoms.positions = atoms_sorted[self.resort].positions
                self.atoms.cell = atoms_sorted.cell

        # we need to move into the folder being extracted in order to
        # let ase check the convergence
        with chdir(folder):
            self.converged = self.read_convergence()
            #NB: this next method call overwrites the identity of self.atoms! We
            #don't like that behavior, so we set it right again afterwards.
            o = self.atoms

            self.set_results(self.atoms)
            self.atoms._calc = self
            E = self.get_potential_energy(atoms=self.atoms)
            F = self.forces
            S = self.atoms.get_stress(False)
            self.atoms.add_property(self.force_name, F)
            self.atoms.add_param(self.virial_name, S*self.atoms.get_volume())
            self.atoms.add_param(self.energy_name, E)

            o.copy_from(self.atoms)
            self.atoms = o
            self.atoms.set_calculator(self)
            
        self.cleanup(folder,clean_level=cleanup)

    def set_atoms(self, atoms):
        #We do *not* do an atoms.copy() like ASE because our atoms objects are hybrids.
        self.atoms = atoms
        
    def cleanup(self, folder, clean_level="default"):
        """Performs cleanup on the folder where the calculation was
        performed. The clean_level determines which files get removed.
        Args:
            folder (str): the folder to be cleaned.
            clean_level (str): the level of cleaning to be done.
        """

        light = ["CHG", "XDATCAR", "DOSCAR", "PCDAT"]
        default =["CHGCAR", "WAVECAR", "IBZKPT", 
                  "EIGENVAL", "DOSCAR", "PCDAT"]
        aggressive = ["vasprun.xml", "OUTCAR", "CONTCAR", "OSZICAR"]

        if clean_level == "light":
            rm_files = light
        elif clean_level == "aggressive":
            rm_files = light + default + aggressive
        else:
            rm_files = light + default
        
        for f in rm_files:
            target = path.join(folder,f)
            if path.isfile(target):
                remove(target)

    def to_dict(self):
        """Writes the current version number of the code being run to a
        dictionary along with the parameters of the code.
        Args:
            folder (str): path to the folder in which the executable was run.
        """
        vasp_dict = {"folder":self.folder.relpace(self.contr,'$control$'),
                     "ran_seed":self.ran_seed,
                     "contr_dir":'$control$', "kwargs": self.kwargs,
                     "args": self.args}

        # run vasp in the root directory in order to determine the
        # version number.
        if hasattr(self,"potcars"):
            potdict = self.potcars.copy()
            with open(path.join(self.contr_dir, "NAME"), "r") as f:
                name = f.readline().strip()
            namehash = str(sha1(name.encode("ASCII")).hexdigest())
            for hid, hpath in paths[namehash].items():
                if potdict["directory"] == hpath:
                    potdict["directory"] = hid
                    break
            vasp_dict["kwargs"]["potcars"] = potdict
        if hasattr(self,"kpoints"):
            vasp_dict["kwargs"]["kpoints"] = self.kpoints
        if self.version is None:
            if self.executable is not None:
                data = execute([self.executable],self.contr_dir,venv=True)
                try:
                    vasp_dict["version"] = data["output"][0].strip().split()[0]
                    self.version = vasp_dict["version"]
                except:
                    msg.warn("Error attempting to get VASP version.")
            else:
                msg.warn("VASP execution path not set in calculator. Could not "
                         "determine VASP version.")
                self.version = ""
        else:
            vasp_dict["version"] = self.version
        # Files that need to be removed after being created by the
        # vasp executable.

        files = ["CHG", "CHGCAR", "WAVECAR", "XDATCAR", "vasprun.xml", "OUTCAR",
                 "IBZKPT", "OSZICAR", "CONTCAR", "EIGENVAL", "DOSCAR", "PCDAT"]

        for f in files:
            if path.isfile(path.join(self.contr_dir,f)):
                remove(path.join(self.contr_dir,f))

        return vasp_dict
