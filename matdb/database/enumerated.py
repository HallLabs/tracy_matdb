"""Database of configurations that is created from an enumerated list of structures.
"""
from .basic import Database
from matdb import msg
from os import path
import numpy as np

class EnumDatabase(Database):
    """Sets up the calculations for a random sampling of structures from
    an enumerated list.

    Args:
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        incar (dict): specify additional settings for the INCAR file (i.e.,
          differing from, or in addition to those in the global set).
        kpoints (dict): specify additional settings for the PRECALC file (i.e.,
          differing from, or in addition to those in the global set).
        execution (dict): specify override parameters for the execution
          templates in this database.
        max_size (int): largest multiple of the primitive cell to be used
          in the enumeration.
        min_size (int): largest multiple of the primitive cell to be used
          in the enumeration.
        n_systems (int): the number of systems to be selected from the enumerated
          list to form the initial database.
        species (list): the atomic species included in the system.
        seed (hashable): a seed to feed to the random number generator.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Database`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
       nconfigs (int): the number of configurations to include in the datababes.
       max_size (int): the largest allowed cell size in the database.
       min_size (int): the smallest allowed cell size in the database.
       knary (int): the number of atomic species in the enumeration.
       species (list): the atomic species in the system

    """
    def __init__(self, atoms=None, root=None, controller=None, parent=None, incar={},
                 kpoints={}, execution={}, sizes=None,
                 nconfigs=None, species=None, name="enum", seed=None):

        self.name = name
        super(EnumDatabase, self).__init__(atoms,incar,kpoints,execution,
                                           path.join(root,self.name),
                                           parent,"E",nconfigs=None)
        self.nconfigs = nconfigs
        if len(sizes)==1 and isinstance(sizes,list):
            self.min_size = 1
            self.max_size = sizes[0]
        elif len(sizes)==2 and isinstance(sizes,list):
            self.min_size = sizes[0]
            self.max_size = sizes[1]
        else:
            raise ValueError("The sizes specified must be a list of 1 or 2 values."
                             "If one value then it must be the largest cell size to include,"
                             "i.e., [32]. If 2 values then it the first value is the smallest "
                             "and the second value is the largest cell size to include, "
                             "i.e., [10,12].")

        self.seed = seed
        self.knary = len(species)
        self.species = species
        self.incar = incar
                
    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        from numpy import count_nonzero as cnz
        from tqdm import tqdm
        from matdb.database.basic import can_cleanup

        # If we haven't made a folder for each config then we aren't ready.
        if len(self.configs) != self.nconfigs:
            return False
        
        done = {}
        for f in tqdm(self.configs.values()):
            done[f] = can_cleanup(f)

        ddata = cnz(done.values())
        N = len(self.configs)        

        if ddata == N:
            return True
        else:
            return False
    
    def cleanup(self):
        """If possible cleans up the DFT runs to ensure that the database ran to
        completion.

        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.

        """

        if not self.ready():
            return False

        if not super(EnumDatabase, self).cleanup():
            return False
        
        return True

    def _build_lattice_file(self,lattice):
        """Creates the lattice.is file that phenum needs in order to perform
        the enumeration.
        
        Args:
            lattice (str): The string that identifies which lattice to use as the
              parent lattice, i.e., 'fcc', 'bcc', 'sc', 'hcp'.
        """

        from os import path
        target = path.join(self.root, "lattice.in")

        # We need to create a dictionary of the arguments to pass into
        # the template.
        settings = {}
        settings["template"] = "lattice.in"
        settings["min_cell_size"] = self.min_size
        settings["max_cell_size"] = self.max_size
        settings["conc_res"] = "F"
        settings["incl_arrows"] = "F"
        settings["k_nary"] = self.knary

        if lattice.lower() == "sc":
            settings["vec_1"] = "1 0 0"
            settings["vec_2"] = "0 1 0"
            settings["vec_3"] = "0 0 1"
            settings["n_basis"] = "1"
            settings["atomic_basis"] = ["0 0 0"]
        elif lattice.lower() == "fcc":
            settings["vec_1"] = "0 0.5 0.5"
            settings["vec_2"] = "0.5 0 0.5"
            settings["vec_3"] = "0.5 0.5 0"
            settings["n_basis"] = "1"
            settings["atomic_basis"] = ["0 0 0"]
        elif lattice.lower() == "bcc":
            settings["vec_1"] = "-0.5 0.5 0.5"
            settings["vec_2"] = "0.5 -0.5 0.5"
            settings["vec_3"] = "0.5 0.5 -0.5"
            settings["n_basis"] = "1"
            settings["atomic_basis"] = ["0 0 0"]
        elif lattice.lower() == "hcp":
            settings["vec_1"] = "1 0 0"
            settings["vec_2"] = ".5 0.866025403784439 0"
            settings["vec_3"] = "0 0 1.6329931618554521"
            settings["n_basis"] = "2"
            settings["atomic_basis"] = ["0 0 0","0.5 0.28867513459  0.81649658093"]
        else:
            raise ValueError("The {} lattice type is unrecognized.".format(lattice))
        
        
        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template("lattice.in")

        with open(target,'w') as f:
            f.write(template.render(**settings))

    def setup(self, rerun=False):

        """Sets up the initial DFT calculations.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """
        from os import chdir, system, getcwd
        from glob import glob
        folders_ok = super(EnumDatabase, self).setup()
        if folders_ok and not rerun:
            return

        #We also don't want to setup again if we have the results already.
        if self.ready():
            return

        if not folders_ok:
            # We need to determine how many configurations from each
            # lattice type to include.
            sub_nconfigs = [self.nconfigs//3]*3
            if sum(sub_nconfigs) != self.nconfigs:
                sub_nconfigs[0] += (self.nconfigs-sum(sub_nconfigs))

            # For each lattice type we need to construct a lattice.in
            # file then run phenum so that each system gets the
            # correct number of configurations.
            current = getcwd()
            chdir(self.root)
            lattices = ["fcc","bcc","hcp"]
            dind = 0
            for i in range(3):
                self._build_lattice_file(lattices[i])
                outfile = "-outfile enum.out.{}".format(lattices[i])
                dist = "all {}".format(sub_nconfigs[i])
                this_sys = "-species {}".format(" ".join(self.species))
                # Perform the enumeration.
                if self.seed is None:
                    system("enumeration.py -enum -distribution {0} {1}".format(dist,outfile))
                else:
                    system("enumeration.py -enum -distribution {0} {1} -seed {2}".format(dist,outfile,self.seed))
                    
                system("rm polya* enum.in")
                # extract the POSCARS
                system("makeStr.py all -input enum.out.{0} {1}".format(lattices[i],this_sys))

                # Now we need to 
                from quippy.atoms import Atoms
                for dposcar in glob("vasp.*"):
                    datoms = Atoms(dposcar,format="POSCAR")
                    self.create(datoms,cid=dind)
                    dind += 1
                    
            chdir(current)

        # Last of all, create the job file to execute the job array.
        self.jobfile(rerun)
            

