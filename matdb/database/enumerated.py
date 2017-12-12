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
       arrows (list): list of arrow restrictions.
       concs (list): a list of the concetration restrictions.
       eps (float): the floating point tolerance
       arrow_res (bool): True if arrows are present.
       conc_res (bool): True if concetrations are being restricted.
       lattice (list): the lattice vectors for the system.
       basis (list): the atomic basis for the system.

    """
    def __init__(self, atoms=None, root=None, controller=None, parent=None, incar={},
                 kpoints={}, execution={}, sizes=None, basis=None, lattice=None, concs=None,
                 arrows = None, eps=None, nconfigs=None, species=None, name="enum", seed=None):

        self.name = name
        super(EnumDatabase, self).__init__(atoms,incar,kpoints,execution,
                                           path.join(root,self.name),
                                           parent,"E",nconfigs=None)
        self.nconfigs = nconfigs
        if eps is None:
            self.eps = 10**(-3)
        else:
            self.eps = eps
            
        #determine the lattice.
        if isinstance(lattice,str):
            if lattice.lower() == "fcc":
                self.lattice = [[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]]
            elif lattice.lower() == "bcc":
                self.lattice = [[0.5,0.5,-0.5],[0.5,-0.5,0.5],[-0.5,0.5,0.5]]
            elif lattice.lower() == "sc":
                self.lattice = [[1,0,0],[0,1,0],[0,0,1]]
            elif lattice.lower() == "hcp":
                self.lattice = [[1,0,0],[.5, 0.866025403784439, 0],[0, 0, 1.6329931618554521]]
            else:
                msg.err("The lattice type {} is unsupported. Please enter your lattice vectors "
                        "as a 3x3 matrix with the vectors as rows in the config file "
                        "(i.e. [a1,a2,a3]).".format(lattice))
        elif isinstance(lattice,list):
            if len(lattice) == 3 and len(lattice[0]) == 3:
                self.lattice = lattice
            else: 
                msg.err("The lattice vectors must be a 3x3 matrix with the vectors as rows.")
        else:
            msg.err("The lattice vectors must either be a string of 'sc', 'fcc', 'hcp', 'bcc', "
                    "or a 3x3 matrix with the vectors a rows, not {}".format(lattice))

        if isinstance(basis,list):
            if len(basis[0]) == 3:
                self.basis = basis
            else:
                msg.err("The atomic basis must be a list of lists that is nx3 where n is "
                        "the number of atoms in the basis.")
        elif basis is None:
            if lattice.lower() != "hcp":
                self.basis = [[0,0,0]]
            else:
                self.basis = [[0,0,0],[0.5,0.28867513459,0.81649658093]]                
        else:
            msg.err("The atomic basis must be a list of lists that is nx3 where n is "
                    "the number of atoms in the basis or left blank.")

        if isinstance(species,list):
            self.knary = len(species)
            if len(self.knary) == 1:
                self.species = [species[0],species[0]]
            else:
                self.species = species
        else:
            msg.err("The species must be a list of atomic elements.")
            

        if concs is not None and len(concs) == len(species):
            if len(concs[0]) == 3:
                self.concs = concs
                self.conc_res = True
            else:
                msg.err("The concetrations must be a nx3 list.")
        elif concs is None:
            self.concs = concs
            self.conc_res = False
        else:
            msg.err("The number of species and the concentrations must be have the "
                    "same length.")

        if arrows is not None and len(arrows) == len(species):
            if isinstance(arrows[0],(int,float)) and arrows[0]<=1:
                self.arrows = arrows
                self.arrow_res = False
            else :
                msg.err("The arrows must be a list of values <= 1.")
        elif arrows is None:
            self.arrows = arrows
            self.arrow_res = False
        else:
            msg.err("The number of species and arrow concentrations must have the "
                    "same length.")
            
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

    def _build_lattice_file(self):
        """Creates the lattice.is file that phenum needs in order to perform
        the enumeration.
        
        """

        from os import path
        target = path.join(self.root, "lattice.in")

        # We need to create a dictionary of the arguments to pass into
        # the template.
        settings = {}
        settings["template"] = "lattice.in"
        settings["min_cell_size"] = self.min_size
        settings["max_cell_size"] = self.max_size
        if self.conc_res:
            settings["conc_res"] = "T"
            if self.arrow_res:
                temp = []
                for i, a in enumerate(self.arrows):
                    temp.append("{0} {1}".format(" ".join([str(j) in self.concs[i]]),a))
            else:
                temp = [" ".join([str(i) for i in j]) for j in self.concs]
                    
            settings["concetrations"] = temp
        else:
            settings["conc_res"] = "F"
        if self.arrow_res:
            settings["incl_arrows"] = "T"
        else:
            settings["incl_arrows"] = "F"
        settings["lattice"] = [" ".join([str(i) for i in j]) for j in self.lattice]
        settings["k_nary"] = self.knary
        settings["atomic_basis"] = [" ".join([str(i) for i in j]) for j in self.basis]
        settings["n_basis"] = len(self.basis)        
        
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
            dind = 0
            self._build_lattice_file(lattices[i])
            outfile = "-outfile enum.out"
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
            

