"""Group of configurations that is created from an enumerated list of structures.
"""
from .basic import Group
from matdb import msg
from os import path, getcwd, chdir, remove, listdir, mkdir
import numpy as np
from six import string_types

class Enumerated(Group):
    """Sets up the calculations for a random sampling of structures from
    an enumerated list.

    Args:
        sizes (list): a list containing the smallest and larges cell sizes to 
            include in the database.
        basis (list): the atomic basis to use for the enumeration.
        lattice (list or str): either the name of the lattice to use 
            ('sc','fcc','bcc', 'hcp') or the atomic basis vectors to use.
        concs (list): the concentrations of each atomic species.
        arrows (list): the maximum number of atoms to displace.
        eps (float): floating point tolerance for comparisons.
        species (list): the atomic species included in the system.
        rseed (hashable): a seed to feed to the random number generator.
        rattle (float): the amount to rattle the atoms by.
        keep_supers (bool): True if the superperiodic cells are to be kept in the
            enumerated list.
        displace (float): the amount to displace atoms with arrows.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
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
       rattle (float): the amount to rattle the atoms in the config by.

    """
    def __init__(self, sizes=None, basis=None, lattice=None, concs=None,
                 arrows = None, eps=None, species=None, name="enum",
                 rattle=None, rseed=None, keep_supers=None, displace=None,
                 **dbargs):
        self.name = name
        dbargs['prefix'] = "E"
        dbargs['cls'] = Enumerated
        if "Enum" not in dbargs['root']:
            mkdir("Enum")
            chdir("Enum")
            dbargs['root'] = path.join(dbargs['root'],"Enum")
        super(Enumerated, self).__init__(**dbargs)
        
        if eps is not None:
            self.eps = 10**(-3)
        else:
            self.eps = eps

        self.rseed = rseed
        self._get_lattice(lattice)
        self._get_basis(basis)
        self._get_species(species)
        self._get_concs(concs)
        self._get_arrows(arrows)
        if rattle is None:
            self.rattle = 0.0
        else:
            self.rattle = float(rattle)
        if keep_supers  is None:
            self.keep_supers = False
        else:
            self.keep_supers = keep_supers
            
        if displace is None:
            self.displace = 0.0
        else:
            self.displace = displace
        
        if sizes is not None and len(sizes)==1 and isinstance(sizes,list):
            self.min_size = 1
            self.max_size = sizes[0]
        elif sizes is not None and len(sizes)==2 and isinstance(sizes,list):
            self.min_size = sizes[0]
            self.max_size = sizes[1]
        elif sizes is not None:
            raise ValueError("The sizes specified must be a list of 1 or 2 values."
                             "If one value then it must be the largest cell size to include,"
                             "i.e., [32]. If 2 values then it the first value is the smallest "
                             "and the second value is the largest cell size to include, "
                             "i.e., [10,12].")
        self.euids = None
        self._load_euids()

    def _get_lattice(self,lattice):
        """Gets the lattice vectors for the system.
        
        Args:
            lattice (str or list): either a string containing the lattice 
                name or a 3x3 list of the vectors as [a1,a2,a3].
        """
        # determine the lattice.
        if isinstance(lattice,string_types):
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
            if len(lattice) == 3 and len(lattice[0]) == 3 and np.linalg.det(lattice) != 0:
                self.lattice = lattice
            else: 
                msg.err("The lattice vectors must be a 3x3 matrix with the vectors as rows "
                        "and the vectors must be linearly independent.")
        elif lattice is not None:
            msg.err("The lattice vectors must either be a string of 'sc', 'fcc', 'hcp', 'bcc', "
                    "or a 3x3 matrix with the vectors a rows.")
        else:
            self.lattice = None

    def _get_basis(self,basis):
        """Determines the atomic basis vectors for the system.
        
        Args:
            basis (list or None): A list of the atomic basis vectors or None,
                if None then the default basis is assumed.
        """
        if basis is not None and isinstance(basis,list):
            if len(basis[0]) == 3:
                self.basis = basis
            else:
                msg.err("The atomic basis must be a list of lists that is nx3 where n is "
                        "the number of atoms in the basis.")
        elif basis is None:
            if self.lattice is not None:
                if isinstance(self.lattice,string_types) and self.lattice.lower() != "hcp":
                    self.basis = [[0,0,0]]
                else:
                    self.basis = [[0,0,0],[0.5,0.28867513459,0.81649658093]]
            else:
                self.basis = [[0,0,0]]
        else:
            msg.err("The atomic basis must be a list of lists that is nx3 where n is "
                    "the number of atoms in the basis or left blank.")

    def _get_species(self,species):
        """Parses the atomic species as passed in by the user.
         
        Args:
            species (list): a list of the names of the atomic species present
                in the system.
        """
        if isinstance(species,list):
            self.knary = len(species)
            if self.knary == 1:
                self.species = [species[0],species[0]]
            else:
                self.species = species
        else:
            msg.err("The species must be a list of atomic elements.")
            

    def _get_concs(self,concs):
        """Parses the atomic concentrations passed in by the user.

        Args:
            concs (list): a list of the atomic concentrations.
        """
        if concs is not None and  len(concs) == len(self.species):
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

    def _get_arrows(self,arrows):
        """Parses the arrow concentrations passed in by the user.

        Args:
            arrows (list): a list of the arrow concentration on each element.
        """
        if arrows is not None and len(arrows) == len(self.species):
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
        
    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        if len(self.sequence) == 0:
            return len(self.atoms_paths) == self.nconfigs
        else:            
            return all(e.ready() for e in self.sequence.values())
    
    @property
    def euid_file(self):
        """Returns the full path to the euid file for this group.
        """
        return path.join(self.root,"euids.pkl")

    def _load_euids(self):
        """Loads the list ef `euid` from the `rset.pkl` file for this
        database group.
        """
        if self.euids is None:
            self.euids = self.load_pkl(self.euid_file)
        return self.euids

    @property
    def atoms_paths(self):
        """Returns a list of full paths to the folders that have `atoms.json` objects
        for the latest result set.
        """
        result = []
        for euid in self.euids:
            folder = self.index[euid]
            target = path.join(folder,"atoms.json")
            if path.isfile(target):
                result.append(folder)

        return result 
        
    def rset(self):
        """Returns a :class:`quippy.AtomsList`, one for each config in the
        latest result set.
        """
        from matdb.database.basic import atoms_from_json
        if len(self.sequence) == 0:
            #Return the configurations from this group; it is at the
            #bottom of the stack
            result = quippy.AtomsList()
            for epath in self.atoms_paths:
                result.append(atoms_from_json)
            return result
        else:
            result = []
            for e in self.sequence.values():
                result.extend(e.rset)
            return result
    
    def _build_lattice_file(self,target):
        """Creates the 'lattice.in' file that phenum needs in order to perform
        the enumeration.

        Args:
            target (str): relative path to where the 'lattice.in' file should be 
                saved.
        """

        from os import path
        target = path.join(target, "lattice.in")

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


    def _setup_configs(self, rerun=False):
        """Sets up the database structure for the enumeration code and creates
        the 'lattice.in' file. Also loops over the enumeration routine
        until the desired number of configurations have been reached
        (this is important for enumerations over small systems where
        the number of systems that are superperiodic is sinificant and
        so the number reported by polya is significantly larger than
        the actual number of unique configs).
        
        Args:
            group (:class:`matdb.database.basic.Group`): An instance of 
                the group class.

        """
        # We need to construct a lattice.in file then run phenum so that
        # each system gets the correct number of configurations.
        current = getcwd()
        chdir(self.root)
        dind = 0
        self._build_lattice_file(self.root)
        # Perform the enumeration, we allow for multiple attempts since the
        # number of configs returned the first time could be to small for
        # enumerations over small systems.
        recurse = 0
        while dind<self.nconfigs and recurse<3:
            dind = self._enumerate(dind)
            recurse += 1
        chdir(current)

        # Last of all, create the job file to execute the job array.
        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.euids,self.euid_file)

    def _enumerate(self,dind):
        """Performs the enumeration using phenum and creates the files in the
        correct folder for each system enumerated.

        Args:
            dind (int): The number of configs found so far.
        """
        from phenum.enumeration import _enum_out
        from phenum.makeStr import _make_structures
        from glob import glob
        
        _enum_out({"input":"enum.in","outfile":"enum.out",
                   "seed":self.rseed if self.rseed is None else self.rseed+dind,
                   "lattice":"lattice.in","distribution":["all",str(self.nconfigs-dind)],
                   "super":self.keep_supers,"sizes":None,"savedist":None,"filter":None,
                   "acceptrate":None})

        remove("enum.in")
        [remove(f) for f in listdir('.') if f.startswith("polya.")]
        # extract the POSCARS
        euids = _make_structures({"structures":None,"input":"enum.out",
                                  "species":self.species,"rattle":self.rattle,
                                  "mink":"t","outfile":"vasp.{}","displace":self.displace}
                                 ,return_euids=True)

        # Now we need to create the folder for each system we've enumerated
        if self.euids is None:
            self.euids = []
        from quippy.atoms import Atoms
        for count, dposcar in enumerate(glob("vasp.*")):
            if euids[count] not in self.euids:
                dind += 1
                datoms = Atoms(dposcar,format="POSCAR")
                self.create(datoms,cid=dind)
                self.index[euids[count]] = self.configs[dind]
                self.euids.append(euids[count])
        [remove(f) for f in listdir('.') if f.startswith("vasp.")]

        return dind
        
    def setup(self, rerun=False):
        """Enumerates the desired number of structures and setups up a folder
        for each one.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.

        """
        super(Enumerated, self).setup(self._setup_configs,rerun)
