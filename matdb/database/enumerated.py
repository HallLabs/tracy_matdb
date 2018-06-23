"""Group of configurations that is created from an enumerated list of structures.
"""
from matdb.database import Group
from matdb import msg
from os import path, getcwd, chdir, remove, listdir, mkdir
import numpy as np
from six import string_types
from matdb.atoms import Atoms, AtomsList
from matdb.utility import copyonce
from glob import glob
from jinja2 import Environment, PackageLoader

class Enumerated(Group):
    """Sets up the calculations for a random sampling of structures from
    an enumerated list.

    Args:
        sizes (list): a list containing the smallest and larges cell sizes to
            include in the database.
        lattice (list or str): either the name of the lattice to use
            ('sc','fcc','bcc', 'hcp') or the atomic basis vectors to use.
        basis (list, optional): the atomic basis to use for the enumeration.
            Defaults to [0,0,0] the origin.
        concs (list, optional): the concentrations of each atomic species.
            Defaults to None.
        arrows (list, optional): the maximum number of atoms to displace.
            Defaults ot None.
        eps (float, optional): floating point tolerance for comparisons.
            Defaults to 1E-3.
        ran_seed (int or float, optional): a seed to feed to the random number generator.
            Defaults to None.
        rattle (float, optional): the amount to rattle the atoms by. Defaults to 0.0.
        keep_supers (bool, optional): True if the superperiodic cells are to be kept
            in the enumerated list. Defaults to False.
        displace (float, optional): the amount to displace atoms with arrows. Defaults
            to 0.0.

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
                 arrows = None, eps=None, name="enum", rattle=None, ran_seed=None,
                 keep_supers=None, displace=None, **dbargs):
        self.name = name
        dbargs['prefix'] = "E"
        dbargs['cls'] = Enumerated
        if "Enum" not in dbargs['root']:
            new_root =path.join(dbargs['root'],"Enum")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Enumerated, self).__init__(**dbargs)

        if eps is None:
            self.eps = 10**(-3)
        else:
            self.eps = eps

        self.ran_seed = ran_seed if ran_seed is not None else self.database.parent.ran_seed
        self._get_lattice(lattice)
        self._get_basis(basis)
        self.species = self.database.parent.species
        self.knary = len(self.species)
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
        else:
            raise ValueError("The sizes specified must be a list of 1 or 2 values."
                             "If one value then it must be the largest cell size to include,"
                             "i.e., [32]. If 2 values then it the first value is the smallest "
                             "and the second value is the largest cell size to include, "
                             "i.e., [10,12].")
        self.euids = None
        self._load_euids()

    def sub_dict(self):
        """Writes the attributes of this instance of the class to a dictionary.
        """
        enum_dict = {}
        enum_dict["sizes"] = [self.min_size,self.max_size]
        enum_dict["basis"] = self.basis
        enum_dict["lattice"] = self.lattice
        enum_dict["concs"] = self.concs
        enum_dict["arrows"] = self.arrows
        enum_dict["eps"] = self.eps
        enum_dict["name"] = self.name
        enum_dict["rattle"] = self.rattle
        enum_dict["ran_seed"] = self.ran_seed
        enum_dict["keep_supers"] = self.keep_supers
        enum_dict["displace"] = self.displace
        return enum_dict

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
                self.lattice_name = "fcc"
            elif lattice.lower() == "bcc":
                self.lattice = [[0.5,0.5,-0.5],[0.5,-0.5,0.5],[-0.5,0.5,0.5]]
                self.lattice_name = "bcc"
            elif lattice.lower() == "sc":
                self.lattice = [[1,0,0],[0,1,0],[0,0,1]]
                self.lattice_name = "sc"
            elif lattice.lower() == "hcp":
                self.lattice = [[1,0,0],[.5, 0.866025403784439, 0],[0, 0, 1.6329931618554521]]
                self.lattice_name = "hcp"
            else: #pragma: no cover
                msg.err("The lattice type {} is unsupported. Please enter your lattice vectors "
                        "as a 3x3 matrix with the vectors as rows in the config file "
                        "(i.e. [a1,a2,a3]).".format(lattice))
        elif isinstance(lattice,list):
            if len(lattice) == 3 and len(lattice[0]) == 3 and np.linalg.det(lattice) != 0:
                self.lattice = lattice
                self.lattice_name = "custom"
            else: #pragma: no cover
                msg.err("The lattice vectors must be a 3x3 matrix with the vectors as rows "
                        "and the vectors must be linearly independent.")
        elif lattice is not None: #pragma: no cover
            msg.err("The lattice vectors must either be a string of 'sc', 'fcc', 'hcp', 'bcc', "
                    "or a 3x3 matrix with the vectors a rows.")
        else:
            self.lattice = None
            self.lattice_name = None

    def _get_basis(self,basis):
        """Determines the atomic basis vectors for the system.

        Args:
            basis (list or None): A list of the atomic basis vectors or None,
                if None then the default basis is assumed.
        """
        if basis is not None and isinstance(basis,list):
            if len(basis[0]) == 3:
                self.basis = basis
            else: #pragma: no cover
                msg.err("The atomic basis must be a list of lists that is nx3 where n is "
                        "the number of atoms in the basis.")
        elif basis is None:
            if self.lattice is not None:
                if self.lattice_name != "hcp":
                    self.basis = [[0,0,0]]
                else:
                    self.basis = [[0,0,0],[0.5,0.28867513459,0.81649658093]]
            else:
                self.basis = [[0,0,0]]
        else: #pragma: no cover
            msg.err("The atomic basis must be a list of lists that is nx3 where n is "
                    "the number of atoms in the basis or left blank.")

    def _get_concs(self,concs):
        """Parses the atomic concentrations passed in by the user.

        Args:
            concs (list): a list of the atomic concentrations.
        """
        if concs is not None and  len(concs) == len(self.species):
            if len(concs[0]) == 3:
                self.concs = concs
                self.conc_res = True
            else: #pragma: no cover
                msg.err("The concetrations must be a nx3 list.")
        elif concs is None:
            self.concs = concs
            self.conc_res = False
        else: #pragma: no cover
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
                self.arrow_res = True
            else: #pragma: no cover
                msg.err("The arrows must be a list of values <= 1.")
        elif arrows is None:
            self.arrows = arrows
            self.arrow_res = False
        else: #pragma: no cover
            msg.err("The number of species and arrow concentrations must have the "
                    "same length.")

    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            return len(self.fitting_configs) == self.nconfigs
        else:
            for e in self.sequence.values():
                temp=e.ready()
            return all(e.ready() for e in self.sequence.values())

    @property
    def euid_file(self):
        """Returns the full path to the euid file for this group.
        """
        return path.join(self.root,"euids.pkl")

    def _load_euids(self):
        """Loads the list of `euid` from the `rset.pkl` file for this
        database group.
        """
        if self.euids is None:
            self.euids = self.load_pkl(self.euid_file)
        return self.euids

    @property
    def fitting_configs(self):
        """Returns a list of full paths to the folders that have `atoms.h5` objects
        for the latest result set.
        """
        result = []
        for euid in self.euids:
            folder = self.index[euid]
            target = path.join(folder,"atoms.h5")
            if path.isfile(target):
                result.append(folder)

        return result

    @property
    def rset(self):
        """Returns a :class:`matdb.atoms.AtomsList`, one for each config in the
        latest result set.
        """
        if len(self.sequence) == 0:
            #Return the configurations from this group; it is at the
            #bottom of the stack
            result = AtomsList()
            for epath in self.fitting_configs:
                result.append(Atoms(path.join(epath,"atoms.h5")))
            return result
        else:
            result = []
            for e in self.sequence.values():
                result.extend(e.rset())
            return result

    def _build_lattice_file(self,target):
        """Creates the 'lattice.in' file that phenum needs in order to perform
        the enumeration.

        Args:
            target (str): relative path to where the 'lattice.in' file should be
                saved.
        """

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
                    temp.append("{0} {1}".format(" ".join([str(j) for j in self.concs[i]]),a))
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
        self._build_lattice_file(self.root)
        dind = 0
        # Perform the enumeration, we allow for multiple attempts since the
        # number of configs returned the first time could be to small for
        # enumerations over small systems.
        chdir(self.root)
        recurse = 0
        while dind<self.nconfigs and recurse<5:
            dind = self._enumerate(dind,recurse,current)
            recurse += 1
        chdir(current)

        # Last of all, create the job file to execute the job array.
        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.euids,self.euid_file)

    def _enumerate(self,dind,recurse,home):
        """Performs the enumeration using phenum and creates the files in the
        correct folder for each system enumerated.

        Args:
            dind (int): The number of configs found so far.
            recurse (int): The number of times we've attempted to find
                a unique set of enumerations over the same range.
            home (str): The home directory.
        """

        _enum_out({"input":"enum.in","outfile":"enum.out",
                   "seed":self.ran_seed if self.ran_seed is None else self.ran_seed+dind+recurse,
                   "lattice":"lattice.in","distribution":["all",str(self.nconfigs-dind)],
                   "super":self.keep_supers,"sizes":None,"savedist":None,"filter":None,
                   "acceptrate":None})

        remove("enum.in")
        [remove(f) for f in listdir('.') if f.startswith("polya.")]
        # extract the POSCARS
        euids = _make_structures({"structures":None, "input":"enum.out",
                                  "species":self.species, "rattle":self.rattle,
                                  "mink":"t", "outfile":"vasp.{}", "displace":self.displace,
                                  "config":"f", "remove_zeros":"f"}, return_euids=True)

        # Now we need to create the folder for each system we've enumerated
        if self.euids is None:
            self.euids = []
        current = getcwd()
        for count, dposcar in enumerate(glob("vasp.*")):
            if euids[count] not in self.euids:
                dind += 1
                datoms = Atoms(dposcar,format="vasp")
                chdir(home)
                self.create(datoms,cid=dind)
                chdir(current)
                copyonce(dposcar,path.join(self.configs[dind],"POSCAR_orig"))
                self.index[str(euids[count].hexdigest())] = self.configs[dind]
                self.euids.append(str(euids[count].hexdigest()))
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

        if len(self.sequence) != 0:
            self.index = {}
            self.euids = []
            for seq in self.sequence.values():
                self.index.update(seq.index)
                self.euids.extend(seq.euids)
        self.save_index()
        self.save_pkl(self.euids,self.euid_file)
