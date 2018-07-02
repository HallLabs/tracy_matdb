"""Group of configurations selected from the prototypes database.
"""

from os import path, getcwd, chdir, remove, listdir, mkdir
import numpy as np
from six import string_types
from glob import glob
from copy import deepcopy
from hashlib import sha1

from phenum.element_data import get_lattice_parameter

from matdb.database import Group
from matdb import msg
from matdb.atoms import Atoms, AtomsList
from matdb.utility import _get_reporoot, chdir

class Prototypes(Group):
    """Constructs a selection of configurations based off the AFLOW
    prototype structures list.

    Args:
        atoms (matdb.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.Database): parent sequence to which this database
          belongs. Could also be another :class:`Hessian`.
        structures (dict): a dictionary stating either the number or the file names
          of the unary, binary, and ternary sized systems.
        ran_seed (int): a random seed to use that would override the global seed.
        permutations (dict): A dictionary of lists that would specify the
          atomic permutations to restrict the configurations to. For
          example to only use A:B  and no B:A permutations the
          dict would be {"binary": [["A","B"]]}, where A and B were
          replaced with the correct atomic species.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
	<<Additional attributes your database group will have>>.

    """
    def __init__(self, name="prototype", structures=None, ran_seed=None, permutations=None,
                 **dbargs):
        self.name = name
        self.seeded = False
        dbargs["prefix"] = "P"
        dbargs["cls"] = Prototypes
        if "Prototypes" not in dbargs['root']:
            from os import mkdir
            new_root =path.join(dbargs['root'],"Prototypes")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Prototypes, self).__init__(**dbargs)

        self.in_structures = structures
        self.ran_seed = ran_seed
        self.permutations = permutations
        self.species = self.database.parent.species

        #Make sure that we override the global calculator default values with
        #those settings that we know are needed for good phonon calculations.
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs["calculator"]:
                calcargs.update(dbargs["calculator"])
                dbargs["calculator"] = calcargs

        # The prototypes are saved into the file prototypes.tar.gz, if
        # this is the first time prototypes has been run we need to unpack it.
        template_root = path.join(_get_reporoot(), "matdb", "templates")
        if not path.isdir(path.join(template_root, "uniqueUnaries")):
            import tarfile
            with chdir(template_root):
                tarf = "prototypes.tar.gz"
                tar = tarfile.open(tarf, "r:gz")
                tar.extractall()
                tar.close()

        # parse the structures to make a list of paths to the source folders for the
        if self.ran_seed is not None:
            import random
            random.seed(self.ran_seed)

        self.puuids = None
        self._load_puuids()
        self.nconfigs = 0

        self.structures = {}
        for k,v in structures.items():
            if k.lower() == "unary":
                cand_path = path.join(template_root, "uniqueUnaries")
            elif k.lower() == "binary":
                cand_path = path.join(template_root, "uniqueBinaries")
            elif k.lower() == "ternary":
                cand_path = path.join(template_root, "uniqueTernaries")
            else: # pragma: no cover
                msg.warn("Must specify the system size, i.e., unary, binary, or "
                         "ternary. {} not recognized".format(k))
                continue
            if isinstance(v,list):
                self.structures[k.lower()] = []
                for prot in v:
                    files = glob("{0}/*{1}*".format(cand_path,prot))
                    if len(files) < 1: # pragma: no cover
                        msg.warn("No prototypes of size {0} matched the string "
                                 "{1}".format(k, prot))
                    else:
                        self.structures[k.lower()].extend(files)
            elif isinstance(v, str) and v == "all":
                files = glob("{0}/*".format(cand_path))
                self.structures[k.lower()] = files
            elif isinstance(v, int):
                from random import shuffle
                files = glob("{0}/*".format(cand_path))
                shuffle(files)
                keep = files[:v]
                self.structures[k.lower()] = keep
            else: #pragma: no cover
                msg.err("Couldn't parse {0} structures for {1} case. Must be either "
                        "a list of file names, 'all', or an int.".format(v, k))

            if self.permutations is not None and  k.lower() in self.permutations.keys():
                self.nconfigs += len(self.structures[k.lower()])*len(self.permutations[k.lower()])
            else:
                if k.lower() == "unary":
                    self.nconfigs += len(self.structures[k.lower()])*3
                elif k.lower() == "binary" or k.lower() == "ternary":
                    self.nconfigs += len(self.structures[k.lower()])*6
                else: #pragma: no cover
                    continue

    @property
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group.
        """
        configs = AtomsList()
        if len(self.sequence) == 0:
            for config in self.config_atoms.values():
                configs.append(config)
        else:
            for seq in self.sequence.values():
                configs.extend(seq.fitting_configs)

        return configs

    def sub_dict(self):
        """Returns a dict needed to initialize the class.
        """
        args = {"structures": self.in_structures, "ran_seed": self.ran_seed,
                "permutations": self.permutations, "prefix": self.prefix, "name": self.name}
        return args

    @property
    def rset(self):
        """Returns all the configurations that are in this database.

        Returns:
            list: of :class:`matdb.atoms.Atoms`
        """
        if len(self.sequence) == 0:
            #We are at the bottom of the stack;
            return self.fitting_configs
        else:
            #Check where we are in the stack. If we are just below the database,
            #then we want to return <<your description of the rset here>>
	    #If we are not, then we must a parameter grid of sequences
            #to select from.
	    return [res for p in self.sequence.values() for res in p.rset]

    def ready(self):
        """Returns True if all the calculations have been completed.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            if len(self.configs) >= 1:
                result = True
                for config in self.configs.values():
                    if not path.isfile(path.join(config, "atoms.h5")):
                        result = False
                        break
            else:
                result = False
            if not result:
                msg.std("{} is not ready. Exiting.".format(self.root), 2)
            return result
        else:
            ready = False
            for p in self.sequence.values():
                if not p.ready():
                    msg.std("{} is not ready. Exiting.".format(p.root), 2)
                    break
            else:
                ready = True
            return ready

    def setup(self, rerun=False):
        """Creates all the configurations for calculation from the prototype
        databases.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.
        """
        super(Prototypes, self).setup(self._setup_configs, rerun)

    def _setup_configs(self, rerun):
        """Loops over the choosen prototype structures and possible
        occupations for those structures, i.e., A:B and B:A, to
        generate all the configurations needed for th calculation.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.
        """
        #We also don't want to setup again if we have the results already.
        if self.ready() and not rerun:
            return

        if self.puuids is None:
            self.puuids = []
        if not self.is_setup():
            from itertools import product
            # Loop over the sizes in the saved structures
            for k, v in self.structures.items():
                perms = self._get_perms(k)
                for fpath, perm in product(v,perms):
                    hash_str = "{0}-{1}".format(fpath.split("/")[-1],"".join(perm))
                    hash_str = str(sha1(hash_str).hexdigest())
                    if hash_str not in self.puuids:
                        self.puuids.append(hash_str)
                        self._correct_poscar(fpath, path.join(self.root, "POSCAR"), perm)
                        datoms = Atoms(path.join(self.root, "POSCAR"),format="vasp")
                        self.create(datoms)
                        remove(path.join(self.root, "POSCAR"))

        # Last of all, create the job file to execute the job array.
        self.save_pkl(self.puuids,self.puuid_file)
        self.jobfile(rerun)

    def _get_perms(self, size):
        """Finds the allowed permutations of the atomic species that need to
        be permuted over.

        Args:
            size (str): Any of "unary", "binary", or "ternary".
        """

        from itertools import permutations
        res = None
        if size == "unary":
            r = 1
            res = self.permutations["unary"] if "unary" in self.permutations.keys() else None
        elif size == "binary":
            r = 2
            res = self.permutations["binary"] if "binary" in self.permutations.keys() else None
        elif size == "ternary":
            r = 3
            res = self.permutations["ternary"] if "ternary" in self.permutations.keys() else None
        else: # pragma: no cover
            msg.err("{} is not a valid prototype size.".format(size))
            return None

        possible_perms = [list(i) for i in permutations(self.species, r=r)]
        perms_copy = deepcopy(possible_perms)
        if res is not None:
            for i in perms_copy:
                if i not in res:
                    possible_perms.remove(i)

        return possible_perms

    def _correct_poscar(self, source, target, species):
        """Corrects the POSCAR so that it has the correct lattice parameter
        and title string for the system to be read into ASE.

        Args:
            source (str): the path to the prototype POSCAR.
            target (str): the path to the output POSCAR.
            species (list): a list of the species to be used for this POSCAR.
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

        # fix the title.
        f_lines[0] = "{0} : {1}".format(" ".join(species), f_lines[0])

        # fix the lattice parameter
        lat_param, ttl = get_lattice_parameter(species, concs, lat_vecs, sum(concs), " ")
        f_lines[1] = "{} \n".format(lat_param)

        with open(target, "w+") as f:
            for line in f_lines:
                f.write(line)

    def _load_puuids(self):
        """Loads the list of `euid` from the `rset.pkl` file for this
        database group.
        """
        if self.puuids is None:
            self.puuids = self.load_pkl(self.puuid_file)
        return self.puuids

    @property
    def puuid_file(self):
        """Returns the full path to the euid file for this group.
        """
        return path.join(self.root,"puuids.pkl")
