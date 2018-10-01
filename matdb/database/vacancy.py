'''Group of matdb.atoms.Atoms objects with atomic vacancies
taken from a seed configuration.
'''

from os import path, mkdir
from itertools import islice, combinations, groupby

import numpy as np  # for use with arrays and random
from scipy.special import comb as choose  # find number of unique combinations

from matdb import msg
from matdb.atoms import Atoms, AtomsList
from matdb.database import Group  # create the vacancies group

class Vacancy(Group):
    '''Vacancy.py: Group to create atomic vacancies from a seed configuration.

    Args:
        name(str): default name Vacancy
        ran_seed (hashable):(=1 default) seed for the random number generator
             for index of vacancies selection.
        vac_per_atom (int < 1): The number of vacancies to include per
             atom in the cell. (i.e. 0.1 would be 1 in every 10 atoms.)
        min_index (int):(default=0) Default choice with the same ran_seed would
             produce the same vacancies in each cell.
        dbargs (dict): dictionary of arguments to be passed to the
            `Group` class.
    .. note:: Additional attributes are also exposed by the super class
          :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
             collection. This is also the name of the folder in which all of
             its calculations will be performed.
        num_atom(int): The number of atoms present in each atoms object.
        num_vac(int): The number of vacancies per cell.
        select_atoms(list): list of lists with indices of atoms to be removed
        unique_perm(int): number of possible combinations

    Returns:
        vacancies(AtomsList): list of atoms objects of length nconfigs with
             unique vacancies for each cell.
    '''

    def __init__(self, ran_seed=None, vac_per_atom=0, min_index=0,
                 name="Vacancy", **dbargs):
        self.name = name
        self.seeded = True
        dbargs['prefix'] = "V"
        dbargs['cls'] = Vacancy
        if "Vacancy" not in dbargs['root']:
            new_root = path.join(dbargs['root'], "Vacancy")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Vacancy, self).__init__(**dbargs)
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs[
                    "calculator"]:
                calcargs.update(dbargs["calculator"])
        self.ran_seed = (ran_seed if ran_seed is not None
                         else self.database.parent.ran_seed)
        self.vac_per_atom = vac_per_atom
        self.min_index = min_index
        self.vuids = None
        self._load_vuids()

    def sub_dict(self):
        """Writes the attributes of this instance of the class to a dictionary.
        """
        vac_dict = {}
        vac_dict["vac_per_atom"] = self.vac_per_atom
        vac_dict["rseed"] = self.ran_seed
        vac_dict["min_index"] = self.min_index
        vac_dict["name"] = self.name
        return vac_dict

    @property
    def vuid_file(self):
        """Returns the full path to the vuid file for this group.
        """
        return path.join(self.root, "vuids.pkl")

    def _load_vuids(self):
        """Loads the list of `vuid` from the `rset.pkl` file for this
        database group.
        """
        if self.vuids is None:
            self.vuids = self.load_pkl(self.vuid_file)
        return self.vuids

    def _get_random_choice(self, select_atoms, num_atoms, num_vac):
        """choose nconfigs lists of unique indices to remove from each atoms
        object config.

        Args:
            select_atoms(list of lists): the selected indices to be removed
                 from each config.
            num_atoms(int): number of atoms present in the cell.
            num_vac(int): number of vacancies to include in each config.
            extension_index(list of lists): selected_atoms for min_index>0.
        Returns:
            select_atoms(list of lists): the indices to remove from each
                 config.
        """
        # Choose random combination options until the length is reached.
        extension_index = []
        while(len(select_atoms) != self.nconfigs+self.min_index):
            # choose randomly from the possible combinations
            choice = sorted(list(np.random.choice(range(num_atoms),
                                                  size=num_vac,
                                                  replace=False)))
            # add the list choice to the selected atoms list of lists
            select_atoms.append(choice)
            # check for duplicates in the list and remove them.
            select_atoms = list(select_atoms for select_atoms, _
                                in groupby(sorted(select_atoms)))
            # for each d.b. only add unique choices.
            if(self.min_index > 0 and len(select_atoms) >
               (self.min_index+len(extension_index))):
                extension_index.append(choice)
        # choose the extension_index containing only unused values.
        if(self.min_index > 0):
            select_atoms = extension_index
        return select_atoms

    def _get_combinations(self, select_atoms, num_atoms, num_vac):
        '''This Approach allows for simple, efficient random iteration of all
        possible vacancies for small cell sizes. Limiting this approach that
        n choose k is less than 1000.

        Args:
            select_atoms(list of lists): the selected indices to be removed
                 from each config.
            num_atoms(int): number of atoms present in the cell.
            num_vac(int): number of vacancies to include in each config.
        Returns:
            select_atoms(list of lists): the indices to remove from each
                 config.
        '''
        atomic_vacancies = range(choose(num_atoms, num_vac, exact=True))
        np.random.shuffle(atomic_vacancies)  # shuffle all possible options
        atomic_vacancies = list(islice(atomic_vacancies, self.min_index,
                                       (self.min_index+self.nconfigs)))
        for i in atomic_vacancies:
            select_atoms.append(list(islice(combinations(range(num_atoms),
                                                         num_vac), i, i+1)))
        return select_atoms

    @property
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group.
        """
        if len(self.sequence) == 0:
            return len(self.atoms_paths()) == self.nconfigs
        else:
            result = []
            for g in self.sequence.values():
                result.extend(g.fitting_configs)
            return result

    def ready(self):
        """Returns True if this database has finished its computations
        and is ready to be used.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            result = len(self.atoms_paths()) == self.nconfigs
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

    def atoms_paths(self):
        """Returns a list of full paths to the folders that have
        `atoms.hdf5` objects for the latest result set.
        """
        result = []
        for vuid in self.vuids:
            folder = self.index[vuid]
            target = path.join(folder, "pre_comp_atoms.h5")
            if path.isfile(target):
                result.append(folder)
        return result

    def rset(self):
        """Returns a :class:`matdb.atoms.AtomsList`, one for each config in the
        latest result set.
        """
        if len(self.sequence) == 0:
            # Return the configurations from this group; it is at the
            # bottom of the stack
            result = AtomsList()
            for epath in self.atoms_paths():
                result.append(Atoms(path.join(epath, 'pre_comp_atoms.h5')))
            return result
        else:
            result = []
            for e in self.sequence.values():
                result.extend(e.rset())
            return result

    def setup(self, rerun=0):
        """Sets up a copy of the atoms object for the current config.

        Args:
            rerun (int): when > 0, recreate job files; if > 1, recreate the
                folders even if they already exist.
        """
        super(Vacancy, self).setup(self._setup_configs, rerun)
        if len(self.sequence) != 0:
            self.index = {}
            self.vuids = []
            for seq in self.sequence.values():
                self.index.update(seq.index)
                self.vuids.extend(seq.vuids)
        self.save_index()
        self.save_pkl(self.vuids, self.vuid_file)

    def _setup_configs(self, rerun=0):
        """Loops over the vacancies routine until the desired number of
        configurations have been reached
        Args:
            group (:class:`matdb.database.basic.Group`): An instance of
                the group class.
            rerun (int): when > 0, recreate job files; if > 1, recreate the
                folders even if they already exist.
        """
        from hashlib import sha1
        vacs, indices = self._get_vacancies()
        if self.vuids is None:
            self.vuids = []
        if(not self.is_setup() or rerun > 1):
            for vac in vacs:
                self.create(vac)
                vuid = sha1(str(indices[len(self.vuids)])).hexdigest()
                self.vuids.append(vuid)
                self.index[vuid] = self.configs[len(self.configs)]
        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.vuids, self.vuid_file)

    def _get_vacancies(self):
        '''Vacancies.py: Group to create atomic vacancies from a seed configuration.

        Args:
            atom_seed (list, str, matdb.atoms.Atoms): The location of the
                 files that will be read into to make the atoms object or an
                 atoms object.
            ran_seed (hashable):(=1 default) seed for the random number
                 generator for index of vacancies selection.
            nconfigs (int): number of cells with vacancies to create.
            vac_per_atom (int < 1): The number of vacancies to include per
                 atom in the cell. (i.e. 0.1 would be 1 in every 10 atoms.)
            min_index (int):(default=0) Default choice with the same ran_seed
                 would produce the same vacancies in each cell.

        .. note:: Additional attributes are also exposed by the super class
              :class:`Group`.

        Attributes:
            name (str): name of this database type relative to the over
                 database collection. This is also the name of the folder
                 in which all of its calculations will be performed.
            num_atom(int): The number of atoms present in each atoms object.
            num_vac(int): The number of vacancies per cell.
            seed_state(tuple, len=4): values 1,3-4 are set by ran_seed after
                 the first call to np.random and do not change, value 2 gives
                 the ith value of a call to random
            select_atoms(list): list of lists with indices of atoms to be
                 removed
            unique_perm(int): number of possible combinations
        Returns:
            vacancies(AtomsList): an list of atoms objects of length nconfigs
                 with unique vacancies for each cell.
        '''
        select_atoms = []  # list of lists with indices of atoms to be removed
        num_atoms = int(len(self.atoms.get_positions()))  # number of atoms
        num_vac = int(num_atoms * self.vac_per_atom)

        np.random.seed(self.ran_seed)  # Set the random seed for reproduction
        if(choose(num_atoms, num_vac) > 1000):
            select_atoms = self._get_random_choice(select_atoms, num_atoms,
                                                   num_vac)
        else:
            select_atoms = self._get_combinations(select_atoms, num_atoms,
                                                  num_vac)
        atom_seed = AtomsList()
        for i in select_atoms:
            local_atoms = self.atoms.copy()
            del local_atoms[i]
            atom_seed.append(local_atoms)
        return atom_seed, select_atoms
