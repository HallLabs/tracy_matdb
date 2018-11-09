'''Group of matdb.atoms.Atoms objects with changed stoichiometries for atomic
positions from a seed configuration.
'''
from os import path, mkdir
import random

import numpy as np

from matdb import msg
from matdb.atoms import AtomsList, Atoms
from matdb.database import Group

class Substitution(Group):
    """Substitution.py: A Group to create substitutions in the stoichiometry from
    a seed configuration.

    Args:
        name(str): Default name Substitution
        stoich (list of lists): each list contains the decimal concentration
             of each element in the system where followed by decimal fraction
             of the number of configs which follow this stoichiometry.
             The decimal concentration in each list as well as the decimal
             fraction of nconfigs across all lists must sum to 1.
        ran_seed (int):seed for the random number generator. To allow for
             reproducibility and extension of databases.
        min_index (int):(default=0) For extension of databases set a nonzero
             starting index for extension calculations to avoid repeat values.
        dbargs (dict): dictionary of arguments to be passed to the
            `Group` class.
    .. note:: Additional attributes are also exposed by the super class
            :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
            collection. This is also the name of the folder in which
            all of its calculations will be performed.
        num_atom(int): The number of atoms present in each atoms object.

    Returns:
        vacancies(AtomsList): list of atoms objects of length nconfigs with
        unique vacancies for each cell.
    """
    def __init__(self, ran_seed=None, stoich=None, min_index=0,
                 name="Substitution", **dbargs):
        self.name = name
        self.seeded = True
        dbargs['prefix'] = "S"
        dbargs['cls'] = Substitution
        if "Substitution" not in dbargs['root']:
            new_root = path.join(dbargs['root'], "Substitution")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Substitution, self).__init__(**dbargs)
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs[
                    "calculator"]:
                calcargs.update(dbargs["calculator"])
        self.ran_seed = (ran_seed if ran_seed is not None
                         else self.database.parent.ran_seed)
        self.stoich = stoich
        self.min_index = min_index
        self.suids = None
        self._load_suids()

    def sub_dict(self):
        """Writes the attributes of this instance of the class to a dictionary.
        """
        sub_dict = {}
        sub_dict["stoich"] = self.stoich
        sub_dict["rseed"] = self.ran_seed
        sub_dict["min_index"] = self.min_index
        sub_dict["name"] = self.name
        return sub_dict

    @property
    def suid_file(self):
        """Returns the full path to the suid file for this group.
        """
        return path.join(self.root, "suids.pkl")

    def _load_suids(self):
        """Loads the list of `suid` from the `rset.pkl` file for this
        database group.
        """
        if self.suids is None:
            self.suids = self.load_pkl(self.suid_file)
        return self.suids

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
        """Returns a list of full paths to the folders that have `atoms.h5`
        objects for the latest result set.
        """
        result = []
        for suid in self.suids:
            folder = self.index[suid]
            target = path.join(folder, "pre_comp_atoms.h5")
            if path.isfile(target):
                result.append(folder)
        return result

    def _set_stoichiometry(self):
        '''Get a list of length nconfigs by num_atoms in atoms with random
        unique stoichiometries according to the stoich specifications.

        Attributes:
            f(factorial function): Numpy's factorial function.
            num_atoms(int): The number of atoms present in the seed object
                passed to the group.
            elements(rp.array): The elements present in the seed config in
                alphabetical order.
            num_ele(int): The number of elements
            combs(set): The unique combinations chosen by random.
            repeats(int): The total number of possible unique permutations is
                num_atoms!/repeats! where repeats is equal to repeat * (number
                of times element was repeated)! for each element in atoms.

        Returns:
            combs(list): a indexed list of each random unique choice of
                stoichiometry, length nconfigs.
        '''
        f = np.math.factorial
        # returns the numbers of atoms present in the cell.
        num_atoms = self.atoms.get_number_of_atoms()

        # get the elements in the object in alphabetical order.
        elements = np.unique(sorted(self.atoms.get_chemical_symbols()))
        num_ele = len(elements)
        # The percent of configs arguments in self.stoich must sum to 1.
        # The stoichiometry for each stoichiometry choice sums to 1.
        if (1.0 != sum(list(map(list, zip(*self.stoich)))[num_ele]) or
            len(self.stoich) !=
            sum([True if 1.0 == sum(self.stoich[i][:(num_ele)])
                 else False for i in range(len(self.stoich))])):
            raise ValueError('The stoichiometry arguments must sum to 1')

        combs, old_combs, tot_index = set(), set(), 0
        # for each choice of stoichiometry.
        for j in range(len(self.stoich)):
            # Percentage of nconfigs chosen for the chosen stoichiometry.
            max_index = _round_to_even(self.stoich[j][num_ele]*self.nconfigs)
            tot_index += max_index + self.min_index
            # Reset atomic number to an empty list.

            atomic_numbers, repeats = [], 1
            # for each element present in atoms
            for i in range(num_ele):
                # Get the number of rep. permutations in permutations(atoms)
                repeats = repeats * f(_round_to_even(num_atoms *
                                                     self.stoich[j][i]))
                # Create a sorted list for this choice of stoichiometry.
                atomic_numbers += ([elements[i]] * _round_to_even(
                    num_atoms * self.stoich[j][i]))

            # If stoich asks for more permutations than exist raise an error.
            if(max_index + self.min_index <= float(f(num_atoms))/(repeats)):
                # Until random.shuffle has created the correct number of sets.
                while len(combs) < tot_index:
                    # Shuffle the atomic numbers list.
                    random.shuffle(atomic_numbers)
                    # Add shuffled list and check for uniqueness.
                    combs.add(tuple(atomic_numbers))
                    if len(combs) < (tot_index + self.min_index - max_index
                                     and self.min_index != 0):
                        old_combs.add(tuple(atomic_numbers))
            else:
                raise ValueError('There are not enough unique configs')

        combs = list(combs.difference(old_combs))
        return combs

    def setup(self, rerun=0):
        """Sets up a copy of the atoms object for the current config.

        Args:
            rerun (int): when > 0, recreate job files; if > 1, recreate the
                folders even if they already exist.
        """
        super(Substitution, self).setup(self._setup_configs, rerun)
        if len(self.sequence) != 0:
            self.index = {}
            self.suids = []
            for seq in self.sequence.values():
                self.index.update(seq.index)
                self.suids.extend(seq.suids)
        self.save_index()
        self.save_pkl(self.suids, self.suid_file)

    def _setup_configs(self, rerun=0):
        """Loops over the substitutions routine until the desired number of
        configurations have been reached

        Args:
            group (:class:`matdb.database.basic.Group`): An instance of
                the group class.
            rerun (int): when > 0, recreate job files; if > 1, recreate the
                folders even if they already exist.
        """
        from hashlib import sha1
        subs = self._get_substitution()
        if self.suids is None:
            self.suids = []
        if(not self.is_setup() or rerun > 1):
            for sub in subs:
                self.create(sub)
                chem_form = sub.get_chemical_formula(mode='reduce')
                suid = sha1(chem_form.encode()).hexdigest()
                self.suids.append(suid)
                self.index[suid] = self.configs[len(self.configs)]
        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.suids, self.suid_file)

    def _get_substitution(self):
        '''
        '''
        np.random.seed(self.ran_seed)  # Set the seed for reproducibility.
        combs = self._set_stoichiometry()
        seed_atoms = AtomsList()

        for i in combs:
            local_atoms = self.atoms.copy()
            local_atoms.set_chemical_symbols(i)
            seed_atoms.append(local_atoms)
        return seed_atoms


def _round_to_even(n):
        '''While unit tests for python --version < 3.2 the built-in round
        function rounds .5 away from zero. The stoichiometry dict relies on
        round to even to compute correct values.
        '''
        if n % 1 == 0.5 and int(n) % 2 == 0:
            return int(n)
        else:
            return int(round(n))
