'''Group of quippy.atoms.Atoms objects with displaced atomic positions from
a seed configuration.
'''
from matdb.atoms import Atoms, AtomsList
from matdb.database import Group
from os import path, mkdir
from matdb import msg
import numpy as np


class Distortions(Group):
    '''Distortions.py: Group to create from a seed configuration and distorts
    the atom positions or displaces the atoms randomly within a normal
    distribution of some standard deviation, std.
    Args:
        name (str): Distortions: 'Dist'
        rattle (float): the standard deviation of the normal distribution
            of atom deviations.
        ran_seed (hashable):(=1 default) seed for the normal distribution
            random number generator in ase.Atoms.rattle()
        volume_factor (float): the volume factor of the repeated cells
            (i.e 1==Same Cell Volume as atom_seed)
        cov_diag=(float): value along the diagonal of the 9x9 covariance
            matrix. Related to the standard deviaton of the lattice vector
            distortions
        min_index (int):(default=0) Default choice with the same ran_seed would
            produce the same choices of volume factor matrices.
        dbargs (dict): dictionary of arguments to be passed to the
            `Group` class.
    .. note:: Additional attributes are also exposed by the super class
            :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
             collection. This is also the name of the folder in which all of
             it's calculations will be performed.
        scaling_matrix(np.array): 3x3 matrix of determinant=volume factor.
        mean(tuple): 3x3 identity matrix, mean of the multivariate normal
             distribution of values in the scaling matrix.
        cov(list of lists): 9x9 diagonal covariance matrix with cav_diag along
             the diagonal.
        rattle (float): the amount to rattle the atoms in the config by.

    Returns:
        distortions (np.n darray): an array of atoms objects of length
             num_cells with distorted atom positions according to the normal
             distribution specified.
    '''
    def __init__(self, rattle=0, ran_seed=None, volume_factor=1.0,
                 cov_diag=1, min_index=0, name="Dist", **dbargs):
        self.name = name
        self.seeded = True
        dbargs['prefix'] = "D"
        dbargs['cls'] = Distortions
        if "Dist" not in dbargs['root']:
            new_root = path.join(dbargs['root'], "Dist")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Distortions, self).__init__(**dbargs)
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs[
                    "calculator"]:
                calcargs.update(dbargs["calculator"])
        self.ran_seed = (ran_seed if ran_seed is not None
                         else self.database.parent.ran_seed)
        self.rattle = float(rattle)
        self.volume_factor = volume_factor
        self.cov_diag = cov_diag
        self.min_index = min_index
        self.duids = None
        self._load_duids()

    def sub_dict(self):
        """Writes the attributes of this instance of the class to a dictionary.
        """
        dist_dict = {}
        dist_dict["rattle"] = self.rattle
        dist_dict["rseed"] = self.ran_seed
        dist_dict["volume_factor"] = self.volume_factor
        dist_dict["cov_diag"] = self.cov_diag
        dist_dict["min_index"] = self.min_index
        dist_dict["name"] = self.name
        return dist_dict

    @property
    def duid_file(self):
        """Returns the full path to the duid file for this group.
        """
        return path.join(self.root, "duids.pkl")

    def _load_duids(self):
        """Loads the list of `duid` from the `rset.pkl` file for this
        database group.
        """
        if self.duids is None:
            self.duids = self.load_pkl(self.duid_file)
        return self.duids

    @property
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group.
        """
        if len(self.sequence) == 0:
            return len(self.atoms_paths) == self.nconfigs
        else:
            result = []
            for g in self.sequence.values():
                result.extend(g.fitting_configs)
            return result

    @property
    def rset(self):
        """Returns a :class:`matdb.atoms.AtomsList`, one for each config in the
        latest result set.
        """
        if len(self.sequence) == 0:
            # Return the configurations from this group; it is at the
            # bottom of the stack
            result = AtomsList()
            for epath in self.fitting_configs:
                result.append(Atoms(path.join(epath, 'atoms.h5')))
            return result
        else:
            result = []
            for e in self.sequence.values():
                result.extend(e.rset())
            return result

    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            result = len(self.atoms_paths) == self.nconfigs
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

    @property
    def atoms_paths(self):
        """Returns a list of full paths to the folders that have `atoms.json` objects
        for the latest result set.
        """
        result = []
        for duid in self.duids:
            folder = self.index[duid]
            target = path.join(folder, "atoms.h5")
            if path.isfile(target):
                result.append(folder)
        return result

    def setup(self, rerun=False):
        """Sets up a copy of the atoms object for the current config.
        Args:
            rerun (bool): when True, recreate the folders even if they
             already exist.
        """
        super(Distortions, self).setup(self._setup_configs, rerun)
        if len(self.sequence) != 0:
            self.index = {}
            self.duids = []
            for seq in self.sequence.values():
                self.index.update(seq.index)
                self.duids.extend(seq.duids)
        self.save_index()
        self.save_pkl(self.duids, self.duid_file)

    def _setup_configs(self, rerun=False):
        """Loops over the distortions routine until the desired number of
        configurations have been reached
        Attributes:
            dists(): length nconfigs AtomsList
        Args:
            group (:class:`matdb.database.basic.Group`): An instance of the
                 group class.
        """
        from hashlib import sha1
        dists = self._get_distortions()

        if self.duids is None:
            self.duids = []
        for dist in dists:
            self.create(dist)
            chem_form = dist.get_chemical_formula(mode='reduce')
            duid = sha1(chem_form).hexdigest()
            self.duids.append(duid)
            self.index[duid] = self.configs[len(self.configs)]
        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.duids, self.duid_file)

    def _get_scaling_matrix(self):
        from numpy.random import multivariate_normal
        from numpy.linalg import det
        from numpy import reshape
        np.random.seed(self.ran_seed)
        mean, n = (1, 0, 0, 0, 1, 0, 0, 0, 1), self.cov_diag
        cov = [[n, 0, 0, 0, 0, 0, 0, 0, 0], [0, n, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, n, 0, 0, 0, 0, 0, 0], [0, 0, 0, n, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, n, 0, 0, 0, 0], [0, 0, 0, 0, 0, n, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, n, 0, 0], [0, 0, 0, 0, 0, 0, 0, n, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, n]]

        scaling_matrix = []
        for i in range(self.min_index, (self.nconfigs+self.min_index)):
            matrix = reshape(multivariate_normal(mean, cov), (3, 3))
            scaling_matrix.append(
                (self.volume_factor/det(matrix))**(1./3)*matrix)
        return scaling_matrix

    def _get_distortions(self):
        """Perform the duplication of the atom_seed and displacement of atom cells.
        Attributes:
            volume_factor (int): the volume factor of the repeated cells
                 (i.e 1==Same Cell Volume as atom_seed)
            cell_choice (ase.Atoms): each repeated atom_seed is rattled and
                 saved to to the distortions array.
        Returns:
            distortions (np.n darray): an array of atoms objects of length
                 num_cells with distorted atom positions according to the
                 normal distribution specified.
        """
        if(self.volume_factor != 1.0):
            scaling_matrix = self._get_scaling_matrix()
        atom_seed = AtomsList()
        for i in scaling_matrix:
            local_atoms = self.atoms.copy()
            local_atoms.set_cell(local_atoms.get_cell()*i)
            if (self.rattle != 0.0):
                local_atoms.rattle(stdev=self.rattle)
            atom_seed.append(local_atoms)
        return atom_seed
