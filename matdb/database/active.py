"""Group of configurations that is created from an enumerated list of structures.
"""
from matdb.database import Group
from matdb import msg
from os import path, getcwd, chdir, remove, listdir, mkdir
import numpy as np
from six import string_types
from glob import glob
from matdb.atoms import AtomsList, Atoms

class Active(Group):
    """Sets up the calculations for a set of configurations that are being
    added to by the active learning approach.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        auids (list): the unique ids for each config in the group
    """
    def __init__(self, name="active", **dbargs):
        self.name = name
        dbargs['prefix'] = "Ac"
        dbargs['cls'] = Active
        dbargs['trainable'] = True
        if "calculator" in dbargs:
            self.calcargs = dbargs["calculator"].copy()
            del dbargs["calculator"]
        if "Active" not in dbargs['root']:
            new_root =path.join(dbargs['root'],"Active")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Active, self).__init__(**dbargs)

        self.auids = None
        self._load_auids()
        self.nconfigs = len(self.auids) if self.auids is not None else 0
        self.last_iteration = None
        cur_iter = len(glob("iter_*.pkl"))
        self.iter_file = path.join(self.root,"iter_{}.pkl".format(cur_iter))
        self._load_last_iter()        
        
    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        return len(self.fitting_configs) == self.nconfigs
    
    @property
    def auid_file(self):
        """Returns the full path to the euid file for this group.
        """
        return path.join(self.root,"auids.pkl")

    def _load_auids(self):
        """Loads the list of `auid` from the `rset.pkl` file for this
        database group.
        """
        if self.auids is None:
            self.auids = self.load_pkl(self.auid_file)

    def _load_last_iter(self):
        """Loads the list of paths from the `iter_{}.pkl` file for this
        database group's last iteration.
        """
        if self.last_iteration is None:
            self.last_iteration = self.load_pkl(self.iter_file)
            
    @property
    def fitting_configs(self):
        """Returns a list of full paths to the folders that have `atoms.json` objects
        for the full result set.
        """
        result = []
        for auid in self.auids:
            folder = self.index[str(auid)]
            target = path.join(folder,"atoms.h5")
            if path.isfile(target):
                result.append(folder)

        return result 
        
    def rset(self):
        """Returns a :class:`matdb.atoms.AtomsList`, one for each config in the
        latest result set.
        """

        #Return the configurations from this group; it is at the
        #bottom of the stack
        result = AtomsList()
        for epath in self.fitting_configs:
            result.append(Atoms(epath))
        return result

    def add_configs(self,new_configs,iteration):
        """Adds the atoms objects in the list to the configs of the active set.

        Args:
            new_configs (list): list of `matdb.atoms.Atoms` objects to be added
                to the active learning set.
        """

        self.new_configs = new_configs
        self.nconfigs += len(new_configs)
        self.iter_file = path.join(self.root,"iter_{}.pkl".format(iteration))
        self.last_iteration = {}

    def _setup_configs(self, rerun=False):
        """Sets up the database structure for the active set and creates a
        folder for the `calculator` to run in for each config.        
        """
        # We need to make sure that none of the configs in this round
        # of the active learning have already been visited by
        # constructing their auids and then verifying that the auid
        # hasn't been visited before.
        dind = len(self.auids)
        iter_ind = 0
        from hashlib import sha1 
        for config in self.new_configs:
            auid = sha1(''.join(tuple(tuple([tuple(i) for i in config.cell]),
                        tuple([tuple(i) for i in config.positions]),
                        tuple(config.get_chemical_symbols()))).encode('utf-8'))
            if auid in self.auids:
                self.nconfigs -= 1
                continue
            
            dind += 1
            self.create(config,cid=dind)
            self.index[auid] = self.configs[dind]
            self.last_iteration[iter_ind] = self.configs[dind]
            self.auids.append(auid)
            iter_ind += 1
            
        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.auids,self.auid_file)
        self.save_pkl(self.last_iteration,self.iter_file)
        
    def setup(self, rerun=False):
        """Enumerates the desired number of structures and setups up a folder
        for each one.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.

        """
        super(Active, self).setup(self._setup_configs,rerun)           
