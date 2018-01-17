"""Group of configurations that is created from an enumerated list of structures.
"""
from .basic import Group
from matdb import msg
from os import path, getcwd, chdir, remove, listdir, mkdir
import numpy as np
from six import string_types
import quippy

class Active(Group):
    """Sets up the calculations for a set of configurations that are being
    added to by the active learning approach.

    Args:
        new_configs (list): list of `quippy.atoms.Atoms` objects to be added
            to the active learning set.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        auids (list): the unique ids for each config in the group
    """
    def __init__(self, new_configs=None, name="active", **dbargs):
        self.name = name
        dbargs['prefix'] = "A"
        dbargs['cls'] = Active
        dbargs['trainable'] = True
        if "Active" not in dbargs['root']:
            new_root =path.join(dbargs['root'],"Active")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Active, self).__init__(**dbargs)

        self.new_configs = new_configs
        
        self.auids = None
        self._load_auids()
        self.nconfigs = len(self.auids) + len(new_configs)
        
    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        return len(self.atoms_paths) == self.nconfigs
    
    @property
    def auid_file(self):
        """Returns the full path to the euid file for this group.
        """
        return path.join(self.root,"auids.pkl")

    def _load_auids(self):
        """Loads the list of `euid` from the `rset.pkl` file for this
        database group.
        """
        if self.auids is None:
            self.auids = self.load_pkl(self.auid_file)
        return self.auids

    @property
    def atoms_paths(self):
        """Returns a list of full paths to the folders that have `atoms.json` objects
        for the latest result set.
        """
        result = []
        for auid in self.auids:
            folder = self.index[str(auid)]
            target = path.join(folder,"atoms.json")
            if path.isfile(target):
                result.append(folder)

        return result 
        
    def rset(self):
        """Returns a :class:`quippy.AtomsList`, one for each config in the
        latest result set.
        """
        from matdb.database.basic import atoms_from_json
        #Return the configurations from this group; it is at the
        #bottom of the stack
        result = quippy.AtomsList()
        for epath in self.atoms_paths:
            result.append(atoms_from_json)
        return result

    def _setup_configs(self, rerun=False):
        """Sets up the database structure for the active set and creates a
        folder for the `calculator` to run in for each config.        
        """
        # We need to make sure that none of the configs in this round
        # of the active learning have already been visited by
        # constructing their auids and then verifying that the auid
        # hasn't been visited before.
        did = len(self.auids)
        for config in self.new_configs:
            auid = hash(tuple([tuple(i) for i in config.cell]),
                        tuple([tuple(i) for i in config.positions]),
                        tuple(config.get_chemical_symbols()))
            if auid in self.auids:
                self.nconfigs -= 1
                continue
            
            dind += 1
            self.create(config,cid=dind)
            self.index[auid] = self.configs[dind]
            self.auids.append(auid)
            
        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.auids,self.auid_file)
        
    def setup(self, rerun=False):
        """Enumerates the desired number of structures and setups up a folder
        for each one.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.

        """
        super(Active, self).setup(self._setup_configs,rerun)           
