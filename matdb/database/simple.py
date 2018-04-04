"""A database of nothing but seed configurations.
"""
from matdb.database import Group
from matdb import msg
from os import path, getcwd, chdir, remove, listdir, mkdir
import numpy as np
from six import string_types
from matdb.atoms import Atoms, AtomsList

class Simple(Group):
    """<<Your description here>>.
    Args:
        atoms (matdb.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.Database): parent sequence to which this database
          belongs. Could also be another :class:`Hessian`.
	  
    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
    """

    def __init__(self, name="simple", **dbargs):
        self.name = name
        self.seeded = True
        dbargs["prefix"] = "S1"
        dbargs["cls"] = Simple
        if "Simple" not in dbargs['root']:
            from os import mkdir
            new_root =path.join(dbargs['root'],"Simple")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Simple, self).__init__(**dbargs)

        self.nconfigs = len(self._seed)
        
        #Make sure that we override the global calculator default values with
        #those settings that we know are needed for good phonon calculations.
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs["calculator"]:
                calcargs.update(dbargs["calculator"])
                self._set_calc_defaults(calcargs)
                dbargs["calculator"] = calcargs            		

    @property
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group. 
        """
        result = []
        for config in self.configs:
            folder = self.index[euid]
            target = path.join(folder,"atoms.h5")
            if path.isfile(target):
                result.append(folder)

        return configs

    @property
    def rset(self):
        """
        Returns:
            list: of :class:`matdb.atoms.Atoms`
        """
        if len(self.sequence) == 0:
            #We are at the bottom of the stack; 
            result = AtomsList()
            for config in self.fitting_configs:
                result.append(Atoms(path.join(config,"atoms.h5")))
            return result
        else:
            #Check where we are in the stack. If we are just below the database,
            #then we want to return <<your description of the rset here>>
	    #If we are not, then we must a parameter grid of sequences
            #to select from.
	    return [p.rset for p in self.sequence.values()]

    def ready(self):
        """Returns True if all the calculations have been completed.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            return len(self.configs) == 1
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
        """<<explanation of setup function.>>

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """
        super(Simple, self).setup(self._setup_configs, rerun)
            
    def _setup_configs(self, rerun):
        """<<explanation of function.>>

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """        
        #We also don't want to setup again if we have the results already.
        if self.ready():
            return

        if not self.is_setup():
            self.create(self.atoms)

        self.jobfile(rerun)
            
    def ready(self):
        """Returns True if all the phonon calculations have been completed, the
        force sets have been created, and the DOS has been calculated.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            return False
        else:
            ready = False
            for p in self.sequence.values():
                if not p.ready():
                    msg.std("{} is not ready. Exiting.".format(p.root), 2)
                    break
            else:
                ready = True
            return ready

    def sub_dict(self):
        """Writes the attributes of this instance of the class to a dictionary.
        """
        return {}
