"""A database of nothing but seed configurations.
"""
from os import path, getcwd, chdir, remove, listdir, mkdir

import numpy as np
from six import string_types

from matdb import msg
from matdb.database import Group
from matdb.atoms import Atoms, AtomsList

class Manual(Group):
    """A basic group that just sets up a calculator for each atoms object
    specified by the user.

    Args:
        name (str): the name of the database group.
        extractable (bool): True if a calculation is to be performed.
        dbargs (dict): a dictionary of arguments for the Group class.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
    """

    def __init__(self, name="manual", extractable=True,
                 **dbargs):
        self.name = name
        self.extractable = extractable
        self.seeded = True
        dbargs["prefix"] = "S1"
        dbargs["cls"] = Manual
        if "Manual" not in dbargs['root']:
            from os import mkdir
            new_root =path.join(dbargs['root'],"Manual")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Manual, self).__init__(**dbargs)

        if not self.extractable:
            self._trainable = False

        self.nconfigs = 1

        #Make sure that we override the global calculator default values with
        #those settings that we know are needed for good phonon calculations.
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs["calculator"]:
                calcargs.update(dbargs["calculator"])
                self._set_calc_defaults(calcargs)
                dbargs["calculator"] = calcargs

    def _set_calc_defaults(self, calcargs):
        """ No implementation for this method.
        """

    @property
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group.
        """
        if len(self.sequence) == 0:
            result = []
            for folder in self.configs.values():
                target = path.join(folder,"atoms.h5")
                if path.isfile(target):
                    result.append(folder)
            return result
        else:
            result = []
            for g in self.sequence.values():
                result.extend(g.fitting_configs)
            return result

    @property
    def rset(self):
        """Returns the reusable set to the next database group.

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
            #then we want to return the atoms objects for all database entries.
	    #If we are not, then we must a parameter grid of sequences
            #to select from.
            result = []
            for g in self.sequence.values():
                result.extend(g.rset)
            return AtomsList(result)

    def setup(self, rerun=False):
        """Creates a folder for each seed configuration.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.
        """
        super(Manual, self).setup(self._setup_configs, rerun)

    def _setup_configs(self, rerun):
        """
        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.
        """
        #We also don't want to setup again if we have the results already.
        if self.ready():
            return

        if not self.is_setup():
            self.create(self.atoms, extractable=self.extractable)

        if self.extractable:
            self.jobfile(rerun)

    def ready(self):
        """Returns True if all the calculations have been completed.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            if not self.extractable and self.is_setup():
                return True
            else:
                #A zero-length sequence can mean we have a set of seeds that
                #were specified, *or* that we have a single seed that is itself
                #an atoms object (instead of a list of atoms objects).
                if (len(self.fitting_configs) == len(self._seed) or
                    (len(self.fitting_configs) == 1 and isinstance(self._seed, Atoms))):
                    return True
                else:
                    return False
        else:
            ready = True
            for p in self.sequence.values():
                if not p.ready():
                    msg.std("{} is not ready. Exiting.".format(p.root), 2)
                    ready = False
                    break
            return ready

    def sub_dict(self):
        """Writes the attributes of this instance of the class to a dictionary.
        """
        return {"name": self.name, "extractable": self.extractable}

    def can_extract(self):
        """Runs post-execution routines to clean-up the calculations. This super class
        implementation only checks that each of the sub-config directories has
        the necessary files needed for extraction of output.
        """
        if not self.extractable:
            return self.is_setup()
        else:
            return super(Manual, self).can_extract()
