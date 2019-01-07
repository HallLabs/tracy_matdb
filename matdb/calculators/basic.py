"""Implements an abstract base class that all `matdb` calculators need to
implement. 
"""
import abc
import re
from os import path
from hashlib import sha1

from matdb.base import abstractstatic
from matdb.utility import recursive_getattr, recursive_setattr, config_specs
from matdb.calculators.utility import paths

_rxhash = re.compile(r"\b[0-9a-f]{32,40}\b", re.I)
"""Regex object to match SHA1 hashes."""

class AsyncCalculator(object):
    """Represents a calculator such as :class:`ase.Calculator` that can be run
    in multiple stages to compute properties (such as in a distributed or HCP
    environment).

    Attributes:
        key (str): short, lower-case name to identify the calculator type.
    """
    key = None
    pathattrs = []
    __metaclass__ = abc.ABCMeta

    def init_calc(self, kwargs):

        if "key" in kwargs:
            self.key = kwargs.pop("key")

        namehash = str(sha1(config_specs["name"].encode("ASCII")).hexdigest())
        for pathattr in self.pathattrs:
            attrval = recursive_getattr(kwargs, pathattr)
            #Test to see if this is a hash reference to a globally-stored
            #calculator absolute path.
            if attrval is not None and _rxhash.match(attrval):
                abspath = paths[namehash][self.key][attrval]
                #Overwrite the value of the hash with the actual absolute path
                #to the directory.
                recursive_setattr(kwargs, pathattr, abspath)

    @property
    def energy_name(self):
        """Returns the name of the energy property that this trainer writes onto
        an atoms object when using its calculator to calculate energy.
        """
        return "{}_energy".format(self.key)

    @property
    def force_name(self):
        """Returns the name of the force property that this trainer writes onto
        an atoms object when using its calculator to calculate force.
        """
        return "{}_force".format(self.key)

    @property
    def stress_name(self):
        """Returns the name of the force property that this trainer writes onto
        an atoms object when using its calculator to calculate force.
        """
        return "{}_stress".format(self.key)

    @property
    def virial_name(self):
        """Returns the name of the virial property that this trainer writes onto
        an atoms object when using its calculator to calculate virial.
        """
        return "{}_virial".format(self.key)
            
    @abc.abstractmethod
    def can_execute(self, folder):
        """Returns `True` if the specified folder is ready to execute using the
        calculator's main executable.

        Args:
            folder (str): path to the folder in which the executable will be
              run.
        """
        pass

    @abc.abstractmethod
    def can_extract(self, folder):
        """Returns True if the specified folder has completed executing and the results
        are available for use.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass

    @abc.abstractmethod
    def is_executing(self, folder):
        """Returns True if the specified folder is in process of executing. This
        means that files/output has been produced to indicate that the process
        started and that the folder is not ready to extract yet.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass

    @abc.abstractmethod
    def create(self, folder, rewrite=False):
        """Creates all necessary input files for the calculator's executable.

        Args:
            folder (str): path to the folder in which to create input files.
        """
        pass

    @abc.abstractmethod
    def extract(self, folder, cleanup="default", asis=False):
        """Extracts results from completed calculations and sets them on the
        :class:`ase.Atoms` object.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass

    @abc.abstractmethod
    def cleanup(self, folder, clean_level="default"):
        """Extracts results from completed calculations and sets them on the
        :class:`ase.Atoms` object.

        Args:
            folder (str): path to the folder in which the executable was run.
            clean_level (str): the level of cleanup to perform.
        """
        pass

    @abc.abstractmethod
    def to_dict(self):
        """Writes the current version number of the code being run to a
        dictionary along with the parameters of the code.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass

class SyncCalculator(object):
    """Represents a calculator such as :class:`ase.Calculator` that can be run
    synchronously stages to compute properties (does not require a distributed
    or HCP environment).

    Attributes:
        key (str): short, lower-case name to identify the calculator type.
    """
    key = None
    pathattrs = []
    
    def init_calc(self, kwargs):
        if "key" in kwargs:
            self.key = kwargs.pop("key")

        namehash = str(sha1(config_specs["name"].encode("ASCII")).hexdigest())        
        #This duplication with the AsyncCalculator hasn't been sufficient to
        #warrant yet another super class, so we just have it again here.
        for pathattr in self.pathattrs:
            attrval = recursive_getattr(kwargs, pathattr)
            #Test to see if this is a hash reference to a globally-stored
            #calculator absolute path.
            if attrval is not None and _rxhash.match(attrval):
                abspath = paths[namehash][self.key][attrval]
                #Overwrite the value of the hash with the actual absolute path
                #to the directory.
                recursive_setattr(kwargs, pathattr, abspath)
            
    @property
    def energy_name(self):
        """Returns the name of the energy property that this trainer writes onto
        an atoms object when using its calculator to calculate energy.
        """
        return "{}_energy".format(self.key)
    @property
    def force_name(self):
        """Returns the name of the force property that this trainer writes onto
        an atoms object when using its calculator to calculate force.
        """
        return "{}_force".format(self.key)
    @property
    def virial_name(self):
        """Returns the name of the virial property that this trainer writes onto
        an atoms object when using its calculator to calculate virial.
        """
        return "{}_virial".format(self.key)
            
    @abc.abstractmethod
    def can_execute(self):
        """Returns `True` if this calculation can calculate properties for the
        specified atoms object.

        Args:
            atoms (matdb.Atoms): config to test executability for.
        """
        pass

    @abc.abstractmethod
    def can_extract(self):
        """Returns True if the specified atoms object has completed executing and the
        results are available for use.

        Args:
            atoms (matdb.Atoms): config to check execution completion for.
        """
        pass

    @abc.abstractmethod
    def is_executing(self):
        """Returns True if the specified config is in process of executing.

        Args:
            atoms (matdb.Atoms): config to check execution for.
        """
        pass

    @abc.abstractmethod
    def create(self, rewrite=False):
        """Creates all necessary input files for the calculator's executable.

        Args:
            folder (str): path to the folder in which to create input files.
        """
        pass

    @abc.abstractmethod
    def to_dict(self):
        """Writes the current version number of the code being run to a
        dictionary along with the parameters of the code.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass
