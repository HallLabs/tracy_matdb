"""Implements an abstract base class that all `matdb` calculators need to
implement. 
"""
import abc
from matdb.base import abstractstatic

class AsyncCalculator(object):
    """Represents a calculator such as :class:`ase.Calculator` that can be run
    in multiple stages to compute properties (such as in a distributed or HCP
    environment).
    """
    __metaclass__ = abc.ABCMeta
    @abstractstatic
    def from_folder(folder):
        """Reconstructs the calculator class attached to an atoms object from a folder.
        """
        pass
    
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
    def can_cleanup(self, folder):
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
        started and that the folder is not ready to cleanup yet.

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
    def cleanup(self, folder):
        """Extracts results from completed calculations and sets them on the
        :class:`ase.Atoms` object.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass

class SyncCalculator(object):
    """Represents a calculator such as :class:`ase.Calculator` that can be run
    synchronously stages to compute properties (does not require a distributed
    or HCP environment).
    """
    @abc.abstractmethod
    def can_execute(self):
        """Returns `True` if this calculation can calculate properties for the
        specified atoms object.

        Args:
            atoms (quippy.Atoms): config to test executability for.
        """
        pass

    @abc.abstractmethod
    def can_cleanup(self):
        """Returns True if the specified atoms object has completed executing and the
        results are available for use.

        Args:
            atoms (quippy.Atoms): config to check execution completion for.
        """
        pass

    @abc.abstractmethod
    def is_executing(self):
        """Returns True if the specified config is in process of executing.

        Args:
            atoms (quippy.Atoms): config to check execution for.
        """
        pass

    @abc.abstractmethod
    def create(self, rewrite=False):
        """Creates all necessary input files for the calculator's executable.

        Args:
            folder (str): path to the folder in which to create input files.
        """
        pass
