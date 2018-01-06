"""Implements an asynchronous calculator for interacting with the AFLOW database
via :class:`aflow.entries.Entry` objects. This allows for asynchronous download
of AFLOW results and makes the overall abstract framework of
:class:`matdb.database.basic.Group` work correctly (since all configurations are
required to have a calculator attached to perform tasks related to cleanup,
checking for execution readiness, etc.

.. todo: This class should ideally inherit from :class:`ase.Calculator`. That
  requires correlating the AFLOW keywords with the physical quantities that the
  calculator can produce.
"""
from .basic import SyncCalculator
from os import path

class AsyncAflow(SyncCalculator):
    """Represents an asynchronous calculator for constructing
    :class:`quippy.Atoms` objects from :class:`aflow.entries.Entry` objects.

    .. note:: Because asynchronous is only supported properly in python 3 and
      this project has to be compatible with both major versions, we don't
      actually support :mod:`asyncio` yet. This means that :meth:`is_executing`
      always returns false because it runs synchronously. This is also why it
      sub-classes :class:`matdb.calculators.basic.SyncCalculator` while the
      documentation talks about async. When async is available in `aflow`, we
      can just change the sub-class and leave most documentation intact.

    Args:
        entry (aflow.entries.Entry): database entry to download data for.
        folder (str): path to the folder where the result will be stored.

    Attributes:
        atoms (quippy.Atoms): atoms object created from the database entry. Is
          `None` until the download is performed.
    """
    def __init__(self, atoms, folder, entry=None, **kwargs):
        self.entry = entry
        self.folder = folder
        self.atoms = atoms
        
    @property
    def entry_file(self):
        """Returns the full path to the pkl file that represents the cached
        :class:`aflow.entries.Entry` object.
        """
        return path.join(folder, "entry.pkl")
        
    def can_execute(self):
        """Returns `True` if the specified folder is ready to execute using the
        calculator's main executable.

        Args:
            folder (str): path to the folder in which the executable will be
              run.
        """
        return self.entry is not None

    def can_cleanup(self):
        """Returns True if the specified folder has completed executing and the results
        are available for use.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        return self.atoms is not None

    def is_executing(self):
        """Returns True if the specified folder is in process of executing. This
        means that files/output has been produced to indicate that the process
        started and that the folder is not ready to cleanup yet.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        return False

    def create(self, rewrite=False):
        """Creates all necessary input files for the calculator's executable.

        Args:
            folder (str): path to the folder in which to create input files.
        """
        if not path.isfile(self.entry_file):
            with open(self.entry_file, "w+") as f:
                pickle.dump(self.entry, f)

    def cleanup(self):
        """Extracts results from completed calculations and sets them on the
        :class:`ase.Atoms` object. This involves executing the actual call to
        the AFLOW database for downloading config information, etc.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass
