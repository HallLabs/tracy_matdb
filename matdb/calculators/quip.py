"""Implements a generalized synchronous calculator interface for using
:class:`quippy.Potential` objects.
"""
from matdb.calculators.basic import SyncCalculator
import quippy

class SyncQuip(quippy.Potential, SyncCalculator):
    """Implements a synchronous `matdb` calculator for QUIP potentials.
    """
    def __init__(self, atoms, folder, calcargs=None, calckw=None):
        self.calcargs = [] if calcargs is None else calcargs
        self.calckw = {} if calckw is None else calckw
        super(SyncQuip, self).__init__(*self.calcargs, **self.calckw)
        self.atoms = atoms
        self.folder = folder
        self.name = "Quip"

    def todict(self):
        return {"calcargs": self.calcargs, "calckw": self.calckw}

    def can_execute(self):
        """Returns `True` if this calculation can calculate properties for the
        specified atoms object.

        Args:
            atoms (quippy.Atoms): config to test executability for.
        """
        return True

    def can_cleanup(self):
        """Returns True if the specified atoms object has completed executing and the
        results are available for use.

        Args:
            atoms (quippy.Atoms): config to check execution completion for.
        """
        return True

    def is_executing(self):
        """Returns True if the specified config is in process of executing.

        Args:
            atoms (quippy.Atoms): config to check execution for.
        """
        return False

    def create(self):
        """Initializes the calculator for the specified atoms object if
        necessary.

        Args:
            atoms (quippy.Atoms): config to initialize for.
        """
        pass
