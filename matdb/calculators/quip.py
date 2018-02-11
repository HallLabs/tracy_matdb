"""Implements a generalized synchronous calculator interface for using
:class:`quippy.Potential` objects.
"""
from matdb.calculators.basic import SyncCalculator
from quippy.atoms import Atoms
import quippy

class SyncQuip(quippy.Potential, SyncCalculator):
    """Implements a synchronous `matdb` calculator for QUIP potentials.
    """
    def __init__(self, atoms, folder, calcargs=None, calckw=None):
        self.calcargs = [] if calcargs is None else calcargs
        self.calckw = {} if calckw is None else calckw
        super(SyncQuip, self).__init__(*self.calcargs, **self.calckw)
        # self.atoms = atoms
        self._convert_atoms(atoms)
        self.folder = folder
        self.name = "Quip"

    def _convert_atoms(self,atoms):
        """Converts an :class:`matdb.atoms.Atoms` object to a
        :class:`quippy.atoms.Atoms` object.

        Args:
            atoms (matdb.atoms.Atoms): the atoms object to 
              perform calculations on.
        """
        props = atoms.properties.copy()
        params = atoms.params.copy()
        del atoms.info['properties']
        del atoms.info['params']
        
        kwargs = {"properties":props, "params":params, "positions":atoms.positions,
                  "numbers":atoms.get_atomic_numbers(),
                  "cell":atoms.get_cell(), "pbc":atoms.get_pbc(),
                  "constraint":atoms.constraints, "info":atoms.info}
        if atoms.calc is not None:
            kwargs["calculator"]=atoms.calc
            kwargs["momenta"]=atoms.get_momenta()
            kwargs["masses"]=atoms.get_masses()
            kwargs["magmons"]=atoms.get_magnetic_moments()
            kwargs["charges"]=get_charges()
        self.atoms = Atoms(**kwargs)
        
    def todict(self):
        return {"calcargs": self.calcargs, "calckw": self.calckw}

    # def calc(self,kwargs):
    #     """Replaces the calc function with one that returns a matdb atoms object.
        
    #     Args:
    #         kwargs (dict): the key work arguments to the :clas:`quippy.Potential` 
    #           calc function.
    #     """
    #     import matdb
    #     temp_A = self._convert_atoms()
    #     super(SyncQuip,self).calc(temp_A,**kwargs)
    #     self.aotms = matdb.atoms.Atoms(temp_A)

    def can_execute(self):
        """Returns `True` if this calculation can calculate properties for the
        specified atoms object.
        """
        return True

    def can_cleanup(self):
        """Returns True if the specified atoms object has completed executing and the
        results are available for use.
        """
        return True

    def is_executing(self):
        """Returns True if the specified config is in process of executing.
        """
        return False

    def create(self):
        """Initializes the calculator for the specified atoms object if
        necessary.
        """
        pass
