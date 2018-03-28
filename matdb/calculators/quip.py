"""Implements a generalized synchronous calculator interface for using
:class:`quippy.Potential` objects.
"""
from matdb.calculators.basic import SyncCalculator
from quippy.atoms import Atoms
import quippy
import numpy as np

class SyncQuip(quippy.Potential, SyncCalculator):
    """Implements a synchronous `matdb` calculator for QUIP potentials.

    Args:
        atoms (quippy.Atoms): configuration to calculate using VASP.
        folder (str): path to the directory where the calculation should take
          place.
        contr_dir (str): The absolute path of the controller's root directory.
        ran_seed (int or float): the random seed to be used for this calculator.
    """
    def __init__(self, atoms, folder, contr_dir, ran_seed, *args, **kwargs):
        self.args = args[0]
        self.kwargs = kwargs
        self.ran_seed = ran_seed
        self.contr_dir = contr_dir
        self.version = None
        super(SyncQuip, self).__init__(*self.args, **self.kwargs)
        self.atoms = atoms
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
        if 'force' in props:
            props['force'] = np.transpose(props['force'])
        info = atoms.info.copy()
        del info["params"]
        del info["properties"]
        
        kwargs = {"properties":props, "params":params, "positions":atoms.positions,
                  "numbers":atoms.get_atomic_numbers(),
                  "cell":atoms.get_cell(), "pbc":atoms.get_pbc(),
                  "constraint":atoms.constraints, "info":info,
                  "n":len(atoms.positions)}
        return Atoms(**kwargs)

    def calc(self,atoms,**kwargs):
        """Replaces the calc function with one that returns a matdb atoms object.
        
        Args:
            atoms (matdb.atoms.Atoms): the atoms object to do calculations on.
            kwargs (dict): the key work arguments to the :clas:`quippy.Potential` 
              calc function.
        """
        from matdb.atoms import Atoms as matAtoms
        temp_A = self._convert_atoms(atoms)
        super(SyncQuip,self).calc(temp_A,**kwargs)
        for key, val in temp_A.params.items():
            if isinstance(val,np.ndarray):
                new_val = np.array(val)
            else:
                new_val = val
            atoms.add_param(key,new_val)

        for key, val in temp_A.properties.items():
            if isinstance(val,np.ndarray):
                new_val = np.array(val)
            else:
                new_val = val
            if key=="force":
                new_val = np.transpose(new_val)
            atoms.add_property(key,new_val)
        if not np.allclose(atoms.positions,temp_A.positions): 
            atoms.positions = temp_A.positions

    def can_execute(self):
        """Returns `True` if this calculation can calculate properties for the
        specified atoms object.
        """
        return True

    def can_extract(self):
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

    def to_dict(self):
        """Writes the current version number of the code being run to a
        dictionary along with the parameters of the code.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        quip_dict = {"folder":self.folder, "ran_seed":self.ran_seed,
                     "contr_dir":self.contr_dir, "kwargs": self.kwargs,
                     "args": self.args}

        # Need to determine how to find the quip version number.
        # quip_dict["version"] = data["output"][0].strip().split()[0]

        return quip_dict
