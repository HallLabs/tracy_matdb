"""Implements a generalized synchronous calculator interface for using
:class:`quippy.Potential` objects.
"""
from matdb.calculators.basic import SyncCalculator
from quippy.atoms import Atoms
import quippy
import numpy as np
from os import path
from matdb.atoms import Atoms as matdbAtoms

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
        self.init_calc(kwargs)
        self.args = args
        self.kwargs = kwargs
        self.ran_seed = ran_seed
        self.contr_dir = contr_dir
        self.version = None
        super(SyncQuip, self).__init__(*self.args, **self.kwargs)
        self.atoms = atoms
        self.folder = folder
        self.name = "Quip"

        if self.key is None:
            #Set default keys; for a GAP potential, we use the name of the
            #file. For other potential types, we just use the type of the
            #potential.
            if self.args[0].lower() == "ip gap":
                pfile = self.kwargs.get("param_filename")
                self.key = path.splitext(pfile)[0]
            else:
                self.key = self.args[0].split()[1].lower()

    def __repr__(self):
        kwargs = ', '.join(map(lambda i: '='.join(i), self.kwargs.items()))
        relroot = self.folder.replace(self.contr_dir, '.')
        return "Quip({}, root={}, {})".format(self.args[0], relroot, kwargs)
                
    def _convert_atoms(self,atoms):
        """Converts an :class:`matdb.atoms.Atoms` object to a
        :class:`quippy.atoms.Atoms` object.

        Args:
            atoms (matdb.atoms.Atoms): the atoms object to 
              perform calculations on.
        """
        if isinstance(atoms, Atoms):
            #No need to do a conversion if we are already a `quippy.Atoms`.
            return atoms
        
        props = atoms.properties.copy()
        params = atoms.params.copy()
        force = None
        if 'force' in props:
            force = props['force'].transpose()
            props['force'] = force

        #Unfortunately, the momenta gets set to zeros by default anytime it is
        #requested from the ASE atoms object. Because of the quippy transpose
        #problem, we don't always transpose. Also, we can't pass it in as a
        #property in the props dict, that raises weird shaping errors in quippy.
        momenta = None
        if "momenta" in props:
            momenta = props["momenta"]
            if momenta.shape[0] != atoms.n:
                momenta = momenta.transpose()
            del props["momenta"]

        info = atoms.info.copy()
        if isinstance(atoms, matdbAtoms):
            del info["params"]
            del info["properties"]
        
        kwargs = {"properties":props, "params":params, "positions":atoms.positions,
                  "numbers":atoms.get_atomic_numbers(), "momenta": momenta,
                  "cell":atoms.get_cell(), "pbc":atoms.get_pbc(),
                  "constraint":atoms.constraints, "info":info,
                  "n":len(atoms.positions)}
        result = Atoms(**kwargs)

        if atoms._calc is not None:
            result.set_calculator(atoms._calc)
        return result

    def set_atoms(self, atoms):
        """Sets the live atoms object for the calculator.
        """
        self.atoms = atoms

    def get_property(self, name, atoms=None, allow_calculation=True, rename=False):
        """Overrides the get_property on ASE atoms object to account for atoms
        conversion between `matdb` and back.
        """
        if atoms is not None:
            atoms = self.atoms

        if not isinstance(atoms, Atoms):
            calcatoms = self._convert_atoms(atoms)
            calcatoms.set_calculator(self)
            calcatoms.set_cutoff(self.cutoff())
            calcatoms.calc_connect()
        else:
            calcatoms = atoms

        #If we already have a quippy.Atoms object, we can call the super method directly
        #and sidestep all the conversions we have going on.
        r = super(SyncQuip, self).get_property(name, calcatoms, allow_calculation)
        if not isinstance(atoms, Atoms):
            self._update_results(atoms, calcatoms, rename)
            
        return r

    def _update_results(self, atoms, calcatoms, rename=False):
        """Updates the parameters and properties on the `atoms` object using
        another one for which calculations were performed.

        Args:
            atoms: passed in atoms object, probably a :class:`ase.Atoms` or a
              :class:`matdb.Atoms`.
            calcatoms (quippy.Atoms): object that the QUIP potential was run
              for.
            rename (bool): when True, include the calculator key as part of the
              quantity names for results.
        """
        #Unfortunately, quippy uses fortran arrays that are
        #transposes. Depending on who calls this function, the atoms object can
        #be ASE, quippy, or `matdb`. We have to perform checks for all relevant
        #arrays that can be transposed.
        for key, val in calcatoms.params.items():
            if isinstance(val,np.ndarray):
                new_val = np.array(val)
            else:
                new_val = val

            if isinstance(atoms, matdbAtoms) and rename:
                atoms.add_param("{}_{}".format(self.key, key), new_val)
            else:
                atoms.params[key] = new_val

        for key, val in calcatoms.properties.items():
            if isinstance(val, np.ndarray):
                new_val = np.array(val)
            else:
                new_val = val
            if key in ["force", "momenta"] and isinstance(atoms, matdbAtoms):
                new_val = np.transpose(new_val)
            if key=="species" and isinstance(atoms, Atoms):
                new_val = list(map(lambda r: ''.join(r), new_val.T))

            if isinstance(atoms, matdbAtoms) and rename:
                atoms.add_property("{}_{}".format(self.key, key), new_val)
            else:
                atoms.add_property(key, new_val)
                
        if not np.allclose(atoms.positions, calcatoms.positions): 
            atoms.positions = calcatoms.positions        
        
    def calc(self, atoms, rename=False, **kwargs):
        """Replaces the calc function with one that returns a matdb atoms object.
        
        Args:
            atoms (matdb.atoms.Atoms): the atoms object to do calculations on.
            kwargs (dict): the key work arguments to the :clas:`quippy.Potential` 
              calc function.
        """
        if not isinstance(atoms, Atoms):
            temp_A = self._convert_atoms(atoms)
            temp_A.set_calculator(self)
            temp_A.set_cutoff(self.cutoff())
            temp_A.calc_connect()
        else:
            temp_A = atoms
            
        super(SyncQuip, self).calc(temp_A, **kwargs)
        if not isinstance(atoms, Atoms):
            self._update_results(atoms, temp_A, rename)

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
