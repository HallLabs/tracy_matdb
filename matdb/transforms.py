"""Often it is desirable to transform an atoms object into another atoms
object. This is a one-to-one mapping. If you are looking to apply a
transformation that creates new configs, you should use a regular
`matdb.database.Group` subclass.

To create your own transform, make an importable function that accepts a single
positional argument, `at` of type :class:`matdb.Atoms` and returns an object of
the same type.
"""
import numpy as np
from ase.build import make_supercell

def conform_supercell(supercell):
    """Conforms the specified supercell to be a 3x3 matrix.

    Args:
        supercell: can be a list of length 3 or 9, or a :class:`numpy.ndarray`
          that has shape (3, 3), (3,), (9,), (3, 1), (1, 3), (9, 1) or (1, 9).
    """
    if supercell is None:
        return
    if isinstance(supercell, (list, tuple)):
        assert len(supercell) == 3 or len(supercell) == 9
        if len(supercell) == 3:
            scell = np.diag(supercell)
        else:
            scell = np.array(supercell).reshape(3, 3)
    elif isinstance(supercell, np.ndarray):
        if supercell.shape == (3, 3):
            scell = supercell
        else:
            scell = conform_supercell(supercell.flatten().tolist())
    return scell

def supercell(at, supercell=None, min_distance=None, max_multiple=None):
    """Transforms the given atoms object into a supercell.

    .. note:: If you don't specify an explicit supercell using `supercell`
      keyword, `matdb` will enumerate all unique supercells up to
      `max_multiple`. The selected supercell will be the *smallest* one that has
      at least `min_distance` between every atom and its periodic image.

    Args:
        supercell (list, tuple, numpy.ndarray): a list of `int` multiples that
          forms the supercell matrix. It should have length 3 or 9. If length 3,
          it is assumed to be a diagonal matrix.
        min_distance (float): for enumerated supercell selection, the minimum
          distance (in Ang.) that should exist between every atom and its
          periodic image in the supercell.
        max_multiple (int): the number of atoms in the supercell is the number
          of atoms in the primitive, times the determinant of the supercell
          matrix. This value specifies the maximum determinant of enumerated
          supercell matrices. `matdb` will keep enumerating until all matrices
          of this determinant have been investigated.
    """
    if supercell is not None:
        if len(supercell) == 3:
            scell = np.identity(supercell)
        elif len(supercell) == 9:
            scell = np.array(supercell).reshape(3, 3)
        else:
            scell = np.array(supercell)
        assert scell.shape == (3, 3)

        return make_supercell(at, scell)
    else:
        raise NotImplementedError("Supercell enumeration still needs a wheel "
                                  "implemented before it can be used here...")        
      
