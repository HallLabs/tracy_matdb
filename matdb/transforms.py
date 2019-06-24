"""Often it is desirable to transform an atom's object into another atom's
object. This is a one-to-one mapping. If you are looking to apply a
transformation that creates new configurations, you should use a regular
`matdb.database.Group` subclass.

To create your own `transform` function, declare an importable function that accepts a single
positional argument, `at` of type :class:`matdb.atoms.Atoms` and returns an object of
the same type.

Copyright (C) 2019  HALL LABS

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

If you have any questions contact: wmorgan@tracy.com
"""
import numpy as np
from collections import namedtuple
from operator import itemgetter

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

def _get_supers(at, sizes):
    """Returns a list of supercells for the atoms object that are optimal for
    the given target sizes.

    Args:
        at (matdb.atoms.Atoms): atoms object to make supercells for.
        sizes (list): list of `int` target cell *sizes*. These should be multiples of
          :attr:`matdb.atoms.Atoms.n` or they will be rounded to the nearest multiple.
    """
    from supercell import get_supers
    _sizes = [s/at.n if s % at.n == 0 else int(round(s/float(at.n))) for s in sizes]
    maxn = max(_sizes)
    _super = get_supers(at, maxn)
    
    #Create a data container that makes it easier to work with the results.
    HNF = namedtuple("HNF", ['hnf', 'rmin', 'rmax', 'pg', 'det', 'size'])
    supers = {}
    for args in zip(*_super):
        #args[-1] is the determinant of the supercell matrix.
        size = at.n*args[-1]
        if size not in supers:
            supers[size] = []
        supers[size].append(HNF(*(args + (size,))))

    #Next, find the closest matches to the desired sizes for the supercells.
    choices = {}
    available = np.array(sorted(supers.keys()))
    for s in sizes:
        if s in supers:
            choices[s] = supers[s]
        else:
            compromise = np.where(available >= s)
            if len(compromise) > 0 and len(compromise[0]) > 0:
                choice = available[compromise[0][0]]
            else:
                choice = max(available)
            choices[s] = supers[choice]

    #Finally, choose the best supercell for each size. Best in this case means
    #it has the largest rmin and the largest point group size.
    result = {}
    for s, hs in choices.items():
        best = next(iter(sorted(hs, key=itemgetter(1, 3), reverse=True)))
        result[s] = best

    return result

def supercell(at, supercell=None, min_distance=None, max_multiple=None):
    """Transforms an atoms object into a supercell.

    .. note:: If you don't specify an explicit supercell using `supercell`
      keyword, `matdb` will enumerate all unique supercells up to
      `max_multiple`. The selected supercell will be the *smallest* one that has
      at least `min_distance` between every atom and its periodic image.

    Args:
        supercell (list, tuple, numpy.ndarray): a list of `int` multiples that
          forms the supercell matrix. It should have length 3 or 9. If length 3,
          it is assumed to be a diagonal matrix.
        min_distance (float): for enumerated supercell selection, the minimum
          distance (in Angstrom) that should exist between every atom and its
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
      
