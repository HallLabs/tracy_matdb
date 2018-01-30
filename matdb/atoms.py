"""Implementation of Atoms object and AtomsList. Borrows from quippy
code for some of the implementation.
"""

import ase
import numpy as np
from copy import deepcopy

class Atoms(ase.Atoms):
    """An implementation of the :class:`ase.Atoms` object that adds the
    additional attributes of params and properties and reads the atoms
    object in from file.

    """

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None,
                 momenta=None, masses=None, magmoms=None, charges=None,
                 scaled_positions=None, cell=None, pbc=None, constraint=None,
                 calculator=None, info=None, n=None, lattice=None,
                 properties=None, params=None, fixed_size=None, set_species=True,
                 fpointer=None, finalise=True,
                 **readargs):

        if "format" in readargs and readargs["format"] is not None:
            super(Atoms, self).__init__(*list(ase.io.iread(symbols,**readargs)))

        else:
            super(Atoms, self).__init__(symbols, positions, numbers,
                                        tags, momenta, masses, magmoms, charges,
                                        scaled_positions, cell, pbc, constraint,
                                        calculator)

        self.info = {"params":{},"properties":{}}
            
        if properties is not None:
            self.properties.update(properties)
        if params is not None:
            self.params.update(params)

        if info is not None:
            self.params.update(info)

        self._initialised = True

    def __getattribute__(self,name):
        if name=="params":
            return self.info["params"]
        elif name=="properties":
            return self.info["properties"]
        elif name self.info["properties"]:
            return self.info["properties"][name]
        else:
            return ase.Atoms.__getattribute__(self,name)

    def _get_info(self):
        """ASE info dictionary

        Entries are actually stored in the params dictionary.
        """
        info = self.info.copy()
        if "params" in info:
            info.pop("params")
        if "properties" in info:
            infe.pop("properties")
        
        return info

    def _set_info(self, value):
        """Set ASE info dictionary.

        Entries are actually stored in tho params dictionary.  Note
        that clearing Atoms.info doesn't empty params,
        """

        self.params.update(value)
    
    def _set_properties(self,value):
        """Gets the properties from the ASE info dictionary.
        """
        self.info["properties"] = value
    
    def _set_params(self,value):
        """Gets the properties from the ASE info dictionary.
        """
        self.info["params"] = value
        
    def __del__(self):
        attributes = list(vars(self))
        for attr in attributes:
            if isinstance(attr,dict):
                self.attr = {}
            else:
                self.attr = None

    def copy_from(self, other):
        """Replace contents of this Atoms object with data from `other`."""

        from ase.spacegroup import Spacegroup
        self.__class__.__del__(self)
        if isinstance(other, Atoms):
            Atoms.__init__(self, n=other.n, lattice=other.lattice,
                                  properties=other.properties, params=other.params)

            self.cutoff = other.cutoff
            self.cutoff_skin = other.cutoff_skin
            self.nneightol = other.nneightol

        elif isinstance(other, ase.Atoms):
            Atoms.__init__(self, n=0, lattice=np.eye(3))

            # copy params/info dicts
            if hasattr(other, 'params'):
                self.params.update(other.params)
            if hasattr(other, 'info'):
                self.params.update(other.info)
                if 'nneightol' in other.info:
                    self.nneightol = other.info['nneightol']
                if 'cutoff' in other.info:
                    self.set_cutoff(other.info['cutoff'],
                                    other.info.get('cutoff_break'))
                if isinstance(other.info.get('spacegroup', None), Spacegroup):
                    self.params['spacegroup'] = other.info['spacegroup'].symbol

            # create extra properties for any non-standard arrays
            standard_ase_arrays = ['positions', 'numbers', 'masses', 'initial_charges',
                                   'momenta', 'tags', 'initial_magmoms' ]

            self.constraints = deepcopy(other.constraints)

        else:
            raise TypeError('can only copy from instances of matdb.Atoms or ase.Atoms')

        # copy any normal (not Fortran) attributes
        for k, v in other.__dict__.iteritems():
            if not k.startswith('_') and k not in self.__dict__:
                self.__dict__[k] = v
        
class AtomsList(list):
    """An AtomsList like object for storing lists of Atoms objects read in
    from file.
    """

    def __init__(self, source=[], format=None, start=None, stop=None, step=None,
                 rename=None, **kwargs):

        self.source = source
        self.format = format
        self._start  = start
        self._stop   = stop
        self._step   = step
        # if the source has a wildcard or would somehow be a list we
        # need to iterate over it here.
        tmp_ar = [ase.io.read(source_file, format, start, stop, step,
                             rename=rename, **kwargs) for source_file in source]
        list.__init__(self, list(iter(tmp_ar)))

    def __getattr__(self, name):
        if name.startswith('__'):
            # don't override any special attributes
            raise AttributeError

        try:
            return self.source.__getattr__(name)
        except AttributeError:
            try:
                seq = [getattr(at, name) for at in iter(self)]
            except AttributeError:
                raise
            if seq == []:
                return None
            elif type(seq[0]) in (FortranArray, np.ndarray):
                return mockNDarray(*seq)
            else:
                return seq

    def __getslice__(self, first, last):
        return self.__getitem__(slice(first,last,None))

    def __getitem__(self, idx):
        if isinstance(idx, list) or isinstance(idx, np.ndarray):
            idx = np.array(idx)
            if idx.dtype.kind not in ('b', 'i'):
                raise IndexError("Array used for fancy indexing must be of type integer or bool")
            if idx.dtype.kind == 'b':
                idx = idx.nonzero()[0]
            res = []
            for i in idx:
                at = list.__getitem__(self,i)
                res.append(at)
        else:
            res = list.__getitem__(self, idx)
        if isinstance(res, list):
            res = AtomsList(res)
        return res

    def iterframes(self, reverse=False):
        if reverse:
            return reversed(self)
        else:
            return iter(self)

    @property
    def random_access(self):
        return True

    def sort(self, cmp=None, key=None, reverse=False, attr=None):
        """
        Sort the AtomsList in place. This is the same as the standard
        :meth:`list.sort` method, except for the additional `attr`
        argument. If this is present then the sorted list will be
        ordered by the :class:`Atoms` attribute `attr`, e.g.::

           al.sort(attr='energy')

        will order the configurations by their `energy` (assuming that
        :attr:`Atoms.params` contains an entry named `energy` for each
        configuration; otherwise an :exc:`AttributError` will be raised).
        """
        if attr is None:
            list.sort(self, cmp, key, reverse)
        else:
            if cmp is not None or key is not None:
                raise ValueError('If attr is present, cmp and key must not be present')
            list.sort(self, key=operator.attrgetter(attr), reverse=reverse)


    def apply(self, func):
        return np.array([func(at) for at in self])
        
