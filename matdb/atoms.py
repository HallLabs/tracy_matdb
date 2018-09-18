"""Implementation of Atoms object and AtomsList. Borrows from quippy
code for some of the implementation.
"""

import ase
from ase.calculators.singlepoint import SinglePointCalculator
from ase.build import make_supercell
import numpy as np
from copy import deepcopy
from itertools import product

import h5py
from ase import io
from six import string_types
import lazy_import
from os import path
from uuid import uuid4
from matdb import msg
from matdb.transforms import conform_supercell

calculators = lazy_import.lazy_module("matdb.calculators")

def _recursively_convert_units(in_dict, split = False):
    """Recursively goes through a dictionary and converts it's units to be
    numpy instead of standard arrays.
    Args:
        in_dict (dict): the input dictionary.
    Returns:
        a copy of the dict with the entries converted to numpy ints,
        floats, and arrays.
    """
    dict_copy = in_dict.copy()
    for key, item in dict_copy.items():
        if isinstance(item,int):
            dict_copy[key] = np.int64(item)
        elif isinstance(item,float):
            dict_copy[key] = np.float64(item)
        elif isinstance(item,dict):
            dict_copy[key] = _recursively_convert_units(item)
        elif isinstance(item,list):
            if split and isinstance(item[0], Atoms):
                atoms_dict = {}
                for i,a in enumerate(item):
                    atoms_dict[str(i)] = a.to_dict()
                dict_copy[key] = _recursively_convert_units(atoms_dict)
            else:
                dict_copy[key] = np.array(item)
        elif item is None: #pragma: no cover I'm not sure this ever
                           #will happen but better safe than sorry.
            del dict_copy[key]
    return dict_copy

def _calc_name_converter(name):
    """Converts the name returned and saved by the ase calculator to the
    matdb calculator instance name.
    """
    name_dict = {"vasp":"Vasp"}
    return name_dict[name] if name in name_dict else name

class Atoms(ase.Atoms):
    """An implementation of the :class:`ase.Atoms` object that adds the
    additional attributes of params and properties.
    .. note:: All arguments are optional. The user only needs to
    specify the symbols or the atomic numbers for the system not both.
    Args:
        symbols (str): The chemical symbols for the system, i.e., 'Si8' or 'CoNb'
        positions (list): The (x,y,z) position of each atom in the cell.
        numbers (list): List of the atomic numbers of the atoms.
        momenta (list): The momenta of each atom.
        masses (list): The masses of each atom.
        charges (list): The charges of each atom.
        cell (list): The lattice vectors for the cell.
        pbc (list): list of bools for the periodic boundary conditions in x y 
          and z. 
        calculator (object): a `matdb` calculator object.
        info (dict): a dictionary containing other info (this will get stored in 
          the params dictionary.
        n (int): the number of atoms in the cell.
        properties (dict): a dictionary of properties where the keys are the property
          names and the values are a list containing the property value for each atom.
        params (dict): a dictionary of parameters that apply to the entire system.
        group_uuid (str): the uuid for the group.
        uuid (str): a uuid4 str for unique identification.
    .. note:: Additional attributes are also exposed by the super class
      :class:`ase.Atoms`.
    Attributes:
        properties (dict): a dictionary of properties where the keys are the property
          names and the values are a list containing the property value for each atom.
        params (dict): a dictionary of parameters that apply to the entire system.
        n (int): the number of atoms in the cell
        calc (object): the `matdb` calculator to be used for calculations.
    """

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None,
                 momenta=None, masses=None, magmoms=None, charges=None,
                 scaled_positions=None, cell=None, pbc=None, constraint=None,
                 calculator=None, info=None, n=None, celldisp=None,
                 properties=None, params=None, fixed_size=None, set_species=True,
                 fpointer=None, finalise=True, group_uuid=None, uuid=None,
                 **readargs):

        if (symbols is not None and not isinstance(symbols,string_types)) or (
                symbols is not None and path.exists(symbols)):
            try:
                self.copy_from(symbols)
            except TypeError:
                self.read(symbols,**readargs)

        else:
            #NB: make sure that we match exactly on the keyword arguments. Do
            #*NOT* use positional arguments for this constructor.
            super(Atoms, self).__init__(symbols=symbols,
                 positions=positions, numbers=numbers,
                 tags=tags, momenta=momenta, masses=masses,
                 magmoms=magmoms, charges=charges,
                 scaled_positions=scaled_positions,
                 cell=cell, pbc=pbc, celldisp=celldisp,
                 constraint=constraint,
                 calculator=calculator, info=info)

        self.n = n if n is not None else len(self.positions)
        if self.calc is None:
            self.calc = calculator if calculator is not None else None
        else:
            self.calc = calculator if calculator is not None else self.calc

        if "params" not in self.info:
            self.info["params"]={}
        if "properties" not in self.info:
            self.info["properties"]={}

        if isinstance(symbols,ase.Atoms):
            for k, v in symbols.arrays.items():
                if k not in ['positions','numbers']:
                    self.add_property(k,v)
                if k in self.info["params"]: # pragma: no cover This
                                             # should never happen
                                             # check in place just in
                                             # case.
                    del self.info["params"][k]
                if k in self.info: # pragma: no cover This should
                                   # never happen check in place just
                                   # in case.
                    del self.info[k]

        if hasattr(self,"calc"):
            if hasattr(self.calc,"results"):
                for k, v in self.calc.results.items():
                    if k != 'force':
                        self.add_param(k,v)
                    else:
                        self.add_property(k,v)                    

        if properties is not None:
            for k, v in properties.items():
                self.add_property(k,v)
                
        if params is not None:
            for k, v in params.items():
                self.add_param(k,v)

        if info is not None:
            for k, v in info.items():
                if k not in ["params","properties"]:
                    self.add_param(k,v)
                
        if self.info is not None:
            for k, v in self.info.items():
                if k not in ["params","properties"]:
                    self.add_param(k,v)
                    del self.info[k]

        if not hasattr(self, "group_uuid"):
            self.group_uuid = group_uuid
    
        self.uuid = uuid if uuid is not None else str(uuid4())
                
        self._initialised = True

    def make_supercell(self, supercell):
        """Returns a new :class:`matdb.Atoms` object that is a supercell of the
        current one.
        """
        scell = conform_supercell(supercell)
        result = make_supercell(self, scell)
        return Atoms(result)
        
    def add_property(self,name,value):
        """Adds an attribute to the class instance.
        Args:
            name (str): the name of the attribute.
            value: the value/values that are associated with the attribute.
        """
        name = str(name)
        if hasattr(self,name) or name in self.info["properties"]:
            self.info["properties"][name] = value
        else:
            self.info["properties"][name]=value

    def add_param(self,name,value):
        """Adds an attribute to the class instance.
        Args:
            name (str): the name of the attribute.
            value: the value/values that are associated with the attribute.
        """
        name = str(name)
        if hasattr(self,name) or name in self.info["params"]:
            self.info["params"][name] = value
        else:
            self.info["params"][name]=value
        
    def rm_param(self,name):
        """Removes a parameter as attribute from the class instance and info dictionary.
        Args:
            name (str): the name of the attribute.
        """
        if name in self.info["params"]:
            del self.info["params"][name]

    def rm_property(self, name):
        """Removes a property as attribute from the class instance and info dictionary.
        Args:
            name (str): the name of the property/attribute.
        """
        if name in self.info["properties"]:
            del self.info["properties"][name]
            
    def __del__(self):
        attributes = list(vars(self))
        for attr in attributes:
            if isinstance(getattr(self,attr),dict):
                self.attr = {}
            else:
                self.attr = None

    def __getattr__(self, name):
        if name in ["params", "properties"]:
            return self.info[name]
        else:
            _dict = object.__getattribute__(self, "__dict__")
            if "info" in _dict:
                info = object.__getattribute__(self, "info")
                if "params" in info and name in info["params"]:
                    return info["params"][name]
                elif "properties" in info and name in info["properties"]:
                    return info["properties"][name]
                else:
                    return object.__getattribute__(self, name)
            else:
                return object.__getattribute__(self, name)

    def __setattr__(self, name, value):
        if name in ["params", "properties"]:
            self.info[name] = value
        else:
            if "info" in object.__getattribute__(self, "__dict__"):
                info = object.__getattribute__(self, "info")
                if "params" in info and name in info["params"]:
                    info["params"][name] = value
                elif "properties" in info and name in info["properties"]:
                    info["properties"][name] = value
                else:
                    return super(Atoms, self).__setattr__(name, value)
            else:
                return super(Atoms, self).__setattr__(name, value)
        
    def copy(self):
        """Returns a copy of this atoms object that has different pointers to
        self, values, etc.
        """
        result = Atoms()
        result.copy_from(self)
        return result
                
    def copy_from(self, other):
        """Replace contents of this Atoms object with data from `other`."""

        from ase.spacegroup import Spacegroup
        self.__class__.__del__(self)
        
        if isinstance(other, Atoms):
            # We need to convert the attributes of the other atoms
            # object so that we can initialize this one properly.
            symbols = other.get_chemical_symbols()
            symbols = ''.join([i+str(symbols.count(i)) for i in set(symbols)])

            magmoms = None
            if hasattr(other, "magnetic_moments") and other.magnetic_moments is not None:
                #Call the get in this try block would setup a new calculator to try and
                #calculate the moments. We are interested in a *copy*, meaning that the
                #quantity should already exist.
                try:
                    magmoms = other.get_magnetic_moment()
                except:
                    pass
            try:
                charges = other.get_charges()
            except:
                charges = None
            try:
                constraint = other.constraint
            except:
                constraint = None
                
            masses = other.get_masses()
            momenta = other.get_momenta()
            info = deepcopy(other.info)
            group_uuid = other.group_uuid
            
            self.__init__(symbols=symbols, positions=other.positions, n=other.n,
                          properties=other.properties, magmoms=magmoms,
                          params=other.params, masses=masses, momenta=momenta,
                          charges=charges, cell=other.cell, pbc=other.pbc,
                          constraint=constraint, info=info, calculator=other.calc,
                          group_uuid = group_uuid)

        elif isinstance(other, ase.Atoms):
            super(Atoms, self).__init__(other)
            if "params" not in self.info:
                self.info["params"]={}
            if "properties" not in self.info:
                self.info["properties"]={}

            # copy info dict
            if hasattr(other, 'info'):
                self.params.update(other.info)
                if 'nneightol' in other.info:
                    self.add_param("nneightol",other.info['nneightol'])
                if 'cutoff' in other.info:
                    self.add_param("cutoff",other.info['cutoff'])
                    self.add_param("cutoff_break",other.info.get('cutoff_break'))

            self.constraints = deepcopy(other.constraints)
            self.group_uuid = None
            self.uuid = str(uuid4())
            self.n = len(other)

        else:
            raise TypeError('can only copy from instances of matdb.Atoms or ase.Atoms')
        
        # copy any normal attributes we've missed
        for k, v in other.__dict__.iteritems(): #pragma: no cover
            if not k.startswith('_') and k not in self.__dict__:
                self.__dict__[k] = v

    def read(self,target="atoms.h5",**kwargs):
        """Reads an atoms object in from file.
        Args:
            target (str): The path to the target file. Default "atoms.h5".
        """

        frmt = target.split('.')[-1]
        if frmt == "h5" or frmt == "hdf5":
            from matdb.io import load_dict_from_h5
            with h5py.File(target,"r") as hf:
                data = load_dict_from_h5(hf)
            if "atom" in data.keys()[0]:
                data = data[data.keys()[0]]
            self.__init__(**data)
            if "calc" in data:
                calc = getattr(calculators, _calc_name_converter(data["calc"]))
                args = data["calc_args"] if "calc_args" in data else None
                kwargs = data["calc_kwargs"] if 'calc_kwargs' in data else None
                if args is not None:
                    if kwargs is not None:
                        calc = calc(self, data["folder"], data["calc_contr_dir"],
                                    data["calc_ran_seed"], *args, **kwargs)
                    else:
                        calc = calc(self, data["folder"], data["calc_contr_dir"],
                                    data["calc_ran_seed"], *args)
                else: #pragma: no cover This case has never come up in
                      #testing, however we wil keep it here to be
                      #verbose.
                    if kwargs is not None:
                        calc = calc(self, data["folder"], data["calc_contr_dir"],
                                    data["calc_ran_seed"], **kwargs)
                    else: 
                        calc = calc(self, data["folder"], data["calc_contr_dir"],
                                    data["calc_ran_seed"])
                self.set_calculator(calc)

        else:
            self.__init__(io.read(target,**kwargs))

    def to_quippy(self):
        """Converts this :class:`matdb.atoms.Atoms` object to a
        :class:`quippy.atoms.Atoms` object.

        Args:
            atoms (matdb.atoms.Atoms): the atoms object to perform calculations
            on.
        """
        import quippy
        props = self.properties.copy()
        params = self.params.copy()
        for k, v in props.items():
            if k == "momenta":
                continue
            if k == "species" and len(v.shape) == 2:
                #Weird string handling inside of quippy means we have to treat
                #this specially...
                props[k] = list(map(lambda r: ''.join(r), v.T))
            elif len(v.shape) == 2 and v.shape[0] != 3:
                props[k] = v.transpose()
                
        #Unfortunately, the momenta gets set to zeros by default anytime it is
        #requested from the ASE atoms object. Because of the quippy transpose
        #problem, we don't always transpose. Also, we can't pass it in as a
        #property in the props dict, that raises weird shaping errors in quippy.
        momenta = None
        if "momenta" in props:
            momenta = props["momenta"]
            if momenta.shape[0] != self.n:
                momenta = momenta.transpose()
            del props["momenta"]

        info = self.info.copy()
        if isinstance(self, Atoms):
            del info["params"]
            del info["properties"]

        kwargs = {"properties":props, "params":params, "positions": self.positions,
                  "numbers": self.get_atomic_numbers(), "momenta": momenta,
                  "cell": self.get_cell(), "pbc": self.get_pbc(),
                  "constraint": self.constraints, "info":info,
                  "n": len(self.positions)}

        return quippy.Atoms(**kwargs)
            
    def S(self, cutoff=5., nmax=8, lmax=8, sigma=0.5, trans_width=0.5,
          average=True, normalize=True):
        """Returns the SOAP vectors for each environment in the specified atoms object.
        """
        from quippy.descriptors import Descriptor as D
        #Convert to quippy atoms so that we can get the SOAP vectors.
        atoms = self.to_quippy()
        descstr = ("soap cutoff={0:.1f} n_max={1:d} l_max={2:d} "
                   "atom_sigma={3:.2f} n_species={6} species_Z={{{4}}} n_Z={6} Z={{{9}}} "
                   "trans_width={5:.2f} normalise={7} average={8}")

        
        Z = np.unique(atoms.get_atomic_numbers())
        savg = 'T' if average else 'F'
        soaps = []

        for Z1, Z2 in product(Z, repeat=2):
            soapstr = descstr.format(cutoff, nmax, lmax, sigma, str(Z1),
                                     trans_width, 1, 'F', savg, str(Z2))
            soapy = D(soapstr)
            msg.info(soapstr, 2)
            atoms.set_cutoff(soapy.cutoff())
            atoms.calc_connect()
            PZ = soapy.calc(atoms)
            if average:
                soaps.append(PZ["descriptor"].flatten())
            else:
                soaps.append(PZ["descriptor"])

        P = np.hstack(soaps)
        if normalize:
            return P/np.linalg.norm(P)
        else:
            return P
            
    def to_dict(self):
        """Converts the contents of a :class:`matdb.atoms.Atoms` object to a
        dictionary so it can be saved to file.
        Args:
            atoms (matdb.atams.Atoms): the atoms object to be converted to 
              a dictionary
        Returns:
            A dictionary containing the relavent parts of an atoms object to 
            be saved.
        """
        import sys
        from matdb import __version__
        
        data = {}
        data["n"] = np.int64(len(self.positions))
        data["pbc"] = np.array(self.pbc)
        data["params"] = _recursively_convert_units(self.params.copy())
        data["properties"] = {}
        for prop, value in self.properties.items():
            if value is not None:
                if prop not in ["pos", "species", "Z", "n_neighb", "map_shift"]:
                    data["properties"].update(_recursively_convert_units({prop:value}))

        data["positions"] = np.array(self.positions)
        data["cell"] = np.array(self.cell)
        if self.calc is not None and not isinstance(self.calc, SinglePointCalculator):
            calc_dict = self.calc.to_dict()
            data["calc"] = self.calc.name
            data["calc_contr_dir"] = calc_dict["contr_dir"]
            if "version" in calc_dict:
                data["calc_version"] = calc_dict["version"] 
            if hasattr(self.calc,"args"):
                data["calc_args"] = self.calc.args
            if "kwargs" in calc_dict:
                data["calc_kwargs"] = _recursively_convert_units(calc_dict["kwargs"])
            if hasattr(self.calc,"folder"):
                data["folder"] = self.calc.folder
            if hasattr(self.calc,"ran_seed"):
                data["calc_ran_seed"] = np.float64(self.calc.ran_seed)
            if hasattr(self.calc, "kpoints") and self.calc.kpoints is not None:
                data["calc_kwargs"]["kpoints"] = _recursively_convert_units(self.calc.kpoints)
            if hasattr(self.calc, "potcars") and self.calc.kpoints is not None:
                data["calc_kwargs"]["potcars"] = _recursively_convert_units(self.calc.potcars)
            
        symbols = self.get_chemical_symbols()
        data["symbols"] = ''.join([i+str(symbols.count(i)) for i in set(symbols)])
        if self.group_uuid is not None:
            data["group_uuid"] = self.group_uuid
        data["uuid"] = self.uuid
        data["python_version"] = sys.version
        data["version"] = np.array(__version__)
        return data

    def write(self,target="atoms.h5",**kwargs):
        """Writes an atoms object to file.
        Args:
            target (str): The path to the target file. Default is "atoms.h5".
        """

        frmt = target.split('.')[-1]
        if frmt == "h5" or frmt == "hdf5":
            from matdb.io import save_dict_to_h5
            with h5py.File(target,"w") as hf:
                data = self.to_dict()
                save_dict_to_h5(hf,data,'/')
        else:
            io.write(target,self,**kwargs)
            
class AtomsList(list):
    """An AtomsList like object for storing lists of Atoms objects read in
    from file.
    """

    def __init__(self, source=[], frmt=None, start=None, stop=None, step=None,
                 **readargs):

        self.source = source
        self.frmt = frmt
        self._start  = start
        self._stop   = stop
        self._step   = step
        # if the source has a wildcard or would somehow be a list we
        # need to iterate over it here.
        tmp_ar = None
        if isinstance(source,list) and len(source)>0:
            if not isinstance(source[0],Atoms):
                for source_file in source:
                    self.read(source_file,**readargs)
            else:
                tmp_ar = source
        elif isinstance(source,list) and len(source) == 0:
            tmp_ar = []
        else:
            if isinstance(source,Atoms):
                tmp_ar = [source]
            else:
                self.read(source,**readargs)

        if tmp_ar is not None:
            list.__init__(self, list(iter(tmp_ar)))

    def __getattr__(self, name):
        if name.startswith('__'):
            # don't override any special attributes
            raise AttributeError
        if name =="get_positions":
            # In order to write out using ase we need to not support
            # this attribute.
            raise AttributeError
        try:
            return self.source.__getattr__(name)
        except AttributeError:
            seq = []
            for at in iter(self):
                try:
                    seq.append(getattr(at, name))
                except AttributeError:
                    pass
            if seq == []: #pragma: no cover
                return None
            else:
                return seq

    def __getslice__(self, first, last):
        return self.__getitem__(slice(first,last,None))

    def __getitem__(self, idx):
        if isinstance(idx, list) or isinstance(idx, np.ndarray):
            idx = np.array(idx)
            if idx.dtype.kind not in ('b', 'i'):
                raise IndexError("Array used for fancy indexing must be of type integer or bool")
            if idx.dtype.kind == 'b': #pragma: no cover
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
        import operator
        if attr is None:
            list.sort(self, cmp, key, reverse)
        else:
            if cmp is not None or key is not None:
                raise ValueError('If attr is present, cmp and key must not be present')
            list.sort(self, key=operator.attrgetter(attr), reverse=reverse)


    def apply(self, func):
        return np.array([func(at) for at in self])
        
    def read(self,target,**kwargs):
        """Reads an atoms object in from file.
        
        Args:
            target (str): The path to the target file.
            kwargs (dict): A dictionary of arguments to pass to the ase 
              read function.
        """
        frmt = target.split('.')[-1]
        if frmt == "h5" or frmt == "hdf5":
            from matdb.io import load_dict_from_h5
            with h5py.File(target,"r") as hf:
                data = load_dict_from_h5(hf)
            # If the data was read in from an hdf5 file written by the
            # AtomsList object then each atom will have a tag with it. We check
            # for this by looking for the word 'atom' inside the first key, if
            # it's present we assume that all the contents of the file are an
            # atoms list. If it's not then we assume this is a single atoms
            # object.

            #NB! It is possible that the atoms list object could be an *empty*
            #AtomsList that was written to disk. In that case, just use an empty
            #list.
            if len(data) == 0:
                atoms = []
            elif "atom" in data.keys()[0]:
                if isinstance(data.values()[0],dict):
                    atoms = [Atoms(**d) for d in data.values()]
                elif isinstance(data.values()[0],string_types):
                    atoms = [Atoms(d) for d in data.values()]
                else: #pragma: no cover
                    msg.err("The data format {} isn't supported for reading AtomLists "
                            "from hdf5 files.".format(type(data.values()[0])))
            else:
                atoms = [Atoms(target,**kwargs)]
            if len(self) >0:
                self.extend(atoms)
            else:
                self.__init__(atoms)
        else:
            atoms = [Atoms(d) for d in io.read(target,index=':',**kwargs)]
            if len(self) >0:
                self.extend(atoms)
            else:
                self.__init__(atoms)
            
    def write(self,target,**kwargs):
        """Writes an atoms object to file.
        Args:
            target (str): The path to the target file.
            kwargs (dict): A dictionary of key word args to pass to the ase 
              write function.
        """

        frmt = target.split('.')[-1]
        if frmt == "h5" or frmt == "hdf5":
            from matdb.io import save_dict_to_h5
            with h5py.File(target,"w") as hf:
                for atom in self:
                    data = atom.to_dict()
                    hf.create_group("atom_{}".format(data["uuid"]))
                    save_dict_to_h5(hf,data,"/atom_{}/".format(data["uuid"]))
        else:
            io.write(target,self)
