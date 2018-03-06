"""Implements a :class:`matdb.databas.basic.Group` for querying the
`aflowlib.org` database. This group is an abstraction of the functionality in
the `aflow` python package. Since that package allows keywords to be defined and
manipulated using overloaded versions of built-in operators, we abstract that
here into a domain-specific language.

Domain-Specific Language
------------------------

.. note:: Keyword names *must* exactly match those in `aflowlib`.

- To describe a keyword, just use its name `enthalpy_formation_atom`.
- To describe a *filter*, enter it as you would in python code, but with each
  element as a list.

.. note:: Operators *must* come from the list described below.

Some example filters:

- Select entries with space-group 216: `["spacegroup_relax", "==", 216]`.
- Select entries for Silicon: `["species", "==", "Si"]`.
- Select band gap greater than 6 eV: `["Egap", ">", 6]`.

You can nest these lists of conditions with logical operators for complex
queries. For example, python code `((Egap > 0) & (Egap < 2)) | ((Egap > 5) &
(Egap < 7))` selects a band gap between 0 and 2 or 5 and 7:

.. code-block:: yaml

  filters:
    - - - ["Egap", ">", 0]
        - "&"
        - ["Egap", "<", 2]
      - "|"
      - - ["Egap", ">", 5]
        - "&"
        - ["Egap", "<", 7]
    - ["species", "==" ,"Si"]
  select:
    - "agl_thermal_conductivity_300K"
  orderby:
    keyword: "agl_thermal_conductivity_300K"
    desc: True

Operators
---------

1. `>` and `<` behave as expected. However, these are overloaded for
   string comparisons in the spirit of the AFLUX endpoint. For example
   `author < "curtarolo"` will match `*curtarolo` and `author >
   "curtarolo" will match `curtarolo*`.
2. `==` behaves as expected for all keywords.
3. `%` allows for string searches. `author % "curtarolo"` matches
   `*curtarolo*`.
4. `~` inverts the filter (equivalent to a boolean `not`).
5. `&` is the logical `and` between two conditions.
6. `|` is the logical `or` between two conditions.

"""
from .basic import Group
import aflow
import operator
from os import path
import pickle

operators = {
    '<': operator.lt,
    '>': operator.gt,
    '==': operator.eq,
    '%': operator.mod,
    '~': operator.inv,
    '&': operator.and_,
    '|': operator.or_
}
"""dict: keys are `str` operators that are possible in AFLOW filters; values are
the corresponding functions in :mod:`operator`.
"""

def kfilter(dsl):
    """Constructs a :class:`aflow.keywords.Keyword` for the given DSL entries.
    
    Args:
        dsl (list): of `str` or other primitive types that defines a filter.
    """
    if any(isinstance(c, list) for c in dsl):
        #This is a nested (aka complex) filter with multiple filters that need
        #to be constructed recursively.
        result = []
        for c in dsl:
            if isinstance(c, list):
                result.append(kfilter(c))
            else:
                result.append(c)
        return kfilter(c)
    else:
        #The list must have at least a keyword name and an operator. Unary
        #operators (such as not) don't need any other arguments.
        assert len(dsl) in [2, 3]
        if len(dsl) == 2:
            #Unary operator case; only two entries.
            op, kw = dsl
            assert op in operators
            kword = get_kw(kw)
            return operators[op](kword)
        
        elif len(dsl) == 3:
            #Has an operator and an argument.
            l, op, right = dsl
            assert op in operators
            left = get_kw(l)
            return operators[op](left, right)

def get_kw(kwstr):
    """Returns the :class:`aflow.keywords.Keyword` instance with the specified
    name.
    """
    return getattr(aflow.K, kwstr)
        
class Aflow(Group):
    """Represents a group of configurations that are downloaded from AFLOW.

    Args:
        catalog (str): one of the catalogs supported on AFLOW: ['icsd', 'lib1',
          'lib2', 'lib3']. Also supports a `list` of catalog names.
        batch_size (int): number of data entries to return per HTTP request.
        filters (list): of DSL-compatible filters; see documentation above.
        select (list): of `str` keyword names to include in the result.
        orderby (dict): with keys `keyword` and `reverse` specifying an optional
          keyword to order the result by.
        exclude (list): of `str` keywords to exclude from the result.
        keywords (dict): keys are keyword obects accessible from `aflow.K`;
          values are desired `str` names in the parameters dictionary of the
          atoms object. See :meth:`aflow.entries.Entry.atoms`.

    Attributes:
        catalog (str): one of the catalogs supported on AFLOW: ['icsd', 'lib1',
          'lib2', 'lib3']. Also supports a `list` of catalog names.
        batch_size (int): number of data entries to return per HTTP request.
        filters (list): of :class:`aflow.keywords.Keyword` for filtering the
          results.
        select (list): of :class:`aflow.keywords.Keyword` to include in the
          result.
        orderby (dict): :class:`aflow.keywords.Keyword` to order the result by.
        reverse (bool): when True, reverse the ordering of the results.
        exclude (list): of :class:`aflow.keywords.Keyword` to exclude from the
          result.
        auids (list): of `str` ids from AFLOW database that are the latest
          result from executing the query.
    """
    def __init__(self, catalog=None, batch_size=100, filters=None, select=None,
                 orderby=None, exclude=None, name="aflow", keywords=None,
                 **dbargs):
        self.name = name
        dbargs["prefix"] = 'A'
        dbargs["calculator"] = {"name": "Aflow"}
        super(Aflow, self).__init__(**dbargs)
        
        self.catalog = catalog
        self.batch_size = batch_size
        
        self.filters = []
        if filters is not None:
            self.filters = list(map(kfilter, filters))
            
        self.select = []
        if select is not None:
            self.select = list(map(get_kw, select))
            
        self.reverse = False
        self.orderby = None
        if orderby is not None:
            self.orderby = orderby["keyword"]
            self.reverse = orderby.get("reverse", False)

        self.exclude = []
        if exclude is not None:
            self.exclude = list(map(get_kw, select))

        self.keywords = {}
        if self.keywords is None:
            self.keywords = {
                aflow.K.energy_cell: "ref_energy",
                aflow.K.forces: "ref_force"
            }

        self.auids = None
        self._load_auids()

    @property
    def auid_file(self):
        """Returns the full path to the auid file for this group.
        """
        return path.join(self.root, "auids.pkl")

    def _load_auids(self):
        """Loads the list of `auid` from the `rset.pkl` file for this database
        group. 
        """
        if self.auids is None:
            self.auids = self.load_pkl(self.auid_file)
        return self.auids

    @property
    def atoms_paths(self):
        """Returns a list of full paths to the folders that have `atoms.json` objects
        for the latest result set.
        """
        result = []
        for auid in self.auids:
            folder = self.index[auid]
            target = path.join(folder, "atoms.h5")
            if path.isfile(target):
                result.append(folder)
    
        return result
    
    @property
    def rset(self):
        """Returns a :class:`matdb.atoms.AtomsList`, one for each config in the
        latest result set.
        """

        from matdb.atoms import Atoms, AtomsList
        result = AtomsList()
        for apath in self.atoms_paths:
            result.append(Atoms(apath))
        return result

    def ready(self):
        """Determines if all the AFLOW configurations specified in the query
        have been downloaded and saved yet.
        """
        #We need to count the number of `atoms.json` files that we have
        #corresponding to the auids in the list.
        return len(self.atoms_paths) == self.nconfigs

    def setup(self, rerun=False):
        """Executes the query against the AFLOW database and downloads the
        configurations specified by the query. Each is created in its own
        folder. However, the `atoms.json` files will not be created until
        :meth:`cleanup` is called.
        """
        super(Aflow, self).setup(self._setup_configs, rerun)

    def _build_query(self):
        """Constructs the :class:`aflow.control.Query` object for requesting
        data from the server.
        """
        result = aflow.search(self.catalog, self.batch_size)
        for f in self.filters:
            result = result.filter(f)
        result = result.select(*self.select).exclude(*self.exclude)
        if self.orderby:
            result = result.orderby(self.orderby, self.reverse)

        if self.nconfigs:
            return result[0:self.nconfigs]
        else:
            return result
    
    def _setup_configs(self, rerun):
        """Sets up the folders for each of the configs retrieved from the AFLOW
        query.

        Args:
            rerun (bool): when True, re-execute the query, even if we have the
              correct number of configurations already.
        """
        from tqdm import tqdm
        if len(self.configs) == self.nconfigs and not rerun:
            return
            
        #Execute the AFLOW query; look at the results and then see which of them
        #need to be created into new folders.
        query = self._build_query()
        auids = []

        try:
            for entry in tqdm(query):
                #No matter what, we need to make sure that we store the auid for
                #the latest query result set.
                auids.append(entry.auid)
                #If the auid is in the index, it has already been downloaded
                #before and has its own folder.
                if entry.auid in self.index:
                    continue

                #This creates the folder and configures the atoms object in the
                #group. However, it does *not* create `atoms.json`, which
                #happens only when cleanup is called.
                atoms = entry.atoms(quippy=True, keywords=self.keywords)                
                cid = self.create(atoms, calcargs={"entry": entry})
                self.index[entry.auid] = self.configs[cid]
        finally:
            self.save_index()

        self.auids = auids
        self.save_pkl(self.auids, self.auid_file)
