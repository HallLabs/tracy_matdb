"""Class for generating and interacting with a database of
configurations generated from liquid-temperature molecular dynamics.
"""
from .basic import Database
class LiquidDatabase(Database):
    """Represents a sub-sampled molecular dynamics run created at a
    specific temperature.


    Args:
        atoms (quippy.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        incar (dict): specify additional settings for the INCAR file (i.e.,
          differing from, or in addition to those in the global set).
        kpoints (dict): specify additional settings for the PRECALC file (i.e.,
          differing from, or in addition to those in the global set).

    .. note:: Additional attributes are also exposed by the super class
      :class:`Database`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
    """
    name = "liquid"
        
    def __init__(self, atoms=None, root=None, parent=None, incar={},
                 kpoints={}, execution={}, nconfigs=None,
                 samplerate=100, strain=3.):
        super(LiquidDatabase, self).__init__(atoms, incar, kpoints, execution,
                                             path.join(root, self.name),
                                             parent, "L", nconfigs=nconfigs)
        self.samplerate = samplerate
        self.nsteps = nconfigs*samplerate
        self.mdpdir = path.join(self.root, "p")
        self.mdmdir = path.join(self.root, "m")
        
        from os import mkdir
        if not path.isdir(self.mdpdir):
            mkdir(self.mdpdir)
        if not path.isdir(self.mdmdir):
            mkdir(self.mdmdir)

        self._update_incar()

    def _update_incar(self):
        """Adds the usual settings for the INCAR file when performing
        MD calculations. They are only added if they weren't already
        specified in the config file.
        """
        usuals = {
            "maxmix": 40, # reuse mixer from one MD step to next          
            "ncore": 4,   # one orbital on 4 cores                        
            "nelmin": 4,  # minimum 4 steps per time step, avoid breaking after 2 steps
            "ibrion": 0,
            "nwrite": 0,
            "lcharg": False,
            # canonic (Nose) MD with XDATCAR updated every N steps
            "nblock": self.samplerate, 
            "smass": 3,
            "potim": 1.
        }
        for k, v in usuals.items():
            if k not in self.incar:
                self.incar[k] = v
