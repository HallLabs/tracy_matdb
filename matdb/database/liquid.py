#Copyright (C) 2019  HALL LABS
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#If you have any questions contact: wmorgan@tracy.com
"""Class for generating and interacting with a group of
configurations generated from liquid-temperature molecular dynamics.
"""
from os import path

from matdb.database import Group

class LiquidGroup(Group):
    """Represents a sub-sampled molecular dynamics run created at a
    liquid transition temperature.

    Args:
        atoms (matdb.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.controller.DatabaseCollection): parent collection
          to which this database belongs.
        incar (dict): specify additional settings for the INCAR file (i.e.,
          differing from, or in addition to those in the global set).
        kpoints (dict): specify additional settings for the PRECALC file (i.e.,
          differing from, or in addition to those in the global set).
        execution (dict): specify override parameters for the execution
          templates in this database.

    .. note:: Additional attributes are also exposed by the super class
      :class:`~matdb.database.Group`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
    """
    def __init__(self, atoms=None, root=None, parent=None, incar={},
                 kpoints={}, execution={}, nconfigs=None, mdbase="md",
                 name="liquid"):
        self.name = name
        super(LiquidGroup, self).__init__(atoms, incar, kpoints, execution,
                                             path.join(root, self.name),
                                             parent, "L", nconfigs=nconfigs)
        self.mdbase = self.parent.databases[mdbase]

    def ready(self):
        """Determines if this database is finished calculating by testing the
        existence of the xyz database file in the root folder.
        """
        return (path.isfile(path.join(self.root, "output.xyz")) and
                len(self.configs) == self._nsuccess)

    def extract(self, cleanup="default"):
        """Generates the XYZ database file for all the sub-configs in this
        liquid database.

        Args:
            cleanup (str): the level of cleanup to perform after extraction.

        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.
        """
        #First, we need to check that the MD is done; then we can subsample it
        #and run the individual DFT calculations.
        if not self.mdbase.ready():
            return

        if not super(LiquidGroup, self).extract(cleanup=cleanup):
            return

        return self.xyz(config_type="liq")

    def setup(self, rerun=False):
        """Sets up a DFT folder for each of the subsampled configurations in the
        base MD database.

        Args:
            rerun (bool): when True, recreate the job file. If the folders
              don't exist yet, they will still be created.
        """
        if not self.mdbase.ready():
            return

        folders_ok = super(LiquidGroup, self).setup()
        if folders_ok and not rerun:
            return

        #We also don't want to setup again if we have the results already.
        if self.ready():
            return

        if not folders_ok:
            from matdb.atoms import Atoms
            from tqdm import tqdm
            with open(self.mdbase.subsamples) as f:
                for line in tqdm(f):
                    config = line.strip()
                    datoms = Atoms(config, format="POSCAR")
                    self.create(datoms)

        # Last of all, create the job file to execute the job array.
        self.jobfile(rerun)
