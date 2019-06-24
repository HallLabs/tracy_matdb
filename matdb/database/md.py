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
"""Class for generating MD runs that can be subsampled to generate
databases of configurations.
"""
from os import path

import numpy as np

from matdb import msg
from matdb.database import Group

class DynamicsGroup(Group):
    """Represents a molecular dynamics run created at a specific
    temperature.

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
      :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.

    """
    def __init__(self, atoms=None, root=None, parent=None, incar={},
                 kpoints={}, execution={}, nsteps=None,
                 samplerate=100, strains=None, tstart=None, tend=None,
                 supercell=None, name="md"):
        self.name = name
        msg.warn("The DM group is only configured for VASP at this time.")
        super(DynamicsGroup, self).__init__(atoms, incar, kpoints, execution,
                                               path.join(root, self.name),
                                               parent, "D", nconfigs=None)
        self.samplerate = samplerate
        self.nsteps = nsteps
        self.strains = [0] if strains is None else strains
        self.tstart = tstart
        self.tend = tend
        self.supercell = supercell

        if supercell is None:
            self.seed = self.atoms.copy()
        else:
            msg.warn("Not Implemnted: At this time specifying a supercell is not "
                    "yet implemented in `matdb` but will be available in "
                    "latter versions. Using seed configuration instead.")
            self.seed = self.atoms.copy()
            
            
        self._update_incar()
        self._update_kpoints()

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
            "potim": 1.,
            "nsw": self.nsteps*self.samplerate
        }
        for k, v in usuals.items():
            if k not in self.incar:
                self.incar[k] = v

        #We also need to add the temperature settings to the INCAR for this
        #particular configuration.
        self.incar["TEBEG"] = self.tstart
        self.incar["TEEND"] = self.tend

    def _update_kpoints(self):
        """Adds the custom k-points settings (i.e., gamma point only) needed for
        MD.
        """
        usuals = {
            "matdb": "gamma"
        }
        for k, v in usuals.items():
            if k not in self.kpoints:
                self.kpoints[k] = v

    def _xdatcar_ok(self):
        """Determines if the MD run is finished calculating by testing the
        existence of the XDATCAR files with correct number of lines.
        """
        for i, folder in self.configs.items():
            target = path.join(folder, "XDATCAR")
            if path.isfile(target):
                #Make sure we have all the steps we need. There should
                #be at least as many entries as requested configs. We
                #also have the lattice vectors, label and atom counts,
                #which adds 5 lines.
                from matdb.utility import linecount
                minlines = self.nsteps * (self.seed.n + 1)
                nlines = linecount(target)
                if nlines < minlines + 5:
                    return False
            else:
                return False

        return True

    @property
    def subsamples(self):
        """Returns the path to the file that has a list of paths to sub-sampled
        POSCARs.
        """
        return path.join(self.root, "subsamples.dat")

    def _parse_md(self, target):
        """Extracts every frame from XDATCAR that is a multiple of `N` and
        constructs a POSCAR for it.

        Args:
            target (str): path to the directory storing `XDATCAR`.
        """
        from os import path
        from subprocess import Popen, PIPE

        header = []
        li = 0
        Natoms = 0
        current = None
        writes = []

        from tqdm import tqdm
        pbar = tqdm(total=self.nsteps)

        try:
            with open(path.join(target, "XDATCAR")) as f:
                for line in f:
                    if li <= 6:
                        if li == 0:
                            header.append("{} ! MD: {{}}")
                        else:
                            header.append(line.strip())
                        li += 1

                    if li == 7:
                        header[0] = header[0].format(header[5])
                        header.remove(header[5])
                        Natoms = sum(map(int, header[5].split()))
                        li += 1

                    if "configuration" in line:
                        vals = line.split()
                        if len(vals) == 3:
                            NL = int(vals[2])
                        else:
                            NL = int(vals[1].split('=')[-1])

                        #Write the POSCAR for this frame.
                        if current is not None:
                            NC = eval(current[0].split()[-1])
                            outfile = path.join(target, "POSCAR.{}".format(NC/self.samplerate))
                            with open(outfile, 'w') as f:
                                f.write('\n'.join(current))
                            writes.append(outfile + '\n')
                            pbar.update(1)

                            if len(current) != Natoms + 7:
                                emsg = "ERROR: MD {} has {} lines."
                                msg.err(emsg.format(NC, len(current)))

                        if NL % self.samplerate == 0:
                            current = list(header)
                            current[0] = current[0].format(NL)
                            current.append("Direct")
                        else:
                            current = None

                    elif current is not None:
                        current.append(line.strip())

            #Handle the last frame in the set.
            if current is not None:
                NC = eval(current[0].split()[-1])
                outfile = path.join(target, "POSCAR.{}".format(NC))
                with open(outfile, 'w') as f:
                    f.write('\n'.join(current))
                writes.append(outfile + '\n')
                pbar.update(1)
                if len(current) != Natoms + 7:
                    emsg = "ERROR: MD {} has {} lines."
                    msg.err(emsg.format(NC, len(current)))
        finally:
            pbar.close()

        return writes

    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        from matdb.utility import linecount
        return linecount(self.subsamples) == self.nsteps * (len(self.strains))

    def extract(self, cleanup="default"):
        """Parses the XDATCAR files to create a list of configurations
        that can be run using high-accuracy DFT.

        Args:
            cleanup (str): the level of cleanup to perform after extraction.

        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.
        """
        #First, we need to check that the MD is done; then we can subsample it
        #and run the individual DFT calculations.
        if not self._xdatcar_ok():
            msg.std("XDATCAR incomplete; can't extract the MD.", 2)
            return False

        subsamples = []
        for i, folder in self.configs.items():
            subsamples.extend(self._parse_md(folder))

        #Write the list of sub-sample file paths to disk.
        with open(self.subsamples, 'w') as f:
            f.writelines(subsamples)

        return len(subsamples) > 0

    def setup(self, rerun=False):
        """Sets up the initial MD run.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.
        """
        folders_ok = super(DynamicsGroup, self).setup()
        if folders_ok and not rerun:
            return

        #We also don't want to setup again if we have the results already.
        if self.ready():
            return

        if not folders_ok:
            for strain in self.strains:
                if strain == 0:
                    self.create(self.seed, sort=True)
                    continue

                #Strain the base cell for the MD run. Make sure we copy the
                #atoms first so that we don't mess up pointer references.
                datoms = self.seed.copy()
                smat = (1.+strain/100.)**(1./3.)
                datoms.cell = datoms.cell*smat
                self.create(datoms, sort=True)

        # Last of all, create the job file to execute the job array.
        self.jobfile(rerun)
