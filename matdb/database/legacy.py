"""Although `matdb` provides plenty of functionality for generating databases,
there are also many existing databases that have interesting data. Even though
we don't necessarily know exactly how they were made, the data can still be
useful. This module provides a simple class that adapts legacy databases to a
format that can be used by the `matdb` fitting machinery.
"""
from glob import glob
from os import path, mkdir

import numpy as np
from tqdm import tqdm

from matdb import msg
from matdb.atoms import AtomsList, Atoms
from matdb.database.utility import dbconfig, split
from matdb.utility import chdir, dbcat, symlink

def _atoms_conform(dbfile, energy, force, virial):
    """Determines whether the specified database conforms to the constraints for
    the database.

    Args:
        dbfile (str): name of the database file to check.
        energy (str): name of the parameter that describes reference energy.
        force (str): name of the parameter that describes reference forces.
        virial (str): name of the parameter that describes reference virial
          tensor.

    Returns:
        tuple: `(params, force)`, where `params` is a dictionary where keys are
        target parameter names and values are source parameter names. If the
        dictionary has no values, then all the required parameters are already
        in the atoms object. `force` is a boolean that is True if `force` must
        be copied to `ref_force` to conform with `matdb` conventions.
    """
    emsg = "Cannot find {0} under parameter name {1}."
    params = {}
    doforce = False
    a = Atoms(dbfile)

    if energy not in a.params:
        raise ValueError(emsg.format("energy", energy))
    #These next two have unit tests and pytest says these are raised, but
    #coverage doesn't catch them for some reason...
    if force not in a.properties:# pragma: no cover
        raise ValueError(emsg.format("force", force))
    if virial not in a.params:# pragma: no cover
        raise ValueError(emsg.format("virial", virial))

    if energy != "ref_energy":
        params["ref_energy"] = energy
    if force != "ref_force":
        doforce = True
    if virial != "ref_virial":
        params["ref_virial"] = virial

    if "config_type" not in a.params:
        params["config_type"] = None

    return params, doforce

class LegacyDatabase(object):
    """Reperesents a database read in from one or more files that were not
    created by `matdb`.

    Args:
        name (str): the name to use for the database.
        root (str): root directory in which all other database sequences for
          the configurations in the same specification will be stored.
        controller (matdb.database.controller.Controller): instance controlling
          multiple configurations.
        splits (dict): keys are split names; values are `float` *training*
          percentages to use.
        folder (str): path to the directory here the database files are stored.
        pattern (str): or a list of `str` that provide :func:`~fnmatch.fnmatch`
          patterns for selecting files in `folder`. All the provided files will
          be combined to form the final database.
        config_type (str): the configuration type to use in labeling the instances in
          the database files.
        energy (str): name of the parameter that describes DFT/reference energy.
        force (str): name of the parameter that describes DFT/reference forces.
        virial (str): name of the parameter that describes DFT/reference virial
          tensor.

    Attributes:
        dbfiles (list): of `str` file paths that were sources for this legacy
          database.
    """
    def __init__(self, name=None, root=None, controller=None, splits=None,
                 folder=None, pattern=None, config_type=None, energy="ref_energy",
                 force="ref_force", virial="ref_virial", limit=None):
        self.name = name
        self.root = path.join(root, self.name)
        if not path.isdir(self.root):
            # from os import mkdir
            mkdir(self.root)

        self.controller = controller
        self.splits = {} if splits is None else splits
        self.folder = folder

        if self.controller is None:
            self.ran_seed = 0
        else:
            self.ran_seed = self.controller.ran_seed

        self._dbfile = path.join(self.root, "legacy-{}.h5".format(limit))
        """str: path to the combined legacy database, with limits included.
        """
        self._dbfull = path.join(self.root, "legacy.h5")
        """str: path to the combined legacy database, *without* limits.
        """
        self.dbfiles = []
        self.config_type = config_type

        # from matdb.database.utility import dbconfig
        config = dbconfig(self._dbfull)
        if path.isfile(self._dbfile) and len(config) > 0:
            self.dbfiles = [db[0] for db in config["sources"]]
            self.config_type = config["config_type"]
            self.folder = folder
        else:
            # from matdb.utility import dbcat
            if not path.isfile(self._dbfull):
                self._create_dbfull(folder, pattern, energy, force, virial, config_type)

            if limit is not None:
                msg.std("Slicing limit subset of full {} db.".format(self.name))
                full = AtomsList(self._dbfull)
                N = np.arange(len(full))
                np.random.shuffle(N)
                ids = N[0:limit]
                part = full[ids]
                part.write(self._dbfile)
                dbcat([self._dbfull], self._dbfile, docat=False, limit=limit,
                      ids=ids)
            else:
                # from matdb.utility import symlink
                symlink(self._dbfile, self._dbfull)

        #The rest of matdb expects each database to have an atoms object that is
        #representative. Just take the first config in the combined database.
        self.atoms = Atoms(self._dbfile)

    def _create_dbfull(self, folder, pattern, energy, force, virial, config_type):
        """Creates the full combined database.
        """
        # from matdb.utility import chdir, dbcat
        # from glob import glob
        # from tqdm import tqdm
        # from os import path

        #NB! There is a subtle bug here: if you try and open a matdb.atoms.Atoms
        #within the context manager of `chdir`, something messes up with the
        #memory sharing in fortran and it dies. This has to be separate.
        with chdir(folder):
            self.dbfiles = glob(pattern)
        rewrites = []

        for dbfile in self.dbfiles:
            #Look at the first configuration in the atoms list to
            #determine if it matches the energy, force, virial and
            #config type parameter names.
            dbpath = path.join(folder, dbfile)
            params, doforce = _atoms_conform(dbpath, energy, force, virial)
            if len(params) > 0 or doforce:
                msg.std("Conforming database file {}.".format(dbpath))
                al = AtomsList(dbpath)
                outpath = path.join(self.root, dbfile.replace(".xyz",".h5"))
                for ai in tqdm(al):
                    for target, source in params.items():
                        if (target == "config_type" and
                            config_type is not None):
                            ai.params[target] = config_type
                        else:
                            ai.add_param(target,ai.params[source])
                            del ai.params[source]
                            if source in ai.info: #pragma: no cover
                                                  #(if things were
                                                  #dane correctly by
                                                  #the atoms object
                                                  #this should never
                                                  #be used. It exists
                                                  #mainly as a
                                                  #safegaurd.
                                msg.warn("The atoms object didn't properly "
                                         "update the parameters of the legacy "
                                         "atoms object.")
                                del ai.info[source]

                    if doforce:
                        ai.add_property("ref_force",ai.properties[force])
                        del ai.properties[force]

                al.write(outpath)

                #Mark this db as non-conforming so that we created a new
                #version of it.
                rewrites.append(dbfile)

                dbcat([dbpath], outpath, docat=False, renames=params,
                      doforce=doforce)

        # We want a single file to hold all of the data for all the atoms in the database.
        all_atoms = AtomsList()
        for dbfile in self.dbfiles:
            if dbfile in rewrites:
                infile = dbfile.replace(".xyz",".h5")
                all_atoms.extend(AtomsList(path.join(self.root, infile)))
            else:
                dbpath = path.join(folder, dbfile)
                all_atoms.extend(AtomsList(dbpath))

        all_atoms.write(self._dbfull)

        #Finally, create the config file.
        # from matdb.utility import dbcat
        with chdir(folder):
            dbcat(self.dbfiles, self._dbfull, config_type=self.config_type, docat=False)

    def train_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used for training.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{}-train.h5".format(split))

    def holdout_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used to validate the potential fit.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{}-holdout.h5".format(split))

    def super_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used to *super* validate the potential fit.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{}-super.h5".format(split))

    def split(self, recalc=0):
        """Splits the database multiple times, one for each `split` setting in
        the database specification.
        """
        # from matdb.database.utility import split

        # Get the AtomsList object
        subconfs = AtomsList(self._dbfile)

        file_targets = {"train": self.train_file, "holdout": self.holdout_file,
                        "super": self.super_file}

        split(subconfs, self.splits, file_targets, self.root, self.ran_seed,
              dbfile=self._dbfile, recalc=recalc)
