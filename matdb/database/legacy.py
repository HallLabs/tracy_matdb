"""Although `matdb` provides plenty of functionality for generating databases,
there are also many existing databases that have interesting data. Even though
we don't necessarily know exactly how they were made, the data can still be
useful. This module provides a simple class that adapts legacy databases to a
format that can be used by the `matdb` fitting machinery.
"""
from os import path
import quippy
import numpy as np
from matdb import msg

def _quick_write(atlist, outpath):
    """Writes the atoms list to file using the
    :class:`~quippy.cinoutput.CInOutputWriter` for optimization.

    Args:
        atlist (quippy.AtomsList): atoms to write to XYZ file.
        outpath (str): full path to the location of the output file to write.
    """
    import quippy.cinoutput as qcio
    from tqdm import tqdm
    out = qcio.CInOutputWriter(outpath)
    try:
        for ai in tqdm(atlist):
            ai.write(out)
    finally:
        out.close()

def _atoms_conform(dbfile, energy, force, virial):
    """Determines whether the specified database conforms to the constraints for
    the database.

    Args:
        dbfile (str): name of the database file to check.
        energy (str): name of the parameter that describes DFT/reference energy.
        force (str): name of the parameter that describes DFT/reference forces.
        virial (str): name of the parameter that describes DFT/reference virial
          tensor.

    Returns:
        tuple: `(params, force)`, where `params` is a dictionary where keys are
        target parameter names and values are source parameter names. If the
        dictionary has no values, then all the required parameters are already
        in the atoms object. `force` is a boolean that is True if `force` must
        be copied to `dft_force` to conform with `matdb` conventions.
    """
    emsg = "Cannot find {0} under parameter name {1}."
    params = {}
    doforce = False
    a = quippy.Atoms(dbfile)
    
    if energy not in a.params:
        raise ValueError(emsg.format("energy", energy))
    #These next two have unit tests and pytest says these are raised, but
    #coverage doesn't catch them for some reason...
    if force not in a.properties:# pragma: no cover
        raise ValueError(emsg.format("force", force))
    if virial not in a.params:# pragma: no cover
        raise ValueError(emsg.format("virial", virial))

    if energy != "dft_energy":
        params["dft_energy"] = energy
    if force != "dft_force":
        doforce = True
    if virial != "dft_virial":
        params["dft_virial"] = virial

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
                 folder=None, pattern=None, config_type=None, energy="dft_energy",
                 force="dft_force", virial="dft_virial", limit=None):
        self.name = name
        self.root = path.join(root, self.name)
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        self.controller = controller
        self.splits = {} if splits is None else splits
        self.folder = folder

        self._dbfile = path.join(self.root, "legacy-{}.xyz".format(limit))
        """str: path to the combined legacy database, with limits included.
        """
        self._dbfull = path.join(self.root, "legacy.xyz")
        """str: path to the combined legacy database, *without* limits.
        """
        self.dbfiles = []
        self.config_type = config_type

        from matdb.utility import dbconfig
        config = dbconfig(self._dbfull)
        if path.isfile(self._dbfile) and len(config) > 0:
            self.dbfiles = [db[0] for db in config["sources"]]
            self.config_type = config["config_type"]
            self.folder = folder
        else:
            from matdb.utility import dbcat
            if not path.isfile(self._dbfull):
                self._create_dbfull(folder, pattern, energy, force, virial, config_type)
            
            if limit is not None:
                msg.std("Slicing limit subset of full {} db.".format(self.name))
                full = quippy.AtomsList(self._dbfull)
                N = np.arange(len(full))
                np.random.shuffle(N)
                ids = N[0:limit]
                part = full[ids]
                part.write(self._dbfile)
                dbcat([self._dbfull], self._dbfile, docat=False, limit=limit,
                      ids=ids)
            else:
                from matdb.utility import symlink
                symlink(self._dbfile, self._dbfull)

        #The rest of matdb expects each database to have an atoms object that is
        #representative. Just take the first config in the combined database.
        self.atoms = quippy.Atoms(self._dbfile)

    def _create_dbfull(self, folder, pattern, energy, force, virial, config_type):
        """Creates the full combined database.
        """
        from matdb.utility import chdir, dbcat
        from glob import glob
        from tqdm import tqdm
        import quippy.cinoutput as qcio
        from os import path
        
        #NB! There is a subtle bug here: if you try and open a quippy.Atoms
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
                al = quippy.AtomsList(dbpath)
                outpath = path.join(self.root, dbfile)
                out = qcio.CInOutputWriter(outpath)
                try:
                    for ai in tqdm(al):
                        for target, source in params.items():
                            if (target == "config_type" and
                                  config_type is not None):
                                ai.params[target] = config_type
                            else:
                                ai.params[target] = ai.params[source]
                                del ai.params[source]

                        if doforce:
                            ai.properties["dft_force"] = ai.properties[force]
                            del ai.properties[force]

                        ai.write(out)

                    #Mark this db as non-conforming so that we created a new
                    #version of it.
                    rewrites.append(dbfile)
                finally:
                    out.close()

                dbcat([dbpath], outpath, docat=False, renames=params,
                      doforce=doforce)

        catdbs = []
        for dbfile in self.dbfiles:
            if dbfile in rewrites:
                catdbs.append(path.join(self.root, dbfile))
            else:
                catdbs.append(path.join(folder, dbfile))

        self.dbfiles = catdbs

        #Finally, concatenate all the files together to form the final
        #database.
        from matdb.utility import dbcat
        with chdir(folder):
            dbcat(self.dbfiles, self._dbfull, config_type=self.config_type)
        
    def train_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used for training.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{0}-train.xyz".format(split))

    def holdout_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used to validate the potential fit.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{0}-holdout.xyz".format(split))

    def super_file(self, split):
        """Returns the full path to the XYZ database file that can be
        used to *super* validate the potential fit.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.root, "{0}-super.xyz".format(split))
                
    def split(self, recalc=0):
        """Splits the database multiple times, one for each `split` setting in
        the database specification.
        """
        for name in self.splits:
            self._split(name, recalc)
    
    def _split(self, name, recalc=0):
        """Splits the total available data in the combined legacy file into a training
        and holdout sets.

        Args:
            name (str): name of the split to perform.
            recalc (int): when non-zero, re-split the data and overwrite any
              existing *.xyz files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level XYZ files, increase this value.
        """
        train_file = self.train_file(name)
        holdout_file = self.holdout_file(name)
        super_file = self.super_file(name)
        if (path.isfile(train_file) and path.isfile(holdout_file)
            and path.isfile(super_file) and recalc <= 0):
            return

        train_perc = self.splits[name]
        
        #Compile a list of all the sub-configurations we can include in the
        #training.
        from cPickle import dump, load
        from tqdm import tqdm
        idfile = path.join(self.root, "{0}-ids.pkl".format(name))

        #Either way, we will have to compile a list of all available atoms in
        #the database files.
        msg.info("Working on split {} for {}.".format(name, self.name))
        subconfs = quippy.AtomsList(self._dbfile)

        if path.isfile(idfile):
            with open(idfile, 'rb') as f:
                data = load(f)

            ids = data["ids"]
            Ntrain = data["Ntrain"]
            Nhold = data["Nhold"]
            Ntot = data["Ntot"]
            Nsuper = data["Nsuper"]
        else:
            Ntot = len(subconfs)
            Ntrain = int(np.ceil(Ntot*train_perc))
            ids = np.arange(Ntot)
            Nhold = int(np.ceil((Ntot-Ntrain)*train_perc))
            Nsuper = Ntot-Ntrain-Nhold
            np.random.shuffle(ids)

            #We need to save these ids so that we don't mess up the statistics on
            #the training and validation sets.
            data = {
                "ids": ids,
                "Ntrain": Ntrain,
                "Nhold": Nhold,
                "Ntot": Ntot,
                "Nsuper": Nsuper
            }
            with open(idfile, 'wb') as f:
                dump(data, f)

        #Only write the minimum necessary files. Use dbcat to create the
        #database version and configuration information. There is duplication
        #here because we also store the ids again. We retain the pkl file above
        #so that we can recreate *exactly* the same split again later.
        from matdb.utility import dbcat
        if not path.isfile(train_file):
            tids = ids[0:Ntrain]
            altrain = subconfs[tids]
            _quick_write(altrain, train_file)
            dbcat([self._dbfile], train_file, docat=False, ids=tids, N=Ntrain)
        if not path.isfile(holdout_file):
            hids = ids[Ntrain:-Nsuper]
            alhold = subconfs[hids]
            _quick_write(alhold, holdout_file)
            dbcat([self._dbfile], holdout_file, docat=False, ids=hids, N=Nhold)
        if not path.isfile(super_file):
            sids = ids[-Nsuper:]
            alsuper = subconfs[sids]
            _quick_write(alsuper, super_file)
            dbcat([self._dbfile], super_file, docat=False, ids=sids, N=Nsuper)
