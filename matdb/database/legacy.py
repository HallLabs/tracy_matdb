"""Although `matdb` provides plenty of functionality for generating databases,
there are also many existing databases that have interesting data. Even though
we don't necessarily know exactly how they were made, the data can still be
useful. This module provides a simple class that adapts legacy databases to a
format that can be used by the `matdb` fitting machinery.
"""
from os import path
import quippy

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
    if force not in a.properties:
        raise ValueError(emsg.format("force", force))
    if virial not in a.params:
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
        parent (matdb.database.controller.Controller): instance controlling
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
        config (dict): key-value pairs describing how the `legacy.xyz` combined
          database file was created.
    """
    def __init__(self, name=None, root=None, controller=None, splits=None,
                 folder=None, pattern=None, config_type=None, energy="dft_energy",
                 force="dft_force", virial="dft_virial"):
        self.name = name
        self.root = path.join(root, self.name)
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        self.controller = controller
        self.splits = {} if splits is None else splits
        self.folder = folder

        self._dbfile = path.join(self.root, "legacy.xyz")
        """str: path to the combined legacy database.
        """
        self._configfile = path.join(self.root, "config.json")
        """str: path to the file that stores configuration information about how
        the :attr:`_dbfile` was created.
        """
        self.config = {}
        self.dbfiles = []
        self.config_type = config_type

        import json
        if path.isfile(self._dbfile) and path.isfile(self._configfile):
            with open(self._configfile) as f:
                self.config = json.load(f)

            self.dbfiles = self.config["dbfiles"]
            self.config_type = self.config["config_type"]
            self.folder = folder
        else:
            from matdb.utility import chdir
            from glob import glob
            from tqdm import tqdm
            import quippy.cinoutput as qcio
            
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
                    al = quippy.AtomsList(dbpath)
                    out = qcio.CInOutputWriter(path.join(self.root, dbfile))
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

            catdbs = []
            for dbfile in self.dbfiles:
                if dbfile in rewrites:
                    catdbs.append(path.join(self.root, dbfile))
                else:
                    catdbs.append(path.join(folder, dbfile))
                    
            self.config["dbfiles"] = catdbs
            self.config["config_type"] = self.config_type
            self.config["folder"] = folder

            #Finally, concatenate all the files together to form the final
            #database.
            from matdb.utility import cat
            with chdir(folder):
                cat(self.config["dbfiles"], self._dbfile)

            with open(self._configfile, 'w') as f:
                json.dump(self.config, f)

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
        idfile = path.join(self.root, "{0}-ids.pkl".format(name))
        if path.isfile(idfile):
            with open(idfile, 'rb') as f:
                data = load(f)
            subconfs = data["subconfs"]
            ids = data["ids"]
            Ntrain = data["Ntrain"]
            Nhold = data["Nhold"]
            Ntot = data["Ntot"]
            Nsuper = data["Nsuper"]
        else:
            Ntot = 0
            subconfs = {}
            for fi, dbfile in enumerate(self.dbfiles):
                subconfs[fi] = path.join(self.config["folder"], dbfile)
                        
            Ntot = len(subconfs)
            Ntrain = int(np.ceil(Ntot*train_perc))
            ids = np.arange(len(subconfs))
            Nhold = int(np.ceil((Ntot-Ntrain)*train_perc))
            Nsuper = Ntot-Ntrain-Nhold
            np.random.shuffle(ids)

            #We need to save these ids so that we don't mess up the statistics on
            #the training and validation sets.
            data = {
                "subconfs": subconfs,
                "ids": ids,
                "Ntrain": Ntrain,
                "Nhold": Nhold,
                "Ntot": Ntot,
                "Nsuper": Nsuper
            }
            with open(idfile, 'wb') as f:
                dump(data, f)
        
        tids = ids[0:Ntrain]
        hids = ids[Ntrain:-Nsuper]
        sids = ids[-Nsuper:]

        def subset(subconfs, idlist, recalc):
            from matdb.io import vasp_to_xyz
            from tqdm import tqdm
            files = []
            for aid in tqdm(idlist):
                if vasp_to_xyz(subconfs[aid], recalc=recalc-1):
                    files.append(path.join(subconfs[aid], "output.xyz"))

            return files

        arrconfs = np.array(subconfs)
        trainfiles = arrconfs[tids]
        holdfiles = arrconfs[hids]
        superfiles = arrconfs[sids]

        from matdb.utility import cat
        cat(trainfiles, train_file)
        cat(holdfiles, holdout_file)
        cat(superfiles, super_file)
