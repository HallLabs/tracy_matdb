"""Contains all the utility functions that belong to the database groups."""

from uuid import uuid4
from cPickle import dump, load
from os import path, rename, remove
import numpy as np
from glob import glob
from tqdm import tqdm

from matdb import msg
from matdb.atoms import AtomsList
# from matdb.utility import dbcat
# from matdb.utility import load_datetime
# from matdb.utility import special_values

def split(atlist, splits, targets, dbdir, ran_seed, dbfile=None, recalc=0,
          nonsplit=None):
    """Splits the :class:`matdb.atoms.AtomsList` multiple times, one for
    each `split` setting in the database specification.

    Args:
        atlsit (AtomsList, or list): the list of :class:`matdb.atams.Atoms` objects
          to be split or a list to the files containing the atoms objects.
        splits (dict): the splits to perform.
        targets (dict): the files to save the splits in, these should
          contain a {} in the name which will be replaced with the split
          name. The dictionary must have the format {"train": file_name,
          "holdout": file_name, "super": file_name}.
        dbdir (str): the root *splits* directory for the database.
        dbfile (str): the _dbfile for a legacy database.
        ran_seed (int or float): the random seed for the splits (i.e. the controllers
          random seed).
        recalc (int): when non-zero, re-split the data and overwrite any
          existing *.h5 files. This parameter decreases as
          rewrites proceed down the stack. To re-calculate
          lower-level h5 files, increase this value.
        nonsplit (AtomsList): a list of atoms to include in the training
          set "as-is" because they cannot be split (they only have meaning
          together).
    """
    from matdb.utility import dbcat

    assert nonsplit is None or isinstance(nonsplit, AtomsList)
    for name, train_perc in splits.items():
        train_file = targets["train"](name)
        holdout_file = targets["holdout"](name)
        super_file = targets["super"](name)
        idfile = path.join(dbdir, "{0}-ids.pkl".format(name))

        if (path.isfile(train_file) and path.isfile(holdout_file)
            and path.isfile(super_file)):
            if recalc <= 0:
                return
            else:
                if path.isfile(idfile):
                    with open(idfile, 'rb') as f:
                        data = load(f)
                new_idfile = path.join(self.root,"{0}_{1}-ids.pkl".format(name,data["uuid"]))
                self.parent.uuids[data["uuid"]] = new_idfile
                for fname in [train_file,holdout_file,super_file]:
                    new_name = fname.replace(name,"{0}_{1}".format(name,data["uuid"]))
                    rename(fname,new_name)
                remove(idfile)

        #Compile a list of all the sub-configurations we can include in the
        #training.
        if not isinstance(atlist,AtomsList):
            subconfs = AtomsList(atlist)
        else:
            subconfs = atlist

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
            Ntot = len(subconfs)
            Ntrain = int(np.ceil(Ntot*train_perc))
            ids = np.arange(Ntot)
            Nhold = int(np.ceil((Ntot-Ntrain)*train_perc))
            Nsuper = Ntot-Ntrain-Nhold
            np.random.shuffle(ids)

            #We need to save these ids so that we don't mess up the statistics on
            #the training and validation sets.
            data = {
                "uuid": str(uuid4()),
                "subconfs": subconfs,
                "ids": ids,
                "Ntrain": Ntrain,
                "Nhold": Nhold,
                "Ntot": Ntot,
                "Nsuper": Nsuper,
                "ran_seed": ran_seed
            }
            with open(idfile, 'wb') as f:
                dump(data, f)

        #Only write the minimum necessary files. Use dbcat to create the
        #database version and configuration information. There is duplication
        #here because we also store the ids again. We retain the pkl file above
        #so that we can recreate *exactly* the same split again later.
        if not path.isfile(train_file):
            tids = ids[0:Ntrain]
            #Make sure that we have some atoms to write in the first place!
            if len(tids) > 0:
                altrain = subconfs[tids]
            else:
                altrain = AtomsList()
            #Add the unsplittable configurations to the training set as-is.
            Nunsplit = 0
            if nonsplit is not None:
                altrain.extend(nonsplit)
                Nunsplit = len(nonsplit)
            altrain.write(train_file)

            if dbfile is not None:
                dbcat([dbfile], train_file, docat=False, ids=tids, N=Ntrain+Nunsplit)
            else:
                dbcat([], train_file, docat=False, ids=tids, N=Ntrain+Nunsplit)
        if not path.isfile(holdout_file):
            hids = ids[Ntrain:-Nsuper]
            alhold = subconfs[hids]
            alhold.write(holdout_file)
            if dbfile is not None:
                dbcat([dbfile], holdout_file, docat=False, ids=hids, N=Nhold)
            else:
                dbcat([], holdout_file, docat=False, ids=hids, N=Nhold)
        if not path.isfile(super_file):
            sids = ids[-Nsuper:]
            alsuper = subconfs[sids]
            alsuper.write(super_file)
            if dbfile is not None:
                dbcat([dbfile], super_file, docat=False, ids=sids, N=Nsuper)
            else:
                dbcat([], super_file, docat=False, ids=sids, N=Nsuper)

def dbconfig(dbfile):
    """Returns the database configuration `dict` of the specified database file.

    Args:
        dbfile (str): path to the database file to get config information for.
    """
    from matdb.utility import load_datetime

    confpath = dbfile + ".json"
    if not path.isfile(confpath):
        return {}

    import json
    with open(confpath) as f:
        config = json.load(f, object_pairs_hook=load_datetime)

    return config


def parse_path(root,seeds,ran_seed=None):
    """Finds the full path to the seed files for this system.
    Args:
        root (str): the root directory for the databes.
        seeds (str or list of str): the seeds for the database that
            need to be parsed and have the root folder found.
        ran_seed (optional): the seed for the random generator if used.

    Returns:
        seed_files (list): a list of the seed files for the database.
    """
    from matdb.utility import special_values
    from itertools import product

    seed_files = []
    for seed in seeds:
        # if there is a '/' in the seed then this is a path to the
        # seed file configuration, otherwise we assume the user put
        # the seed file in the `seed` directory in root.
        if len(seed.split("/")) > 1:
            this_seeds = []
            seed_path = root
            for segment in seed.split("/"):
                if "*" in segment:
                    if len(this_seeds) >=1:
                        this_level = []
                        for t_path in res:
                            this_level.extend(glob(path.join(seed_path,t_path,segment)))
                    else:
                        this_level = glob(path.join(seed_path,segment))
                else:
                    this_level = [segment]
                if len(this_seeds) >= 1:
                    this_seeds.extend([path.join(*i) for i in product(this_seeds,this_level)])
                else:
                    this_seeds.extend(this_level)

        else:
            seed_path = path.join(root,"seed")
            to_parse = seed
            if "*" in to_parse:
                this_seeds = glob(path.join(seed_path,to_parse))
            else:
                this_seeds = [to_parse]

        for ts in this_seeds:
            t_seed = path.join(seed_path,ts)
            if path.isfile(t_seed):
                seed_files.append(t_seed)
            else:
                msg.err("The seed file {} could not be found.".format(t_seed))

        return seed_files
