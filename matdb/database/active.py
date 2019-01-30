"""Group of configurations that is created during an mtp active learning iteration.
"""
from hashlib import sha1
from glob import glob
from os import path, getcwd, chdir, remove, listdir, mkdir

import numpy as np
from six import string_types
from tqdm import tqdm

from matdb import msg
from matdb.database import Group
from matdb.atoms import AtomsList, Atoms

class Active(Group):
    """Sets up the calculations for a set of configurations that are being
    added to by the active learning approach.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        auids (list): the unique ids for each config in the group
    """
    def __init__(self, name="active", **dbargs):
        self.name = name
        dbargs['prefix'] = "Ac"
        dbargs['cls'] = Active
        dbargs['trainable'] = True
        if "calculator" in dbargs:
            self.calcargs = dbargs["calculator"].copy()
            del dbargs["calculator"]
        if "Active" not in dbargs['root']:
            new_root =path.join(dbargs['root'],"Active")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        super(Active, self).__init__(**dbargs)

        self.auids = None
        self._load_auids()
        self.nconfigs = len(self.auids) if self.auids is not None else 0
        self.last_iteration = None
        cur_iter = len(glob(path.join(self.root,"iter_*.pkl")))
        self.iter_file = path.join(self.root,"iter_{}.pkl".format(cur_iter))
        self._load_last_iter()

    def ready(self):
        """Returns True if this database has finished its computations and is
        ready to be used.
        """
        return len(self.fitting_configs) == self.nconfigs

    def can_extract(self):
        """Runs post-execution routines to clean-up the calculations. 
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            if (len(self.configs) != self.nconfigs and
                self.nconfigs is not None):
                #We need to have at least one folder for each config;
                #otherwise we aren't ready to go.
                return False

            result = False
            for f, a in zip(list(self.configs.keys()), list(self.config_atoms.keys())):
                if not self.config_atoms.get(a).calc.can_extract(self.configs.get(f)):
                    msg.std("Config {} not ready for extraction.".format(self.configs.get(f)), 2)
                    # remove the not-ready-for-extraction item from the dictionaries
                    self.configs.pop(f)
                    self.config_atoms.pop(a)
                    # continue processing the rest.
                    continue
                else:
                    result = True
            return result
        else:
            return all(group.can_extract() for group in self.sequence.values())

    @property
    def auid_file(self):
        """Returns the full path to the euid file for this group.
        """
        return path.join(self.root,"auids.pkl")

    def _load_auids(self):
        """Loads the list of `auid` from the `rset.pkl` file for this
        database group.
        """
        if self.auids is None:
            self.auids = self.load_pkl(self.auid_file)

    def _load_last_iter(self):
        """Loads the list of paths from the `iter_{}.pkl` file for this
        database group's last iteration.
        """
        self.last_iteration = self.load_pkl(self.iter_file)

    @property
    def last_config_atoms(self):
        """Returns the atoms objects from the last iteration.
        """

        last_atoms = {}
        last_count = 0
        if self.last_iteration is None or (self.last_iteration is not None and
                                           len(self.last_iteration) == 0):
            return None
        else:
            for i in range(len(self.config_atoms)-len(self.last_iteration),
                           len(self.config_atoms)):
                last_atoms[last_count] = self.config_atoms[list(self.config_atoms.keys())[i]]
                last_count += 1

        return last_atoms

    @property
    def fitting_configs(self):
        """Returns a list of full paths to the folders that have `atoms.json` objects
        for the full result set.
        """
        result = []
        for auid in self.auids:
            folder = self.index[str(auid)]
            target = path.join(folder,"atoms.h5")
            if path.isfile(target):
                result.append(target)

        return result

    @property
    def rset(self):
        """Returns a :class:`matdb.atoms.AtomsList`, one for each config in the
        latest result set.
        """

        #Return the configurations from this group; it is at the
        #bottom of the stack
        result = AtomsList()
        for epath in self.fitting_configs:
            result.append(Atoms(epath))
        return result

    def add_configs(self,new_configs,iteration):
        """Adds the atoms objects in the list to the configs of the active set.

        Args:
            new_configs (list): list of `matdb.atoms.Atoms` objects to be added
                to the active learning set.
        """

        self.new_configs = new_configs
        self.nconfigs += len(new_configs)
        self.iter_file = path.join(self.root,"iter_{}.pkl".format(iteration))
        self.last_iteration = {}

    def _setup_configs(self, rerun=False):
        """Sets up the database structure for the active set and creates a
        folder for the `calculator` to run in for each config.
        """
        # We need to make sure that none of the configs in this round
        # of the active learning have already been visited by
        # constructing their auids and then verifying that the auid
        # hasn't been visited before.
        dind = len(self.auids) if self.auids is not None else 0
        iter_ind = 0
        pbar = tqdm(total=len(self.new_configs))
        for config in self.new_configs:
            auid = sha1(''.join(map(str, (tuple([tuple(i) for i in config.cell]),
                        tuple([tuple(i) for i in config.positions]),
                        tuple(config.get_chemical_symbols())))).encode('utf-8'))
            if self.auids is not None:
                if auid.hexdigest() in self.auids:
                    self.nconfigs -= 1
                    continue
            else:
                self.auids = []

            dind += 1
            self.create(config,cid=dind)
            self.index[auid.hexdigest()] = self.configs[dind]
            self.last_iteration[iter_ind] = self.configs[dind]
            self.auids.append(auid.hexdigest())
            iter_ind += 1
            pbar.update(1)
        pbar.close()

        self.jobfile(rerun)
        self.save_index()
        self.save_pkl(self.auids,self.auid_file)
        self.save_pkl(self.last_iteration,self.iter_file)

    def setup(self, rerun=False):
        """Enumerates the desired number of structures and setups up a folder
        for each one.

        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist.

        """
        super(Active, self).setup(self._setup_configs,rerun)

    def is_executing(self):
        """Returns True if the most recently added configurations are in
        process of being executed.
        """

        is_executing = False
        if len(self.configs) != self.nconfigs:
            # There have been configurations added to the database
            # that haven't been setup yet.
            return False
        else:
            # We only want to know if the last iteration's files are
            # being executed, the old iterations don't matter.
            for config in self.last_iteration.values():
                try:
                    atoms = Atoms(path.join(config, "pre_comp_atoms.h5"))
                except:
                    atoms = Atoms(path.join(config, "atoms.h5"))
                is_executing = atoms.calc.is_executing(config)
                if is_executing:
                    break

        return is_executing
    
    def execute(self, dryrun=False, recovery=False, env_vars=None):
        """Submits the job script for each of the folders in this
        database if they are ready to run.
        Args:
            dryrun (bool): when True, simulate the submission without
              actually submitting.
            recovery (bool): when True, submit the script for running recovery
              jobs.
            env_vars (dict): of environment variables to set before calling the
              execution. The variables will be set back after execution.
        Returns:
            bool: True if the submission generated a job id
            (considered successful).
        """

        self._expand_sequence()
        if len(self.sequence) == 0:
            jobfile = "recovery.sh" if recovery else "jobfile.sh"
            if not path.isfile(path.join(self.root, jobfile)):
                return False

            if not recovery:
                if not all(a.calc.can_execute(self.configs[i])
                           for i, a in self.config_atoms.items()):
                    return False

                #We also need to check that we haven't already submitted this
                #job. Check to see if it is executing.
                if any(a.calc.is_executing(self.configs[i])
                       for i, a in self.config_atoms.items()):
                    return False

                #Make sure that the calculation isn't complete.
                if self.last_config_atoms is not None:
                    if any(a.calc.can_extract(self.last_iteration[i])
                           for i, a in self.last_config_atoms.items()):
                      return False

            # We must have what we need to execute. Compile the command and
            # submit.
            from matdb.utility import execute
            cargs = [self.database.parent.shell_command, jobfile]
            if dryrun:
                from matdb.msg import okay
                okay("Executed {} in {}".format(' '.join(cargs), self.root))
                return True
            else:
                xres = execute(cargs, self.root, env_vars=env_vars)

            if len(xres["output"]) > 0 and "Submitted" in xres["output"][0]:
                from matdb.msg import okay
                okay("{}: {}".format(self.root, xres["output"][0].strip()))
                return True
            else:
                return False

        else:
            already_executed = []
            for group in self.sequence.values():
                already_executed.append(group.execute(dryrun=dryrun,
                                                      recovery=recovery,
                                                      env_vars=env_vars))
            return all(already_executed)
