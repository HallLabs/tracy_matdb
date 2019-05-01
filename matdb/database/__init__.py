"""Exposes classes and functions for interacting with the database
folders via a simple configuration file.
"""
import abc
import collections
from contextlib import contextmanager
from datetime import datetime
from glob import glob
import pickle
from os import path, mkdir, makedirs, sys, rename, remove
from uuid import uuid4

import ase.db
from collections import OrderedDict
from hashlib import sha1
import h5py
import json
import numpy as np
from importlib import import_module
import re
import six
from tqdm import tqdm

from matdb import __version__, msg
from matdb.atoms import Atoms, AtomsList, _recursively_convert_units
from matdb.database.legacy import LegacyDatabase
from matdb.database.utility import parse_path, split
from matdb.fitting.controller import TController
from matdb.io import read, save_dict_to_h5
from matdb.msg import okay, verbosity
from matdb.utility import (chdir, ParameterGrid, convert_dict_to_str,
                            import_fqdn, is_uuid4, _set_config_paths)
from matdb import calculators

class Group(object):
    """Represents a collection of material configurations (varying in
    structure and composition) from which a machine learning model can
    be created. Includes logic for generating the DFT directories that
    need to be run as well as extracting the relevant data from such
    computations.
    Args:
        execution (dict): key-value pairs of settings for the supercomputer job
          array batch file.
        root (str): full path to the root directory that this database will live
          in.
        parent (matdb.database.controller.Database): the database that this
          group of calculations belong to.
        prefix (str): sub-sampled configurations will be stored using integer
          ids after this prefix; for example `S.1`, `S.2`, etc.
        nconfigs (int): number of displaced configurations to create.
        config_type (str): the type of configuration.
        calculator (dict): a dictionary containing the information for
          the calculator object.
        trainable (bool): True if the groups configs will be used for traning.
        pgrid (ParamaterGrid): The ParameterGrid for the database.
        seeds (list, str, matdb.atoms.Atoms): The location of the files that will be
          read into to make the atoms object or an atoms object.
        cls (subclass): the subclass of :class:`Group`.
        # override (dict): a dictionary of with uuids or paths as the
        #   keys and a dictionary containing parameter: value pairs for
        #   parameters that need to be adjusted.
        transforms (dict): a dictionary of transformations to apply to the
          seeds of the database before calculations are performed. Format is
          {"name": {"args": dict of keyword args}}, where the "name" keyword
          is the fully qualified path to the function.

    Attributes:
        atoms (matdb.atoms.Atoms): a single atomic configuration from
          which many others may be derived using MD, phonon
          displacements, etc.
        configs (dict): keys are integer identifiers for the particular
          configuration; values are paths (relative to the base atoms root
          directory) in which calculations are performed.
        root (str): full path to the root directory that this database will live
          in.
        database (matdb.database.controller.Database): parent database
          to which this group belongs.
        prefix (str): sub-sampled configurations will be stored using integer
          ids after this prefix; for example `S.1`, `S.2`, etc.
        nconfigs (int): number of displaced configurations to create.
        pgrid (ParamaterGrid): The ParameterGrid for the database.
        grpargs (dict): default arguments to construct the new groups; will
          be overridden by any parameter grid specs.
        splittable (bool): when True, this Group can be split into training,
          holdout and super sets because its configs are independent; otherwise,
          the configs are kept together and used *only* in the training set.
    """
    seeded = False
    splittable = True

    def __init__(self, cls=None, root=None, parent=None, prefix='S', pgrid=None,
                 nconfigs=None, calculator=None, seeds=None,
                 config_type=None, execution=None, trainable=False, #override=None,
                 rec_bin=None, transforms = None):
        if isinstance(parent, Database):
            #Because we allow the user to override the name of the group, we
            #have to expand our root to use that name over here. However, for
            #recursively nested groups, we use the root passed in because it
            #already includes the relevant suffixes, etc.
            self.root = path.join(root, "{0}.{1}".format(parent.name, self.name))
        else:
            self.root = root

        if not path.isdir(self.root):
            mkdir(self.root)

        self.index = {}
        self._read_index()

        self.cls = cls
        self.parent = parent
        self.execution = execution if execution is not None else {}
        self.atoms = None
        self.transforms = transforms if transforms is not None else {}

        self._trainable = trainable
        self.is_seed_listed = None
        if seeds is not None:
            self.is_seed_listed = isinstance(seeds, (list, six.string_types))

        self.rec_bin = rec_bin
        self.seeds = None
        self.pgrid = pgrid
        self._seed = seeds
        """list, str, matdb.atoms.Atoms: The location of the files that will be
        read into to make the atoms object or an atoms object. This is the
        parameter that was passed as the seed for the group constructor; for
        recursive constructions, it starts as a list of seed patterns, which
        gets expanded into individual seeds, which may then be coupled to
        parameter grid specs.
        """
        self.grpargs = dict(parent=self, prefix=prefix, nconfigs=nconfigs,
                            trainable=trainable, execution=execution,
                            config_type=config_type, calculator=calculator)

        self.sequence = OrderedDict()
        self.calc = None
        self.calcargs = self.database.calculator.copy()
        if calculator is not None:
            self.calcargs.update(calculator)
        if "name" in self.calcargs:
            self.calc = getattr(calculators, self.calcargs["name"])

        self.prefix = prefix
        self.nconfigs = nconfigs
        self.config_type = config_type

        self._nsuccess = 0
        """int: number of configurations whose output files were successfully
        converted to XYZ format. Should be equal to :attr:`nconfigs` if the
        database is complete.
        """
        self._db_name = self.database.name

        #Try and load existing folders that match the prefix into the configs
        #list.
        self.configs = {}
        self.config_atoms = {}
        self._rx_folder = re.compile(r"{}\.\d+".format(prefix))

        from os import getcwd
        with chdir(self.root):
            for folder in glob("{}.*".format(prefix)):
                #We aren't interested in files, or folders that just happen to match
                #the naming convention.
                if not path.isdir(folder):
                    continue
                if not self._rx_folder.match(folder):
                    continue

                cid = int(folder.split('.')[1])
                self.configs[cid] = path.join(self.root, folder)
                if path.isfile(path.join(folder,"atoms.h5")):
                    self.config_atoms[cid] = Atoms(path.join(self.configs[cid],"atoms.h5"))
                elif path.isfile(path.join(folder,"pre_comp_atoms.h5")):
                    self.config_atoms[cid] = Atoms(path.join(self.configs[cid],"pre_comp_atoms.h5"))
                else:
                    msg.warn("No config atoms available for {}.".format(self.configs[cid]))

        if path.isfile(path.join(
                self.root,"{}_{}_uuid.txt".format(self._db_name,self.prefix))):
            with open(path.join(
                    self.root,"{}_{}_uuid.txt".format(self._db_name,self.prefix)),"r") as f:
                uid = f.readline().strip()
                time_stamp = f.readline().strip()
        else:
            uid = str(uuid4())
            time_stamp = str(datetime.now())
            with open(path.join(
                    self.root,"{}_{}_uuid.txt".format(self._db_name,self.prefix)),"w") as f:
                f.write("{0} \n {1}".format(str(uid),str(time_stamp)))

        self.time_stamp = time_stamp
        self.uuid = uid
        self.database.parent.uuids[str(self.uuid)] = self

        # TODO: This entire code chunk needs to be moved to a class
        # method. It doesn't run properly as part of the class setup
        # because not all of the objects it searches for have been
        # loaded yet.
        
        # self.override = override.copy() if override is not None else {}
        # if bool(override):
        #     for k,v in override:
        #         obj_ins = self.database.controller.find(k)
        #         if isinstance(obj_ins,Atoms):
        #             if "calc" in v:
        #                 if self.rec_bin is not None:
        #                     #construct the path for the object in the recycling bin
        #                     new_path = path.join(self.rec_bin.root,
        #                                          "{}-atoms.h5".format(obj_ins.uuid))
        #                     obj_ins.write(new_path)
        #                     self.database.parent.uuids[obj_ins.uuid] = new_path
        #                     if path.isfile(path.join(obj_ins.calc.folder,"atoms.h5")):
        #                         # from os import remove
        #                         remove(path.join(obj_ins.calc.folder,"atoms.h5"))
        #                     # Make a new uuid for the new atoms object
        #                     # and overwrite the uuid file.
        #                     obj_ins.uuid = str(uuid4())
        #                     obj_ins.time_stamp = str(datetime.now())
        #                     with open(path.join(
        #                             obj_ins.root,
        #                             "{}_{}_uuid.txt".format(obj_ins._db_name,
        #                                                     obj_ins.prefix)),"w+") as f:
        #                         f.write("{0} \n {1}".format(obj_ins.uuid,obj_ins.time_stamp))
        #                     self.database.parent.uuids[obj_ins.uuid] = obj_ins

        #                 lcargs = None
        #                 if "name" in v["calc"]:
        #                     new_calc = getattr(calculators, v["calc"]["name"])
        #                     if "Tracy" in v["calc"]["name"]:
        #                         lcargs = self._tracy_setup(calcargs = v["calc"]["calcargs"])
        #                 else:
        #                     new_calc = obj_ins.calc

        #                 if lcargs is None:
        #                     lcargs = self.calcargs.copy()
        #                     lcargs.update(v["calc"]["calcargs"])
        #                     del lcargs["name"]

        #                 calc = new_calc(atoms, obj_ins.calc.folder, obj_ins.calc.contr_dir,
        #                                 obj_ins.calc.ran_seed, **lcargs)
        #                 obj_ins.set_calculator(calc)
        #         elif path.isfile(obj_ins):
        #             #If the object is a file path then it's pointing
        #             #to an old instance of a class objcet that has
        #             #been saved to file.
        #             msg.warn("Can't update object, it has already been overwritten.")
        #         else:
        #             args = obj_ins.to_dict()
        #             if self.rec_bin is not None:
        #                 for atm in self.fitting_configs():
        #                     # from os import rename
        #                     atms = Atoms(atm)
        #                     new_atm = path.join(selg.rec_bin.root,
        #                                         "{}-atoms.h5".format(atms.uuid))
        #                     self.database.parent.uuids[atms.uuid] = new_atm
        #                     rename(atm,new_atm)
        #             obj_ins.save_pkl(obj_ins.to_dict(),"{}.pkl".format(self.uuid))
        #             self.database.parent.uuids[obj_ins.uuid] = path.join(obj_ins.root,
        #                                                                  "{}.pkl".format(
        #                                                                      self.uuid))
        #             args.update(v)
        #             obj_ins.__init__(**args)
        #             obj_ins.uuid = str(uuid4())
        #             obj_ins.time_stamp = str(datetime.now())
        #             with open(path.join(obj_ins.root,
        #                                 "{}_{}_uuid.txt".format(obj_ins._db_name,
        #                                                         obj_ins.prefix)),"w+") as f:
        #                 f.write("{0} \n {1}".format(obj_ins.uuid,obj_ins.time_stamp))
        #             self.database.parent.uuids[obj_ins.uuid] = obj_ins
        #             obj_ins.setup(rerun=1)

    @property
    def key(self):
        """Returns the directory name for the group (i.e.,
        :func:`os.path.basename`).
        """
        return path.basename(self.root)

    @property
    def calculator(self):
        """Returns a representative calculator for the group. This will be the calculator
        attached to one of the :attr:`config_atoms` in the expanded group.
        """
        if len(self.config_atoms) == 0:
            self._expand_sequence()

        result = None
        if len(self.sequence) > 0:
            result = next(iter(self.sequence.values())).calculator
        else:
            try:
                result = next(iter(self.config_atoms.values())).get_calculator()
            except: #pragma: no cover
                raise
        return result

    @property
    def iconfigs(self):
        """Provides a generator over all the configurations in this group and its
        children.
        """
        self._expand_sequence()
        if len(self.sequence) > 0:
            for group in self.sequence.values():
                for atoms in group.iconfigs:
                    yield atoms
        else:
            for atoms in self.config_atoms.values():
                yield atoms

    def _expand_sequence(self):
        """Recursively expands the nested groups to populate :attr:`sequence`.
        """
        self._expand_seeds(self._seed)
        if self.seeds is not None:
            for seedname, at_seed in self.seeds.items():
                seed_root = path.join(self.root, seedname)
                if not path.isdir(seed_root):
                    mkdir(seed_root)

                clsargs = self.grpargs.copy()
                clsargs["pgrid"] = self.pgrid
                if self.pgrid is None or len(self.pgrid) == 0:
                    clsargs.update(self.pgrid.params)
                clsargs["root"] = seed_root
                clsargs["seeds"] = at_seed
                if isinstance(at_seed,Atoms):
                    self.atoms = at_seed
                if self.cls is None: #pragma: no cover
                    msg.err("The Group must have a class to have seeds.")
                self.sequence[seedname] = self.cls(**clsargs)
        else:
            if self.pgrid is not None and len(self.pgrid) > 0:
                for pkey in self.pgrid:
                    this_root = path.join(self.root, pkey)
                    if not path.isdir(this_root):
                        mkdir(this_root)

                    clsargs = self.grpargs.copy()
                    clsargs.update(self.pgrid[pkey])
                    clsargs["root"] = this_root
                    clsargs["seeds"] = self._seed
                    if isinstance(self._seed,Atoms):
                        self.atoms = self._seed
                    if self.cls is None: #pragma: no cover
                        msg.err("The Group must have a class to have a parameter grid.")
                    self.sequence[pkey] = self.cls(**clsargs)
            else:
                self.atoms = self._seed

    def _expand_seeds(self, seeds):
        """Expands explicitly listed seed wildcard patterns to populate the
        :attr:`seeds` dict.
        Args:
            seeds (list, str, matdb.atoms.Atoms): The location of the files that will be
              read into to make the atoms object or an atoms object.
        """
        if self.is_seed_listed:
            self.seeds = OrderedDict()
            for atomspec in seeds:
                fmt, pattern = atomspec.split(':')
                if fmt == "POSCAR":
                    fmt = "vasp"
                for apath in self.database.parent.relpaths([pattern]):
                    self.seeds[path.basename(apath)] = Atoms(apath, format=fmt)

        # At this point, not sure how to setup a multi-steps database with more than 2 steps
        # has the same name(this is required because of the prev property). 
        # So, keep this uncovered for the unit test.
        elif seeds is None and self.seeded and self.prev is not None: #pragma: no cover
            self.seeds = OrderedDict()
            for n_seeds, a in enumerate(self.prev.rset):
                seedname = "seed-{}".format(n_seeds)
                #NB! The previous rset may be a dict with an "atoms" key and
                #additional keys to pass to the group constructor. We copy it to
                #make sure the original rset doesn't modify. For normal cases
                #where it is simply an Atoms object, the copy performs the same
                #function.
                self.seeds[seedname] = a.copy()

    @property
    def database(self):
        """Returns the parent :class:`matdb.database.controller.Database` instance for
        this group, irrespective of how deep it is in the recursive stack.
        """
        if isinstance(self.parent, Database):
            return self.parent
        elif isinstance(self.parent, Group):
            return self.parent.database
        else: #pragma: no cover
            return None

    @property
    def trainable(self):
        """Determines if the group configs should be used for training.
        """
        if self.dep is None:
            return True
        else:
            return self._trainable

    @abc.abstractproperty
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this group.
        """
        pass #pragma: no cover

    @abc.abstractmethod
    def sub_dict(self):
        """Returns a dictionary of the parameters passed in to create this
        group.
        """
        pass #pragma: no cover

    def to_dict(self,include_time_stamp=True):
        """Returns a dictionary of the parameters passed into the group instance.
        """
        # from matdb import __version__
        # import sys

        kw_dict = self.grpargs.copy()
        args_dict = {"root": self.root, "version":__version__,
                     "python_version":sys.version}
        if include_time_stamp:
            args_dict["datetime"] = str(datetime.now())

        kw_dict.update(args_dict)
        if "parent" in kw_dict:
            del kw_dict["parent"]

        kw_dict.update(self.sub_dict())
        return kw_dict

    @property
    def rset_file(self):
        """Returns the full path to the `rset.h5` file for this group.
        """
        return path.join(self.root, "rset.h5")

    @abc.abstractproperty
    def rset(self):
        """Saves the rset for the group and all sequences of the group.
        """
        pass #pragma: no cover

    def hash_group(self):
        """Hashes the rset and the parameters of the  for the group.
        """
        # from matdb.utility import convert_dict_to_str
        hash_str = convert_dict_to_str(self.to_dict(include_time_stamp=False))
        for atom in self.rset:
            temp_atom = Atoms(atom)
            hash_str += convert_dict_to_str(temp_atom.to_dict())

        return str(sha1(hash_str.encode()).hexdigest())

    def load_pkl(self, file_name, rel_path=None):
        """Loads a pickled obj from the specified location on the path.

        Args:
            file_name (str): the file name to be save too.
            rel_path (str): the relative path from self.root that the file will
              be saved to.
        """
        f_path = path.join(self.root, rel_path, file_name) \
                 if rel_path is not None else path.join(self.root,file_name)

        #Initialze the result as an empty list instead of None
        result = []
        #If the pkl file is not exist or is empty, don't bother to load it
        if path.isfile(f_path) and path.getsize(f_path) > 0:
            with open(f_path,"rb") as f:
                result = pickle.load(f)

        return result

    def save_pkl(self, obj, file_name, rel_path=None):
        """Saves the obj passed to the correct location on the path.

        Args:
            obj (dict): The dictionary to be written to file.
            file_name (str): the file name to be save too.
            rel_path (str): the relative path from self.root that the file will
              be saved to.
        """
        f_path = path.join(self.root, rel_path, file_name) \
                 if rel_path is not None else path.join(self.root,file_name)

        with open(f_path,"wb+") as f:
            pickle.dump(obj,f)

    def save_index(self):
        """Writes the unique index for each of the configs to file along with
        the relative path to the atoms.json file
        """
        with open(path.join(self.root,"index.json"),"w+") as f:
            json.dump(self.index,f)

    def _read_index(self):
        """Reads in the index from the index.json file if it exists.
        """
        if path.isfile(path.join(self.root,"index.json")):
            with open(path.join(self.root,"index.json"),"r") as f:
                self.index = json.load(f)

    @property
    def prev(self):
        """Finds the previous group in the database.
        """
        keylist = list(self.database.steps.keys())
        for i, v in enumerate(keylist):
            if v == self.name and i!=0:
                return self.database.steps[keylist[i-1]]

    @property
    def dep(self):
        """Finds the next, or dependent, group in the databes.
        """
        keylist = list(self.database.steps.keys())
        for i, v in enumerate(keylist):
            if v == self.name and i!=(len(keylist)-1):
                return self.database.steps[keylist[i+1]]

    def is_executing(self):
        """Returns True if the database DFT calculations are in process of being
        executed.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            is_executing = False
            for i, atoms in self.config_atoms.items():
                is_executing = atoms.calc.is_executing(self.configs[i])
                if is_executing:
                    break
        else:
            executing = [group.is_executing() for group in self.sequence.values()]
            is_executing = all(executing)

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
                if any(a.calc.can_extract(self.configs[i])
                       for i, a in self.config_atoms.items()):
                    return False

            # We must have what we need to execute. Compile the command and
            # submit.
            from matdb.utility import execute

            shell_command = self.database.parent.shell_command
            # We suport 'bash' and 'sbatch' shell commands, if it's neighter one 
            # of them, default to 'bash' 
            if shell_command not in ['bash', 'sbatch']:
                shell_command = 'bash' 
            cargs = [shell_command, jobfile]

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

    def recover(self, rerun=0):
        """Compiles a list of all DFT runs that didn't complete and compiles the
        `failures` file. Creates a jobfile to re-run the failed
        folders only.

        Args:
            rerun (int): when > 0, recreate the jobfile even if it
              already exists.
        """

        self._expand_sequence()

        if len(self.sequence) == 0:
            detail = self.status(False)
            failed = [k for k, v in detail["done"].items() if not v]
            identity = "{0}|{1}".format(self._db_name, self.name)
            xpath = path.join(self.root, "failures")

            if len(failed) > 0:
                #Only write a failures file if we had failures.
                with open(xpath, 'w') as f:
                    f.write('\n'.join(failed))

                imsg = "{0}: queued {1:d} configs for recovery."
                msg.info(imsg.format(identity, len(failed)))
            else:
                msg.okay("{0}: no failures.".format(identity))

            #Only create a jobfile if there were actually failures
            if len(failed) > 0:
                self.jobfile(rerun, recovery=True)
            else:
                #Delete any existing recovery files from previous failures.
                from os import remove
                jobfile = path.join(self.root, "recovery.sh")
                if path.isfile(jobfile):
                    remove(jobfile)
                if path.isfile(xpath):
                    remove(xpath)
        else:
            for group in self.sequence.values():
                group.recover(rerun=rerun)

    def jobfile(self, rerun=0, recovery=False):
        """Creates the job array file to run each of the sub-configurations in
        this database.
        Args:
            rerun (int): when > 0, recreate the jobfile even if it
              already exists.
            recovery (bool): when True, configure the jobfile to run
              recovery jobs for those that have previously failed. This uses a
              different template and execution path.
        """

        self._expand_sequence()

        if len(self.sequence) == 0:
            if recovery:
                from matdb.utility import linecount
                target = path.join(self.root, "recovery.sh")
                xpath = path.join(self.root, "failures")
                asize = linecount(xpath)
            else:
                target = path.join(self.root, "jobfile.sh")
                xpath = path.join(self.root, "{}.".format(self.prefix))
                asize = len(self.configs)

            if (path.isfile(target) and rerun == 0) or asize == 0:
                return

            # We use the global execution parameters and then any updates
            # locally. We need to add the execution directory (including prefix) and
            # the number of jobs in the array.
            settings = self.database.execution.copy()
            settings.update(self.execution.items())

            settings["execution_path"] = xpath
            settings["array_size"] = asize

            if "array_limit" in settings and asize < settings["array_limit"]:
                del settings["array_limit"]

            from jinja2 import Environment, PackageLoader
            env = Environment(loader=PackageLoader('matdb', 'templates'))
            if recovery:
                template = env.get_template(settings["template"].replace("array", "recovery"))
            else:
                template = env.get_template(settings["template"])

            # By default, don't do timeout on QE calculation
            if "exec_time_out_minutes" not in settings:
                settings["exec_time_out_minutes"] = 0

            with open(target, 'w') as f:
                f.write(template.render(**settings))
        else:
            for group in self.sequence.values():
                group.jobfile(rerun=rerun, recovery=recovery)

    # def _tracy_setup(self, calcargs=None):
    #     """Extracts the needed information from the group that needs to be
    #     passed to the Tracy calculator.

    #     Args:
    #         calcargs (dict): the updates to the global calculation arguments.
    #     """
    #     calcinput = self.calcargs.copy()
    #     if calcargs is not None:
    #         calcinput.update(calcargs)
    #     del calcinput["name"]

    #     tracy = {}
    #     exec_settings = self.Database.execution.copy()
    #     if self.prev is not None and self.seeded:
    #         tracy["group_preds"] = self.prev.uuid

    #     if "eCommerce" in exec_settings:
    #         tracy["ecommerce"] = exec_settings["eCommerce"]

    #     if "contract_predecessors" in exec_settings:
    #         tracy["contract_preds"] = exec_settings["contract_predecessors"]

    #     if "priority" in exec_settings:
    #         tracy["contract_priority"] = exec_settings["priority"]

    #     keys = ["time", "flops", "minimum_ram", "minimum_mem", "ncores", "network_latency",
    #             "role"]
    #     for key in keys:
    #         if key not in exec_settings.keys():
    #             raise ValueError("{0} must be set by the user.")

    #     tracy["max_time"] = int(exec_settings["time"])*360
    #     tracy["min_flops"] = int(exec_settings["flops"])
    #     tracy["min_ram"] = int(exec_settings["minimum_ram"])
    #     tracy["min_mem"] = int(exec_settings["minimum_mem"])
    #     tracy["ncores"] = int(exec_settings["ncores"])
    #     tracy["max_net_lat"] = int(exec_settings["network_latency"])
    #     tracy["role"] = exec_settings["role"]
    #     if "notifications" in exec_settings:
    #         tracy["notifications"] = exec_settings["notifications"]

    #     return {"calcargs": calcinput, "tracy": tracy}

    def create(self, atoms, cid=None, rewrite=False, sort=None, calcargs=None,
               extractable=True):

        """Creates a folder within this group to calculate properties.
        Args:
            atoms (matdb.atoms.Atoms): atomic configuration to run.
            cid (int): integer configuration id; if not specified, defaults to
              the next available integer.
            rewrite (bool): when True, overwrite any existing files with the
              latest settings.
            sort (bool): when True, sort the atoms by type so that
              supercell writes work correctly.
            calcargs (dict): additional config-specific arguments for the
              calculator that will be created and attached to the atoms object.
            extractable (bool): True if the creation needs to write files for
              later calculation and extraction. If False then we just write
              an atoms.h5 object for later use as a seed for another group.

        Returns:
            int: new integer configuration id if one was auto-assigned.
        """
        from matdb.utility import import_fqdn

        if len(self.sequence)==0:
            if cid is None:
                cid = len(self.configs) + 1

            from os import path, mkdir
            target = path.join(self.root, "{}.{}".format(self.prefix, cid))
            if not path.isdir(target):
                mkdir(target)

            if path.isfile(path.join(target,"uuid.txt")):
                with open(path.join(target,"uuid.txt"),"r") as f:
                    uid = f.readline().strip()
                    time_stamp = f.readline().strip()
            else:
                uid = str(uuid4())
                time_stamp = str(datetime.now())
                with open(path.join(target,"uuid.txt"),"w+") as f:
                    f.write("{0} \n {1}".format(uid,time_stamp))

            if path.isfile(path.join(target,"atoms.h5")) and rewrite:
                from os import rename
                back_up_path = back_up_path.replace('/','.')+'.atoms.h5'
                new_path = path.join(self.rec_bin.root,"{}-atoms.h5".format(uid))
                rename(path.join(target,'atoms.h5'),new_path)
                self.database.parent.uuids[uid] = new_path
                uid = str(uuid4())

            trans_atoms = Atoms()
            trans_atoms.copy_from(atoms)
            for name, func_args in self.transforms.items():
                modojb, func = import_fqdn(name)
                trans_atoms = func(trans_atoms, **func_args)

            trans_atoms.uuid = uid
            trans_atoms.time_stamp = time_stamp
            trans_atoms.group_uuid = self.uuid
            if "Tracy" in self.calcargs["name"]: #pragma: no cover
                lcargs = self._tracy_setup(calcargs)
            else:
                lcargs = self.calcargs.copy()
                del lcargs["name"]
                if calcargs is not None:
                    lcargs.update(calcargs)

            calc = self.calc(trans_atoms, target, self.database.parent.root,
                             self.database.parent.ran_seed, **lcargs)
            calc.create()
            trans_atoms.set_calculator(calc)
            # Turns out we need this for some seedless configurations.
            if extractable:
                trans_atoms.write(path.join(target,"pre_comp_atoms.h5"))
            else:
                trans_atoms.write(path.join(target,"atoms.h5"))
            #Finally, store the configuration for this folder.
            self.configs[cid] = target
            self.config_atoms[cid] = trans_atoms

            self.database.parent.uuids[uid] = trans_atoms
            return cid
        else:
            for group in self.sequence.values():
                group.create(atoms, cid=cid, rewrite=rewrite, sort=sort)

    def ready(self):
        """Determines if this database has been completely initialized *and* has
        all its computations' results ready.
        .. note:: This method should be overloaded by a sub-class.
        Raises:
            NotImplementedError: this method is intended to be overloaded by a
            sub-class.
        """
        raise NotImplementedError("Method `ready` must be overloaded by a "
                                  "sub-class.")

    def is_setup(self):
        """Determines if all the necessary folders for sub-configurations of the seed
        atomic configuration exist.
        """
        self._expand_sequence()

        if len(self.sequence) == 0:
            #Test to see if we have already set the database up.
            confok = False
            #If there is zero new configs, we still need to create an empty iter_*.pkl 
            #file in the database folder
            if self.nconfigs != 0 and (len(self.configs) == self.nconfigs or
                len(self.configs) > 0 and self.nconfigs is None):
                imsg = "The {} database has already been setup."
                msg.info(imsg.format(self.name), 2)
                confok = True

            #Don't run setup if the program is currently executing.
            xok = False
            if self.is_executing():
                xok = True

            result = confok or xok
        else:
            already_setup = [group.is_setup() for group in self.sequence.values()]
            result = all(already_setup)

        return result

    def setup(self, db_setup, rerun=0):
        """Performs the setup of the database using the `db_setup` function
        passed in by the specific group instance.

        Args:
            db_setup (function): a function that will perform the setup for each
                group independently.
            rerun (int): values > 0 cause rerunning of the setup at different granularity.
        """
        if self.prev is None or self.prev.can_extract():
            #Before we attempt to setup the folders, we first need to construct
            #the recursive sequence groups. This cannot happen until the
            #previous group in the database is ready, which is why it happens
            #here rather than in __init__.
            self._expand_sequence()
            if len(self.sequence) == 0:
                ok = self.is_setup()
                if ok and rerun == 0:
                    return
                db_setup(rerun)
                with open(path.join(self.root,"compute.pkl"),"wb+") as f:
                    pickle.dump({"date": datetime.now(),"uuid":self.uuid},f)
            else:
                pbar = tqdm(total=len(self.sequence))
                for group in self.sequence.values():
                    group.setup(rerun=rerun)
                    pbar.update(1)
                pbar.close()
        else:
            return False

    def status(self, print_msg=True):
        """Returns a status message for statistics of sub-configuration
        execution.
        Args:
            print_msg (bool): when True, return a text message with aggregate status
              information; otherwise, return a dict of the numbers involved
        """
        from numpy import count_nonzero as cnz
        self._expand_sequence()
        ready = {}
        done = {}

        summaries = {}
        if len(self.sequence) > 0:
            for group in tqdm(self.sequence.values()):
                subsummary = group.status(False)
                summaries[group.key] = subsummary

        allatoms = list(self.config_atoms.values())
        allconfigs = list(self.configs.values())
        for f, a in zip(allconfigs,allatoms):
            ready[f] = a.calc.can_execute(f)
            done[f] = a.calc.can_extract(f)

        N = len(self.configs)
        is_busy = self.is_executing()

        for groupname, summary in summaries.items():
            ready.update({"{}.{}".format(groupname, k): v
                          for k, v in summary["ready"].items()})
            done.update({"{}.{}".format(groupname, k): v
                          for k, v in summary["done"].items()})
            N += summary["stats"]["N"]
            is_busy = is_busy and summary["busy"]

        rdata, ddata = cnz(list(ready.values())), cnz(list(done.values()))
        rmsg = "ready to execute {}/{};".format(rdata, N)
        dmsg = "finished executing {}/{};".format(ddata, N)
        busy = " busy executing..." if is_busy else ""

        if print_msg:
            return "{} {}{}".format(rmsg, dmsg, busy)
        else:
            result = {
                "ready": ready,
                "done": done,
                "stats": {
                    "ready": rdata,
                    "done": ddata,
                    "N": N
                },
                "busy": is_busy
            }
            return result

    def can_extract(self):
        """Runs post-execution routines to clean-up the calculations. This super class
        implementation only checks that each of the sub-config directories has
        the necessary files needed for extraction of output.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            if (len(self.configs) != self.nconfigs and
                self.nconfigs is not None):
                #We need to have at least one folder for each config;
                #otherwise we aren't ready to go.
                return False

            result = False
            for f, a in zip(self.configs.values(), self.config_atoms.values()):
                if not a.calc.can_extract(f):
                    msg.std("Config {} not ready for extraction.".format(f), 2)
                    break
            else:
                result = True
            return result
        else:
            return all(group.can_extract() for group in self.sequence.values())

    def tarball(self, filename="output.tar.gz"):
        """Creates a zipped tar archive that contains each of the specified
        files in sub-sampled configurations' output folders.

        Args:
            filename (str): name of the zipped archive to create.
        """
        if len(self.sequence) == 0:
            parts = []
            for fname in self.calc.tarball:
                parts.append("{}.*/{}".format(self.prefix, fname))

            targs = ["tar", "-cvzf", filename, ' '.join(parts)]
            # from matdb.utility import execute
            execute(targs, self.root)
        else:
            for group in self.sequence.values():
                group.tarball(filename)

    def extract(self, cleanup="default", asis=False):
        """Creates a hdf5 file for each atoms object in the group.
        Args:
            cleanup (str): the level of cleanup to perform after
              extraction.
            asis (bool): when True, read the results even if the calculation
              didn't converge to required accuracy.
        """        
        self._expand_sequence()
        if len(self.sequence) == 0 and self.can_extract():
            for cid, folder in self.configs.items():
                if path.isfile(path.join(folder, "atoms.h5")):
                    #We don't need to recreate the atoms objects if they already
                    #exist.
                    continue
                from os import remove
                atoms = self.config_atoms[cid]
                #We only write the atoms.h5 if extraction returns successful.
                if atoms.calc.extract(folder, cleanup=cleanup, asis=asis):
                    atoms.write(path.join(folder, "atoms.h5"))
                    
                # For debugging, we really don't want to remove these yet;
                # otherwise it is a *big pain* to recreate them.
                # if path.isfile(path.join(folder, "pre_comp_atoms.h5")):
                #     remove(path.join(folder, "pre_comp_atoms.h5"))
            return self.can_extract()
        elif len(self.sequence) >0:
            pbar = tqdm(total=len(self.sequence))
            cleaned = []
            for group in self.sequence.values():
                cleaned.append(group.extract(cleanup=cleanup, asis=asis))
                pbar.update(1)
            cleaned = all(cleaned)
            pbar.close()
            if cleaned:
                atoms = self.rset
                atoms_dict = {"atom_{}".format(f.uuid): f.to_dict() for f in atoms}
                from matdb.io import save_dict_to_h5
                with h5py.File(path.join(self.root,"rset.h5"),"w") as hf:
                    save_dict_to_h5(hf,atoms_dict,'/')
            return cleaned
        else:
            return self.can_extract()

    def finalize(self):
        """Returns a dictionary of the uuid, hash, parameters, and the rset
        for the group.
        """

        from matdb.atoms import _recursively_convert_units

        final_dict = self.to_dict()
        final_dict["hash"] = self.hash_group()
        final_dict["uuid"] = self.uuid

        rset = self.rset
        atoms = AtomsList(rset)

        # The rset exists for a set of parameters. We want to know
        # which sets of parameters went into making the rset.
        msg.info("Finalizing {}".format(self.name))
        pbar = tqdm(total=len(atoms))
        params_dict = {}
        uuids = []
        for atm in atoms:
            this_uuid = atm.group_uuid
            if this_uuid not in uuids:
                # We need to get the parameters for this instance of
                # the group and it's key to store the parameters under
                uuids.append(this_uuid)
                group_inst = self.database.parent.find(this_uuid)
                key = group_inst.key
                params = group_inst.to_dict(include_time_stamp=True)
                params[key] = group_inst.grpargs
                params_dict[key] = _recursively_convert_units(params)
            pbar.update(1)

        return final_dict

def _conform_atoms(atoms, ekey, fkey, vkey, hesskey):
    """Renames the parameters and properties in the specified atoms object to
    conform to the energy, force, virial and hessian fitting names.

    Args:
        atoms (matdb.Atoms): configuration to conform values to.
        ekey (str): existing energy parameter key name.
        fkey (str): existing force parameter key name.
        vkey (str): existing virial parameter key name.
        hesskey (str): existing hessian parameter key name.

    Returns:

    matdb.Atoms: new atoms object with the quantities renamed.
    """
    #We need to rename the parameters and properties of the individual atoms
    #objects to match the refkey and global choice of "ref_energy",
    #"ref_force" and "ref_virial". *NB* for some of the fitting
    #configs (for example Hessian fitting), a config may only have
    #an eigenvalue/eigenvector pair and no energy, force or virial
    #information.
    ati = Atoms(path.join(atoms, "atoms.h5"))
    if ekey in ati.params:
        energy = ati.params[ekey]
        ati.params["ref_energy"] = energy
        del ati.params[ekey]
    if fkey in ati.properties:
        force = ati.properties[fkey]
        ati.properties["ref_force"] = force
        del ati.properties[fkey]
    if vkey in ati.params:
        virial = ati.params[vkey]
        ati.params["ref_virial"] = virial
        del ati.params[vkey]

    #There may be many hessian parameters depending on whether
    #custom parameters are used per-eigenvalue.
    for pname in list(ati.params.keys()):
        if hesskey in pname:
            eigval = ati.params[pname]
            repkey = pname.replace(hesskey, "ref_hessian")
            ati.params[repkey] = eigval
            del ati.params[pname]
    for pname in list(ati.properties.keys()):
        if hesskey in pname:
            eigvec = ati.properties[pname]
            repkey = pname.replace(hesskey, "ref_hessian")
            ati.properties[repkey] = eigvec
            del ati.properties[pname]

    return ati

class Database(object):
    """Represents a Database of groups (all inheriting from :class:`Group`) that
    are all related be the atomic configuration that they model.
    .. note:: See the list of attributes below which correspond to the sections
      in the YAML database specification file.
    Args:
        name (str): name of the configuration that this database sequence is
          operating for.
        root (str): root directory in which all other database sequences for
          the configurations in the same specification will be stored.
        parent (Controller): instance controlling multiple configurations.
        steps (list): of `dict` describing the kinds of sub-configuration
          database steps to setup.
        splits (dict): keys are split names; values are `float` *training*
          percentages to use.
        ran_seed (int): random seed to use for splits in this database.

    Attributes:
        title (str): title for the alloy system that these databases work with.
        config (str): name of the configuration that this database sequence is
          running for.
        species (list): of `str` element names in the alloy system.
        potcars (dict): keys are lower-case element names; values are *suffixes*
          for the various pseudo-potentials that ship with VASP.
        root (str): root directory in which all other database directories will
          be stored. Defaults to the current directory.
        plotdir (str): path to the directory to store plots in for this
          database. Defaults to the parent controller's plot directory.
        incar (dict): keys are valid settings in the INCAR file that should be
          used as defaults for *all* database calculations (with their
          corresponding values).
        kpoints (dict): keys are valid settings for the Mueller k-point PRECALC
          file that should be used as default for *all* database calculations.
        steps (OrderedDict): keys are step names (e.g. `dft`, `calibration`,
          etc.); values are the corresponding class instances.
        parent (Controller): instance controlling multiple configurations.
        rec_bin (RecycleBin): instance of the recycling bin database.
    """
    def __init__(self, name, root, parent, steps, splits, ran_seed):
        self.name = name
        self.config = name.split('.')[0]
        self.root = root
        self.splits = {} if splits is None else splits
        self.ran_seed = ran_seed
        self.splitroot = path.join(root, "splits", name)

        from os import mkdir, makedirs
        if not path.isdir(self.root):
            mkdir(self.root)
        if not path.isdir(self.splitroot):
            makedirs(self.splitroot)

        self.rec_bin = RecycleBin(parent,root,splits)

        parrefs = ["species", "execution", "plotdir", "calculator"]
        for ref in parrefs:
            setattr(self, ref, getattr(parent, ref))
        self.parent = parent

        self._settings = steps
        """dict: with keys and values describing the kinds of step databases to setup.
        """

        # from os import mkdir
        from matdb.utility import ParameterGrid
        self.steps = OrderedDict()
        for dbspec in steps:
            if isinstance(dbspec, six.string_types):
                #This is a reference to an existing database instance that was
                #defined previously.
                instance = self.parent[dbspec]
                self.steps[instance.name] = instance
                continue

            modname, clsname = dbspec["type"].split('.')
            fqdn = "matdb.database.{}".format(modname)
            module = import_module(fqdn)
            if not hasattr(module, clsname):# pragma: no cover
                #We haven't implemented this database type yet, just skip the
                #initialization for now.
                msg.warn("The {0} group has not been implemented yet.".format(clsname))
                continue

            cls = getattr(module, clsname)

            #Make a copy of the original dictionary so that we don't mess up the
            #pointers; then add in the keyword arguments that are missing.
            cpspec = dbspec.copy()
            del cpspec["type"]

            cpspec["pgrid"] = ParameterGrid(cpspec.copy())
            for k in list(cpspec.keys()):
                if "suffix" in k:
                    del cpspec[k]
                elif "*" == k[-1]:
                    cpspec[k[:-1]] = None
                    del cpspec[k]

            cpspec["root"] = self.root
            cpspec["parent"] = self
            cpspec["rec_bin"] = self.rec_bin
            #Handle the special cases where settings are specified uniquely for
            #each of the configurations separately.
            for k in list(cpspec.keys()):
                if isinstance(cpspec[k], dict) and self.config in cpspec[k]:
                    cpspec[k] = cpspec[k][self.config]

            instance = cls(**cpspec)
            self.steps[instance.name] = instance

        if path.isfile(path.join(self.root,"{}_uuid.txt".format(self.name))):
            with open(path.join(self.root,"{}_uuid.txt".format(self.name)),"r") as f:
                uid = f.readline().strip()
                time_stamp = f.readline().strip()
        else:
            uid = str(uuid4())
            time_stamp = str(datetime.now())
            with open(path.join(self.root,"{}_uuid.txt".format(self.name)),"w+") as f:
                f.write("{0} \n {1}".format(uid,time_stamp))

        self.uuid = uid
        self.time_stamp = time_stamp
        self.parent.uuids[uid] = self

    def hash_db(self):
        """Creates the hash for the database.
        """
        hashes = ''
        for name, step in self.steps.items():
            hashes += ' '
            hashes += step.hash_group()

        return str(sha1(hashes.encode()).hexdigest())

    @property
    def iconfigs(self):
        """Returns a generator over all the configurations in all sub-steps of
        this database.
        """
        for dbname, db in self.isteps():
            for config in db.isteps:
                yield config

    @property
    def isteps(self):
        """Returns a generator over steps in this sequence. The generator yields
        the next step *only* if the previous one is already finished (i.e., the
        `ready()` method returns True.
        """
        previous = None
        for dbname, db in self.steps.items():
            if previous is None or previous[1].ready():
                previous = (dbname, db)
                yield previous
            else:
                raise StopIteration()

    def recover(self, rerun=0):
        """Runs recovery on this database to determine which configs failed and
        then create a jobfile to requeue them for compute.
        Args:
            rerun (int): when > 0, recreate the jobfile even if it
              already exists.
        """
        for dbname, db in self.steps.items():
            db.recover(rerun)

    def status(self, busy=False):
        """Prints a status message for each of the databases relative
        to VASP execution status.
        Args:
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
        """
        from matdb.msg import verbosity
        for dbname, db in self.isteps:
            if not busy:
                busy2 = verbosity < 2 if verbosity is not None else True
                imsg = "{}:{} => {}".format(self.name, dbname, db.status(busy2))
                msg.info(imsg)
            else:
                detail = db.status(False)
                running = [k for k, v in detail["done"].items() if not v]
                for config in running:
                    msg.std(config.replace(self.root, ""))

        msg.blank(level=1)

    def execute(self, recovery=False, env_vars=None, dryrun=False):
        """Submits job array files for any of the databases that are ready to
        execute, but which haven't been submitted yet.
        Args:
            recovery (bool): when True, submit the script for running recovery
              jobs.
            env_vars (dict): of environment variables to set before calling the
              execution. The variables will be set back after execution.
            dryrun (bool): when True, don't submit any jobs for execution, just
              say what would happen.
        """
        ready = True
        for dbname, db in self.isteps:
            if ready:
                ready = (db.ready() or db.execute(recovery=recovery,
                                                  env_vars=env_vars,
                                                  dryrun=dryrun))
            if not ready:
                imsg = ("Group {}.{} is not ready to execute yet, or is "
                        "already executing. Done.")
                msg.info(imsg.format(self.name, dbname))
                break
        msg.blank()

    def train_file(self, split):
        """Returns the full path to the h5 database file that can be
        used for training.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.splitroot, "{}-train.h5".format(split))

    def holdout_file(self, split):
        """Returns the full path to the h5 database file that can be
        used to validate the potential fit.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.splitroot, "{}-holdout.h5".format(split))

    def super_file(self, split):
        """Returns the full path to the h5 database file that can be
        used to *super* validate the potential fit.

        Args:
            split (str): name of the split to use.
        """
        return path.join(self.splitroot, "{}-super.h5".format(split))

    def split(self, recalc=0):
        """Splits the database multiple times, one for each `split` setting in
        the database specification.
        """
        subconfs = AtomsList()
        nonsplit = AtomsList()
        for dbname, db in self.isteps:
            if not db.trainable:
                continue

            ekey, fkey, vkey = ["{}_{}".format(db.calculator.key, q)
                                for q in ["energy", "force", "virial"]]
            hesskey = "{}_hessian".format(db.calculator.key)
            for atconf in db.fitting_configs:
                ati = _conform_atoms(atconf, ekey, fkey, vkey, hesskey)
                if db.splittable:
                    subconfs.append(ati)
                else:
                    nonsplit.append(ati)

        file_targets = {"train": self.train_file, "holdout": self.holdout_file,
                        "super": self.super_file}
        split(subconfs, self.splits, file_targets, self.splitroot,
                        self.ran_seed, recalc=recalc, nonsplit=nonsplit)
        
    def extract(self, cleanup="default", asis=False):
        """Runs the extract methods of each database in the collection, in the
        correct order.

        Args:
            cleanup (str): the level of cleanup to perform after extraction.
            asis (bool): when True, read the results even if the calculation
              didn't converge to required accuracy.
        """
        for dbname, db in self.isteps:
            if not db.extract(cleanup=cleanup, asis=asis):
                imsg = "Group {}:{} is not ready yet. Done."
                msg.info(imsg.format(self.name, dbname), 2)
                break
        msg.blank()

    def setup(self, rerun=0):
        """Sets up the database collection by generating the POTCAR file and
        initializing any databases that haven't already been initialized.
        .. note:: The db setup functions are setup to only execute once, and then
           only if their dependencies have completed their calculations. This
           method can, therefore, be safely called repeatedly between different
           terminal sessions.

        Args:
            rerun (int): when > 0, recreate the folders even if they
              already exist. Higher levels redo more of the work.
        """
        for dbname, db in self.isteps:
            msg.info("Setting up database {}:{}".format(self.name, dbname))
            db.setup(rerun)
        msg.blank()

    def to_dict(self):
        """Returns a dictionary of the database parameters and settings.
        """
        # from matdb.utility import __version__
        # from os import sys
        data = {"version":__version__,"python_version":sys.version,"name":self.name,
                "root":self.root,"steps":self._settings,"uuid":self.uuid}
        return data

    def finalize(self):
        """Finalizes the database to a dictionary that can be saved to an h5 file.
        """
        from matdb.atoms import _recursively_convert_units

        final_dict = self.to_dict()
        final_dict["uuid"] = self.uuid

        if not isinstance(self,RecycleBin):
            final_dict["hash"] = self.hash_db()

            for dbname, db in self.isteps:
                final_dict["dbname"] = db.finalize()

            for name, train_perc in self.splits.items():
                for f in glob(path.join(self.root,"splits",self.name,"{0}*-ids.pkl".format(name))):
                    id_file = open(f,'rb')
                    final_dict[f.split(".pkl")[0]] = _recursively_convert_units(pickle.load(id_file), split=True)
        else:
            final_dict["hash"] = self.hash_bin()

        # for name, trani_perc in self.splits.items():
        #     for f in glob(path.join(self.root,"splits",self.name,"{0}*-ids.pkl".format(name))):
        #         id_file = open(f,'r')
        #         final_dict[f.split(".pkl")[0]] = _recursively_convert_units(pickle.load(id_file))

        return final_dict

# from Wiley's comments, no need to cover this class as of M1(2019-06)
class RecycleBin(Database): #pragma: no cover
    """A database of past calculations to be stored for later use.
    Args:
        parent (Controller): instance controlling multiple configurations.
        root (str): root directory in which the 'RecycleBin' folder will
          be placed.
        splits (dict): keys are split names; values are `float` *training*
          percentages to use.
    """

    def __init__(self,parent,root,splits):
        """Sets up the RecycleBin database.
        """
        self.parent = parent
        self.root = path.join(root,"RecycleBin")
        self.splits = {} if splits is None else splits
        if not path.isdir(self.root):
            from os import mkdir
            mkdir(self.root)

        self.steps = {"rec_bin":self}
        self.trainable = True
        self.uuid = str(uuid4())
        self.time_stamp = str(datetime.now())
        self.parent.uuids[self.uuid] = self

    def to_dict(self):
        return {}

    @property
    def rset(self):
        """Returns a list of all the atoms object files in the RecycleBin."""
        from glob import glob
        return glob(path.join(self.root,"*.h5"))

    def hash_bin(self):
        """Returns a hash of the atoms objcets in the recycle bin.
        """
        from matdb.utility import convert_dict_to_str

        hash_str = ""
        for atom in self.rset:
            temp_atom = Atoms(atom)
            hash_str += convert_dict_to_str(temp_atom.to_dict())

        return str(sha1(hash_str.encode()).hexdigest())

    def setup(self):
        pass

    def extract(self):
        pass

    @property
    def isteps(self):
        pass

    def recover(self, rerun=0):
        pass

    def split(self, recalc=0, dfilter=None):
        """Splits the total available data in the recycle bin.
        Args:
            recalc (int): when non-zero, re-split the data and overwrite any
              existing *.h5 files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level h5 files, increase this value.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        from matdb.utility import chdir

        with chdir(self.root):
            super(RecycleBin,self).split(recalc=recalc, dfilter=dfilter)

    def status(self, busy=False):
        pass
        """Prints a status message for each of the databases relative
        to VASP execution status.
        Args:
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
        """
        msg.std("Ready")

    def execute(self, recovery=False, env_vars=None):
        pass

class Controller(object):
    """Implements methods for tying a configuration dictionary (in
    YAML format) to instances of various databases.
    Args:
        config (str): name of the YML file (without the .yml) that
          specifies all information for constructing the set of databases.
        tmpdir (str): path to a temporary directory to use for the
          database. This is for unit testing purposes.

    Attributes:
        specs (dict): the parsed settings from the YAML configuration file.
        collections (dict): keys are configuration names listed in attribute
          `configs` of the YAML file. Values are the :class:`DatabaseCollection`
          instances.
        legacy (dict): keys are legacy database names; values are
          :class:`~matdb.database.legacy.LegacyDatabase`.
        plotdir (str): path to the directory to store plots in for all
          databases.
        venv (str): name of a virtual environment to activate for plotting
          potentials after fitting.
    """
    def __init__(self, config, tmpdir=None):
        from matdb.utility import relpath, check_deps
        from matdb.calculators.utility import set_paths

        self.versions = check_deps()
        self.config = path.expanduser(path.abspath(config))
        if path.isabs(config):
            root, config = path.split(config)
        else:
            root, config = path.dirname(self.config), config
        self.specs = read(root, config)

        #We allow the user to specify paths relative the matdb repo.
        self.root = relpath(path.expanduser(self.specs["root"]))
        if tmpdir is not None:
            self.root = tmpdir
        name = self.specs["title"].strip().replace(' ', '_')
        with open(path.join(self.root, "NAME"), "w+") as f:
            f.write(name)        
        _set_config_paths(name, self.root)
        set_paths(self.specs)

        self.shell_command = "sbatch"
        self.plotdir = path.join(self.root, "plots")
        self.kpathdir = path.join(self.root, "kpaths")
        self.title = self.specs["title"]
        self.legacy = {}
        self.collections = {}
        self.uuids = {}
        self.species = sorted([s for s in self.specs["species"]])
        if "shell_command" in self.specs:
            self.shell_command = self.specs["shell_command"]

        import random
        self.ran_seed = self.specs.get("random_seed", 0)
        random.seed(self.ran_seed)

        self.execution = self.specs.get("execution", {})
        self.calculator = self.specs.get("calculator", {})

        self.venv = self.specs.get("venv")

        # We need to split out the databases by user-given name to create
        # the sequences.
        from matdb.database.legacy import LegacyDatabase
        for dbspec in self.specs.get("databases", []):
            if dbspec.get("legacy", False):
                cpspec = dbspec.copy()
                cpspec["root"] = self.root
                cpspec["controller"] = self
                cpspec["splits"] = self.specs.get("splits")
                #We allow the user to specify the folder relative to repository
                #root; this is mainly for unit tests.
                cpspec["folder"] = relpath(cpspec["folder"])
                del cpspec["legacy"]
                self.legacy[cpspec["name"]] = LegacyDatabase(**cpspec)
            else:
                dbname = dbspec["name"]
                steps = dbspec["steps"]
                #We use the global random seed by default if one isn't specified
                #explicitly for the database.
                db = Database(dbname, self.root, self,
                              steps, self.specs.get("splits"),
                              self.specs.get("random_seed", self.ran_seed))
                self.collections[dbname] = db

        from os import mkdir
        if not path.isdir(self.plotdir):
            mkdir(self.plotdir)
        if not path.isdir(self.kpathdir):
            mkdir(self.kpathdir)

        #If the controller is going to train any potentials, we also need to
        self.trainers = None
        if "fitting" in self.specs:
            from matdb.fitting.controller import TController
            tdict = self.specs["fitting"].copy()
            tdict["root"] = self.root
            tdict["db"] = self
            self.trainers = TController(**tdict)

    def ifiltered(self, dfilter=None):
        """Returns a *filtered* generator over databases in the collections of this object.

        Args:
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.
        """
        from fnmatch import fnmatch
        for name, dbi in self.collections.items():
            if dfilter is None or any(fnmatch(name, d) for d in dfilter):
                yield (name, dbi)

    def relpaths(self, pattern):
        """Finds the relative paths for the seed configurations within the databases that
        match to the pattern.
        Args:
            pattern (str): the pattern to match.
        """
        from matdb.database.utility import parse_path
        return parse_path(self.root,pattern,ran_seed=self.ran_seed)

    def find(self, pattern):
        """Finds a list of :class:`matdb.database.basic.Group` that match the given
        pattern. The pattern is formed using `group.dbname[[.seed].params]`. `*`
        can be used as a wildcard for any portion of the '.'-separated path.
        .. note:: Actually, an :func:`~fnmatch.fnmatch` pattern can be used.
        Args: pattern (str): fnmatch pattern that follows the convention of the
        DB key. Alternatively the pattern can be a uuid. Examples:

            Get all the dynamical matrix databases for the `Pd`
            configuration. The example assumes that the database name is
            `phonon` and that it includes a dynamical matrix step.
            >>> Pd = Controller("Pd.yml")
            >>> Pd.find("DynMatrix.phonon.Pd.*")
            Get all the database sequences for liquids across all configurations in
            the database.
            >>> CdWO4 = Controller("CdWO4.yml")
            >>> CdWO4.find("*.liquid*")
            >>> CdWO4.find(uuid)
        """
        from matdb.utility import is_uuid4
        if is_uuid4(pattern):
            return self.uuids[pattern]

        from fnmatch import fnmatch
        if pattern.count('/') == 3:
            groupname, dbname, seed, params = pattern.split('/')
        elif pattern.count('/') == 2:
            groupname, dbname, seed = pattern.split('/')
            params = None
        elif pattern.count('/') == 1:
            groupname, dbname = pattern.split('/')
            params = None
            seed = None
        else:
            #We must be searching legacy databases or the pattern given is to
            #match a *database* and not a group.
            dbname = pattern
            groupname, params, seed = None, None, None

        colls = [v for k, v in self.collections.items() if fnmatch(k, dbname)]
        colls.extend([li for ln, li in self.legacy.items() if fnmatch(ln, pattern)])
        result = []
        for dbi in colls:
            if isinstance(dbi, LegacyDatabase) or groupname is None:
                result.append(dbi)
                continue
            groups = [groupi for groupn, groupi in dbi.steps.items()
                      if fnmatch(groupn, groupname)]
            for group in groups:
                group._expand_sequence()
                if len(group.sequence) > 0 and seed is not None:
                    seeds = [si for sn, si in group.sequence.items()
                             if fnmatch(sn, seed)]
                    for seedi in seeds:

                        seedi._expand_sequence()
                        if len(seedi.sequence) > 0 and params is not None:
                            result.extend([si for sn, si in seedi.sequence.items()
                                           if fnmatch(sn, params)])
                        else:
                            result.append(seedi)
                else:
                    result.append(group)

        if groupname == '*':
            #Add all the possible legacy databases.
            result.extend([li for ln, li in self.legacy.items()
                           if fnmatch(ln, groupname)])

        return result

    def steps(self):
        """Compiles a list of all steps in this set of databases.
        """
        result = []
        for db_name, db in self.collections.items():
            for group_name, group in db.steps.items():
                if len(group.sequence) > 0:
                    for seed_name, seed in group.sequence.items():
                        if len(seed.sequence) > 0:
                            for param_name, param in seed.sequence.items():
                                result.append("{0}/{1}/{2}/{3}".format(group_name,db_name,
                                                                           seed_name,param_name))
                        else:
                            result.append("{0}/{1}/{2}".format(group_name,db_name,
                                                                   seed_name))
                else:
                    result.append("{0}/{1}".format(group_name,db_name))

        return sorted(result)

    def sequences(self):
        """Compiles a list of all sequences in this set of databases.
        """
        result = []
        for db_name, db in self.collections.items():
            for group_name, group in db.steps.items():
                if len(group.sequence) > 0:
                    for seed_name, seed in group.sequence.items():
                        if len(seed.sequence) > 0:
                            for param_name, param in seed.sequence.items():
                                result.append("{0}/{1}".format(seed_name,param_name))
                        else:
                            result.append("{0}".format(seed_name))

        return sorted(result)

    def __getitem__(self, key):
        """Returns the database object associated with the given key. This is
        necessary because of the hierarchy of objects needed to implement
        sequence repitition via the `ParamaterGrid`.
        """
        if key.count('/') == 3:
            group, dbname, seed, params = key.split('/')
        elif key.count('/') == 2:
            group, dbname, seed = key.split('/')
            params = None
        else:
            group, dbname = key.split('/')
            seed = None
            params = None

        coll = self.collections[dbname]
        if group.lower() in coll.steps:
            step = coll.steps[group.lower()]
            if seed is not None and seed in step.sequence:
                seq = step.sequence[seed]
                if params is not None and params in seq.sequence:
                    return seq.sequence[params]
                else:
                    return seq
            else:
                return step
        else: #pragma: no cover
            msg.err("The group name {0} could not be found in the steps of "
                    "the database {1}".format(group.lower(),coll.steps.values()))

    def setup(self, rerun=0, dfilter=None):
        """Sets up each of configuration's databases.

        Args:
            rerun (int): when > 0, recreate the folders even if they
              already exist.
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases sequences are returned.
        """
        try:
            for dbname, dbi in self.ifiltered(dfilter):
                dbi.setup(rerun)

        finally:
            # Save the paths to the setup atoms objects for the Tracy
            # extract method.
            uuid_paths = {}
            for matdb_id, matdb_obj in self.uuids.items():
                if isinstance(matdb_obj, Atoms):
                    uuid_paths[matdb_id] = matdb_obj.calc.folder
            with open(path.join(self.root,"atoms_paths.json"),"w+") as f:
                json.dump(uuid_paths, f)
            if path.isfile(path.join(self.root, "user_cred.json")): #pragma: no cover
                remove(path.join(self.root, "user_cred.json"))

    def extract(self, dfilter=None, cleanup="default", asis=False):
        """Runs extract on each of the configuration's databases.

        Args:
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.
            cleanup (str): the level of cleanup to perform after extraction.
            asis (bool): when True, read the results even if the calculation
              didn't converge to required accuracy.
        """
        for dbname, dbi in self.ifiltered(dfilter):
            msg.std("Extracting calculation output from {}.".format(dbname), 2)
            dbi.extract(cleanup=cleanup, asis=asis)

    def execute(self, recovery=False, dfilter=None, env_vars=None, dryrun=False):
        """Submits job array scripts for each database collection.

        Args:
            recovery (bool): when True, submit the script for running recovery
              jobs.
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.
            env_vars (dict): of environment variables to set before calling the
              execution. The variables will be set back after execution.
            dryrun (bool): when True, don't submit any jobs for execution, just
              say what would happen.
        """
        for dbname, dbi in self.ifiltered(dfilter):
            dbi.execute(recovery, env_vars=env_vars, dryrun=dryrun)

    def recover(self, rerun=False, dfilter=None):
        """Runs recovery on this database to determine which configs failed and
        then create a jobfile to requeue them for compute.

        Args:
            rerun (bool): when True, recreate the jobfile even if it
              already exists.
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.
        """
        for dbname, dbi in self.ifiltered(dfilter):
            dbi.recover(rerun)

    def status(self, busy=False, dfilter=None):
        """Prints status messages for each of the configuration's databases.

        Args:
            busy (bool): when True, print a list of the configurations that are
              still busy being computed in DFT.
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.

        """
        for dbname, dbi in self.ifiltered(dfilter):
            dbi.status(busy)

    def split(self, recalc=0, dfilter=None):
        """Splits the total available data in all databases into a training and holdout
        set.
        Args:
            recalc (int): when non-zero, re-split the data and overwrite any
              existing *.h5 files. This parameter decreases as
              rewrites proceed down the stack. To re-calculate
              lower-level h5 files, increase this value.
            dfilter (list): of `str` patterns to match against *database sequence*
              names. This limits which databases sequences are returned.
        """
        for dbname, dbi in self.ifiltered(dfilter):
            dbi.split(recalc)

    def hash_dbs(self, dfilter=None):
        """Hashes the databases into a single hash.

        Args:
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.
        """
        hash_all = ''
        for dbname, seq in self.ifiltered(dfilter):
            hash_all += ' '
            hash_all += seq.hash_db()

        hash_all += ' '
        hash_all += seq.rec_bin.hash_bin()

        hash_all = str(sha1(hash_all.encode()).hexdigest())
        with open(path.join(self.root,"hash.txt"),"w+") as f:
            f.write("{0} \n {1}".format(hash_all,str(datetime.now())))

        return hash_all

    def verify_hash(self, hash_cand, dfilter=None):
        """Verifies that the the candidate hash matches this matdb controller's databases.

        Args:
            hash_cand (str): The candidate hash to check.
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.

        Returns:
            True if the hash matches the databases.
        """
        return hash_cand == self.hash_dbs(dfilter=dfilter)

    def finalize(self, dfilter=None):
        """Creates the finalized version of the databases that were used for
        fitting the potential. Stored in final_{matdb.version}.h5.
        Args:
            dfilter (list): of `str` patterns to match against *database*
              names. This limits which databases are returned.
        """

        # from matdb.io import save_dict_to_h5
        # from matdb import __version__

        final_dict = self.versions.copy()
        final_dict["hash"] = self.hash_dbs(dfilter=dfilter)
        final_dict["yml_file"] = self.specs.copy()
        for dbname, seq in self.ifiltered(dfilter):
            final_dict[dbname] = seq.finalize()

        if (dfilter is not None and "rec_bin" in dfilter) or dfilter is None or dfilter=="*":
            final_dict["rec_bin"] = seq.rec_bin.finalize()

        str_ver = []
        for item in __version__:
            str_ver.append(str(item))
        mdb_ver = ".".join(str_ver)
        target = path.join(self.root,"final_{}.h5".format(mdb_ver))

        with h5py.File(target,"w") as hf:
            save_dict_to_h5(hf,final_dict,'/')
