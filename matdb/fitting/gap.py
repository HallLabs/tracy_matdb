"""Implements classes and methods for performing a GAP fit over a
database defined in the YAML specification file.
"""
from os import path
from matdb import msg
from collections import OrderedDict
import numpy as np
from tqdm import tqdm

from .basic import Trainer
import lazy_import
calculators = lazy_import.lazy_module("matdb.calculators")

from matdb.calculators import Quip
from matdb.atoms import AtomsList
from matdb.database.hessian import Hessian

def update_nbody(settings):
    """Adds the usual n-body settings to the specified dictionary. This function
    *alters* the `settings` parameter.

    Args:
        settings (dict): key-value pairs that are accepted by the GAP
          `distance_Nb` potential type.
    """
    usuals = {
        "compact_clusters": True,
        "sparse_method": "uniform",
        "covariance_type": "ARD_SE",
        "theta_uniform": 1.0,
        "cutoff_transition_width": 1.0
    }
    for k, v in usuals.items():
        if k not in settings:
            settings[k] = v

def update_soap(settings):
    """Adds the usual soap settings to the specified dictionary. This function
    *alters* the `settings` parameter.

    Args:
        settings (dict): key-value pairs that are accepted by the GAP
          `soap` potential type.
    """
    usuals = {
        "cutoff_transition_width": 1.0,
        "atom_sigma": 0.5,
        "zeta": 2,
        "sparse_method": 'cur_points',
        "covariance_type": 'dot_product'
    }
    for k, v in usuals.items():
        if k not in settings:
            settings[k] = v

def dict_to_str(settings, spacer=""):
    """Converts the specified dictionary of QUIP-compatible settings into a
    single string.

    Args:
        settings (dict): key-value pairs that are recognized settings by QUIP
          for descriptors, teach_sparse, eval, etc.
    """
    result = []
    for k, v in settings.items():
        if isinstance(v, (list, set, tuple)):
            value = '{' + ' '.join(map(str, v)) + '}'
        else:
            value = v

        result.append("{}={}".format(k, value))

    joiner = " \\\n  {}".format(spacer)
    return joiner.join(result)

def _n_neighbors(atoms, cutoff):
    """Returns the *average* number of nearest or next-nearest neighbors
    depending on the specified reference cutoff.

    Args:
        atoms (quippy.Atoms): seed configuration to use for getting neighbors.
        cutoff (float): reference cutoff radius around each atom.
    """
    acopy = atoms.copy()
    acopy.set_cutoff(cutoff)
    acopy.calc_connect()
    neighs = []
    for i in range(acopy.n):
        neighs.append(acopy.n_neighbours(i+1))
    
    return np.mean(neighs)

def _rescale_2body(atoms, settings):
    """Calculates the appropriate scaling factor for the RMS error of a 2-body
    potential. 

    Args:
        atoms (quippy.Atoms): seed configuration to use for getting neighbors.
        settings (dict): key-value pairs used to configure the 2-body potential.
    """
    n = _n_neighbors(atoms, settings["cutoff"])
    msg.okay("2-body: scaling by {}.".format(n))
    return n

def _rescale_3body(atoms, settings):
    """Calculates the appropriate scaling factor for the RMS error of a 3-body
    potential. 

    Args:
        atoms (quippy.Atoms): seed configuration to use for getting neighbors.
        settings (dict): key-value pairs used to configure the 3-body potential.
    """
    #First, get the number of nearest neighbors in the cutoff.
    n = _n_neighbors(atoms, settings["cutoff"])
    result = n*(n-1)/2.
    msg.okay("3-body: scaling by {} (from {} neighbors).".format(result, n))
    return result

class GAP(Trainer):
    """Implements a simple wrapper around the GAP training functionality for
    creating GAP potentials.

    Args:
        nb (int): number of bodies to include in the model. Usually one of 2, 3
          or 4 for `distance_Nb` calculations. Use -1 for many-body.
        controller (matdb.fitting.controller.Controller): fitting controller
          provides access to previous fitting steps and training/validation data.
        sigmas (dict): of `float` specifies the expected variance in the energy,
          force, virial, and hessian respectively. It should have length
          4. If a `dict`, then include keys for those config types that should
          be overridden and use "__" for the default sigma values.
        split (str): name of the split specification to use for training.
        execution (dict): settings needed to configure the jobfile for running
          the GAP fit.
        dbs (list): of `str` patterns from the database that should be included
          in the training and validation.
        ogaps (list): of `str` trainer FQNs for previous gaps that were fitted,
          whose parameters should be included in the final GAP IP.
        n_random (int): number of random configs to "rattle" using the
          Hessian group's seed configurations. These are added to the training
          database *without* any energy/force/virial information.
        teach_sparse (dict): key-value pairs for the `teach_sparse` executable.
        gap (dict): parameters for the n-body potential.

    Attributes:
        name (str): name of the folder in which all the fitting for this trainer
          takes place; defaults to `{n}b`.
        root (str): root directory that this trainer operates in.
        e0 (float): reference "zero" for the potential energy surface.
        gap (dict): potential key-value pairs for configuring that particular
          GAP.
        seeds (list): of tuples with `(atoms, sequence)`, where `atoms` is the
          the seed configuration :class:`quippy.Atoms` objects and `seq` is the
          database sequence that was derived from it. The atoms objects may have
          additional with hessian information attached if `hessfit=True`.
        params (dict): key-value pairs that are parameters for the model
          fitting.
    """
    def __init__(self, nb=None, sigmas=None, ogaps=None, teach_sparse=None,
                 gap=None, n_random=0, **trainargs):
        self.name = "{}b".format(nb) if nb > 0 else "mb"
        super(GAP, self).__init__(**trainargs)
        self.nb = nb
        self.e0 = self.controller.e0
        self.sigmas = sigmas
        self.gap = {} if gap is None else gap.copy()
        self.params.update(self.gap)
        self.n_random = n_random
        self.params["n_random"] = n_random

        #We have to run the descriptor string to populate the parameters
        #correctly.
        self.desc_str()
        if sigmas is not None:
            self.params.update({"sigma_{}".format(k): v for k, v in sigmas.items()})
        
        self.ogaps = [] if ogaps is None else [controller[gap] for gap in ogaps]
        self.gp_file = "{}.xml".format(self.name)
        self.teach_sparse = {} if teach_sparse is None else teach_sparse
        
        #Configure the fitting directory for this particular set of
        #potentials. This way, we can get separate directories if the number of
        #n-body or soap potentials changes.
        from os import mkdir
        if not path.isdir(self.root):
            mkdir(self.root)

        #If we are doing hessian training, copy the seed configuration and
        #updates its hessian attributes. Calculating delta uses seed
        #configurations anyway.
        self.seeds = []
        for dbi in self.dbs:
            for group in dbi.steps.values():
                if not group.trainable:
                    continue
                group._expand_sequence()
                self.seeds.append((group.atoms.copy(), group))

    def get_calculator(self):
        """Returns an instance of :class:`ase.Calculator` using the latest
        fitted GAP potential in this trainer.
        """
        if path.isfile(self.gp_file):
            return Quip("IP GAP", param_filename=self.gp_file)

    def ready(self):
        return path.isfile(path.join(self.root, self.gp_file))
        
    @property
    def delta_cache(self):
        """Returns the path the `delta` parameter file for each of the GAPs.
        """
        return path.join(self.root, "delta.dat")
        
    def _get_delta(self):
        """Estimates the value to use for the `delta` parameter for the next GAP
        in the list.

        Returns:
            float: value for the `teach_sparse` delta parameter of the next GAP
            in the sequence.
        """
        #We cache the `delta` value for each GAP in a file so that we don't have
        #to keep re-calculating them.
        if path.isfile(self.delta_cache):
            with open(self.delta_cache) as f:
                delta = float(f.read().strip())
            return delta

        if len(self.ogaps) == 0:
            al = AtomsList(self._trainfile)
            rms = np.std(al.dft_energy)
            msg.info("2B: {0:.3f} data variance.".format(rms), 2)
        else:
            ens = self.ogaps[-1].validate("ref", force=False, virial=False)
            rms = np.std(ens["e_ref"]-ens["e_ip"])
            msg.info("{0}: {1:.3f} RMS error.".format(self.latest, rms), 2)
        
        #We use the first seed configuration in the database and calculate how
        #many nearest neighbors it has within the specified cutoff of the n-body
        #configuration. SOAP doesn't rescale because the descriptors are already
        #per-atom.
        atoms = self.seeds[0][0]
        scalers = {
            2: _rescale_2body,
            3: _rescale_3body,
            -1: lambda a, s: 1.
        }
        scaling = scalers[self.nb](atoms, self.gap)
            
        #Now, we scale the expected error by the number of descriptors per atom
        #in this particular GAP.
        delta = rms/scaling
            
        #Now that we have updated the delta, save it to file.
        with open(self.delta_cache, 'w') as f:
            f.write("{0:.8f}".format(delta))

        return delta
                
    def desc_str(self, spacer="  "):
        """Returns the descriptor string for this training object, if it
        exists.
        """
        pdict = self.gap.copy()
        if "delta" not in pdict:
            pdict["delta"] = self._get_delta()        
        
        #Add default parameters for the parameters.
        if self.nb < 0:
            update_soap(pdict)
        else:
            update_nbody(pdict)
            if "order" not in pdict:
                pdict["order"] = self.nb
        #Overwrite our parameter dictionary for the model. This is needed by
        #other methods and functions that want to do plotting, etc.
        self.params.update(pdict)
                
        settings = []
        settings.append("soap" if self.nb < 0 else "distance_Nb")
        settings.append(dict_to_str(pdict, spacer="  "))
        settings.append("add_species")
        settings.append("")
        settings.append("")

        joiner = " \\\n  {0}".format(spacer)
        return joiner.join(settings)

    def _create_xyz(self):
        """Creates the training.xyz file that `teach_sparse` needs to use.
        """
        import quippy
        target = path.join(self.root, "train.xyz")
        if not path.isfile(target):
            al = AtomsList(self._trainfile)
            ol = quippy.AtomsList()
            for at in tqdm(al):
                #Empty dictionaries in info breaks the quippy fortran
                #implementation. Delete these entries out.
                if len(at.info["params"]) == 0:
                    del at.info["params"]
                if len(at.info["properties"]) == 0:
                    del at.info["properties"]
                
                ai = quippy.Atoms()
                ai.copy_from(at)
                if "properties" in at.info:
                    #We also need to copy the properties in our info onto the
                    #properties of the quippy.Atoms object.
                    ai.arrays.update(at.info["properties"])
                if "params" in ai.params:
                    del ai.params["params"]
                    
                ol.append(ai)
            ol.write(target)

    @property
    def sparse_file(self):
        """Returns the path to the sparse point file that includes random
        configurations to include in the training set, but which have no
        energy/force/virial information.
        """
        return path.join(self.root, "sparse_points.h5")
        
    def _sparse_points(self, recalc=False):
        """Creates random sparse points by rattling atoms for hessian-based training.

        Args:
            recalc (bool): when True, recalculate the sparse points even if the
              file already exists.
        """
        if self.n_random == 0:
            return
        if path.isfile(self.sparse_file) and not recalc:
            return

        #Determine how many random configs we need to generate.
        hessians = [seed for seed, db in self.seeds
                    if isinstance(db, Hessian)]
        n_ratio = int(np.ceil(float(self.n_random)/len(hessians)))        
        N = n_ratio * len(hessians)
        msg.info("Generating {0:d} random configs as sparse points.".format(N))
        
        result = AtomsList()
        for hseed in hessians:
            for i in range(n_ratio):
                atRand = hseed.copy()
                p = atRand.get_positions()
                atRand.set_positions(p  + 0.1*2*(np.random.random_sample(p.shape)))
                result.append(atRand)

        result.write(self.sparse_file)

    def command(self):
        """Returns the `teach_sparse` command that is needed to train the GAP
        potentials specified by this object.

        .. note:: This method also configures the directory that the command
          will run in so that it has the relevant files.
        """
        #Generate any random sparse points and the seed XYZ training file.
        self.compile()
        self._sparse_points()
        self._create_xyz()

        template = ("teach_sparse at_file={train_file} \\\n"
                    "  gap={{ \\\n"
                    "    {gap} \\\n"
                    "  }} \\\n"
                    "  {teach_sparse}")
        
        #Compile the GAP string specifying the potential. Check to see if we
        #have other GAPs that were referenced.
        gaplist = []
        for igap in self.ogaps:
            gaplist.append(igap.desc_str())
        gaplist.append(self.desc_str())
        gapstr = ':'.join(gaplist)

        #e0 needs to be handled correctly for each pure species.
        if isinstance(self.e0, dict):
            e0s = ["{0}:{1:.6f}".format(*i) for i in self.e0.items()]
            e0 = '{' + ':'.join(e0s) + '}'
        else:
            e0 = "{0:.6f}".format(self.e0)

        #sigmas can be speciifed per config_type or just as a list of default
        #ones.
        custom = []
        tsattrs = {
            "e0": e0,
            "gp_file": self.gp_file
        }
        #Override any of the default values using user-specified ones.
        tsattrs.update(self.teach_sparse)
            
        if isinstance(self.sigmas, (list, tuple)):
            tsattrs["default_sigma"] = self.sigmas
        else:
            for k, sigs in self.sigmas.items():
                if k == "__":
                    tsattrs["default_sigma"] = sigs
                else:
                    svals = ':'.join(["{0:.3f}".format(s) for s in sigs])
                    custom.append(k + ':' + svals)

        #Add in the names of the energy, force and virial parameters
        #from dft.
        tsattrs["energy_parameter_name"] = "ref_energy"
        tsattrs["force_parameter_name"] = "ref_force"
        tsattrs["virial_parameter_name"] = "ref_virial"
        tsattrs["hessian_parameter_name"] = "ref_hessian"
        tsattrs["config_type_parameter_name"] = "config_type"
        tsattrs["sigma_parameter_name"] = "csigma"
        
        if len(custom) > 0:
            #If we have custom sigmas, add them in; make sure we have specified
            #which parameter to use for custom_types.
            tsattrs["config_type_sigma"] = '{' + ':'.join(custom) + '}'

        tsstr = dict_to_str(tsattrs)                
        fields = {
            "train_file": "train.xyz",
            "gap": gapstr,
            "teach_sparse": tsstr
        }

        return template.format(**fields)

    def status(self, printed=True):
        """Returns or prints the current status of the GAP training.

        Args:
            printed (bool): when True, print the status to screen; otherwise,
              return a dictionary with the relevant quantities.
        """
        # Our interest is in knowing which GAP model is the latest (if any) and
        # whether the job script has been created for the next one in the
        # sequence.
        result = {
            "trained": path.isfile(self.gp_file),
            "file": self.gp_file,
            "jobfile": path.isfile(self._jobfile)
        }
        
        if printed:
            fqn = "{}.{}".format(self.parent.name, self.name)
            msg.info("{} => Model ready: {}".format(fqn, result["trained"]))
            x = "exists" if result["jobfile"] else "does not exist"
            msg.info("{} => Next jobfile '{}' {}".format(fqn, self._jobfile, x))
            msg.blank()
        else:
            return result
