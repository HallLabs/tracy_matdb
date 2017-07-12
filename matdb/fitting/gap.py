"""Implements classes and methods for performing a GAP fit over a
database defined in the YAML specification file.
"""
from os import path
from matdb import msg
from collections import OrderedDict
import numpy as np
import quippy

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
        "theta_uniform": 1.0
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

def dict_to_str(settings):
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
    return ' '.join(result)

def gpfile_name(gaps, prefix="gp", hessfit=False):
    """Returns the expected gp file name for the specified set of GAPs.

    Args:
        gaps (list): of `int` body sizes for the potentials. SOAP has value 10
          by arbitrary definition.
        hessfit (bool): when True, the fit is being performed using Hessian
          fitting.

    Returns:
        str: gp_{nbody_list}(_{soap})?.xml if valid potentials are in the list;
        otherwise None.
    """
    soap = "_soap" if 10 in gaps else ""
    nbody = 'b'.join([str(k) for k in gaps if k != 10]) + 'b'
    hess = "h" if hessfit else ""
    if len(soap) > 0 or len(nbody) > 1:
        return "{}{}_{}{}.xml".format(prefix, hess, nbody, soap)
    else:
        return None

def _extract_pots(settings):
    """Determines the configuration settings for each of the GAPs or other
    potential pieces specified in the `settings` dictionary.

    Args:
        settings (dict): keys are GAP types (e.g., 2-body, 3-body, SOAP); values
          are `dict` with key-value pairs for the `teach_sparse` command.

    Returns:

        tuple: `(gaps, others)` where gaps is a :class:`collections.OrderedDict`
        of GAP names in the obvious order (2b, 3b, Nb, SOAP) with keys being
        integer ids representing the n-body order (SOAP: 10) and values being
        the string key in `settings`; others represents any other settings' keys
        in the dictionary that are not recognized as parts of potentials.

    """
    from operator import itemgetter
    gaps = []
    others = []
    
    #Before we set the other attributes, let's first see what potentials are
    #present in the settings dictonary.
    from fnmatch import fnmatch
    for k in settings:
        if k.lower() == "soap":
            #The 10 here is arbitrary; we will never create a 10-body
            #potential using distance_Nb.
            gaps.append((10, "soap"))
        elif fnmatch(k.lower(), "*body"):
            n = int(k[0])
            gaps.append((n, k.lower()))
        else:
            others.append(k)

    ogap = OrderedDict(sorted(gaps, key=itemgetter(0)))
    return (ogap, others)

def _get_max_gap(root, hessfit=False):
    """Returns the list of maximally executed GAPs as part of the iterative
    training scheme.

    Args:
        root (str): root directory in which to search for potential files.
        hessfit (bool): when True, look for the maximally executed GAP that was
          trained using hessian fitting.

    Returns:
        list: of `int` n-body codes (SOAP = 10) defining the progress so far in
        the overall fitting scheme.
    """
    from matdb.utility import chdir
    from glob import glob
    nbmax = 0
    with chdir(root):
        pattern = "gph*.xml" if hessfit else "gp*.xml"
        nmax = 0
        soap = False
        
        for potname in glob(pattern):
            potstr = potname[:-4]
            parts = potstr.split('_')
            if len(parts) > 2:
                _, nbstr, soapstr = parts
            else:
                _, nbstr = parts
                soapstr = ""

            #It is possible that we are doing a soap-only fit, in which case the
            #n-body integer will be missing.
            nbmax = nbstr.rstrip('b').split('b')[-1]
            if nbmax != '':
                nbmax = int(nbmax)
            else:
                nbmax = 0
                
            if nbmax > nmax:
                nmax = nbmax
            if soapstr == "soap":
                soap = True

    if nbmax > 0:
        xgaps = list(range(2, nmax+1))
    else:
        xgaps = []
        
    if soap:
        xgaps.append(10)
    return xgaps

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
    for i in acopy.indices:
        for neighb in acopy.connect[i]:
            dist = neighb.distance
            deltain = [abs(dist-s) < 1e-5 for s in neighs]
            if not any(deltain):
                neighs.append(dist)
    shells = sorted(neighs)

    #Determine whether we should use nearest or next-nearest neighbors.
    if len(shells) > 4:
        ncut = shells[1]
    else:
        ncut = shells[0]
    ncut += 1e-4
    
    acopy.set_cutoff(ncut)
    acopy.calc_connect()    
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

class GAPTrainer(object):
    """Implements a simple wrapper around the GAP training functionality that
    can integrate with the methods of the `database` sub-package.

    .. note:: Creating GAP potentials is a multi-fitting procedure. For a
      `2 + 3 + ...` potential, we first need to fit the 2-body, and then (based
      on its errors), configure the fitting for the 2 + 3 potential. The results
      from the previous fit are used to configure the fitting that includes the
      next, higher-order potential type. As a result, the fits take longer as
      more complexity is included *and* the previous, lower-body terms have to
      be refit as part of the total many-body potential.

    Args:
        gap (list): of `str` that are either `*body` for the n-body potentials
          or `soap`. This determines which potentials are included as part of
          the fit.
        db (matdb.database.controller.Controller): database controller whose
          data will be used to train the potentials.
        gp_file (str): name of the output GAP potential file. Defaults to
          `gp_nbody_soap.xml` where `nbody` is a list of the n-body potentials
          included  and `soap` is included if it was part of the fit.
        e0 (float): reference "zero" for the potential energy surface.
        default_sigma (list): of `float` specifies the expected variance in the
          energy, force, virial, and hessian respectively. It should have length
          4. Defaults to `[0.001, 0.01, 1.0, 1.0]`.
        sparse_jitter (float): jitter added to the matrices before inversion to
          help improve the condition number.
        split (float): percentage of the available data to use for
          training (vs. validation).
        root (str): root directory for the entire alloy database
          (i.e., corresponding to a single YAML specification file).
        execution (dict): settings needed to configure the jobfile for running
          the GAP fit.
        configs (list): of `str` configuration ids from the database
          that should be included in the training and validation.
        kwargs (dict): additional parameters accepted by `teach_sparse` or a
          `*body` or `soap` dict that defines the settings for each particular
          GAP.

    Attributes:
        ogap (collections.OrderedDict): keys are integer n-body values
          (SOAP=10), values are the ["2-body", "3-body", ..., "soap"] string
          values.
        gap (list): of `str` that are either `*body` for the n-body potentials
          or `soap`. This determines which potentials are included as part of
          the fit.
        db (matdb.database.controller.Controller): database controller whose
          data will be used to train the potentials.
        gp_file (str): name of the output GAP potential file. Defaults to
          `gp_nbody_soap.xml` where `nbody` is a list of the n-body potentials
          included  and `soap` is included if it was part of the fit.
        name (str): name of the folder in which all the fitting for this trainer
          takes place; defaults to :attr:`gp_file` without the ".xml"
          extension.
        root (str): root directory that this trainer operates in.
        e0 (float): reference "zero" for the potential energy surface.
        default_sigma (list): of `float` specifies the expected variance in the
          energy, force, virial, and hessian respectively. It should have length
          4. Defaults to `[0.001, 0.01, 0.1, 0.1]`.
        sparse_jitter (float): jitter added to the matrices before inversion to
          help improve the condition number.
        split (float): percentage of the available data to use for
          training (vs. validation).
        execution (dict): settings needed to configure the jobfile for running
          the GAP fit.
        teach_sparse (dict): key-value pairs that are passed directly to
          `teach_sparse` that are not in ["e0", "default_sigma",
          "sparse_jitter"].
        gaps (dict): keys are one of ["2body", "3body", ..., "soap"] describing
          the kind of potential. Values are the actual potential key-value
          pairs for configuring that particular GAP.
        seeds (dict): keys are database names in the controller; values are
          copies of the seed configuration :class:`quippy.Atoms` objects in each
          of those databases, with hessian and force information attached.
        state (list): of `int` GAP n-body codes in :attr:`gap` indicating the current
          progress/state in the multi-chain fitting process. This value iterates
          over each of the GAPs in turn, adding one extra one each time the
          previous set has finished executing.
        hessfit (bool): when True, the trainer is performing a hessian-based fit
          over the configurations (as opposed to using the entire configuration
          library).
        latest (str): name of the latest, complete potential file that can be
          used to make predictions for configurations.
        configs (list): of `str` configuration ids from the database
          that should be included in the training and validation.

    """
    def __init__(self, gap=None, db=None, gp_file=None, e0=0., default_sigma=None,
                 sparse_jitter=1e-12, split=0.5, root=None,
                 execution={}, configs=None, **kwargs):
        self.db = db
        if db is None:
            raise ValueError("Cannot train a GAP without a database of "
                             "configurations!")

        ogap, others = _extract_pots(kwargs)
        self.gap = ogap.values() if gap is None else gap
        self.ogap = ogap
        
        #Check that if they specified a list of potentials, that it is ordered
        #correctly and contains all the pieces.
        assert self.gap == ogap.values()
        
        #There are certain required parameters that must have default values. We
        #set these as named keyword arguments in the constructor. All other
        #teach_sparse parameters are contained in kwargs and were separated from
        #the regular potential definitions by _extract_pots() above.
        self.e0 = e0
        if default_sigma is None:
            self.default_sigma = [0.001, 0.01, 0.1, 0.1]
        else:
            self.default_sigma = default_sigma
        self.sparse_jitter = sparse_jitter
        self.split = split
        self.teach_sparse = {k: kwargs[k] for k in others}
        self.execution = execution.copy()
        
        self.configs = ['*'] if configs is None else configs
        if not isinstance(self.configs, (list, set, tuple)):
            self.configs = [self.configs]
        
        #Now, we can finally setup the GAP fitting parameters.
        self.gaps = OrderedDict()
        for kgap in self.gap:
            cdict = kwargs[kgap].copy()
            if kgap == "soap":
                update_soap(cdict)
            else:
                update_nbody(cdict)
            self.gaps[kgap] = cdict
            
        self.hessfit = "hessian_delta" in self.teach_sparse
        if gp_file is None:
            self.gp_file = gpfile_name(ogap.keys(), hessfit=self.hessfit)
        else:
            self.gp_file = gp_file

        #Configure the fitting directory for this particular set of
        #potentials. This way, we can get separate directories if the number of
        #n-body or soap potentials changes.
        from os import mkdir
        self.name = self.gp_file.replace(".xml", "")
        self.root = path.join(root, self.name)
        if not path.isdir(self.root):
            mkdir(self.root)

        #If we are doing hessian training, copy the seed configuration and
        #updates its hessian attributes.
        self.seeds = {}

        #Finally, determine where we are in the execution of the potential
        #fitting chain. `teach_sparse` only writes out a potential file when the
        #fitting is complete, before which it doesn't put out any files.
        self.state = _get_max_gap(self.root, self.hessfit)
        self.latest = gpfile_name(self.state, hessfit=self.hessfit)

        #We need to increment the state; the object assumes that state is what
        #we are *aiming for*, whereas now it only has what has finished
        #executing.
        if self.state is None:
            self.state = min(ogap.keys())
        else:
            diff = np.setdiff1d(ogap.keys(), self.state)
            if len(diff) > 0:
                n_next = min(diff)
                self.state.append(n_next)

        self._calculator = None
        """ase.Calculator: GAP potential for the latest fitted file in the
        sequence (i.e., based on :attr:`latest`).
        """
            
    @property
    def calculator(self):
        """Returns an instance of :class:`ase.Calculator` using the latest
        fitted GAP potential in this trainer.
        """
        if self._calculator is None:
            from quippy.potential import Potential
            from matdb.utility import chdir
            from os import getcwd
            with chdir(self.root):
                self._calculator = Potential("IP GAP", param_filename=self.latest)
        return self._calculator
            
    @property
    def delta_cache(self):
        """Returns the path the `delta` parameter file for each of the GAPs.
        """
        return path.join(self.root, "deltas.json")

    @property
    def validation_list(self):
        """Returns a :class:`quippy.AtomsList` of configurations that can be
        used for potential validation. The atoms list will depend on the type of
        fit (Hessian or not).
        """
        if self.hessfit:
            if len(self.seeds) == 0:
                self._get_hessians()
                
            #Use only those configurations derived from the seed configurations
            #that were actually included.
            result = quippy.AtomsList()
            for name in self.seeds:
                target = self.db.collections[name].holdout_file
                atoms = quippy.AtomsList(target)
                result.extend(atoms)

            return result
        else:
            #Return the global validation list from the whole database.
            return quippy.AtomsList(self.db.holdout_file)
    
    def validate(self):
        """Validates the latest calculator in this training sequence against the
        `holdout.xyz` configurations in the controller database.

        Returns:
            dict: with keys ["e_dft", "e_gap"] where the values are
            DFT-calculated and GAP-calculated predictions for the energies.
        """
        #We need to split to get validation data. If the split has already been
        #done as part of a different training run, then it won't be done a
        #second time.
        self.db.split(self.split)

        from tqdm import tqdm
        al = self.validation_list
        for a in tqdm(al):
            a.set_cutoff(self.calculator.cutoff())
            a.calc_connect()
            self.calculator.calc(a, energy=True, force=True)
    
        e_dft = np.array(al.dft_energy)
        e_gap = np.array(al.energy)
        f_dft = np.array(al.dft_force)
        f_gap = np.array(al.force)

        return {
            "e_dft": e_dft,
            "e_gap": e_gap,
            "f_dft": f_dft,
            "f_gap": f_gap
        }
    
    def _get_delta(self, igap):
        """Estimates the value to use for the `delta` parameter for the next GAP
        in the list.

        Args:
            igap (int): *index* of the n-body code for the GAP (SOAP=10) in
              :attr:`latest`.

        Returns:
            float: value for the `teach_sparse` delta parameter of the next GAP
            in the sequence.
        """
        #We cache the list of previous `delta` values for each GAP in a JSON
        #file so that we don't have to keep re-calculating them.
        import json
        if path.isfile(self.delta_cache):
            with open(self.delta_cache) as f:
                deltas = json.load(f)
        else:
            deltas = {}

        #See if we already have this delta in our list.
        ngap = self.state[igap]
        dname = self.ogap[ngap]
        if dname in deltas:
            return deltas[dname]
            
        #:attr:`latest` has the name of the GAP file with the latest potential
        #in the sequence. Use it to get predictions and errors for a set of
        #configurations.
        if igap > 0:
            ens = self.validate()
            rms = np.std(ens["e_dft"]-ens["e_gap"])
            msg.info("{0}: {1:.3f} RMS error.".format(self.latest, rms))
        else:
            self.db.split(self.split)
            
            al = quippy.AtomsList(self.db.train_file)
            ens = []
            for a in al:
                ens.append(a.dft_energy)
            rms = np.std(ens)
            msg.info("2B: {0:.3f} data variance.".format(rms))
        
        #We use the first seed configuration in the database and calculate how
        #many nearest neighbors it has within the specified cutoff of the n-body
        #configuration. SOAP doesn't rescale because the descriptors are already
        #per-atom.
        for name, dbset in self.db.collections.items():
            atoms = dbset.atoms
            break

        if "2body" == dname:
            scaling =  _rescale_2body(atoms, self.gaps["2body"])
        if "3body" == dname:
            scaling =  _rescale_3body(atoms, self.gaps["3body"])
        if "soap" == dname:
            scaling = 1.
            
        #Now, we scale the expected error by the number of descriptors per atom
        #in this particular GAP. Remember, the last element in state is the
        #*next* one in the sequence. But we need to rescale the one *before*
        #that, i.e., -2.
        deltas[dname] = rms/scaling
            
        #Now that we have updated the deltas for the latest GAP, re-serialize
        #the dict.
        with open(self.delta_cache, 'w') as f:
            json.dump(deltas, f)

        return deltas[dname]
            
    def _get_hessians(self):
        """Extracts the hessian matrix and assigns its eigenvalues and
        eigenvectors to relevant properties on copies of the seed configuration
        atoms objects.
        """
        from matdb.database.phonon import PhononBase
        from fnmatch import fnmatch
        
        hname = "hessian{}"        
        for name, dbset in self.db.collections.items():
            if self.configs is not None:
                #Skip all those seed configurations that are not explicitly
                #included in the training spec.
                if not any([fnmatch(name, p) for p in self.configs]):
                    continue
            
            self.seeds[name] = dbset.atoms.copy()
            dmatrix = dbset.databases[PhononBase.name].dmatrix
            eigvecs, eigvals = dmatrix["eigvecs"], dmatrix["eigvals"]
            atc = self.seeds[name]
            
            ni = 0
            for i, e in enumerate(eigvals):
                if not np.isclose(e, 0.0):
                    ni = ni + 1
                    hsingle = hname.format(ni)
                    atc.add_property(hsingle, 0.0, n_cols=3)
                    H = np.reshape(eigvecs[:,i], (atc.n, 3)).T
                    setattr(atc, hsingle, H)
                    atc.params.set_value(hsingle, e)
            if ni > 0:
                atc.params.set_value("n_hessian", ni)

            #Also add the zero forces. TODO: get the actual forces from the
            #relaxed DFT calculation.
            atc.add_property("force", 0.0, n_cols=3)
            
    def _write_hessian_trainxyz(self, recalc=False):
        """Writes the h_train.xyz file that is needed for the Hessian-based
        training program.

        Args:
            recalc (bool): when True, re-create the xyz file with the seed
              configurations and hessians.
        """
        self._get_hessians()
        outpath = path.join(self.root, self.train_file)
        if path.isfile(outpath) and not recalc:
            return

        import quippy.cinoutput as qcio
        try:
            out = qcio.CInOutputWriter(outpath)
            for a in self.seeds.values():
                a.write(out)
        finally:
            out.close()
                
    def _desc_str(self, ptype, **kwargs):
        """Returns the descriptor string for this training object, if it
        exists.

        Args:
            ptype (str): one of ["soap", "2body", "3body", "4body"].
            kwargs (dict): additional key-value pairs to include in the
              dictionary before compiling the string.
        """
        if ptype in self.gaps:
            pdict = self.gaps[ptype].copy()
            pdict.update(kwargs)

            if pdict.get("sparse_method") == "file":
                pdict["sparse_file"] = "soap_sparse_points.dat"
            
            settings = []
            settings.append("soap" if ptype == "soap" else "distance_Nb")
            settings.append(dict_to_str(pdict))
            settings.append("add_species")
            return ' '.join(settings)
                
    def _soap_sparse_points(self, recalc=False):
        """Calculates the hessian sparse points for hessian-based
        training.

        Args:
            recalc (bool): when True, recalculate the sparse points even if the
              file already exists.
        """
        spfile = path.join(self.root, "soap_sparse_points.dat")
        if path.isfile(spfile) and not recalc:
            return

        #Determine how many sparse points we need to generate.
        n_sparse = self.gaps["soap"].get("n_sparse")
        n_ratio = int(np.ceil(float(n_sparse)/len(self.seeds)))
        #Reset n_sparse to be this whole number.
        n_sparse = n_ratio * len(self.seeds)
        
        desc = quippy.Descriptor(self._desc_str("soap"))

        #Run each of the seed configurations once until we hit the ratio
        #limit.
        spoints = {}
        for name, atc in self.seeds.items():
            seed_pts = []
            n_seed = int(np.ceil(float(n_ratio)/atc.n))
            for ni in range(n_seed):
                atRand = atc.copy()
                quippy.randomise(atRand.pos, 0.2)
                atRand.set_cutoff(desc.cutoff())
                atRand.calc_connect()

                n_descriptors, n_cross = desc.descriptor_sizes(atRand)
                d = quippy.fzeros((desc.dimensions(), n_descriptors))
                desc.calc(atRand, d)
                d = np.array(d)
                
                #We need to add some ones to the SOAP descriptor; they represent
                #the relative weights of each of the sparse points. We weight
                #them all equally.
                seed_pts.append(np.vstack((np.ones(atRand.n), d)))
                
            spoints[name] = np.hstack(seed_pts)[:,0:n_ratio].T

        assert sum([v.shape[0] for v in spoints.values()]) == n_sparse
        allpoints = np.hstack([v.flatten() for v in spoints.values()])
        np.savetxt(spfile, allpoints)

    @property
    def train_file(self):
        """Returns the name of the XYZ training file that will be passed to the
        teach_sparse command.
        """
        if "hessian_delta" in self.teach_sparse:
            return "h_train.xyz"
        else:
            return "train.xyz"

    @property
    def jobfile(self):
        """Returns the name of the jobfile for the current training iteration.
        """
        if self.state is None:
            return "jobfile.sh"
        else:
            nmax = max(self.state)
            if nmax == 10:
                return "jobfile_soap.sh"
            else:
                return "jobfile_{}.sh".format(nmax)
        
    def execute(self, dryrun=False):
        """Submits the job script for the currently configured potential
        training.

        Args:
            dryrun (bool): when True, simulate the submission without
              actually submitting.

        Returns:
            bool: True if the submission generated a job id
            (considered successful).
        """
        if not path.isfile(path.join(self.root, self.jobfile)):
            return False

        #Make sure we have the `train.xyz` or `h_train.xyz` file; it is the only
        #dependency that we need for teach_sparse.
        trainfile = path.join(self.root, self.train_file)
            
        if not path.isfile(trainfile):
            from matdb.msg import std
            std("*train.xyz missing in {}; can't execute.".format(self.root))
            return False
        
        # We must have what we need to execute. Compile the command
        # and submit.
        from matdb.utility import execute
        cargs = ["sbatch", self.jobfile]
        if dryrun:
            from matdb.msg import okay
            okay("Executed {} in {}".format(' '.join(cargs), self.root))
            return True
        else:
            xres = execute(cargs, self.root)

        if len(xres["output"]) > 0 and "Submitted" in xres["output"][0]:
            from matdb.msg import okay
            okay("{}: {}".format(self.root, xres["output"][0].strip()))
            return True
        else:
            return False
            
    def write_jobfile(self):
        """Sets up the job file for the next GAP in the potential that needs to
        be executed. This includes creating folders, training XYZ files and
        setting up any other dependencies required by `teach_sparse`.

        .. note:: Calling the method assumes that you are using a newly
          constructed :class:`GAPTrainer`. The initialization of the class
          determines where in the fitting sequence of a `2 + 3 + ...` GAP
          potential we are. If you call this repeatedly with the same class in
          memory, it won't progress.

        .. note:: This method also runs :meth:`command`, which configures the
          directory that the executable will run in.
        """
        target = path.join(self.root, self.jobfile)
        if path.isfile(target):
            return
        
        # We use the global execution parameters and then any updates
        # locally. We need to add the execution directory (including prefix) and
        # the number of jobs in the array.
        settings = self.execution.copy()
        settings["execution_path"] = self.root
        settings["exec_path"] = self.command()

        from jinja2 import Environment, PackageLoader
        env = Environment(loader=PackageLoader('matdb', 'templates'))
        template = env.get_template(settings["template"])
        with open(target, 'w') as f:
            f.write(template.render(**settings))

    def command(self):
        """Returns the `teach_sparse` command that is needed to train the GAP
        potentials specified by this object.

        .. note:: This method also configures the directory that the command
          will run in so that it has the relevant files.
        """
        #First, tell the DB to split over training and holdout.
        if self.hessfit:
            #Generate the sparse points file and the seed configuration XYZ
            #training file.
            self._write_hessian_trainxyz()
            if ("soap" in self.gaps and
                self.gaps["soap"].get("sparse_method") == "file"):
                self._soap_sparse_points()
        else:
            #The split command writes the XYZ files using a new random subset of
            #the available configurations.
            self.db.split(self.split)
        
        template = "teach_sparse at_file={train_file} gap={{{gap}}} {teach_sparse}"
        #Compile the GAP string specifying the potential.
        gaplist = []
        for i, igap in enumerate(self.state):
            kgap = self.ogap[igap]
            delta = self._get_delta(i)
            if igap < 10:
                gaplist.append(self._desc_str(kgap, delta=delta, order=igap))
            else:
                gaplist.append(self._desc_str(kgap, delta=delta))
            
        gapstr = ':'.join(gaplist)

        tsattrs = {
            "e0": self.e0,
            "default_sigma": self.default_sigma,
            "sparse_jitter": self.sparse_jitter,
            "gp_file": gpfile_name(self.state, hessfit=self.hessfit)
        }
        tsattrs.update(self.teach_sparse)
        tsstr = dict_to_str(tsattrs)
                
        fields = {
            "train_file": self.train_file,
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
            "latest": self.latest,
            "jobfile": path.isfile(path.join(self.root, self.jobfile))
        }
        
        if printed:
            msg.info("Latest GAP model: {}".format(result["latest"]))
            x = "exists" if result["jobfile"] else "does not exist"
            msg.info("Next jobfile '{}' {}".format(self.jobfile, x))
        else:
            return result
