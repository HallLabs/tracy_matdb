"""Implements a hessian database group by extracting eigenvectors and
eigenvalues of the force constants (aka Hessian matrix) for each of the
configurations in `seed`.
"""
from os import path, remove
import numpy as np
from operator import itemgetter
from matdb.database import Group
from matdb.utility import execute, chdir
from phonopy import file_IO
from phonopy.api_phonopy import Phonopy
from phonopy.structure.atoms import Atoms as PhonopyAtoms
from phonopy.cui.phonopy_argparse import get_parser
from phonopy.cui.settings import PhonopyConfParser
from matdb import msg
from matdb.atoms import Atoms, AtomsList
from matdb.transforms import conform_supercell

def roll_fc(hessian):
    """Rolls the specified hessian into the `phonopy` force constants format.
    
    Args:
        hessian (numpy.ndarray): of shape `n_atoms * 3`.
    """
    n = hessian.shape[0]/3
    result = np.zeros((n, n, 3, 3), dtype='double')
    
    for i in range(n):
        for j in range(n):
            result[i, j] = hessian[i*3:(i+1)*3, j*3:(j+1)*3]

    return result

def unroll_fc(fc):
    """Unroll's the phonopy force constants matrix into the Hessian.
    """
    Nc = int(np.sqrt(np.product(fc.shape)))
    result = np.zeros((Nc, Nc), dtype="double")
    num_atom = Nc/3
    for i in range(num_atom):
        for j in range(num_atom):
            result[i*3:(i+1)*3, j*3:(j+1)*3] = fc[i, j]
  
    return result

def phonopy_to_matdb(patoms):
    """Converts a :class:`phonopy.structure.atoms.Atoms` to
    :class:`matdb.atoms.Atoms`. See also :func:`matdb_to_phonopy`.
    """
    #We hard-code the pbc here because phonons only make sense for solids.
    return Atoms(numbers=patoms.get_atomic_numbers(),
                 positions=patoms.get_positions(),
                 magmoms=patoms.get_magnetic_moments(),
                 cell=patoms.get_cell(),
                 pbc=[True, True, True])

def matdb_to_phonopy(matoms):
    """Converts a :class:`matdb.atoms.Atoms` to a
    :class:`phonopy.structure.atoms.Atoms`. See also :func:`phonopy_to_matdb`.
    """
    return PhonopyAtoms(numbers=matoms.get_atomic_numbers(),
                        positions=matoms.positions,
                        masses=matoms.get_masses(),
                        cell=matoms.cell, pbc=matoms.pbc)

class Hessian(Group):
    """Sets up the displacement calculations needed to construct the Hessian
    matrix.
    Args:
        atoms (matdb.atoms.Atoms): seed configuration that will be
          displaced to generate the database.
        root (str): path to the folder where the database directories will
          be stored.
        parent (matdb.database.Database): parent sequence to which this database
          belongs. Could also be another :class:`Hessian`.
        phonons (dict): specifying additional settings for `phonopy`
          configuration files (i.e., differing from, or in addition to those in
          the global set).
        bandmesh (list): of `int`; number of splits in each reciprocal
          lattice vector according to Monkhorst-Pack scheme. Used for
          calculating the phonon bands.
        dosmesh (list): of `int`; number of splits in each reciprocal
          lattice vector according to Monkhorst-Pack scheme. Used for
          calculating the phonon density-of-states.
        tolerance (float): maximum difference in integrated total DOS that may
          exist for a calculation to be considered "good enough". This is used
          when selecting the "best" calculation from a parameter grid.
        dfpt (bool): when True, calculate the force constants using Density
          Functional Perturbation Theory.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
        name (str): name of this database type relative to the over database
          collection. This is also the name of the folder in which all of its
          calculations will be performed.
        supercell (list): of `int`; number of cells in each direction for
          generating the supercell.
        phonodir (str): directory in which all `phonopy` executions take place.
        bandmesh (list): mesh for calculating the phonon bands.
        dosmesh (list): mesh for calculating the phonon density-of-states.
    """
    splittable = False
    
    def __init__(self, phonopy={}, name="hessian", bandmesh=None,
                 dosmesh=None, tolerance=0.1, dfpt=False, **dbargs):
        self.name = name
        self.seeded = True
        dbargs["prefix"] = "H"
        dbargs["cls"] = Hessian
        if "Hessian" not in dbargs['root']:
            from os import mkdir
            new_root =path.join(dbargs['root'],"Hessian")
            if not path.isdir(new_root):
                mkdir(new_root)
            dbargs['root'] = new_root
        
        #Make sure that we override the global calculator default values with
        #those settings that we know are needed for good phonon calculations. We
        #also need to assign the parent here *before* the super.__init__ so that
        #the recursive database lookup works.
        self.dfpt = dfpt
        self.parent = dbargs["parent"]
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs["calculator"]:
                calcargs.update(dbargs["calculator"])
        self._set_calc_defaults(calcargs)
        dbargs["calculator"] = calcargs

        #We only initialize now because we wanted to first update the dbargs for
        #the calculator, etc.
        super(Hessian, self).__init__(**dbargs)
        
        if "dim" in phonopy:
            self.supercell = phonopy["dim"]

            #Make sure that the supercell matrix has positive determinant for
            #phonopy; if it doesn't correct it.
            scell = conform_supercell(self.supercell)
            det = np.linalg.det(scell)
            if det < 0:
                self.supercell = list(np.array(self.supercell)*-1)
        else:
            self.supercell = None
            
        self.bandmesh = bandmesh
        self.dosmesh = dosmesh
        self.tolerance = tolerance
        self.phonodir = path.join(self.root, "phonopy")
        self.phonocache = path.join(self.root, "phoncache")
        self._kpath = None
        """tuple: Special point path in k-space. First term is a list of special point
        labels; second is the list of points corresponding to those labels.
        """

        self._bands = None
        """dict: keys are ['q', 'w', 'path', 'Q'], values are the distance along
        the special path (scalar), phonon frequencies at that distance (vector,
        one component for each frequency), q-positions of the special points
        along the paths, and their corresponding distances.
        """

        self._bzsample = None
        """tuple: returned by :attr:`BZ_sample`."""
        
        self._H = self.load_pkl(self.H_file)
        """np.array: the Hessian matrix, whether it was derived from
        DFPT or from frozen phonon calculations.
        """
        # Only place these directories if we're at the bottom of the stack.
        if self.pgrid is None or (self.pgrid is not None and len(self.pgrid) ==0):
            from os import mkdir
            if not path.isdir(self.phonodir):
                mkdir(self.phonodir)
            if not path.isdir(self.phonocache):
                mkdir(self.phonocache)

    @property
    def fitting_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group. This list includes a single *duplicated* configuration for each
        of the eigenvalue/eigenvector combinations of the Hessian matrix.
        """
        if len(self.sequence) == 0:
            if self.ready():
                return self.config_atoms.values()
            else:
                return AtomsList()
        else:
            result = AtomsList()
            for g in self.sequence.values():
                result.extend(g.fitting_configs)
            return result

    def _hessian_configs(self):
        """Returns a :class:`matdb.atoms.AtomsList` for all configs in this
        group. This list includes a single *duplicated* configuration for each
        of the eigenvalue/eigenvector combinations of the Hessian matrix.

        .. note:: This assumes that the group has actual displacements and is
          not a parent group in the recursive structure.
        """
        configs = AtomsList()

        #Start with a configuration that has energy, force and virial
        #information from the VASP computation. We just grab the first of the
        #config_atoms from this sequence or its children.
        atBase = next(self.iconfigs)
        atcalc = atBase.get_calculator()
        
        #Make sure that energy, force and virial information was found. 
        assert getattr(atBase, atcalc.energy_name) < 0
        assert getattr(atBase, atcalc.force_name).shape == (atBase.n, 3)
        assert getattr(atBase, atcalc.virial_name).shape == (3, 3)
        configs.append(atBase)

        #Now, make a copy of the base atoms object; this object only has the
        #lattice and positions. We will add the eigenvalue and eigenvectors
        #*individually* because they each need different scaling for sigma in
        #the GAP fit.
        atEmpty = atBase.copy()
        for k in list(atBase.params.keys()):
            if "energy" in k or "virial" in k:
                atEmpty.rm_param(k)
        for k in list(atBase.properties.keys()):
            if "force" in k:
                atEmpty.rm_property(k)

        hname = "{}_hessian1".format(self.calc.key)
        #NB: make sure you transpose the eigenvectors matrix before doing the
        #zip!
        evals, evecs = np.linalg.eigh(self.H)
        natoms = len(evals)/3
        l0 = np.max(evals)-np.min(evals)
        sigma0 = 0.01
        lscaling = 0.005

        for l, v in zip(*(evals, evecs.T)):
            #The eigenvalues should all be positive. There may be some really
            #small ones that are essentially zero, but slightly negative.
            if np.abs(l) < 1e-5 or l < 0:
                continue
                    
            #Add this eigenvector to its own configuration.
            atc = atEmpty.copy()
            Hi = np.reshape(v, (natoms, 3))
            atc.properties[hname] =  Hi
                    
            #Same thing for the eigenvalue.
            atc.add_param(hname, l)

            #We want the small eigenvalues to have a weighting of sigma0 and the
            #largest eigenvalue to have a sigma of 10% of its value.
            c = (lscaling*l0-sigma0)/l0**2
            #atc.add_param("{}_hessian_csigma".format(self.calc.key), sigma0 + c*l**2)
            atc.add_param("{}_hessian_csigma".format(self.calc.key), sigma0)
            atc.add_param("n_{}_hessian".format(self.calc.key), 1)
            configs.append(atc)

        return configs

    def sub_dict(self):
        """Returns a dict needed to initialize the class.
        """
        args = {"phonopy":{"dim":self.supercell},"name":self.name,
                "bandmesh":self.bandmesh,"dosmesh":self.dosmesh,
                "tolerance":self.tolerance,"dfpt":self.dfpt}
        return args

    @property
    def BZ_sample(self):
        """Returns a full sampling of the BZ by calculating frequencies at every
        unique k-point as sampled on the :attr:`dosmesh` grid.
        
        Returns:
            tuple: `(q-points, weights, eigenvalues)` where the `q-points` are
            the unique points after symmetry reduction; `weights` are the
            corresponding weights for each point; `eigenvalues` is a
            :class:`numpy.ndarray` of frequencies (in THz) for each point.
        """
        #This is a little convoluted because of how the phonopy API works. We
        #want to get a full sampling of the BZ, but with symmetry, and then
        #compare the eigenvalues at every point.
        if self._bzsample is None:
            with chdir(self.phonodir):
                atoms = matdb_to_phonopy(self.atoms)
                phonpy = Phonopy(atoms, np.array(self.supercell).reshape(3,3))
                phonpy.set_force_constants(roll_fc(self.H))
                phonpy._set_dynamical_matrix()       

                #Phonopy requires full settings to compute the unique grid and
                #eigenvalues. We spoof the command-line parser.
                parser = get_parser()
                (options, args) = parser.parse_args()
                option_list = parser.option_list
                options.mesh_numbers = ' '.join(map(str, self.dosmesh))
                phonopy_conf = PhonopyConfParser(options=options, option_list=option_list)
                settings = phonopy_conf.get_settings()

                #Next, set the mesh on the phonopy object and ask it to reduce and
                #calculate frequencies.
                mesh = settings.get_mesh()
                phonpy.set_mesh(*mesh)
                self._bzsample = phonpy.get_mesh()[0:3]

        return self._bzsample
    
    def compare(self, other):
        """Compares this Hessian group's calculated bands with another Hessian
        group that is taken to the "right" answer.
        Args:
            other (Hessian): correct bands to compare to.
        Returns:
            tuple: `(abs err, rms err)` for the eigenvalues at every point in the
            sampling of the BZ.
        """
        dosmesh = self.dosmesh
        self.dosmesh = other.dosmesh
        sgrid, sweights, sfreqs = self.BZ_sample
        ogrid, oweights, ofreqs = other.BZ_sample

        #Make sure both are running on the same grid. This should always be true
        #unless the primitive is different, or if the grid spacing is different.
        assert np.allclose(sgrid, ogrid)
        
        return (np.mean(np.abs(sfreqs-ofreqs)), np.std(sfreqs-ofreqs))
    
    def _best_bands(self):
        """Returns the name of the band collection that has the smallest *converged*
        phonon bands. This is accomplished by assuming that the largest supercell is
        the "correct" answer, and comparing the total DOS. If the comparitive error
        is within `tolerance`, then it is acceptable. The smallest acceptable
        supercell's key is returned.
        Returns:
            str: the key in the group's sequence that has the smallest acceptable
            supercell size.
        """
        #Find the cell size and DOS for each calculation in the sequence.
        sizes = {k: np.linalg.det(np.reshape(np.array(d.supercell), (3, 3)))
                 for k, d in self.sequence.items()}
        dos = {k: np.loadtxt(d.dos_file)
               for k, d in self.sequence.items()}

        #Find the calculation with the largest cell size and grab its DOS.
        maxkey, maxval = max(sizes.items(), key=itemgetter(1))
        maxdos = dos[maxkey]

        ok = {}
        for k, d in self.sequence.items():
            if k == maxkey:
                continue
            assert dos[k].shape == maxdos.shape
            diff = np.sum(np.abs(dos[k][:,1]-maxdos[:,1]))
            if diff < self.tolerance:
                ok[k] = sizes[k]

        #Now, choose the supercell with the smallest cell size, if one
        #exists. Otherwise warn the user that either the tolerance was too low, or
        #that the calculation may not be converged.
        if len(ok) > 0:
            minkey, minval = min(ok.items(), key=itemgetter(1))
        else:
            msg.warn("Hessian calculation may not be converged. Your tolerance "
                     "may be too high. Returning the largest supercell by default.")
            minkey = maxkey

        return minkey

    @property
    def rset(self):
        """Constructs the Hessian matrix for the *best* convergence parameters
        in this group and it's possible sub-sequences.
        Returns:
            list: of :class:`matdb.atoms.Atoms`; each atoms object will have a `H`
            matrix in its info dictionary.
        """
        if len(self.sequence) == 0:
            #We are at the bottom of the stack; attach the hessian matrix
            #to the atoms object it corresponds to.
            self.atoms.info["H"] = self.H
            result = AtomsList()
            result.append(self.atoms)
            return result
        else:
            #Check where we are in the stack. If we are just below the database,
            #then we want to return a list of hessian matrices and atoms
            #objects. If we are not, then we must a parameter grid of sequences
            #to select from.
            if isinstance(self.parent, Hessian):
                #We have many dynamical matrices to choose from. We need to decide
                #what "best" means and then return that one.
                bestkey = self._best_bands()
                return self.sequence[bestkey].rset
            else:
                result = AtomsList()
                for p in self.sequence.values():
                    result.extend(p.rset)
                return result
            
    def _set_calc_defaults(self, calcargs):
        """Sets the default calculator parameters for phonon calculations based
        on the calculator specified in `calcargs`.
        .. warning:: This method mutates the `calcargs` dictionary.
        Args:
            calcargs (dict): the "calculator" dictionary that is part of the
              group arguments for db group.
        """
        from matdb.calculators import get_calculator_module
        try:
            mod = get_calculator_module(calcargs)
            if mod is not None and hasattr(mod, "phonon_defaults"):
                call = getattr(mod, "phonon_defaults")
                call(calcargs, self.dfpt)
        except:
            pass

    @property
    def dos_file(self):
        """Returns the full path to file that contains the total DOS for the
        phonons.
        """
        return path.join(self.phonodir, "total_dos.dat")
    
    def ready(self):
        """Returns True if all the phonon calculations have been completed, the
        force sets have been created, and the DOS has been calculated.
        """
        self._expand_sequence()
        if len(self.sequence) == 0:
            #If the DOS has been calculated, then all the other steps must have
            #completed correctly.
            result = path.isfile(self.dos_file)
            if not result:
                msg.std("{} is not ready. Exiting.".format(self.root), 2)
            return result
        else:
            ready = False
            for p in self.sequence.values():
                if not p.ready():
                    msg.std("{} is not ready. Exiting.".format(p.root), 2)
                    break
            else:
                ready = True
            return ready

    @property
    def H_file(self):
        """Returns the full path to the hessian matrix pickle file.
        """
        return path.join(self.root, "H.pkl")

    @property
    def H(self):
        """Returns the Hessian matrix extracted from the frozen phonon calculations.
        """
        if self._H is not None:
            return self._H
        
        #Otherwise, we need to calculate it from scratch. This depends on
        #whether we are using DFPT or frozen phonons.
        if not self.dfpt:
            self.calc_forcesets()
            dim = ' '.join(map(str, self.supercell))
            xargs = ["phonopy", '--dim="{}"'.format(dim), "--writefc"]
            execute(xargs, self.phonodir, venv=True)

            if not path.isfile(path.join(self.phonodir, "FORCE_CONSTANTS")):
                msg.err("Cannot convert FORCE_SETS to FORCE_CONSTANTS")
                msg.err(''.join(xargs["output"]))
        else:
            self.calc_fc()
        
        with chdir(self.phonodir):
            result = file_IO.parse_FORCE_CONSTANTS()
            
        self._H = unroll_fc(result)        
        self.save_pkl(self._H, self.H_file)
        return self._H
    
    @property
    def bands(self):
        """Returns the DFT-accurate phonon bands in a format that can be
        consumed by the `matdb` interfaces.
        Returns:
            dict: keys are ['q', 'w', 'path', 'Q'], values are the distance along
            the special path (scalar), phonon frequencies at that distance (vector,
            one component for each frequency), q-positions of the special points
            along the paths, and their corresponding distances.
        """
        if self._bands is None:
            from matdb.phonons import from_yaml
            byaml = path.join(self.phonodir, "band.yaml")
            self._bands = from_yaml(byaml)
        return self._bands

    def get_kpath(self):
        """Returns the special path in k-space for the seed configuration of this group.
        Returns:
            tuple: special path in k-space. First term is a list of special
            point labels; second is the list of points corresponding to those
            labels.
        """
        if self._kpath is None:
            from matdb.kpoints import parsed_kpath
            self._kpath = parsed_kpath(self.atoms)
            
        return self._kpath    
    
    @property
    def kpath(self):
        """Returns the special path in k-space for the seed configuration of this
        database.
        Returns:
            tuple: special path in k-space. First term is a list of special
            point labels; second is the list of points corresponding to those
            labels.
        """
        from matdb.database import Database
        if isinstance(self.parent, Database):
            return self.get_kpath()
        elif isinstance(self.parent, Hessian):
            return self.parent.get_kpath()
    
    def calc_bands(self, recalc=False):
        """Calculates the bands at the special points given by seekpath.
        
        Args:
            recalc (bool): when True, recalculate the DOS, even if the
              file already exists.
        """
        self._expand_sequence()
        bandfile = path.join(self.phonodir, "band.yaml")
        if not recalc and path.isfile(bandfile):
            return
        
        from matdb.phonons import _calc_bands
        _calc_bands(self.atoms, self.H, supercell=self.supercell,
                    outfile=bandfile, grid=self.bandmesh)
    
    def calc_DOS(self, recalc=False):
        """Calculates the *total* density of states.
        Args:
            recalc (bool): when True, recalculate the DOS, even if the
              file already exists.
        """
        dosfile = path.join(self.phonodir, "mesh.yaml")
        if not recalc and path.isfile(dosfile):
            return

        #Make sure we have calculated the force sets already.
        self.calc_forcesets(recalc)
        settings = {
            "ATOM_NAME": ' '.join(self.database.parent.species),
            "DIM": ' '.join(map(str, self.supercell)),
            "MP": ' '.join(map(str, self.dosmesh))
        }
        with open(path.join(self.phonodir, "dos.conf"), 'w') as f:
            for k, v in settings.items():
                f.write("{} = {}\n".format(k, v))

        sargs = ["phonopy", "-p", "dos.conf", "-s"]
        xres = execute(sargs, self.phonodir, venv=True)
        #Make sure that phonopy actually produced files; otherwise show the output
        #(phonopy doesn't write to stderr, only stdout).
        if not path.isfile(dosfile): #pragma: no cover
            msg.std(''.join(xres["error"]))
            msg.err("could not calculate the DOS; see errors.")

    def calc_fc(self, recalc=False):
        """Extracts the force constants from a DFPT Hessian matrix.
        """
        fcfile = path.join(self.phonodir, "FORCE_CONSTANTS")
        if not recalc and path.isfile(fcfile):
            return

        from matdb.calculators import get_calculator_module
        mod = get_calculator_module(self.calcargs)
        call = getattr(mod, "extract_force_constants")
        xres = call(self.configs, self.phonodir)
        
        #Make sure that phonopy actually produced files; otherwise show the
        #output (phonopy doesn't write to stderr, only stdout).
        if not path.isfile(fcfile): #pragma: no cover
            msg.std(''.join(xres["error"]))
            msg.err("could not calculate the force constants from DFPT.")
            
    def calc_forcesets(self, recalc=False):
        """Extracts the force sets from the displacement calculations.
        Args:
            recalc (bool): when True, recalculate the force sets, even if the
              file already exists.
        """
        fsets = path.join(self.phonodir, "FORCE_SETS")
        if not recalc and path.isfile(fsets):
            return
        
        from matdb.calculators import get_calculator_module
        mod = get_calculator_module(self.calcargs)
        call = getattr(mod, "extract_force_sets")
        xres = call(self.configs, self.phonodir)

        #Make sure that phonopy actually produced files; otherwise show the output
        #(phonopy doesn't write to stderr, only stdout).
        if not path.isfile(fsets):#pragma: no cover
            msg.std(''.join(xres["output"]))
            msg.err("Couldn't create the FORCE_SETS in {}.".format(self.phonodir))

    def setup(self, rerun=0):
        """Displaces the seed configuration preparatory to calculating the force
        sets for phonon spectra.
        Args:
            rerun (int): when > 1, recreate the folders even if they
              already exist. If > 0, then recreate the jobfile.
        """        
        super(Hessian, self).setup(self._setup_configs, rerun)
            
    def _setup_configs(self, rerun):
        """Displaces the seed configuration preparatory to calculating the force
        sets for phonon spectra.

        .. note:: This method *appears* to be VASP-specific. However, the
          configurations that are generated by `phonopy` as VASP `POSCAR` files
          are turned into :class:`matdb.atoms.Atoms` objects before they are
          passed to the super class that sets up the actual calculations. So, it
          is quite general.

        Args:
            rerun (int): when > 1, recreate the folders even if they
              already exist. If > 0, recreate the jobfile.
        """        
        #We also don't want to setup again if we have the results already.
        if self.ready() and rerun == 0:
            return

        if not self.is_setup() or rerun > 1:
            from ase.io import write
            write(path.join(self.phonodir, "POSCAR"), self.atoms, "vasp")        
            scell = ' '.join(map(str, self.supercell))
            sargs = ["phonopy", "-d", '--dim="{}"'.format(scell)]
            pres = execute(sargs, self.phonodir, venv=True)

            #Make sure that phonopy produced the supercell. If it didn't, it
            #should have printed an error to *stdout* because it doesn't know
            #about stderr...
            if not path.isfile(path.join(self.phonodir, "SPOSCAR")):
                msg.err('\n'.join(pres["output"]))
            
            from matdb.atoms import Atoms
            if not self.dfpt:
                #Frozen phonons, create a config execution folder for each of
                #the displacements.
                from glob import glob
                with chdir(self.phonodir):
                    for dposcar in glob("POSCAR-*"):
                        dind = int(dposcar.split('-')[1])
                        datoms = Atoms(dposcar, format="vasp")
                        self.create(datoms)
            else:
                #Pull the perfect supercell and set it up for executing with
                #DFPT parameters.
                with chdir(self.phonodir):
                    datoms = Atoms("SPOSCAR", format="vasp")
                    self.create(datoms)

        # Last of all, create the job file to execute the job array.
        self.jobfile(rerun)

    def extract(self, recalc=False, cleanup="default"):
        """Runs post-DFT execution routines to calculate the force-sets and the
        density of states.
        Args:
            recalc (bool): when True, redo any calculations that use the DFT
              outputs to find other quantities.
            cleanup (str): the level of cleanup to perform after extraction.
        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.
        """
        if not super(Hessian, self).extract(cleanup=cleanup):
            return

        if len(self.configs) > 0:
            if not self.dfpt:
                self.calc_forcesets(recalc)
            else:
                self.calc_fc(recalc)
                
            self.calc_DOS(recalc)
            
        return self.ready()
