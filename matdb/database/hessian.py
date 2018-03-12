"""Implements a hessian database group by extracting eigenvectors and
eigenvalues of the force constants (aka Hessian matrix) for each of the
configurations in `seed`.
"""
from os import path, remove
import numpy as np
from operator import itemgetter
from .basic import Group
from matdb.utility import execute, chdir
from phonopy import file_IO
        
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
    def __init__(self, phonopy={}, name="dynmatrix", bandmesh=None,
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
        super(Hessian, self).__init__(**dbargs)
        
        #Make sure that we override the global calculator default values with
        #those settings that we know are needed for good phonon calculations.
        calcargs = self.database.calculator.copy()
        if "calculator" in dbargs:
            if dbargs["calculator"] is not None and "name" in dbargs["calculator"]:
                calcargs.update(dbargs["calculator"])
                self._set_calc_defaults(calcargs)
                dbargs["calculator"] = calcargs            

        if "dim" in phonopy:
            self.supercell = phonopy["dim"]
        else:
            self.supercell = None
            
        self.dfpt = dfpt
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
        configs = quippy.AtomsList()

        #Start with the configuration that has energy, force and virial
        #information from the VASP computation. This is the result of relaxing
        #the pure structure before we created our supercells. Get this from the
        #relax step in the previous group.
        #TODO: implement the relaxation Group...
        atBase = self.atoms.copy()

        #Make sure that energy, force and virial information was found. 
        assert atBase.energy < 0
        assert atBase.force.T.shape[1] == 3
        assert atBase.virial.T.shape == (3, 3)
        configs.append(atBase)

        #Now, make a copy of the base atoms object, delete the force, energy and
        #virial information so that we have just the lattice and positions. We
        #will add the eigenvalue and eigenvectors *individually* because they
        #each need different scaling for sigma in the GAP fit.
        atEmpty = atBase.copy()
        for k in atEmpty.params:
            if "energy" in k or "virial" in k:
                del atEmpty.params[k]
        for k in atEmpty.properties:
            if "force" in k:
                del atEmpty.properties[k]

        hname = "{}_hessian1".format(self.calc.key)
        #NB: make sure you transpose the eigenvectors matrix before doing the
        #zip!
        evals, evecs = np.linalg.eigh(self.H)
        natoms = len(evals)/3
        for l, v in zip(*(evals, evecs.T)):
            #The eigenvalues should all be positive. There may be some really
            #small ones that are essentially zero, but slightly negative.
            if np.abs(l) < 1e-5 or l < 0:
                continue
                    
            #Add this eigenvector to its own configuration.
            atc = atEmpty.copy()
            Hi = np.reshape(v, (natoms, 3))
            atc.arrays[hname] = Hi
                    
            #Same thing for the eigenvalue.
            atc.params.set_value(hname, l)
            
            #This custom scaling reweights by eigenvalue so that larger
            #eigenvalues get fitted more closely. The 0.1 is our "default_sigma"
            #for hessian.
            atc.params.set_value("hessian_csigma", 0.1/np.sqrt(l))
            atc.params.set_value("n_hessian", 1)
            configs.append(atc)

        return configs
                
    def sub_dict(self):
        """Returns a dict needed to initialize the class.
        """
        args = {"phonopy":{"dim":self.supercell},"name":self.name,
                "bandmesh":self.bandmesh,"dosmesh":self.dosmesh,
                "tolerance":self.tolerance,"dfpt":self.dfpt}
        return args
                
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
            msg.warn("DynMatrix calculation may not be converged. Your tolerance "
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
            return [self.atoms]
        else:
            #Check where we are in the stack. If we are just below the database,
            #then we want to return a list of hessian matrices and atoms
            #objects. If we are not, then we must a parameter grid of sequences
            #to select from.
            if isinstance(self.parent, DynMatrix):
                #We have many dynamical matrices to choose from. We need to decide
                #what "best" means and then return that one.
                bestkey = self._best_bands()
                return self.sequence[bestkey].rset
            else:
                return [p.rset for p in self.sequence.values()]
            
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
        if len(self.sequence) == 0:
            #If the DOS has been calculated, then all the other steps must have
            #completed correctly.
            return path.isfile(self.dos_file)
        else:
            return all(p.ready() for p in self.sequence.values())

    @property
    def H_file(self):
        """Returns the full path to the hessian matrix pickle file.
        """
        return path.join(self.root, "H.pkl")

    @property
    def H(self):
        """Returns the Hessian matrix extracted from the frozen phonon calculations at
        the gamma point in the BZ.
        """
        if self._H is not None:
            return self._H
        
        #Otherwise, we need to calculate it from scratch. This depends on
        #whether we are using DFPT or frozen phonons.


        if not self.dfpt:
            dim = ' '.join(map(str, self.supercell))
            xargs = ["phonopy", '--dim="{}"'.format(dim), "--writefc"]
            execute(xargs, self.phonodir, venv=True)    
        
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
            self._kpath = parsed_kpath
            
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
        from matdb.database import Database, Group
        if isinstance(self.parent, Database):
            return self.get_kpath()
        elif isinstance(self.parent, DynMatrix):
            return self.parent.get_kpath()
    
    def calc_bands(self, recalc=False):
        """Calculates the bands at the special points given by seekpath.
        
        Args:
            recalc (bool): when True, recalculate the DOS, even if the
              file already exists.
        """
        bandfile = path.join(self.phonodir, "band.yaml")
        if not recalc and path.isfile(bandfile):
            return
        
        from matdb.phonons import _calc_bands
        _calc_bands(self.atoms, self.H, supercell=self.supercell,
                    outfile=self.bandfile, grid=self.bandmesh)
    
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
            msg.std(''.join(xres["error"]))
            msg.err("Couldn't create the FORCE_SETS.")

    def setup(self, rerun=False):
        """Displaces the seed configuration preparatory to calculating the force
        sets for phonon spectra.
        Args:
            rerun (bool): when True, recreate the folders even if they
              already exist. 
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
            rerun (bool): when True, recreate the folders even if they
              already exist. 
        """        
        #We also don't want to setup again if we have the results already.
        if self.ready():
            return

        if not self.is_setup():
            from ase.io import write

            write(path.join(self.phonodir, "POSCAR"), self.atoms, "vasp")        
            scell = ' '.join(map(str, self.supercell))
            sargs = ["phonopy", "-d", '--dim="{}"'.format(scell)]
            pres = execute(sargs, self.phonodir, venv=True)

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

    def cleanup(self, recalc=False):
        """Runs post-DFT execution routines to calculate the force-sets and the
        density of states.
        Args:
            recalc (bool): when True, redo any calculations that use the DFT
              outputs to find other quantities.
        Returns:
           bool: True if the database is ready; this means that any other
           databases that rely on its outputs can be run.
        """
        if not super(Hessian, self).cleanup():
            return

        if self.dfpt:
            self.calc_forcesets(recalc)
        else:
            self.calc_fc(recalc)
            
        self.calc_DOS(recalc)
        return self.ready()
