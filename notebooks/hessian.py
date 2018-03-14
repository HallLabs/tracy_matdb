"""Implements a hessian database by extracting eigenvectors and eigenvalues of
the dynamical matrix for each of the configurations in `seed`.
"""
from phonopy.structure.atoms import Atoms as PhonopyAtoms
from phonopy.api_phonopy import Phonopy
from phonopy import file_IO
from phonopy.structure.cells import get_supercell
from phonopy.phonon.modulation import Modulation
import numpy as np
import quippy
import ase
from os import path

from .basic import Group
from matdb.utility import chdir

def phonopy_to_ase(patoms):
    """Converts a :class:`phonopy.structure.atoms.Atoms` to :class:`ase.Atoms`.
    """
    return ase.Atoms(symbols=patoms.get_chemical_symbols(),
                     positions=patoms.get_positions(),
                     magmoms=patoms.get_magnetic_moments(),
                     cell=patoms.get_cell())

class HessianSupercell(object):
    """Represents a supercell configuration generator that has a set of eigenvalues
    and eigenvectors compatible with Hessian fitting.

    Args:
        primitive (ase.Atoms): primitive configuration to create the supercell
          from.
        supercell (np.array): of `int` supercell matrix.
        folder (str): path to the `phonopy` folder where `FORCE_SETS` and
          displacements are kept; or where `vasprun.xml` is located when the
          full Hessian is calculated using VASP.

    Attributes:
        phonopy (phonopy.Phonopy): phonopy class for generating dynamical
          matrix, eigenvalues and eigenvectors.
    """
    def __init__(self, primitive, supercell, folder, vasp=False):
        self.atoms = PhonopyAtoms(symbols=primitive.get_chemical_symbols(),
                                  positions=primitive.positions,
                                  masses=primitive.get_masses(),
                                  cell=primitive.cell, pbc=primitive.pbc)

        self.supercell = supercell
        self.folder = folder
        if not vasp:
            self.phonopy = Phonopy(self.atoms, supercell)
            self._get_dynmatrix()
            self.primitive = self.phonopy._dynamical_matrix.get_primitive()
        else:
            from matdb.io import vasp_to_xyz
            self.phonopy = Phonopy(self.atoms, supercell)
            if path.isfile("FORCE_CONSTANTS"):
                fc = file_IO.parse_FORCE_CONSTANTS(filename="FORCE_CONSTANTS")
                self.phonopy.set_force_constants(fc)
                self.phonopy._set_dynamical_matrix()

            self.primitive = primitive
            vasp_to_xyz(folder)
            self.parent = quippy.Atoms(path.join(folder, "output.xyz"))
        
    def diagonalize(self, q):
        """Diagonalizes the dynamical matrix at `q`.

        Args:
            q (numpy.array): q-point to diagonalize at.

        Returns:
            tuple: `(eigvals, eigvecs)`, where the eigenvalues are in units of
            `THz` and both eigenvalues and eigenvectors are in the *primitive*
            cell.
        """
        return self.phonopy.get_frequencies_with_eigenvectors(q)

    def _get_phase_factor(self, modulation, argument):
        """Returns a phase factor corresponding to the given modulation.

        Args:
            modulation (numpy.ndarray): list of displacements, one for each
              atom.
            argument (float): phase angle (in degrees) to manipulate the
              displacement by.
        """
        u = np.ravel(modulation)
        index_max_elem = np.argmax(abs(u))
        max_elem = u[index_max_elem]
        phase_for_zero = max_elem / abs(max_elem)
        phase_factor = np.exp(1j * np.pi * argument / 180) / phase_for_zero

        return phase_factor

    def _map_eigvec_supercell(self, supercell, eigvec):
        """Maps the specified eigenvector to the target supercell.

        Args:
            supercell (phonopy.structure.cells.Supercell): the target supercell
              that the atoms will be displaced in.
            eigvec (numpy.ndarray): target eigenvector from the *primitive* cell.
        """
        s2u_map = supercell.get_supercell_to_unitcell_map()
        u2u_map = supercell.get_unitcell_to_unitcell_map()
        s2uu_map = [u2u_map[x] for x in s2u_map]
        result = []
        
        for i in range(supercell.get_number_of_atoms()):
            eig_index = s2uu_map[i] * 3
            ej = eigvec[eig_index:eig_index + 3]
            result.append(ej)

        return np.array(result)
    
    def _get_displacements(self, supercell, eigvec, q, amplitude, argument):
        """Returns the vector of displacements for a single eigenvector.

        .. warning:: The supercell *must* be comensurate with the selected
          q-vector.

        Args:
            supercell (phonopy.structure.cells.Supercell): the target supercell
              that the atoms will be displaced in.
            eigvec (numpy.ndarray): target eigenvector from the *primitive* cell.
        """
        m = supercell.get_masses()
        spos = supercell.get_scaled_positions()
        dim = supercell.get_supercell_matrix()
        coefs = np.exp(2j * np.pi * np.dot(np.dot(spos, dim.T), q)) / np.sqrt(m)
        meigvec = self._map_eigvec_supercell(supercell, eigvec)
        u = (meigvec.T*coefs).T
        u = np.array(u) / np.sqrt(len(m))
        phase_factor = self._get_phase_factor(u, argument)
        u *= phase_factor * amplitude
            
        return u

    def iterate(self, method="hessian", supercell=None, q=None, nrattle=0):
        """Returns a list of possible configurations, one for each eigenmode in the
        system, that has `supercell` is compatible with the specified `q`
        vector.

        Args:
            method (str): desired method for computing the eigenvalues and
              eigenvectors. Possible choices are :meth:`hessian` or
              :meth:`dynmat`.
            supercell (numpy.ndarray): supercell matrix to use in generating the
              configs.
            q (numpy.ndarray): q-vector that the resulting supercell should be
              compatible with.
            nrattle (int): number of additional, empty configs to include via
              :meth:`quippy.Atoms.rattle`.
        """
        if method == "hessian":
            dmd = self.hessian()
        elif method == "vasp_hessian":
            dmd = self.vasp_hessian()
        else:
            dmd = self.dynmat(supercell, q)
            
        hname = "hessian1"
        seed = quippy.Atoms()
        seed.copy_from(dmd["template"])
        result = quippy.AtomsList()
        result.append(dmd["template"])

        #Delete the force, energy and virial information from the copy, so that
        #we don't duplicate it all over.
        del seed.params["dft_energy"]
        del seed.params["dft_virial"]
        del seed.properties["dft_force"]
        
        for l, v in zip(dmd["eigvals"], dmd["eigvecs"]):
            atc = seed.copy()
            
            #Add this eigenvector to its own configuration.
            atc.add_property(hname, 0.0, n_cols=3)
            H = np.reshape(v.real, (atc.n, 3)).T
            setattr(atc, hname, H)
            
            #Same thing for the eigenvalue. Then save it to the group folder
            #structure.                
            atc.params.set_value(hname, l)
            atc.params.set_value("n_hessian", 1)
            #atc.add_property("force", 0.0, n_cols=3)
            result.append(atc)

        for i in range(nrattle):
            atRand = seed.copy()
            quippy.randomise(atRand.pos, 0.2)
            result.append(atRand)

        # atz = seed.copy()
        # atz.add_property("dft_force", 0.0, n_cols=3)
        #atz.params.set_value("dft_energy", 
            
        return result

    def vasp_hessian(self):
        """Extracts the hessian from `vasprun.xml`.
        """
        import xml.etree.ElementTree as ET
        tree = ET.parse('vasprun.xml')
        root = tree.getroot()
        calc = root.getchildren()[-2]
        dynm = calc.find("dynmat")
        hessian = dynm.getchildren()[0]
        H = []
        for v in hessian.getchildren():
            H.append(map(float, v.text.split()))

        H = -np.array(H)
        eigvals, eigvecs = np.linalg.eigh(H)
        Na = self.primitive.n*int(np.linalg.det(self.supercell))

        result = {
            "template": self.parent,
            "eigvals": [],
            "eigvecs": [],
            "hessian": H
        }
        
        for i, l in enumerate(eigvals):
            if np.abs(l) > 1e-3:
                result["eigvals"].append(l)
                result["eigvecs"].append(eigvecs[:,i].reshape(Na, 3))

        return result
        
    def hessian(self):
        """Returns the non-zero eigenvalues and their corresponding eigenvectors.
        """
        supercell = self.phonopy.get_supercell()
        Na = supercell.get_number_of_atoms()
        Nf = Na*3
        H = self.phonopy._dynamical_matrix._force_constants.reshape((Nf,Nf))
        eigvals, eigvecs = np.linalg.eigh(H)

        result = {
            "template": phonopy_to_ase(supercell),
            "eigvals": [],
            "eigvecs": []
        }
        
        for i, l in enumerate(eigvals):
            if np.abs(l) > 0.1:
                result["eigvals"].append(l)
                result["eigvecs"].append(eigvecs[:,i].reshape(Na, 3))

        return result
    
    def dynmat(self, supercell, q=None, cutoff=0.1):
        """Returns the non-zero eigenvalues and their corresponding eigenvectors
        for the specified supercell.

        Args:
            supercell (numpy.ndarray): supercell matrix to use in generating the
              configs.
            q (numpy.ndarray): q-vector that the resulting supercell should be
              compatible with.
            cutoff (float): minimum value an eigenvalue should have before it is
              included in the set.
        """
        #We need to determine the supercell matrix that is compatible with the
        #given `q` and has `N` atoms.
        scell = get_supercell(self.primitive, supercell)
        eigvals, eigvecs = self.diagonalize(q)

        result = {
            "template": phonopy_to_ase(scell),
            "eigvals": [],
            "eigvecs": []
        }
        for i, l in enumerate(eigvals):
            if np.abs(l) > 0.1:
                meigvec = self._map_eigvec_supercell(scell, eigvecs[:,i])
                result["eigvals"].append(l)
                result["eigvecs"].append(meigvec)

        return result
    
    def _get_dynmatrix(self):
        """Extracts the force constants from `FORCE_SETS` and constructs the
        dynamical matrix for the calculation.
        """
        with chdir(self.folder):
            force_sets = file_IO.parse_FORCE_SETS()
            self.phonopy.set_displacement_dataset(force_sets)
            self.phonopy.produce_force_constants(
                calculate_full_force_constants=True,
                computation_algorithm="svd")
            self.phonopy._set_dynamical_matrix()    
            
class Hessian(Group):
    """Represents a collection of Hessian eigenvectors and eigenvalues
    associated with dynamical matrices calculated for local minima in the phase
    space.

    Args:
        dynmat (dict): dynamical matrix with calculated eigenvectors and
          eigenvals.
    """
    def __init__(self, name="hessian", **dbargs):
        self.name = name
        self.seeded = True
        dbargs["prefix"] = 'H'
        super(Hessian, self).__init__(trainable=True, **dbargs)

    @property
    def rset(self):
        """Returns a list of folder names, each of which has an `atoms.json`
        file with configuration and hessian eigenvector/eigenvalue data
        attached.
        """
        if len(self.sequence) == 0:
            #Just return this group's hessian configurations; it is at the
            #bottom of the stack.
            return self.configs.values()
        else:
            result = []
            for h in self.sequence.values():
                result.extend(h.rset)
            return result

    def ready(self):
        """Determines if all the AFLOW configurations specified in the query
        have been downloaded and saved yet.
        """
        if isinstance(self.seed, dict) and "atoms" in self.seed:
            return len(self.configs) > 0
        else:
            return all(h.ready() for h in self.sequence.values())

    def setup(self, rerun=False):
        """Diagonalizes the dynamical matrix for each seed configuration and creates a
        set of new configurations that includes each eigenvector and eigenvalue.
        """
        super(Hessian, self).setup(self._setup_configs, rerun)
        
    def _setup_configs(self, rerun):
        """Sets up the folders for each of the eigenvalue/eigenvector sub-configs
        calculated from the Hessian matrix. Extracts the hessian matrix and
        assigns its eigenvalues and eigenvectors to relevant properties on
        copies of the seed configuration atoms objects.

        Args:
            rerun (bool): when True, re-diagonalize all the seed
              configurations.
        """
        dmatrix = self.seed["dynmat"]
        eigvecs, eigvals = dmatrix["eigvecs"], dmatrix["eigvals"]
        
        hname = "hessian"
        for i, e in enumerate(eigvals):
            #Ignore zero eigenvalues, not useful.
            if not np.isclose(e, 0.0):
                atc = self.seed["atoms"].copy()
                #Add this eigenvector to its own configuration.
                atc.add_property(hname, 0.0, n_cols=3)
                H = np.reshape(eigvecs[:,i], (atc.n, 3)).T
                setattr(atc, hname, H)
                #Same thing for the eigenvalue. Then save it to the group folder
                #structure.                
                atc.params.set_value(hname, e)
                atc.params.set_value("n_hessian", 1)
                #The seed should already have force information attached from
                #the dynamical matrix step. atc.add_property("force", 0.0, n_cols=3)
                self.create(atc)
