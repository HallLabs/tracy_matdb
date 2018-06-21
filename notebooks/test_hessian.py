"""Tests the logic to construct Hessian training databases.
"""
import pytest
import numpy as np
import quippy
from os import path

from matdb.database.hessian import HessianSupercell
from matdb.utility import relpath

def test_eigvecs_displace():
    """Tests that our eigenvector extraction from the dynamical matrix are
    correct.
    """
    #First, we setup a supercell calculator based on NbP seed.
    target = relpath("./tests/data/hessian")
    smatrix = np.array([-2,  0, -1, -1, -1,  0,  0,  2,  2]).reshape(3,3)
    s0 = quippy.Atoms(path.join(target, "POSCAR"), format="vasp")
    hs = HessianSupercell(s0, smatrix, target)

    #Next, calculate the eigenvectors. We are interested in band 3 for this unit
    #test.
    q = np.array([0.3500000, 0.1500000, 0.0500000])
    #We put 4 here because that is what goes in the phonopy input file. But
    #python is off by 1. We pick the amplitude and phase randomly.
    band = 4 - 1 
    ampl = 16.0492756 
    phase = 53.2673576
    evals, evecs = hs.diagonalize(q)

    #Test the eigenvector of the *primitive* dynamical matrix at that q-point.
    modev = np.array([ 0.38721282+0.j, 0.25880307-0.02083032j,
                       0.20992689-0.07299886j, 0.0205681 +0.38666616j,
                       -0.00705372+0.25954417j, -0.06174484+0.2135081j,
                       -0.43758703+0.07162384j, -0.15599019-0.00706084j,
                       0.00309682-0.11068152j, 0.04827883-0.4407738j,
                       -0.01533681-0.15539491j, -0.11036077+0.00897167j])
    band_eig = evecs[:,band]
    assert np.allclose(band_eig, modev)

    #Next, test the expansion of the eigenvector into the *supercell*. We have a
    #24-atom supercell, so we need an eigenvector of length 24. We need to
    #construct an explicit supercell object to test the internals of the class.
    from phonopy.structure.cells import get_supercell
    scell = get_supercell(hs.primitive, smatrix)
    model = np.load(path.join(target, "super-eigv.npy"))
    assert np.allclose(hs._map_eigvec_supercell(scell, band_eig), model)

    model_disp = np.load(path.join(target, "super_disp.npy"))
    calc_disp = hs._get_displacements(scell, band_eig, q, ampl, phase)
    assert np.allclose(calc_disp, model_disp)
