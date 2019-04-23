# -*- coding: utf-8 -*-
"""Tests the utility functions for the databases.
"""
import pytest
import numpy as np

def test_dbconfig():
    """Tests the db config function.
    """
    from matdb.database.utility import dbconfig
    from os import path
    import json
    from matdb.utility import load_datetime

    dbfile = "./tests/database/files/p-50-2.h5"
    confpath = dbfile + ".json"
    assert path.isfile(confpath)
    
    with open(confpath) as f:
         config = json.load(f, object_pairs_hook=load_datetime)
         #assertIsNotNone(config)
         assert config!= None

def test_swap_column():
    """Tests the swap column function.
    """
    from matdb.database.utility import swap_column

    hnf = np.array([[0,1,0],[1,1,0],[1,0,1]])
    b = np.array([[1,2,0],[3,1,2],[0,0,1]])

    new_hnf, new_b = swap_column(hnf, b, 0)

    assert np.allclose(new_hnf, [[1,0,0],[1,1,0],[0,1,1]])
    assert np.allclose(new_b, [[2,1,0],[1,3,2],[0,0,1]])

def test_hnf():
    """Tests the conversion of integer matrices to hnf form.
    """
    from matdb.database.utility import hermite_normal_form

    n = [[-1, -1, -2], [2, 6, -1], [-7, 0, 5]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [66, 85, 111]]
    true_b = [[-8, -10, -13], [3, 4, 5], [2, 3, 4]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[-8, 5, -1], [3, -5, 5], [-5, -9, 9]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [52, 137, 208]]
    true_b = [[-5, -13, -20], [-9, -24, -37], [-6, -16, -25]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[6, -6, -5], [5, 9, 7], [2, -6, -6]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [68, 34, 96]]
    true_b = [[-2, -1, -3], [47, 24, 67], [-59, -30, -84]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[0, 7, -5], [10, 5, 8], [2, -6, 0]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [296, 416, 462]]
    true_b = [[52, 73, 81], [-32, -45, -50], [-45, -63, -70]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[10, -1, -5], [-1, 7, 2], [-2, -5, 0]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [2, 3, 0], [2, 1, 3]]
    true_b = [[14, 12, 11], [-6, -5, -5], [29, 25, 23]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[1, 0, 0], [3, 2, 0], [1, -1, -2]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [1, 2, 0], [0, 1, 2]]
    true_b = [[1, 0, 0],[-1, 1, 0], [1, -1, -1]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    with pytest.raises(ValueError):
        hermite_normal_form([[1,0,0],[1,0,0],[0,0,0]])


def test_make_primitive():
    """Tests the make_primitive routine.
    """

    from matdb.atoms import Atoms
    from matdb.database.utility import make_primitive
    from phenum.grouptheory import _is_equiv_lattice

    atm = Atoms(cell=[[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
                positions= [[0, 0, 0]], symbols="Pd")

    new_vecs, unique_pos, unique_types, hnf = make_primitive(atm)

    assert np.allclose(new_vecs, atm.cell)
    assert np.allclose(unique_pos, unique_pos)
    assert np.allclose(hnf, np.identity(3))
    assert unique_types[0] == "Pd"

    atm = Atoms(cell=[[0, 0, -1], [0, 1, 0], [1, -0.5, 0.5]],
                positions= [[0, 0, 0], [0.5, 0, -0.5], [0, 0.5, -0.5], [0.5, 0.5, 0]],
                numbers= [41, 41, 41, 41])

    new_vecs, unique_pos, unique_types, hnf = make_primitive(atm)
    new_lat = np.matmul(np.transpose(new_vecs), hnf)
    assert _is_equiv_lattice(new_lat, np.transpose(atm.cell), 1E-3)
    assert np.allclose(new_vecs, [[0.5, 0, -0.5], [0, 0.5, -0.5], [0.5, 0.5, 0]])
    assert np.allclose(unique_pos, [[0,0,0]])
    assert np.allclose(np.linalg.det(hnf), 4)
    assert unique_types[0] == "Nb"

    atm = Atoms(cell=[[0.5, 0.5, 0], [0, 0.5, 0.5], [1.5, -1, 1.5]],
                positions= [[0, 0, 0], [0.5, 0, 0.5], [1, 0, 1], [1.5, 0, 1.5]],
                symbols = "Pd3Ag")

    new_vecs, unique_pos, unique_types, hnf = make_primitive(atm)
    new_lat = np.matmul(new_vecs, hnf)
    assert _is_equiv_lattice(new_lat, atm.cell, 1E-3)
    assert np.allclose(new_vecs, atm.cell)
    assert np.allclose(unique_pos, atm.positions)

    atm = Atoms(cell=[[0, 0, -1], [0, 1, 0], [1, -0.5, 0.5]],
                positions= [[0, 0, 0], [0.5, 0, -0.5], [0, 0.5, -0.5], [0.5, 0.5, 0]],
                symbols= "AlPdAlPd")

    new_vecs, unique_pos, unique_types, hnf = make_primitive(atm)
    new_lat = np.matmul(np.transpose(new_vecs), hnf)
    assert _is_equiv_lattice(new_lat, np.transpose(atm.cell), 1E-3)
    assert np.allclose(new_vecs, [[0, 0.5, -0.5], [0, 0, -1], [1, -0.5, 0.5]])
    assert np.allclose(unique_pos, [[0,0,0],[0.5, 0, -0.5]])
    assert np.allclose(np.linalg.det(hnf), 2)
    assert "Al" in unique_types
    assert "Pd" in unique_types


    atm = Atoms(cell=[[1, 0, 0], [0.5, 0.8660254, 0], [0, 0, 1.6329932]],
                positions= [[0, 0, 0], [0.5, 0.2886751, 0.8164966]],
                symbols= "Al2")

    new_vecs, unique_pos, unique_types, hnf = make_primitive(atm)
    assert np.allclose(new_vecs, atm.cell)
    assert np.allclose(unique_pos, atm.positions)
    assert np.allclose(np.linalg.det(hnf), 1)

    
    with pytest.raises(ValueError):
        atm = Atoms(cell=[[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
                    positions= [[0, 0, 0]])

        stuff = make_primitive(atm)

def test_decompress():
    """Tests that the decompression algorithm works.
    """

    from matdb.atoms import Atoms
    from matdb.database.utility import make_primitive, decompress
    from phenum.grouptheory import _is_equiv_lattice

    atm = Atoms(cell=[[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
                positions= [[0, 0, 0]], symbols="Pd")

    new_vecs, unique_pos, unique_types, hnf = make_primitive(atm)
    unique_types = 1
    hnf_vec = [hnf[0][0], hnf[1][0], hnf[1][1], hnf[2][0], hnf[2][1], hnf[2][2]]
    hnf_int = "".join(["{0}0".format(i+1) for i in hnf_vec])
    lat_vecs, new_basis, new_types = decompress(new_vecs, unique_pos, unique_types, hnf_int)

    assert _is_equiv_lattice(atm.cell, lat_vecs, 1E-3)
    assert np.allclose(atm.positions, new_basis)
    assert new_types == [1]
    

    atm = Atoms(cell=[[0, 0, -1], [0, 1, 0], [1, -0.5, 0.5]],
                positions= [[0, 0, 0], [0.5, 0, -0.5], [0, 0.5, -0.5], [0.5, 0.5, 0]],
                numbers= [41, 41, 41, 41])

    new_vecs, unique_pos, unique_types, hnf = make_primitive(atm)
    unique_types = 1
    hnf_vec = [hnf[0][0], hnf[1][0], hnf[1][1], hnf[2][0], hnf[2][1], hnf[2][2]]
    hnf_int = "".join(["{0}0".format(i+1) for i in hnf_vec])
    lat_vecs, new_basis, new_types = decompress(new_vecs, unique_pos, unique_types, hnf_int)
    
    assert _is_equiv_lattice(np.transpose(atm.cell), np.transpose(lat_vecs), 1E-3)
    assert np.allclose([[0.0, 0.0, 0.0], [0.5, -0.5, -1.0], [0.0, 0.0, -1.0],
                        [-0.5, 0.5, -1.0]], new_basis)
    assert new_types == [1, 1, 1, 1]

def test_get_hnf_int():
    """Tests the conversion of ints to hnfs.
    """
    from matdb.database.utility import _get_hnf_from_int
    
    hnf_int = 120102004030110

    hnf = [[11,0,0],[0,19,0],[3,2,10]]

    out_hnf = _get_hnf_from_int(hnf_int)
    assert np.allclose(hnf, out_hnf)
