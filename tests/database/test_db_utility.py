# -*- coding: utf-8 -*-
"""Tests the utility functions for the databases.
"""
import pytest
import numpy as np
from os import mkdir, path, symlink, remove, rename, listdir
from glob import glob

from matdb.database.active import Active
from matdb.utility import relpath, copyonce
from matdb.utility import _get_reporoot
from matdb.fitting.controller import TController
from matdb.database import Database, Controller
from matdb.database.utility import parse_path, split

@pytest.fixture()
def Act(tmpdir):

    target = relpath("./tests/AgPd/matdb.yml")
    dbdir = str(tmpdir.join("active_db"))
    mkdir(dbdir)
    copyonce(target, path.join(dbdir, "matdb.yml"))
    target = path.join(dbdir,"matdb")

    seed_root = path.join(dbdir, "seed")
    if not path.isdir(seed_root):
        mkdir(seed_root)

    for i in range(1, 4):
        cfg_target = path.join(seed_root, "Pd{0}".format(i))
        cfg_source = path.join(_get_reporoot(), "tests", "database", "files", "Pd", "POSCAR{0}".format(i))
        copyonce(cfg_source, cfg_target)

    cntrl = Controller(target, dbdir)
    db = Database("active", dbdir, cntrl, [{"type":"active.Active"}], {}, 0)
    tcntrl = TController(db=db, root=dbdir, fits={})
    dbargs = {"root": dbdir, "parent": db,
              "calculator": tcntrl.db.calculator}
    result = Active(**dbargs)

    return result

def test_split(Act):
    """Tests the split function.
    """
    from matdb.atoms import Atoms, AtomsList

    at1 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43],info={"rand":10})
    at2 = Atoms("S6",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75]],
                 cell=[6.43,5.43,4.43],info={"rand":10})
    at3 = Atoms("CNi",positions=[[0,0,0],[0.5,0.5,0.5]], info={"rand":8})
    at4 = Atoms("CoV",positions=[[0,0,0],[0.25,0.5,0.25]], info={"rand":8})

    at1.add_param("vasp_energy", 25361.504084423999)
    at2.add_param("vasp_energy", 25362.504084423999)
    at3.add_param("vasp_energy", 25363.504084423999)
    at4.add_param("vasp_energy", 25364.504084423999)

    al = [at4,at2,at1,at3]

    splits = {'A': 0.4, 'B': 0.2}
    splitroot = path.join(Act.root, "split")
    mkdir(splitroot)
    file_targets = {"train": Act.database.train_file, "holdout": Act.database.holdout_file,
                    "super": Act.database.super_file}
    # the split directory should be empty
    assert len(listdir(splitroot)) == 0

    # split with an empty splits
    split(al, {}, file_targets, splitroot, ran_seed=1, recalc=1)

    # the split directory should still be empty
    assert len(listdir(splitroot)) == 0

    # split with an empty atom list
    split([], splits, file_targets, splitroot, ran_seed=1, recalc=1)

    # the split directory should have 2 entries
    assert len(listdir(splitroot)) == 2

    # split 
    split(al, splits, file_targets, splitroot, ran_seed=1)
    # now we should have all the splitted files
    assert path.exists(path.join(splitroot, Act.database.train_file('A')))
    assert path.exists(path.join(splitroot, Act.database.holdout_file('A')))
    assert path.exists(path.join(splitroot, Act.database.super_file('A')))
    assert path.exists(path.join(splitroot, Act.database.train_file('B')))
    assert path.exists(path.join(splitroot, Act.database.holdout_file('B')))
    assert path.exists(path.join(splitroot, Act.database.super_file('B')))
    # the split directory should have 2 entries
    assert len(listdir(splitroot)) == 2

    # split again with recalc
    split(al, splits, file_targets, splitroot, ran_seed=1, recalc=1)
    assert glob(path.join(Act.database.root, "splits", "active", "A_*-train.h5")) 
    assert glob(path.join(Act.database.root, "splits", "active", "A_*-holdout.h5")) 
    assert glob(path.join(Act.database.root, "splits", "active", "A_*-super.h5")) 
    assert glob(path.join(Act.database.root, "splits", "active", "B_*-train.h5")) 
    assert glob(path.join(Act.database.root, "splits", "active", "B_*-holdout.h5")) 
    assert glob(path.join(Act.database.root, "splits", "active", "B_*-super.h5")) 
    # the split directory should have 2 entries
    assert len(listdir(splitroot)) == 2

def test_parse_path(Act):
    """Tests the parse_path function.
    """
    parsed=parse_path(Act.database.root,["Pd*"])  
    assert parsed is not None 
    assert len(parsed) == 3

    parsed=parse_path(Act.database.root,["/Pd/Act/*"])  
    assert parsed is not None 
    assert len(parsed) >= 3

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
