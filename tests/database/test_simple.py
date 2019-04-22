"""Tests the simple group interface.
"""
import pytest
from os import mkdir, path, symlink, remove
import numpy as np
import six

from matdb.database.simple import Manual
from matdb.utility import relpath, compare_tree, copyonce
from matdb.utility import _get_reporoot
from matdb.fitting.controller import TController
from matdb.database import Database, Controller
from matdb.atoms import Atoms

@pytest.fixture()
def Pd(tmpdir):
    target = relpath("./tests/Pd/matdb.yml")
    dbdir = str(tmpdir.join("manual_db"))
    mkdir(dbdir)
    copyonce(target, path.join(dbdir, "matdb.yml"))

    seed_root = path.join(dbdir, "seed")
    if not path.isdir(seed_root):
        mkdir(seed_root)

    for i in range(1, 4):
        cfg_target = path.join(seed_root, "Pd{0}".format(i))
        cfg_source = path.join(_get_reporoot(), "tests", "database", "files", "Pd", "POSCAR{0}".format(i))
        copyonce(cfg_source, cfg_target)

    target = path.join(dbdir,"matdb")
    cntrl = Controller(target, dbdir)

    return cntrl

def test_init(Pd):
    dbdir = Pd.root
    db = Database("manual", dbdir, Pd, [{"type":"simple.Manual"}], {}, 0)
    dbargs = {"root": dbdir, "parent": db, "calculator": Pd.calculator}
    mdb = Manual(extractable=False, **dbargs)
    
    assert mdb is not None
    assert mdb.nconfigs == 1 
    
def test_setup(Pd):
    """Tetsts the setup of the simple.Manual database.
    """

    Pd.setup()

    mdb = Pd.collections['phonon'].steps['manual']

    assert Pd is not None
    assert mdb is not None
    assert mdb.is_setup()
    assert len(mdb.sequence) == 3
    assert len(mdb.sequence['Pd1'].configs) == 1

    folders = {
        "__files__": ["phonon_S1_uuid.txt"],
        "Pd1": {
            "__files__": ["phonon_S1_uuid.txt", "jobfile.sh", "compute.pkl"],
            "S1.1": {
                "__files__": ["INCAR", "PRECALC", "POSCAR", "POTCAR", "pre_comp_atoms.h5",
                              "uuid.txt", "KPOINTS", "ase-sort.dat"]
            }
        },
        "Pd2": {
            "__files__": ["phonon_S1_uuid.txt", "jobfile.sh", "compute.pkl"],
            "S1.1": {
                "__files__": ["INCAR", "PRECALC", "POSCAR", "POTCAR", "pre_comp_atoms.h5",
                              "uuid.txt", "KPOINTS", "ase-sort.dat"]
            }
        },
        "Pd3": {
            "__files__": ["phonon_S1_uuid.txt", "jobfile.sh", "compute.pkl"],
            "S1.1": {
                "__files__": ["INCAR", "PRECALC", "POSCAR", "POTCAR", "pre_comp_atoms.h5",
                              "uuid.txt", "KPOINTS", "ase-sort.dat"]
            }
        }
    }

    dbfolder = mdb.root
    compare_tree(dbfolder,folders)    

    assert mdb.is_setup()
    assert not mdb.ready()

    # We need to fake some VASP output so that we can cleanup the
    # database and get the rset
    src = relpath("./tests/data/Pd/complete/OUTCAR__DynMatrix_phonon_Pd_dim-2.00")
    dbfolder = mdb.root
    for j in range(1,4):
        dest = path.join(dbfolder,"Pd{}".format(j),"S1.1", "OUTCAR")
        symlink(src,dest)
            
    dbfolder = mdb.root
    for j in range(1,4):
        src = path.join(dbfolder,"Pd{}".format(j),"S1.1", "POSCAR")
        dest = path.join(dbfolder,"Pd{}".format(j),"S1.1", "CONTCAR")
        symlink(src,dest)

    #remove(path.join(mdb.root, "Pd1", "S1.1", "pre_comp_atoms.h5"))
    mdb.extract()
    assert len(mdb.config_atoms) == 3
    assert len(mdb.configs) == 3
    assert len(mdb.last_config_atoms) == 3
    assert len(mdb.rset) == 3
    assert not mdb.is_executing()

    remove(path.join(mdb.root, "Pd1", "S1.1", "OUTCAR"))
    dest = path.join(mdb.root, "Pd1", "S1.1", "OUTCAR")
    src = relpath("./tests/data/Pd/basic_fail/S.4/OUTCAR")
    symlink(src,dest)
    assert mdb.is_executing()

    #mdb.rset
    #mdb.fitting_configs

