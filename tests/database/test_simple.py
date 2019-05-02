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

@pytest.fixture()
def Pd_not_extractable(tmpdir):
    target = relpath("./tests/Pd/matdb_not_extractable.yml")
    dbdir = str(tmpdir.join("manual_not_extractale_db"))
    mkdir(dbdir)
    copyonce(target, path.join(dbdir, "matdb.yml"))

    seed_root = path.join(dbdir, "seed")
    if not path.isdir(seed_root):
        mkdir(seed_root)

    cfg_target = path.join(seed_root, "Pd")
    cfg_source = path.join(_get_reporoot(), "tests", "database", "files", "Pd", "POSCAR1")
    copyonce(cfg_source, cfg_target)

    target = path.join(dbdir,"matdb")
    cntrl = Controller(target, dbdir)

    # these 3 lines are just added to cover dbargs["calculator"] not None in __init__ 
    db = Database("phonon", dbdir, cntrl, [{"type":"simple.Manual"}], {}, 0)
    dbargs = {"root": dbdir, "parent": db, "calculator": cntrl.calculator}
    mdb = Manual(extractable=False, **dbargs)

    return cntrl

@pytest.fixture()
def Pd_no_seeds(tmpdir):
    target = relpath("./tests/Pd/matdb_no_seeds.yml")
    dbdir = str(tmpdir.join("manual_no_seeds_db"))
    mkdir(dbdir)
    copyonce(target, path.join(dbdir, "matdb.yml"))

    target = path.join(dbdir,"matdb")
    cntrl = Controller(target, dbdir)

    return cntrl

def test_not_extractable(Pd_not_extractable):
    """ test not extractable 
    """
    mPd = Pd_not_extractable
    mPd.setup()

    mdb = mPd.collections['phonon'].steps['manual']
    assert mdb is not None
    assert mdb.nconfigs == 1 
    assert mdb.sub_dict() == {'extractable': False, 'name': 'manual'}
    assert not mdb.extractable
    assert not mdb._trainable

    assert mdb.is_setup()
    assert len(mdb.sequence) == 1
    assert len(mdb.sequence['Pd'].configs) == 1

    assert mdb.ready()
    assert mdb.can_extract()

    folders = {
        "__files__": ["phonon_S1_uuid.txt"],
        "Pd": {
            "__files__": ["phonon_S1_uuid.txt", "compute.pkl"],
            "S1.1": {
                "__files__": ["INCAR", "PRECALC", "POSCAR", "POTCAR", "atoms.h5",
                              "uuid.txt", "KPOINTS", "ase-sort.dat"]
            }
        }
    }

    dbfolder = mdb.root
    compare_tree(dbfolder,folders)    

    mdb.tarball()
    assert path.isfile(path.join(dbfolder, "Pd", "output.tar.gz"))

def test_all(Pd):
    """Tetsts setup/extract/ready of the simple.Manual database.
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

    assert not mdb.ready()

    mdb.extract()
    assert len(mdb.sequence) == 3
    assert len(mdb.sequence['Pd1'].config_atoms) == 1
    assert len(mdb.sequence['Pd2'].config_atoms) == 1
    assert len(mdb.sequence['Pd3'].config_atoms) == 1
    assert len(mdb.sequence['Pd1'].configs) == 1
    assert len(mdb.sequence['Pd2'].configs) == 1
    assert len(mdb.sequence['Pd3'].configs) == 1

    assert len(mdb.fitting_configs) == 3
    assert len(mdb.rset) == 3
    assert not mdb.is_executing()

    assert mdb.ready()

    # setup again(with rerun=True) on an already ready database
    Pd.setup(rerun=True)
    assert mdb.is_setup()
    assert mdb.ready()
    assert len(mdb.sequence) == 3
    assert len(mdb.sequence['Pd1'].config_atoms) == 1
    assert len(mdb.fitting_configs) == 3
    assert len(mdb.rset) == 3

def test_corner_cases_in_Group(Pd):
    """Tetsts some corner cases in class Group.
    """

    Pd.setup()
    mdb = Pd.collections['phonon'].steps['manual']
    dbfolder = mdb.root

    # create a file named `S1.4`
    open(path.join(dbfolder,"Pd1","S1.4"), 'a').close()
    # create a folder named `S1.A`
    mkdir(path.join(dbfolder,"Pd1","S1.A"))

    # the file `S1.1` and folder `S1.A` should be ignored
    assert mdb.is_setup()
    assert len(mdb.sequence) == 3
    assert len(mdb.sequence['Pd1'].configs) == 1

    # create a empty folder named `S1.5` 
    mkdir(path.join(dbfolder,"Pd1","S1.5"))
    # shold cause a warning message and is_setup() should return False
    assert not mdb.is_setup()

    for atom in mdb.iconfigs:
        assert atom is not None
        
    assert mdb.key == 'phonon.manual'
    assert 'rset.h5' in mdb.rset_file 

    mdb.jobfile()
    assert path.isfile(path.join(dbfolder, "Pd1", "jobfile.sh"))
    assert mdb.execute(dryrun=True)

def test_corner_cases_in_Group_no_seeds(Pd_no_seeds):
    """Tetsts some corner cases in class Group.
    """
    Pd_no_seeds.setup()
    mdb = Pd_no_seeds.collections['phonon'].steps['manual']
    dbfolder = mdb.root

    assert not mdb.is_setup()
    assert mdb.seeded
    assert mdb._seed is None
