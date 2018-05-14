"""Tests the  group interface.
"""
import pytest
# from matdb.database.distortions import Distortions
from os import mkdir, path, symlink, remove
# import numpy as np


@pytest.fixture()
def AlMg(tmpdir):
    """Test the functionality of the distortions group for a small seed
    configuration where the number of possible random combinations is limited
    and repition is more likely.

    Attributes:
        atoms_seed(matdb.atoms): The seed atoms object configuration for db
            generation.
    """
    from matdb.utility import relpath
    from matdb.database import Controller
    from shutil import copy

    target = relpath("./tests/AlMg/matdb")
    dbdir = str(tmpdir.join("almg_db"))
    mkdir(dbdir)

    # We need to copy the POSCAR from the testing directory to tempdir.
    POSCAR = relpath("./tests/AlMg/Al6Mg4")
    mkdir(path.join(dbdir, "seed"))
    copy(POSCAR, path.join(dbdir, "seed", "Al6Mg4"))
    if path.isfile("matdb.yml"):
        remove("matdb.yml")
    symlink("{}.yml".format(target), "matdb.yml")

    result = Controller(target, dbdir)
    return result


def test_AlMg_setup(AlMg):
    """Test the setup of the substitutions database.
    """
    assert not AlMg.collections[
        'distortions'].steps['Dist'].is_setup()

    AlMg.setup()

    dbs = "dist/distortions/Al6Mg4"
    '''
    folders = {
        "__files__": ["compute.pkl", "suids.pkl", "jobfile.sh", "index.json"],
        "S.1": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.2": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.3": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.4": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.5": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.6": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.7": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.8": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.9": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.10": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.11": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.12": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.13": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.14": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "S.15": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        }
    }

    from matdb.utility import compare_tree
    dbfolder = path.join(AgCu.root, dbs)
    compare_tree(dbfolder, folders)

    assert AgCu.collections['substitution'][
        'substitution'].steps['sub'].is_setup()

    # test the suid and index creation for the entire database.
    assert path.isfile(path.join(AgCu.root, "sub/substitution/suids.pkl"))
    assert path.isfile(path.join(AgCu.root, "sub/substitution/index.json"))

    sub = AgCu.collections['substitution']['substitution'].steps['sub']
    assert len(sub.index) == 15
    assert len(sub.suids) == 15

    # assert not enum.ready()
    assert not sub.ready()

    # We need to fake some VASP output so that we can cleanup the
    # database and get the rset

    src = relpath(
        "./tests/data/Pd/complete/OUTCAR__DynMatrix_phonon_Pd_dim-2.00")

    dbfolder = path.join(AgCu.root, dbs)
    for j in range(1, 16):
        dest = path.join(dbfolder, "S.{}".format(j), "OUTCAR")
        symlink(src, dest)

    dbfolder = path.join(AgCu.root, dbs)
    for j in range(1, 16):
        src = path.join(dbfolder, "S.{}".format(j), "POSCAR")
        dest = path.join(dbfolder, "S.{}".format(j), "CONTCAR")
        symlink(src, dest)

    sub.cleanup()
    assert len(sub.atoms_paths()) == 15
    assert len(sub.rset()) == 15
    '''
