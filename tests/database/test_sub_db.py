"""Tests the substitution group interface.
"""
import pytest
from matdb.database.substitution import _round_to_even
from matdb.utility import relpath
from os import mkdir, path, symlink, remove
import numpy as np


@pytest.fixture()
def AgCu(tmpdir):
    """Test the functionality of the substitution group for a small seed
    configuration where the number of possible random combinations is limited
    and repition is more likely.

    Attributes:
        atoms_seed(matdb.atoms): The seed atoms objech configuration for db
            generation.
    """
    from matdb.utility import relpath
    from matdb.database import Controller
    from shutil import copy

    target = relpath("./tests/AgCu/matdb")
    dbdir = str(tmpdir.join("agcu_db"))
    mkdir(dbdir)

    # We need to copy the POSCAR from the testing directory to tempdir.
    POSCAR = relpath("./tests/AgCu/Ag1Cu5")
    mkdir(path.join(dbdir, "seed"))
    copy(POSCAR, path.join(dbdir, "seed", "Ag1Cu5"))
    if path.isfile("matdb.yml"):
        remove("matdb.yml")
    symlink("{}.yml".format(target), "matdb.yml")

    result = Controller(target, dbdir)
    return result


def test_AgCu_setup(AgCu):
    """Test the setup of the substitutions database.
    """
    assert not AgCu.collections[
        'substitution'].steps['Substitution'].is_setup()

    AgCu.setup()

    dbs = "Substitution/substitution/Ag1Cu5"

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

    assert AgCu.collections[
        'substitution'].steps['Substitution'].is_setup()

    # test the suid and index creation for the entire database.
    assert path.isfile(path.join(AgCu.root,
                                 "Substitution/substitution/suids.pkl"))
    assert path.isfile(path.join(AgCu.root,
                                 "Substitution/substitution/index.json"))

    sub = AgCu.collections[
        'substitution'].steps['Substitution']
    assert len(sub.index) == 15
    assert len(sub.suids) == 15
    assert not sub.ready()

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

    assert len(sub.atoms_paths()) == 0


def test_functions(AgCu):
    """Tests the specific functions of the substitutions database
    """
    # test the rounding function to mimic round behavior in python3.
    rounded = [0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3,
               3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7]
    vals = np.arange(0, 6.6, 0.3)
    for i in range(len(vals)):
        assert _round_to_even(vals[i]) == rounded[i]


def test_initial(AgCu):
    """Tests the substitutions specific functions of the database.
    """
    sub = AgCu.collections[
        'substitution'].steps['Substitution']

    sub.__init__(root=sub.root, ran_seed=12, seeds=["vasp:Ag1Cu5"],
                 nconfigs=15,
                 stoich=[[0.5, 0.5, 0.5], [0.7, 0.3, 0.3], [0.4, 0.6, 0.2]],
                 min_index=5, parent=sub.parent)

    assert sub.ran_seed == 12
    assert sub.stoich == [[0.5, 0.5, 0.5], [0.7, 0.3, 0.3], [0.4, 0.6, 0.2]]
    assert sub.min_index == 5


def test_set_stoichiometry(AgCu):
    sub = AgCu.collections['substitution'].steps['Substitution']
    new_pgrid = sub.pgrid
    new_pgrid.params.update({'stoich': [[0.2, 0.8, 0.2], [0.4, 0.6, 0.3],
                                        [0.5, 0.5, 0.5]],
                             'nconfigs': 20, 'min_index': 0})
    sub. __init__(root=sub.root, seeds=["vasp:Ag1Cu5"], parent=sub.parent,
                  pgrid=new_pgrid)
    sub.setup()

    stoich = [[[0.2, 0.3, 0.5], [0.5, 0.5, 0.4], [0.9, 0.1, 0.1]],
              [[0.5, 4.2, 0.5], [0.5, 0.4, 0.4], [0.9, 0.1, 0.1]],
              [[0.5, 0.5, 0.7], [0.5, 0.5, 0.4], [0.9, 0.1, 0.1]],
              [[0.5, 0.5, 1.5], [0.5, 0.5, 0.4], [0.9, 0.1, 0.1]],
              [[0.5, 0.5, 0.5], [0.5, 0.5, 1.4], [0.9, 0.1, 0.1]],
              [[0.8, 0.3, 0.5], [0.5, 0.5, 0.4], [0.9, 0.1, 0.1]],
              [[0.1, 0.9, 0.8], [0.5, 0.5, 0.4], [0.9, 0.1, 0.1]],
              [[0.8, 0.2, 0.1], [0.5, 0.5, 0.1], [0.9, 0.1, 0.8]],
              [[0.5, 0.5, 0.1], [0.6, 0.4, 0.1], [0.9, 0.1, 0.8]]]
    for i in range(8):
        with pytest.raises(ValueError):
            new_pgrid.params.update({'stoich': stoich[i], 'nconfigs': 20})
            sub.__init__(root=sub.root, seeds=["vasp:Ag1Cu5"],
                         parent=sub.parent, pgrid=new_pgrid)
            sub.setup()

