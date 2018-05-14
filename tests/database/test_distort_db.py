"""Tests the  group interface.
"""
import pytest
from os import mkdir, path, symlink, remove
from matdb.database import distortion
from matdb.utility import relpath

@pytest.fixture()
def AlMg(tmpdir):
    """Test the functionality of the distortion group for a small seed
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
    """Test the setup of the distortion database.
    """
    assert not AlMg.collections[
        'distortion'].steps['Distortion'].is_setup()

    AlMg.setup()

    dbs = "Distortion/distortion/Al6Mg4"

    folders = {
        "__files__": ["compute.pkl", "duids.pkl", "jobfile.sh", "index.json"],
        "D.1": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.2": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.3": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.4": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.5": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.6": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.7": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.8": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.9": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.10": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.11": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.12": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.13": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.14": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "D.15": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        }
    }

    from matdb.utility import compare_tree
    dbfolder = path.join(AlMg.root, dbs)
    compare_tree(dbfolder, folders)

    assert AlMg.collections['distortion'].steps[
        'Distortion'].is_setup()

    # test the duid and index creation for the entire database.
    assert path.isfile(path.join(AlMg.root,
                                 "Distortion/distortion/duids.pkl"))
    assert path.isfile(path.join(AlMg.root,
                                 "Distortion/distortion/index.json"))

    dist = AlMg.collections['distortion'].steps['Distortion']
    assert len(dist.duids) == 50

    # assert not enum.ready()
    assert not dist.ready()

    # We need to fake some VASP output so that we can cleanup the
    # database and get the rset

    src = relpath(
        "./tests/data/Pd/complete/OUTCAR__DynMatrix_phonon_Pd_dim-2.00")

    dbfolder = path.join(AlMg.root, dbs)
    for j in range(1, 51):
        dest = path.join(dbfolder, "D.{}".format(j), "OUTCAR")
        symlink(src, dest)

    dbfolder = path.join(AlMg.root, dbs)
    for j in range(1, 51):
        src = path.join(dbfolder, "D.{}".format(j), "POSCAR")
        dest = path.join(dbfolder, "D.{}".format(j), "CONTCAR")
        symlink(src, dest)
