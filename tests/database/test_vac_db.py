"""Tests the vacancy group interface.
"""
import pytest
from matdb.utility import relpath
from os import mkdir, path, symlink, remove

@pytest.fixture()
def AlTi(tmpdir):
    """Test the functionality of the vacancy group for a small seed
    configuration where the number of possible random combinations is limited
    and repition is more likely.

    Attributes:
        atoms_seed(matdb.atoms): The seed atoms objech configuration for db
            generation.
    """
    from matdb.utility import relpath
    from matdb.database import Controller
    from shutil import copy

    target = relpath("./tests/AlTi/matdb")
    dbdir = str(tmpdir.join("alti_db"))
    mkdir(dbdir)

    # We need to copy the POSCAR from the testing directory to tempdir.
    POSCAR = relpath("./tests/AlTi/Al14Ti6")
    mkdir(path.join(dbdir, "seed"))
    copy(POSCAR, path.join(dbdir, "seed", "Al14Ti6"))
    if path.isfile("matdb.yml"):
        remove("matdb.yml")
    symlink("{}.yml".format(target), "matdb.yml")

    result = Controller(target, dbdir)
    return result


def test_AlTi_setup(AlTi):
    """Test the setup of the vacancy database.
    """
    assert not AlTi.collections['vacancy'].steps[
        'Vacancy'].is_setup()

    AlTi.setup()

    dbs = "Vacancy/vacancy.Vacancy/Al14Ti6"

    folders = {
        "__files__": ["compute.pkl", "vuids.pkl", "jobfile.sh", "index.json"],
        "V.1": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "V.2": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "V.3": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "V.4": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "V.5": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "V.6": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        },
        "V.7": {
            "__files__": ["INCAR", "PRECALC", "POSCAR"]
        }
    }

    from matdb.utility import compare_tree
    dbfolder = path.join(AlTi.root, dbs)
    compare_tree(dbfolder, folders)

    assert AlTi.collections[
        'vacancy'].steps['Vacancy'].is_setup()

    # test the vuid and index creation for the entire database.
    assert path.isfile(path.join(AlTi.root,
                                 "Vacancy/vacancy.Vacancy/vuids.pkl"))
    assert path.isfile(path.join(AlTi.root,
                                 "Vacancy/vacancy.Vacancy/index.json"))

    vac = AlTi.collections[
        'vacancy'].steps['Vacancy']
    assert len(vac.index) == 50
    assert len(vac.vuids) == 50
    assert vac.ready()

    # We need to fake some VASP output so that we can cleanup the
    # database and get the rset

    src = relpath(
        "./tests/data/Pd/complete/OUTCAR__DynMatrix_phonon_Pd_dim-27.00")

    dbfolder = path.join(AlTi.root, dbs)
    for j in range(1, 51):
        dest = path.join(dbfolder, "V.{}".format(j), "OUTCAR")
        symlink(src, dest)

    dbfolder = path.join(AlTi.root, dbs)
    for j in range(1, 51):
        src = path.join(dbfolder, "V.{}".format(j), "POSCAR")
        dest = path.join(dbfolder, "V.{}".format(j), "CONTCAR")
        symlink(src, dest)

    assert len(vac.atoms_paths()) == 50
    assert len(vac.rset()) == 50
