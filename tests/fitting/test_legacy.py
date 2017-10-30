"""Tests fitting using legacy databases.
"""
import pytest

@pytest.fixture()
def legDB(tmpdir):
    from matdb.utility import relpath, reporoot
    from matdb.database.controller import Controller
    from os import mkdir

    target = relpath("./tests/legacy/matdb.yaml")
    dbdir = str(tmpdir.join("legacy"))
    mkdir(dbdir)

    return Controller(target, tmpdir=dbdir)

def test_jobfile(legDB):
    """Tests creation of jobfiles for legacy databases.
    """
    legDB.trainers.jobfiles()

    legdbs = legDB.trainers.find("2b-A.*")
    for tb in legdbs:
        assert path.isfile(path.join(tb.root, "jobfile.sh"))
        assert path.isfile(path.join(tb.root, "train.xyz"))

