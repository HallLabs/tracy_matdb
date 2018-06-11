"""Tests fitting using legacy databases.
"""
import pytest

@pytest.fixture()
def mtpDB(tmpdir):
    from matdb.utility import relpath, reporoot
    from matdb.database import Controller
    from os import mkdir

    target = relpath("./tests/AgPd/matdb")
    dbdir = str(tmpdir.join("mtp"))
    mkdir(dbdir)

    return Controller(target, tmpdir=dbdir)

def test_init(mtpDB):
    """Tests the initialization of the mtp potenital fitter.
    """
    mtp = mtpDB.trainers.fits['AgPt'].sequences['AgPt'].steps['mtp']

    assert not mtp.ready()
