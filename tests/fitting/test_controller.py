"""Tests methods of the controlling objects for the fitting.
"""
import pytest
from os import path

def test_jobfile(Pd_db):
    """Tests the writing of the jobfile for the fit, which includes constructing
    the command.
    """
    #Test the finder; we don't do this in a separate test because the tests are
    #independent and it takes long to setup the trainer grids.
    _2bA = Pd_db.trainers.find("2b-A*.*")
    model = ["2b-A-0", "2b-A-1", "2b-A-2", "2b-A-3", "2b-A-4", "2b-A-5"]
    _2b = [s.parent.name for s in _2bA]
    assert _2b == model

    _2bB = Pd_db.trainers.find("2b-B.*")
    
    Pd_db.trainers.jobfiles()
    for tb in _2bA:
        assert path.isfile(path.join(tb.root, "jobfile.sh"))
        assert path.isfile(path.join(tb.root, "train.xyz"))
    for tb in _2bB:
        assert path.isfile(path.join(tb.root, "jobfile.sh"))
        assert path.isfile(path.join(tb.root, "train.xyz"))    

@pytest.mark.skip()
def test_split(Pd_db):
    """Tests the splitting of available data for the Pd database.
    """
    pass
