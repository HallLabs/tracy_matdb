"""Tests the legacy database interface.
"""
import pytest
from matdb.database.phonon import DynMatrix, Calibration, Modulation
from matdb.utility import relpath
from os import mkdir, path, symlink
import quippy
import numpy as np

@pytest.fixture()
def Pd(tmpdir):
    from matdb.utility import relpath
    from matdb.database.controller import Controller
    from os import mkdir, symlink, remove

    target = relpath("./tests/Pd/matdb")
    dbdir = str(tmpdir.join("pd_db"))
    mkdir(dbdir)
    
    #We need to copy the POSCAR over from the testing directory to the temporary
    #one.
    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")
    copy(POSCAR, dbdir)
    
    result = Controller(target, dbdir)
    return result

@pytest.fixture
def dynmatdb(Pd):
    """Returns a DynMatrix database for the Pd dynmatrix database.
    """
    Pd.setup()
        
    return Pd["Pd.phonon-16"].steps['dynmatrix']

def test_setup_and_cleanup(dynmatdb):
    """Tests the setup and cleanup funtions that have been overwritten
    from the basic Database.
    """
    
    dynmatdb.setup()
    dynmatdb.setup(rerun=True)

    dynmatdb.cleanup()

def test_dmatrix(dynmatdb):
    """Tests the dmatrix extraction.
    """

    result = {'eigvals': np.array([ 0.,  0., -0.]),
              'dynmat': np.array([[ 0.+0.j,  0.+0.j,  0.+0.j],
                                  [-0.+0.j,  0.+0.j,  0.+0.j],
                                  [-0.+0.j, -0.+0.j, -0.+0.j]]),
              'eigvecs': np.array([[ 1.+0.j,  0.+0.j,  0.+0.j],
                                   [ 0.+0.j,  1.+0.j,  0.+0.j],
                                   [ 0.+0.j,  0.+0.j,  1.+0.j]])}
    src = relpath("./tests/data/Pd/dynmatrix/FORCE_SETS__Pd.phonon-16")
    target = path.join(dynmatdb.root,"phonopy","FORCE_SETS")
    symlink(src,target)

    temp = dynmatdb.dmatrix
    assert np.all(result['eigvals'] == temp['eigvals'])
    assert np.all(result['dynmat'] == temp['dynmat'])
    assert np.all(result['eigvecs'] == temp['eigvecs'])

    dynmatdb.dmatrix
