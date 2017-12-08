"""Tests the legacy database interface.
"""
import pytest
from matdb.database.phonon import DynMatrix, Calibration, Modulation
from matdb.utility import relpath
from os import mkdir, path, symlink
import quippy
import numpy as np

def compare_dicts(dict1,dict2):
    """Compares two dictionaries to see if they are the same.
    """

    if dict1.keys() != dict2.keys():
        return False

    for key in dict1:
        if not np.allclose(dict1[key],dict2[key]):
            return False

    return True

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
    assert compare_dicts(result,temp)

    dynmatdb.dmatrix

def test_forcecalc(dynmatdb):
    """Tests the force calculation for the dynamical matrix.
    """
    src = relpath("./tests/data/Pd/complete/vasprun.xml__Pd.phonon-16")
    target = path.join(dynmatdb.root,"W.1","vasprun.xml")
    symlink(src,target)
    dynmatdb.calc_forcesets()
    assert path.isfile(path.join(dynmatdb.root,"phonopy","FORCE_SETS"))
    
    with open(path.join(dynmatdb.root,"phonopy","FORCE_SETS"),'r') as f:
        temp1 = f.read()
    with open(relpath("./tests/data/Pd/dynmatrix/FORCE_SETS__Pd.phonon-16"),'r') as f:
        temp2 = f.read()

    assert temp1==temp2
    
    
def test_calc_bands(dynmatdb):
    """Tests a band calculation and the construction of the band.yaml file.
    """

    src = relpath("./tests/data/Pd/dynmatrix/FORCE_SETS__Pd.phonon-16")
    target = path.join(dynmatdb.root,"phonopy","FORCE_SETS")
    symlink(src,target)

    dynmatdb.calc_bands()

    assert path.isfile(path.join(dynmatdb.root,"phonopy","band.yaml"))
    
    from matdb.phonons import from_yaml
    temp1 = from_yaml(relpath("./tests/data/Pd/dynmatrix/band.yml__Pd.phonon-16"))
    temp2 = from_yaml(path.join(dynmatdb.root,"phonopy","band.yaml"))
    assert compare_dicts(temp1,temp2)
    assert compare_dicts(temp1,dynmatdb.bands)
    dynmatdb.calc_bands()

# def test_DOS(dynmatdb):
#     """Tests of the DOS calculation.
#     """

    
