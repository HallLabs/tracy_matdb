"""Tests the database descriptor functions.
"""
import pytest
from matdb.utility import relpath
from matdb.database import basic
import matdb.database.descriptors as ds
from pandas.util.testing import assert_frame_equal
import pandas as pd

@pytest.fixture()
def atoms():
    """An atoms object for the descriptors to work on.
    """

    POSCAR = relpath("./tests/Pd/POSCAR")
    from quippy.atoms import Atoms
    datoms = Atoms(POSCAR,format="POSCAR")

    return datoms

def test_d2b(atoms):
    """Tests the histogram of the 2-body terms."""
        
    result = ds.d2b([atoms],6)
    test_path = relpath("./tests/data/output/d2b_out1.csv")
    data = pd.read_csv(test_path,sep="\t")
    assert_frame_equal(data,result)

    result = ds.d2b([atoms],3)
    test_path = relpath("./tests/data/output/d2b_out2.csv")
    data = pd.read_csv(test_path,sep="\t")
    assert_frame_equal(data,result)

def test_d3b(atoms):
    """Tests the histogram of the 2-body terms."""
    
    result = ds.d3b([atoms],6)
    test_path = relpath("./tests/data/output/d3b_out1.csv")
    data = pd.read_csv(test_path,sep="\t")
    assert_frame_equal(data,result)

    result = ds.d3b([atoms],3)
    test_path = relpath("./tests/data/output/d3b_out2.csv")
    data = pd.read_csv(test_path,sep="\t")
    assert_frame_equal(data,result)
    
def test_dNb(atoms):
    """Tests the histogram of the 2-body terms."""
    
    result = ds.dNb([atoms],{2:6,3:4})
    test_path = relpath("./tests/data/output/dNb_out1.csv")
    data = pd.read_csv(test_path,sep="\t")
    assert_frame_equal(data,result.reset_index(drop=True))

    result = ds.dNb([atoms],{2:3,3:6})
    test_path = relpath("./tests/data/output/dNb_out2.csv")
    data = pd.read_csv(test_path,sep="\t")
    assert_frame_equal(data,result.reset_index(drop=True))
