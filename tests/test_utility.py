"""Tests the utility functions.
"""
import pytest
from os import path
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

def test_execute():
    """Tests the execution via shell subprocess in a different folder.
    """
    from matdb.utility import execute, reporoot
    sargs = ["pwd"]
    target = path.join(reporoot, "matdb")
    xres = execute(sargs, target, nlines=1)
    assert xres["output"][0].decode("UTF-8").strip() == target

    sargs = ["cat dummy-file"]
    target = path.join(reporoot, "matdb")
    xres = execute(sargs, target, nlines=1)
    assert len(xres["error"]) > 0

def test_cat(tmpdir):
    """Tests concatenation of multiple files.
    """
    from matdb.utility import cat, execute, reporoot
    files = [path.join(reporoot, "tests", "files", f)
             for f in ["A.txt", "B.txt"]]
    outfile = str(tmpdir.join("cat_C.txt"))
    cat(files, outfile)

    sargs = ["diff", "C.txt", outfile]
    xres = execute(sargs, path.join(reporoot, "tests/files"))
    assert len(xres["output"]) == 0

def test_symlink(tmpdir):
    """Tests symbolic linking of a file.
    """
    from matdb.utility import reporoot, symlink, execute
    target = path.join(reporoot, "__init__.py")
    source = str(tmpdir.join("symlink_init"))
    symlink(source, target)

    from os import readlink
    result = readlink(source)
    assert result == target

    dirsource = str(tmpdir.join("dummy-dir"))
    assert symlink(dirsource, reporoot) is None

def test_slicer():
    """Tests the slicer of objects.
    """

    from matdb.utility import slicer

    obj = range(100)

    slicer(obj,[1])
    slicer(obj,(1))
    assert [1, 2, 3, 6, 7] == slicer(obj,(1,4,6,8))
    assert [96, 97, 98, 99] == slicer(obj,(-4,None))
    assert [0, 1, 2, 3, 4, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99] == slicer(obj,(0,5,90,None))

def test_py_execute():
    """Tests the execution of python modules.
    """

    from matdb.utility import _py_execute

    assert np.allclose(_py_execute("numpy","linspace","(1,2)"),
                       np.array([ 1.        ,  1.02040816,  1.04081633,  1.06122449,  1.08163265,
                                  1.10204082,  1.12244898,  1.14285714,  1.16326531,  1.18367347,
                                  1.20408163,  1.2244898 ,  1.24489796,  1.26530612,  1.28571429,
                                  1.30612245,  1.32653061,  1.34693878,  1.36734694,  1.3877551 ,
                                  1.40816327,  1.42857143,  1.44897959,  1.46938776,  1.48979592,
                                  1.51020408,  1.53061224,  1.55102041,  1.57142857,  1.59183673,
                                  1.6122449 ,  1.63265306,  1.65306122,  1.67346939,  1.69387755,
                                  1.71428571,  1.73469388,  1.75510204,  1.7755102 ,  1.79591837,
                                  1.81632653,  1.83673469,  1.85714286,  1.87755102,  1.89795918,
                                  1.91836735,  1.93877551,  1.95918367,  1.97959184,  2.        ]
                       )) 

def test_special_values():
    """Tests the special values for matdb.
    """

    from matdb.utility import special_values
    import scipy

    assert np.allclose(special_values("linspace(1,3)"),
                       np.array([ 1.        ,  1.04081633,  1.08163265,  1.12244898,  1.16326531,
                                  1.20408163,  1.24489796,  1.28571429,  1.32653061,  1.36734694,
                                  1.40816327,  1.44897959,  1.48979592,  1.53061224,  1.57142857,
                                  1.6122449 ,  1.65306122,  1.69387755,  1.73469388,  1.7755102 ,
                                  1.81632653,  1.85714286,  1.89795918,  1.93877551,  1.97959184,
                                  2.02040816,  2.06122449,  2.10204082,  2.14285714,  2.18367347,
                                  2.2244898 ,  2.26530612,  2.30612245,  2.34693878,  2.3877551 ,
                                  2.42857143,  2.46938776,  2.51020408,  2.55102041,  2.59183673,
                                  2.63265306,  2.67346939,  2.71428571,  2.75510204,  2.79591837,
                                  2.83673469,  2.87755102,  2.91836735,  2.95918367,  3.        ]
                       ))

    assert special_values(["linspace(1,3)"]) == ["linspace(1,3)"]
    assert special_values("random:choice(10)",seed=1) == 5
    assert special_values("dist:alpha(0,10)") == "dist:alpha(0,10)"
    assert isinstance(special_values("distr:alpha(0,10)"),
                      scipy.stats._distn_infrastructure.rv_frozen)


def test_flatten_dict():
    """Tests the flattening of a dictionary.
    """
    from matdb.utility import flatten_dict
    
    dict1 = {"a":[1,0],"b":[2,0],"c":1.6}
    dict2 = flatten_dict(dict1)
    
    assert compare_dicts(dict1,dict2)

    dict1 = {"a":{"d":[1,0],"e":0},"b":[2,0],"c":2.5}
    dict2 = flatten_dict(dict1)
    dict1 = {"d":[1,0],"e":0,"b":[2,0],"c":2.5}

def test_special_functions():
    """Tests the special function evaluation.
    """

    from matdb.utility import special_functions

    assert np.allclose([1,8],special_functions("linalg:det",[[1,0,0,0,1,0,0,0,1],[2,2,2]]))
    assert np.allclose([3,6],special_functions("numpy:sum",[[1,0,0,0,1,0,0,0,1],[2,2,2]]))
    assert not special_functions(None,None)
    assert not special_functions("stuff",None)
    assert not special_functions("stuff",[1,2])
