# -*- coding: utf-8 -*-
"""Tests the utility functions.
"""
import pytest
from os import path, remove
import numpy as np
import six

def compare_dicts(dict1,dict2):
    """Compares two dictionaries to see if they are the same.
    """

    if dict1.keys() != dict2.keys():
        return False

    for key in dict1:
        print("H2")
        if not np.allclose(dict1[key],dict2[key]):
            return False

    return True

def compare_nested_dicts(dict1,dict2):
    """Compares two dictionaries to see if they are the same.
    """

    if sorted(dict1.keys()) != sorted(dict2.keys()):
        return False

    for key in dict1:
        if isinstance(dict1[key],dict):
            res = compare_nested_dicts(dict1[key],dict2[key])
            if not res:
                return False
            else:
                continue
        if not isinstance(dict1[key],six.string_types) and not np.allclose(dict1[key],dict2[key]):
            return False
        elif isinstance(dict1[key],six.string_types) and not dict1[key] == dict2[key]:
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
    
    symlink(source, target)
    result = readlink(source)
    assert result == target

    dirsource = str(tmpdir.join("dummy-dir"))
    from os import mkdir
    mkdir(dirsource)
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

def test_special_functions():
    """Tests the special function evaluation.
    """

    from matdb.utility import special_functions

    assert np.allclose(1,special_functions({"func":"linalg:det","reshape":(3,3)},[1,0,0,0,1,0,0,0,1]))
    assert np.allclose(3,special_functions("numpy:sum",[1,0,0,0,1,0,0,0,1]))
    with pytest.raises(ValueError):
        special_functions(None,None)

def test_is_number():
    """Tests the is_number function.
    """
    from matdb.utility import is_number
    import sys
    
    assert is_number('1.6')
    if sys.version_info >= (3, 0):    
        assert is_number('٥')
    else:
        assert is_number(u'٥')
    assert not is_number('str')


def test_ParameterGrid():
    """Tests the creation of a ParamaterGrid and it's functionality.
    """
    from matdb.utility import ParameterGrid

    db_info = {'phonopy': {'dim*': [[2, 0, 0, 0, 2, 0, 0, 0, 2], [0, 1, 1, 1, 0, 1, 1, 1, 0],
                                    [0, 2, 2, 2, 0, 2, 2, 2, 0], [-2, 2, 2, 2, -2, 2, 2, 2, -2],
                                    [-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 0, 0, 0, 3, 0, 0, 0, 3]],
                           'dim_suffix': {'func':'linalg:det','reshape':[3,3]}},
               'parent':  'controller', 'dosmesh*': [[10, 10, 10],[12,12,12]],
               'dosmesh_suffix*':[30,"tt"],'atoms': 'Pd',
               'root': '/root/codes/matdb/tests/Pd', 'bandmesh': [13, 13, 13]}

    pgrid = ParameterGrid(db_info)

    assert len(pgrid) == 12
    assert compare_nested_dicts(pgrid['dos-tt-dim-16.00'],
                         {'phonopy': {'dim': [0, 2, 2, 2, 0, 2, 2, 2, 0]},
                          'dosmesh': [12, 12, 12], 'bandmesh': [13, 13, 13]})
    pgrid.pop('dos-tt-dim-16.00')
    assert not ('dos-tt-dim-16.00' in pgrid)
    assert not pgrid == ParameterGrid(db_info)

def test_get_grid():
    """Tests the get_grid method.
    """
    from matdb.utility import get_grid
    from numpy import array
    db_info = {"lattice*": ["fcc", "bcc", "hcp"], 
          "calculator": {"encut*": [700, 800, 900]},
          "normal": [1, 2, 3],
          "double": {"single": {"dog*": [range(9), list(np.diag([9,2,1])), list(np.diag([3,2,4]))]}},
          "encut_suffix*": [70,80,90],
          "dog_suffix": {"func": "linalg:det",
                         "reshape": [3 ,3]},
          "lattice_suffix": "{}"
        }

    test = get_grid(db_info)
    model = {'enc-70-dog-0.00-lat-bcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-0.00-lat-bcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-0.00-lat-bcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-18.00-lat-bcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-18.00-lat-bcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-18.00-lat-bcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-24.00-lat-bcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-24.00-lat-bcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-24.00-lat-bcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-0.00-lat-fcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-0.00-lat-fcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-0.00-lat-fcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-18.00-lat-fcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-18.00-lat-fcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-18.00-lat-fcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-24.00-lat-fcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-24.00-lat-fcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-24.00-lat-fcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-0.00-lat-hcp': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-80-dog-0.00-lat-hcp': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-90-dog-0.00-lat-hcp': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-70-dog-18.00-lat-hcp': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-80-dog-18.00-lat-hcp': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-90-dog-18.00-lat-hcp': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-70-dog-24.00-lat-hcp': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-80-dog-24.00-lat-hcp': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-90-dog-24.00-lat-hcp': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]}}

    assert compare_nested_dicts(test,model)
    
def test_hdf5_in_out():
    """Tests the writing of dictionaries to hdf5 and reading back out.
    """

    import h5py
    from matdb.utility import load_dict_from_h5, save_dict_to_h5

    dict_1_in = {"a":{"B":np.int64(1),"C":np.int64(3),"D":{"temp":np.array([10,11,12])}},
                 "n":np.array([3,2]),"t":np.int64(5)}
    dict_2_in = {"a":np.int64(10),"b":np.array([1,2,10])}

    hf = h5py.File("temp.h5","w")
    save_dict_to_h5(hf,dict_1_in,"/")
    hf.close()

    hf = h5py.File("temp.h5","r")
    out = load_dict_from_h5(hf)
    hf.close()
    assert compare_nested_dicts(dict_1_in,out)
    remove("temp.h5")
    
    hf = h5py.File("temp.h5","w")
    save_dict_to_h5(hf,dict_2_in,"/")
    hf.close()

    hf = h5py.File("temp.h5","r")
    out = load_dict_from_h5(hf)
    hf.close()
    assert compare_dicts(dict_2_in,out)
    remove("temp.h5")

    hf = h5py.File("temp.h5","w")

    with pytest.raises(ValueError):
        save_dict_to_h5(hf,{"a":2},"/")
    hf.close()
    remove("temp.h5")

def test_linecount():
    """Tests the linecount method in utility.
    """
    from matdb.utility import linecount

    assert linecount("temp") == 0

def test_safeupdate():
    """Tests the safe_update method in utility.
    """

    from matdb.atoms import Atoms
    from matdb.utility import safe_update
    
    al = Atoms("Si",positions=[[0,0,0]])
    al.add_param("energy",None)
    kv = {"positions":[[0.5,0.5,0.5]],"energy":10}
    safe_update(al,kv)
    
    assert np.allclose(al.positions,[[0,0,0]])
    assert al.energy == 10

def test_objupdate():
    """Tests the obj_update method in utility.
    """

    from matdb.atoms import Atoms
    from matdb.utility import obj_update
    
    al = Atoms("Si",positions=[[0,0,0]])
    k = "positions"
    
    al = obj_update(al,k,[[0.5,0.5,0.5]])

    assert np.allclose(al.positions,[[0.5,0.5,0.5]])

def test_copyonce():
    """Tests the copyonce method in utility.
    """

    from matdb.utility import copyonce, touch
    from os import remove, path

    touch("temp1.txt")
    copyonce("temp1.txt","temp2.txt")

    assert path.isfile("temp2.txt")

    remove("temp1.txt")
    remove("temp2.txt")

def test_which():
    """Tests th which method in utility.
    """

    from matdb.utility import which

    assert which('/bin/rm')=='/bin/rm'

def test_parse_date():
    """Tests the date parser.
    """

    from matdb.utility import parse_date
    from dateutil import parser

    dates = ['10-2-1987','10-4-1894']
    temp = parse_date(dates)
    assert len(temp)==2

    for i in range(2):
        assert parser.parse(dates[i]) == temp[i]

    with pytest.raises(ValueError):
        parse_date((10,2,1987))
