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

def test_execution():
    """Tests the execution via shell subprocess in a different folder.
    """
    from matdb.utility import execute, reporoot
    sargs = ["pwd"]
    target = path.join(reporoot, "matdb")
    xres = execute(sargs, target, nlines=1, env_vars={"VASP_PP_PATH":"~/."})
    assert xres["output"][0].strip() == target

    sargs = ["cat dummy-file"]
    target = path.join(reporoot, "matdb")
    xres = execute(sargs, target, nlines=1)
    assert len(xres["error"]) > 0

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
    assert np.allclose(special_values("[0,5,10,12]"),[1, 2, 3, 4, 5, 11])

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

def test_get_suffix():
    """Tests uncovored lines in get_suffix.
    """
    from matdb.utility import get_suffix
    
    d = {"A":10,"B":20}
    k = "A"
    index=1
    values = 10
    assert get_suffix(d,k,index,values) == '-2'

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
    assert pgrid.__eq__(set(pgrid))

    assert len(pgrid) == 12
    assert compare_nested_dicts(pgrid['dos-tt-dim-16.00'],
                         {'phonopy': {'dim': [0, 2, 2, 2, 0, 2, 2, 2, 0]},
                          'dosmesh': [12, 12, 12], 'bandmesh': [13, 13, 13]})
    pgrid.add('dos-tt-dim-16.00',10)
    assert pgrid['dos-tt-dim-16.00'] != 10

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

def test_linecount():
    """Tests the linecount method in utility.
    """
    from matdb.utility import linecount

    assert linecount("temp") == 0

    with open("temp.txt","w+") as f:
        f.write("This is a test file. \n It is rather boring. \n. Thanks for listening.")
    assert linecount("temp.txt") == 3
    remove("temp.txt")

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

    al = Atoms("Si",positions=[[0,0,0]])
    c = Atoms("C",positions=[[0,0,0]])
    k = "Si.positions"

    temp = [{"Si":al},{"c":c}]
    temp = obj_update(temp,k,[[0.5,0.5,0.5]])

    assert np.allclose(temp[0]["Si"].positions,[[0.5,0.5,0.5]])

    al = Atoms("Si",positions=[[0,0,0]])
    temp = {"Si":[[0,0,0]]}
    k = "Si"
    temp = obj_update(temp,k,[[0.5,0.5,0.5]],copy=False)
    assert np.allclose(temp["Si"],[[0.5,0.5,0.5]])

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

def test_is_uuid4():
    """Tests the is_uuid4 funciton.
    """
    from matdb.utility import is_uuid4
    
    assert is_uuid4('4b602114-858d-455d-8152-27a2683af17e')
    assert is_uuid4('4b602114858d455d815227a2683af17e')
    assert not is_uuid4(0)
    assert not is_uuid4('DynMatrix/phonon/Pd/dim-2.00')
    assert not is_uuid4('a0b1c2d3e4f5ghijklmnopqrstuvwxyz')

def test_redirect_stdout():
    """Tests the redirection of stdout.
    """
    from matdb.utility import redirect_stdout
    
    with open("temp.txt",'w') as f:
        with redirect_stdout(f):
            print("Wow")

    with open("temp.txt",'r') as f:
        temp = f.readline().strip()

    assert temp=="Wow"
    remove("temp.txt")

    with pytest.raises(IOError):
        with open("temp.txt",'r') as f:
            with redirect_stdout(f):
                print("Wow")
        
def test_execute():
    """Tests missing lines of execute.
    """
    from matdb.utility import execute

    # test early breaking on stderr and stdout.
    temp=execute(('which','python'),'.',nlines=0)
    temp=execute(('python','enum.x'),'.',nlines=0)
    
    temp=execute(('which','python'),'.',env_vars={"POTENTIALS_DIR":'1'})
    assert "bin/python\n" in temp['output'][0]
    assert temp['error'] == []

def test_load_datetime():
    """Tests missing cases from load datetime pairs.
    """
    from matdb.utility import load_datetime

    data=[[1,10]]
    out = load_datetime(data)

    assert out[1] == 10

    data = [[1,[10,12]]]
    out = load_datetime(data)
    assert out[1] == [10,12]
            

def test_dbcate():
    """Tests missing lines in dbcat.
    """
    from matdb.utility import dbcat
    # Test to make just that 'temp2.txt.json' doesn't get written if
    # files can't be cated.
    dbcat(['temp1.txt'],'temp2.txt',sources=["temp3.txt"])

    assert not path.isfile('temp2.txt.json')
    remove("temp2.txt")
    
    dbcat(['temp1.txt'],'temp2.txt')

    assert path.isfile('temp2.txt')
    remove("temp2.txt")

def test_getattrs():
    """Tetsts the getting of attributes from a chain of attributes.
    """
    from matdb.utility import getattrs
    obj = {"a":{"b":20}}
    assert 20 == getattrs(obj,'a.b')

    from matdb.atoms import Atoms
    at = Atoms("C4")
    assert np.allclose(np.array([[ 0.,  0.,  0.],
                                 [ 0.,  0.,  0.], [ 0.,  0.,  0.]]),
                       getattrs(at,'cell'))
    
def test_chdir(tmpdir):
    """Tests the chdir context manager.
    """

    from matdb.utility import chdir
    from os import getcwd, mkdir
    
    target = str(tmpdir.join("chdir"))
    if not path.isdir(target):
        mkdir(target)
    
    with chdir(target):
        assert getcwd() == target

    cur_dir = getcwd()
    try:
        with chdir(str(tmpdir.join('not_chdir'))):
            l = 1
    except:
        l = 1
    assert getcwd() == cur_dir

def test_dict_update():
    """Tests dictionary update.
    """

    from matdb.utility import dict_update
    
    a = {"one":1, "two":2, "four":{"b":2}}
    b = {"two":2.01, "three":3, "four":{"a":1}}

    dict_update(a,b)
    out = {"one":1, "two":2, "three":3, "four":{"a":1, "b":2}}

    assert compare_nested_dicts(a,out)

def test_rel_path():
    """Tests the relative path. 
    """

    from matdb.utility import relpath

    temp = relpath("./tests")
    assert "matdb/tests" in temp

def test_compare_tree(tmpdir):
    """Tests the folder comparison method.
    """

    from os import mkdir
    from matdb.utility import touch, compare_tree
    test_dir = str(tmpdir.join("comp_tree"))
    mkdir(test_dir)

    touch(path.join(test_dir,"compute.pkl"))
    touch(path.join(test_dir,"jobfile.sh"))
    mkdir(path.join(test_dir,"phonopy"))
    touch(path.join(test_dir,"phonopy","POSCAR"))
    
    folders = {
        "__files__": ["compute.pkl","jobfile.sh"],
        "phonopy": {
            "__files__": ["POSCAR"]
        },
    }

    compare_tree(test_dir,folders)

def test_pgrid():
    """Tests the paramater grid generation.
    """

    from matdb.utility import pgrid

    opts = {"dim*":[1],"tum":[2],"leg":3}
    ignore = ["leg"]

    grid, keys = pgrid(opts,ignore)

    assert keys == ["dim","tum"]
    assert grid == [(1,[2])]

def test_dict_to_str():
    """Tests convertion of dict to str.
    """
    from matdb.utility import convert_dict_to_str
    
    in_dict = {"a":1, "b":{"c":2}}
    out = convert_dict_to_str(in_dict)

    assert out == "'a':'1';'b':''c':'2';';"

def test_check_deps():
    """Tests the dependency checker.
    """

    from matdb.utility import required_packages, check_deps

    assert ["ase", "beautifulsoup4", "certifi", "chardet", "cycler", "h5py", "html5lib",
            "idna", "matplotlib", "mpld3", "numpy", "phenum", "phonopy", "pyparsing",
            "python-dateutil", "pytz", "PyYAML", "requests", "subprocess32", "termcolor", 
            "tqdm", "urllib3", "webencodings", "seekpath"] == required_packages()

    res = check_deps()

    print(res)
    assert res["numpy"] is not None
