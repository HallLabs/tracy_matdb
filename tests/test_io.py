"""Tests the I/O functions for template hierarchy parsing in matdb.
"""
import pytest
from matdb.io import read
from matdb.utility import reporoot, relpath
from matdb.atoms import AtomsList
from matdb.io import cfg_to_xyz
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

def test_cfg(tmpdir):
    """Tests conversion of MTP's CFG format to XYZ.
    """
    target = str(tmpdir.join("cfg.xyz"))
    model = AtomsList(relpath("tests/files/io_convert/atoms.xyz"))
    conv = cfg_to_xyz(relpath("tests/files/io_convert/atoms.cfg"), target,
               species=[46, 47])

    for a, b in zip(model, conv):
        assert np.allclose(a.get_positions(), b.get_positions())
        assert np.allclose(a.get_forces(), b.calc.results["forces"])
        assert np.allclose(a.get_stress(), b.calc.results["stress"], atol=1e-4, rtol=1e-3)
        assert a.get_total_energy() == b.calc.results["energy"]
        
def test_schema():
    """Tests the template read for the whole schema to make sure that contexts
    work as expected.
    """
    varfull = read(reporoot, "tests/io/schema")
    #The full dictionary is large and unwieldly and gets tested below in
    #separate pieces. For now, we just test that the context hookups worked
    #correctly. Edges uses a relative path specifier with `../` syntax, so it is
    #most likely to mess up.
    model = [
        {'targets': 'B',
         'name': 'AtoB',
         'doc': 'Connects A to B.',
         'sources': 'A',
         'properties': [
             {'dtype': 'int',
              'doc': 'some integer.',
              'example': 3,
              'name': 'value'}]
        }]
    assert varfull["edges"] == model

def test_raises_errors():
    """Makes sure errors are raised where appropriate.
    """
    with pytest.raises(ValueError):
        template = read(reporoot, "tests/dummy")

def test_corner_cases():
    """Tests the corner cases in the template reader.
    """
    varcorner = read(reporoot, "tests/io/corner")
    model = {
        'first': {'a': 0, 'b': 1},
        'second': [{'a': 0, 'b': 1}]
    }
    assert varcorner == model
        
def test_read_types():
    """Tests the recursive read of a directory of edge template files with
    sub-directories, relative paths, etc.
    """
    varkind = read(reporoot, "tests/io/types/varkind")
    model = {
        'name': 'VarKind',
        'properties': [
            {'name': 'tolerance',
             'example': 0.2,
             'dtype': 'float',
             'doc': 'heard of epsilon.'}
        ],
        'doc': 'Simple custom type for unit testing.'
    }
    assert varkind == model
        
def test_read_edges():
    """Tests the recursive read of a directory of edge template files with
    sub-directories, relative paths, etc.
    """
    AtoB = read(reporoot, "tests/io/edges/atob")
    model = {
        'targets': 'B',
        'properties': [
            {'dtype': 'int',
             'doc': 'some integer.',
             'name': 'value',
             'example': 3}],
        'doc': 'Connects A to B.',
        'name': 'AtoB',
        'sources': 'A'
    }
    assert AtoB == model
        
def test_read_verts():
    """Tests the recursive read of a directory of vertex template files with
    sub-directories, relative paths, etc.
    """
    A = read(reporoot, "tests/io/verts/a")
    modelA = {
        'name': 'A',
        'properties': [
            {'name': 'kind',
             'dtype': 'VarKind',
             'example': 'kind(1.8)',
             'keytype': 'float',
             'doc': 'testing how kind the vertex is.'
            }],
        'doc': 'Simple vertex for unit tests.'
    }
    assert A == modelA

    B = read(reporoot, "tests/io/verts/b")
    modelB = {
        'vtype': 'logic',
        'doc': 'Simple schema for unit tests.  ',
        'properties': [
            {'example': 3,
             'doc': 'some integer.',
             'dtype': 'int',
             'name': 'value'},
            {'example': "'cc'",
             'doc': 'b.text',
             'dtype': 'str',
             'name': 'text'}
        ],
        'name': 'B'}
    assert B == modelB
    
def test_hdf5_in_out():
    """Tests the writing of dictionaries to hdf5 and reading back out.
    """

    import h5py
    from matdb.io import load_dict_from_h5, save_dict_to_h5
    from os import remove

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

def test_vasp_xyz(tmpdir):
    """Tests vasp to xyz function.
    """
    from matdb.io import vasp_to_xyz
    from os import remove

    target = "vasp.xyz"
    model = vasp_to_xyz(relpath("tests/files/io_convert/"), target,
                        properties=["species", "pos", "z"],
                        parameters=["energy", "virial"], config_type="temp")

    assert model

    model = vasp_to_xyz(relpath("tests/files/io_convert/"), target, config_type="temp")
    assert model

    remove(relpath("tests/files/io_convert/vasp.xyz"))
    remove(relpath("tests/files/io_convert/vasp.xyz.idx"))   

def test_atoms_to_cfg(tmpdir):
    """Tests the writing of an atoms object to cfg format. 
    """
    from matdb.io import atoms_to_cfg
    from matdb.atoms import Atoms
    from os import path, remove
    from ase.calculators.singlepoint import SinglePointCalculator
    from matdb.calculators import Vasp

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    target = "train.cfg"

    atoms_to_cfg(atSi, target)

    assert path.isfile(target)

    remove(target)
    
    energy = 10.0
    forces = [[1,1,1],[2,2,2],[3,3,3],[4,4,4],[5,5,5],[6,6,6],[7,7,7],[8,8,8]]
    stress = [0,1,2,3,4,5]

    atSi.calc = SinglePointCalculator(atSi, energy=energy, forces=np.array(forces),
                                          stress=stress)

    atSi.add_property("energy", energy)
    atSi.add_property("stress", stress)
    atSi.add_param("force", forces)
    atoms_to_cfg(atSi, target)

    assert path.isfile(target)
    
    remove(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    type_map = {0:1}
    
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp",
                                                 "versions": {"Si":"05Jan2001"}}, "xc":"pbe"}
    atSi.calc = Vasp(atSi, str(tmpdir.join("Vasp")) , str(tmpdir), 0, **kwargs)

    atoms_to_cfg(atSi, target, type_map=type_map, config_id="test")
    assert path.isfile(target)
    
    remove(target)
