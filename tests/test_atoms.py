# -*- coding: utf-8 -*-
"""Tests the atoms object and related functions.
"""
import pytest
from os import path, remove, mkdir
import numpy as np
import six

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

def test_recursive_convert():
    """Tests the recursive unit conversion.
    """
    from matdb.atoms import _recursively_convert_units
    
    dict_in = {"a":10,"b":0.34,"c":[1,2,3]}
    dict_out = {"a":np.int64(10),"b":np.float64(0.34),"c":np.array([1,2,3])}

    test = _recursively_convert_units(dict_in)
    
    assert compare_nested_dicts(test,dict_out)
    
    dict_in = {"a":10,"b":{"B":0.34,"C":np.array([1,4]),"D":{"E":111,"F":[0,2]}},
               "c":{"P":[1,2,3],"I":0}}
    dict_out = {"a":np.int64(10),"b":{"B":np.float64(0.34),"C":np.array([1,4]),
                                      "D":{"E":np.int64(111),"F":np.array([0,2])}},
               "c":{"P":np.array([1,2,3]),"I":np.int64(0)}}

    test = _recursively_convert_units(dict_in)
    
    assert compare_nested_dicts(test,dict_out)
    
    dict_in = {"a":10,"b":0.34,"c":[1,2,3]}
    dict_out = {"a":np.float64(10),"b":np.int64(0.34),"c":np.array([1,2,3])}

    test = _recursively_convert_units(dict_in)
    
    assert not compare_nested_dicts(test,dict_out)

def test_hdf5(tmpdir):
    """Tests whether an atoms object with calculated parameters can be saved to
    JSON and then restored.
    """
    from matdb.calculators import Quip
    from matdb.atoms import Atoms
    target = str(tmpdir.join("to_hdf5"))
    if not path.isdir(target):
        mkdir(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    potSW = Quip(atSi, target, ["IP SW"])
    atSi.set_calculator(potSW)
    potSW.calc(atSi, energy=True, force=True, virial=True)
    atSi.properties["rand"] = np.random.randint(0, 100, 8)
    atSi.write(target=path.join(target,"temp.h5"))
    atR = Atoms()
    atR.read(target=path.join(target,"temp.h5"))

    assert atR.energy == atSi.energy
    assert isinstance(atR, Atoms)
    assert np.allclose(atR.force, atSi.force)
    assert np.allclose(atR.virial, atSi.virial)
    assert np.allclose(atR.properties["rand"], atSi.properties["rand"])
    assert np.allclose(atR.positions, atSi.positions)
    remove(path.join(target,"temp.h5"))
