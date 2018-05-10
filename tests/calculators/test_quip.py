"""Tests the vasp calculator in matdb.
"""
import pytest
from matdb.atoms import Atoms
#from matdb.calculators import Quip
Quip = pytest.importorskip('matdb.calculators.Quip')
import six
import numpy as np
from ase import Atoms as aseAtoms

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
        if not isinstance(dict1[key],(six.string_types,list)) and not np.allclose(dict1[key],dict2[key]):
            return False
        elif isinstance(dict1[key],(six.string_types,list)) and not dict1[key] == dict2[key]:
            return False

    return True

def test_quip_steup():
    """Tests the initialization of the quip calculator.
    """

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    potSW = Quip(atSi, '.', '.', 0, ["IP SW"])

    potSW.create()
    
    assert potSW.name == "Quip"
    assert potSW.ran_seed == 0

def test_methods():
    """Tests various Quip methods.
    """
    
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    potSW = Quip(atSi, '.', '.', 0, ["IP SW"])

    assert potSW.can_execute()
    assert potSW.can_extract()
    assert not potSW.is_executing()

    quip_dict = {"folder":'.', "ran_seed":0, "contr_dir":'.', "kwargs": {}, 
                 "args": ["IP SW"]}
    out = potSW.to_dict()
    assert compare_nested_dicts(out, quip_dict)

def test_convert_types():
    """Tests of the conversion to ase atoms types.
    """

    atSi = Atoms("Si",positions=[[0,0,0]],
                 cell=[5.43,5.43,5.43])
    atSi.add_property("force",[[1,2,3]])
    potSW = Quip(atSi, '.', '.', 0, ["IP SW"])
    new_atoms = potSW._convert_atoms(atSi)

    assert isinstance(new_atoms, aseAtoms)
    assert new_atoms.n == 1
    assert np.allclose(new_atoms.positions,atSi.positions)
    assert np.allclose(new_atoms.properties["force"],np.transpose(atSi.force))

def test_calc():
    """Tests the calc method.
    """
    atSi = Atoms("Si",positions=[[0,0,0]],
                 cell=[5.43,5.43,5.43])

    atSi.add_property("force",[[1,2,3]])
    atSi.add_param("energy",12345)
    atSi.add_param("virial",[1,2,3,4,5,6,7,8,9])
    potSW = Quip(atSi, '.', '.', 0, ["IP SW"])
    atSi = Atoms("Si",positions=[[0.1,0.1,0.1]],
                 cell=[5.43,5.43,5.43])
    atSi.add_property("force",[[1,2,3]])
    atSi.add_param("energy",12345)
    atSi.add_param("virial",[1,2,3,4,5,6,7,8,9])
    potSW.calc(atSi)

    assert 1==1
