"""Tests the vasp calculator in matdb.
"""
import pytest
from matdb.atoms import Atoms
from matdb.calculators import Aflow
import six
import numpy as np
from os import path
from types import NoneType

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
        if not isinstance(dict1[key],(six.string_types,list,NoneType)) and not np.allclose(dict1[key],dict2[key]):
            return False
        elif isinstance(dict1[key],(six.string_types,list,NoneType)) and not dict1[key] == dict2[key]:
            return False

    return True

def test_init(paper):
    """Tests the initialize function and some other subroutines.
    """
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    entry = paper[0]
    kwargs  = {}
    args = []

    calc = Aflow(atSi, '.', '.', 0, entry=entry, *args, **kwargs)

    assert calc.entry_file == "./entry.pkl"
    assert calc.can_execute()
    assert calc.can_extract()
    assert not calc.is_executing()
    calc.extract()

def test_to_dict():
    """Tests the to dict method.
    """
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    kwargs  = {}
    args = []

    calc = Aflow(atSi, '.', '.', 0, *args, **kwargs)
    out = {"folder":'.', "ran_seed":0, "contr_dir":'.', "kwargs":{"entry":None}, "args":[]}
    assert compare_nested_dicts(calc.to_dict(),out)

def test_calc(paper,tmpdir):
    """Tests the create method.
    """
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    entry = paper[0]
    kwargs  = {}
    args = []

    calc = Aflow(atSi, str(tmpdir), '.', 0, entry=entry, *args, **kwargs)

    calc.create()
    assert path.isfile(path.join(tmpdir,calc.entry_file))
