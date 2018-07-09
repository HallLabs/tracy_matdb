"""Tests the vasp calculator in matdb.
"""
import pytest
from os import path, mkdir, remove
import six
import numpy as np
import json

from matdb.atoms import Atoms
from matdb.utility import reporoot, relpath, symlink, touch
from matdb.calculators import TracyQE
from matdb.exceptions import VersionError

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
        if not isinstance(dict1[key],(six.string_types, list)) and not np.allclose(dict1[key],dict2[key]):
            return False
        elif isinstance(dict1[key],six.string_types) and not dict1[key] == dict2[key]:
            return False
        elif isinstance(dict1[key],list):
            if not len(dict1[key]) == len(dict2[key]):
                return False
            else:
                return all([dict1[key][i] == dict2[key][i] for i in range(len(dict1[key]))])
            
    return True

def test_get_calc_mod():
    """Tests the __init__ modules get_calculator_module.
    """

    from matdb.calculators import get_calculator_module
    from types import ModuleType
    
    mod = get_calculator_module({"name":"TracyQE"})
    assert isinstance(mod, ModuleType)

def test_tracy_qe_setup(tmpdir):
    """Tests TracyQE calculator initialization.
    """

    target = str(tmpdir.join("TracyQE"))
    mkdir(target)
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])

    kwargs = {"calcargs":{"potcars": {"directory":"./tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
                          "kpoints":{"method": "MP", "divisions": (3, 3, 3)},
                          "input_data": {"control":{"calculation": "relax", "prefix": "test"}}},
              "tracy": {"max_time": 10, "min_flops": 10, "min_ram": 10, "min_mem": 10,
                        "ncores": 1, "role": "Enumerated", "notifications": "stuff"}}

    with open(path.join(target, "contract.txt"), "w+") as f:
        f.write("12345")
    with open(path.join(target, "post_print.txt"), "w+") as f:
        f.write("abcde")
    authenticate = {"user": "user", "password": "password"}
    with open(str(tmpdir.join("user_cred.json")), "w+") as f:
        json.dump(authenticate,f)
    calc = TracyQE(atm, target, '.', 0, **kwargs)

    assert isinstance(calc, TracyQE)
    assert calc.parameters["input_data"]["control"]["calculation"] == "relax"
    assert "kpts" in calc.parameters
    assert calc.out_file == "test"
    assert calc.group_preds is None
    assert calc.sys_specs["max_time"] == 10
    assert calc.contract_id == "12345"
    assert calc.after_print == "abcde"
    remove(str(tmpdir.join("user_cred.json")))

def test_Tracy(tmpdir):
    """Tests the writing of the input files for the TracyQE calculator.
    """

    target = str(tmpdir.join("TracyQE"))
    mkdir(target)
    atm = Atoms("Al",positions=[[0, 0, 0]])#, cell=[[1.0.0],[0,1,0],[0,0,1]])

    kwargs = {"calcargs":{"potcars": {"directory":"./tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
                          "versions": {"Al": ["2.0.1", ".5.1"]}},
                          "kpoints":{"method": "MP", "divisions": (3, 3, 3)},
                          "input_data": {"control":{"calculation": "relax", "prefix": "test"}}},
              "tracy": {"max_time": 10, "min_flops": 10, "min_ram": 10, "min_mem": 10,
                        "ncores": 1, "role": "Enumerated", "notifications": "stuff"}}

    authenticate = {"user": "user", "password": "password"}
    with open(str(tmpdir.join("user_cred.json")), "w+") as f:
        json.dump(authenticate,f)
    calc = TracyQE(atm, target, '.', 0, **kwargs)

    with pytest.raises(ValueError):
        compressed = calc._compress_struct(calc.atoms)

    atm = Atoms("Al",positions=[[0, 0, 0]], cell=[[1,0,0],[0,1,0],[0,0,1]])
    calc = TracyQE(atm, target, '.', 0, **kwargs)
    compressed = calc._compress_struct(calc.atoms)

    assert np.allclose(compressed["a"], atm.cell)
    assert np.allclose(compressed["b"], [[0,0,0]])
    assert np.allclose(compressed["t"], 1)
    assert np.allclose(compressed["h"], 201020101020)

    calc.is_executing(target)
    calc.can_extract(target)

    source = calc._get_source()
    assert source == 3

    tmp_source = calc.role
    calc.role = "dump"
    with pytest.raises(ValueError):
        calc._get_source()
    remove(str(tmpdir.join("user_cred.json")))

def test_TracyQE(tmpdir):
    """Tests the Tracy QE methods."""
    
    target = str(tmpdir.join("TracyQE"))
    mkdir(target)
    atm = Atoms("Al",positions=[[0, 0, 0]], cell=[[1,0,0],[0,1,0],[0,0,1]])

    kwargs = {"calcargs":{"potcars": {"directory":"./tests",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
                          "versio2ns": {"Al": ["2.0.1", ".5.1"]}},
                          "kpoints":{"method": "MP", "divisions": (3, 3, 3)},
                          "input_data": {"control":{"calculation": "relax", "prefix": "test"}}},
              "tracy": {"max_time": 10, "min_flops": 10, "min_ram": 10, "min_mem": 10,
                        "ncores": 1, "role": "Enumerated", "notifications": "stuff",
                        "group_preds": "abcd", "contract_preds":"1234", "notifications": "stuff"}}

    authenticate = {"user": "user", "password": "password"}
    with open(str(tmpdir.join("user_cred.json")), "w+") as f:
        json.dump(authenticate,f)
    calc = TracyQE(atm, target, '.', 0, **kwargs)

    calc.write_input(atm)

    assert path.isfile(path.join(target, "espresso.pwi"))
    assert path.isfile(path.join(target, "submission.json"))
    
    assert calc.can_execute(target)

    calc.write_input(atm)

    calc.extract(target)

    Qe_dict = calc.to_dict()

    assert "kwargs" in Qe_dict
    assert Qe_dict["folder"] == target
    assert Qe_dict["ran_seed"] == 0
    remove(str(tmpdir.join("user_cred.json")))

def test_array_to_int():
    """Tests conversion of arrays to ints."""
    atm = Atoms("Al",positions=[[0, 0, 0]], cell=[[1,0,0],[0,1,0],[0,0,1]])

    kwargs = {"calcargs":{"potcars": {"directory":"./tests",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
                          "versions": {"Al": ["2.0.1", ".5.1"]}},
                          "kpoints":{"method": "MP", "divisions": (3, 3, 3)},
                          "input_data": {"control":{"calculation": "relax", "prefix": "test"}}},
              "tracy": {"max_time": 10, "min_flops": 10, "min_ram": 10, "min_mem": 10,
                        "ncores": 1, "role": "Enumerated", "notifications": "stuff",
                        "group_preds": "abcd", "contract_preds":"1234"}}

    authenticate = {"user": "user", "password": "password"}
    with open(path.join('.',"user_cred.json"), "w+") as f:
        json.dump(authenticate,f)
    calc = TracyQE(atm, '.', '.', 0, **kwargs)

    in_put = [1,2,3,4,5]
    out_put = 12345

    test_out = calc._intarray_to_int(in_put)

    assert out_put == test_out

    in_put = [1,2,3,4,5]
    out_put = 2030405060

    test_out = calc._intarray_to_int(in_put, pad=True)

    assert out_put == test_out
    remove(str(tmpdir.join("user_cred.json")))
