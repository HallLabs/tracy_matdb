"""Tests the I/O functions for template hierarchy parsing in matdb.
"""
import pytest
from matdb.atoms import Atoms
from matdb.utility import reporoot, relpath
from matdb.calculators import Vasp
from os import path, mkdir, remove
import six
import numpy as np

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

def test_get_calc_mod():
    """Tests the __init__ modules get_calculator_module.
    """

    from matdb.calculators import get_calculator_module
    from types import ModuleType
    
    mod = get_calculator_module({"name":"Vasp"})
    assert isinstance(mod, ModuleType)

def test_vasp_setup(tmpdir):
    """Tests Vasp calculator initialization.
    """

    target = str(tmpdir.join("Vasp"))        
    atm = Atoms("Si",positions=[[0,0,0]])

    with pytest.raises(ValueError):
        kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp"}}
        calc = Vasp(atm, '.', 'def', 0, **kwargs)

    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp"}, "xc":"pbe"}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    assert calc.potcars["xc"] == "pbe"
    
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    assert calc.kwargs["xc"] == "pbe"    


def test_write_potcar(tmpdir):
    """Tests the writing of the POTCAR file. 
    """
    target = str(tmpdir.join("Vasp"))        
    atm = Atoms("Si",positions=[[0,0,0]])
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    calc._write_potcar()
    assert path.isfile(path.join(calc.contr_dir,"POTCAR"))
    assert path.isfile(path.join(calc.folder,"POTCAR"))
    remove(path.join(calc.contr_dir,"POTCAR"))
    
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp",
                                                 "xc":"pbe", "version":"28June2018"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)
    calc._write_potcar()
    assert path.isfile(path.join(calc.contr_dir,"POTCAR"))
    calc._write_potcar()
    assert path.isfile(path.join(calc.folder,"POTCAR"))

def test_write_input(tmpdir):
    """Tests the writing of the input files. 
    """
    target = str(tmpdir.join("Vasp"))        
    atm = Atoms("Si",positions=[[0,0,0]])
    kwargs = {"kpoints":{"method":"mueller","mindistance":50},
              "potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    calc.write_input(atm,directory=target)
    assert path.isfile(path.join(calc.folder,"POTCAR"))
    assert path.isfile(path.join(calc.folder,"INCAR"))
    assert path.isfile(path.join(calc.folder,"POSCAR"))
    assert path.isfile(path.join(calc.folder,"KPOINTS"))
    
    kwargs = {"potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    calc.write_input(atm,directory=target)
    assert path.isfile(path.join(calc.folder,"POTCAR"))
    assert path.isfile(path.join(calc.folder,"INCAR"))
    assert path.isfile(path.join(calc.folder,"POSCAR"))
    assert path.isfile(path.join(calc.folder,"KPOINTS"))

    calc.create()
    assert path.isfile(path.join(calc.folder,"POTCAR"))
    assert path.isfile(path.join(calc.folder,"INCAR"))
    assert path.isfile(path.join(calc.folder,"POSCAR"))
    assert path.isfile(path.join(calc.folder,"KPOINTS"))

def test_can_execute(tmpdir):
    """Tests the can_execute method.
    """

    target = str(tmpdir.join("Vasp"))        
    atm = Atoms("Si",positions=[[0,0,0]])
    kwargs = {"kpoints":{"method":"mueller","mindistance":50},
              "potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    assert not calc.can_execute(target)
    calc.write_input(atm,directory=target)

    assert not calc.can_execute("def")
    assert calc.can_execute(target)

def test_can_extract(tmpdir):
    """Tests the can_extract method. We'll also test is_executing at the
    same time.
    """

    from matdb.utility import symlink, relpath

    target = str(tmpdir.join("Vasp"))        
    atm = Atoms("Si",positions=[[0,0,0]])
    kwargs = {"kpoints":{"method":"mueller","mindistance":50},
              "potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    assert not calc.can_extract("def")
    assert not calc.can_extract(target)
    symlink(path.join(calc.folder,"OUTCAR"),
            relpath("tests/files/VASP/OUTCAR_incomplete"))
    assert not calc.can_extract(target)
    assert calc.is_executing(target)
    symlink(path.join(calc.folder,"OUTCAR"),
            relpath("tests/files/VASP/OUTCAR_complete"))
    assert calc.can_extract(target)
    assert not calc.is_executing(target)

def test_extract(tmpdir):
    """Tests the extract method and cleanup method.
    """
    from matdb.utility import symlink, relpath, touch

    target = str(tmpdir.join("Vasp"))        
    atm = Atoms("Si",positions=[[0,0,0]],cell=[1,1,1])
    kwargs = {"ibrion":5, "nsw": 1, "kpoints":{"method":"mueller","mindistance":50},
              "potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)
    
    calc.write_input(atm, target)    
    symlink(path.join(calc.folder,"OUTCAR"),
            relpath("tests/files/VASP/OUTCAR_complete"))
    symlink(path.join(calc.folder,"CONTCAR"),
            path.join(calc.folder,"POSCAR"))
    calc.extract(target)

    assert hasattr(calc.atoms,"vasp_force")
    assert hasattr(calc.atoms,"vasp_stress")
    assert hasattr(calc.atoms,"vasp_energy")
    assert calc.atoms.vasp_energy is not None
    assert calc.atoms.vasp_stress is not None
    assert calc.atoms.vasp_force is not None

    touch(path.join(calc.folder,"CHG"))
    calc.cleanup(target,clean_level="light")
    assert not path.isfile(path.join(calc.folder,"CHG"))
    touch(path.join(calc.folder,"CHGCAR"))
    calc.cleanup(target)
    assert not path.isfile(path.join(calc.folder,"CHGCAR"))
    touch(path.join(calc.folder,"vasprun.xml"))
    calc.cleanup(target,clean_level="aggressive")
    assert not path.isfile(path.join(calc.folder,"OUTCAR"))
    assert not path.isfile(path.join(calc.folder,"vasprun.xml"))

def test_to_dict(tmpdir):
    """Tests the calculator to_dict method.
    """

    target = str(tmpdir.join("Vasp"))        
    atm = Atoms("Si",positions=[[0,0,0]],cell=[1,1,1])
    kwargs = {"ibrion":5, "xc":"pbe", "nsw": 1, "kpoints":{"method":"mueller","mindistance":50},
              "potcars": {"directory":"./tests/vasp", "xc":"pbe"}}
    calc = Vasp(atm, target, str(tmpdir), 0, **kwargs)

    calc_dict = calc.to_dict()

    out = {"folder":target, "ran_seed":0, "contr_dir":str(tmpdir),
           "kwargs": kwargs, "args": (), "version": "vasp.4.6.35"}

    assert compare_nested_dicts(calc_dict,out)

    calc_dict = calc.to_dict()
    assert compare_nested_dicts(calc_dict,out)

def test_phonon_defaults():
    """Tests the default phonon settings.
    """

    from matdb.calculators.vasp import phonon_defaults

    d = {"name":"Vasp"}
    phonon_defaults(d)

    out = {"name":"Vasp", "encut": 500, "ediff": '1.0e-08', "ialgo": 38, "ismear": 0,
           "lreal": False, "addgrid": True, "lwave": False, "lcharg": False,
           "ibrion":-1}
    assert compare_nested_dicts(d,out)

    d = {"name":"Vasp", "encut":200, "ediff":'1.0e-01'}
    phonon_defaults(d, dfpt=True)

    out = {"name":"Vasp", "encut": 200, "ediff": '1.0e-01', "ialgo": 38, "ismear": 0,
           "lreal": False, "addgrid": True, "lwave": False, "lcharg": False,
           "ibrion":8}
    assert compare_nested_dicts(d,out)

def test_extract_force_sets(tmpdir):
    """Tests the extract_force_sets and extract_farce_constants functions.
    """

    from matdb.calculators.vasp import extract_force_sets, extract_force_constants
    from matdb.utility import relpath, symlink

    phonon_dir = str(tmpdir.join("phononpy"))
    mkdir(phonon_dir)
    configs = {"1":str(tmpdir)}

    res = extract_force_sets(configs,phonon_dir)
    assert compare_nested_dicts(res,{"error":""})

    symlink(path.join(str(tmpdir),"vasprun.xml"),
            relpath("tests/files/VASP/vasprun.xml_complete"))

    res = extract_force_sets(configs,phonon_dir)
    assert res["error"] == []

    res = extract_force_constants(configs,phonon_dir)
    assert res["error"] == []
