"""Tests the vasp calculator in matdb.
"""
import pytest
from os import path, mkdir, remove
import six
import numpy as np

from matdb.atoms import Atoms
from matdb.utility import reporoot, relpath, symlink, touch
from matdb.calculators import Qe
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
    
    mod = get_calculator_module({"name":"Qe"})
    assert isinstance(mod, ModuleType)

def test_qe_setup(tmpdir):
    """Tests Vasp calculator initialization.
    """

    target = str(tmpdir.join("Qe"))        
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])

    kwargs = {"potcars": {"directory":"./tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "MP", "divisions": (3, 3, 3)},
              "input_data": {"control":{"calculation": "relax", "prefix": "test"}}}

    calc = Qe(atm, target, 'def', 0, **kwargs)

    assert isinstance(calc, Qe)
    assert calc.parameters["input_data"]["control"]["calculation"] == "relax"
    assert "kpts" in calc.parameters
    assert calc.out_file == "test"

    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1},
              "input_data": {"control":{"calculation": "relax", "prefix":"test",
                                        "outdir": "temp"}}}

    calc = Qe(atm, target, '.', 0, **kwargs)
    calc.write_input(atm)

    assert calc.parameters["koffset"] == 1
    assert "kspacing" in calc.parameters
    assert calc.out_file == "temp/test"
    
    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1}}

    calc = Qe(atm, target, '.', 0, **kwargs)
    calc.write_input(atm)
    assert "input_data" not in calc.in_kwargs
    

def test_chekc_potcar(tmpdir):
    """Tests the checking of the potential files.
    """
    target = str(tmpdir.join("Qe"))        
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])
    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1},
              "input_data": {"control":{"calculation": "relax"}},
              "output": "test.xml"}

    calc = Qe(atm, target, '.', 0, **kwargs)

    calc.potcars["versions"]["Al"][1] = ".4.1"
    with pytest.raises(VersionError):
        calc._check_potcars()
    
    calc.potcars["versions"]["Al"][0] = "1.4.1"
    with pytest.raises(VersionError):
        calc._check_potcars()

    del calc.potcars["directory"]
    with pytest.raises(IOError):
        calc._check_potcars()

def test_write_input_can_execute(tmpdir):
    """Tests the writing of the input files and the tests of the can_execute method. 
    """
    target = str(tmpdir.join("Qe"))     
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])
    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1},
              "input_data": {"control":{"calculation": "relax"}},
              "output": "test.xml"}
    calc = Qe(atm, target, '.', 0, **kwargs)

    assert not calc.can_execute(target)

    calc.create()
    assert path.isfile(path.join(target, "espresso.pwi"))
    assert calc.can_execute(target)
    
    calc.folder = path.join(target,"Qe_2")
    calc.write_input(atm)
    assert path.isfile(path.join(calc.folder, "espresso.pwi"))

    assert not calc.can_execute(path.join(target,"Qe_3"))

def test_can_extract(tmpdir):
    """Tests the can_extract method. We'll also test is_executing at the
    same time.
    """

    target = str(tmpdir.join("Qe"))     
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])
    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1},
              "input_data": {"control":{"calculation": "relax", "prefix": "test"}}}
    calc = Qe(atm, target, '.', 0, **kwargs)
    calc.create()

    assert not calc.can_extract("def")
    assert not calc.can_extract(target)
    assert not calc.is_executing(target)

    touch(path.join(calc.folder, "CRASH"))
    assert not calc.can_extract(target)
    remove(path.join(calc.folder, "CRASH"))
    
    symlink(path.join(calc.folder,"pwscf.xml"),
            relpath("tests/qe/complete.xml"))
    mkdir(path.join(calc.folder, "pwscf.save"))
    
    assert calc.can_extract(target)
    remove(path.join(calc.folder, "test.xml"))
    mkdir(path.join(calc.folder, "pwscf.save"))
    
    symlink(path.join(calc.folder,"pwscf.xml"),
            relpath("tests/qe/fail.xml"))
    assert not calc.can_extract(target)
    
def test_extract(tmpdir):
    """Tests the extract method and cleanup method.
    """
    target = str(tmpdir.join("Qe"))     
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])
    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1},
              "input_data": {"control":{"calculation": "relax", "prefix": "test"}}}
    calc = Qe(atm, target, '.', 0, **kwargs)
    calc.create()

    symlink(path.join(calc.folder,"test.xml"),
            relpath("tests/qe/complete.xml"))
    mkdir(path.join(calc.folder, "test.save"))
    calc.extract(target)

    assert hasattr(calc.atoms,"qe_force")
    assert hasattr(calc.atoms,"qe_stress")
    assert hasattr(calc.atoms,"qe_energy")
    assert calc.atoms.qe_energy is not None
    assert calc.atoms.qe_stress is not None
    assert calc.atoms.qe_force is not None

    touch(path.join(calc.folder, "test.save", "paw.txt"))
    calc.cleanup(target,clean_level="light")
    assert not path.isfile(path.join(calc.folder,"test.save", "paw.txt"))
    touch(path.join(calc.folder, "test.save", "charge-density.dat"))
    calc.cleanup(target)
    assert not path.isfile(path.join(calc.folder, "test.save", "charge-density.dat"))
    calc.cleanup(target,clean_level="aggressive")
    assert not path.isfile(path.join(calc.folder,"test.xml"))
    assert not path.isdir(path.join(calc.folder,"test.save"))

def test_to_dict(tmpdir):
    """Tests the calculator to_dict method.
    """

    target = str(tmpdir.join("Qe"))     
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])
    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1},
              "input_data": {"control":{"calculation": "relax", "prefix": "test"}}}
    calc = Qe(atm, target, '.', 0, **kwargs)
    symlink(path.join(calc.folder,"test.xml"),
            relpath("tests/qe/complete.xml"))
    calc.extract(target)

    calc_dict = calc.to_dict()

    out = {"folder":target, "ran_seed":0, "contr_dir": '.',
           "kwargs": kwargs, "args": (), "version": "6.2 (svn rev. 14038)"}

    assert compare_nested_dicts(calc_dict,out)

def test_read(tmpdir):
    """Tests the read function of the QE calculator.
    """
    
    target = str(tmpdir.join("Qe"))     
    atm = Atoms("AlPd",positions=[[0, 0, 0],[0.5, 0.5, 0.5]])
    kwargs = {"potcars": {"directory":"~/codes/matdb/tests/qe",
                          "potentials": {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
                                         "Pd": "Pd_ONCV_PBE-1.0.upf"},
                          "versions": {"Al": ["2.0.1", ".5.1"], "Pd": ["2.0.1", "2.1.1"]}},
              "kpoints":{"method": "kspacing", "spacing": 0.1, "offset": 1},
              "input_data": {"control":{"calculation": "relax", "prefix": "test"}}}
    calc = Qe(atm, target, '.', 0, **kwargs)
    symlink(path.join(calc.folder,"test.xml"),
            relpath("tests/qe/complete.xml"))

    output = calc._read(path.join(calc.folder,"test.xml"))

    assert output["convergence"] == 4.068079462655824e-7
    assert np.allclose(output["atoms"], [0,0,0])
    assert np.allclose(output["cell"], [[-3.75,0,3.75],[0,3.75,3.75],[-3.75,3.75,0]])
    assert output["etot"] == -1.975055613913407e1
    assert np.allclose(output["forces"], [0,0,0])
    assert np.allclose(output["stress"], [[1.578434139006113e-4, -1.219727444046192e-19,
                                           -9.486769009248164e-20],
                                          [-1.490777987167569e-19, 1.578434139006113e-4,
                                           9.486769009248164e-20],
                                          [-6.776263578034403e-20, 1.219727444046192e-19,
                                           1.578434139006113e-4]])

    assert calc.version == '6.2 (svn rev. 14038)'
