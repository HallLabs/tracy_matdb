# -*- coding: utf-8 -*-
"""Tests the atoms object and related functions.
"""
import pytest
from os import path, remove, mkdir
import numpy as np
import six

def globals_setup(new_root):
    """Sets up the globals for the calculator instance.
    """
    from matdb.io import read
    from matdb.utility import relpath
    from matdb.calculators.utility import paths, set_paths

    target = relpath("./tests/AgPd/matdb")
    config = path.expanduser(path.abspath(target))
    if path.isabs(config):
        root, config = path.split(config)
    else:
        root, config = path.dirname(config), config
        
    configyml = read(root, config)
    configyml["root"] = new_root

    set_paths(configyml)

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
    from matdb.calculators import Vasp
    from matdb.atoms import Atoms
    from matdb.utility import _set_config_paths

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("to_hdf5"))
    globals_setup(target)
    if not path.isdir(target):
        mkdir(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp",
                                                 "versions":{"Si": '05Jan2001'}}, "xc":"pbe"}    
    potSW = Vasp(atSi, target, str(tmpdir), 0, **kwargs)
    atSi.set_calculator(potSW)
    atSi.add_property("vasp_force", [[-18057.59589857, -18057.59589857, -18057.59589857],
         [ -2997.55626529,  -2997.55626529,  -2997.55626529],
         [  3044.17916471,   3044.17916471, -34130.71583118],
         [ 18969.1757571 ,  18969.1757571 ,  11159.15815145],
         [  3044.17916471, -34130.71583118,   3044.17916471],
         [ 18969.1757571 ,  11159.15815145,  18969.1757571 ],
         [-34130.71583118,   3044.17916471,   3044.17916471],
         [ 11159.15815145,  18969.1757571 ,  18969.1757571 ]])
    atSi.add_param("vasp_energy", 25360.504084423999)
    atSi.add_param("vasp_virial", [[ 33538.34327189,  11045.88697112,  11045.88697112],
         [ 11045.88697112,  33538.34327189,  11045.88697112],
         [ 11045.88697112,  11045.88697112,  33538.34327189]])
    atSi.properties["rand"] = np.random.randint(0, 100, 8)
    atSi.write(target=path.join(target,"temp.h5"))
    atR = Atoms()
    atR.read(target=path.join(target,"temp.h5"))

    # check that the calculator got transfered properly.
    assert type(atSi.calc) == type(atR.calc)
    assert atSi.calc.args == atR.calc.args
    assert atSi.calc.kwargs == atR.calc.kwargs

    # check that the other properties got transfered properly.
    assert atR.vasp_energy == atSi.vasp_energy
    assert isinstance(atR, Atoms)
    assert np.allclose(atR.vasp_force, atSi.vasp_force)
    assert np.allclose(atR.vasp_virial, atSi.vasp_virial)
    assert np.allclose(atR.properties["rand"], atSi.properties["rand"])
    assert np.allclose(atR.positions, atSi.positions)
    remove(path.join(target,"temp.h5"))

def test_remote_read(tmpdir):
    """Tests the reading in of a atoms.h5 file from another directory."""
    from matdb.atoms import Atoms
    from matdb.utility import _set_config_paths, reporoot
    from os import path

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("to_hdf5"))
    globals_setup(target)
    at = Atoms()
    at.read(target=path.join(reporoot, "tests", "files", "test.h5"))
    assert isinstance(at, Atoms)

def test_Atoms_creation(tmpdir):
    """Tests the initialization of the atoms objcet.
    """
    from matdb.atoms import Atoms
    from ase.atoms import Atoms as aseAtoms
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    atR = Atoms(atSi)
    assert atR==atSi
    
    atSi = aseAtoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    atSi.info['nneightol'] = 10
    atSi.info['cutoff'] = 5
    atSi.info['cutoff_break'] = 10
    atR = Atoms(atSi)

    assert np.allclose(atR.positions,atSi.positions)
    assert np.allclose(atR.cell,atSi.cell)
    assert hasattr(atR,'nneightol')
    assert hasattr(atR,'cutoff')
    assert hasattr(atR,'cutoff_break')
    
    from matdb.utility import _set_config_paths
    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("make_atoms"))
    globals_setup(target)
    if not path.isdir(target):
        mkdir(target)
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    atSi.write(target=path.join(target,"temp.xyz"))
    atR = Atoms(path.join(target,"temp.xyz"))
    
    assert np.allclose(atR.positions,atSi.positions)
    assert np.allclose(atR.cell,atSi.cell)
    remove(path.join(target,"temp.xyz"))

def test_Atoms_attributes(tmpdir):
    """Tests the some of the attributes of an atoms object.
    """
    from matdb.calculators import Vasp
    from matdb.atoms import Atoms
    from matdb.utility import _set_config_paths, reporoot
    from os import path

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("atom_attributes"))
    globals_setup(target)
    if not path.isdir(target):
        mkdir(target)

    at1 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp",
                                                 "versions":{"Si": '05Jan2001'}}, "xc":"pbe"}    
    potSW = Vasp(at1, target, str(tmpdir), 0, **kwargs)
    at1.set_calculator(potSW)
    at1.add_property("vasp_force", [[-18057.59589857, -18057.59589857, -18057.59589857],
         [ -2997.55626529,  -2997.55626529,  -2997.55626529],
         [  3044.17916471,   3044.17916471, -34130.71583118],
         [ 18969.1757571 ,  18969.1757571 ,  11159.15815145],
         [  3044.17916471, -34130.71583118,   3044.17916471],
         [ 18969.1757571 ,  11159.15815145,  18969.1757571 ],
         [-34130.71583118,   3044.17916471,   3044.17916471],
         [ 11159.15815145,  18969.1757571 ,  18969.1757571 ]])
    at1.add_param("vasp_energy", 1234)

    assert at1.get_energy() == 1234

    at1.rm_param("vasp_energy")
    assert not "vasp_energy" in at1.info["params"]
    at1.rm_property("vasp_force")
    assert not "vasp_force" in at1.info["properties"]

    at1.params = {"vasp_energy": 1234}
    assert at1.vasp_energy == 1234

    at1.properties["vasp_force"] = [[-18057.59589857, -18057.59589857, -18057.59589857],
         [ -2997.55626529,  -2997.55626529,  -2997.55626529],
         [  3044.17916471,   3044.17916471, -34130.71583118],
         [ 18969.1757571 ,  18969.1757571 ,  11159.15815145],
         [  3044.17916471, -34130.71583118,   3044.17916471],
         [ 18969.1757571 ,  11159.15815145,  18969.1757571 ],
         [-34130.71583118,   3044.17916471,   3044.17916471],
         [ 11159.15815145,  18969.1757571 ,  18969.1757571 ]]
    assert "vasp_force" in at1.info["properties"]
    
def test_AtomsList_creation(tmpdir):
    """Tests the creation of the AtomsList object. 
    """
    from matdb.calculators import Vasp
    from matdb.atoms import Atoms, AtomsList
    from matdb.utility import _set_config_paths

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("make_AtomsList"))
    globals_setup(target)

    if not path.isdir(target):
        mkdir(target)
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp",
                                                 "versions":{"Si": '05Jan2001'}}, "xc":"pbe"}    

        
    at1 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43],info={"rand":10})
    potSW = Vasp(at1, target, str(tmpdir), 0, **kwargs)
    at1.set_calculator(potSW)
    at1.add_property("vasp_force", [[-18057.59589857, -18057.59589857, -18057.59589857],
         [ -2997.55626529,  -2997.55626529,  -2997.55626529],
         [  3044.17916471,   3044.17916471, -34130.71583118],
         [ 18969.1757571 ,  18969.1757571 ,  11159.15815145],
         [  3044.17916471, -34130.71583118,   3044.17916471],
         [ 18969.1757571 ,  11159.15815145,  18969.1757571 ],
         [-34130.71583118,   3044.17916471,   3044.17916471],
         [ 11159.15815145,  18969.1757571 ,  18969.1757571 ]])
    at1.add_param("vasp_energy", 1234)
    
    at2 = Atoms("S6",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75]],
                 cell=[6.43,5.43,4.43],info={"rand":10},calculator=potSW)
    kwargs = {"kpoints":{"rmin":50}, "potcars": {"directory":"./tests/vasp",
                                                 "versions":{"S": '06Sep2000'}}, "xc":"pbe"}    
    potSW = Vasp(at2, target, str(tmpdir), 0, **kwargs)
    at2.set_calculator(potSW)
    at2.add_property("vasp_force", np.random.randint(0, 100, (6,3)))
    at2.add_param("vasp_energy", 4321)
    
    al1 = AtomsList([at1,at2])
    
    assert len(al1) == 2

    at1.write(target=path.join(target,"temp1.h5"))
    at2.write(target=path.join(target,"temp2.h5"))

    al2 = AtomsList([path.join(target,"temp1.h5"),path.join(target,"temp2.h5")])
    assert len(al2) == 2
    assert isinstance(al2[0],Atoms)

    empty_list = AtomsList([])
    assert len(empty_list) == 0

    al3 = AtomsList(at1)
    assert al3[0] == at1

    al4 = AtomsList(path.join(target,"temp1.h5"))
    assert len(al4) == 1
    assert isinstance(al4[0],Atoms)
    assert al4[0].vasp_energy == at1.vasp_energy

def test_AtomsList_attributes():
    """Tests the atoms lists attributes.
    """
    from matdb.atoms import Atoms, AtomsList

    at1 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43],info={"rand":10})
    
    at2 = Atoms("S6",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75]],
                 cell=[6.43,5.43,4.43],info={"rand":10})

    at3 = Atoms("CNi",positions=[[0,0,0],[0.5,0.5,0.5]], info={"rand":8})
    at4 = Atoms(info={"rand":7})
    at4.copy_from(at3)
    
    al1 = AtomsList([at1,at2,at3,at4])

    alpos = al1.positions
    assert np.allclose(alpos[0],at1.positions)
    assert np.allclose(alpos[1],at2.positions)
    assert np.allclose(alpos[2],at3.positions)
    assert np.allclose(alpos[3],at4.positions)

    with pytest.raises(AttributeError):
        al1.__getattr__('__dict__')
    assert al1.energy is None

    alslice = al1[0:2]
    assert len(alslice) == 2
    assert alslice[0] == at1
    assert alslice[1] == at2

    alitems = al1[[0,1,3]]
    assert len(alitems) == 3
    assert alitems[0] == at1
    assert alitems[1] == at2
    assert alitems[2] == at4

    with pytest.raises(IndexError):
        al1[[0,2.3]]

    alitems = al1[[True,True,False,True]]
    assert len(alitems) == 3

    for i in al1.iterframes():
        assert isinstance(i,Atoms)
    
    for i in al1.iterframes(reverse=True):
        assert isinstance(i,Atoms)

    assert al1.random_access

    def get_pos(atoms):
        return atoms.positions

    alpos=al1.apply(get_pos)
    assert np.allclose(alpos[0],at1.positions)
    assert np.allclose(alpos[1],at2.positions)
    assert np.allclose(alpos[2],at3.positions)
    assert np.allclose(alpos[3],at4.positions)

    al1.sort(reverse=True)
    al1.sort(attr='rand')

    with pytest.raises(ValueError):
        al1.sort(attr='positions',key=2)

def test_AtomsList_io(tmpdir):
    """Tests the AtomsList writing and reading from file.
    """
    from matdb.atoms import Atoms, AtomsList

    at1 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43],info={"rand":10})
    
    at2 = Atoms("S6",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75]],
                 cell=[6.43,5.43,4.43],info={"rand":10})

    at3 = Atoms("CNi",positions=[[0,0,0],[0.5,0.5,0.5]])
    at4 = Atoms()
    at4.copy_from(at3)
    
    al1 = AtomsList([at1,at2,at3,at4])
    
    target = str(tmpdir.join("atomList_to_hdf5"))
    if not path.isdir(target):
        mkdir(target)

    al1.write(path.join(target,"temp.h5"))

    aR = AtomsList()
    aR.read(path.join(target,"temp.h5"))

    assert len(aR) == len(al1)

    alpos = aR.positions
    assert any([np.allclose(alpos[i],at1.positions) for i in range(4) if
                len(alpos[i])==len(at1.positions)])
    assert any([np.allclose(alpos[i],at2.positions) for i in range(4) if
                len(alpos[i])==len(at2.positions)])
    assert any([np.allclose(alpos[i],at3.positions) for i in range(4) if
                len(alpos[i])==len(at3.positions)])
    assert any([np.allclose(alpos[i],at4.positions) for i in range(4) if
                len(alpos[i])==len(at4.positions)])

    al1.write(path.join(target,"temp.xyz"))

    aR = AtomsList()
    aR.read(path.join(target,"temp.xyz"))

    assert len(aR) == len(al1)

    aR.read(path.join(target,"temp.xyz"))
    assert len(aR) == 2*len(al1)


    # Test reading in of a single atoms object.

    aR1 = Atoms(path.join(target,"temp.h5"))
    assert isinstance(aR1,Atoms)
    assert any([np.allclose(alpos[i],at1.positions) for i in range(4) if
                len(alpos[i])==len(at1.positions)])

def test_ase_atoms_conversion(tmpdir):
    """Tests the conversion of an ase atoms objcet to a 'matdb.atoms.Atoms' object. 
    """

    from matdb.atoms import Atoms as matAtoms
    from ase.atoms import Atoms
    from matdb.utility import _set_config_paths

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("from_ase"))
    globals_setup(target)
    if not path.isdir(target):
        mkdir(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    aR = matAtoms(atSi)

    assert np.allclose(aR.positions, atSi.positions)
    assert aR.calc == atSi.calc 
    
def test_to_dict(tmpdir):
    """Tests the conversion of atoms to dictionaries.
    """
    
    from matdb.calculators import Vasp
    from matdb.atoms import Atoms as Atoms
    from matdb.utility import _set_config_paths

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("atoms_dict"))
    globals_setup(target)

    if not path.isdir(target):
        mkdir(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    kwargs = {"encut":400, "kpoints": {"rmin": 50},
              "potcars":{"xc": "pbe", "directory": "./tests/vasp", "versions": {"Si": "05Jan2001"}}}
    
    calc = Vasp(atSi, target, '.', 0, **kwargs)

    atSi.set_calculator(calc)
    atSi.group_uuid = "123456"
    Sidict = atSi.to_dict()

    
    assert "calc" in Sidict
    assert "calc_kwargs" in Sidict
    assert Sidict["calc_kwargs"]["encut"] == 400
    assert Sidict["group_uuid"] == "123456"
    assert "potcars" in Sidict["calc_kwargs"]
    assert "kpoints" in Sidict["calc_kwargs"]

def test_read_atoms(tmpdir):
    """Tests the reading of atoms objects from files.
    """
    
    from matdb.calculators import Vasp
    from matdb.atoms import Atoms as Atoms
    from matdb.utility import _set_config_paths

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("read_atoms"))
    globals_setup(target)

    if not path.isdir(target):
        mkdir(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    kwargs = {"encut":400, "kpoints": {"rmin": 50},
              "potcars":{"xc": "pbe", "directory": "./tests/vasp", "versions": {"Si": "05Jan2001"}}}

    calc = Vasp(atSi, target, '.', 0, **kwargs)
    atSi.set_calculator(calc)

    atSi.set_calculator(calc)
    atSi.add_property("vasp_force", [[-18057.59589857, -18057.59589857, -18057.59589857],
         [ -2997.55626529,  -2997.55626529,  -2997.55626529],
         [  3044.17916471,   3044.17916471, -34130.71583118],
         [ 18969.1757571 ,  18969.1757571 ,  11159.15815145],
         [  3044.17916471, -34130.71583118,   3044.17916471],
         [ 18969.1757571 ,  11159.15815145,  18969.1757571 ],
         [-34130.71583118,   3044.17916471,   3044.17916471],
         [ 11159.15815145,  18969.1757571 ,  18969.1757571 ]])
    atSi.add_param("vasp_energy", 25360.504084423999)
    atSi.add_param("vasp_virial", [[ 33538.34327189,  11045.88697112,  11045.88697112],
         [ 11045.88697112,  33538.34327189,  11045.88697112],
         [ 11045.88697112,  11045.88697112,  33538.34327189]])
    atSi.group_uuid = "123456"

    temp = path.join(target,"temp.h5")
    atSi.write(temp)

    atR = Atoms(temp)
    
    assert atR.calc.name == "Vasp"
    assert hasattr(atR.calc,"potcars")
    assert atR.calc.kwargs["encut"] == 400
    assert np.allclose(atR.positions,[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]])
    
def test_reading_multiple_files(tmpdir):
    """Tests the reading in of multiple atoms objects to an AtomsList.
    """
    
    from matdb.calculators import Vasp
    from matdb.atoms import Atoms as Atoms, AtomsList
    from matdb.io import save_dict_to_h5
    import h5py
    from matdb.utility import _set_config_paths

    _set_config_paths("AgPd_Enumerated", str(tmpdir))
    target = str(tmpdir.join("read_atoms2"))
    globals_setup(target)

    if not path.isdir(target):
        mkdir(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    kwargs = {"encut":400, "kpoints": {"rmin": 50},
              "potcars":{"xc": "pbe", "directory": "./tests/vasp", "versions": {"Si": "05Jan2001"}}}
    
    calc = Vasp(atSi, target, '.', 0, **kwargs)
    atSi.set_calculator(calc)

    temp = path.join(target,"temp.h5")
    atSi.write(temp)    

    atSi2 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[6.43,6.43,6.43])

    kwargs = {"encut":600, "kpoints": {"rmin": 50},
              "potcars":{"xc": "pbe", "directory": "./tests/vasp", "versions": {"Si": "05Jan2001"}}}
    
    calc = Vasp(atSi2, target, '.', 0, **kwargs)
    atSi2.set_calculator(calc)
    temp2 = path.join(target,"temp2.h5")
    atSi2.write(temp2)

    atRL = AtomsList([temp,temp2])

    assert len(atRL) == 2
    assert atRL[0].calc.kwargs["encut"] != atRL[1].calc.kwargs["encut"]
    assert atRL[1].calc.kwargs["encut"] in [400,600]
    assert atRL[0].calc.kwargs["encut"] in [400,600]

    atom_dict = {"atom_1":temp, "atom_2": temp2}

    temp3 = path.join(target,"temp3.h5")
    with h5py.File(temp3,"w") as hf:
        save_dict_to_h5(hf,atom_dict,'/')

    atRL = AtomsList(temp3)

    assert len(atRL) == 2
    assert atRL[0].calc.kwargs["encut"] != atRL[1].calc.kwargs["encut"]
    assert atRL[1].calc.kwargs["encut"] in [400,600]
    assert atRL[0].calc.kwargs["encut"] in [400,600]
