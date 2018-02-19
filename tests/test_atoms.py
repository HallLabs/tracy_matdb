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

    # check that the calculator got transfered properly.
    assert type(atSi.calc) == type(atR.calc)
    assert atSi.calc.args == atR.calc.args
    assert atSi.calc.kwargs == atR.calc.kwargs

    # check that the other properties got transfered properly.
    assert atR.energy == atSi.energy
    assert isinstance(atR, Atoms)
    assert np.allclose(atR.force, atSi.force)
    assert np.allclose(atR.virial, atSi.virial)
    assert np.allclose(atR.properties["rand"], atSi.properties["rand"])
    assert np.allclose(atR.positions, atSi.positions)
    remove(path.join(target,"temp.h5"))

def test_Atoms_creation(tmpdir):
    """Tests the initialization of the atoms objcet.
    """
    from matdb.calculators import Quip
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
    
    target = str(tmpdir.join("make_atoms"))
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

    potSW = Quip(atSi, target, ["IP SW"])
    potSW.results["energy"] = 1234
    potSW.results["force"] = np.random.randint(0, 100, (8,3))
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43],info={"rand":10},calculator=potSW)

    assert "rand" in atSi.params
    assert atSi.params['energy'] == potSW.results['energy']
    assert np.allclose(atSi.properties['force'], potSW.results['force'])

    atSi.add_property('force',np.random.randint(0, 100, (8,3)))
    assert not np.allclose(atSi.properties['force'], potSW.results['force'])

def test_AtomsList_creation(tmpdir):
    """Tests the creation of the AtomsList object. 
    """
    from matdb.calculators import Quip
    from matdb.atoms import Atoms, AtomsList

    target = str(tmpdir.join("make_AtomsList"))
    if not path.isdir(target):
        mkdir(target)

    at1 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43],info={"rand":10})
    potSW = Quip(at1, target, ["IP SW"])
    potSW.results["energy"] = 1234
    potSW.results["force"] = np.random.randint(0, 100, (8,3))
    at1 = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43],info={"rand":10},calculator=potSW)
    
    at2 = Atoms("S6",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75]],
                 cell=[6.43,5.43,4.43],info={"rand":10},calculator=potSW)
    potSW = Quip(at2, target, ["IP SW"])
    potSW.results["energy"] = 4321
    potSW.results["force"] = np.random.randint(0, 100, (6,3))
    at2 = Atoms("S6",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75]],
                 cell=[6.43,5.43,4.43],info={"rand":10},calculator=potSW)
    
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
    assert al4[0].energy == at1.energy

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

    at3 = Atoms("CNi",positions=[[0,0,0],[0.5,0.5,0.5]])
    at4 = Atoms()
    at4.copy_from(at3)
    
    al1 = AtomsList([at1,at2,at3,at4])

    alpos = al1.positions
    assert np.allclose(alpos[0],at1.positions)
    assert np.allclose(alpos[1],at2.positions)
    assert np.allclose(alpos[2],at3.positions)
    assert np.allclose(alpos[3],at4.positions)

    with pytest.raises(AttributeError):
        al1.__getattr__('__dict__')
    with pytest.raises(AttributeError):
        al1.energy

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
    al1.sort(attr='calc')

    with pytest.raises(ValueError):
        al1.sort(attr='positions',cmp=2)

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
    assert np.allclose(alpos[0],at1.positions)
    assert np.allclose(alpos[1],at2.positions)
    assert np.allclose(alpos[2],at3.positions)
    assert np.allclose(alpos[3],at4.positions)

    al1.write(path.join(target,"temp.xyz"))

    aR = AtomsList()
    aR.read(path.join(target,"temp.xyz"))

    assert len(aR) == len(al1)
