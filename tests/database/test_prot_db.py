"""Tests the enumeration group interface.
"""
import pytest
from matdb.database.prototype import Prototypes
from matdb.utility import relpath
from os import mkdir, path, symlink, remove
import quippy
import numpy as np
import six

def compare_nested_dicts(dict1,dict2):
    """Compares two dictionaries to see if they are the same.
    """

    print(dict1)
    print(dict2)
    if sorted(dict1.keys()) != sorted(dict2.keys()):
        print("H1", dict1.keys(), dict2.keys())
        return False

    for key in dict1:
        if isinstance(dict1[key],dict):
            res = compare_nested_dicts(dict1[key],dict2[key])
            if not res:
                print("H2", key, dict1[key], dict2[key])
                return False
            else:
                continue
        print("H3",dict1[key], dict2[key], key)
        if isinstance(dict1[key], list) and not dict1[key] == dict2[key]:
            return False
            return False
        elif isinstance(dict1[key],six.string_types) and not dict1[key] == dict2[key]:
            print("H4", key, dict1[key], dict2[key])
            return False

    return True

@pytest.fixture()
def CoNiTi(tmpdir):
    from matdb.utility import relpath
    from matdb.database import Controller
    from os import mkdir, symlink, remove

    target = relpath("./tests/files/CoNiTi")
    dbdir = str(tmpdir.join("coniti_db"))
    mkdir(dbdir)
    
    result = Controller(target, dbdir)
    return result

def test_setup(CoNiTi):
    """Tetsts the setup of the enumerated database.
    """

    assert not CoNiTi.collections['prototype']['prototype'].steps['prototype'].is_setup()

    CoNiTi.setup()

    db = "Prototypes/prototype/ord-1"

    folders = {
        "__files__": ["compute.pkl", "jobfile.sh", "prototype_P_uuid.txt", "puuids.pkl"]}
    for i in range(1,170):
        folders["P.{0}".format(i)] = {"__files__": ["INCAR", "PRECALC", "uuid.txt",
                                                    "POSCAR", "pre_comp_atoms.h5"]}

    dbfolder = path.join(CoNiTi.root,db)
    from matdb.utility import compare_tree
    compare_tree(dbfolder,folders)

    assert CoNiTi.collections['prototype']['prototype'].steps['prototype'].is_setup()

    prot = CoNiTi.collections['prototype']['prototype'].steps['prototype']
    assert len(prot.sequence['ord-1'].puuids) == 169

    assert not prot.ready()

    # We need to create fake atoms.h5 objects so that the system will
    # think that VASP has run and the caculations have been extracted.
    dbfolder = path.join(CoNiTi.root,db)
    for j in range(1,170):
        src = path.join(dbfolder,"P.{}".format(j),"pre_comp_atoms.h5")
        dest = path.join(dbfolder,"P.{}".format(j),"atoms.h5")
        symlink(src,dest)

    assert len(prot.rset) == 169
    assert len(prot.fitting_configs) == 169

    # We run the setup one more time to ensure quick returns
    assert prot.ready()
    prot._setup_configs(False)

def test_to_dict(CoNiTi):
    """Tests the to dict method.
    """
    from matdb import __version__
    
    prot = CoNiTi.collections['prototype']['prototype'].steps['prototype']

    out = prot.to_dict()
    print(out)
    model = {'name': 'prototype', 'trainable': False, 'prefix': 'P',
             'calculator': {'kpoints': {'method': 'mueller', 'mindistance': 40},
                            'pp': 'pbe', 'nsw': 1, 'name': 'Vasp',
                            'potcars': {'directory': './tests/vasp', 'xc': 'PBE'}},
             'permutations': {'ternary': [['Co', 'Ni', 'Ti']]}, 'ran_seed': 10,
             'config_type': None,
             'version': __version__, 'override': {}, 'execution': None, 'root': prot.root,
             'structures': {'unary': 'all',
                         'binary': ['b210_', 'b20_', 'b211_', 'b212_', 'b215_',
                                    'b216_', 'b217_', 'b218_', 'b219_', 'b221_'],
                         'ternary': 10}, 'nconfigs': None}
    # some things we can't know before hand like the python version
    # and the datetime stamp for these entries we simply need to
    # verify that they are present.
    assert "python_version" in out
    assert "datetime" in out
    model["python_version"] = out["python_version"]
    model["datetime"] = out["datetime"]
    assert compare_nested_dicts(out, model)
