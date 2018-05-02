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
        print(dict1[key], dict2[key], key)
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

    db = "Prototypes/prototype/"

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
    assert len(prot.puuids) == 169

    assert not prot.ready()
    # We need to fake some VASP output so that we can cleanup the
    # database and get the rset

    src = relpath("./tests/data/Pd/complete/OUTCAR__DynMatrix_phonon_Pd_dim-2.00")
    dbfolder = path.join(CoNiTi.root,db)
    for j in range(1,170):
        dest = path.join(dbfolder,"P.{}".format(j),"OUTCAR")
        symlink(src,dest)
            
        src = path.join(dbfolder,"P.{}".format(j),"POSCAR")
        dest = path.join(dbfolder,"P.{}".format(j),"CONTCAR")
        symlink(src,dest)

    prot.extract()
    assert len(prot.rset) == 169
    assert len(prot.fitting_configs) == 169

def test_to_dict(CoNiTi):
    """Tests the to dict method.
    """
    from matdb import __version__
    
    prot = CoNiTi.collections['prototype']['prototype'].steps['prototype']

    out = prot.to_dict()
    print(out)
    model = {'name': 'prototype', 'calculator': None, 'trainable': False, 'prefix': 'P',
     'order': {'ternary': [['Co', 'Ni', 'Ti']]}, 'ran_seed': 10, 'config_type': None,
     'version': __version__, 'override': {}, 'execution': None, 'root': prot.root,
     'structs': {'unary': 'all', 'binary': 10, 'ternary': 10}, 'nconfigs': None}
    # some things we can't know before hand like the python version
    # and the datetime stamp for these entries we simply need to
    # verify that they are present.
    assert "python_version" in out
    assert "datetime" in out
    model["python_version"] = out["python_version"]
    model["datetime"] = out["datetime"]
    assert compare_nested_dicts(out, model)
