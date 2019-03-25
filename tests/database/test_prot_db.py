"""Tests the enumeration group interface.
"""
import pytest
from matdb.database.prototype import Prototypes
from matdb.utility import relpath
from os import mkdir, path, symlink, remove
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
        if isinstance(dict1[key], list) and not dict1[key] == dict2[key]:
            return False
            return False
        elif isinstance(dict1[key],six.string_types) and not dict1[key] == dict2[key]:
            return False

    return True

@pytest.fixture()
def CoNiTi(tmpdir):
    from matdb.utility import relpath, copyonce
    from matdb.database import Controller
    from os import mkdir, symlink, remove

    target = relpath("./tests/files/CoNiTi.yml")
    dbdir = str(tmpdir.join("coniti_db"))
    mkdir(dbdir)
    copyonce(target, path.join(dbdir, "matdb.yml"))
    target = path.join(dbdir,"matdb")
    
    result = Controller(target, dbdir)
    return result

def test_setup(CoNiTi):
    """Tetsts the setup of the enumerated database.
    """

    assert not CoNiTi.collections['prototype'].steps['prototype'].is_setup()

    CoNiTi.setup()

    db = "Prototypes/prototype.prototype/per-1"

    folders = {
        "__files__": ["compute.pkl", "jobfile.sh", "prototype_P_uuid.txt", "puuids.pkl"]}
    for i in range(1,251):
        folders["P.{0}".format(i)] = {"__files__": ["INCAR", "PRECALC", "uuid.txt",
                                                    "POSCAR", "pre_comp_atoms.h5"]}

    dbfolder = path.join(CoNiTi.root,db)
    from matdb.utility import compare_tree
    compare_tree(dbfolder,folders)

    assert CoNiTi.collections['prototype'].steps['prototype'].is_setup()

    prot = CoNiTi.collections['prototype'].steps['prototype']
    assert len(prot.sequence['per-1'].puuids) == 250

    assert not prot.ready()

    # We need to create fake atoms.h5 objects so that the system will
    # think that VASP has run and the caculations have been extracted.
    dbfolder = path.join(CoNiTi.root,db)
    for j in range(1,251):
        src = path.join(dbfolder,"P.{}".format(j),"pre_comp_atoms.h5")
        dest = path.join(dbfolder,"P.{}".format(j),"atoms.h5")
        symlink(src,dest)

    assert len(prot.rset) == 250
    assert len(prot.fitting_configs) == 250

    # We run the setup one more time to ensure quick returns
    assert prot.ready()
    prot._setup_configs(False)

def test_to_dict(CoNiTi):
    """Tests the to dict method.
    """
    from matdb import __version__
    
    prot = CoNiTi.collections['prototype'].steps['prototype']

    out = prot.to_dict()
    model = {'name': 'prototype', 'trainable': False, 'prefix': 'P',
             'calculator': {'kpoints': {'method': 'mueller', 'mindistance': 40},
                            'pp': 'pbe', 'nsw': 1, 'name': 'Vasp',
                            'potcars': {'directory': './tests/vasp', 'xc': 'PBE',
                                        'versions': {'Co': '02Aug2007',
                                                     'Ni': '02Aug2007',
                                                     'Ti': '08Apr2002'}}},
             'permutations': {'ternary': [['Co', 'Ni', 'Ti']]}, 'ran_seed': 10,
             'config_type': None,
             'version': __version__, 'execution': None, 'root': prot.root,
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
