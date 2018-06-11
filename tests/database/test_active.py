"""Tests the enumeration group interface.
"""
import pytest
from os import mkdir, path, symlink, remove
import numpy as np
import six

from matdb.database.active import Active
from matdb.utility import relpath, compare_tree
from matdb.fitting.controller import TController
from matdb.database import Database, Controller
from matdb.atoms import Atoms

@pytest.fixture()
def Act(tmpdir):

    target = relpath("./tests/AgPd/matdb")
    dbdir = str(tmpdir.join("active_db"))
    mkdir(dbdir)

    cntrl = Controller(target, dbdir)
    db = Database("active", dbdir, cntrl, [{"type":"active.Active"}], {}, 0)
    tcntrl = TController(db=db, root=dbdir, fits={})
    dbargs = {"root": dbdir, "parent": db,
              "calculator": tcntrl.db.calculator}
    result = Active(**dbargs)
    return result

def test_all_active(Act):
    """Tetsts the setup of the Active database.
    """
    assert Act.last_iteration is None
    assert Act.last_config_atoms is None

    configs = []
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    configs.append(atSi)
    atSi = Atoms("Si6",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75]],
                 cell=[5.5, 5.5, 5.5])
    configs.append(atSi)
    atSi = Atoms("Si4",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25]],
                 cell=[5.4,5.4,5.4])
    configs.append(atSi)

    Act.add_configs(configs, 1)

    assert Act.iter_file == path.join(Act.root, "iter_1.pkl")
    assert Act.nconfigs == 3
    assert not Act.is_executing()

    Act.setup()

    assert len(Act.last_iteration) == 3
    
    folders = {
        "__files__": ["compute.pkl", "auids.pkl", "jobfile.sh", "index.json", "iter_1.pkl",
                      "active_Ac_uuid.txt"],
        "Ac.1": {
            "__files__": ["INCAR", "PRECALC", "POSCAR", "POTCAR", "pre_comp_atoms.h5",
                          "uuid.txt"]
        },
        "Ac.2": {
            "__files__": ["INCAR", "PRECALC", "POSCAR", "POTCAR", "pre_comp_atoms.h5",
                          "uuid.txt"]
        },
        "Ac.3": {
            "__files__": ["INCAR", "PRECALC", "POSCAR", "POTCAR", "pre_comp_atoms.h5",
                          "uuid.txt"]
        }
    }

    dbfolder = Act.root
    compare_tree(dbfolder,folders)    

    assert Act.is_setup()

    assert not Act.ready()
    # We need to fake some VASP output so that we can cleanup the
    # database and get the rset

    src = relpath("./tests/data/Pd/complete/OUTCAR__DynMatrix_phonon_Pd_dim-2.00")
    dbfolder = Act.root
    for j in range(1,4):
        dest = path.join(dbfolder,"Ac.{}".format(j),"OUTCAR")
        symlink(src,dest)
            
    dbfolder = Act.root
    for j in range(1,4):
        src = path.join(dbfolder,"Ac.{}".format(j),"POSCAR")
        dest = path.join(dbfolder,"Ac.{}".format(j),"CONTCAR")
        symlink(src,dest)

    remove(path.join(Act.root, "Ac.1", "pre_comp_atoms.h5"))
    Act.extract()
    assert len(Act.config_atoms) == 3
    assert len(Act.configs) == 3
    assert len(Act.last_config_atoms) == 3
    assert len(Act.rset) == 3
    assert not Act.is_executing()

    remove(path.join(Act.root, "Ac.1", "OUTCAR"))
    dest = path.join(Act.root, "Ac.1", "OUTCAR")
    src = relpath("./tests/data/Pd/basic_fail/S.4/OUTCAR")
    symlink(src,dest)
    assert Act.is_executing()

    # test the addition of a second set of new configs

    # first config is to test if the exclusion of duplicates works.
    configs = []
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])
    configs.append(atSi)

    
    atSi = Atoms("Si",positions=[[0,0,0]],
                 cell=[3,3,3])
    configs.append(atSi)

    atSi = Atoms("Si2",positions=[[0,0,0], [0.5,0.75,0.25]],
                 cell=[3.43,3.43,3.43])
    configs.append(atSi)

    Act.add_configs(configs, 2)

    assert Act.iter_file == path.join(Act.root, "iter_2.pkl")
    assert Act.nconfigs == 6
    assert not Act.is_executing()

    Act.setup()
    
    assert Act.nconfigs == 5    
    assert len(Act.last_config_atoms) == 2

    Act.iter_file = path.join(Act.root, "iter_1.pkl")
    Act._load_last_iter()
    assert len(Act.last_iteration) == 3
