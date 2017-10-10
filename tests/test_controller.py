"""Tests the controller and database collection objects methods
directly. These tests rely on the `./tests/Si` directory, which has the model
outputs that temporary directories will be compared to.
"""
import numpy as np
from os import path
def test_repeater_multi():
    """Tests the `niterations` functionality on simple Pd.
    """
    from matdb.utility import relpath
    from matdb.database.controller import Controller
    target = relpath("./tests/Pd/matdb.yaml")
    control = Controller(target)
    control.setup()

    #The matdb.yml file specifies the following databases:
    dbs = ["Pd.phonon-{}".format(i) for i in (2, 4, 16, 32, 54)]
    #Each one should have a folder for: ["dynmatrix", "modulations"]
    #On the first go, the modulations folder will be empty because the DFT
    #calculations haven't been performed yet. However, dynmatrix should have DFT
    #folders ready to go.
    folders = {
        "dynmatrix": {
            "__files__": ["INCAR", "PRECALC"],
            "phonopy": {
                "__files__": ["POSCAR", "POSCAR-001", "disp.yaml", "phonopy_disp.yaml"]
            },
            "phoncache": {},
            "W.1": {
                "__files__": ["INCAR", "POSCAR", "POTCAR", "PRECALC", "KPOINTS"]
            }
        }
    }

    from matdb.utility import compare_tree
    for db in dbs:
        dbfolder = path.join(control.root, db)
        compare_tree(dbfolder, folders)
    
# def test_split():
#     """Tests the splitting logic and that the ids from original
#     randomization are saved correctly.
#     """
#     from cPickle import load
#     with open("PdAg50/ids.pkl", 'rb') as f:
#         d = load(f)
#     assert d["Nsuper"]+d["Nhold"]+d["Ntrain"] == d["Ntot"]
#     assert np.round(d["Nsuper"]/float(d["Nhold"]), 2) == 0.33
#     assert np.round((d["Nsuper"] + d["Nhold"])/float(d["Ntrain"])) == 0.33
