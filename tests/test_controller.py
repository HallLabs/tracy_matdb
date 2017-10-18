"""Tests the controller and database collection objects methods
directly. These tests rely on the `./tests/Si` directory, which has the model
outputs that temporary directories will be compared to.
"""
import numpy as np
from os import path
import pytest
from matdb.utility import reporoot
    
def _mimic_vasp(folder, xroot):
    """Copies a `vasprun.xml` and `OUTCAR ` output files from the given folder into
    the execution directory to mimic what VASP would have done.

    Args:
        folder (str): path to the folder where the model files are stored.
        xroot (str): path to the root folder where the config steps are stored.
    """
    from matdb.utility import chdir
    from glob import glob
    from os import path
    from matdb.utility import symlink

    files = ["vasprun.xml", "OUTCAR"]
    
    with chdir(folder):
        for vaspfile in files:
            pattern = vaspfile + "__*"
            for dft in glob(pattern):
                name, config = dft.split("__")
                xpath = path.join(xroot, config, "dynmatrix", "W.1")
                #We want to make some testing assertions to ensure that the
                #stubs ran correctly.
                assert path.isfile(path.join(xpath, "CONTCAR"))
                assert path.isfile(path.join(xpath, ".matdb.module"))
                target = path.join(xpath, name)
                symlink(target, path.join(folder, dft))

@pytest.fixture()
def Pd(tmpdir):
    from matdb.utility import relpath
    from matdb.database.controller import Controller
    from os import mkdir

    target = relpath("./tests/Pd/matdb.yaml")
    dbdir = str(tmpdir.join("pd_db"))
    mkdir(dbdir)
    
    #We need to copy the POSCAR over from the testing directory to the temporary
    #one.
    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")
    copy(POSCAR, dbdir)
    
    return Controller(target, dbdir)

def test_find(Pd):
    """Tests the find function with pattern matching.
    """
    steps = Pd.find("Pd.phonon*.dynmatrix")
    model = ['Pd.phonon-2', 'Pd.phonon-16', 'Pd.phonon-32', 'Pd.phonon-4', 'Pd.phonon-54']
    assert model == [s.parent.name for s in steps]

    steps = Pd.find("Pd.phonon-*2.dynmatrix")
    model = ['Pd.phonon-2', 'Pd.phonon-32']
    assert model == [s.parent.name for s in steps]

def test_repeater_multi(Pd):
    """Tests the `niterations` functionality on simple Pd.
    """
    Pd.setup()

    modelroot = path.join(Pd.root, "Pd.phonon-2", "dynmatrix")
    assert Pd["Pd.phonon-2.dynmatrix"].root == modelroot
    
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
        dbfolder = path.join(Pd.root, db)
        compare_tree(dbfolder, folders)

    #Test the status, we should have some folder ready to execute.
    Pd.status()

    #Test execution command; this uses stubs for all of the commands that would
    #be executed.
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"})
    folder = path.join(reporoot, "tests", "data", "Pd", "recover")
    _mimic_vasp(folder, Pd.root)

    #Now that we have vasp files, we can queue recovery. The recover
    #`vasprun.xml` files that we linked to are not complete for all the
    #structures (on purpose).
    Pd.recover()
    recoveries = ["Pd.phonon-16.dynmatrix",
                  "Pd.phonon-32.dynmatrix",
                  "Pd.phonon-54.dynmatrix"]
    okay = ["Pd.phonon-2.dynmatrix",
            "Pd.phonon-4.dynmatrix"]
    for rkey in recoveries:
        assert path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
    for rkey in okay:
        assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))

    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"}, recovery=True)
    folder = path.join(reporoot, "tests", "data", "Pd", "rexecute")
    _mimic_vasp(folder, Pd.root)   

    #Now that we have recovered vasp files, we can queue a *second*
    #recovery. All but one of the `vasprun.xml` files that we linked to are
    #complete for all the structures. This test that we can recover and execute
    #multiple times in order.
    Pd.recover()
    recoveries = ["Pd.phonon-32.dynmatrix"]
    okay = ["Pd.phonon-2.dynmatrix",
            "Pd.phonon-16.dynmatrix",
            "Pd.phonon-4.dynmatrix",
            "Pd.phonon-54.dynmatrix"]
    for rkey in recoveries:
        assert path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
    for rkey in okay:
        assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))

    #Now copy the final files, which are all good.
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"}, recovery=True)
    folder = path.join(reporoot, "tests", "data", "Pd", "complete")
    _mimic_vasp(folder, Pd.root)

    Pd.recover()
    okay = ["Pd.phonon-2.dynmatrix",
            "Pd.phonon-16.dynmatrix",
            "Pd.phonon-4.dynmatrix",
            "Pd.phonon-54.dynmatrix",
            "Pd.phonon-32.dynmatrix"]
    for rkey in okay:
        assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
    
    #We are finally ready to cleanup the phonon database.
    Pd.cleanup()
    
    for rkey in okay:
        assert path.isfile(path.join(Pd[rkey].root, "phonopy", "FORCE_SETS"))
        assert path.isfile(path.join(Pd[rkey].root, "phonopy", "total_dos.dat"))
    
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
