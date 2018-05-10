"""Tests the controller and database collection objects methods
directly. These tests rely on the `./tests/Pd` directory, which has the model
outputs that temporary directories will be compared to.
"""
import numpy as np
from os import path
import pytest
from matdb.utility import reporoot, relpath
import six

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
                xpath = path.join(xroot, path.join(*config.split("_")), "W.1")
                #We want to make some testing assertions to ensure that the
                #stubs ran correctly.
                assert path.isfile(path.join(xpath, "CONTCAR"))
                assert path.isfile(path.join(xpath, ".matdb.module"))
                target = path.join(xpath, name)
                symlink(target, path.join(folder, dft))

@pytest.fixture()
def Pd(tmpdir):
    from matdb.utility import relpath
    from matdb.database import Controller
    from os import mkdir, symlink, remove, path

    target = relpath("./tests/Pd/matdb")
    dbdir = str(tmpdir.join("pd_db"))
    mkdir(dbdir)
    
    #We need to copy the POSCAR over from the testing directory to the temporary
    #one.
    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")
    mkdir(path.join(dbdir,"seed"))
    copy(POSCAR, path.join(dbdir,"seed","Pd"))
    if path.isfile("matdb.yml"):
        remove("matdb.yml")
    symlink("{}.yml".format(target),"matdb.yml")
    
    result = Controller("matdb", dbdir)
    remove("matdb.yml")
    result = Controller(target, dbdir)
    return result

@pytest.fixture()
def dynPd(Pd):
    Pd.setup()
    
    #First, we need to copy the FORCE_SETS and total_dos.dat files so that we
    #don't have to recompile those (they are tested elsewhere).
    from matdb.utility import symlink    
    troot = path.join(reporoot, "tests", "data", "Pd", "hessian")
    files = ["FORCE_SETS", "total_dos.dat", "mesh.yaml"]
    for seq in Pd.find("hessian/phonon/Pd/dim-*"):
        for filename in files:
            target = path.join(seq.root, "phonopy", filename)
            source = path.join(troot, "{0}__{1}".format(filename, seq.parent.name))
            symlink(target, source)

    return Pd

# # @pytest.mark.skip()
# def test_Pd_phonplot(dynPd, tmpdir):
#     """Tests the plotting of phonon bands for supercell convergence test in Pd.
#     """
#     from matdb.plotting.comparative import band_plot
#     dbs = dynPd.find("Pd.phonon-*.hessian")
#     target = str(tmpdir.join("Pd.phonon-convergence.pdf"))
#     args = {
#         "dim": 3,
#         "save": target
#     }
#     band_plot(dbs, **args)
#     assert path.isfile(target)
    
def test_Pd_setup(Pd):
    """Makes sure the initial folders were setup according to the spec.
    """
    Pd.setup()
    modelroot = path.join(Pd.root, "Manual","phonon","Pd")
    assert Pd["Manual/phonon/Pd/"].root == modelroot
    
    #The matdb.yml file specifies the following database:
    dbs = ["Manual/phonon/Pd/"]
    #Each one should have a folder for: ["hessian", "modulations"]
    #On the first go, the modulations folder will be empty because the DFT
    #calculations haven't been performed yet. However, hessian should have DFT
    #folders ready to go.
    folders = {
        "__files__": ["compute.pkl","jobfile.sh"],
        "S1.1": {
            "__files__": ["INCAR", "POSCAR", "POTCAR", "PRECALC", "KPOINTS"]
        }
    }

    from matdb.utility import compare_tree
    for db in dbs:
        dbfolder = path.join(Pd.root, db)
        compare_tree(dbfolder, folders)
    
def test_steps(Pd):
    """Tests compilation of all steps in the database.
    """
    assert Pd.steps() == ['manual/phonon']
    Pd.setup()
    steps = sorted(['manual/phonon/Pd/dim-2.00', 'manual/phonon/Pd/dim-4.00',
                    'manual/phonon/Pd/dim-16.00', 'manual/phonon/Pd/dim-27.00',
                    'manual/phonon/Pd/dim-32.00','manual/phonon/Pd/dim-8.00'])
    assert Pd.steps() == steps
    
    seqs = sorted(['Pd/dim-2.00', 'Pd/dim-16.00', 'Pd/dim-32.00',
                   'Pd/dim-4.00', 'Pd/dim-27.00','Pd/dim-8.00'])
    assert Pd.sequences() == seqs

def test_find(Pd):
    """Tests the find function with pattern matching.
    """
    Pd.setup()
    steps = Pd.find("maual/phonon/Pd/dim-*")
    model = ['manual', 'manual', 'manual','manual','manual','manual']
    assert model == [s.parent.name for s in steps]
    model = [path.join(Pd.root,'Manual/phonon/Pd/dim-2.00'),
             path.join(Pd.root,'Manual/phonon/Pd/dim-4.00'),
             path.join(Pd.root,'Manual/phonon/Pd/dim-8.00'),
             path.join(Pd.root,'Manual/phonon/Pd/dim-16.00'),
             path.join(Pd.root,'Manual/phonon/Pd/dim-27.00'),
             path.join(Pd.root,'Manual/phonon/Pd/dim-32.00')]
    assert sorted(model) == sorted([s.root for s in steps])

    steps = Pd.find("manual/phonon/Pd")
    model = ['manual']
    assert model == [s.parent.name for s in steps]
    model = [path.join(Pd.root,'Manual/phonon/Pd')]
    assert model == [s.root for s in steps]

    # test uuid finding.
    assert all([Pd.find(s.uuid)==s for s in steps])
    
# @pytest.mark.skip()    
def test_Pd_hessian(Pd):
    """Tests the `niterations` functionality and some of the standard
    methods of the class on simple Pd.

    """
    Pd.setup()
    
    #Test the status, we should have some folder ready to execute.
    Pd.status()
    Pd.status(busy=True)

    #Test execution command; this uses stubs for all of the commands that would
    #be executed.
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"})
    folder = path.join(reporoot, "tests", "data", "Pd", "recover")
    _mimic_vasp(folder, Pd.root)

    #Now that we have vasp files, we can queue recovery. The recover
    #`vasprun.xml` files that we linked to are not complete for all the
    #structures (on purpose).
    Pd.recover()
    recoveries = ["hessian/phonon/Pd/dim-16.00",
                  "hessian/phonon/Pd/dim-32.00",
                  "hessian/phonon/Pd/dim-27.00"]
    okay = ["hessian/phonon/Pd/dim-2.00",
            "hessian/phonon/Pd/dim-4.00"]
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
    recoveries = ["hessian/phonon/Pd/dim-32.00"]
    okay = ["hessian/phonon/Pd/dim-2.00",
            "hessian/phonon/Pd/dim-16.00",
            "hessian/phonon/Pd/dim-4.00",
            "hessian/phonon/Pd/dim-27.00"]
    for rkey in recoveries:
        assert path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
    for rkey in okay:
        assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))

    #Now copy the final files, which are all good.
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"}, recovery=True)
    folder = path.join(reporoot, "tests", "data", "Pd", "complete")
    _mimic_vasp(folder, Pd.root)

    Pd.recover()
    okay = ["hessian/phonon/Pd/dim-2.00",
            "hessian/phonon/Pd/dim-16.00",
            "hessian/phonon/Pd/dim-4.00",
            "hessian/phonon/Pd/dim-27.00",
            "hessian/phonon/Pd/dim-32.00"]
    for rkey in okay:
        assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
    
    #We are finally ready to cleanup the phonon database.
    Pd.cleanup()
    
    for rkey in okay:
        assert path.isfile(path.join(Pd[rkey].root, "phonopy", "FORCE_SETS"))
        assert path.isfile(path.join(Pd[rkey].root, "phonopy", "total_dos.dat"))
