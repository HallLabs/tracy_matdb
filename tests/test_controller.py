"""Tests the controller and database collection objects methods
directly. These tests rely on the `./tests/Pd` directory, which has the model
outputs that temporary directories will be compared to.
"""
import numpy as np
from os import path
import pytest
from matdb.utility import reporoot, relpath
from matdb.atoms import AtomsList
import six

def _mimic_vasp(folder, xroot, prefix="W.1"):
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
                xpath = path.join(xroot, path.join(*config.split("_")), prefix)
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
def Pd_copy(tmpdir):
    from matdb.utility import relpath
    from matdb.database import Controller
    from os import path, remove, mkdir

    target = relpath("./tests/Pd/matdb_copy")
    dbdir = str(tmpdir.join("pd_db_copy"))

    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")

    if path.isfile("matdb_copy.yml"):
        remove("matdb_copy.yml")
        
    result = Controller(target, dbdir)

    mkdir(path.join(dbdir,"seed"))
    copy(POSCAR, path.join(dbdir, "seed", "Pd"))
    result = Controller(target, dbdir)
    return result

@pytest.fixture()
def Pd_split(tmpdir):
    from matdb.utility import relpath
    from matdb.database import Controller
    from os import mkdir, path

    target = relpath("./tests/Pd/matdb_split")
    dbdir = str(tmpdir.join("pd_db_splits"))
    mkdir(dbdir)

    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")
    mkdir(path.join(dbdir,"seed"))
    copy(POSCAR, path.join(dbdir,"seed","Pd-1"))
    copy(POSCAR, path.join(dbdir,"seed","Pd-2"))
    copy(POSCAR, path.join(dbdir,"seed","Pd-3"))
    copy(POSCAR, path.join(dbdir,"seed","Pd-4"))
    copy(POSCAR, path.join(dbdir,"seed","Pd-5"))

    result = Controller(target, dbdir)
    return result
    
    
@pytest.fixture()
def Pd_2(tmpdir):
    from matdb.utility import relpath
    from matdb.database import Controller
    from os import mkdir, symlink, remove, path

    target = relpath("./tests/Pd/matdb")
    dbdir = str(tmpdir.join("pd_2_db"))
    if not path.isdir(dbdir):
        mkdir(dbdir)

    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")
    if not path.isdir(path.join(dbdir,"seed")):
        mkdir(path.join(dbdir,"seed"))
    copy(POSCAR, path.join(dbdir,"seed","Pd"))
    if path.isfile("matdb.yml"):
        remove("matdb.yml")
    symlink("{}.yml".format(target),"matdb.yml")

    result = Controller("matdb",dbdir)
    remove("matdb.yml")
    result = Controller(target, dbdir)
    return result

@pytest.fixture()
def AgPd(tmpdir):
    from matdb.utility import relpath
    from matdb.database import Controller
    from os import mkdir, symlink, remove, path

    target = relpath("./tests/AgPd/matdb_manual")
    dbdir = str(tmpdir.join("agpd_db"))
    mkdir(dbdir)

    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")
    mkdir(path.join(dbdir,"seed"))
    copy(POSCAR,
         path.join(dbdir,"seed","Pd"))
    if path.isfile("matdb.yml"):
        remove("matdb.yml")
    symlink("{}.yml".format(target),"matdb.yml")

    result = Controller("matdb",dbdir)
    remove("matdb.yml")
    result = Controller(target,dbdir)
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

#@pytest.mark.skip()
def test_Pd_setup(Pd, Pd_copy):
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

    #Now we will test some of the border cases of the database __init__ method
    Pd_copy.setup()

    db = "Manual/phonon/Pd"
    dbfolder = path.join(Pd_copy.root, db)
    compare_tree(dbfolder, folders)

#@pytest.mark.skip()
def test_steps(Pd):
    """Tests compilation of all steps in the database.
    """
    assert Pd.steps() == ['manual/phonon']
    Pd.setup()
    steps = sorted(['manual/phonon/Pd'])
    assert Pd.steps() == steps
    
    seqs = sorted(['Pd'])
    assert Pd.sequences() == seqs
#@pytest.mark.skip()
def test_find(Pd):
    """Tests the find function and the __getitem__ method with pattern matching.
    """
    Pd.setup()
    steps = Pd.find("manual/phonon")
    model = ['phonon']
    assert model == [s.parent.name for s in steps]
    model = [path.join(Pd.root,'Manual/phonon')]
    assert sorted(model) == sorted([s.root for s in steps])
   
    steps = Pd.find("*/phonon")
    model = ['phonon']
    assert model == [s.parent.name for s in steps]
    model = [path.join(Pd.root,'Manual/phonon')]
    assert model == [s.root for s in steps]

    steps = Pd.find("manual/phonon/Pd")
    model = ['manual']
    assert model == [s.parent.name for s in steps]
    model = [path.join(Pd.root,'Manual/phonon/Pd')]
    assert model == [s.root for s in steps]

    steps = Pd.find("phonon")
    model = ['phonon']
    assert model == [s.name for s in steps]
    model = [Pd.root]
    assert model == [s.root for s in steps]

    steps = Pd.find('*')
    model = ['phonon']
    assert model == [s.name for s in steps]
    model = [Pd.root]
    assert model == [s.root for s in steps]

    steps = Pd.find("manual/phonon/Pd/S1.1")
    model = [('phonon','manual')]
    assert model == [(s.parent.parent.name,s.name) for s in steps]
    model = [path.join(Pd.root,'Manual/phonon/Pd')]
    assert model == [s.root for s in steps]
                             
    # test uuid finding.
    assert all([Pd.find(s.uuid)==s for s in steps])

    # test the __getitem__ method
    model = 'phonon'
    modelroot = path.join(Pd.root,'Manual/phonon')
    group = Pd["manual/phonon"]
    assert group.parent.name == model
    assert group.root == modelroot

    model = 'manual'
    modelroot = path.join(Pd.root,'Manual/phonon/Pd')
    group = Pd["manual/phonon/Pd"]
    assert group.parent.name == model
    assert group.root == modelroot

    group = Pd["manual/phonon/Pd/S1.1"]
    assert group.parent.name == model
    assert group.root == modelroot

    group = Pd["enumeration/phonon"]
    assert group == None

#@pytest.mark.skip()
def test_execute(Pd, capsys):
    """Tests the execute and extract methods 
    """
    from os import path
    from matdb.utility import relpath, chdir

    Pd.status()
    output = capsys.readouterr()
    status = "ready to execute 0/0; finished executing 0/0;"
    assert status in output.out

    # test to ensure execute prints error message if ran before setup
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID":"1"})
    output = capsys.readouterr()
    status = "Group phonon.manual is not ready to execute yet, or is already executing. Done."
    assert status in output.out

    # Run setup and status to make sure the staus output is correct
    Pd.setup()
    Pd.status()
    output = capsys.readouterr()
    status = "ready to execute 1/1; finished executing 0/1;"
    assert status in output.out

    # Execute the jobfile and test an incomplete OUTCAR to check the status
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID":"1"})
    folder = path.join(reporoot, "tests", "data", "Pd", "manual_recover")
    _mimic_vasp(folder, Pd.root,"S1.1")

    Pd.status(True)
    busy_status = "Pd./Manual/phonon/Pd/S1.1"
    output = capsys.readouterr()
    assert busy_status in output.out

    #now use a complete OUTCAR and test the status again
    folder = path.join(reporoot, "tests", "data", "Pd", "manual")
    _mimic_vasp(folder,Pd.root,"S1.1")
    
    Pd.status()
    output = capsys.readouterr()
    status = "ready to execute 1/1; finished executing 1/1;"
    assert status in output.out

    # Run exctract and test to see if the status is correct
    Pd.extract()
    Pd.status()
    output = capsys.readouterr()
    status = "ready to execute 1/1; finished executing 1/1;"
    assert status in output.out

    # Run extract again to make sure the atoms.h5 files are no rewritten
    Pd.extract()

#@pytest.mark.skip()
def test_recovery(Pd):
    """Tests the rerun on unfinshed jobs
    """
    from os import path
    from matdb.utility import symlink, chdir
    from glob import glob

    Pd.setup()
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID":"1"})

    files = ["vasprun.xml", "OUTCAR"]
    folder = path.join(reporoot, "tests", "data", "Pd", "manual_recover")
    with chdir(folder):
        for vaspfile in files:
            pattern = vaspfile + "__*"
            for dft in glob(pattern):
                name, config = dft.split("__")
                xpath = path.join(Pd.root, path.join(*config.split("_")), "S1.1")
                target = path.join(xpath, name)
                symlink(target, path.join(folder, dft))

    Pd.extract()
    
    Pd.recover(True)
    assert path.isfile(path.join(Pd.root,"Manual","phonon","Pd","recovery.sh"))
    assert path.isfile(path.join(Pd.root,"Manual","phonon","Pd","failures"))

    folder = path.join(reporoot, "tests", "data", "Pd", "manual")
    _mimic_vasp(folder,Pd.root,"S1.1")
    Pd.recover(True)
    assert not path.isfile(path.join(Pd.root,"Manual","phonon","Pd","recovery.sh"))
    assert not path.isfile(path.join(Pd.root,"Manual","phonon","Pd","failures"))

#@pytest.mark.skip()
def test_hash(Pd,Pd_2):
    """Tests the hash_dbs and verify_hash methods
    """
    from os import path, chdir
    Pd.setup()
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID":"1"})
    folder = path.join(reporoot, "tests", "data", "Pd", "manual")
    _mimic_vasp(folder,Pd.root,"S1.1")
    db_hash = Pd.hash_dbs()
    assert Pd.verify_hash(db_hash)
    
    # test to make sure a change in the databas produces a different hash
    Pd_2.setup()
    Pd_2.execute(env_vars={"SLURM_ARRAY_TASK_ID":"1"})
    _mimic_vasp(folder,Pd_2.root,"S1.1")
    assert Pd_2.hash_dbs != db_hash

#@pytest.mark.skip()
def test_finalize(Pd):
    """ Test the finalize function in the controller module
    """
    from os import path
    Pd.setup()
    Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID":"1"})
    folder = path.join(reporoot,"tests","data","Pd","manual")
    _mimic_vasp(folder,Pd.root,"S1.1")

    Pd.extract()
    Pd.split()
    Pd.finalize()

    from matdb.io import load_dict_from_h5
    from matdb import __version__
    import h5py

    str_ver = []
    for item in __version__:
        str_ver.append(str(item))
    mdb_ver = ".".join(str_ver)
    target = path.join(Pd.root, "final_{}.h5".format(mdb_ver))
    with h5py.File(target, "r") as hf:
        loaded_final = load_dict_from_h5(hf)
    assert path.isfile(target)

def test_split(Pd_split):
    """ Test the split function in the controller object
    """
    from os import path
    Pd_split.setup()
    Pd_split.execute(env_vars={"SLURM_ARRAY_TASK_ID":"1"})
    folder = path.join(reporoot,"tests","data","Pd","manual_split")
    _mimic_vasp(folder,Pd_split.root,"S1.1")

    Pd_split.extract()
    Pd_split.split()

    for dbname, db in Pd_split.collections.items():
        for s, p in db.splits.items():
            tfile = path.join(db.train_file(s).format(s))
            hfile = path.join(db.holdout_file(s).format(s))
            sfile = path.join(db.super_file(s).format(s))

            tal = AtomsList(tfile)
            hal = AtomsList(hfile)
            sal = AtomsList(sfile)

            assert len(tal) == int(np.ceil(5*p))
            assert len(hal) == int(np.ceil((5-len(tal))*p))
            assert len(sal) == 5-len(tal)-len(hal)
    
@pytest.mark.skip()    
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
