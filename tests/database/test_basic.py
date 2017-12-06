"""Tests the basic database class structure from which all others
inherit.
"""
import pytest
from matdb.utility import relpath
from matdb.database import basic

@pytest.fixture(scope="module", autouse=True)
def Pd_f(request):
    """Returns a :class:`matdb.database.basic.Database` using Pd as a
    seed configuration, the linked to database tests the failurs or False
    staments for the basic Database.
    """
    from matdb.database.basic import Database

    POSCAR = relpath("./tests/Pd/POSCAR")
    from quippy.atoms import Atoms
    datoms = Atoms(POSCAR,format="POSCAR")
    incar = {"prec": 'a',"encut": 400, "isym": 0, "lwave": False, "lreal": 'auto',
             "ediff": '1e-5',"ismear": 1,"sigma": 0.1}
    kpoints = {"mindistance": 20}
    execution = {"template": 'run_array_ml.sh',"time": 4,"ntasks": 1,"nodes": 1,
                 "mem_per_cpu": 4,"job_name": 'Pd DB',"partition": 'physics',
                 "array_limit": 50, "modules_load": ['mkl/11.2.0'],"exec_path": 'vasp'}
    parent = None
    root = relpath("./tests/data/Pd/basic_fail")

    Pd_db = Database(datoms,incar,kpoints,execution,root,parent,nconfigs=3,prefix="N")
    Pd_db = Database(datoms,incar,kpoints,execution,root,parent,nconfigs=3)

    return Pd_db

@pytest.fixture(scope="module", autouse=True)
def Pd_p(request):
    """Returns a :class:`matdb.database.basic.Database` using Pd as a
    seed configuration.
    """
    from matdb.database.basic import Database

    POSCAR = relpath("./tests/Pd/POSCAR")
    from quippy.atoms import Atoms
    datoms = Atoms(POSCAR,format="POSCAR")
    incar = {"prec": 'a',"encut": 400, "isym": 0, "lwave": False, "lreal": 'auto',
             "ediff": '1e-5',"ismear": 1,"sigma": 0.1}
    kpoints = {"mindistance": 20}
    execution = {"template": 'run_array_ml.sh',"time": 4,"ntasks": 1,"nodes": 1,
                 "mem_per_cpu": 4,"job_name": 'Pd DB',"partition": 'physics',
                 "array_limit": 50, "modules_load": ['mkl/11.2.0'],"exec_path": 'vasp'}
    parent = None
    root = relpath("./tests/data/Pd/basic_pass")

    Pd_db = Database(datoms,incar,kpoints,execution,root,parent,nconfigs=3)

    return Pd_db

def test_can_execute():
    """Tests if a folder is ready to be executed.
    """

    assert not basic.can_execute("temp")

    PATH = relpath("./tests/data/Pd/basic_fail/S.1")
    assert basic.can_execute(PATH)

    PATH = relpath("./tests/data/Pd/basic_fail/S.2")
    assert not basic.can_execute(PATH)
    
def test_can_cleanup(Pd_f):
    """Tests if a folder is ready to be executed.
    """

    assert not basic.can_cleanup("temp")

    PATH = relpath("./tests/data/Pd/basic_fail/S.1")
    assert not basic.can_cleanup(PATH)

    PATH = relpath("./tests/data/Pd/basic_fail/S.2")
    assert basic.can_cleanup(PATH)

    PATH = relpath("./tests/data/Pd/basic_fail/S.3")
    assert basic.can_cleanup(PATH)
    
    PATH = relpath("./tests/data/Pd/basic_fail/S.4")
    assert not basic.can_cleanup(PATH)

    assert not Pd_f.cleanup()
    
def test_execution_fails(Pd_f):
    """ Tests of the basic Database execution related subroutines that should 
    return False or don't require a completed database.
    """
    from os import remove, path
    
    PATH = relpath("./tests/data/Pd/basic_fail/")
    jobfile = relpath("./tests/data/Pd/basic_fail/jobfile.sh")
    if path.isfile(jobfile):
        remove(jobfile)
    
    assert Pd_f.is_executing()        
    assert not Pd_f.execute()

    with open(jobfile,"w+") as f:
        f.write('/n')

    assert not Pd_f.execute()
    remove(jobfile)
        
def test_execution(Pd_p):
    """ Tests of the basic Database execution related subroutines.
    """
    from os import remove
    
    outfile = relpath("./tests/data/Pd/basic_pass/S.1/OUTCAR")
    with open(outfile,"w+") as f:
        f.write('/n')

    assert not Pd_p.execute()        
    remove(outfile)

    assert not Pd_p.execute()    
    assert Pd_p.execute(dryrun=True)
