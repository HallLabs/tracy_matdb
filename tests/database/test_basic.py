"""Tests the basic database class structure from which all others
inherit.
"""
import pytest
from os import path, mkdir
from matdb.utility import relpath
from matdb.database import basic
import quippy
import numpy as np
import ase
from matdb.atoms import Atoms

def test_json(tmpdir):
    """Tests whether an atoms object with calculated parameters can be saved to
    JSON and then restored.
    """
    from matdb.database.basic import atoms_to_json, atoms_from_json
    from matdb.calculators import Quip
    target = str(tmpdir.join("database_tojson"))
    if not path.isdir(target):
        mkdir(target)

    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])#quippy.diamond(5.43,14)
    potSW = Quip(atSi, target, ["IP SW"])
    atSi.set_calculator(potSW)
    # potSW.calc(atSi, energy=True, force=True, virial=True)
    atSi.properties["rand"] = np.random.randint(0, 100, 8)
    atSi.properties["rand"] = np.random.randint(0, 100, 8)
    atSi.add_param("hessian_1",np.random.random())
    atSi.add_param("energy",-34.67998604386673)
    atSi.add_param("virial",[[0.05319905934507262, 0.0, 0.0], [0.0, 0.05319905934507262, 0.0], [0.0, 0.0, 0.05319905934507262]])
    atSi.add_property("force",[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    
    
    atoms_to_json(atSi, target)
    atR = atoms_from_json(target)

    assert atR.energy == atSi.energy
    assert isinstance(atR, Atoms)
    assert np.allclose(atR.force, atSi.force)
    assert np.allclose(atR.virial, atSi.virial)
    assert atR.hessian_1 == atSi.hessian_1
    assert np.allclose(atR.properties["rand"], atSi.properties["rand"])
    assert np.allclose(atR.positions, atSi.positions)

# @pytest.fixture(scope="module", autouse=True)
# def Pd_f(request):

#     """Returns a :class:`matdb.database.basic.Database` using Pd as a
#     seed configuration, the linked to database tests the failurs or False
#     staments for the basic Database.
#     """
#     from matdb.database.basic import Group

#     POSCAR = relpath("./tests/Pd/POSCAR")
#     from quippy.atoms import Atoms
#     datoms = Atoms(POSCAR,format="POSCAR")
#     incar = {"prec": 'a',"encut": 400, "isym": 0, "lwave": False, "lreal": 'auto',
#              "ediff": '1e-5',"ismear": 1,"sigma": 0.1}
#     kpoints = {"mindistance": 20}
#     execution = {"template": 'run_array_ml.sh',"time": 4,"ntasks": 1,"nodes": 1,
#                  "mem_per_cpu": 4,"job_name": 'Pd DB',"partition": 'physics',
#                  "array_limit": 50, "modules_load": ['mkl/11.2.0'],"exec_path": 'vasp'}
#     parent = None
#     root = relpath("./tests/data/Pd/basic_fail")

#     calcargs = incar.copy()
#     calcargs["name"] = "Vasp"
#     calcargs["kpoints"] = kpoints    
#     Pd_db = Group(root=root, parent=parent, execution=execution,
#                   nconfigs=3, prefix="N", calculator=calcargs)
#     return Pd_db

# @pytest.fixture(scope="module", autouse=True)
# def Pd_p(request):
#     """Returns a :class:`matdb.database.basic.Group` using Pd as a
#     seed configuration.
#     """
#     from matdb.database.basic import Group
#     from quippy.atoms import Atoms
#     incar = {"prec": 'a',"encut": 400, "isym": 0, "lwave": False, "lreal": 'auto',
#              "ediff": '1e-5',"ismear": 1,"sigma": 0.1}
#     kpoints = {"mindistance": 20}
#     execution = {"template": 'run_array_ml.sh',"time": 4,"ntasks": 1,"nodes": 1,
#                  "mem_per_cpu": 4,"job_name": 'Pd DB',"partition": 'physics',
#                  "array_limit": 50, "modules_load": ['mkl/11.2.0'],"exec_path": 'vasp'}
#     parent = None
#     root = relpath("./tests/data/Pd/basic_pass")

#     calculator=incar.copy()
#     calculator["name"] = "Vasp"
#     calculator["kpoints"] = kpoints
#     Pd_db = Group(root=root, parent=parent, seed="Pd", execution=execution,
#                   nconfigs=3, calculator=calculator)
#     return Pd_db

# def test_can_execute():
#     """Tests if a folder is ready to be executed.
#     """

#     assert not basic.can_execute("temp")

#     PATH = relpath("./tests/data/Pd/basic_fail/S.1")
#     assert basic.can_execute(PATH)

#     PATH = relpath("./tests/data/Pd/basic_fail/S.2")
#     assert not basic.can_execute(PATH)
    
# def test_can_cleanup(Pd_f):
#     """Tests if a folder is ready to be executed.
#     """

#     assert not basic.can_cleanup("temp")

#     PATH = relpath("./tests/data/Pd/basic_fail/S.1")
#     assert not basic.can_cleanup(PATH)

#     PATH = relpath("./tests/data/Pd/basic_fail/S.2")
#     assert basic.can_cleanup(PATH)

#     PATH = relpath("./tests/data/Pd/basic_fail/S.3")
#     assert basic.can_cleanup(PATH)
    
#     PATH = relpath("./tests/data/Pd/basic_fail/S.4")
#     assert not basic.can_cleanup(PATH)

#     assert not Pd_f.cleanup()
    
# def test_execution_fails(Pd_f):
#     """ Tests of the basic Group execution related subroutines that should 
#     return False or don't require a completed database.
#     """
#     from os import remove, path
    
#     PATH = relpath("./tests/data/Pd/basic_fail/")
#     jobfile = relpath("./tests/data/Pd/basic_fail/jobfile.sh")
#     if path.isfile(jobfile):
#         remove(jobfile)
    
#     assert Pd_f.is_executing()        
#     assert not Pd_f.execute()

#     with open(jobfile,"w+") as f:
#         f.write('/n')

#     assert not Pd_f.execute()
#     remove(jobfile)
        
# def test_execution(Pd_p):
#     """ Tests of the basic Group execution related subroutines.
#     """
#     from os import remove
    
#     outfile = relpath("./tests/data/Pd/basic_pass/S.1/OUTCAR")
#     with open(outfile,"w+") as f:
#         f.write('/n')

#     assert not Pd_p.execute()        
#     remove(outfile)

#     assert not Pd_p.execute()    
#     assert Pd_p.execute(dryrun=True)

# def test_jobfile(Pd_p):
#     """Tests the jobfile creation corner cases.
#     """

#     Pd_p.jobfile()

# # Needs parent to perform
# #   Pd_p.jobfile(recovery=True)
# #
# # def test_recover(Pd_p):
# #     """Tests the recovery function.
# #     """
# #     Pd_p.recover()
# #
# # def test_create(Pd_p):
# #     """Tests the create function's corner cases.
# #     """
# #
# #     POSCAR = relpath("./tests/Pd/POSCAR")
# #     from quippy.atoms import Atoms
# #     datoms = Atoms(POSCAR,format="POSCAR")
# #     Pd_p.create(datoms)

# def test_setup(Pd_f):
#     """Tests the setup funciton's corner cases.
#     """
#     from os import remove
#     outfile = relpath("./tests/data/Pd/basic_pass/S.1/OUTCAR")
#     with open(outfile,"w+") as f:
#         f.write('/n')

#     assert Pd_f.setup()        
#     remove(outfile)

# def test_status(Pd_p):
#     """Tests the status of the database.
#     """

#     assert "ready to execute 3/3; finished executing 0/3;" == Pd_p.status()

#     status = {'ready': {'/root/codes/matdb/tests/data/Pd/basic_pass/S.1': True, '/root/codes/matdb/tests/data/Pd/basic_pass/S.2': True, '/root/codes/matdb/tests/data/Pd/basic_pass/S.3': True}, 'stats': {'ready': 3, 'done': 0, 'N': 3}, 'busy': False, 'done': {'/root/codes/matdb/tests/data/Pd/basic_pass/S.1': False, '/root/codes/matdb/tests/data/Pd/basic_pass/S.2': False, '/root/codes/matdb/tests/data/Pd/basic_pass/S.3': False}}
#     assert status == Pd_p.status(print_msg=False)

# def test_xyz(Pd_p):
#     """Tests the coversion to xyz files.
#     """
#     from os import symlink, remove, path
#     src = relpath("./tests/data/Pd/complete/vasprun.xml__Pd.phonon-16")
#     for i in range(1,4):
#         dest = relpath("./tests/data/Pd/basic_pass/S.{}/vasprun.xml".format(i))
#         symlink(src,dest)
    
#     assert Pd_p.xyz(combine=True)
#     xyzf = relpath("./tests/data/Pd/basic_pass/output.xyz")
#     assert path.isfile(xyzf)
#     remove(xyzf)
#     for i in range(1,4):
#         dest = relpath("./tests/data/Pd/basic_pass/S.{}/vasprun.xml".format(i))
#         remove(dest)
#         dest = relpath("./tests/data/Pd/basic_pass/S.{}/output.xyz".format(i))
#         remove(dest)
#         dest = relpath("./tests/data/Pd/basic_pass/S.{}/output.xyz.idx".format(i))
#         remove(dest)
    
# def test_tarball(Pd_p):
#     """Tests the tar ball creation.
#     """

#     from os import symlink, remove, path
#     src = relpath("./tests/data/Pd/complete/OUTCAR__Pd.phonon-16")
#     for i in range(1,4):
#         dest = relpath("./tests/data/Pd/basic_pass/S.{}/OUTCAR".format(i))
#         symlink(src,dest)

#     Pd_p.tarball()
#     target = relpath("./tests/data/Pd/basic_pass/output.tar.gz")
#     assert path.isfile(target)

#     remove(target)
#     for i in range(1,4):
#         dest = relpath("./tests/data/Pd/basic_pass/S.{}/OUTCAR".format(i))
#         remove(dest)
    
