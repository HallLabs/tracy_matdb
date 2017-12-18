"""Tests the controller and database collection objects methods
directly. These tests rely on the `./tests/Pd` directory, which has the model
outputs that temporary directories will be compared to.
"""
import numpy as np
from os import path
import pytest
from matdb.utility import reporoot, relpath
    
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

# @pytest.fixture()
# def Pd(tmpdir):
#     from matdb.utility import relpath
#     from matdb.database.controller import Controller
#     from os import mkdir, symlink, remove

#     target = relpath("./tests/Pd/matdb")
#     dbdir = str(tmpdir.join("pd_db"))
#     mkdir(dbdir)
    
#     #We need to copy the POSCAR over from the testing directory to the temporary
#     #one.
#     from shutil import copy
#     POSCAR = relpath("./tests/Pd/POSCAR")
#     copy(POSCAR, dbdir)
#     symlink("{}.yml".format(target),"matdb.yml")
    
#     result = Controller("matdb", dbdir)
#     remove("matdb.yml")
#     result = Controller(target, dbdir)
#     return result

# @pytest.fixture()
# def dynPd(Pd):
#     Pd.setup()
    
#     #First, we need to copy the FORCE_SETS and total_dos.dat files so that we
#     #don't have to recompile those (they are tested elsewhere).
#     from matdb.utility import symlink    
#     troot = path.join(reporoot, "tests", "data", "Pd", "dynmatrix")
#     files = ["FORCE_SETS", "total_dos.dat", "mesh.yaml"]
#     for seq in Pd.find("Pd.phonon-*.dynmatrix"):
#         for filename in files:
#             target = path.join(seq.root, "phonopy", filename)
#             source = path.join(troot, "{0}__{1}".format(filename, seq.parent.name))
#             symlink(target, source)

#     return Pd

# # @pytest.mark.skip()
# def test_Pd_phonplot(dynPd, tmpdir):
#     """Tests the plotting of phonon bands for supercell convergence test in Pd.
#     """
#     from matdb.plotting.comparative import band_plot
#     dbs = dynPd.find("Pd.phonon-*.dynmatrix")
#     target = str(tmpdir.join("Pd.phonon-convergence.pdf"))
#     args = {
#         "dim": 3,
#         "save": target
#     }
#     band_plot(dbs, **args)
#     assert path.isfile(target)
    
# def test_Pd_setup(Pd):
#     """Makes sure the initial folders were setup according to the spec.
#     """
#     Pd.setup()
#     modelroot = path.join(Pd.root, "Pd.phonon-2", "dynmatrix")
#     assert Pd["Pd.phonon-2.dynmatrix"].root == modelroot
    
#     #The matdb.yml file specifies the following databases:
#     dbs = ["Pd.phonon-{}".format(i) for i in (2, 4, 16, 32, 54)]
#     #Each one should have a folder for: ["dynmatrix", "modulations"]
#     #On the first go, the modulations folder will be empty because the DFT
#     #calculations haven't been performed yet. However, dynmatrix should have DFT
#     #folders ready to go.
#     folders = {
#         "dynmatrix": {
#             "__files__": ["INCAR", "PRECALC"],
#             "phonopy": {
#                 "__files__": ["POSCAR", "POSCAR-001", "disp.yaml", "phonopy_disp.yaml"]
#             },
#             "phoncache": {},
#             "W.1": {
#                 "__files__": ["INCAR", "POSCAR", "POTCAR", "PRECALC", "KPOINTS"]
#             }
#         }
#     }

#     from matdb.utility import compare_tree
#     for db in dbs:
#         dbfolder = path.join(Pd.root, db)
#         compare_tree(dbfolder, folders)
    
# def test_steps(Pd):
#     """Tests compilation of all steps in the database.
#     """
#     steps = ['Pd.modulate.dynmatrix', 'Pd.modulate.modulations',                          
#              'Pd.phonon-16.dynmatrix', 'Pd.phonon-2.dynmatrix',
#              'Pd.phonon-32.dynmatrix', 'Pd.phonon-4.dynmatrix',
#              'Pd.phonon-54.dynmatrix'] 
    
#     seqs = sorted(['Pd.phonon-2', 'Pd.phonon-16', 'Pd.phonon-32',
#                    'Pd.phonon-4', 'Pd.phonon-54', 'Pd.modulate'])
#     assert Pd.steps() == steps
#     assert Pd.sequences() == seqs

# def test_find(Pd):
#     """Tests the find function with pattern matching.
#     """
#     assert Pd['Pd.modulate.dynmatrix'] is Pd['Pd.phonon-16.dynmatrix']
    
#     steps = Pd.find("Pd.phonon*.dynmatrix")
#     model = ['Pd.phonon-2', 'Pd.phonon-16', 'Pd.phonon-32', 'Pd.phonon-4', 'Pd.phonon-54']
#     assert model == [s.parent.name for s in steps]

#     steps = Pd.find("Pd.phonon-*2.dynmatrix")
#     model = ['Pd.phonon-2', 'Pd.phonon-32']
#     assert model == [s.parent.name for s in steps]
    
# def test_Pd_modulation(dynPd):
#     """Tests generation of modulated configurations along phonon modes.
#     """
#     from matdb.utility import symlink
#     #Call setup again since we have force_sets and DOS available.
#     dynPd.cleanup()
#     dynPd.setup()

#     folders = {
#         "modulations": {
#             "__files__": ["INCAR", "PRECALC", "jobfile.sh"]
#         }
#     }

#     files = ["INCAR", "POSCAR", "POTCAR", "PRECALC", "KPOINTS"]
#     for i in range(1, 6):
#         key = "M.{0:d}".format(i)
#         folders["modulations"][key] = {}
#         folders["modulations"][key]["__files__"] = files

#     from matdb.utility import compare_tree
#     dbfolder = path.join(dynPd.root, "Pd.modulate")
#     compare_tree(dbfolder, folders)    
#     folder = path.join(reporoot, "tests", "data", "Pd", "modulations")
#     dest = dynPd["Pd.modulate"].steps["modulations"].root
#     for i in range(1,6):
#         symlink(path.join(dest,"M.{}".format(i),"OUTCAR"),path.join(folder,"OUTCAR__M.{}".format(i)))
#         symlink(path.join(dest,"M.{}".format(i),"output.xyz"),path.join(folder,"output.xyz__M.{}".format(i)))
#     dynPd["Pd.modulate"].steps["modulations"].cleanup()
#     dynPd["Pd.modulate"].steps["modulations"].ready()
#     dynPd["Pd.modulate"].steps["modulations"].setup()
#     # dynPd["Pd.modulate"].steps["modulations"].setup(rerun=True)
    
# # @pytest.mark.skip()    
# def test_Pd_dynmatrix(Pd):
#     """Tests the `niterations` functionality and some of the standard
#     methods of the class on simple Pd.

#     """
#     Pd.setup()
        
#     #Test the status, we should have some folder ready to execute.
#     Pd.status()
#     Pd.status(busy=True)

#     #Test execution command; this uses stubs for all of the commands that would
#     #be executed.
#     Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"})
#     folder = path.join(reporoot, "tests", "data", "Pd", "recover")
#     _mimic_vasp(folder, Pd.root)

#     #Now that we have vasp files, we can queue recovery. The recover
#     #`vasprun.xml` files that we linked to are not complete for all the
#     #structures (on purpose).
#     Pd.recover()
#     recoveries = ["Pd.phonon-16.dynmatrix",
#                   "Pd.phonon-32.dynmatrix",
#                   "Pd.phonon-54.dynmatrix"]
#     okay = ["Pd.phonon-2.dynmatrix",
#             "Pd.phonon-4.dynmatrix"]
#     for rkey in recoveries:
#         assert path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
#     for rkey in okay:
#         assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))

#     Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"}, recovery=True)
#     folder = path.join(reporoot, "tests", "data", "Pd", "rexecute")
#     _mimic_vasp(folder, Pd.root)   

#     #Now that we have recovered vasp files, we can queue a *second*
#     #recovery. All but one of the `vasprun.xml` files that we linked to are
#     #complete for all the structures. This test that we can recover and execute
#     #multiple times in order.
#     Pd.recover()
#     recoveries = ["Pd.phonon-32.dynmatrix"]
#     okay = ["Pd.phonon-2.dynmatrix",
#             "Pd.phonon-16.dynmatrix",
#             "Pd.phonon-4.dynmatrix",
#             "Pd.phonon-54.dynmatrix"]
#     for rkey in recoveries:
#         assert path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
#     for rkey in okay:
#         assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))

#     #Now copy the final files, which are all good.
#     Pd.execute(env_vars={"SLURM_ARRAY_TASK_ID": "1"}, recovery=True)
#     folder = path.join(reporoot, "tests", "data", "Pd", "complete")
#     _mimic_vasp(folder, Pd.root)

#     Pd.recover()
#     okay = ["Pd.phonon-2.dynmatrix",
#             "Pd.phonon-16.dynmatrix",
#             "Pd.phonon-4.dynmatrix",
#             "Pd.phonon-54.dynmatrix",
#             "Pd.phonon-32.dynmatrix"]
#     for rkey in okay:
#         assert not path.isfile(path.join(Pd[rkey].root, "recovery.sh"))
    
#     #We are finally ready to cleanup the phonon database.
#     Pd.cleanup()
    
#     for rkey in okay:
#         assert path.isfile(path.join(Pd[rkey].root, "phonopy", "FORCE_SETS"))
#         assert path.isfile(path.join(Pd[rkey].root, "phonopy", "total_dos.dat"))
    
# def test_split(Pd):
#     """Tests the splitting logic and that the ids from original
#     randomization are saved correctly.
#     """
#     from os import remove, path, symlink
    
#     Pd.setup()
#     Pd["Pd.phonon-4"]._split("A")
#     remove(path.join(Pd["Pd.phonon-4"].root,"A-super.xyz"))
#     Pd["Pd.phonon-4"]._split("A")

#     POSCAR = relpath("./tests/Pd/POSCAR")
#     output = path.join(Pd["Pd.phonon-4"].root,"output.xyz")
#     symlink(POSCAR,output)
#     list(Pd["Pd.phonon-4"].isteps)[0][1].configs = {"1":Pd["Pd.phonon-4"].root,"2":''}
#     remove(path.join(Pd["Pd.phonon-4"].root,"A-super.xyz"))
#     remove(path.join(Pd["Pd.phonon-4"].root,"A-ids.pkl"))
#     Pd["Pd.phonon-4"]._split("A")
#     remove(output)

#     Pd.split()

# def test_Repeater(Pd):
#     """Tests the initialization of the squence repeater.
#     """

#     from matdb.database.controller import Repeater
#     root = Pd["Pd.phonon-4"].root
#     POSCAR = relpath("./tests/Pd/POSCAR")
#     temp = Repeater("temp",POSCAR,root,Pd,[{'phonopy': {'dim': [2, 2, 2]}, 'kpoints':
#                                              {'mindistance': 30}, 'dosmesh': [10, 10, 10],
#                                              'type': 'phonon.DynMatrix',
#                                              'bandmesh': [13, 13, 13]}],
#                      niterations=[{'phonopy.dim': [0, 1, 1, 1, 0, 1, 1, 1, 0]}])

def test_db_paramgrid():
    """Tests the construction of the parametr grid for a given database.
    """

    from matdb.utility import flatten_dict
    from matdb.database.controller import db_pgrid
    
    db_info = {'phonopy': {'dim*': [[2, 2, 2], [0, 1, 1, 1, 0, 1, 1, 1, 0],
                                    [0, 2, 2, 2, 0, 2, 2, 2, 0], [-2, 2, 2, 2, -2, 2, 2, 2, -2],
                                    [-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3, 3]],
                           'dim_suffix': 'linalg:det'},
               'parent':  'controller', 'dosmesh*': [[10, 10, 10],[12,12,12]], 'atoms': 'Pd',
               'root': '/root/codes/matdb/tests/Pd', 'bandmesh': [13, 13, 13]}

    test = db_pgrid(flatten_dict(db_info))
    actual = ([([2, 2, 2], [10, 10, 10]), ([2, 2, 2], [12, 12, 12]),
                   ([0, 1, 1, 1, 0, 1, 1, 1, 0], [10, 10, 10]),
                   ([0, 1, 1, 1, 0, 1, 1, 1, 0], [12, 12, 12]),
                   ([0, 2, 2, 2, 0, 2, 2, 2, 0], [10, 10, 10]),
                   ([0, 2, 2, 2, 0, 2, 2, 2, 0], [12, 12, 12]),
                   ([-2, 2, 2, 2, -2, 2, 2, 2, -2], [10, 10, 10]),
                   ([-2, 2, 2, 2, -2, 2, 2, 2, -2], [12, 12, 12]),
                   ([-1, 1, 1, 1, -1, 1, 1, 1, -1], [10, 10, 10]),
                   ([-1, 1, 1, 1, -1, 1, 1, 1, -1], [12, 12, 12]),
                   ([3, 3, 3], [10, 10, 10]),  ([3, 3, 3], [12, 12, 12])],
                  ['dim', 'dosmesh'],
                  [(8.0, 1.0), (8.0, 2.0), (2.0, 1.0), (2.0, 2.0), (16.0, 1.0), (16.0, 2.0),
                   (32.0, 1.0), (32.0, 2.0), (4.0, 1.0), (4.0, 2.0), (27.0, 1.0), (27.0, 2.0)])

    assert test==actual

    db_info = {'phonopy': {'dim*': [[2, 2, 2], [0, 1, 1, 1, 0, 1, 1, 1, 0],
                                    [0, 2, 2, 2, 0, 2, 2, 2, 0], [-2, 2, 2, 2, -2, 2, 2, 2, -2],
                                    [-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3, 3]],
                           'dim_suffix': 'linalg:det'},
               'parent':  'controller', 'dosmesh*': [[10, 10, 10],[12,12,12]],
               'dosmesh_suffix*':[30,36],'atoms': 'Pd',
               'root': '/root/codes/matdb/tests/Pd', 'bandmesh': [13, 13, 13]}
    
    test = db_pgrid(flatten_dict(db_info))
    actual = ([([2, 2, 2], [10, 10, 10]), ([2, 2, 2], [12, 12, 12]),
                   ([0, 1, 1, 1, 0, 1, 1, 1, 0], [10, 10, 10]),
                   ([0, 1, 1, 1, 0, 1, 1, 1, 0], [12, 12, 12]),
                   ([0, 2, 2, 2, 0, 2, 2, 2, 0], [10, 10, 10]),
                   ([0, 2, 2, 2, 0, 2, 2, 2, 0], [12, 12, 12]),
                   ([-2, 2, 2, 2, -2, 2, 2, 2, -2], [10, 10, 10]),
                   ([-2, 2, 2, 2, -2, 2, 2, 2, -2], [12, 12, 12]),
                   ([-1, 1, 1, 1, -1, 1, 1, 1, -1], [10, 10, 10]),
                   ([-1, 1, 1, 1, -1, 1, 1, 1, -1], [12, 12, 12]),
                   ([3, 3, 3], [10, 10, 10]),  ([3, 3, 3], [12, 12, 12])],
                  ['dim', 'dosmesh'],
                  [(8.0, 30.0), (8.0, 36.0), (2.0, 30.0), (2.0, 36.0), (16.0, 30.0),
                   (16.0, 36.0), (32.0, 30.0), (32.0, 36.0), (4.0, 30.0), (4.0, 36.0),
                   (27.0, 30.0), (27.0, 36.0)])

    assert test==actual
    
    db_info = {'phonopy': {'dim*': [[2, 2, 2], [0, 1, 1, 1, 0, 1, 1, 1, 0],
                                    [0, 2, 2, 2, 0, 2, 2, 2, 0], [-2, 2, 2, 2, -2, 2, 2, 2, -2],
                                    [-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3, 3]],
                           'dim_suffix': 'linalg:det'},
               'parent':  'controller', 'dosmesh*': [[10, 10, 10],[12,12,12]],
               'dosmesh_suffix':'random','atoms': 'Pd',
               'root': '/root/codes/matdb/tests/Pd', 'bandmesh': [13, 13, 13]}

    test = db_pgrid(flatten_dict(db_info),ignore_=["seed"])
    actual = ([([2, 2, 2], [10, 10, 10]), ([2, 2, 2], [12, 12, 12]),
                   ([0, 1, 1, 1, 0, 1, 1, 1, 0], [10, 10, 10]),
                   ([0, 1, 1, 1, 0, 1, 1, 1, 0], [12, 12, 12]),
                   ([0, 2, 2, 2, 0, 2, 2, 2, 0], [10, 10, 10]),
                   ([0, 2, 2, 2, 0, 2, 2, 2, 0], [12, 12, 12]),
                   ([-2, 2, 2, 2, -2, 2, 2, 2, -2], [10, 10, 10]),
                   ([-2, 2, 2, 2, -2, 2, 2, 2, -2], [12, 12, 12]),
                   ([-1, 1, 1, 1, -1, 1, 1, 1, -1], [10, 10, 10]),
                   ([-1, 1, 1, 1, -1, 1, 1, 1, -1], [12, 12, 12]),
                   ([3, 3, 3], [10, 10, 10]),  ([3, 3, 3], [12, 12, 12])],
                  ['dim', 'dosmesh'],
                  [(8.0, 1.0), (8.0, 2.0), (2.0, 1.0), (2.0, 2.0), (16.0, 1.0), (16.0, 2.0),
                   (32.0, 1.0), (32.0, 2.0), (4.0, 1.0), (4.0, 2.0), (27.0, 1.0), (27.0, 2.0)])

    assert test==actual

def test_ParameterGrid():
    """Tests the creation of a ParamaterGrid and it's functionality.
    """

    from matdb.database.controller import ParameterGrid
    db_info = {'phonopy': {'dim*': [[2, 2, 2], [0, 1, 1, 1, 0, 1, 1, 1, 0],
                                    [0, 2, 2, 2, 0, 2, 2, 2, 0], [-2, 2, 2, 2, -2, 2, 2, 2, -2],
                                    [-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 3, 3]],
                           'dim_suffix': 'linalg:det'},
               'parent':  'controller', 'dosmesh*': [[10, 10, 10],[12,12,12]],
               'dosmesh_suffix*':[30,"tt"],'atoms': 'Pd',
               'root': '/root/codes/matdb/tests/Pd', 'bandmesh': [13, 13, 13]}

    pgrid = ParameterGrid(db_info)

    for i in pgrid:
        folder = pgrid.to_str(i)
        assert i==pgrid.from_str(folder)

    pgrid.add((8.0,30),())
    assert len(pgrid) == 12
    assert pgrid[(8.0,30)] == ([2, 2, 2], [10, 10, 10])
    pgrid.pop((8.0,30))
    assert not ((8.0,30) in pgrid)
    assert not pgrid == ParameterGrid(db_info)
    assert pgrid == [(8.0, 'tt'), (2.0, 30), (2.0, 'tt'), (16.0, 30), (16.0, 'tt'), (32.0, 30), (32.0, 'tt'), (4.0, 30), (4.0, 'tt'), (27.0, 30), (27.0, 'tt')]
