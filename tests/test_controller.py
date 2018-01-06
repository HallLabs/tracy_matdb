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
                xpath = path.join(xroot, config, "dynmatrix", "W.1")
                #We want to make some testing assertions to ensure that the
                #stubs ran correctly.
                assert path.isfile(path.join(xpath, "CONTCAR"))
                assert path.isfile(path.join(xpath, ".matdb.module"))
                target = path.join(xpath, name)
                symlink(target, path.join(folder, dft))

def compare_dicts(dict1,dict2):
    """Compares two dictionaries to see if they are the same.
    """

    if sorted(dict1.keys()) != sorted(dict2.keys()):
        # print(dict1.keys(),dict2.keys())
        print(len(dict1.keys()),len(dict2.keys()))
        print([(a,b) for a, b in zip(dict1.keys(),dict2.keys()) if a!=b])
        return False

    for key in dict1:
        if isinstance(dict1[key],dict):
            res = compare_dicts(dict1[key],dict2[key])
            if not res:
                return False
            else:
                continue
        if not isinstance(dict1[key],six.string_types) and not np.allclose(dict1[key],dict2[key]):
            print(dict1[key],dict2[key])
            return False
        elif isinstance(dict1[key],six.string_types) and not dict1[key] == dict2[key]:
            print(dict1[key],dict2[key])
            return False

    return True
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

def test_ParameterGrid():
    """Tests the creation of a ParamaterGrid and it's functionality.
    """
    from matdb.database.controller import ParameterGrid

    db_info = {'phonopy': {'dim*': [[2, 0, 0, 0, 2, 0, 0, 0, 2], [0, 1, 1, 1, 0, 1, 1, 1, 0],
                                    [0, 2, 2, 2, 0, 2, 2, 2, 0], [-2, 2, 2, 2, -2, 2, 2, 2, -2],
                                    [-1, 1, 1, 1, -1, 1, 1, 1, -1], [3, 0, 0, 0, 3, 0, 0, 0, 3]],
                           'dim_suffix': {'func':'linalg:det','reshape':[3,3]}},
               'parent':  'controller', 'dosmesh*': [[10, 10, 10],[12,12,12]],
               'dosmesh_suffix*':[30,"tt"],'atoms': 'Pd',
               'root': '/root/codes/matdb/tests/Pd', 'bandmesh': [13, 13, 13]}

    pgrid = ParameterGrid(db_info)

    assert len(pgrid) == 12
    assert compare_dicts(pgrid['dos-tt-dim-16.00'],
                         {'phonopy': {'dim': [0, 2, 2, 2, 0, 2, 2, 2, 0]},
                          'dosmesh': [12, 12, 12], 'bandmesh': [13, 13, 13]})
    pgrid.pop('dos-tt-dim-16.00')
    assert not ('dos-tt-dim-16.00' in pgrid)
    assert not pgrid == ParameterGrid(db_info)


def test_get_grid():
    """Tests the get_grid method.
    """
    from matdb.database.controller import get_grid
    from numpy import array
    db_info = {"lattice*": ["fcc", "bcc", "hcp"], 
          "calculator": {"encut*": [700, 800, 900]},
          "normal": [1, 2, 3],
          "double": {"single": {"dog*": [range(9), list(np.diag([9,2,1])), list(np.diag([3,2,4]))]}},
          "encut_suffix*": [70,80,90],
          "dog_suffix": {"func": "linalg:det",
                         "reshape": [3 ,3]},
          "lattice_suffix": "{}"
        }

    test = get_grid(db_info)
    model = {'enc-70-dog-0.00-lat-bcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-0.00-lat-bcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-0.00-lat-bcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-18.00-lat-bcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-18.00-lat-bcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-18.00-lat-bcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-24.00-lat-bcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-24.00-lat-bcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-24.00-lat-bcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'bcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-0.00-lat-fcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-0.00-lat-fcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-0.00-lat-fcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-18.00-lat-fcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-18.00-lat-fcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-18.00-lat-fcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-24.00-lat-fcc': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-80-dog-24.00-lat-fcc': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-90-dog-24.00-lat-fcc': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'fcc',
  'normal': [1, 2, 3]},
 'enc-70-dog-0.00-lat-hcp': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-80-dog-0.00-lat-hcp': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-90-dog-0.00-lat-hcp': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [0, 1, 2, 3, 4, 5, 6, 7, 8]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-70-dog-18.00-lat-hcp': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-80-dog-18.00-lat-hcp': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-90-dog-18.00-lat-hcp': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([9, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 1])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-70-dog-24.00-lat-hcp': {'calculator': {'encut': 700},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-80-dog-24.00-lat-hcp': {'calculator': {'encut': 800},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]},
 'enc-90-dog-24.00-lat-hcp': {'calculator': {'encut': 900},
  'double': {'single': {'dog': [array([3, 0, 0]),
     array([0, 2, 0]),
     array([0, 0, 4])]}},
  'lattice': 'hcp',
  'normal': [1, 2, 3]}}

    assert compare_dicts(test,model)
