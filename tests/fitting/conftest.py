"""Provides fixtures for a pre-configured database that includes output files so
that a database is ready for testing the trainers.
"""
import pytest
from os import path

@pytest.fixture()
def Pd_db(tmpdir):
    from matdb.utility import relpath, reporoot
    from matdb.database.controller import Controller
    from os import mkdir

    target = relpath("./tests/Pd/matdb")
    dbdir = str(tmpdir.join("pd_db"))
    mkdir(dbdir)
    
    #We need to copy the POSCAR over from the testing directory to the temporary
    #one.
    from shutil import copy
    POSCAR = relpath("./tests/Pd/POSCAR")
    copy(POSCAR, dbdir)
    
    Pd = Controller(target, dbdir)
    Pd.setup()
    
    #First, we need to copy the FORCE_SETS and total_dos.dat files so that we
    #don't have to recompile those (they are tested elsewhere).
    from matdb.utility import symlink    
    troot = path.join(reporoot, "tests", "data", "Pd", "dynmatrix")
    files = ["FORCE_SETS", "total_dos.dat", "mesh.yaml"]
    for seq in Pd.find("Pd.phonon-*.dynmatrix"):
        for filename in files:
            target = path.join(seq.root, "phonopy", filename)
            source = path.join(troot, "{0}__{1}".format(filename, seq.parent.name))
            symlink(target, source)

    Pd.cleanup()
    Pd.setup()

    files = ["OUTCAR", "output.xyz"]
    seq = Pd["Pd.modulate.modulations"]
    for i in range(1, 6):
        key = "M.{0:d}".format(i)
        for filename in files:
            target = path.join(seq.root, key, filename)
            source = path.join(troot, "{0}__{1}".format(filename, key))
            symlink(target, source)        
        
    return Pd
