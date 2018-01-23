"""Tests the enumeration group interface.
"""
import pytest
from matdb.database.enumerated import Enumerated
from matdb.utility import relpath
from os import mkdir, path, symlink
import quippy
import numpy as np

@pytest.fixture()
def AgPd(tmpdir):
    from matdb.utility import relpath
    from matdb.database.controller import Controller
    from os import mkdir, symlink, remove

    target = relpath("./tests/AgPd/matdb")
    dbdir = str(tmpdir.join("agpd_db"))
    mkdir(dbdir)
    
    result = Controller(target, dbdir)
    return result

def test_setup(AgPd):
    """Tetsts the setup of the enumerated database.
    """

    assert not AgPd.collections['enumerated']['enumerated'].steps['enum'].is_setup()
    
    AgPd.setup()

    dbs = ["Enum/enumerated/lat-{}".format(i) for i in (1,2)]

    folders = {
        "__files__": ["compute.pkl","euids.pkl","jobfile.sh","enum.out",
                      "lattice.in","index.json"],
        "E.1": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.2": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.3": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.4": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.5": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.6": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.7": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.8": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.9": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        },
        "E.10": {
            "__files__": ["INCAR","PRECALC","POSCAR"]
        }
    }

    from matdb.utility import compare_tree
    for db in dbs:
        dbfolder = path.join(AgPd.root,db)
        compare_tree(dbfolder,folders)    

    assert AgPd.collections['enumerated']['enumerated'].steps['enum'].is_setup()

    # test the euid and index creation for the entire database.
    assert path.isfile(path.join(AgPd.root,"Enum/enumerated/euids.pkl"))
    assert path.isfile(path.join(AgPd.root,"Enum/enumerated/index.json"))

    enum = AgPd.collections['enumerated']['enumerated'].steps['enum']
    assert len(enum.index) == 20
    assert len(enum.euids) == 20

def test_functions(AgPd):
    """Tests the enumerated specific functions of the database.
    """

    enum = AgPd.collections['enumerated']['enumerated'].steps['enum']

    #test the lattice settings for the system
    enum._get_lattice("sc")
    assert np.allclose(enum.lattice,[[1,0,0],[0, 1, 0],[0, 0, 1]])
    assert enum.lattice_name == "sc"
    enum._get_lattice([[1,0,0],[0,1.5,0],[0,0.5,2]])
    assert np.allclose(enum.lattice,[[1,0,0],[0,1.5,0],[0,0.5,2]])
    assert enum.lattice_name == "custom"
    enum._get_lattice("hcp")
    assert np.allclose(enum.lattice,[[1,0,0],[.5, 0.866025403784439, 0],[0, 0, 1.6329931618554521]])
    assert enum.lattice_name == "hcp"

    #test the basis setting function
    enum._get_basis(None)
    assert np.allclose(enum.basis, [[0,0,0],[0.5,0.28867513459,0.81649658093]])
    enum._get_basis([[1,2,3],[4,5,6]])
    assert np.allclose(enum.basis, [[1,2,3],[4,5,6]])

    # test the concentrations setting function
    assert not enum.conc_res
    assert enum.concs is None
    enum._get_concs([[0,2,6],[0,2,6]])
    assert enum.conc_res
    assert np.allclose(enum.concs,[[0,2,6],[0,2,6]])

    # test the arrow concentration setting function
    assert not enum.arrow_res
    assert enum.arrows is None
    enum._get_arrows([1,1])
    assert np.allclose(enum.arrows,[1,1])
    assert enum.arrow_res

