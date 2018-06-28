"""Tests the legacy database interface.
"""
import pytest
from matdb.database.legacy import LegacyDatabase
from matdb.utility import relpath
from os import mkdir, path
from matdb.atoms import AtomsList, Atoms
import numpy as np

@pytest.fixture
def phondb(tmpdir):
    """Returns a legacy database for the split-up phonon-50 AgPd database.
    """
    splits = {
        "A": .5,
        "B": .75
    }
    root = str(tmpdir.join("legacy"))
    if not path.isdir(root):
        mkdir(root)
    folder = relpath("./tests/data/legacy")

    return LegacyDatabase("AgPd-50", root, None, splits, folder, "p-50-*.xyz",
                          "ph", energy="dft_energy", force="dft_force", virial="dft_virial")

@pytest.fixture
def rendb(tmpdir):
    """Returns a legacy database for the *renamed* phonon-50 AgPd database.
    """
    splits = {
        "0": .3,
        "1": .6
    }
    root = str(tmpdir.join("legacy"))
    if not path.isdir(root):
        mkdir(root)
    folder = relpath("./tests/data/legacy")

    return LegacyDatabase("R-50", root, None, splits, folder, "r-50-*.xyz",
                          "re", energy="energy", force="force", virial="virial",
                          limit=80)

def test_split(phondb):
    """Tests splitting of available data from a legacy database.
    """
    supers = {}
    
    phondb.split()
    for s, p in phondb.splits.items():
        tfile = path.join(phondb.train_file(s).format(s))
        hfile = path.join(phondb.holdout_file(s).format(s))
        sfile = path.join(phondb.super_file(s).format(s))

        tal = AtomsList(tfile)
        hal = AtomsList(hfile)
        sal = AtomsList(sfile)
        supers[s] = sal

        assert len(tal) == int(np.ceil(150*p))
        assert len(hal) == int(np.ceil((150-len(tal))*p))
        assert len(sal) == 150-len(tal)-len(hal)

        assert path.isfile(path.join(phondb.root, "{}-ids.pkl".format(s)))

    #Now, make sure that we return quickly if split is called again.
    phondb.split()

    #Remove one of the files so that we can trigger reading from existing ids.
    from os import remove
    for s, p in phondb.splits.items():
        sfile = path.join(phondb.super_file(s).format(s))
        remove(sfile)
        
    phondb.split()

    for s, p in phondb.splits.items():
        sfile = path.join(phondb.super_file(s).format(s))
        sal = AtomsList(sfile)
        assert sal == supers[s]

def test_errors(tmpdir, rendb):
    """Tests raising of exceptions for badly named parameters.
    """
    root = str(tmpdir.join("legacy"))
    if not path.isdir(root):
        mkdir(root)
    folder = relpath("./tests/data/legacy")
    splits = None
    with pytest.raises(ValueError):
        LegacyDatabase("errors", root, None, splits, folder, "r-50-*.xyz")
    with pytest.raises(ValueError):
        LegacyDatabase("errors", root, None, splits, folder, "r-50-*.xyz",
                       energy="dft_energy")
    with pytest.raises(ValueError):
        LegacyDatabase("errors", root, None, splits, folder, "r-50-*.xyz",
                       energy="dft_energy", force="dft_force")
    
def test_merge(phondb):
    """Tests merger of the databases into a single DB.
    """
    combined = AtomsList(phondb._dbfile)
    assert len(combined) == 150
    assert  combined[20].params["config_type"] == "ph"

    root = path.dirname(phondb.root)
    newdb = LegacyDatabase("AgPd-50", root, None, phondb.splits, phondb.folder,
                           "p-50-*.xyz", "ph")
    assert newdb.atoms == phondb.atoms
    
def test_rename(rendb):
    """Tests renaming of properties to meet `matdb` conventions.
    """
    first = Atoms(rendb._dbfile)
    assert "ref_energy" in first.params
    assert "ref_force" in first.properties
    assert "ref_virial" in first.params
    assert  first.params["config_type"] == "re"

    al = AtomsList(rendb._dbfile)
    assert len(al) == 80

def test_ran_seed(tmpdir):
    """Tests that the random seed gets set by the controller properly.
    """

    from matdb.utility import relpath
    from matdb.database import Controller

    target = relpath("./tests/general/ran_seed")
    dbdir = str(tmpdir.join("general"))
    mkdir(dbdir)
    cntrl = Controller(target, dbdir)
    
    splits = {
        "A": .5,
        "B": .75
    }
    root = str(tmpdir.join("legacy"))
    if not path.isdir(root):
        mkdir(root)
    folder = relpath("./tests/data/legacy")
    
    result = LegacyDatabase("AgPd-50", root, cntrl, splits, folder, "p-50-*.xyz",
                            "ph")

    assert result.ran_seed == 100
