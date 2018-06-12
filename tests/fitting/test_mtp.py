"""Tests fitting using legacy databases.
"""
import pytest

@pytest.fixture()
def mtpfit(tmpdir):
    from matdb.utility import relpath, reporoot
    from matdb.database import Controller
    from os import mkdir

    target = relpath("./tests/AgPd/matdb")
    dbdir = str(tmpdir.join("mtp"))
    mkdir(dbdir)

    cntrl = Controller(target, tmpdir=dbdir) 
    mtp = cntrl.trainers.fits['AgPd'].sequences['AgPd'].steps['mtp']

    return mtp

def test_prot_to_cfg():
    """Tests the conversion of a prototype structure to a cfg file.
    """
    from matdb.fitting.mtp import _prot_to_cfg
    from matdb.utility import _get_reporoot, touch
    from os import path, remove

    template_root = path.join(_get_reporoot(), "matdb", "templates")
    if not path.isdir(path.join(template_root, "uniqueUnaries")):
        with chdir(template_root):
            tarf = "prototypes.tar.gz"
            tar = tarfile.open(tarf, "r:gz")
            tar.extractall()
            tar.close()
    source = path.join(template_root, "uniqueUnaries", "A10_POSCAR.orig")
    species = ["Al"]
    relax_file = "relax.cfg"
    type_map = {0:0}
    root = "."
    min_atoms = 1
    max_atoms = None

    touch(relax_file)
    _prot_to_cfg(source, species, relax_file, type_map, root, min_atoms, max_atoms)

    assert path.isfile(relax_file)

    lc = 0
    with open(relax_file, "r") as f:
        for line in f:
            lc += 1

    print(lc)
    assert lc == 12
    remove(relax_file)

def test_init(mtpfit):

    """Tests the initialization of the mtp potenital fitter.
    """

    assert not mtpfit.ready()
