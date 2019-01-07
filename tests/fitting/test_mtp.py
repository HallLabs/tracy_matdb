"""Tests fitting using legacy databases.
"""
import pytest
import json

import warnings
warnings.filterwarnings("ignore")

@pytest.fixture()
def mtpdb(tmpdir):
    from matdb.utility import relpath, reporoot, copyonce
    from matdb.database import Controller
    from os import mkdir, path

    target = relpath("./tests/mtp/CoWV.yml")
    dbdir = str(tmpdir.join("mlp_tests"))
    mkdir(dbdir)
    copyonce(target, path.join(dbdir, "matdb.yml"))
    target = path.join(dbdir,"matdb")

    cntrl = Controller(target, tmpdir=dbdir)

    return cntrl

def test_prot_to_cfg():
    """Tests the conversion of a prototype structure to a cfg file.
    """
    from matdb.fitting.mtp import _prot_to_cfg
    from matdb.utility import _get_reporoot, touch, chdir
    from os import path, remove
    import tarfile

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

    assert lc == 12

    source = path.join(template_root, "uniqueBinaries", "b31_POSCAR.orig")
    species = ["Al", "Cu"]
    relax_file = "relax.cfg"
    type_map = {0:0, 1:1}
    root = "."
    min_atoms = 1
    max_atoms = 4
    _prot_to_cfg(source, species, relax_file, type_map, root, min_atoms, max_atoms)
    
    lc = 0
    with open(relax_file, "r") as f:
        for line in f:
            lc += 1

    assert lc == 12

    remove(relax_file)

def test_create_to_relax():
    """Tests the creation of the to_relax_file.
    """
    from matdb.fitting.mtp import create_to_relax
    from matdb.utility import _get_reporoot
    from os import remove, path
    import shutil

    template_root = path.join(_get_reporoot(), "matdb", "templates")
    if path.isdir(path.join(template_root, "uniqueUnaries")):
        shutil.rmtree(path.join(template_root, "uniqueUnaries"), ignore_errors=True)
        shutil.rmtree(path.join(template_root, "uniqueBinaries"), ignore_errors=True)
        shutil.rmtree(path.join(template_root, "uniqueTernaries"), ignore_errors=True)
        
    
    setup_args = {}
    target = "to-relax.cfg"
    setup_args["phenum_args"] = {"config":"t", "species": ["Al", "Cu"],  "structures": "all",
                          "debug": False, "example": False, "displace":0.0, 
                          "mink":True, "outfile":target, "verbose": None, 
                          "rattle":0.0, "mapping": None}
    setup_args["species"] = ["Al", "Cu"]
    setup_args["crystals"] = ["fcc"]
    setup_args["min_atoms"] = 1
    setup_args["max_atoms"] = 2
    setup_args["root"] = "."

    create_to_relax(setup_args)

    struct_count = 0
    with open(target,"r") as f:
        for line in f:
            if "BEGIN_CFG" in line:
                struct_count += 1

    assert struct_count == 4
    remove(target)

    setup_args["crystals"] = ["fcc", "prototypes"]
    create_to_relax(setup_args)

    struct_count = 0
    with open(target,"r") as f:
        for line in f:
            if "BEGIN_CFG" in line:
                struct_count += 1
    
    assert struct_count == 138
    remove(target)

def test_is_float():
    """Tests the is_float function.
    """

    from matdb.fitting.mtp import _is_float

    assert _is_float("0.0")
    assert _is_float(1.0)
    assert _is_float(1)
    assert _is_float("1")
    assert not _is_float("stuff")

def test_relax_ini(mtpdb):
    """Tests the setting of the relax_ini settings.
    """

    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    relaxargs = {}
    mtpfit._set_relax_ini(relaxargs)
    assert mtpfit.relax_ini["calc_efs"] == "TRUE"
    assert mtpfit.relax_ini["efs_ignore"] == "FALSE"
    assert mtpfit.relax_ini["active_learn"] == "TRUE"
    assert mtpfit.relax_ini["fit_setting"] == "FALSE"
    assert mtpfit.relax_ini["site_weight"] == "0.0"
    assert mtpfit.relax_ini["energy_weight"] == "1.0"
    assert mtpfit.relax_ini["force_weight"] == "0.001"
    assert mtpfit.relax_ini["stress_weight"] == "0.0001"
    assert mtpfit.relax_ini["extrap_threshold"] == "2.0"
    assert mtpfit.relax_ini["threshold_break"] == "10.0"

    relaxargs["calc-efs"] = "FALSE"
    relaxargs["efs-ignore"] = "TRUE"
    relaxargs["active-learn"] = "FALSE"
    relaxargs["fit"] = "TRUE"
    relaxargs["site-weight"] = "2.0"
    relaxargs["energy-weight"] = "2.0"
    relaxargs["force-weight"] = "0.1"
    relaxargs["stress-weight"] = "0.1"
    relaxargs["threshold"] = "0.1"
    relaxargs["threshold-break"] = "0.1"
    mtpfit._set_relax_ini(relaxargs)
    assert mtpfit.relax_ini["calc_efs"] == "FALSE"
    assert mtpfit.relax_ini["efs_ignore"] == "TRUE"
    assert mtpfit.relax_ini["active_learn"] == "FALSE"
    assert mtpfit.relax_ini["fit_setting"] == "TRUE"
    assert mtpfit.relax_ini["site_weight"] == "2.0"
    assert mtpfit.relax_ini["energy_weight"] == "2.0"
    assert mtpfit.relax_ini["force_weight"] == "0.1"
    assert mtpfit.relax_ini["stress_weight"] == "0.1"
    assert mtpfit.relax_ini["extrap_threshold"] == "0.1"
    assert mtpfit.relax_ini["threshold_break"] == "0.1"

    relaxargs["calc-efs"] = "1.0"
    relaxargs["stress-weight"] = "stuff"
    mtpfit._set_relax_ini(relaxargs)
    assert mtpfit.relax_ini["calc_efs"] == "TRUE"
    assert mtpfit.relax_ini["stress_weight"] == "0.0001"

def test_set_local_attributes(mtpdb):
    """Tests the initialization of the mtp potenital fitter.
    """

    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    mtpargs = {}
    mtpfit._set_local_attributes(mtpargs)
    assert mtpfit.relax_args == {}
    assert mtpfit.train_args == {}
    assert mtpfit.relax_min_atoms == 1
    assert mtpfit.relax_max_atoms is None
    assert mtpfit.grade_args == {}
    assert mtpfit.select_args == {}
    assert mtpfit.use_mpi
    assert mtpfit.ran_seed == 0
    assert not mtpfit.use_unrelaxed
    assert mtpfit.crystals_to_relax == ["sc", "fcc", "bcc", "hcp", "prototypes"]

    mtpargs["relax_ini"] = {"threshold_break": "1.0"}
    mtpargs["relax"] = {"stress-tolerance": 0.1}
    mtpargs["ran_seed"] = 1
    mtpargs["train"] = {"max-iter": 10}
    mtpargs["smallest_relax_cell"] = 2
    mtpargs["largest_relax_cell"] = 3
    mtpargs["calc-grade"] = {"nbh-weight": 10}
    mtpargs["select"] = {"nbh-weight": 10}
    mtpargs["use_mpi"] = "False"
    mtpargs["use_unrelaxed"] = "True"
    mtpargs["crystals_to_relax"] = ["sc", "bcc"]
    mtpfit._set_local_attributes(mtpargs)
    assert mtpfit.relax_args == {"stress-tolerance": 0.1}
    assert mtpfit.train_args == {"max-iter": 10}
    assert mtpfit.relax_min_atoms == 2
    assert mtpfit.relax_max_atoms == 3
    assert mtpfit.grade_args == {"nbh-weight": 10}
    assert mtpfit.select_args == {"nbh-weight": 10}
    assert not mtpfit.use_mpi
    assert mtpfit.ran_seed == 1
    assert mtpfit.use_unrelaxed
    assert mtpfit.crystals_to_relax == ["sc", "bcc"]

    mtpargs["smallest_relax_cell"] = 4
    mtpargs["largest_relax_cell"] = 3
    with pytest.raises(ValueError):
        mtpfit._set_local_attributes(mtpargs)

    mtpargs["smallest_relax_cell"] = 4
    mtpargs["largest_relax_cell"] = -1
    with pytest.raises(ValueError):
        mtpfit._set_local_attributes(mtpargs)

    mtpargs["smallest_relax_cell"] = -1
    mtpargs["largest_relax_cell"] = None
    with pytest.raises(ValueError):
        mtpfit._set_local_attributes(mtpargs)
    
def test_init(mtpdb):
    """Tests the __init__ function.
    """
    from os import path

    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    with open(path.join(mtpfit.root, "status.txt"), "w+") as f:
        f.write("done 1 0")

    assert mtpfit.iter_status is None
        
    mtpargs = {}
    mtpfit.__init__(controller=mtpfit.controller, dbs=mtpfit._dbs,
                    execution=mtpfit.execution, split=mtpfit.split, root=mtpfit.root,
                    parent=mtpfit.parent, dbfilter=mtpfit.dbfilter, **mtpargs)

    assert mtpfit.iter_status == "done"
    assert not mtpfit.ready()


def test_make_input_files(mtpdb):
    """Tests the creation of various input files.
    """
    from os import path, mkdir
    from matdb.utility import touch, _get_reporoot, copyonce
    from matdb.atoms import Atoms

    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    touch(path.join(mtpfit.root, "to-relax.cfg"))
    mtpfit._make_pot_initial()

    assert path.isfile(path.join(mtpfit.root, "pot.mtp"))

    with open(path.join(mtpfit.root, "pot.mtp"), "r") as f:
        for line in f:
            assert not "{{" in line and not "}}" in line

    mtpfit._make_relax_ini()
    assert path.isfile(path.join(mtpfit.root, "relax.ini"))
    with open(path.join(mtpfit.root, "relax.ini"), "r") as f:
        for line in f:
            assert not "{{" in line and not "}}" in line

    mtpfit._setup_to_relax_cfg()
    assert path.isfile(path.join(mtpfit.root, "to_relax.json"))

    with open(path.join(mtpfit.root, "to_relax.json")) as f:
        to_relax_dict = json.load(f)

    assert to_relax_dict["phenum_args"] == {"config":"t", "species":mtpfit.species,
                                            "structures": "all", "debug": False, 
                                            "example": False, "displace":0.0, "mink":True,
                                            "outfile":path.join(mtpfit.root,"to-relax.cfg"),
                                            "verbose": None, "rattle":0.0, 
                                            "mapping": None}
    assert to_relax_dict["crystals"] == mtpfit.crystals_to_relax
    assert to_relax_dict["species"] == mtpfit.species
    assert to_relax_dict["min_atoms"] == mtpfit.relax_min_atoms
    assert to_relax_dict["max_atoms"] == mtpfit.relax_max_atoms
    assert to_relax_dict["root"] == mtpfit.root
    
    mtpfit.species = ["Co", "Ni", "Ti", "W", "Li"]

    assert mtpfit._setup_to_relax_cfg() is None

def test_train_setup(mtpdb):
    """Tests the mtp training setup.
    """
    #test make_train_cfg
    #first build the database directory.
    from matdb.utility import copyonce, _get_reporoot
    from os import path, mkdir, remove
    from matdb.atoms import Atoms

    seed_root = path.join(mtpdb.root, "seed")
    mkdir(seed_root)
    templates = path.join(_get_reporoot(), "tests", "mtp", "training")
    for i in range(1,11):
        target = path.join(seed_root, "vasp.{0}".format(i))
        source = path.join(templates, "POSCAR{0}".format(i))
        copyonce(source, target)    
    
    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    cntrl_root = mtpdb.root
    db_root = path.join(cntrl_root, "Manual", "test.manual")
    mtpdb.setup()

    
    files = ["OUTCAR{}", "CONTCAR{}"]
    for i in range(1,11):
        target = path.join(db_root, "vasp.{0}".format(i), "S1.1")
        for f in files:
            copyonce(path.join(templates, f.format(i)), path.join(target, f.format('')))
    mtpdb.extract()
    mtpfit._make_train_cfg(1)
    assert path.isfile(path.join(mtpfit.root, "train.cfg"))

    new_configs = []
    for i in range(1,11):
        source = path.join(templates, "POSCAR{0}".format(i))
        new_configs.append(Atoms(source))

    mtpfit.active.add_configs(new_configs, 2)
    mtpfit.active.setup()
    act_root = path.join(mtpdb.root, "Active", "active.active")
    files = ["OUTCAR{}", "CONTCAR{}"]
    for i in range(1,11):
        target = path.join(act_root, "Ac.{0}".format(i))
        for f in files:
            copyonce(path.join(templates, f.format(i)), path.join(target, f.format('')))

    mtpfit.active.last_iteration = None
    mtpfit._make_train_cfg(2)
    assert path.isfile(path.join(mtpfit.root, "train.cfg"))

    mtpfit.active.last_iteration = None
    remove(mtpfit.active.iter_file)
    with pytest.raises(IOError):
        mtpfit._make_train_cfg(2)

def test_templates(mtpdb):
    """Tests the creation of the execution templates.
    """

    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    #train template
    template = "mpirun -n 72 mlp train pot.mtp train.cfg > training.txt"

    assert mtpfit._train_template() == template

    mtpfit.use_mpi = False
    mtpfit.train_args = {"curr-pot-name":"name1", "valid-cfgs": "name2.cfg",
                         "bfgs-conv-tol":0.001}

    template = "mlp train pot.mtp train.cfg --bfgs-conv-tol=0.001 > training.txt"
    assert mtpfit._train_template() == template
    
    #calc-grade template
    template = "mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg"
    assert mtpfit._calc_grade_template() == template

    mtpfit.grade_args = {"force-weight":10, "mvs-filename":"name1"}

    template_parts = ["mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg",
                      "--force-weight=10"]
    mtp_out = mtpfit._calc_grade_template()
    for part in template_parts:
        assert part in mtp_out

    mtpfit.use_mpi = True
    #relax template
    template_parts = ["mpirun", "-n 72", "mlp relax relax.ini",
                      "--cfg-filename=to-relax.cfg", "--save-relaxed=relaxed.cfg",
                      "--log=relax_log.txt", "--save-unrelaxed=unrelaxed.cfg"]
    mtp_out = mtpfit._relax_template()
    for part in template_parts:
        assert part in mtp_out

    mtpfit.use_mpi = False
    mtpfit.relax_args = {"log": "name1", "bfgs-wolfe_c1": True, "max-step": 10}
    template_parts = ["mlp relax relax.ini", "--cfg-filename=to-relax.cfg",
                      "--save-relaxed=relaxed.cfg", "--log=relax_log.txt",
                      "--save-unrelaxed=unrelaxed.cfg", "--max-step=10",
                      "--bfgs-wolfe_c1"]
    mtp_out = mtpfit._relax_template()
    for part in template_parts:
        assert part in mtp_out

    #select step
    template_parts = ["mlp select-add pot.mtp train.cfg candidate.cfg new_training.cfg",
                      "--force-weight=10"]
    mtpfit.select_args = {"mvs-filename":"name1", "force-weight":10}
    mtp_out = mtpfit._select_template()
    for part in template_parts:
        assert part in mtp_out

def test_update_split(mtpdb):
    """Tests the updating of the split params.
    """
    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    mtpfit.split = None
    mtpfit._update_split_params()
    assert mtpfit.params["split"] == 1

def test_status(mtpdb):
    """Tests the status function.
    """

    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    mtpfit.status()

    result = mtpfit.status(printed=False)

    assert not result["trained"]
    assert not result["jobfile"] 

def test_command_functions(mtpdb):
    """All the commands are iterative, so we'll have to test them one
    after the other.
    """

    #first build the database directory and the initial database of
    #structures.
    from matdb.utility import copyonce, _get_reporoot, touch
    from os import path, mkdir, remove, rename
    from matdb.atoms import Atoms

    seed_root = path.join(mtpdb.root, "seed")
    mkdir(seed_root)
    templates = path.join(_get_reporoot(), "tests", "mtp", "training")
    for i in range(1,11):
        target = path.join(seed_root, "vasp.{0}".format(i))
        source = path.join(templates, "POSCAR{0}".format(i))
        copyonce(source, target)    
    
    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    cntrl_root = mtpdb.root
    db_root = path.join(cntrl_root, "Manual", "test.manual")
    mtpdb.setup()

    files = ["OUTCAR{}", "CONTCAR{}"]
    for i in range(1,11):
        target = path.join(db_root, "vasp.{0}".format(i), "S1.1")
        for f in files:
            copyonce(path.join(templates, f.format(i)), path.join(target, f.format('')))
    mtpdb.extract()

    touch(path.join(mtpfit.root, "new_training.cfg"))
    touch(path.join(mtpfit.root, "relaxed.cfg"))

    # With the initial database setup we can now run the mtp commands.

    # First test the initial training setup script
    cmd_template = mtpfit.command()
    target = path.join(mtpfit.root, "status.txt")
    bc = path.join(mtpfit.root, "status.txt_train")
    copyonce(target,bc)
    assert path.isfile(target)

    with open(target, "r") as f:
        line = f.read()

    assert line.strip() == "relax_setup 1 0"
    assert path.isfile(path.join(mtpfit.root, "relax.ini"))
    assert path.isfile(path.join(mtpfit.root, "pot.mtp"))
    assert cmd_template == mtpfit._train_template

    # Now we setup for the test of the relaxation_setup
    copyonce(path.join(mtpfit.root, "pot.mtp"), path.join(mtpfit.root, "Trained.mtp_"))
    cmd_template = mtpfit.command()

    assert cmd_template == "matdb_mtp_to_relax.py"
    with open(target, "r") as f:
        line = f.read()
    assert line.strip() == "relax 1 0"

    rename(bc, target)
    copyonce(target,bc)
    copyonce(path.join(mtpfit.root, "pot.mtp"), path.join(mtpfit.root, "Trained.mtp_"))
    touch(path.join(mtpfit.root, "unrelaxed.cfg"))
    mtpfit.use_unrelaxed = True
    cmd_template = mtpfit.command()

    assert cmd_template == mtpfit._relax_template()
    with open(target, "r") as f:
        line = f.read()
    assert line.strip() == "select 1 0"

    mtpfit.use_unrelaxed = False
    rename(bc, target)
    copyonce(path.join(mtpfit.root, "pot.mtp"), path.join(mtpfit.root, "Trained.mtp_"))
    touch(path.join(mtpfit.root, "candidate.cfg"))
    cmd_template = mtpfit.command()
    assert cmd_template == mtpfit._relax_template()
    with open(target, "r") as f:
        line = f.read()
    assert line.strip() == "select 1 0"    

    # Now we can test the select step
    for i in range(72):
        touch(path.join(mtpfit.root, "candidate.cfg_{0}".format(i)))

    cmd_template = mtpfit.command()
    assert path.isfile(path.join(mtpfit.root, "candidate.cfg"))
    with open(target, "r") as f:
        line = f.read()
    assert line.strip() == "add 1 0"    

#This test doesn't work because something is wrong with the jinja template import
@pytest.mark.skip()
def test_command_functions2(mtpdb):
    # Finally we test the add step
    from matdb.utility import copyonce, _get_reporoot, touch
    from os import path, mkdir, remove, rename
    from matdb.atoms import Atoms

    mtpfit = mtpdb.trainers.fits['CoWV_mtp'].sequences['CoWV_mtp'].steps['mtp']
    target = path.join(mtpfit.root, "status.txt")
    with open(target, "w+") as f:
        f.write("add 1 0")
    touch(path.join(mtpfit.root, "relaxed.cfg"))
    touch(path.join(mtpfit.root, "new_training.cfg"))
    cmd_template = mtpfit.command()
    assert cmd_template == ''

    assert len(mtpfit.active.last_iteration) == 10
    assert path.isfile(path.join(mtpfit.root, "relaxed.cfg_iter_1"))
    assert path.isfile(path.join(mtpfit.root, "new_training.cfg_iter_1"))
    
    templates = path.join(_get_reporoot(), "tests", "mtp", "training")
    act_root = path.join(mtpdb.root, "Active", "active.active")
    files = ["OUTCAR{}", "CONTCAR{}"]
    for i in range(1,11):
        target = path.join(act_root, "Ac.{0}".format(i))
        for f in files:
            copyonce(path.join(templates, f.format(i)), path.join(target, f.format('')))

    target = path.join(mtpfit.root, "status.txt")
    with open(target, "r") as f:
        line = f.read()
    assert line.strip() == "done 1 0 10"    
    cmd_template = mtpfit.command()
    assert cmd_template == mtpfit._train_template
