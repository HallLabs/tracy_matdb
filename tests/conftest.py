import pytest

@pytest.fixture(scope="session", autouse=True)
def stubs(request, tmpdir_factory):
    #We need to hack $PATH so that `vasp` and `module` point to the unit testing
    #stubs that we created.
    from os import path, pathsep, environ
    from matdb.utility import which, symlink, touch
    stubs = {
        "vasp": "matdb_vasp.py",
        "module": "matdb_module.py",
        "sbatch": "matdb_sbatch.py",
        "getKPoints": "matdb_getkpoints.py"
    }

    #Validate the existence of installed stub scripts in the python path.
    binpaths = {}
    for name, xstub in stubs.items():
        binpath = which(xstub)
        if binpath is None:
            emsg = "Cannot find `{0}` stub for unit testing.".format(name)
            raise EnvironmentError(emsg)
        binpaths[name] = binpath

    #Create a new temporary directory as part of the unit test; create symlinks
    #in the directory to the stubs and hack the path to look in that directory
    #first.
    stubpath = str(tmpdir_factory.mktemp("stubs"))
    for name, xstub in stubs.items():
        lpath = path.join(stubpath, name)
        symlink(lpath, binpaths[name])
        
    environ["PATH"] = stubpath + pathsep + environ["PATH"]
    from matdb.utility import execute, touch
    xres = execute(["./module", "load", "mkl/*"], stubpath)
    assert path.isfile(path.join(stubpath, ".matdb.module"))
    xres = execute(["./sbatch", "-c", "pwd"], stubpath)
    assert xres["output"] == []
    # assert xres["error"][0].strip() == "Failed to submit"
    symlink(stubpath+"/sbatch.sh",stubpath+"/sbatch")
    xres = execute(["./sbatch", stubpath+"/sbatch.sh"], stubpath)
    temp = xres["output"][-1].strip().split()[0:3]
    assert ' '.join(temp) == "Submitted batch job"
    # assert xres["error"] == []

    touch(path.join(stubpath, "PRECALC"))
    xres = execute(["./getKPoints"], stubpath)
    assert path.isfile(path.join(stubpath, "KPOINTS"))
    xres = execute(["./vasp"], stubpath)
    assert path.isfile(path.join(stubpath, "CONTCAR"))
    
    def restore():
        paths = environ["PATH"].split(pathsep)
        #Make sure that we actually added our temporary directory to the path.
        assert paths[0] == stubpath
        environ["PATH"] = pathsep.join(paths[1:])
    request.addfinalizer(restore)
    
    return stubpath

@pytest.fixture
def paper():
    """Returns a search query that mimics the one shown in the paper
    for electrically-insulating heat sinks.
    """
    import aflow
    import aflow.keywords as kw
    result = aflow.search(batch_size=20
        ).select(kw.agl_thermal_conductivity_300K
        ).filter(kw.Egap > 6).orderby(kw.agl_thermal_conductivity_300K, True)

    #Let's pre-fill the responses from the saved JSON files so that the tests
    #run faster *and* so that the results are predictable.
    import json
    result._N = 912
    
    n = -1
    with open("tests/files/aflow/data.json") as f:
        response = json.loads(f.read())
        result.responses[n] = response

    return result[0:20]
