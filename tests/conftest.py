import pytest

@pytest.fixture(scope="session", autouse=True)
def stubs(request, tmpdir_factory):
    #We need to hack $PATH so that `vasp` and `module` point to the unit testing
    #stubs that we created.
    from os import path, pathsep, environ
    from matdb.utility import which, symlink
    stubs = {
        "vasp": "matdb_vasp.py",
        "module": "matdb_module.py",
        "sbatch": "bash",
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
    xres = execute(["module", "load", "mkl/*"], stubpath)
    assert path.isfile(path.join(stubpath, ".matdb.module"))
    xres = execute(["sbatch", "-c", "pwd"], stubpath)
    assert xres["output"][0].strip() == stubpath

    touch(path.join(stubpath, "PRECALC"))
    xres = execute(["getKPoints"], stubpath)
    assert path.isfile(path.join(stubpath, "KPOINTS"))
    xres = execute(["vasp"], stubpath)
    assert path.isfile(path.join(stubpath, "CONTCAR"))
    
    def restore():
        paths = environ["PATH"].split(pathsep)
        #Make sure that we actually added our temporary directory to the path.
        assert paths[0] == stubpath
        environ["PATH"] = pathsep.join(paths[1:])
    request.addfinalizer(restore)
    
    return stubpath