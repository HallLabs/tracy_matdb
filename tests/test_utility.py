"""Tests the utility functions.
"""
import pytest
from os import path
def test_execute():
    """Tests the execution via shell subprocess in a different folder.
    """
    from matdb.utility import execute, reporoot
    sargs = ["pwd"]
    target = path.join(reporoot, "matdb")
    xres = execute(sargs, target, nlines=1)
    assert xres["output"][0].decode("UTF-8").strip() == target

    sargs = ["cat dummy-file"]
    target = path.join(reporoot, "matdb")
    xres = execute(sargs, target, nlines=1)
    assert len(xres["error"]) > 0

def test_cat(tmpdir):
    """Tests concatenation of multiple files.
    """
    from matdb.utility import cat, execute, reporoot
    files = [path.join(reporoot, "tests", "files", f)
             for f in ["A.txt", "B.txt"]]
    outfile = str(tmpdir.join("cat_C.txt"))
    cat(files, outfile)

    sargs = ["diff", "C.txt", outfile]
    xres = execute(sargs, path.join(reporoot, "tests/files"))
    assert len(xres["output"]) == 0

def test_symlink(tmpdir):
    """Tests symbolic linking of a file.
    """
    from matdb.utility import reporoot, symlink, execute
    target = path.join(reporoot, "__init__.py")
    source = str(tmpdir.join("symlink_init"))
    symlink(source, target)

    from os import readlink
    result = readlink(source)
    assert result == target

    dirsource = str(tmpdir.join("dummy-dir"))
    assert symlink(dirsource, reporoot) is None
