"""Tests the vasp calculator in matdb.
"""
import pytest
from os import path, mkdir, remove
import six
import numpy as np
from hashlib import sha1
from matdb.io import read
from matdb.utility import relpath

def test_get_calc_hashes():
    """Tests the get_calculator_hashes function."""
    from matdb.calculators.utility import get_calculator_hashes
    key = "matdb"
    target = relpath("./tests/AgPd/matdb")
    config = path.expanduser(path.abspath(target))
    if path.isabs(config):
        root, config = path.split(config)
    else:
        root, config = path.dirname(config), config
        
    configyml = read(root, config)
    bc = ""
    cp = {}
    get_calculator_hashes(key, configyml, bc, cp)

    assert "Vasp" in cp
    assert len(cp["Vasp"]) == 2
    for k ,v in cp["Vasp"].items():
        assert "tests/vasp" in v

def test_set_paths():
    """Tests the setting of the global paths."""
    from matdb.calculators.utility import paths, set_paths

    assert paths == {}
    target = relpath("./tests/AgPd/matdb")
    config = path.expanduser(path.abspath(target))
    if path.isabs(config):
        root, config = path.split(config)
    else:
        root, config = path.dirname(config), config
        
    configyml = read(root, config)

    set_paths(configyml)

    name = configyml["title"].strip().replace(' ', '_')
    namehash = str(sha1(name.encode("ASCII")).hexdigest())
    assert namehash in paths
    assert "Vasp" in paths[namehash]
    for k ,v in paths[namehash]["Vasp"].items():
        assert "tests/vasp" in v    
