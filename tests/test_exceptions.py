"""Tests the custom matdb exceptions.
"""
import pytest

def test_versions():
    """Tests the versions exception.
    """
    from matdb.exceptions import VersionError

    with pytest.raises(VersionError):
        raise VersionError("Test of version error.")

def test_species():
    """Tests the Species exception.
    """
    from matdb.exceptions import SpeciesError

    with pytest.raises(SpeciesError):
        raise SpeciesError("Test of version error.")

def test_logic():
    """Tests the logic errors.
    """
    from matdb.exceptions import LogicError

    with pytest.raises(LogicError):
        raise LogicError("Test of version error.")

def test_mlp():
    """Tests the logic errors.
    """
    from matdb.exceptions import MlpError

    with pytest.raises(MlpError):
        raise MlpError("Test of version error.")
