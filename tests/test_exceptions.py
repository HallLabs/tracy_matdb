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
    """Tests the versions exception.
    """
    from matdb.exceptions import VersionError

    with pytest.raises(SpeciesError):
        raise SpeciesError("Test of version error.")
