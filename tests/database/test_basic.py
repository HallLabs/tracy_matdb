"""Tests the basic database class structure from which all others
inherit.
"""
import pytest
@pytest.fixture(scope="module", autouse=True)
def Pd(request):
    """Returns a :class:`matdb.database.basic.Database` using Pd as a
    seed configuration.
    """
    from matdb.database.basic import Database
    
