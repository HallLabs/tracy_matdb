""" Tests the vegard interpolation in data.py
"""
import pytest
from matdb.data import vegard

def test_vegard():
    assert vegard(["H"],[3]) == 3.75
    assert vegard(["He","Li","Be"],[2,2,2]) == 3.11666666666666666667
    assert vegard(["B","C","N"],[1,3,5]) == 4.4038888888888888
