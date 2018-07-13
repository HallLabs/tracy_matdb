# -*- coding: utf-8 -*-
"""Tests the descriptors functions.
"""
import pytest
from matdb.atoms import Atoms

def test_soap():
    """Tests the soap fingerprint.
    """
    from matdb.descriptors import soap
    atSi = Atoms("Si8",positions=[[0,0,0],[0.25,0.25,0.25],[0.5,0.5,0],[0.75,0.75,0.25],
                                  [0.5,0,0.5],[0.75,0.25,0.75],[0,0.5,0.5],[0.25,0.75,0.75]],
                 cell=[5.43,5.43,5.43])

    assert soap(atSi) == []
