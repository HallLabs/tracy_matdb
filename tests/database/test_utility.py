# -*- coding: utf-8 -*-
"""Tests the utility functions for the databases.
"""
import pytest
import numpy as np

def test_swap_column():
    """Tests the swap column function.
    """
    from matdb.database.utility import swap_column

    hnf = np.array([[0,1,0],[1,1,0],[1,0,1]])
    b = np.array([[1,2,0],[3,1,2],[0,0,1]])

    new_hnf, new_b = swap_column(hnf, b, 0)

    assert np.allclose(new_hnf, [[1,0,0],[1,1,0],[0,1,1]])
    assert np.allclose(new_b, [[2,1,0],[1,3,2],[0,0,1]])

def test_hnf():
    """Tests the conversion of integer matrices to hnf form.
    """
    from matdb.database.utility import hermite_normal_form

    n = [[-1, -1, -2], [2, 6, -1], [-7, 0, 5]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [66, 85, 111]]
    true_b = [[-8, -10, -13], [3, 4, 5], [2, 3, 4]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[-8, 5, -1], [3, -5, 5], [-5, -9, 9]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [52, 137, 208]]
    true_b = [[-5, -13, -20], [-9, -24, -37], [-6, -16, -25]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[6, -6, -5], [5, 9, 7], [2, -6, -6]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [68, 34, 96]]
    true_b = [[-2, -1, -3], [47, 24, 67], [-59, -30, -84]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[0, 7, -5], [10, 5, 8], [2, -6, 0]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [0, 1, 0], [296, 416, 462]]
    true_b = [[52, 73, 81], [-32, -45, -50], [-45, -63, -70]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[10, -1, -5], [-1, 7, 2], [-2, -5, 0]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [2, 3, 0], [2, 1, 3]]
    true_b = [[14, 12, 11], [-6, -5, -5], [29, 25, 23]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    n = [[1, 0, 0], [3, 2, 0], [1, -1, -2]]
    
    hnf, b = hermite_normal_form(n)

    true_hnf = [[1, 0, 0], [1, 2, 0], [0, 1, 2]]
    true_b = [[1, 0, 0],[-1, 1, 0], [1, -1, -1]]

    assert np.allclose(hnf, true_hnf)
    assert np.allclose(b, true_b)

    with pytest.raises(ValueError):
        hermite_normal_form([[1,0,0],[1,0,0],[0,0,0]])
