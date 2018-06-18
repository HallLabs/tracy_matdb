"""Tests the matdb.kpoints module functions."""
import numpy as np
def test_qpoint_finder():
    try:
        from matdb.kpoints import find_qpoints
    except ImportError:
        return
    from matdb.atoms import Atoms
    Si = Atoms(cell=np.identity(3), positions=[[0,0,0]], numbers=[14])
    supercell = (2,3,1)
    model = np.array([[ 0.        ,  0.        ,  0.        ],
                 [ 0.        ,  0.33333333,  0.        ],
                 [ 0.        , -0.33333333,  0.        ],
                 [ 0.5       ,  0.        ,  0.        ],
                 [ 0.5       ,  0.33333333,  0.        ],
                 [-0.5       , -0.33333333,  0.        ]])
    result = find_qpoints(Si, supercell)
    assert np.allclose(model, result)

def test_kpath():
    """Performs a simple test to extract the kpath for an alloy.
    """
    from matdb.atoms import Atoms
    from matdb.kpoints import parsed_kpath
    from matdb.utility import relpath
    filepath = relpath("tests/files/POSCAR-AgPd-50")
    at0 = Atoms(filepath, format="vasp")
    model = {'band': [[0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.5],
                      [0.25, 0.25, 0.25],
                      [0.0, 0.5, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.5, 0.5, -0.5],
                      [0.3126058432665583, 0.6873941567334416, -0.3126058432665583],
                      [0.0, 0.0, 0.0],
                      [-0.12521168653311662, 0.12521168653311662, 0.5],
                      [0.5, 0.5, -0.5]],
             'labels': ['\\Gamma',
                        'X',
                        'P',
                        'N',
                        '\\Gamma',
                        'M',
                        'S|S_0',
                        '\\Gamma|X',
                        'R|G',
                        'M']
    }
        
    labels, band = parsed_kpath(at0)
    assert labels == model["labels"]
    assert band == model["band"]
