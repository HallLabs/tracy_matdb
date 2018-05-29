"""This contains the algorithms needed to decompress a compressed crystal.
"""

import numpy as np

def decompress(prim, basis, types, hnf_vals):
    """Decompresses the crystal back into it's original form.

    Args:
        prim (list): the primitive lattice vectors as rows of a matrix.
        basis (list): the atomic basis vectors as rows of a matrix.
        types (list): list of integers for the atomic species.
        hnf_vals (list): integer hnf entries.

    Returns:
        The new crystal lattice vectors, atomic basis and atomic types.
    """

    hnf = [[hnf[0], 0, 0], [hnf[1], hnf[2], 0], [hnf[3], hnf[4], hnf[5]]]
    lat_vecs = np.transpose(np.matmul(np.transpose(prim), hnf))

    vol_fact = hnf[0]*hnf[2]*hnf[5]

    new_basis = []
    new_types = []
    prim = np.array(prim)
    for a in range(hnf_vals[0]):
        for b in range(hnf_vals[2]):
            for c in range(hnf_vals[5]):
                #calculate the vector that will point to a new atom in
                #the basis by taking a linear combination of the
                #primitive cell vectors.
                add_vec = prim[0]*a + prim[1]*b + prim[c]*c
                for old_t, old_b in zip(types, basis):
                    new_b = list(np.array(old_b)+add_vec)
                    new_basis.append(bring_into_cell(new_b, lat_vecs))
                    new_types.append(old_t)

    if vol_fact*len(basis) != len(new_basis):
        raise ValueError("Error occured in decompression.")
                    
    return lat_vecs, new_basis, new_types
                    
def bring_into_cell(point, vecs):
    """Brings the point into the cell defined by vecs.
    
    Args:
        point (list): the points listed in cartesian coordinates.
        vecs (list): the lattice vectors of the crystal as rows of a matrix.

    Returns:
        the point translated into the unit cell still in cartesion coordinates.
    """
    
    cart_to_latt = np.linalg.inv(np.transpose(vecs))
    lat_to_cart = np.transpose(vecs)

    eps = 1E-3
    point = np.matmul(cart_to_latt, point)
    while any(point-eps > 1) or any(point<-eps):
        for i, p in enumerate(point):
            if p >= 1-eps:
                point[i] -= 1
            elif p < -eps:
                point[i] += 1

    point = list(np.matmul(lat_to_cart, point))

    return point
