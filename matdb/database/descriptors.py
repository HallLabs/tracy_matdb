"""Methods for analyzing coverage in descriptor space for 2+3+... GAPs. It is
impossible to produce a good 2-body and 3-body potential if the data provided to
GAP is sparse between those parts of descriptor space where Hessian information
is available. To remedy this, we construct heuristics to determine when we have
reached the threshold of descriptor coverage to provide good behavior for all
pieces of the overall GAP potential.
"""
import quippy
import numpy as np
import pandas as pd

def d2b(atoms, cutoff):
    """Returns histogram data for the specified *list* of atoms within the given
    cutoff.

    Args:
        atoms (list): A list of atoms objects
        cutof (float): The distance to cutoff the contributions at.

    Returns:
        pandas.dataframe: of ["N","D","T"] as the columns.
    """
    D2b = quippy.Descriptor("distance_Nb order=2 cutoff={0:.2f}".format(cutoff))
    rs = []
    for ai in atoms:
        ai.set_cutoff(cutoff)
        ai.calc_connect()
        r2b = D2b.calc(ai)
        rs.append(r2b["descriptor"].flatten())

    R2b = np.hstack(rs)
    result = np.histogram(R2b)
    stacked = np.vstack((result[0], result[1][:-1])).T
    result = pd.DataFrame(stacked, columns=["N", "D"])
    result["T"] = "d2"
    return result

def d3b(atoms, cutoff):
    """Returns histogram data for the specified *list* of atoms within the given
    cutoff.

    Args:
        atoms (list): A list of atoms objects
        cutof (float): The distance to cutoff the contributions at.

    Returns:
        pandas.dataframe: of ["N","D","T"] as the columns.
    """
    D3b = quippy.Descriptor("distance_Nb order=3 cutoff={0:.2f}".format(cutoff))
    ds = []
    for ai in atoms:
        ai.set_cutoff(cutoff)
        ai.calc_connect()
        r3b = D3b.calc(ai)
        ds.append(r3b["descriptor"])

    R3b = np.vstack(ds)
    Ns = []
    Ds = []
    Ts = []
    for i in range(3):
        hist = (np.histogram(R3b[:,i]))
        Ns.append(hist[0])
        Ds.append(hist[1][:-1])
        Ts.append(np.array(["d3{0:d}".format(i) for j in range(len(hist[0]))]))

    contents = {
        "N": np.hstack(Ns),
        "D": np.hstack(Ds),
        "T": np.hstack(Ts)
    }

    return pd.DataFrame(contents)

def dNb(atoms, cutoffs):
    """Constructs the descriptors for all the specified N-body descriptor
    cutoffs and the specified atoms list or object.

    Args:
        atoms (quippy.AtomsList): of atoms to calculate descriptors for.
        cutoffs (dict): keys are `int` number of bodies in the descriptors;
          values are the cutoffs in Angstroms.
    """
    dmap = {
        2: d2b,
        3: d3b
    }
    dfs = []
    for n, cutoff in cutoffs.items():
        dfs.append(dmap[n](atoms, cutoff))

    return pd.concat(dfs)
