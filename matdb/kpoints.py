"""To automatically produce the phonon band structure plots, we need to
have an appropriate special path in the BZ. `matdb` interfaces with
`materialscloud.org` to extract their recommended path in the BZ.
"""
from os import path, getcwd, chdir, system
import numpy as np

def kpath(atoms):
    """Returns a list of the special k-points in the BZ at which the
    phonons should be sampled.

    Args:
        atoms (matdb.Atoms): structure to get the path for.

    Returns:

        tuple: `(path, points)`, where `path` is a list of *names* of the special points
        and `points` is a dict with the same names as keys and :class:`numpy.ndarray` as
        values.
    """
    from seekpath.hpkot import get_path
    s = (atoms.cell, atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    kdict = get_path(s)

    names = [kdict["path"][0][0]]
    for ki in range(len(kdict["path"]))[1:]:
        s0, e0 = kdict["path"][ki-1]
        s1, e1 = kdict["path"][ki]

        if e0 == s1:
            names.append(s1)
        else:
            names.append((e0, s1))
    names.append(kdict["path"][-1][1])
    
    pts = {k: np.array(v) for k, v in kdict["point_coords"].items()}
    return (names, pts)   

def _gamma_only(target, atoms):
    """Saves a KPOINTS file with gamma-point only in the specified
    folder.

    Args:
        target (str): path to the directory to create the KPOINTS file
          in.
    """
    kpoints = path.join(target, "KPOINTS")
    with open(kpoints, 'w') as f:
        f.write("""Gamma-point only
1 ! one k-point
rec ! in units of the reciprocal lattice vector
0 0 0 1 ! 3 coordinates and weight
""")

def _mueller(target,atoms,mindistance=None):
    """Gets the Mueller style k-points from Mueller's server.
    
    Args:
        target(str): path to the directory to create the KPOINTS file.
        rmin (float): the cutoff for the k-point density by Mueller's
            metric.
    """

    precalc = path.join(target,"PRECALC")
    with open(precalc,"w+") as f:
        if mindistance is not None:
            f.write("INCLUDEGAMMA=AUTO \n"
                    "MINDISTANCE={}".format(mindistance))
        else:
            raise ValueError("'mindistiance' must be provided for Mueller k-points.")
    cur_dir = getcwd()
    chdir(target)
    system("getKPoints")
    chdir(cur_dir)
       
def custom(target, key, atoms=None):
    """Creates the KPOINT file with the specified custom configuration
    in `target`.

    Args:
        target (str): path to the directory to create the KPOINTS file
          in.
        key (str): one of ['gamma'], specifies which custom KPOINTS
          file to generate (see note above).
        atoms (quippy.Atoms): atoms object to generate KPOINTS for.
    """
    select = {
        "gamma": _gamma_only,
        "mueller": _mueller
    }
    if key["method"] in select:
        method_args = key.copy()
        del method_args["method"]
        select[key["method"]](target, atoms,**method_args)
    else:
        emsg = "'{}' is not a valid key for custom k-points."
        raise ValueError(emsg.format(key))
