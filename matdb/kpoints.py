"""To automatically produce the phonon band structure plots, we need to
have an appropriate special path in the BZ. `matdb` interfaces with
`materialscloud.org` to extract their recommended path in the BZ.
"""
from os import path, getcwd, chdir, system
import numpy as np

def find_qpoints(atoms, supercell):
    """Finds the `q` points that would be "pinned down" by a calculation of the
    given supercell.

    Args:
        atoms (matdb.Atoms): *primitive* atoms object.
        supercell (numpy.ndarray): of shape `(3, 3)` or a list/tuple of shape 3
          or 9.
    """
    from kgridgen import kpointgeneration
    A = np.asarray(atoms.cell, order='F')
    B = np.asarray(atoms.make_supercell(supercell).cell, order='F')
    n = int(np.linalg.det(B)/np.linalg.det(A))    
    qpoints = np.zeros((n, 3), order='F')
    kpointgeneration.findqpointsinzone(A, B, n, qpoints)
    return qpoints

def parsed_kpath(atoms):
    """Gets the special path in the BZ for the structure with the specified
    atoms object and then parses the results into the format required by the package
    machinery.
    Args:
        atoms (:class:`matdb.atoms.Atoms`): a matdb `Atoms` object.
    Returns:
        tuple: result of querying the materialscloud.org special path
        service. First term is a list of special point labels; second is the
        list of points corresponding to those labels.
    """
    ktup = kpath(atoms)
    band = []
    labels = []
    names, points = ktup

    def fix_gamma(s):
        return r"\Gamma" if s == "GAMMA" else s
    
    for name in names:
        #Unfortunately, the web service that returns the path names returns
        #GAMMA for \Gamma, so that the labels need to be fixed.
        if isinstance(name, tuple):
            key = name[0]
            labels.append("{}|{}".format(*map(fix_gamma, name)))
        else:
            key = name
            labels.append(fix_gamma(name))

        band.append(points[key].tolist())

    return (labels, band)

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
