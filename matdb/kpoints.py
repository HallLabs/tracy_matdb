"""To automatically produce the phonon band structure plots, we need to
have an appropriate special path in the BZ. `matdb` interfaces with
`materialscloud.org` to extract their recommended path in the BZ.
"""
from os import path
import numpy as np

def kpath(poscar):
    """Returns a list of the special k-points in the BZ at which the
    phonons should be sampled.

    Args:
        poscar (str): full path to the POSCAR file that contains the
          structural information.

    Returns:

        tuple: `(path, points)`, where `path` is a list of *names* of the special points
        and `points` is a dict with the same names as keys and :class:`numpy.ndarray` as
        values.
    """
    import requests
    files = { "structurefile": open(path.abspath(path.expanduser(poscar)), 'rb') }
    data = { "fileformat": "vasp" }
    url = "http://www.materialscloud.org/tools/seekpath/process_structure/"
    r = requests.post(url, data=data, files=files)

    from bs4 import BeautifulSoup
    parsed_html = BeautifulSoup(r.text, "html5lib")

    from json import loads
    kptcode = parsed_html.body.find("code", attrs={"id": "rawcodedata"}).get_text()
    kptcode = kptcode.replace(u'\xa0', u' ')
    kdict = loads(kptcode)

    names = [kdict[u"path"][0][0]]
    for ki in range(len(kdict[u"path"]))[1:]:
        s0, e0 = kdict[u"path"][ki-1]
        s1, e1 = kdict[u"path"][ki]

        if e0 == s1:
            names.append(s1)
        else:
            names.append((e0, s1))
    names.append(kdict[u"path"][-1][1])
    
    pts = {k: np.array(v) for k, v in kdict[u"kpoints"].items()}
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
        "gamma": _gamma_only
    }
    if key in select:
        select[key](target, atoms)
    else:
        emsg = "'{}' is not a valid key for custom k-points."
        raise ValueError(emsg.format(key))
