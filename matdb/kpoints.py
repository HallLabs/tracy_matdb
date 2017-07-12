"""Functions for converting k-points between different units,
requesting special paths in the BZ and so on.
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
