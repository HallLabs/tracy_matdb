"""Tests the generation of interactive plots using the custom matdb
plugin. Basically, we just create a plot with two points, but each of those
points has a set of three random plots generated to add detail to the point.
"""
import pytest
import numpy as np
from os import path, mkdir
import matplotlib.pyplot as plt

def test_pointdetail(tmpdir):
    """Tests creation of a point-detail plot of a sine-wave.
    """
    from matdb.plotting.matd3 import PointDetailImage
    folder = str(tmpdir.join("plots"))
    mkdir(folder)
    
    x = np.linspace(0, np.pi, 25)
    y = np.sin(x)
    image = PointDetailImage(x, y, folder=folder, index=0, name="sine")
    assert image.url == "sine-0.png"
    assert path.isfile(path.join(folder, image.filename))

    from matdb.utility import relpath
    img2 = PointDetailImage(x, y, index=0, name="sine", base64=True)
    with open(relpath("tests/plotting/sine-base64.dat")) as f:
        model = f.read()
    assert img2.url == model
    plt.close('all')
    
def test_html(tmpdir):
    """Tests the creation of the HTML package folder with custom point detail
    plots.
    """
    #Generate two sets of three images, sin, cos and tan.
    from matdb.plotting.matd3 import PointDetailImage as PDI
    folder = str(tmpdir.join("html"))
    mkdir(folder)

    def plot_tuple(imgtype):
        x0 = x[imgtype];
        x1 = x0 + np.pi/4.
        ufunc = getattr(np, imgtype)
        return [(x0, ufunc(x0)), (x1, ufunc(x1))]
    
    x = {
        "sin": np.linspace(0, np.pi, 25),
        "cos": np.linspace(np.pi/4, 3*np.pi/4, 25),
        "tan": np.linspace(-np.pi/4, 0, 25)
    }
    xy = {k: plot_tuple(k) for k in x}
    images = {k: [PDI(*v, folder=folder, index=vi, name=k)
                  for vi, v in enumerate(xy[k])]
              for k in x}

    data = {
        "david": {
            "location": (0.25, 0.25),
            "index": 0
        },
        "goliath": {
            "location": (0.5, 0.5),
            "index": 0
        }
    }
    for k, detail in images.items():
        data["david"][k] = detail[0]
        data["goliath"][k] = detail[1]

    from matdb.plotting.matd3 import html
    titles = {"sin": "Sine", "cos": "Cosine", "tan": "Tangent"}
    subplot_kw = {
        "title": "David vs. Goliath",
        "xlabel": "Strength",
        "ylabel": "Beauty"
    }
    plot_kw = {
        "alpha": 0.3
    }
    html(data, folder, titles=titles, subplot_kw=subplot_kw, plot_kw=plot_kw)
    #Make sure all the files were created. Since it is a UI, it has to be tested
    #manually.
    ifiles = ["sin-0.png", "sin-1.png", "tan-0.png",
              "tan-1.png", "cos-0.png", "cos-1.png"]
    for ifile in ifiles:
        assert path.isfile(path.join(folder, ifile))

    assert path.isfile(path.join(folder, "imagepoint.js"))
    assert path.isfile(path.join(folder, "imagepoint.css"))
    assert path.isfile(path.join(folder, "index.html"))
    plt.close('all')    
