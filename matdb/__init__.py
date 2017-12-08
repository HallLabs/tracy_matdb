__version__ = [1, 0, 4]
from matdb import base
from matdb import msg
import matplotlib
import os

if os.system != "nt": #pragma: no cover
    matplotlib.use("Agg")
