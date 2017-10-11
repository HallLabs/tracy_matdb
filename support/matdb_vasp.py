"""This is a stub for vasp to enable unit testing. It copies a `vasprun.xml`
output file from a location specified by a `.matdb.json` file in the directory
where it runs to that directory.

The `.matdb.json` file should be created by the `module` stub that mimics
`module load` and `module unload` capabilities.
"""
from os import path

if not path.isfile(".matdb.json"):
    raise FileNotFoundError("Cannot find .matdb.json file required by stub.")

import json
with open(".matdb.json") as f:
    data = json.load(f)

if "vasprun.xml" in data:
    from shutil import copy
    current = path.abspath('.')
    copy(data["vasprun.xml"], current)
