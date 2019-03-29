#!/usr/bin/python

"""This is a unit testing stub for the `getKPoints` executable that constructs
`KPOINTS` files for the DFT runs. It just creates a `KPOINTS` file with a
comment to trick `matdb` into thinking that the folder is ready to execute. This
avoids the network latency associated with the requests so that tests run
quickly.
"""
from os import path
if not path.isfile("PRECALC"):
    raise EnvironmentError("Cannot find PRECALC file to initiate KPOINTS.")

from numpy.random import randint, random
N = randint(1, 10)
lines = [
    "MATDB Dummy KPOINTS file",
    str(N)
]
for i in range(N):
    lines.append(str(random()))
    
with open("KPOINTS", 'w') as f:
    f.write('\n'.join(lines))    
