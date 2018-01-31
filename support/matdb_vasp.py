"""This is a stub for vasp to enable unit testing. It creates empty
files that mimic the ones VASP actually creates when it executes. We
don't include the OUTCAR because we want to be able to simulate VASP
failing as well, if the OUTCAR is needed then we simply copy an
existing OUTCAR from `tests/data`. We're copying the POSCAR to the
CONTCAR so that the cleanup methods can grab the 'relaxed' atomic
positions.
"""
from matdb.utility import touch, symlink
from os import path
files = ["CONTCAR", "WAVECAR", "CHGCAR"]
for fname in files:
    touch(fname)

if path.isfile("POSCAR"):
    symlink("CONTCAR","POSCAR")
