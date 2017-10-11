"""This is a stub for vasp to enable unit testing. It creates empty files that
mimic the ones VASP actually creates when it executes.
"""
from matdb.utility import touch
files = ["CONTCAR", "WAVECAR", "CHGCAR"]
for fname in files:
    touch(fname)
