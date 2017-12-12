"""This is a stub for the sbatch script on the supercomputer.
"""
from matdb.utility import touch
from os import getcwd
import sys

if len(sys.argv) >1 and sys.argv[1][-3:] == ".sh":
    from matdb.utility import execute
    execute(["bash",sys.argv[1]],getcwd())
    sys.stdout.write("Submitted batch job {}\n".format(abs(hash(sys.argv[1]))))
else:
    sys.stderr.write("Failed to submit\n")


