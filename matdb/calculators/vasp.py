"""Implements a `matdb` compatible subclass of the
:class:`ase.calculators.vasp.Vasp` calculator.
"""
from ase.calculators.vasp import Vasp
from os import path, stat
import mmap
from matdb.calculators.basic import AsyncCalculator
from matdb.utility import chdir
from matdb import msg

class AsyncVasp(AsyncCalculator, Vasp):
    """Represents a calculator that can compute material properties with VASP,
    but which can do so asynchronously.

    Args:
        atoms (quippy.Atoms): configuration to calculate properties for.
        folder (str): path to the directory to run this configuration in.
    """
    def __init__(self, atoms, folder, *args, **kwargs):
        super(Vasp, self).__init__(*args, **kwargs)
        
        self.folder = folder
        self.initialize(atoms)
        with chdir(folder):
            self.write_input(atoms)    

    def can_execute(self, folder):
        """Returns True if the specified folder is ready to execute VASP
        in.
        """
        if not path.isdir(folder):
            return False
        
        required = ["INCAR", "POSCAR", "KPOINTS", "POTCAR"]
        present = {}
        for rfile in required:
            target = path.join(folder, rfile)
            sizeok = stat(target).st_size > 25
            present[rfile] = path.isfile(target) and sizeok

        if not all(present.values()):
            for f, ok in present.items():
                if not ok:
                    msg.info("{} not present for VASP execution.".format(f), 2)
        return all(present.values())

    def can_cleanup(self, folder):
        """Returns True if the specified VASP folder has completed
        executing and the results are available for use.
        """
        if not path.isdir(folder):
            return False
    
        #If we can extract a final total energy from the OUTCAR file, we
        #consider the calculation to be finished.
        outcar = path.join(folder, "OUTCAR")
        if not path.isfile(outcar):
            return False

        line = None
        with open(outcar, 'r') as f:
            # memory-map the file, size 0 means whole file
            m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)  
            i = m.rfind('free  energy')
            if i > 0:
                # seek to the location and get the rest of the line.
                m.seek(i)
                line = m.readline()

        if line is not None:
            return "TOTEN" in line or "Error" in line
        else:
            return False
