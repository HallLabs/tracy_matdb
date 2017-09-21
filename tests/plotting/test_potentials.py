"""Tests the generalized plotting of all potential-related plots that are
helpful in deciding how to proceed with a GAP fit.

.. note:: The :func:`~matdb.plotting.potentials.generate` function has the
  capability to produce plots for *all* the plot types. So, we just test it
  directly and ignore tests for the others since they are tested indirectly.
"""
import pytest
import quippy
from os import path
def test_generate(tmpdir):
    """Tests generation of all plots supported by the potential plotter.
    """
    from matdb.plotting.potentials import generate
    folder = tmpdir.join("plotting")
    #We will plot some figures for a binary system using a simple Lennard-Jones
    #potential. Since we are doing binary, we expect 3 dimer plots and 4 trimer
    #plots. So, we end up with the following expected file names:
    files = ["energy-0.png", "force-0.png", "virial-0.png", "evsv-0.png",
             "AgPd-0.png", "AgAg-0.png", "PdPd-0.png", "Ag2Pd.png" "AgPd2.png",
             "Ag3.png", "Pd3.png"]
    plots = "efvodt"

    #Now, grab a test atoms list for AgPd that has lots of variation, then get a
    #LJ potential and set its parameters to be suitable for AgPd.
    from matdb.utility import relpath
    atoms = quippy.AtomsList(relpath("tests/files/AgPd.xyz"))
    pot = quippy.Potential("IP LJ", param_filename=relpath("tests/files/LJ-AgPd.xml"))    
    generate(plots, pot, atoms, str(folder))

    for fname in files:
        assert path.isfile(str(folder.join(fname)))    
