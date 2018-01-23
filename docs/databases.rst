`matdb` Database Types
======================

Depending on the needs of the potential, different configuration types
need to be generated for training. Each configuration type is referred
to loosely as a "Group" and `matdb` supports multiple group types:

1. Liquid from high-temperature sub-sampled molecular dynamics.
2. Phonon modulated structures by displacing atoms along eigenvectors
   of the dynamical matrix.
3. Enumerated structures selected from an enumerated list of
   symmetrically distinct configurations.
4. AFLOW which allows the querying of the `aflowlib.org` database.
   
.. toctree::
   :maxdepth: 1
   :caption: Modules:

   database/phonon.rst
   database/liquid.rst
   database/enuerated.rst
   database/aflow.rst
