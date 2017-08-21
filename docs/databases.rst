`matdb` Database Types
======================

Depending on the needs of the potential, different configuration types
need to be generated for training. Each configuration type is referred
to loosely as a "database" and `matdb` supports multiple database
types:

1. Liquid from high-temperature sub-sampled molecular dynamics.
2. Phonon modulated structures by displacing atoms along eigenvectors
   of the dynamical matrix.

.. toctree::
   :maxdepth: 1
   :caption: Modules:

   database/phonon.rst
   database/liquid.rst
