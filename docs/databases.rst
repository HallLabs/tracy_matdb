`matdb` Database Types
======================

Depending on the needs of the potential, different configuration types
need to be generated for training. Each configuration type is referred
to loosely as a "Group" and `matdb` supports multiple group types:

1. Enumerated structures selected from an enumerated list of
   symmetrically distinct configurations.
2. AFLOW which allows the querying of the `aflowlib.org` database.
3. Active which contains structures being used by `mtp` active learning.
4. Distortion which distorts the atoms within the seed configurations.
5. Substitution which changes the stoichiometry of the seed
   configuration.
6. Vacancy which introduces vacancies into the atomic basis of the
   seed configuration.
   
.. toctree::
   :maxdepth: 1
   :caption: Modules:

   database/enumerated.rst
   database/aflow.rst
   database/active.rst
   ..
      database/distribution.rst
      database/substitution.rst
      database/vacancy.rst
