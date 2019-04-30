Vacancy Database Documentation
==============================

Example yml File Entry
**********************

.. code-block:: yaml

   steps:
     - type: "vacancy.Vacancy"
       seeds: ["POSCAR:PdAg25","POSCAR:PdAgp0","POSCAR:PdAg75"]
       vac_per_atom: 0.2

Vacancy Database
****************

.. automodule:: matdb.database.vacancy
   :synopsis: Database sub-classes for implementing databases of
              site vacancies in lattices.
   :members:      

