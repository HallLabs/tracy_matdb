Manual Database Documentation
=============================

To use this database the user must place every desired configuration
of atoms they want in the database into the `SEED` directory of the
`root` directory specified in the `matdb.yml` file.

Example yml File Entry
**********************

.. code-block:: yaml

  steps:
  - type: 'simple.Manual'
    seeds: ["POSCAR:PdAg25","POSCAR:PdAgp0","POSCAR:PdAg75"]
   
Manual Database
***************

.. automodule:: matdb.database.simple
   :synopsis: Database sub-classes for implementing databases of
              user provided configurations.
   :members:      

