Distortion Database Documentation
=================================

Example yml File Entry
**********************

.. code-block:: yaml

   steps:
     - type: "distortion.Distortion"
       seeds: ["POSCAR:PdAg25","POSCAR:PdAgp0","POSCAR:PdAg75"]
       rattle: 0.1
       volume_factor: 1.1
       cov_diag: [1, 2, 3]

Distortion Database
*******************

.. automodule:: matdb.database.distortion
   :synopsis: Database sub-classes for implementing databases of
              lattice and atomic distortions.
   :members:      

