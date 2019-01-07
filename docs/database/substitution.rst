Substitution Database Documentation
===================================

Example yml File Entry
**********************

.. code-block:: yaml

   steps:
     - type: "substitution.Substitution"
       seeds: ["POSCAR:PdAg25"]
       stoich:
         - [0.25, 0.75, 2]
	 - [0.4, 0.6, 5]
	 - [0.5, 0.5, 3]

This is an example of a binary system for which we want the seed
configuration to be altered so as to produce 2 configurations with a
25% A and 75% B, 5 with 40% A and 60% B, and 3 with 50% A and 50% B
resulting in a total of 10 configurations from the initial seed
configuration. If multiple seeds are provided then this would yield
10*(the number of seed configurations).

.. automodule:: matdb.database.substitution
   :synopsis: Database sub-classes for implementing databases of
              lattice with atomic substitutions throughout the
              lattice.
   :members:      

