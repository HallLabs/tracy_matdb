Calculators
===========

YML File Specifications
-----------------------

Within the `matdb.yml` file the only keywords that can be used for
each calculator are the `name` and `exec_path` fields. The `name`
field is should simply supply the name of the calculator desired, most
of the time the given name, if otherwise it will be specified in that
calculators documentation. The `exec_path` is optional and provides
the path to the executable for the calculator, i.e., `vasp.x` or
`pw.x`. If `exec_path` is not supplied then it is assumed that the
executable is on your system path. All other fields are calculator
dependent and the individual calculators documentation should be
consulted.

Supported Calculators
---------------------

`matdb` was designed such that it can be extended to interface with
any `ase` calculator. Details on haw to add a specific calculator to
`matdb` can be found below. 

Calculators currently supported are: (**!!Is AFLOW a calculator?? More of a library..**)

.. toctree::
   :maxdepth: 1

   calculators/qe.rst
   calculators/vasp.rst
   calculators/aflow.rst


.. toctree::
   :caption: Calculator Utility and other modules
   :maxdepth: 1

   calculators/basic.rst
   calculators/utility.rst
