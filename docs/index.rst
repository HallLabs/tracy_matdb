.. Tracy Science documentation master file, created by
   sphinx-quickstart on Mon Apr  1 13:53:39 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Materials Database Generator API Documentation
==============================================

`matdb` is a python package for generating databases of configurations
needed for the building of interatomic potentials. While we are
working toward a fully-automated approach, in the meantime databases
can be constructed by building an appropriate :doc:`matdbyml`. The
following links document the internal API and give some examples.

.. toctree::
   :maxdepth: 1
   :caption: Input files and high-level control:
	     
   matdbyml.rst
   workflow.rst

.. toctree::
   :maxdepth: 1
   :caption: Control Scripts:
	     
   scripts.rst
   execution.rst

.. toctree::
   :maxdepth: 1
   :caption: Sub-packages and Modules:

   calculator.rst
   databases.rst
   fitting.rst
   kpoints.rst
   io.rst
   data.rst
   base.rst
   atoms.rst
   descriptors.rst
   queries.rst
   utility.rst
   transforms.rst
   exceptions.rst
   logs.rst
   msg.rst
   phonons.rst
   plotting.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
