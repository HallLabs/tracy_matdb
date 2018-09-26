Quantum Espresso
=================================

Quantum Espresso Input
----------------------

The Quantum Espresso (QE) calculator can be used by supplying the
`calculator` `name` field with 'qe' rather than the full name. Input
fields for this calculator are `input_data`, `kpoints`, and
`potcars`. The `input_data` field most of the parameters that will be
written into the QE controll file and needs to be broken up into
subfields `control`, `system`, `electrons`, `ions`, `cell`
`atomic_species`, `constraints`, `occupations`, `atomic_forces`. Valid
entries within these subfields can be found in the QE `manual`_.

.. note:: The `cell_parameters`, `k_points`, and `atomic_positions`
          fields mentioned in the QE `manual`_ are not used because
          they are provided by other entries into the calculator.

Other fields needed for the QE calculator are `kpoints` and `potcars`.
	  
K-points
^^^^^^^^

The `kpoints` field is used to determine the grid of k-points to be
used by QE. For QE the field takes `method` and `offset` entries in
addition to parameters based off of which method is used. Methods
supported for QE include:

1) `MP` for Monkhorst-Pack k-point grids. This method requires the
   additional `divisions` parameter be specified, the divisions takes
   3 integers which indicate how aften to split the reciprocal cell in
   each direction to form the k-point grid.

.. code-block:: yml

   kpoints:
     method: 'MP'
     divisions: 3 3 3

2) `kspacing` which is used to indicate the space desired between
   k-points. This method takes a single additional entry of `spacing`.

.. code-block:: yml

   kpoints:
     method: 'kspacing'
     spacing: 0.1

The `offset` field is optional for this calculator. If you use it
specify the offset vector you would like, if not supplied then no
offset will be used.
     
POTCARS
^^^^^^^

This section is used to construct the potentials section of the
control file. For QE the fields that need to be supplied are
`directory` which provides the path to where the potential files are
located on your system, `versions` and `potentials`. The `potentials`
section requires the name of the potential file for each element be
supplied. The `versions` fields requires the version numbers for the
specific potentials being used are supplied. For example the versions
from the potential `Al.pbe-n-kjpaw_psl.1.0.0.UPF`::

  <UPF version="2.0.1">
    <PP_INFO>
      Generated using "atomic" code by A. Dal Corso  v.5.1
      Author: ADC
      Generation date: 10Oct2014
      Pseudopotential type: PAW
      Element: Al
      Functional:  SLA  PW   PBX  PBC

Are '2.0.1' and '.5.1'. Different potential files store the versions
in slightly different places in the front matter but two are always
supplied and both should be put in the `YML` file.

A full example of the `matdb.yml` file section for a QE calculator is
below:

.. code-block:: yml
   
   calculator:
     name: 'qe'
     exec_path: '/home/tracy/q-e-qe-6.3/bin/pw.x'
     input_data:
       control:
         calculation: 'relax'
	 prefix: 'test'
	 outdir: 'output'
     kpoints:
       method: 'MP'
       divisions: 3 3 3
       offset: 0.5 0.5 0.5
     potcars:
       directory: '/home/tracy/pslibrary/pz/PSEUDOPOTENTIALS/'
       potentials:
         Al: 'Al.pbe-n-kjpaw_psl.1.0.0.UPF'
	 Pd: 'Pd_ONCV_PBE-1.0.upf'
       versions:
         Al: '2.0.1' '.5.1'
	 Pd: '2.0.1' '2.1.1'

.. automodule:: matdb.calculators.qe
   :synopsis: Quantum Espresso subclass for interacting with the `ase`
              Quantum Espresso calculator for job file creation.
   :members:      

.. _manual: https://www.quantum-espresso.org/Doc/INPUT_PW.html
