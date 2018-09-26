VASP
====

VASP Input
----------

The VASP calculator can take as input anything that could be put into
a `VASP` `INCAR` file. For full documentation of options please see
the `VASP` `manual`_.

Additionally the `matdb.yml` file alse takes the `kpoints` and
`potcar` fields as required input for the VASP calculator.

K-points
^^^^^^^^

The `kpoints` field is used to determine the k-point grid that will be
used for calculations. This field takes a `method` and any options
that need to be specified for that method. Currently supported methods
are:

1) `mueller` this method querries the k-point server [kps_ref]_ setup
   by Tim Mueller. To use this method you must have the `getKPoints`
   on your system path, the script can be downloaded `here`_. The only
   option accepted by this method is `mindistance` for which we
   recomend a value of ~50, the larger the number the denser the
   k-point grid that will be returned. An example of this option is:

.. code-block:: yaml

   kpoints:
     method: 'mueller'
     mindistance: 50

2) `gamma` this method will create a `KPOINTS` file with only a single
   k-point in it centered at gamma, i.e., the origin.

.. code-block:: yaml

   kpoints:
     method: 'gamma'

3) `kspacing` this option doesn't actually get written into the
   `kpoints` keyword block since it is placed into the `VASP` `INCAR`
   file. To use it simly specify the distance between k-points desired.

.. code-block:: yaml

   kspacing: 0.01

POTCAR
^^^^^^
   
This field given `matdb` the information needed to construct the
`POTCAR` file. The `potcars` section needs to contain the `directory`
field which contains the path to the VASP POTCARs you have, the `xc`
field which determines the type of `POTCAR` to use, i.e., "PBE" or
"GGA", the `versions` field and the `setups` field. The `versions`
field must supply a version for each chemical species in the system,
the versions are specified by the date found in the first line of the
`POTCAR` file::

    PAW_PBE Al 04Jan2001

The `setups` field specifies if a special `POTCAR` needs to be
used. For example VASP has the `Ag`, `Ag_pv`, `Ag_GW`, and `Ag_sv_GW`
potentials for Ag. Depending on the type and accuracy of the
calculation you are performing you would specify which of these you
potentials you wanted to use by supplying `` Ag: '_pv' `` in the
`setups` section. A complete example:

.. code-block:: yml
   
   potcars:
     directory: './tests/vasp'
   xc: 'PBE'
     versions:
       Ag: '02Apr2005'
       Pd: '28Jan2005'
     setups:
       Ag: '_pv'


A complete example of the VASP calculor `yml` section:

.. code-block:: yml

   calculator:
     name: Vasp
     nsw: 1
     pp: 'pbe'
     kpoints:
       method: 'mueller'
       mindistance: 30
     potcars:
       directory: './tests/vasp'
       xc: 'PBE'
       versions:
         Ag: '02Apr2005'
	 Pd: '04Jan2005'
       setups:
         Ag: '_pv'
		
VASP calculator
---------------

.. automodule:: matdb.calculators.vasp
   :synopsis: VASP subclass for interacting with the `ase` VASP
              calculator for job file creation.
   :members:      

      
.. _manual: https://cms.mpi.univie.ac.at/wiki/index.php/Category:INCAR
.. [kps_ref] Pandu Wisesa, Kyle A. McGill, and Tim Mueller Physcis Review B 93, 155109 (2016)
.. _here: http://muellergroup.jhu.edu/K-Points.html
