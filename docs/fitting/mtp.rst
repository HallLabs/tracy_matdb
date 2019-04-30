Fitting for MTP Potentials
==========================

MTP (Moment Tensor Potential) use a polynomial basis to expand the
interatomic potential and offer fast calculation times for relaxing
cells and calculating the energies and forces of an atomic
configuration. The pontentials are also adaptable in terms of
accuracy, in that they can be expanded to be more accurate but then
they take a little longer to run.

.. warning:: MTP fits require unrelaxed calculations to be
             effecetive. Please ensure the database you are using does
             not contain relaxations.

	     
YML file specifications
-----------------------

To fit an MTP potential the use `mtp.MTP` in the `type` field of the
fitting step. The MTP section of the yml file requires the `species`
section be provided and can take `relax_ini`, `relax`, `ran_seed`,
`train`, `smallest_relax_cell`, `largest_relax_cell`, `calc_grade`,
`select`, `use_mpi`, `use_unrelaxed`, and `crystals_to_relax`.

relax_ini
*********

The `relax_ini` section can contain the following fields and values:

- **calc-efs**: (bool) If `TRUE` then `mlp` will calculate the energy,
   force, and stress of the configuration of atoms. Input can be
   either `TRUE` or `FALSE`, default value `TRUE`.
- **efs-ignore**: (bool) If `TRUE` then energy, force and stress
   calculations will be ignored in fitting situations. Input can be
   either `TRUE` or `FALSE`, default value `FALSE`.
- **active-learn**: (bool) Turns on/off the selection, i.e., active
   learning engine. Default value `TRUE`.
- **fit**: (bool) Turns fitting on/off. Input can be either `TRUE` or
   `FALSE`, default value `FALSE`.
- **site-weight**: (float) Sets the weights for the site energy
   calculation in the selection procedure. Default value 0.0.
- **energy-weight**: (float) Sets the weights for the energy
   calculation in the selection procedure. Default value 1.0.
- **force-weight**: (float) Sets the weights for the force calculation
   in the selection procedure. Default value 0.001.
- **stress-weight**: (float) Sets the weights for the force calculation
   in the selection procedure. Default value 0.0001.
- **threshold**: (float) Sets the maximum allowed extrapolation level
   for a configuration, i.e., anything larger than this will be
   included in the list of structures for the active learning to
   select from. Default of 2.0.
- **threshold-break**: (float) Sets the breaking threshold for the
    code. If a structure reaches this extrapolation grade then the
    code will break and not try to relax it any further. Default of
    10.0.

crystals_to_relax
*****************

The crystals to relax field takes a list of crystal structures to
include in the `to_relax.cfg` file. Valid options to include are `sc`,
`bcc`, `fcc`, `hcp`, and `prototypes` where `sc` stands for simple
cubic, `bcc` stands for body-centered cubic, `fcc` stands for
face-centered cubic, `hcp` stands for hexagonal close packed, and
`prototypes` is a list of structures that have been found to exist in
nature. The default behavior is to include all of these in the
`to_relax.cfg` file.

train
*****

The `train` section of the yml file can defines the options and values
that will get passed to the `mlp` `train` command when training the
interatomic potential. Valid options are:

- **energy-weight**: (float) The weight of energies in the
   fitting. Default is 1.0.
- **force-weight**: (float) The weight of forces in the
   fitting. Default is 0.01.
- **stress-weight**: (float) The weight of stresses in the
   fitting. Default is 0.001.
- **scale-by-force**: (float) If given value is >0 then configurations
   near equilibrium (with roughly force < the given value) get more
   weight. Default is 0.0.
- **max-iter**: (integer) The maximal number of BFGS (**!!explain!!**)
   iterations. Default is 1000.
- **trained-pot-name**: (string) The filename for trained potential to
   be saved to. Default is `Trained.mtp_`.
- **bfgs-conv-tol**: (float) The stopping criterion for optimization in
   the BFGS algorithm. Default is 1e-8.
- **weighting**: (string) Determines how to weight configuration wtih
   different sizes relative to each other. Options are `molecules`,
   `structures`, or `vibrations`, default is `vibrations`
- **init-params**: (string) Determins how to initialize parameters if a
   potential was not pre-fitted. Default is random. Options include
   `random` and `same` - this is when interaction of all species is
   the same (more accurate fit, but longer optimization), default is
   `random`.
- **skip-preinit**: If present then skip the 75 iterations done when
    params are not given.

calc_grade
**********

The `calc_grade` section of the yml file contains the options and
values that will get passed to the `mlp` `calc-grade` command for
calculating the grade of each configuration in the active set and
determining the best set of active configurations. Valid options are:

- **init-threshold**: (float) Sets the initial threshold to
   1+value. Default is 1e-5.
- **select-threshold**: (float) Sets the selection threshold to the value,
   i.e., nothing will be selected unless it has at least this
   extrapolation grade. Default is 1.1.
- **swap-threshold**: (float) Sets the swap threshold to the
   value. Default is 1.0000001.
- **energy-weight**: (float) Sets the weight for energy
   equation. Default is 1.
- **force-weight**: (float) Sets the weight for force
   equations. Default is 0.
- **stress-weight**: (float) Sets the weight for stress
   equations. Default is 0.
- **nbh-weight**: (float) Sets the weight for site energy
   equations. Default is 0.

relax
*****

The `relax` section of the yml file contains the options for the
relaxation step in the fitting iteration. The inputs will be passed to
the `mlp` `relax` command. Valid options are:

- **pressure**: (float) Sets an external pressure, in GPa. Default 0.
- **iteration_limit**: (integer) Sets the maximum number of BFGS
   iterations to take.
- **min_dist**: (float) Sets the minimum allowed distance between
   atoms. If atoms get any closer than this the relaxation will be
   terminated. Expressed in Angstroms.
- **force-tolerance**: (float) Sets the force tolerance for relaxation
   termination. If the forces on the atoms drop below this number then
   the relaxation will end. Zero disables atomic
   relaxations. Expressed in eV/Angstrom.
- **stress-tolerance**: (float) Sets the stress tolerance for
   relaxation tcermination. If the stress on the system drops below
   this number then relaxation will stop. Zero disables atomic
   relaxations. Expressed in GPa.
- **max-step**: (float) Sets the maximum allowed displacement of atoms
   and lattice vectors. Expressed in Angrstroms.
- **min-step**: (float) Sets the minimum atomic and cell vector
   displacement for relaxation. If the displacementens are smaller
   than this the relaxation stops, i.e., the cell is already
   relaxed. Expressed in Angrstroms.
- **bfgs-wolfe_c1**: If specified the bfgs wolfe code 1 will be used
   for relaxation.
- **bfgs-wolfe_c2**: If specified the bfgs wolfe code 2 will be used
   for relaxations.

select
******

The `select` section of the yml file contains options that will be
passed to the `mlp` `select` command. This command is used for the
active learning step and will decide which new structures from those
which failed to relax will be should be included in the next training
iteration. Valid options are:

- **init-threshold**: (float) Sets the initial threshold. Default 1E-5.
- **select-threshold**: (float) Sets the selection threshold. Default 1.1.
- **swap-threshold**: (float) Sets the threshold for configuration
   swapping. Default 1.0000001.
- **energy-weight**: (float) Sets the weight for the energy
   equation. Default 1.
- **force-weight**: (float) Sets the weight for the force
   equations. Default 0.
- **stress-weight**: (float) Sets the weight for the stress
   equations. Default 0.
- **nbh-weight**: (float) Sets the weight for the site energy
   equation. Default 0.
- **selection-limit**: (integer) Sets the largest number of structures
   that can be added to the training set in this iteration. Default 0
   (disables this function).
- **weighting**: (string) Sets the way of weighting the functional for
   better fitting. Options are `vibrations`, `molecules`,
   `structures`. Default is `vibrations`.

Other
*****

- **smallest_relax_cell**: (integer) The smallest cell size to include in the
   `to_relax.cfg` file. Values must be positive integers smaller than
   `largest_relax_cell`, default value 1.
- **largest_relax_cell**: (integer) The largest cell size to include in the
   `to_relax.cfg` file. Values must be positive integers larger that
   `smallest_relax_cell`, default value depends on the number of
   species in the system.
- **ran_seed**: (float) random seed to be used for MTP related calculations
   when a random seed is needed.
- **species**: (list of strings) This field takes a list of the atomic species in the
   system.
- **use_unrelaxed**: (bool) If true then if the `unrelaxed.cfg` file exists
   (during the relaxation step of the worflow, see below,
   configurations that can't be relaxed are saved to the
   `unrelaxed.cfg` file) then use it to start the relaxation step from
   rather than `to_relax.cfg`. Can take either `True` or `False`,
   default `False`.
- **use_mpi**: (bool) If TRUE then the `mtp` commands written to the
   jobfiles that can be executed using `mpi` will be. Can take `True`
   or `False`, default value `True`.


Example YML section
*******************

Example of the `fits` portion of the `fitting` section of the yml file
with `mtp`. For full documentation of the `fitting` section see
:doc:`../fitting`.

.. code-block:: yaml

  fits:
    - name: "AgPd"
      dbs: ["*"]
      steps:
        - type: "mtp.MTP"
          select: 
            selection-limit: 200
	  relax:
	    iteration_limit: 500
	  calc_grade:
	    selec-threshopld: 5.0
	  train:
	    max-iter: 500
	  crystals_to_relax:
	    - bcc
	    - fcc
	    - prototypes
	  larges_relax_cell: 8
	  use_mpi: False
          species:
            - "Ag"
            - "Pd"
   
MTP Potential Fitting Workflow
------------------------------

MTP uses an active learning framework to learn the potential landscape
for a material. Since this is obviously an interative process it
requires a slightly different workflow than other fitting
methods. Depending on which step in the process is being used the
number of cores and amound of memory to use specified in the `yml`
file's `execution` section should be adjusted.

The steps, once the initial database is complete, are as follows:

First Iteration
***************

1. Run `matdb_train.py cntrl -t` sets up the files needed for an
   `mtp` model to be trained.
2. Run `matdb_train.py cntrl -x` submits the computations to the queue
   (can be run on multiple processors).
3. Run `matdb_trani.py cntrl -t` sets up the files to build a list of
   structures to test the potential against.
4. Run `matdb_train.py cntrl -x` submits the computation to the queue (single core job).
5. Run `matdb_train.py cntrl -t` sets up files to perform a test
   relaxation step on the structure list created in steps 3 and 4.
6. Run `matdb_trani.py cntrl -x` submits the relaxation calculations
   to the queue (can be run with multiple processors).
7. Run `matdb_train.py cntrl -t` sets up the files needed to select
   new structures for the potential, i.e., the active learning
   algorithm selects structures that will improve the potential.
8. Run `matdb_train.py cntrl -x` to submit the job to the queue
   (single core process).
9. Run `matdb_train.py cntrl -t` adss the selected configurations to
   the :doc:`../database/active` for calculation with a first principles
   code and submits them to the job queue.

Second Iteration
****************

The second iteration can be started once the first principles
calculations are complete.

1. Run `matdb_build.py cntrl -e` to extract data from the new
   configurations in the :doc:`../database/active`.
2. Run `matdb_train.py cntrl -t` sets up the files to update the `mtp`
   potential.
3. Run `matdb_train.py cntrl -x` submits the computations to the queue
   (can be run on multiple processors).
4. Run `matdb_train.py cntrl -t` sets up files to perform a test
   relaxation step on the structure list created in steps 3 and 4 of
   the first iteration.
5. Run `matdb_trani.py cntrl -x` submits the relaxation calculations
   to the queue (can be run with multiple processors).
6. Run `matdb_train.py cntrl -t` sets up the files needed to select
   new structures for the potential, i.e., the active learning
   algorithm selects structures that will improve the potential.
7. Run `matdb_train.py cntrl -x` to submit the job to the queue
   (single core process).
8. Run `matdb_train.py cntrl -t` adss the selected configurations to
   the :doc:`../database/active` for calculation with a first principles
   code and submits them to the job queue.

The process for the second iteration is repeated until either the
potential is able to relax all structures in the list of test
configurations.
   
.. note:: It is possible to change the list of configurations that the
          potential is trained against at any point, however doing so
          will require that the list be re-computed and add those
          steps back into the workflow for that iteration.

MTP Modules
-----------

.. automodule:: matdb.fitting.mtp
   :synopsis: Trainers for creating MTP potentials automatically.
   :members:      

