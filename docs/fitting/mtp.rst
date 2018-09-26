Fitting for MTP Potentials
==========================

MTP (Moment Tensor Potentials) use a polynomial basis to expand the
interatomic potential and offer fast calculation times for relaxing
cells and calculating the energies and forces of an atomic
configuration. The pontentials are also adaptable in terms of
accuracy, in that they can be expanded to be more accurate but then
they take a little longer to run.

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
   the :doc:`database/active` for calculation with a first principles
   code and submits them to the job queue.

Second Iteration
****************

The second iteration can be started once the first principles
calculations are complete.

1. Run `matdb_build.py cntrl -e` to extract data from the new
   configurations in the :doc:`database/active`.
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
   the :doc:`database/active` for calculation with a first principles
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
