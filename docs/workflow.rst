Workflow
========

Database Construction Workflow
------------------------------

After you have familiarized yourself with :doc:`matdbyml`, you are
ready to initiate a construction workflow, throughout the rest of this
documentation it is assumed that a `yml` file name `cntrl.yml` is
being used though any name ending in `yml` will work. The workflow
uses the :doc:`scripts/build` to grant high-level access to sequences
of API calls.

The workflow for `matdb` database creation typically follows these steps:

1. Run `matdb_build.py cntrl -s` to run the setup of all VASP
   folders that need to be computed.
2. Run `matdb_build.py cntrl --status` to see whether the folders
   are setup correctly.
3. Run `matdb_build.py cntrl -x` to queue the jobs with the job
   scheduler.
4. Run `matdb_build.py cntrl -e` to extract the relevant
   parameters for the next round of calculations.

Since it is possible to have some databases depond on the output of
others for seed configurations it these steps are repeated until all
databases have been calculated and cleaned up.

Potential Construction Workflow
-------------------------------

Once all desired databases have been computed it is possible to train
an interatomic potential. The pontetial training as controlled by the
same `yml` file that the database creation was written from.

The workflow for `matdb` potential traininig looks like:

1. Run `matdb_train.py cntrl -t` to setup the files needed for the
   training.
2. Run `matdb_train.py cntrl -x` submit the jobs for training to
   the job scheduler.

For some of the interatomic potentials, such as :doc:`fitting/mtp`
which employs active learning, it is neccessary to permform the above
steps iteratively. Details of the exact requirements for each fitting
method are available in the documentation of the fitter. Fitting
methods currently supported are:

1. :doc:`fitting/mtp`
