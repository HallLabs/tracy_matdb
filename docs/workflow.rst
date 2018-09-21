Database Construction Workflow
==============================

After you have familiarized yourself with :doc:`matdbyml`, you are
ready to initiate a construction workflow. The workflow uses the
:doc:`scripts/build` to grant high-level access to sequences of API
calls.

The workflow for `matdb` typically follows these steps:

1. Run `matdb_build.py matdb.yaml -s` to run the setup of all VASP
   folders that need to be computed.
2. Run `matdb_build.py matdb.yaml --status` to see whether the folders
   are setup correctly.
3. Run `matdb_build.py matdb.yaml -x` to queue the jobs with the job
   scheduler.
4. Run `matdb_build.py matdb.yaml -e` to extract the relevant
   parameters for the next round of calculations.

These steps are repeated until all databases have been calculated and
cleaned up. 
