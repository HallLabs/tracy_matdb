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
4. Run `matdb_build.py matdb.yaml -c` to cleanup the VASP calculations
   and extract the relevant parameters for the next round of
   calculations.

These steps are repeated until all databases have been calculated and
cleaned up. For example, after running this the first time, the
`PhononDFT` calculations will be finished and cleaned up. These are
needed by `PhononDatabase` to created modulated configurations, so we
run all four steps again. This can be repeated until the database is
ready.
