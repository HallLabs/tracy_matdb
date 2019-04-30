SBATCH
======

`sbatch` is used to submit jobs to a `slurm` style (**!!maybe a hyperlink to `slurm` style??**) workload manager on
many supercomputers. There may be differences in implementation of
slurm templates for different supercomputers but the following fields
should be applicaple to most use cases.

yml File specifications
-----------------------

There are a number of fields that must be provided for population to
the `sbatch` template. These fields are:

- **template**: (string) The name of the template file located in
  `matdb/templates` to populate.
- **time**: (integer) Time, in hours, that job will be allowed to run for.
- **ntasks**: (integer) The number of tasks that are being performed, i.e., how
  many cpus to request.
- **nodes**: (integer) The number of nodes to request.
- **mem_per_cpu**: (integer) The amounnt of memory, in giga-bytes, to
  request on each cpu used.
- **exec_path**: (string) Path to the executable that will be called
  by `slurm`.
- **job_name**: (string) The job name for human reference.

Additional opptional parameters include:
  
- **partition**: (string) The name of a special partition to assign
  the jobs to.
- **array_limit**: (integer) Only to be used if the template supports
  job arrays. Sets the maximum number of jobs that can run at a time.
- **modules_load**: (list) Loads the specified modules into the session.
- **modules_unload**: (list) Unloads the specified modules into the session.
- **preamble**: (string) Anything that needs to be done immediately
  before execution of the script.
- **options**: (list) Any other `slurm` options that the user wishes
  to include.

For the preamble, you can use multi-line strings in YAML like:

.. code-block:: yaml

   preamble: >
     first line
     second line
     and so on...


Templates included in `matdb` for `sbatch` are:

1) `run_array_ml.sh`: used to run job arrays, such as the training of
   databases where multiple calculations can be happening
   simultaneously.
2) `run_single_ml.sh`: used to run a single job that won't make uso fo
   job arrays.


Example
-------

.. code-block:: yaml

   execution:
     template: 'run_array_ml.sh'
     time: 4
     ntasks: 1
     nodes: 1
     mem_per_cpu: 4
     job_name: 'AgPd DB'
     partition: 'physics'
     array_limit: 50
     modules_load: ['mkl/11.2.0']
     exec_path: 'vasp'

