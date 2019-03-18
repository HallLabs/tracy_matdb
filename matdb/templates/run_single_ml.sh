#!/bin/bash

#SBATCH --time={{ time }}:00:00   # walltime
#SBATCH --ntasks={{ ntasks }}   # number of processor cores (i.e. tasks)
#SBATCH --nodes={{ nodes }}   # number of nodes
#SBATCH --mem-per-cpu={{ mem_per_cpu  }}   # memory per CPU core
#SBATCH -J "{{ job_name }}"   # job name
{%- if options %}
{%- for option in options %}
#SBATCH {{option}}
{%- endfor %}
{%- endif %}
{%- if partition %}
#SBATCH --partition={{ partition }}
{%- endif %}

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
cd {{ execution_path }}

{%- if modules_unload %}
{%- for modname in modules_unload %}
module unload {{ modname }}
{%- endfor %}
{%- endif %}
{%- if modules_load %}
{%- for modname in modules_load %}
module load {{ modname }}
{%- endfor %}
{%- endif %}

{%- if preamble %}
{{preamble}}
{%- endif %}
# Get the path to the executable; should be on user's path after the modules have been loaded.
{{ exec_path }}

