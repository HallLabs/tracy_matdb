"""This module contains the functions needed to submit jobs in different environments.
"""

from matdb.utility import execute

def sbatch(jobfile, folder, kwargs):
    """Submits jobs to a supercomputer with sbatch.

    Args:
       jobfile (str): the name of jobfile that needs to be submitted.
       folder (str): the path to the jobfile.
       kwargs (dict): other arguments to pass to the execute function.
    """

    execute(["sbatch",jobfile],**kwargs)

def tracy():
    """Submits jobs for the tracey workflow.
    """

    raise NotImplementedError()
