"""This module contains the functions needed to submit jobs in different environments.

Copyright (C) 2019  HALL LABS

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

If you have any questions contact: wmorgan@tracy.com
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
