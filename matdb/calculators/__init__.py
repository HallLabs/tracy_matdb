#Copyright (C) 2019  HALL LABS
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#If you have any questions contact: wmorgan@tracy.com

# the local name for the modules are used as the "name" instance attribute of the correspending class 
from .vasp import AsyncVasp as Vasp
from .aflux import AsyncAflow as Aflow
from .qe import AsyncQe as Qe 

from matdb import msg

def get_calc_class(name):
    """Gets the class definition objects for the calculator that has the
    specified name.
    """
    globs = globals()
    try:
        for k, v in globs.items():
            if k.lower() == name.lower() and k != name.lower():
                target = v
                break
    except KeyError: # pragma: no cover
        msg.err("Cannot import calculator {}. ".format(name) +
                "Does not exist at package level.")

    return target

def get_calculator_module(calcargs):
    """Returns the module corresponding to the calculator mentioned in
    `calcargs`.

    Args:
        calcargs (dict): the "calculator" dictionary that is part of the
          arguments for a db group.
    """
    from matdb import calculators
    from inspect import getmodule
    cls = getattr(calculators, calcargs["name"])
    mod = None

    try:
        mod = getmodule(cls)
    except: #pragma: no cover
        pass

    return mod
