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
import warnings
class abstractstatic(staticmethod):
    __slots__ = ()
    def __init__(self, function):
        super(abstractstatic, self).__init__(function)
        function.__isabstractmethod__ = True
    __isabstractmethod__ = True
    
def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func

def exhandler(function, parser):
    """If "examples" was specified in 'args', the specified function
    is called and the application exits.

    :arg function: the function that prints the examples.
    :arg parser: the initialized instance of the parser that has the
      additional, script-specific parameters.
    """
    args = vars(bparser.parse_known_args()[0])
    if args["examples"]:
        function()
        return
    if args["verbose"]:
        from matdb.msg import set_verbosity
        set_verbosity(args["verbose"])

    args.update(vars(parser.parse_known_args()[0]))
    return args

def _common_parser():
    """Returns a parser with common command-line options for all the scripts
    in the matdb suite.
    """
    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--examples", action="store_true",
                        help="See detailed help and examples for this script.")
    parser.add_argument("--verbose", default=0, type=int,
                        help="See verbose output as the script runs.")
    parser.add_argument('--action', nargs=1, choices=['save','print'], default='print',
                        help="Specify what to do with the output (print or save)")
    parser.add_argument("--debug", action="store_true",
                        help="Print verbose calculation information for debugging.")

    return parser

bparser = _common_parser()
testmode = False
"""bool: when True, the package is operating in unit test mode, which changes
how plotting is handled.
"""
def set_testmode(testing):
    """Sets the package testing mode.
    """
    global testmode
    testmode = testing

debug = False
"""bool: when True, loggers throughout the `matdb` system are set to debug level
to produce more verbose logging output. Can also be `int` corresponding to level
in :mod:`logging`.
"""
def set_debug(debugging):
    """Sets the package debug mode.
    """
    global debug
    debug = debugging
