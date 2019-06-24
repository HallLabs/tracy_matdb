"""Contains all the custom defined exceptions that matdb can throw.

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

class Error(Exception):
    """Base class for exceptions."""
    pass

class VersionError(Error):
    """Exception raised when version conditions aren't met.
    
    Attributes:
        message (str): Explanation of the error.
    """
    def __init__(self, message):
        self.message = message
    
class SpeciesError(Error):
    """Exception raised when atomic species don't match.
    
    Attributes:
        message (str): Explanation of the error.
    """
    def __init__(self, message):
        self.message = message
    
    
class LogicError(Error):
    """Exception raised when internal logic failed, most likely a bug!
    
    Attributes:
        message (str): Explanation of the error.
    """
    def __init__(self, message):
        self.message = message
    
class MlpError(Error):
    """Exception raised when the mlp code failed to produce the correct output.
    
    Attributes:
        message (str): Explanation of the error.
    """
    def __init__(self, message):
        self.message = message
