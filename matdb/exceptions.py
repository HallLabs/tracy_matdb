"""Contains all the custom defined exceptions that matdb can throw.
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
    
