"""Class for generating and interacting with a database of
configurations generated from liquid-temperature molecular dynamics.
"""
from .basic import Database
class LiquidDatabase(Database):
    """Represents a sub-sampled molecular dynamics run created at a
    specific temperature.
    """
    def __init__(self, atoms):
        pass
