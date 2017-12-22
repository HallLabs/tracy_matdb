"""Implements an abstract base class that all `matdb` calculators need to
implement. 
"""
import abc

class AsyncCalculator(object):
    """Represents a calculator such as :class:`ase.Calculator` that can be run
    in multiple stages to compute properties (such as in a distributed or HCP
    environment).
    """
    @abc.abstractmethod
    def can_execute(self, folder):
        """Returns `True` if the specified folder is ready to execute using the
        calculator's main executable.

        Args:
            folder (str): path to the folder in which the executable will be
              run.
        """
        pass

    @abc.abstractmethod
    def can_cleanup(self, folder):
        """Returns True if the specified folder has completed executing and the results
        are available for use.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass

    @abc.abstractmethod
    def is_executing(self, folder):
        """Returns True if the specified folder is in process of executing. This
        means that files/output has been produced to indicate that the process
        started and that the folder is not ready to cleanup yet.

        Args:
            folder (str): path to the folder in which the executable was run.
        """
        pass

    @abc.abstractmethod
    def create(self, folder, rewrite=False):
        """Creates all necessary input files for the calculator's executable.

        Args:
            folder (str): path to the folder in which to create input files.
        """
        pass

    @abc.abstractmethod
    def convert(self, folder, filename="output.xyz", fmt="xyz",
                properties=["species", "pos", "z", "ref_force"],
                parameters=["ref_energy", "ref_virial"],
                recalc=False):
        """Converts the output from the calculation to a new format.

        Args:
            folder (str): path to the folder in which the executable was run.
            filename (str): name of the file/object to create; this is created in
              each sub-sampled configurations directory/group.
            properties (list): of `str` *atom* property names (such as position,
              force, Z, etc.) to include in the output.
            parameters (list): of `str` *configuration* paramater names (such as
              energy, stress, etc.).
            recalc (bool): when True, re-create the files/objects, even if they
              already exist. 
        """
        pass

class SyncCalculator(object):
    """Represents a calculator such as :class:`ase.Calculator` that can be run
    synchronously stages to compute properties (does not require a distributed
    or HCP environment).
    """
    @abc.abstractmethod
    def can_execute(self, atoms):
        """Returns `True` if this calculation can calculate properties for the
        specified atoms object.

        Args:
            atoms (quippy.Atoms): config to test executability for.
        """
        pass

    @abc.abstractmethod
    def can_cleanup(self, atoms):
        """Returns True if the specified atoms object has completed executing and the
        results are available for use.

        Args:
            atoms (quippy.Atoms): config to check execution completion for.
        """
        pass

    @abc.abstractmethod
    def is_executing(self, atoms):
        """Returns True if the specified config is in process of executing.

        Args:
            atoms (quippy.Atoms): config to check execution for.
        """
        pass

    @abc.abstractmethod
    def execute(self, atoms):
        """Executes the calculator for the specified configuration.

        Args:
            atoms (quippy.Atoms): config to execute calculations for.
        """
        pass

    @abc.abstractmethod
    def create(self, atoms):
        """Initializes the calculator for the specified atoms object if
        necessary.

        Args:
            atoms (quippy.Atoms): config to initialize for.
        """
        pass

    @abc.abstractmethod
    def xyz(self, atoms, filename="output.xyz",
            properties=["species", "pos", "z", "ref_force"],
            parameters=["ref_energy", "ref_virial"],
            recalc=False):
        """Generates an extended XYZ file for the specified configuration.

        Args:
            atoms (quippy.Atoms): configuration to generate XYZ file for.
            filename (str): name of the XYZ file to create; this is created in
              each sub-sampled configurations directory.
            properties (list): of `str` *atom* property names (such as position,
              force, Z, etc.) to include in the XYZ file.
            parameters (list): of `str` *configuration* paramater names (such as
              energy, stress, etc.).
            recalc (bool): when True, re-create the XYZ files, even if they already
              exist. 
        """
        pass
