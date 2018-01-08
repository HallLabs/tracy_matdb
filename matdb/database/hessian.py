"""Implements a hessian database by extracting eigenvectors and eigenvalues of
the dynamical matrix for each of the configurations in `seed`.
"""
from .basic import Group
class Hessian(Group):
    """Represents a collection of Hessian eigenvectors and eigenvalues
    associated with dynamical matrices calculated for local minima in the phase
    space.

    Args:
        dynmat (dict): dynamical matrix with calculated eigenvectors and
          eigenvals.
    """
    def __init__(self, name="hessian", **dbargs):
        self.name = name
        self.seeded = True
        dbargs["prefix"] = 'H'
        super(Hessian, self).__init__(**dbargs)

    @property
    def rset(self):
        """Returns a list of folder names, each of which has an `atoms.json`
        file with configuration and hessian eigenvector/eigenvalue data
        attached.
        """
        if len(self.sequence) == 0:
            #Just return this group's hessian configurations; it is at the
            #bottom of the stack.
            return self.configs.values()
        else:
            result = []
            for h in self.sequence.values():
                result.extend(h.rset)
            return result

    def ready(self):
        """Determines if all the AFLOW configurations specified in the query
        have been downloaded and saved yet.
        """
        if isinstance(self.seed, dict) and "atoms" in self.seed:
            return len(self.configs) > 0
        else:
            return all(h.ready() for h in self.sequence.values())

    def setup(self, rerun=False):
        """Diagonalizes the dynamical matrix for each seed configuration and creates a
        set of new configurations that includes each eigenvector and eigenvalue.
        """
        super(Hessian, self).setup(self._setup_configs, rerun)
        
    def _setup_configs(self, rerun):
        """Sets up the folders for each of the eigenvalue/eigenvector sub-configs
        calculated from the Hessian matrix. Extracts the hessian matrix and
        assigns its eigenvalues and eigenvectors to relevant properties on
        copies of the seed configuration atoms objects.

        Args:
            rerun (bool): when True, re-diagonalize all the seed
              configurations.
        """
        dmatrix = self.seed["dynmat"]
        eigvecs, eigvals = dmatrix["eigvecs"], dmatrix["eigvals"]
        
        hname = "hessian"
        for i, e in enumerate(eigvals):
            #Ignore zero eigenvalues, not useful.
            if not np.isclose(e, 0.0):
                atc = self.seed["atoms"].copy()
                #Add this eigenvector to its own configuration.
                atc.add_property(hname, 0.0, n_cols=3)
                H = np.reshape(eigvecs[:,i], (atc.n, 3)).T
                setattr(atc, hname, H)
                #Same thing for the eigenvalue. Then save it to the group folder
                #structure.                
                atc.params.set_value(hname, e)
                atc.params.set_value("n_hessian", 1)
                #The seed should already have force information attached from
                #the dynamical matrix step. atc.add_property("force", 0.0, n_cols=3)
                self.create(atc)
