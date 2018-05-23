"""Implements the classes needed for the Tracy calculations, i.e.,
those calculations that will be added to the Tracy compute queue.
"""

from matdb.database.utility import make_primitive
from matdb.descriptors import soap
from matdb.calculators import Qe

class Tracy(object):
    """Represents a calculator that will be submitted to the Tracy queue.

    Args:
        atoms (matdb.Atoms): the atoms object for this calculator instance.
        role (str): the role of the user, i.e., 'Scientist'.
    Attributes:
        role (str): The role of the user, i.e., "Cheif Scientist".
    """

    def __init__(self, atoms, role=None, notifications=None,
                 exec_dict = None, group_preds=None, contract_preds=None):

        self.atoms = atoms
        self.contract_type = None
        self.input_dict = None
        self.ecommerce_priority = 1
        

    def _compress_struct(self, atoms):
        """Compresess the input atoms object so that it is ready to be sent to
        the queue.
        
        Args:
            atoms (matdb.Atoms): the atoms object for the calculaton.
        
        Returns:
            A dictionary of the compressed structure that the
            decompression algorithm can unpack.
        """
        a_vecs, pos, types, hnf = make_primitive(self.atoms)
        type_map = {k:v for v, k in enumerate(np.unique(types))}

        result = {"a": [list(a) for a in a_vecs],
                  "b": [list(b) for b in pos],
                  "s": list(np.unique(types)),
                  "t": [type_map[i] for i in types]
                  "h": hnf}
        return result

    def _get_source():
        """Determines the priority of the source submitting the calculations.
        """
        source_dict = {"Chief Scientist": 10,
                       "Scientist": 9,
                       "Intern": 8,
                       "ASPEN": 7,
                       "Enumerated": 3}

        return source_dict[self.role]    

    def submit(self):

        """Submits the job to the authetication code before being sent to the queue.
        """
        package = {}
        package["MatDB ID"] = self.atoms.uuid
        package["MatDB Group ID"] = self.atoms.group_uuid
        # package["Source ID"] = Get this from Josh's scripit
        package["Source"] = def._get_source()
        package["Contract Type"] = self.contract_type
        package["Input Dictionary"] = self.input_dict
        package["Input Dictionary"]["cryst"] = self._compress_struct(self.atoms)
        package["Before fingerprint"] = soap(self.atoms)
        package["Contract Priority"] = 1#self._get_contract_priority()
        package["eCommerce Priority"] =
        # package["Date Ready"] = None
        job_reqs = self._get_job_reqs()
        for k,v in job_needs.values():
            package[k] = v

        # if :
        #     package["Group Predecessors"] =
        #     package["Contract Predecossors"] = 

    def get_authentication(self):


class Tracy_QE(Tracy, Qe):

        """Represents a DFT calculation that will be submitted to the Tracy queue.

    Args:

    Attributes:
        input_data (dict): A dictionary of the keywords and args needed to perform
          the calculations, using the ASE format.
    
    """

    def __init__(self, atoms, folder, contr_dir, ran_seed, *args, **kwargs):
        
        self.contract_type = 1
        self.input_dict = kwargs["input_data"]
        del kwargs["input_data"]
        if "tprnfor" in self.input_dict:
            self.input_dict["Control"]["tprnfor"] = self.input_dict["tprnfor"]
            del self.input_dict["tprnfor"]
        else:
            self.input_dict["Control"]["tprnfor"] = True

        if "tstress" in self.input_dict:
            self.input_dict["Control"]["tstress"] = self.input_dict["tstress"]
            del self.input_dict["tprnfor"]
        else:
            self.input_dict["Control"]["tstress"] = True

        super(Tracy_QE, Tracy).__init__(**kwargs)
        
            
        

    def _check_potcar(self):
        """We don't construct the potetial files on the user's end so we don't
        actuall want to do any checking.
        """
        pass

    def create(self, atoms):
        """We don't want to write any input files either.
        """
        pass

    def 
