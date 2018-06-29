"""Implements the classes needed for the Tracy calculations, i.e.,
those calculations that will be added to the Tracy compute queue.
"""

from matdb.calculators import Qe

class Tracy(object):
    """Represents a calculator that will be submitted to the Tracy queue.

    Args:
        calculation_type (str)
    """

    def __init__(self, contract_type=None, matdb_id=None, source_id=None, source=None,
                 matdb_group_id=None, input_dict=None, pre_print=None, post_print=None,
                 notifications=None, contract_priority=None, ecommerce_priority=None,
                 max_time=None, min_flops=None, min_ram=None, min_mem_storage=None,
                 n_cores=None, max_net_latency=None, group_preds=None, contract_preds=None):
        


    def submit(self):

    def get_authentication(self):


class Tracy_QE(Tracy, Qe):
    """Represents a DFT calculation that will be submitted to the Tracy queue.

    Args:

    Attributes:
    
    """

    def __init__(self, atoms, folder, contr_dir, ran_seed, *args, **kwargs):
        
