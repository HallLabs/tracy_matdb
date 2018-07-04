"""Implements the classes needed for the Tracy calculations, i.e.,
those calculations that will be added to the Tracy compute queue.
"""
from datetime import datetime
from os import path
from random import seed, uniform
import json
import abc
import numpy as np
import binascii

from matdb.database.utility import make_primitive
from matdb.descriptors import soap
from matdb.calculators.basic import AsyncCalculator
from matdb.calculators import Qe
from matdb.utility impor _get_reporoot

def _get_struct_info(struct, local_label):
    """Gets the structure info for the atomic species.

    Args:
        struct (str): The structure for which data is wanted.
        local_label (int): The local, i.e., this instance of matdb 
          label for the structure.
    
    Returns:
        Dictionary of structure data.
    """

    struct_data = ""
    with open(path.join(_get_reporoot(), "matdb", "templates",
                        "struct_enum.out.bin"), "rb") as f:
        for line in f:
            struct_data += _text_from_bits(line)

    data_dict = json.loads(struct_data)

    for s in temp["structs"]["details"]:
        if s["title"] == struct:
            detalis = s.copy()
            details["title"] = local_label
            details["fileName"] = details["fileName"].replace(struct,str(local_label))

    return details

def _text_from_bits(bits, encoding='utf-8', errors='surrogatepass'):
    """Converts the data in bits to texte.

    Args:
        bits (binary): binary data to convert.

    Returns:
        The text inside the binary.
    """
    n = int(bits, 2)
    return int2bytes(n).decode(encoding, errors)

def _int2bytes(i):
    """Converts the integer to bytes.
    
    Args:
        i (int): integer to convert.
    """
    hex_string = '%x' % i
    n = len(hex_string)
    return binascii.unhexlify(hex_string.zfill(n + (n & 1)))        

class Tracy(AsyncCalculator):
    """Represents a calculator that will be submitted to the Tracy queue.

    Args:
    Attributes:
        role (str): The role of the user, i.e., "Cheif Scientist".
    """
    key = "tracy"
    
    def __init__(self, folder, role=None, notifications=None,
                 group_preds=None, contract_preds=None, ecommerce=None,
                 contract_priority=None, max_time=None, min_flops=None,
                 min_mem=None, ncores=None, max_net_lat=None, min_ram=None,
                 user=None, password=None):

        # If this folder has already been submitted then there will be
        # a file containing the contract number in it. We want to read
        # this in so we can check it's status later.
        if user is None or password is None
            raise ValueError("A user name and password must be specified to use a "
                             "Tracy calculotar.")

        self.user = user
        self.password = password
        self.folder = folder
        if path.isfile(path.join(folder, "contract.txt")):
            with open(path.join(folder, "contract.txt"), "r") as f:        
                self.contract_id = f.readline().strip()
        else:
            self.contract_id = None
            
        if path.isfile(path.join(folder, "post_print.txt")):
            with open(path.join(folder, "post_print.txt"), "r") as f:
                self.after_print = f.readline().strip()
        else:
            self.after_print = None

        # All these will be set by either input or latter methods
        self.role = role
        self.type_map = None
        self.contract_type = None
        self.input_dict = None
        self.ecommerce = ecommerce if ecommerce is not None else 1
        self.group_preds = group_preds
        self.contract_preds = contract_preds
        self.contract_priority = contract_priority if contract_priority is not None else 1
        self.sys_specs = {"max_time": max_time, "min_flops": min_flops, "min_mem": min_mem,
                          "ncores": ncores, "min_ram": min_ram, "max_net_lat": max_net_lat}
        self.notifications = notifications

    def _compress_struct(self, atoms):
        """Compresess the input atoms object so that it is ready to be sent to
        the queue.
        
        Args:
            atoms (matdb.Atoms): the atoms object for the calculaton.
        
        Returns:
            A dictionary of the compressed structure that the
            decompression algorithm can unpack.
        """
        if np.allclose(0, atoms.cell):
            raise ValueError("The Atoms object must contian cell vectors in order "
                             "to be compressed.")
        
        a_vecs, pos, types, hnf = make_primitive(self.atoms)
        if self.type_map is None:
            self.type_map = {}
            for v, k in enumerate(types):
                self.type_map[k] = v+1
            

        result = {"a": [list(a) for a in a_vecs],
                  "b": [list(b) for b in pos],
                  "t": self._intarray_to_int([self.type_map[i] for i in types]),
                  "h": self._intarray_to_int([hnf[0][0], hnf[1][0],
                                              hnf[1][1], hnf[2][0],
                                              hnf[2][1], hnf[2][2]], pad=True)}
        return result

    @staticmethod
    def _intarray_to_int(int_array, pad=False):
        """Converts an integer into a single integer.

        Args:
            int_array (list): the integer array to convert.
            pad (bool): if True pad the values with zeros.
        """

        if pad:
            frm_str = "{0}0"
            shift = 1
        else:
            frm_str = "{0}"
            shift = 0

        full_int = "".join([frm_str.format(i+shift) for i in int_array])

        return int(full_int)

    @abc.abstractmethod
    def get_input_dict(self): #pragma: no cover
        """Constructs the input dictionary from the files written.
        """
        pass
    
    def _get_source(self):
        """Determines the priority of the source submitting the calculations.
        """
        source_dict = {"Chief Scientist": 10,
                       "Scientist": 9,
                       "Intern": 8,
                       "ASPEN": 7,
                       "Enumerated": 3}

        if self.role in source_dict.keys():
            return source_dict[self.role]
        else:
            raise ValueError("The source {0} is not recognized. Cannot assign "
                             "priority.".format(self.role))
    

    def prep_submit(self, folder):
        """Submits the job to the authetication code before being sent to the queue.
        """
        self.get_input_dict()
        package = {}
        package["matDbId"] = self.atoms.uuid
        package["matDbGroupId"] = self.atoms.group_uuid
        # package["Source ID"] = Get this from Josh's scripit
        package["source"] = self._get_source()
        package["contractType"] = self.contract_type
        package["input"] = self.input_dict
        package["input"]["cryst"] = self._compress_struct(self.atoms)
        package["beforeFingerprint"] = soap(self.atoms)
        package["contractPriority"] = self.contract_priority
        package["eCommercePriority"] = self.ecommerce
        package["maximumProcessingTime"] = self.sys_specs["max_time"]
        package["minimumFlops"] = self.sys_specs["min_flops"]
        package["minimumRam"] = self.sys_specs["min_ram"]
        package["minimumStorage"] = self.sys_specs["min_mem"]
        package["numberOfCores"] = self.sys_specs["ncores"]
        package["maximumNetworkLatency"] = self.sys_specs["max_net_lat"]
        package["user"] = self.user
        package["password"] = self.password

        #Move this to the wrapper script that reads the json.
        # if self.can_execute(self.folder):
        #     fmt = '%Y%m%d%H%M%S'
        #     package["Date Ready"] = datetime.now().strftime(fmt)
            
        if self.group_preds is not None:
            package["groupPredecessors"] = self.group_preds
        if self.contract_preds is not None:
            package["contractPredecossors"] = self.contract_preds

        if self.notifications is not None:
            package["notifications"] = self.notifications
            
        target = path.join(folder, "submission.json")

        with open(target, "w+") as f:
            json.dump(package, f)
        

    def can_extract(self, folder):
        """Returns `True` if the calculation has been completed and the data
        is ready for extraction.
    
        Args:
            folder (str): the directory to the folder.
        """
        # Here we need to do a query of the endpoint to determine if
        # the calculation has been completed.
        pass

    def is_executing(self, folder):
        """Returns `True` if the calculation is being performed on the queue.
        """
        # Here we need to da a querry of the endpoint to determin if
        # the calculation has been started.
        pass
        
class Tracy_QE(Tracy, Qe):
    """Represents a DFT calculation that will be submitted to the Tracy queue.

    Args:
        atoms (matdb.Atoms): configuration to calculate using QE.
        folder (str): path to the directory where the calculation should take
          place.
        contr_dir (str): The absolute path of the controller's root directory.
        ran_seed (int or float): the random seed to be used for this calculator.
        kwargs (dict): a dictionary of the kwargs to be passed to the queue.

    Attributes:
        input_data (dict): A dictionary of the keywords and args needed to perform
          the calculations, using the ASE format.
    
    """

    key = "tracy_qe"
        
    def __init__(self, atoms, folder, contr_dir, ran_seed, *args, **kwargs):

        self.in_kwargs = kwargs.copy()
        self.contract_type = 1
        self.QE_input = kwargs["calcargs"]
        self.tracy_input = kwargs["tracy"]
        self.ran_seed = ran_seed
        self.contr_dir = contr_dir
        self.folder = folder
        self.atoms = atoms

        if self.ran_seed is not None:
            seed(self.ran_seed)            

        Qe.__init__(self, atoms, folder, contr_dir, ran_seed, **self.QE_input)
        Tracy.__init__(self, folder, **self.tracy_input)

    def _check_potcars(self):
        """We don't construct the potetial files on the user's end so we don't
        actuall want to do any checking.
        """
        pass

    def write_input(self, atoms):
        """Writes the input to a folder.
        
        Args:
            atoms (Atoms): an instance of :class:`matdb.atoms.Atoms`.
        """

        Qe.write_input(self, atoms)
        self.prep_submit(self.folder)

    def get_input_dict(self):
        """Reads in the input file to a dictionary.
        """
        self.input_dict = {}
        self.type_map = {}
        skip_until_next_key = False
        k_val = None
        with open(path.join(self.folder, "espresso.pwi"), "r") as f:
            for line in f:
                if not line.strip() == '':
                    if "&" in line:
                        key = line.strip()[1:]
                        skip_until_next_key == False
                    elif line.strip().lower() == "atomic_species":
                        key = "atomic_species"
                        skip_until_next_key == False
                    elif line.strip().split()[0].lower() == "k_points":
                        key = "k_points"
                        k_val = line.strip().split()[1]
                        skip_until_next_key == False
                    elif line.strip() != '' and line.strip().split()[0].lower() in ["cell_paramaters", "atomic_positions"]:
                        skip_until_next_key = True
                    elif not skip_until_next_key:
                        if "=" in line:
                            sub_key, val = line.strip().split("=")
                            sub_key = sub_key.strip()
                            val = val.strip()
                            if sub_key == "pseudo_dir":
                                continue
                            self.input_dict[key][sub_key] = val
                        elif line.strip() == "/":
                            continue
                        elif k_val is not None and key == "k_points":
                            self.input_dict[key][k_val] = line.strip()
                        elif key == "atomic_species":
                            species = line.strip().split()[0]
                            if species not in self.type_map.keys():
                                self.type_map[species] = len(self.type_map.keys()) + 1
                            
                            self.input_dict[key][self.type_map[species]] = uniform(0, 100)
                        else: #pragma: no cover
                            msg.warn("Could no process line {0} of file "
                                     "{1}".format(line, path.join(self.folder, "espresso.pwi")))
                            
                    if key not in self.input_dict.keys():
                        self.input_dict[key] = {}

       self.input_dict["structure_data"] = self._get_data()

    def _get_data(self):
        """Uses the QE input to construct the dictionary of information needed
        for the cloud compute code.
        """
        results = {"numEntries": len(self.type_map.keys()),
                   "details":[]}

        for struct, val in self.type_map.items():
            results["details"].append(_get_struct_info(struct, val))

        return results

    def can_execute(self, folder):
        """Returns True if the specified file is ready to be submitted to the queue.
        
        Args:
            folder (str): the path to the folder.
        """
        qe_ready = Qe.can_execute(self, folder)
        sub_ready = path.isfile(path.join(folder,"submission.json"))
        return qe_ready and sub_ready
                    
    def extract(self, folder):
        """Extracts the results from the returned data package.
        """
        # This will be handled externally for now.
        pass
        
    def to_dict(self):
        """Converts the arguments of the calculation to a dictionary.
        """
        results = {"kwargs": self.in_kwargs, "folder": self.folder,
                   "ran_seed": self.ran_seed, "contr_dir": self.contr_dir}
        return results        
