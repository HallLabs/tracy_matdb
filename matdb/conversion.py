"""Implements conversion between file formats for :class:`matdb.AtomsList`
objects.
"""
from tqdm import tqdm
from os import path
from matdb import msg
def to_xyz(atomslist, outfile, overwrite=False):
    """Converts the specified list of atomic configurations to the extended XYZ
    format. 

    .. note:: This function requires the use of :mod:`quippy`.

    Args:
        atomslist (matdb.AtomsList): list of configurations to write to XYZ.
        outfile (str): path to the output file to write with the
          configurations.
        overwrite (bool): when True, overwrite the output file if it already
          exists; otherwise, issue a warning and do nothing.
    """
    import quippy
    if path.isfile(outfile):
        msg.warn("The output file {} already".format(outfile) +
                 " exists; aborting XYZ conversion.")
        return
    
    ol = quippy.AtomsList()
    for at in tqdm(atomslist):
        #Empty dictionaries in info breaks the quippy fortran
        #implementation. Delete these entries out.
        if len(at.info["params"]) == 0:
            del at.info["params"]
        if len(at.info["properties"]) == 0:
            del at.info["properties"]
        
        ai = quippy.Atoms()
        ai.copy_from(at)
        if "properties" in at.info:
            #We also need to copy the properties in our info onto the
            #properties of the quippy.Atoms object.
            ai.arrays.update(at.info["properties"])
        if "params" in ai.params:
            del ai.params["params"]
            
        ol.append(ai)
    ol.write(outfile)
