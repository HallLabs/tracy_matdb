"""Functions for interacting with various output formats.
"""
def vasp_to_xyz(folder, outfile="output.xyz", recalc=False,
                properties=["species", "pos", "z", "dft_force"],
                parameters=["dft_energy", "dft_virial"]):
    """Creates an extended XYZ file for the calculated structure in
    OUTCAR for the given folder.

    Args:
        folder (str): path to the folder to convert.
        outfile (str): name of the XYZ file to create. The file will
          be created in the same folder as the original if no absolute
          path is given.
        recalc (bool): when True, re-convert the OUTCAR file, even if
          the target XYZ file already exists.
    """
    from os import path, stat
    if not path.isabs(outfile):
        #Convert to absolute path if one wasn't given.
        outfile = path.join(folder, outfile)
        
    if (path.isfile(outfile)
        and stat(outfile).st_size > 100
        and not recalc):
        return True
        
    p = ','.join(properties)
    P = ','.join(parameters)
    renames = [("energy", "dft_energy"), ("force", "dft_force"),
               ("virial", "dft_virial")]
    sargs = ["convert.py", "-I", "OUTCAR", "-p", p, "-P", P, "-f", "xyz"]
    for s, d in renames:
        sargs.append("-n")
        sargs.append(s)
        sargs.append(d)
        
    sargs.extend(["-o", outfile, "OUTCAR"])

    from matdb.utility import execute
    execute(sargs, folder)
    return path.isfile(outfile) and stat(outfile).st_size > 100
