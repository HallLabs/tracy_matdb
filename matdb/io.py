"""Functions for interacting with various output formats.
"""
def vasp_to_xyz(folder, outfile="output.xyz", recalc=0,
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
        and recalc <= 0):
        return True
        
    p = ','.join(properties)
    P = ','.join(parameters)
    renames = [("energy", "dft_energy"), ("force", "dft_force"),
               ("virial", "dft_virial")]
    sargs = ["convert.py", "-I", "vasprun.xml", "-p", p, "-P", P, "-f", "xyz"]
    for s, d in renames:
        sargs.append("-n")
        sargs.append(s)
        sargs.append(d)
        
    sargs.extend(["-o", outfile, "vasprun.xml"])

    from matdb.utility import execute
    execute(sargs, folder, errignore="OMP_STACKSIZE")
    return path.isfile(outfile) and stat(outfile).st_size > 100

import yaml
from os import path
from matdb.utility import chdir
def is_link(obj):
    """Determines whether the specified object is a link according to the `matdb`
    templating specification.
    """
    result = False
    if isinstance(obj, str):
        if len(obj) > 0:
            result = obj[0] == ":"
    return result

def _unpack_obj(context, obj, lcontext=None):
    """Unpacks each item of the specified object recursively so that all
    dictionary values are visited and all list items are also visited.

    .. warning:: `obj` will be mutated if any value it considers turns out to be
      a link (according to :func:`is_link`). In that case, the file descriptor
      will be placed by the actual contents of the YAML file that the link
      points to.

    Args:
        context (str): path to the root folder where the yaml file is
          located. Needed for relative paths of file links.
        lcontext (dict): local context for the items in `obj`. Keys are the
          names of keys in `obj`; values are relative folder paths that should
          be used as the context for reads within that item.
    """
    if isinstance(obj, dict):
        result = obj
        for k, o in obj.items():
            ncontext = context
            #If the template specifies a relative context for this item,
            #then switch out the context for all of its children.
            if lcontext is not None and k in lcontext:
                with chdir(context):
                    ncontext = path.abspath(lcontext[k])
            
            if is_link(o):
                result[k] = read(ncontext, o)
            else:
                result[k] = _unpack_obj(ncontext, o)
    elif isinstance(obj, (list, set, tuple)):
        result = []
        for o in obj:
            if is_link(o):
                result.append(read(context, o))
            else:
                result.append(_unpack_obj(context, o))
    else:
        result = obj
                
    return result
                
def read(context, yfile):
    """Reads in the specified YAML file, following any additional file
    directives to compile a full representation of the template hierarchy for
    the root file.

    Args:
        context (str): path to the root folder where the yaml file is
          located. Needed for relative paths of file links.
        yfile (str): name of the template YAML file *relative* to
        `context`. Should *not* include the `.yaml` or `.yml` extension.
    """
    with chdir(context):
        if yfile[0] == ":":
            root = path.abspath(yfile[1:])
        else:
            root = path.abspath(yfile)

    if path.isfile(root + ".yml"):
        target = root + ".yml"
    else:
        emsg = ("The specified template file '{}' was not found relative "
                "to the given context directory ('{}'). Note that all files"
                " should use the `.yml` extension, *not* `.yaml`.")
        raise ValueError(emsg.format(yfile, context))
        
    with open(target, 'r') as stream:
        result = yaml.load(stream)

    #Determine the new context for recursive file links within the values of
    #this file.
    ncontext = path.dirname(target)

    #The specification allows for a "local" context that describes folder
    #locations for specific items within the template.
    lcontext = None
    if isinstance(result, dict) and "context" in result:
        lcontext = result["context"]
        del result["context"]
    
    #The unpacking command will mutate the values in result so that file links
    #are expanded to be full-fledged python objects from their YAML files.
    _unpack_obj(ncontext, result, lcontext)
    return result
