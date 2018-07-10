"""Functions for interacting with various output formats.
"""
from os import path

import ase
from ase.calculators.singlepoint import SinglePointCalculator
import h5py
import numpy as np
import re
import yaml

# from matdb.atoms import Atoms, AtomsList
from matdb.atoms import Atoms, AtomsList
from matdb.utility import chdir, execute

_rxcfg = re.compile(r"[a-z\s:\n]+", re.I)

_cfg_pos = ["cartes_x", "cartes_y", "cartes_z"]
"""list: of `str` label names in the CFG file corresponding to the x, y and z
positions respectively.
"""
_cfg_force = ["fx", "fy", "fz"]
"""list of `str` label names in the CFG file corresponding to the x, y and z
forces respectively.
"""

def order_stress(xx=None, yy=None, zz=None, yz=None, xz=None, xy=None):
    """Ensures that the six unique components of the stress tensor are presented
    in the correct order assumed by ASE.
    """
    return (xx, yy, zz, yz, xz, xy)

def symmetrize(xx=None, yy=None, zz=None, yz=None, xz=None, xy=None):
    """Returns the 3x3 stress matrix from the specified components.
    .. note:: the components are: xx yy zz xy yz zx.
    """
    from numpy import array
    return array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def atoms_to_cfg(atm, target, config_id=None, type_map=None):
    """Converts an :class:`matdb.atoms.Atoms` object to a cfg file.

    Args:
        atm (matdb.atoms.Atoms): the atoms object.
        target (str): path to the cfg file to be written.
        config_id (str): the config id for the atoms object.
        type_map (dict): a type map to match the species to a larger
          system than present in the system.
    """

    chem_syms = atm.get_chemical_symbols()
    pos = atm.positions

    if hasattr(atm, "calc") and atm.calc is not None and hasattr(atm.calc, "key"):
        calc_name = "{0}_".format(atm.calc.key)
    else:
        calc_name = ""

    local_map = {}
    for i, specs in enumerate(np.unique(chem_syms)):
        if type_map is None:
            local_map[specs] = i
        else:
            local_map[specs] =  type_map[i]

    with open(target, "w+") as f:
        f.write("BEGIN_CFG\n Size\n")
        f.write("    {}\n".format(len(pos)))
        f.write(" SuperCell\n")
        # Then write out the lattice vectors.
        lat_vecs = atm.cell
        for i in range(3):
            f.write("   {}\n".format("      ".join(
                ["{0: .6f}".format(j) for j in lat_vecs[i]])))

        f.write("  ")

        if "{0}force".format(calc_name) in atm.properties:
            force = atm.properties["{0}force".format(calc_name)]
            f.write(" AtomData:  id type       cartes_x      cartes_y      "
                    "cartes_z           fx          fy          fz\n")
        else:
            f.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z\n")
        iAt = 0
        for type, loc in zip(chem_syms, pos):
            out_lab = local_map[type]
            if "{0}force".format(calc_name) in atm.properties:
                f.write("             {0}    {1}       "
                        "{2}    {3}\n".format(iAt+1, out_lab,
                                              "  ".join(["{0: .8f}".format(i) for i in loc]),
                                              "    ".join(["{0: .8f}".format(i) for i in force[iAt]])))
            else:
                f.write("             {0}    {1}       "
                        "{2}\n".format(iAt+1, out_lab,
                                       "  ".join(["{0: .8f}".format(i) for i in loc])))
            iAt += 1

        if "{0}energy".format(calc_name) in atm.params:
            f.write("  Energy\n")
            f.write("        {0}\n".format(atm.params["{0}energy".format(calc_name)]))

        if "{0}virial".format(calc_name) in atm.params:
            virial = [atm.params["vasp_virial"][0][0], atm.params["vasp_virial"][1][1],
                      atm.params["vasp_virial"][2][2], atm.params["vasp_virial"][2][1],
                      atm.params["vasp_virial"][2][0], atm.params["vasp_virial"][0][1]]
            f.write(" Stress:   xx          yy          zz          yz          xz          xy\n")
            f.write("            {0}\n".format("    ".join(["{0: .8f}".format(i) for i in virial])))

        if config_id is None:
            conf_id = "{0}_{1}".format("".join(np.unique(chem_syms)),len(pos))
        else:
            conf_id = config_id
        f.write(" Feature   conf_id  {}\n".format(conf_id))
        f.write("END_CFG\n\n")

def _cfgd_to_atoms(cfgd, species=None):
    """Converts a CFG dictionary to an atoms object.
    Args:
        cfgd (dict): of a single config extracted by :func:`cfg_to_xyz`.
        species (list): of element names corresponding to the integer species in
          the CFG dictionary.
    """
    from matdb.atoms import Atoms

    lattice = np.array(cfgd["Supercell"]["vals"])
    natoms = cfgd["Size"]["vals"][0][0]
    stressdict = cfgd["PlusStress"]
    stress = {k: v for k, v in zip(stressdict["cols"], stressdict["vals"][0])}
    energy = cfgd["Energy"]["vals"][0][0]

    alabels = cfgd["AtomData"]["cols"]
    positions = []
    forces = []
    types = []

    for entry in cfgd["AtomData"]["vals"]:
        vals = {k: v for k, v in zip(alabels, entry)}
        pos = []
        force = []
        types.append(species[vals["type"]])
        for poslabel in _cfg_pos:
            pos.append(vals[poslabel])
        positions.append(pos)

        for flabel in _cfg_force:
            force.append(vals[flabel])
        forces.append(force)

    aseatoms = ase.Atoms(numbers=types, positions=np.array(positions),
                         cell=lattice)
    aseatoms.calc = SinglePointCalculator(
                                    aseatoms, energy=energy,
                                    forces=np.array(forces),
                                    stress=order_stress(**stress))
    aseatoms.get_total_energy()
    aseatoms.get_forces()
    aseatoms.get_stress()

    result = Atoms()
    result.copy_from(aseatoms)
    result.pbc = True

    prefix = ""
    if "EFS_by" in cfgd["features"]:
        prefix = "{}_".format(cfgd["features"]["EFS_by"][0].lower())

    force_name, energy_name, virial_name = ["{}{}".format(prefix, l) for l in
                                            ["force", "energy", "virial"]]
    result.properties[force_name] = aseatoms.calc.results["forces"].T
    result.params[energy_name] = energy
    result.params[virial_name] = symmetrize(*aseatoms.calc.results["stress"])

    assert result.n == natoms

    return result

def cfg_to_xyz(cfgfile, outfile="output.xyz", config_type=None, species=None):
    """Converts MTP's CFG forrmat to XYZ.
    .. note:: Multiple frames in the CFG file will be converted to multiple
      frames in the XYZ file.
    Args:
        cfgfile (str): path to the file to convert.
        config_type (str): name of the config_type to assign to each
          configuration.
        species (list): of element names corresponding to the integer species in
          the CFG dictionary.
    """

    from matdb.atoms import AtomsList

    configs = []
    cfgd = None
    with open(cfgfile) as f:
        for line in f:
            if line.strip() == '':
                continue

            if "BEGIN_CFG" in line:
                cfgd = {"features": {}}
            elif isinstance(cfgd, dict) and "END_CFG" not in line:
                if _rxcfg.match(line.strip()):
                    if ':' in line:
                        raw = line.strip().split()
                        label = raw[0].rstrip(':')
                        cols = raw[1:]
                        cfgd[label] = {
                            "cols": cols,
                            "vals": []
                        }
                    else:
                        label = line.strip()
                        cfgd[label] = {
                            "vals": []
                        }

                    if "Feature" in label:
                        fvals = label.split()
                        feature = fvals[1]
                        values = fvals[2:]
                        cfgd["features"][feature] = values
                        del cfgd[label]
                else:
                    parsed = map(eval, line.strip().split())
                    cfgd[label]["vals"].append(parsed)
            elif "END_CFG" in line:
                if cfgd is not None:
                    configs.append(cfgd)
                cfgd = None

    result = AtomsList()
    for cfg in configs:
        atoms = _cfgd_to_atoms(cfg, species)
        result.append(atoms)

    dirname = path.dirname(cfgfile)
    result.write(path.join(dirname, outfile))
    return result

def vasp_to_xyz(folder, outfile="output.xyz", recalc=0,
                properties=["species", "pos", "z", "dft_force"],
                parameters=["dft_energy", "dft_virial"],
                config_type=None):
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

    from matdb.atoms import Atoms
    # from os import path, stat

    if not path.isabs(outfile):
        #Convert to absolute path if one wasn't given.
        outfile = path.join(folder, outfile)

    if (path.isfile(outfile)
        and stat(outfile).st_size > 100
        and recalc <= 0):
        return True

    p = ','.join(properties)
    P = ','.join(parameters)
    renames = [("energy", "vasp_energy"), ("force", "vasp_force"),
               ("virial", "vasp_virial")]
    sargs = ["convert.py", "-I", "vasprun.xml", "-p", p, "-P", P, "-f", "xyz"]
    for s, d in renames:
        sargs.append("-n")
        sargs.append(s)
        sargs.append(d)

    sargs.extend(["-o", outfile, "vasprun.xml"])

    from matdb.utility import execute
    execute(sargs, folder, errignore="OMP_STACKSIZE")

    if config_type is not None:
        #We need to load the XYZ file, add the config_type parameter and then
        #save it again. The -e parameter of convert.py should save this, but in
        #our experiments, it doesn't :(.
        a = Atoms(outfile)
        a.params["config_type"] = config_type
        a.write(outfile)

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

def save_dict_to_h5(h5file, dic, path='/'):
    """Saves a nested dictionary to an open hdf5 file.
    Args:
        h5file (file object): the h5 file to be saved to.
        dic (dict): the dictionary to save.
        path (str, optional): the path within the h5 file that the dict will be
            saved to. Default is '/'.
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes, tuple, float, int)):
            h5file[path + key] = item
        elif isinstance(item, np.ndarray):
            if item.ndim==1 and isinstance(item[0], np.ndarray): #pragma: no cover
                # this chunk of code is only ever used by the unit testing suite.
                dt = h5py.special_dtype(vlen=np.float64)
                h5file.create_dataset(path+key, data=item, dtype=dt)
            else:
                h5file[path + key] = item
        elif isinstance(item, dict):
            save_dict_to_h5(h5file, item, path + key + '/')
        elif isinstance(item, list):
            for i in range(len(item)):
                try:
                    save_dict_to_h5(h5file, item[i], path + key + '/')
                except AttributeError:
                    h5file[path + key + '/' + str(i)] = item
        elif item is None:
            continue

        else:
            raise ValueError('Cannot save %s type'%type(item))

def load_dict_from_h5(h5file, path='/'):
    """Reads an open hdf5 file into a dictionary.
    Args:
        h5file (file object): the h5 file to be read.
        path (str, optional): the path within the h5 file presently being
            read. Default is '/'.
    Returns:
        ans (dict): a dictionary containing the contents of the h5 file.
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = load_dict_from_h5(h5file, path + key + '/')
    return ans
