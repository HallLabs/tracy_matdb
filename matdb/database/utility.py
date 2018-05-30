"""Contains all the utility functions that belong to the database groups."""

from uuid import uuid4
from cPickle import dump, load
from os import path, rename, remove
import numpy as np
from glob import glob
from tqdm import tqdm
from copy import deepcopy
from itertools import combinations

from phenum.symmetry import bring_into_cell, _does_mapping_exist, _get_transformations
from phenum.grouptheory import _find_minmax_indices

from matdb import msg
from matdb.atoms import AtomsList
from matdb.utility import dbcat

def split(atlist, splits, targets, dbdir, ran_seed, dbfile=None, recalc=0):
    """Splits the :class:`matdb.atoms.AtomsList` multiple times, one for
    each `split` setting in the database specification.

    Args:
        atlsit (AtomsList, or list): the list of :class:`matdb.atams.Atoms` objects
          to be split or a list to the files containing the atoms objects.
        splits (dict): the splits to perform.
        targets (dict): the files to save the splits in, these should 
          contain a {} in the name which will be replaced with the split
          name. The dictionary must have the format {"train": file_name, 
          "holdout": file_name, "super": file_name}.
        dbdir (str): the root *splits* directory for the database.
        dbfile (str): the _dbfile for a legacy database.
        ran_seed (int or float): the random seed for the splits (i.e. the controllers
          random seed).
        recalc (int): when non-zero, re-split the data and overwrite any
          existing *.h5 files. This parameter decreases as
          rewrites proceed down the stack. To re-calculate
          lower-level h5 files, increase this value.
    """
    for name, train_perc in splits.items():
        train_file = targets["train"](name)
        holdout_file = targets["holdout"](name)
        super_file = targets["super"](name)
        idfile = path.join(dbdir, "{0}-ids.pkl".format(name))
        
        if (path.isfile(train_file) and path.isfile(holdout_file)
            and path.isfile(super_file)):
            if recalc <= 0:
                return
            else:
                if path.isfile(idfile):
                    with open(idfile, 'rb') as f:
                        data = load(f)
                new_idfile = path.join(self.root,"{0}_{1}-ids.pkl".format(name,data["uuid"]))
                self.parent.uuids[data["uuid"]] = new_idfile
                for fname in [train_file,holdout_file,super_file]:
                    new_name = fname.replace(name,"{0}_{1}".format(name,data["uuid"]))
                    rename(fname,new_name)
                remove(idfile)

        #Compile a list of all the sub-configurations we can include in the
        #training.
        if not isinstance(atlist,AtomsList):
            subconfs = AtomsList(atlist)
        else:
            subconfs = atlist

        if path.isfile(idfile):
            with open(idfile, 'rb') as f:
                data = load(f)
            subconfs = data["subconfs"]
            ids = data["ids"]
            Ntrain = data["Ntrain"]
            Nhold = data["Nhold"]
            Ntot = data["Ntot"]
            Nsuper = data["Nsuper"]
        else:
            Ntot = len(subconfs)
            Ntrain = int(np.ceil(Ntot*train_perc))
            ids = np.arange(Ntot)
            Nhold = int(np.ceil((Ntot-Ntrain)*train_perc))
            Nsuper = Ntot-Ntrain-Nhold
            np.random.shuffle(ids)

            #We need to save these ids so that we don't mess up the statistics on
            #the training and validation sets.
            data = {
                "uuid": str(uuid4()),
                "subconfs": subconfs,
                "ids": ids,
                "Ntrain": Ntrain,
                "Nhold": Nhold,
                "Ntot": Ntot,
                "Nsuper": Nsuper,
                "ran_seed": ran_seed
            }
            with open(idfile, 'wb') as f:
                dump(data, f)

        #Only write the minimum necessary files. Use dbcat to create the
        #database version and configuration information. There is duplication
        #here because we also store the ids again. We retain the pkl file above
        #so that we can recreate *exactly* the same split again later.
        if not path.isfile(train_file):
            tids = ids[0:Ntrain]
            altrain = subconfs[tids]
            altrain.write(train_file)
            if dbfile is not None:
                dbcat([dbfile], train_file, docat=False, ids=tids, N=Ntrain)
            else:
                dbcat([], train_file, docat=False, ids=tids, N=Ntrain)
        if not path.isfile(holdout_file):
            hids = ids[Ntrain:-Nsuper]
            alhold = subconfs[hids]
            alhold.write(holdout_file)
            if dbfile is not None:
                dbcat([dbfile], holdout_file, docat=False, ids=hids, N=Nhold)
            else:
                dbcat([], holdout_file, docat=False, ids=hids, N=Nhold)
        if not path.isfile(super_file):
            sids = ids[-Nsuper:]
            alsuper = subconfs[sids]
            alsuper.write(super_file)
            if dbfile is not None:
                dbcat([dbfile], super_file, docat=False, ids=sids, N=Nsuper)
            else:
                dbcat([], super_file, docat=False, ids=sids, N=Nsuper)

def dbconfig(dbfile):
    """Returns the database configuration `dict` of the specified database file.

    Args:
        dbfile (str): path to the database file to get config information for.
    """
    from matdb.utility import load_datetime
    
    confpath = dbfile + ".json"
    if not path.isfile(confpath):
        return {}
    
    import json
    with open(confpath) as f:
        config = json.load(f, object_pairs_hook=load_datetime)

    return config


def parse_path(root,seeds,ran_seed=None):
    """Finds the full path to the seed files for this system.
    Args:
        root (str): the root directory for the databes.
        seeds (str or list of str): the seeds for the database that 
            need to be parsed and have the root folder found.
        ran_seed (optional): the seed for the random generator if used.
    
    Returns:
        seed_files (list): a list of the seed files for the database.
    """
    from matdb.utility import special_values
    from itertools import product
    
    seed_files = []
    for seed in seeds:
        # if there is a '/' in the seed then this is a path to the
        # seed file configuration, otherwise we assume the user put
        # the seed file in the `seed` directory in root.
        if len(seed.split("/")) > 1:
            this_seeds = []
            seed_path = root
            for segment in seed.split("/"):
                if "*" in segment:
                    if len(this_seeds) >=1:
                        this_level = []
                        for t_path in res:
                            this_level.extend(glob(path.join(seed_path,t_path,segment)))
                    else:
                        this_level = glob(path.join(seed_path,segment))
                else:
                    this_level = [segment]
                if len(this_seeds) >= 1:
                    this_seeds.extend([path.join(*i) for i in product(this_seeds,this_level)])
                else:
                    this_seeds.extend(this_level)                    

        else:
            seed_path = path.join(root,"seed")
            to_parse = seed
            if "*" in to_parse:
                this_seeds = glob(path.join(seed_path,to_parse))
            else:
                this_seeds = [to_parse]

        for ts in this_seeds:
            t_seed = path.join(seed_path,ts)
            if path.isfile(t_seed):
                seed_files.append(t_seed)
            else:
                msg.err("The seed file {} could not be found.".format(t_seed))

        return seed_files

def make_primitive(atm, eps=None):
    """Finds the primitive cell of a given crystal and the HNF needed to
    convert the primitive to the current crystal.

    Args:
        atm (matdb.Atoms): an atoms object.
        esp (float): floating point tollerance.

    Returns:
        The lattice and atomic basis vectors of the primitive cell along
        with the HNF needed to map the primitive cell to the original cell.
    """
    #extract data from the atoms object to be used.
    a_vecs = atm.cell
    atom_pos = atm.positions
    atom_type = None
    new_vecs = None
    unique_pos = None
    unique_types = None
    try:
        atom_type = atm.get_chemical_symbols()
        atom_type = [i for i in atom_type if i != "X"]
        # mapping = {k:v for v, k in enumerate(np.unique(temp))}
        # atom_type = [mapping[i] for i in types]
    except: #pragma: no cover
        atom_type = atm.numbers
        # mapping = {k:v for v, k in enumerate(np.unique(temp))}
        # atom_type = [mapping[i] for i in types]

    if atom_type is None or len(atom_type) == 0:
        raise ValueError("The atoms object doesn't contain species information. "
                         "The primitive cell can't be found without the species information.")
    if eps is None:
        eps = 1E-3

    #Armed with the data we can now make the cell primitive.
    num_atoms = len(atom_pos)
    
    latt_to_cart, cart_to_latt = _get_transformations(a_vecs)

    #Ensure that all the basis atoms are in the unit cell.
    atom_pos = np.array([bring_into_cell(pos, cart_to_latt, latt_to_cart, eps) for pos in atom_pos])

    #If the cell isn't primitive then there will be lattice vectors
    #inside the cell. If a lattice vector is inside the unit cell then
    #it will bring any atom inside the cell back unto itself via a
    #fractional transaltion. Such a translation must exist for each
    #atom of the same type.
    fracts = []
    for i_atom, a_type in enumerate(atom_type):
        #Only proceed for atoms of the same type as the first type of
        #atom that aren't the first atom.
        if not a_type == atom_type[0] or i_atom == 0:
            continue

        fract = atom_pos[i_atom] - atom_pos[0]
        fract = np.array(bring_into_cell(fract, cart_to_latt, latt_to_cart, eps))
        for j_atom, this_type in enumerate(atom_type):
            #Find the new location of the atom after the fractional
            #translation.
            new_pos = fract + atom_pos[j_atom]
            new_pos = np.array(bring_into_cell(new_pos, cart_to_latt, latt_to_cart, eps))
            mapped = _does_mapping_exist(new_pos, this_type, atom_pos, atom_type, eps)
            if not mapped:
                break

        #If the loop was successfull (mapped==True) then this
        #translation takes all atoms to another of the same type so it
        #is a valid fractional translation and should be kept.
        if mapped:
            fracts.append(list(fract))

    #If the lattice isn't primitive then extra lattice points, i.e.,
    #fractional rotations, were found above.
    if len(fracts) > 0:
        #Collect all lattices points, i.e., potential new primitive vectors.
        lattice_points = deepcopy(fracts)
        lattice_points.extend(a_vecs)

        #Consider all possible triplets of the lattices points to see
        #if the form a set of primitive basis vectors. A triplet will
        #form a basis set if all lattice points will be integer
        #combinations of the triplet.
        for new_vecs in combinations(lattice_points, 3):
            try:
                cart_to_latt = np.linalg.inv(np.transpose(new_vecs))
            except: #pragma: no cover
                continue
            for point in lattice_points:
                vec = np.matmul(cart_to_latt, point)

                #Check if the new vecs are all integers. If not then
                #this lattice point wasn't preserved, move on to the
                #next one.
                mapped = True
                if not np.allclose(vec, np.rint(vec), rtol=eps): #pragma: no cover
                    # I could never get this flag to fire.
                    mapped = False
                    msg.err("Reached portion of code that has never been tested. "
                            "Please submit following to developers for testing: "
                            "a_vecs: {0}, atom_pos: {1}, atom_types: "
                            "{2}".format(a_vecs, atom_pos, atom_type))
                    break

            #If all lattice points were mapped then we've found a valid new basis and can exit.
            if mapped:
                break 

        if not mapped: #pragma: no cover
            raise LogicError("Error in make_primitive. Valid basis not found.")

        #Bring all the atoms into the new cell.
        atom_pos = np.array([bring_into_cell(pos, cart_to_latt, np.transpose(new_vecs), eps)
                             for pos in atom_pos])
        # Check for redundant atoms in the basis. If any are present then remove them.
        unique_pos = []
        unique_types = []
        removed = np.zeros(len(atom_type))
        for i_atom, this_type in enumerate(atom_type):
            pos = atom_pos[i_atom]
            mapped = False
            for j in range(i_atom+1, len(atom_type)):
                if (atom_type[j] == this_type) and (np.allclose(pos, atom_pos[j], rtol=eps)):
                    mapped = True
                    removed[i_atom] = i_atom+1
                    
            if not mapped:
                unique_pos.append(pos)
                unique_types.append(this_type)

    if new_vecs is None:
        new_vecs = a_vecs
    if unique_pos is None:
        unique_pos = atom_pos
    if unique_types is None:
        unique_types = atom_type
    #Now that we have the new atomic and lattice basis we need to
    #define the HNF that mapps the primitive to the supercell.
    n = np.rint(np.matmul(np.linalg.inv(np.transpose(new_vecs)), np.transpose(a_vecs)))
    hnf, b = hermite_normal_form(n)
    hnf = hnf.astype(int)
    return (new_vecs, unique_pos, unique_types, hnf)

def hermite_normal_form(n):
    """Converts an integer matrix to hermite normal form.
    
    Args:
        n (list): a 2D list of integers.

    Returns:
        The integer matrix converted to Hermite Normal Form and 
        the matrix needed to transform the input matrix to the HNF.
    """

    if np.equal(np.linalg.det(n), 0):
        raise ValueError("The input matrix N is singular in hermite_normal_form.")

    hnf = np.array(deepcopy(n))
    b = np.identity(3)

    # Keep doing colum operations until all elements in row 1 are zero
    # except the one on the diagonal.
    count = 0
    while (list(hnf[0]).count(0) != 2):
        count += 1
        min_index, max_index = _find_minmax_indices(hnf[0])
        minm = hnf[0, min_index]
        multiple = hnf[0, max_index]//minm

        hnf[:,max_index] = hnf[:, max_index] - multiple*hnf[:, min_index]
        b[:,max_index] = b[:, max_index] - multiple*b[:, min_index]

        if not np.allclose(np.matmul(n,b), hnf): #pragma: no cover
            raise LogicError("COLS1: Transformation matrix failed in hermite_normal_form")
        
    if hnf[0,0] == 0:
        hnf, b = swap_column(hnf, b, 0)
    if abs(hnf[0,0] < 0):
        hnf[:,0] = -hnf[:,0]
        b[:,0] = -b[:,0]
    if not np.allclose(np.matmul(n,b),hnf): #pragma: no cover
        raise LogicError("COLSWAP1: Transformation matrix failed in hermite_normal_form")

    #Now work on getting hnf[1][2]==0.
    while not hnf[1,2] == 0:
        if hnf[1,1] == 0:
            temp_col = deepcopy(hnf[:,1])
            hnf[:,1] = hnf[:,2]
            hnf[:,2] = temp_col
            
            temp_col = deepcopy(b[:,1])
            b[:,1] = b[:,2]
            b[:,2] = temp_col
        if hnf[1,2] == 0:
            continue

        if (abs(hnf[1,2])<abs(hnf[1,1])):
            max_idx = 1
            min_idx = 2
        else:
            max_idx = 2
            min_idx = 1

        multiple = hnf[1,max_idx]//hnf[1,min_idx]
        hnf[:, max_idx] = hnf[:,max_idx] - multiple*hnf[:, min_idx]
        b[:, max_idx] = b[:,max_idx] - multiple*b[:, min_idx]
        if not np.allclose(np.matmul(n,b), hnf): #pragma: no cover
            raise LogicError("COLS2: Transformation matrix failed in hermite_normal_form")

    if hnf[1,1]<0:
        hnf[:,1] = -hnf[:,1]
        b[:,1] = -b[:,1]
    if not np.allclose(np.matmul(n,b),hnf): #pragma: no cover
        raise LogicError("COLSWAP2: Transformation matrix failed in hermite_normal_form")

    if hnf[2,2]<0:
        hnf[:,2] = -hnf[:,2]
        b[:,2] = -b[:,2]

    if not (hnf[0,1]==0 and hnf[0,2]==0 and hnf[1,2]==0): #pragma: no cover
        print(hnf)
        raise LogicError("hermite_normal_form not lower triangular.")

    if not np.allclose(np.matmul(n,b),hnf): #pragma: no cover
        raise LogicError("End Part 1: Transformation matrix failed in hermite_normal_form")

    #Now that the matrix is lower triangular we need to make sure that
    #the off diagonal elemnts are less than the diagonal elements.

    while (hnf[1,1] <= hnf[1,0] or hnf[1,0]<0):
        multiple = -1
        if hnf[1,1] <= hnf[1,0]:
            multiple = 1
        hnf[:,0] = hnf[:,0] - multiple*hnf[:,1]
        b[:,0] = b[:,0] - multiple*b[:,1]

    for j in [0,1]:
        while (hnf[2,2] <= hnf[2,j] or hnf[2,j]<0):
            multiple = -1
            if hnf[2,2] <= hnf[2,j]:
                multiple = 1
            hnf[:,j] = hnf[:,j] - multiple*hnf[:,2]
            b[:,j] = b[:,j] - multiple*b[:,2]
        
    if not np.allclose(np.matmul(n,b),hnf): #pragma: no cover
        raise LogicError("End: Transformation matrix failed in hermite_normal_form")
        
    if not (hnf[0,1]==0 and hnf[0,2]==0 and hnf[1,2]==0): #pragma: no cover
        raise LogicError("END: hermite_normal_form not lower triangular.")

    if (hnf[1,2]<0 and hnf[2,1]<0 or hnf[2,1]<0): #pragma: no cover
        raise LogicError("END: negative off diagonals (hermite_normal_form).")

    if (hnf[1,0]>hnf[1,1] or hnf[2,0]>hnf[2,2] or hnf[2,1]>hnf[2,2]): #pragma: no cover
        raise LogicError("END: off diagonals larger than diagonals (hermite_normal_form).")
    
    return hnf, b

def swap_column(hnf, b, row):
    """Swaps the non-zero element in the designated row for both hnf and b
    matrices.

    Args:
        hnf (numpy.ndarray): an integer matrix.
        b (numpy.ndarray): an integer matrix.
        row (int): the row that the swap will be centered on.
    
    Returns:
        The hnf and b matrices with their columns swapped so that hnf[row,row] is 
        non-zero.
    """

    min_idx, max_idx = _find_minmax_indices(abs(hnf[row,row:]))
    max_idx += row

    temp_col = deepcopy(hnf[:,row])
    hnf[:,row] = hnf[:,max_idx]
    hnf[:,max_idx] = temp_col

    temp_col = deepcopy(b[:,row])
    b[:,row] = b[:,max_idx]
    b[:,max_idx] = temp_col

    return hnf, b

def decompress(prim, basis, types, hnf_vals):
    """Decompresses the crystal back into it's original form.

    Args:
        prim (list): the primitive lattice vectors as rows of a matrix.
        basis (list): the atomic basis vectors as rows of a matrix.
        types (list): list of integers for the atomic species.
        hnf_vals (list): integer hnf entries.

    Returns:
        The new crystal lattice vectors, atomic basis and atomic types.
    """

    hnf = [[hnf_vals[0], 0, 0], [hnf_vals[1], hnf_vals[2], 0],
           [hnf_vals[3], hnf_vals[4], hnf_vals[5]]]
    lat_vecs = np.transpose(np.matmul(np.transpose(prim), hnf))

    vol_fact = hnf_vals[0]*hnf_vals[2]*hnf_vals[5]

    latt_to_cart, cart_to_latt = _get_transformations(np.transpose(lat_vecs))
    eps = 1E-3
    new_basis = []
    new_types = []
    prim = np.array(prim)
    for a in range(hnf_vals[0]):
        for b in range(hnf_vals[2]):
            for c in range(hnf_vals[5]):
                #calculate the vector that will point to a new atom in
                #the basis by taking a linear combination of the
                #primitive cell vectors.
                add_vec = prim[0]*a + prim[1]*b + prim[2]*c
                for old_t, old_b in zip(types, basis):
                    new_b = list(np.array(old_b)+add_vec)
                    new_basis.append(bring_into_cell(new_b, latt_to_cart, cart_to_latt, eps))
                    new_types.append(old_t)

    if vol_fact*len(basis) != len(new_basis): #pragma: no cover
        raise ValueError("Error occured in decompression.")
                    
    return lat_vecs, new_basis, new_types
