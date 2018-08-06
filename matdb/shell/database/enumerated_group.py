
import os

from numpy import zeros as np_zeros

from ruamel.yaml.util import load_yaml_guess_indent


def group_options(name):
    # < Update DocString > (2)
    """Sets up the calculations for a random sampling of structures from
    an enumerated list.

    Args:
        sizes (list): a list containing the smallest and larges cell sizes to
            include in the database.
        lattice (list or str): either the name of the lattice to use
            ('sc','fcc','bcc', 'hcp') or the atomic basis vectors to use.
        basis (list, optional): the atomic basis to use for the enumeration.
            Defaults to [0,0,0] the origin.
        concs (list, optional): the concentrations of each atomic species.
            Defaults to None.
        arrows (list, optional): the maximum number of atoms to displace.
            Defaults ot None.
        eps (float, optional): floating point tolerance for comparisons.
            Defaults to 1E-3.
        ran_seed (int or float, optional): a seed to feed to the random number generator.
            Defaults to None.
        rattle (float, optional): the amount to rattle the atoms by. Defaults to 0.0.
        keep_supers (bool, optional): True if the superperiodic cells are to be kept
            in the enumerated list. Defaults to False.
        displace (float, optional): the amount to displace atoms with arrows. Defaults
            to 0.0.

    .. note:: Additional attributes are also exposed by the super class
      :class:`Group`.

    Attributes:
       max_size (int): the largest allowed cell size in the database.
       min_size (int): the smallest allowed cell size in the database.
       knary (int): the number of atomic species in the enumeration.
       species (list): the atomic species in the system
       arrows (list): list of arrow restrictions.
       concs (list): a list of the concetration restrictions.
       eps (float): the floating point tolerance
       arrow_res (bool): True if arrows are present.
       conc_res (bool): True if concetrations are being restricted.
       lattice (list): the lattice vectors for the system.
       basis (list): the atomic basis for the system.
       rattle (float): the amount to rattle the atoms in the config by.

    """
    # For parameters which need the current user selections config['parameter']
    try:
        config, ind, bsi = load_yaml_guess_indent(open('../matdb.yml'))
    except(OSError, IOError):
        config, ind, bsi = load_yaml_guess_indent(open('input.yml'))
    parameters = {}
    # ----------------------------------------------------------------
    if(name == 'setup'):
        parameters = {
                      'sizes': [], 'lattice': [], 'basis': [],
                      'concs': [], 'arrows': [], 'eps': '',
                      'ran_seed': '', 'rattle': '', 'keep_supers': '',
                      'displace': ''
                      }

    # Range parameter
    elif(name == 'sizes'):
        exempt_message = """ An upper limit of 32 (atoms/cell) will sufficiently
            search over almost every system. However, there are some rare cases
            where a system has more than that number of atoms."""
        range_message = """
            A list containing the smallest and larges cell sizes to
            include in the database.

            The sizes specified must be a list of 1 or 2 values.
            If one value then it must be the largest cell size to include,
            i.e., [32]. If 2 values then it the first value is the smallest
            and the second value is the largest cell size to include,
            i.e., [10,12].

            Upper limit: 32
            Lower limit: 0
            """
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 32,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 2,
                           # (bool) Are higher values allowed?
                           'top_exempt': True,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': True})
    elif(name == 'lattice'):
        exempt_message = """ The coordinate system you have entered is not recognized
            as one of the standard crystal structures. Please enter the desired
            lattice vectors as a list of lists.
            eg. [[1,1,0],
                 [1,0,1],
                 [0,1,1]]
        """
        dict_message = """
            The standard crystal systems are simple cubic (sc), face-centered cubic (fcc),
            body-centered cubic (bcc), and hexagonal close packed (hcp). Please choose one
            of the standard structures or select custom lattice vectors (clv)
            """
        key_names = ['fcc', 'bcc', 'hcp', 'sc']
        defaults = [[]]
        parameters.update({
                           'type':'dictionary',
                           'dict_message': dict_message,
                           'key_names':key_names,
                           'defaults':defaults
                           })
    elif(name == 'basis'):
        options = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
        'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',
        'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
        'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
        'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba',
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
        'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
        'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
        'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
        'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg',
        'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

        choice_message = """ Please list the atomic basis to use for the enumeration.
                The symbols should be in list form, example (ethanol): ['C','H','O']
                Defaults to [0,0,0] the origin. """

        parameters.update({'options': options,  # List of avaliable options.
                           # Are values outside of options allowed?
                           'exempt': False,
                          'exempt_message': exempt_message,
                           'choice_message': choice_message,
                           # (bool) True if accepts multiple arguments
                           'multiple': True})

    elif(name == 'arrows'):
        exempt_message = """"""
        range_message = """
                arrows (list, optional): the maximum number of atoms to displace.
                    Defaults to None.
                The values of 'arrows' are fractional values that will tell us the percent
                fraction of a given element in your system that is allowed to be
                displaced by the 'displace' parameter.
                 
                """
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 1,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 0.05,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': True})

    elif(name == 'concs'):
        list_message = """
    List the concentration of each element as a fraction of 1.
    For example (ethanol): basis = ['C','H','O'],
    concs = [.22, .66, .11]
    """
        # Names for different argument types in a subList
        # (i.e for seeds "calc_type", and "relPath" to a potential)
        arg_names = ['length', 'numerator', 'denominator']
        # The format that exists between argument names.
        # (i.e. for seeds "vasp:Al6Mg4" the first arg_format would be ":")
        # (i.e. for ['name1', name2] the arg_format would be ",")
        arg_format = ['None', ","]
        # If arg min is 1 a default value is expected.
        arg_min = [1, 1, 1]
        # For cases where a limited number of arguments are allowed.
        arg_max = [len(config['species']), 2, 1]
        parameters.update({'type': 'list',
                           'list_message': list_message,
                           'arg_names': arg_names,
                           'arg_format': arg_format,
                           'arg_min': arg_min,
                           'arg_max': arg_max
                           })

    elif(name == 'length'):
        exempt_message = """ ---- """
        range_message = """ Inputing the number of concentration restrictions on the basis atoms.
                Min = 1 (Only one species is restricted)
                Max = number of species (all species have a restriction)

            """
        parameters.update({'type': 'range',
                           'bottom': 1,  # (Float) Lowest available value.
                           'top': len(config['species']),  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 1,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})

    elif(name == 'denominator'):
        exempt_message = """ ---- """
        range_message = """ Input the denominator, or base for the fractional restriction.
                For example:


            """
        parameters.update({'type': 'range',
                           'bottom': 1,  # (Float) Lowest available value.
                           'top': len(config['species']),  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 1,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})

    elif(name == 'numerator'):
        exempt_message = """ ---- """
        range_message = """ Range - Give the numerator(s) for the concentration restriction.

            """
        parameters.update({'type': 'range',
                           'bottom': 1,  # (Float) Lowest available value.
                           'top': 100,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 1,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': True})
    elif(name == 'eps'):
       exempt_message = """ You are choosing a custom floating point tolerance
       different than the default value (default: 1E-3). Be sure you know
       what you're doing before submitting the change.
       """
       range_message = """ (float, optional): floating point tolerance for comparisons.
       Defaults to 1E-3.
       """
       parameters.update({'type': 'range',
       'bottom': .0001,  # (Float) Lowest available value.
       'top': .01,  # (Float) Highest available value.
       # (Float) step for arange tab completion range list.
       'step': .001,
       # (bool) Are higher values allowed?
       'top_exempt': True,
       # (bool) Are lower values allowed?
       'bottom_exempt': True,
       'exempt_message': exempt_message,
       'range_message': range_message,
       # (bool) True if accepts multiple arguments.
       'multiple': False})

    elif(name == 'ran_seed'):
       options = ['integer', 'Random']
       exempt_message = ""
       choice_message = """
       The ran_seed parameter is set to enable both pseudo random config
       generation and scientific reproducibility. Setting the seed to the same value
       with the same parameters will result in the same set of atomic configurations.
       """
       parameters.update({'options': options,  # List of avaliable options.
       # Are values outside of options allowed?
       'exempt': True,
       'exempt_message': exempt_message,
       'choice_message': choice_message,
       # (bool) True if accepts multiple arguments
       'multiple': False})

    elif(name == 'rattle'):
       exempt_message = """
       WARNING: rattle values greater than 1 lead to larger atom positon
       distortions than desired for most cases. To change update the rattle option."""
       range_message = """
       The rattle parameter defines the standord deviation of the atoms in each
       configuration from the lattice vectors. The defined range is 0-1, values
       greater than 1 are allowed with a warning message."""
       parameters.update({'bottom': 0,  # (Float) Lowest available value.
       'top': 1,  # (Float) Highest available value.
       # (Float) step for arange tab completion range list.
       'step': 0.01,
       # (bool) Are higher values allowed?
       'top_exempt': True,
       # (bool) Are lower values allowed?
       'bottom_exempt': False,
       'exempt_message': exempt_message,
       'range_message': range_message,
       # (bool) True if accepts multiple arguments.
       'multiple': False})

    elif(name == 'keep_supers'):
        # < Code to generate function options. >
        exempt_message = """"""
        choice_message = """ (bool, optional): True if the superperiodic cells are to be kept
            in the enumerated list. Defaults to False.

        """
        parameters.update({'options': [],  # List of avaliable options.
                           # Are values outside of options allowed?
                           'exempt': False,
                          'exempt_message': exempt_message,
                           'choice_message': choice_message,
                           # (bool) True if accepts multiple arguments
                           'multiple': False})

    elif(name == 'displace'):
        exempt_message = """ """
        range_message = """ displace (float, optional): the amount to displace atoms with arrows. Defaults
            to 0.0.
            Typically the max value you'll want to give for displacemenet is .1, or
            10%. The range of acceptable values is [0,1]

        """
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 1,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 0.01,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})
    # -------------------------------------------------------------------------
    return parameters
