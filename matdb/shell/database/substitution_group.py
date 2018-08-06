from ruamel.yaml.util import load_yaml_guess_indent


def group_options(name):
    """A Group to create substitutions in the stoichiometry from
    a seed configuration.
    Args:
        name(str): Default name Substitution
        stoich (list of lists): each list contains the decimal concentration
             of each element in the system where followed by decimal fraction
             of the number of configs which follow this stoichiometry.
             The decimal concentration in each list as well as the decimal
             fraction of nconfigs across all lists must sum to 1.
        ran_seed (int):seed for the random number generator. To allow for
             reproducibility and extension of databases.
        min_index (int):(default=0) For extension of databases set a nonzero
             starting index for extension calculations to avoid repeat values.
        dbargs (dict): dictionary of arguments to be passed to the
            `Group` class.
    """
    try:
        config, ind, bsi = load_yaml_guess_indent(open('matdb.yml'))
    except(OSError, IOError):
        config, ind, bsi = load_yaml_guess_indent(open('input.yml'))
    parameters = {}

    if(name == 'setup'):
        # subList, choice, choice
        parameters = {'stoich': [], 'ran_seed': '', 'min_index': ''}
    # -------------------------------------------------------------------------
    # SubList parameter
    elif(name == 'stoich'):
        list_message = """
    Specify the stoichiometry for the database. The stoichiometry takes the
form of a list of lists. of dimensions (number of elements + 1 by number
of desired stoichiometries.) For example, for Al Ti with 3 unique
stoichiometries: [0.5, 0.5, 0.1] says that 0.1 percent of the database will
contain 50% Al and 50% Ti. The stoichiometry values for each list must sum
ta 1 and the percentage value, the last one, must sum to 1 accross all lists.
(i.e. [[0.3, 0.7, 0.5], [0.4, 0.6, 0.5]] here 0.5 and 0.5 sum to 1.)"""
        arg_names = ['length', 'stoichiometry', 'concentration']
        arg_format = ["None", ","]
        # If arg min is 1 a default value is expected.
        arg_min = [1, 1, 1]
        # For cases where a limited number of arguments are allowed.
        arg_max = [None, None, 1]
        # 2 iff a list is 2 dimensional.
        arg_dim = 2
        parameters.update({'type': 'list',
                           'list_message': list_message,
                           'arg_names': arg_names,
                           'arg_format': arg_format,
                           'arg_min': arg_min,
                           'arg_max': arg_max,
                           'arg_dim': arg_dim})
    # ------------------------------------------------------------------------
    # Range parameter
    elif(name == 'length'):
        exempt_message = ""
        range_message = """
    Specify the number of unique configurations or stoichiometries to use."""
        parameters.update({'type': 'range',
                           'bottom': 1,  # (Float) Lowest available value.
                           'top': 10,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 1,
                           # (bool) Are higher values allowed?
                           'top_exempt': True,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})
    # ------------------------------------------------------------------------
    # Range parameter
    elif(name == 'stoichiometry'):
        range_message = """
    The stoichiometry must be given as at least 1 float that sums to 1.0."""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 1,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 0.1,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': "",
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': True})
    # ------------------------------------------------------------------------
    # Range parameter
    elif(name == 'concentration'):
        range_message = """
    The concentration must be given as a single float less than 1.0."""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 1,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 0.1,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': "",
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})
    # -------------------------------------------------------------------------
    # Range parameter
    elif(name == 'ran_seed'):
        range_message = """
    Choose an integer to set the database seed for reproducibility."""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 100,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 1,
                           # (bool) Are higher values allowed?
                           'top_exempt': True,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': "",
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})
    # -------------------------------------------------------------------------
    # Range parameter
    elif(name == 'min_index'):
        exempt_message = """
    WARNING: Using a nonzero min_index value. This is primarily used to extend
an existing database with the same seed value."""
        range_message = """
    For a set seed value a nonzero index begins picking database
configurations without duplication from that seed index. A zero
min_index with the same seed and nconfigs as a previous database will produce
the same results."""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 1,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 0.1,
                           # (bool) Are higher values allowed?
                           'top_exempt': True,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})
    # -------------------------------------------------------------------------
    return parameters
