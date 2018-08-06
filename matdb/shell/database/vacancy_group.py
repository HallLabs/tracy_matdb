from ruamel.yaml.util import load_yaml_guess_indent


def group_options(name):
    """Group to create atomic vacancies from a seed configuration.
    Args:
        name(str): default name Vacancy
        ran_seed (hashable):(=1 default) seed for the random number generator
             for index of vacancies selection.
        vac_per_atom (int < 1): The number of vacancies to include per
             atom in the cell. (i.e. 0.1 would be 1 in every 10 atoms.)
        min_index (int):(default=0) Default choice with the same ran_seed would
             produce the same vacancies in each cell.
    """
    try:
        config, ind, bsi = load_yaml_guess_indent(open('matdb.yml'))
    except(OSError, IOError):
        config, ind, bsi = load_yaml_guess_indent(open('input.yml'))
    parameters = {}
    if(name == 'setup'):
        parameters = {'ran_seed': '', 'vac_per_atom': '', 'min_index': ''}
    # ------------------------------------------------------------------------
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
    # Range parameter
    elif(name == 'vac_per_atom'):
        range_message = """
    This message is printed if no value is selected or if the value is outside
the allowed range."""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 1,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 0.05,
                           # (bool) Are higher values allowed?
                           'top_exempt': False,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': "",
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})
    # -------------------------------------------------------------------------
    return parameters
