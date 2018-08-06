from ruamel.yaml.util import load_yaml_guess_indent


def group_options(name):
    """The vasp calculator is a a `matdb` compatible subclass of the
    :class:`ase.calculators.vasp.Vasp` calculator.
    Args:
        param1(Choice): what param1 is and why it uses the choice function.
        param2(Range): what param2 is and why it uses the range function.
        .
        .
        .
        continue for each value avaliable in the group.
    """
    # For parameters which need the current user selections config['parameter']
    try:
        config, ind, bsi = load_yaml_guess_indent(open('matdb.yml'))
    except(OSError, IOError):
        config, ind, bsi = load_yaml_guess_indent(open('input.yml'))
    parameters = {}
    if(name == 'setup'):
        parameters = {'param1': '', 'param2': '', 'param3': '', 'param4': '',
                      'param5': '', 'param6': ''}
    # <For each parameter in the group choose the appropriate template and
    # update the corresponding variables..> (3)
    # -------------------------------------------------------------------------
    # Choice parameter.
    elif(name == 'param1'):
        # < Code to generate function options. >
        exempt_message = """
    Parameter 1 message to be printed if exempt is True and a value not in
options is selected."""
        choice_message = """
    The message to be printed if an incorrect or no entry is made. This
message should explain the function of the value and any relevant information
about the available options."""
        parameters.update({'options': [],  # List of avaliable options.
                           # Are values outside of options allowed?
                           'exempt': True,
                          'exempt_message': exempt_message,
                           'choice_message': choice_message,
                           # (bool) True if accepts multiple arguments
                           'multiple': True})
    # ------------------------------------------------------------------------
    # Range parameter
    elif(name == 'param2'):
        exempt_message = """
    This message is printed if either top or bottom exempt is True and the
user's value falls outside the specified range but is still allowed."""
        range_message = """
    This message is printed if no value is selected or if the value is outside
the allowed range."""
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
                           'multiple': True})
    # -------------------------------------------------------------------------
    # Directory parameter
    elif(name == 'param3'):
        exempt_message = """
    This message is printed if the user input is not an existing directory and
the exempt parameter is true."""
        directory_message = """
    This message describes the use(s) of this parameter and prompts the user
for a relative path to a directory. This message is printed, when quiet is
False, whenever calculator setup begins."""
        parameters.update({'type': 'directory',
                           'root': '../../calculator',  # Path from cwd.
                           'exempt': False,
                           'exempt_message': exempt_message,
                           'directory_message': directory_message})
    # -------------------------------------------------------------------------
    # File parameter
    elif(name == 'param4'):
        file_exempt_message = """
    This message is printed if the user input is not an existing file and
the file_exempt parameter is true."""
        type_exempt_message = """
    This message is printed if the user input is not of the correct file type
and the type_exempt parameter is true."""
        file_message = """
    This message describes the use(s) of this parameter and prompts the user
for a relative path to a file. This message is printed, when quiet is
False, whenever calculator setup begins."""
        parameters.update({'type': 'file',
                           'root': '../../calculator',  # Path from cwd.
                           'file_type': '.py',  # How the desired file ends.
                           'type_exempt': False,
                           'file_exempt': False,
                           'type_exempt_message': type_exempt_message,
                           'file_exempt_message': file_exempt_message,
                           'file_message': file_message})
    # -------------------------------------------------------------------------
    # SubDict parameter
    elif(name == 'param5'):
        dict_message = """
    This message is printed whenever group setup is called if quiet is True. It
is also printed upon selecting the sub dictionary."""
        # Names of the keys in this sub dictionary.
        key_names = ['name1', 'name2', 'name3', 'name4']
        # True is a default value should be used if no value is chosen.
        defaults = [True, True, False, False]
        parameters.update({'type': 'dictionary',
                           'dict_message': dict_message,
                           'key_names': key_names,
                           'defaults': defaults})
    # -------------------------------------------------------------------------
    # SubList parameter
    elif(name == 'param6'):
        list_message = """
    This message is printed whenever group setup is called if quiet is True. It
is also printed upon selecting the sub list."""
        # Names for different argument types in a subList
        # (i.e for seeds "calc_type", and "relPath" to a potential)
        arg_names = ['name1', 'name2', 'name3']
        # The format that exists between argument names.
        # (i.e. for seeds "vasp:Al6Mg4" the first arg_format would be ":")
        # (i.e. for ['name1', name2] the arg_format would be ",")
        arg_format = [":", ","]
        # If arg min is 1 a default value is expected.
        arg_min = [1, 1, 0]
        # For cases where a limited number of arguments are allowed.
        arg_max = [None, None, 1]
        parameters.update({'type': 'list',
                           'list_message': list_message,
                           'arg_names': arg_names,
                           'arg_format': arg_format,
                           'arg_min': arg_min,
                           'arg_max': arg_max})
    # -------------------------------------------------------------------------
    return parameters
