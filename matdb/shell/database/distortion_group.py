from ruamel.yaml.util import load_yaml_guess_indent
import os


def group_options(name):
    """The distortion group takes a seed input configuration and can distort
lattice vectors, cell volume, and atom positions.
    Args:
        name(Default): distortion is name of this group.
        type(Default): distortion.Distorton is the group type.
        seeds(Choice): list of leed configurations, takes multiple arguments
            and must be a list of paths to each seed from the setwd or project
            root directory. The calculator is then used to choose the format.
        rattle(Range): The standord deviation of the atom positions from the
            lattice positions.
        nconfigs(Range): The number of unique configurations to create from
            the sead.
        ran_seed(Choice): either an integer seed or left blank.
        volume_factor(Range): The cell volume parameter adjusts the lattice
            vectors to decrease the overall volume of the unit cell.
        cov_diag(Range): the value on the diagonal of a covariance matrix in
            angstroms. Effects the lattice vectors not atom positions.
        min_index(Choice): Should be 0 or length of a previous database
            calculation with the same setup parameters.
    """
    try:
        config, ind, bsi = load_yaml_guess_indent(open('matdb.yml'))
    except(OSError, IOError):
        config, ind, bsi = load_yaml_guess_indent(open('input.yml'))
    parameters = {}

    if(name == 'setup'):
        parameters = {'seeds': [], 'rattle': '', 'nconfigs': '',
                      'ran_seed': '', 'volume_factor': '',
                      'cov_diag': '', 'min_index': ''}
    # Choice parameter.
    elif(name == 'seeds'):
        if(not config['root'] == ''):
            options = os.listdir('../../../' + config['root'])
        else:
            options = ['path from matdb to seed config.']
        exempt_message = """
    WARNING: Using a seed not found in the root directory. The seed must
include the path from the project or root directory. To change the project
directory update setwd in the MATDB context, to change the seed update seeds.
"""
        choice_message = """
    The seeds parameter is a list of seed configurations for the database
calculations. The seeds must be located in the project directory or have a
path specifying their location from that directory."""
        parameters.update({'type': 'choice',
                           'options': options,  # List of avaliable options.
                           # Are values outside of options allowed?
                           'exempt': True,
                          'exempt_message': exempt_message,
                           'choice_message': choice_message,
                           # (bool) True if accepts multiple arguments
                           'multiple': True})
    # Range parameter
    elif(name == 'rattle'):
        exempt_message = """
    WARNING: rattle values greater than 1 lead to larger atom positon
distortions than desired for most cases. To change update the rattle option."""
        range_message = """
    The rattle parameter defines the standord deviation of the atoms in each
configuration from the lattice vectors. The defined range is 0-1, values
greater than 1 are allowed with a warning message."""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
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
    # Range parameter
    elif(name == 'nconfigs'):
        exempt_message = """
    WARNING: Set to create more than 10,000 configs if this is correct
continue, otherwise update nconfigs to the desired value."""
        range_message = """
    The nconfigs parameter is the number of atomic configurations to create
with unique distortions from the the seed configuration. The defined range is
1-10,000 with higher value available with a warning message."""
        parameters.update({'type': 'range',
                           'bottom': 1,  # (Float) Lowest available value.
                           'top': 10000,  # (Float) Highest available value.
                           # (Float) step for arange tab completion range list.
                           'step': 5,
                           # (bool) Are higher values allowed?
                           'top_exempt': True,
                           # (bool) Are lower values allowed?
                           'bottom_exempt': False,
                           'exempt_message': exempt_message,
                           'range_message': range_message,
                           # (bool) True if accepts multiple arguments.
                           'multiple': False})
    # Range parameter
    elif(name == 'volume_factor'):
        exempt_message = """
    WARNING: A volume_factor greater than 10 will result in a larger unit cell
than desired for most applications. To change update volume_factor."""
        range_message = """
    The volume_factor parameter changes the unit cell volume by changing the
lattice vectors. The defined range is 0-10 with values greater than 10 allowed
with a warning message."""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
                           'top': 10,  # (Float) Highest available value.
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
    # Range parameter
    elif(name == 'cov_diag'):
        exempt_message = """
    WARNING: cov_diag values greater than 1 lead to larger lattice vector
distortions than desired for most applications."""
        range_message = """
    The cov_diag parameter sets the diagonal of a covariance matrix in
angstroms. The lattice vectors are distorted according to this matrixi"""
        parameters.update({'type': 'range',
                           'bottom': 0,  # (Float) Lowest available value.
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
    # Choice parameter.
    if(name == 'ran_seed'):
        options = ['integer', 'Random']
        exempt_message = ""
        choice_message = """
    The ran_seed parameter is set to enable both pseudo random config
generation and scientific reproducibility. Setting the seed to the same value
with the same parameters will result in the same set of atomic configurations.
"""
        parameters.update({'type': 'choice',
                           'options': options,  # List of avaliable options.
                           # Are values outside of options allowed?
                           'exempt': True,
                          'exempt_message': exempt_message,
                           'choice_message': choice_message,
                           # (bool) True if accepts multiple arguments
                           'multiple': False})
    # Choice parameter.
    if(name == 'min_index'):
        options = [0, 'length of previous calculations.']
        exempt_message = """
    WARNING: A non zero min_index is being used. This should only be done to
extend an existing database too avoid duplication."""
        choice_message = """
    The min_index parameter sets the atomic configurations to begin from the
nth pseudo random result from ran_seed. This allows for the extentsion of an
existing database with the same parameters without fear of duplication."""
        parameters.update({'type': 'choice',
                           'options': options,  # List of avaliable options.
                           # Are values outside of options allowed?
                           'exempt': True,
                          'exempt_message': exempt_message,
                           'choice_message': choice_message,
                           # (bool) True if accepts multiple arguments
                           'multiple': False})
    return parameters
