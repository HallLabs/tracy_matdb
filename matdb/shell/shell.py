import readline
import os


class InputError(Exception):
    """
    """
    pass


class Range(object):
    """Provide a range of reasonable parameters for a matdb group input.
    reason and consequences for outside values.
    """
    def __init__(self, bottom, top, bottom_exempt=False, top_exempt=False,
                 exempt_message=None, range_message=None):
        """Args:
            bottom(float): The bottom available or tested value in the
                range.
            top(float): The top value in the range.
            exempt(bool): True if a value outside the range makes sense
                in some cases.
            exempt_message(str): The warning message for out of bounds
                parameter choice. None leads to default message.
            range_message(str): reason or explanation for selected range.
                default message if None.
        """
        self.bottom = bottom
        self.top = top
        self.top_exempt = top_exempt
        self.bottom_exempt = bottom_exempt
        if(exempt_message):
            self.exempt_message = exempt_message
        else:
            self.exempt_message = """WARNING: Using value in untested range."""
        if(range_message):
            self.range_message = range_message
        else:
            self.range_message = """The list of ranges: {} - {} are those
        that MATDB has been setup to handle.""".format(bottom, top)

    def validate(self, value):
        """assert that the user value is within the specified range.
        Args:
            value(input): User input from the shell. Tab completion for
                available range. Incorrect input prompts an explanation
                for why only a certain rasge is reasonable.
            exempt(bool): Where appropriate exempt as True allows a value
                not in the range is processed.
        Return:
            bool: True or False, a False value will prompt an AssertionError.
        """
        value = float(value)
        if(self.top_exempt or self.bottom_exempt
           or value >= self.bottom and value < self.top):
            try:
                assert value >= self.bottom and value < self.top
                return True
            except AssertionError:
                raise InputError(self.range_message)
        else:
            print(self.exempt_message)

    def message(self):
        """Print the range message on tab completion without text.
        """
        print(self.range_message)


class Choice(object):
    """Represents one of a few possible choices that a person could enter for a
    given parameter. This is a simple validation class. The possible or
    sensible parameters must be updated and passed in.
    """
    def __init__(self, options, exempt=None, exempt_message=None,
                 choice_message=None):
        """Read in the available options for this location in the yaml file.
        Args:
            options(list): list of tested values which have been shown to
                perform as expected for the given calculation.
        """
        self.options = options
        if(exempt):
            self.exempt = exempt
        else:
            self.exempt = False
        if(exempt_message):
            self.exempt_message = exempt_message
        else:
            self.exempt_message = """WARNING: Using untested choice."""
        if(choice_message):
            self.choice_message = choice_message
        else:
            self.choice_message = """The list of choices: {} are those
        that MATDB has been setup to handle.""".format(options)

    def validate(self, value):
        """assert that the user value is one of the available options.
        Args:
            value(input): User input from the shell. Tab completion for
                available options. Incorrect input prompts an explanation
                for why only certain options are reasonable.
            exempt(bool): Where appropriate exempt as True allows a value
                not in "options" to be processed.
        Return:
            bool: True or False, a False value will prompt an AssertionError.
        """
        try:
            value = float(value)
        except ValueError:
            value = str(value)
        if(self.exempt is False or value in self.options):
            try:
                assert value in self.options
                return True
            except AssertionError:
                raise InputError(self.choice_message)
        else:
            print(self.exempt_message)
            return True


class Directory(object):
    """Represents a relative directory object. This simple validation class
    accepts the root directory for the yaml argument as root and asserts that
    the specified string is a relative path to a directory from the location.
    """
    def __init__(self, root, exempt=False, exempt_message=None,
                 directory_message=None):
        """Read in the available options for this location in the yaml file.
        Args:
            options(list): list of tested values which have been shown to
                perform as expected for the given calculation.
        """
        if os.path.isdir(root):
            self.root = root
        else:
            self.root = os.getcwd()
            raise InputError('Invalid root directory, using cwd().')
        self.exempt = exempt
        if(exempt_message):
            self.exempt_message = exempt_message
        else:
            self.exempt_message = """WARNING: Using unknown directory."""
        if(directory_message):
            self.directory_message = directory_message
        else:
            self.directory_message = """Specify a relative path to a directory
            from {}.""".format(self.root)

    def validate(self, value):
        """assert that the user value is one of the available options.
        Args:
            value(input): User input from the shell. Tab completion for
                available options. Incorrect input prompts an explanation
                for why only certain options are reasonable.
            exempt(bool): Where appropriate exempt as True allows a value
                not in "options" to be processed.
        Return:
            bool: True or False, a False value will prompt an AssertionError.
        """
        if(self.exempt is False or os.path.isdir(self.root + value)):
            try:
                assert os.path.isdir(self.root + value)
                return True
            except AssertionError:
                raise InputError(self.directory_message)
        else:
            print(self.exempt_message)
            return True


class File(object):
    """Represents a relative path to a file. This simple validation class
    accepts the root directory for the yaml argument as root and asserts that
    the specified string is a relative path to a file from the location.
    """
    def __init__(self, root, file_type='all', type_exempt=False,
                 file_exempt=False, type_exempt_message=None,
                 file_exempt_message=None, file_message=None):
        """Read in the available options for this location in the yaml file.
        Args:
            options(list): list of tested values which have been shown to
                perform as expected for the given calculation.
        """
        self.file_type = file_type
        if os.path.isdir(root):
            self.root = root
        else:
            self.root = os.getcwd()
            raise InputError('Invalid root directory, using cwd().')
        self.type_exempt = type_exempt
        if(type_exempt_message):
            self.type_exempt_message = type_exempt_message
        else:
            self.type_exempt_message = """WARNING: Using unexpected file type.
            """
        self.file_exempt = file_exempt
        if(file_exempt_message):
            self.file_exempt_message = file_exempt_message
        else:
            self.file_exempt_message = """WARNING: Using unknown file path."""
        if(file_message):
            self.file_message = file_message
        else:
            self.file_message = """Specify a relative path to a file
            from {}.""".format(self.root)

    def validate(self, value):
        """assert that the user value is one of the available options.
        Args:
            value(input): User input from the shell. Tab completion for
                available options. Incorrect input prompts an explanation
                for why only certain options are reasonable.
            exempt(bool): Where appropriate exempt as True allows a value
                not in "options" to be processed.
        Return:
            bool: True or False, a False value will prompt an AssertionError.
        """
        if(not self.type_exempt or value.endswith(self.file_type)):
            if(not self.filexempt or os.path.isfile(self.root + value)):
                try:
                    assert os.path.isfile(self.root + value)
                    return True
                except AssertionError:
                    raise InputError(self.file_message)
            else:
                print(self.file_exempt_message)
                return True
        else:
            print(self.type_exempt_message)
            return True


class History_log():
    """Save or load the history file when saving, exiting, or loading
    a project.
    """
    def __init__(self):
        """
        """

    def save(self):
        """save the history log to a file and clear the history
        to avoid duplication.
        """
        readline.write_history_file('{}.log'.format('matdb'))
        readline.clear_history()
        fname = "input.yml"
        return fname

    def load(self):
        """load the history file after clearing any existing history from
        previous projects.
        """
        readline.clear_history()
        readline.read_history_file('{}.log'.format('matdb'))
        fname = 'matdb.yml'
        return fname
