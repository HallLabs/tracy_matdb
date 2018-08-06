#!/usr/bin/env python3.5
"""The Script for the global matdb_shell commands. This script is called to
edit the matdb.yml file for different projects and interact with the
building, calculator, fitting, platting, and analysis defaults."""
from shell import Choice
from subprocess import check_output as check
from ruamel.yaml.util import load_yaml_guess_indent  # For parsing yaml args.
from ase.data import chemical_symbols  # List of Element symbols.
from cmd2 import Cmd, with_category  # Subclassed to build the shell.
import cmd2  # Shell building package.
import ruamel.yaml  # Efficient reading of the input.yml.
import os  # For saving and loading databases.
import readline  # Read to history log.
import functools  # Auto complete the file path.
import argparse  # Parse multiline_commands.
import time


class MatdbShell(Cmd):
    """An Interactive Shell Interface to MATDB. A Subclass of Cmd2. The
    matdb_shell is the base context for the functionality of the shell.
    Subclasses are accessed through workon command.
    """
    def __init__(self):
        """Initialize global options and check for a project yaml file.
        vars:
            homedir(str): File path to the matdb/shell directory.
            loader(func): Load and save the current project state by storing
                the matdb.yml and log files to a project directory.
            multiline_commands(list): List of Cmd functions which receive
                multiple commandline entries.
            fname(str): The name of the project yaml file (always matdb.yml).
            config(yml object): Stores the current matdb.yml settings.
            ind(yml object): Stores the current matdb.yml indent.
            bsi(yml object): Stores matdb.yml Block Sequence Indent.
            contexts(list): List of available sub contexts.
            venvs(list): List of available virtualenvs.
        """
        # Call to Cmd class, allows for ipy embed terminal functionality.
        self.matdb_log = 'matdb.log'
        if (os.path.isfile(self.matdb_log)):
            pass
        else:
            os.system('touch {}'.format(self.matdb_log))
        super(MatdbShell, self).__init__(
            use_ipython=True, persistent_history_file=self.matdb_log)
        self.homedir = os.getcwd()
        self.setting_up = False
        self.multiline_commands = ['species']

        try:
            self.fname = 'matdb.yml'
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))
        except (IOError, OSError):
            self.fname = "input.yml"
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))
        self.contexts = ['database', 'calculator', 'plotting',
                         'fitting', 'analysis']
        # Read in the available virtualenvs
        self.venvs = check("./venvs").decode("utf-8")[:-1].split()

        database_dscrp = """
    database_builder(str): The database builder context controls the species,
lattice types, lattice defects, size, etc... found in the matdb database.
        """
        calculator_dscrp = """
    calculator(str): contols the calculator that will be used. The available
calculators are updated in the context.
        """
        fitting_dscrp = """
    fitting(str): controls the defaults for fitting matdb results. These can
be viewed and edited from within the context.
        """
        analysis_dscrp = """
    analysis(str): The analysis context controls which tests are performed and
displayed. The user may also view, and modify the available datastructures.
        """
        plotting_dscrp = """
    plotting(str): The plotting context controls the plot specifics. Defaults
are updated as the other contexts are changed to include sensible plotting
options. These can be viewed and overridden from within the context.
        """
        self.cont_help = {'database': database_dscrp,
                          'calculator': calculator_dscrp,
                          'plotting': plotting_dscrp, 'fitting': fitting_dscrp,
                          'analysis': analysis_dscrp}

    def preloop(self):
        message = """
    The matdb shell allows for interactive editing of the matdb.yml file. All
of the functionality of the code can be accessed via this shell. type > "help"
to get a list of available functions in the current namespace or type > help
<function> to get specfic information about the use of each function. >
"gettingStarted" has several examples and will give further instructions about
the functionality of the shell.

type "set quiet true" to remove help messages.
        """
        self.pfeedback(message)

    @with_category('MATDB')
    def do_setup(self, args):
        """
    Interactive list of suggested tasks. The numbered list displays the state
of each function.
    If the value has not been updated it is printed in bold. After completing
a task the updated list will be displayed until the user exits.
        """
        matdb_context = ['title', 'species', 'root',
                         'venv']
        defaults = []
        for task in matdb_context:
            if(self.config[task] == ''):
                if(task == 'venv'):
                    task = 'virtual_env'
                if(task == 'root'):
                    task = 'setwd'
                do = (self.colorize('Set {} parameter. '
                                    .format(task), 'bold'))
                defaults.append((task, do))
            else:
                val = self.config[task]
                if(task == 'venv'):
                    task = 'virtual_env'
                if(task == 'root'):
                    task = 'setwd'
                do = ("Update {} = {}".format(task, val))
                defaults.append((task, do))
        for context in self.contexts:
            try:
                self.config[context]
                do = "Return to update {} context.".format(context)
            except KeyError:
                do = self.colorize("Enter {} context.".format(context), 'bold')
            defaults.append((context, do))
        defaults.append(('quit', self.colorize('Quit to manual setup.',
                                               'bold')))
        task = self.select(defaults, "Choose a task: ")
        check = ('type y to continue to the {} context.'.format(task) +
                 ' Or press any key to choose another task: ')
        if(task == 'quit'):
            self.setting_up = False
        else:
            self.setting_up = True
            if(task in self.contexts):
                self.do_help('workon')
                self.poutput(self.cont_help[task])
                choice = input(check)
                if(choice == 'y'):
                    self.do_workon(task)
                else:
                    self.do_setup('\n')
            else:
                self.poutput(getattr(MatdbShell, 'do_' + task).__doc__)
                choice = input(
                    'type y to continue or any other key to return: ')
                if(not choice == 'y'):
                    self.do_setup('\n')
                else:
                    if(task == 'setwd'):
                        os.chdir('../../')
                    self.poutput('\nStarting task')
                    self.poutput(self.colorize(
                        "TYPE: {} <your {}>\n".format(task, task), 'bold'))

    @with_category('MATDB: Project Options')
    def do_setwd(self, args):
        """
    Set the working Directory, project files will be stored in this location.
See help for load_db and save_db for more information on how to save current
projects or continue saved ones.
        """
        try:
            assert os.path.isdir(args)
            self.poutput("Changed Working Directory")
            self._write_to_file(str(args), parameter='root')
        except AssertionError:
            self.poutput("{} is an invalid directory, or navigate."
                         .format(args))
        if(self.setting_up is True):
            self.do_setup('\n')
        os.chdir(self.homedir)

    complete_setwd = functools.partialmethod(cmd2.Cmd.path_complete,
                                             dir_only=True)

    @with_category('MATDB: Project Options')
    def do_title(self, args, parameter='title'):
        """
    Create a title for the current matdb project. e.g. Al-Ti Database. The
title can be changed at any time from the (matdb) context.
    Args:
title(str): The chosen title for the database calculation
        """
        self._write_to_file(str(args), parameter='title')
        if(self.setting_up is True):
            self.do_setup('\n')

    def complete_title(self, text, line, start_index, end_index):
        """Default uses the global then current title for tab completion
        """
        if text:
            return [ele for ele in chemical_symbols if
                    ele.startswith(text)]
        else:
            return [self.config['title']]

    species_parser = argparse.ArgumentParser()
    species_parser.add_argument('elements', nargs='+',
                                help='Elements in system.')

    @cmd2.with_argparser(species_parser)
    @with_category('MATDB: Project Options')
    def do_species(self, args):
        """
    Choose from a list of available elements those elements which will
be included in the current project. (i.e Al Ti)
        """
        species = []
        allowed_elements = Choice(chemical_symbols)
        for element in args.elements:
            allowed_elements.validate(element)
            species.append(element)
        self._write_to_file(species, parameter='species')
        if(self.setting_up is True):
            self.do_setup('\n')

    def complete_species(self, text, line, start_index, end_index):
        """Default uses the global then current title for tab completion
        """
        if text:
            return [ele for ele in chemical_symbols if
                    ele.startswith(text)]
        else:
            return chemical_symbols

    @with_category('MATDB: Project Options')
    def do_virtual_env(self, args):
        """
    Choose a virtual environment to use for running MATDB. tab gives a list of
available virtualenvs.
        """
        # Choice Object, virtenvs must be in choices.
        virtenvs = Choice(self.venvs)
        virtenvs.validate(args)
        # Write validated env argument to the yaml.
        self._write_to_file(args, parameter='venv')
        if(self.setting_up is True):
            self.do_setup('\n')

    def complete_virtual_env(self, text, line, start_index, end_index):
        """Default uses the global then current title for tab completion
        """
        if text:
            return [venv for venv in self.venvs if
                    venv.startswith(text)]
        else:
            return (self.venvs)

    @with_category('MATDB')
    def do_workon(self, args):
        """
    To view and change matdb defaults in any of the available contexts matdb
covers. > workon "context" open a shell where the defaults for the selected
context can be viewed and changed and the available data structures can be
examined.
        """
        import importlib
        try:
            context = importlib.import_module(args, 'ContextShell')
            prompt = context.ContextShell()
            message = ('\n Starting matdb {}...\n'.format(args))
            prompt.prompt = '(matdb/{}) > '.format(args)
            readline.write_history_file(self.matdb_log)
            readline.clear_history()
            prompt.cmdloop(message)
            readline.read_history_file(self.matdb_log)
        except (IOError, OSError):
            self.poutput("{} is not an available option. ".format(args) +
                         "Choose from the available")
            self.poutput("contexts: {}".format(self.contexts))
            self.poutput("or type >(help workon) to get more information on " +
                         "available options")
        if(self.setting_up is True):
            self.do_setup('\n')

    def complete_workon(self, text, line, start_index, end_index):
        """
    Return the available contexts from a list or display contexts that
startswith the typed letters.
        """
        if text:
            return [cont for cont in self.contexts if
                    cont.startswith(text)]
        else:
            return self.contexts

    @with_category('MATDB')
    def do_saveAndClose_db(self, proj_name):
        """
    Save the matdb.yml and history files to a directory. provide a database
name or a path to a directory and name.
        """
        """To save the .yml and log files for scientific reproducibility and
        begin working with a new database or open an existing one.

        save_db(str, abs path): specify the absolute path where the database
            file should be stored.
        name_db(str): Specify a name for the constructed database and
            framework.
        load_db(str, path): To work on an existing database specify the
            absolute path for the existing database.
        name_existing(str): If the existing database has been stored in the
            default directory it can be loaded by speciying the name.
        """
        os.chdir(self.homedir)
        if(os.path.isfile('matdb.yml')):
            saved = []
            for files in os.listdir('./projects'):
                saved.append(files)
            if(not os.path.isdir(proj_name)):
                proj_dir = "./projects/{}".format(proj_name)
            else:
                proj_dir = proj_name
            try:
                os.mkdir(proj_dir)
            except(OSError, IOError):
                msg = ("the database {} already exists.)".format(proj_dir) +
                       "try Using a different name or file path.")
                self.perror(msg, exception_type=OSError)
                self.pfeedback("Type, project name, to save by name.")
                self.pfeedback("Project must NOT exist in: {}".format(saved))
                self.pfeedback("or specify a path to a directory.")
            readline.write_history_file(self.matdb_log)
            os.system("mv matdb.yml matdb.log {}".format(proj_dir))
            readline.clear_history()
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open("input.yml"))
            for i in self.contexts:
                if(os.path.isfile('{}/{}.log'.format(i, i))):
                    os.system("mv {}/{}.log {}".format(i, i, proj_dir))
            self.poutput("saving project {}".format(proj_name))
        else:
            self.perror("No changes to be saved. No changes have been made.")

    complete_saveAndCloe_db = functools.partialmethod(cmd2.Cmd.path_complete,
                                                      dir_only=True)

    @with_category('MATDB')
    def do_load_db(self, args, remove=False):
        """
    load a saved database. If the project was saved to the default location
give the project name. Otherwise navigate to the correct project directory.
        """
        if(os.path.isfile('matdb.yml')):
            remove = input("WARNING: Unsaved changes. go to save_db to save" +
                           "the current project or type(y) to continue and" +
                           "overwrite unsaved changes: ")
            if(remove == 'y'):
                remove = True
        if(remove is True or not os.path.isfile('matdb.yml')):
            os.chdir(self.homedir)
            if(not os.path.isdir(args)):
                proj_dir = "./projects/{}".format(args)
            else:
                proj_dir = args
            for i in self.contexts:
                if(os.path.isfile('{}/{}.log'.format(proj_dir, i))):
                    os.system("mv {}/{}.log {}/".format(proj_dir, i, i))
            os.system("mv {}/matdb.* ./".format(proj_dir))
            os.system("rm -r {}".format(proj_dir))
            readline.clear_history()
            readline.read_history_file(self.matdb_log)
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))
            self.poutput("{} database loaded.".format(args))

    def complete_load_db(self, text, line, start_index, end_index):
        """
Tab completion for saving and opening databases. Selec a name from the saved
list or provide a path to the project directory.
        """
        saved = []
        for files in os.listdir('./projects'):
            saved.append(files)
        if(text):
            return [name for name in saved
                    if name.startswith(text)]
        else:
            return saved

    def _write_to_file(self, value, parameter=None):
        """Write the change made to the matdb.yml file.
        """
        self.config[parameter] = value
        ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                    indent=self.ind,
                                    block_seq_indent=self.bsi)
        self.poutput("\nThe {}: {} has been set.\n".format(parameter,
                                                           self.config[
                                                               parameter]))
        time.sleep(2)


if __name__ == '__main__':
    prompt = MatdbShell()
    message = ('Starting database shell matdb...\n')
    prompt.prompt = '(matdb) > '
    prompt.cmdloop(message)
