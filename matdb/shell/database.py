#!/usr/bin/env python3.5
"""The Script for the database matdb_shell commands. This script is called to
edit the matdb.yml database context arguments for different projects and
interact with the plotting global defaults"""
from shell import Choice, Range, Directory, File
from ruamel.yaml.util import load_yaml_guess_indent  # For parsing yaml args.
from cmd2 import Cmd, with_category  # Subclassed to build the shell.
import cmd2
import ruamel.yaml
import functools
import readline
import os  # For saving and loading databases.
import numpy as np
import sys
import importlib
sys.path.insert(0, './database')


class ContextShell(Cmd):
    """An Interactive Shell Interface to matdb/database_builder.
    The database_builder subclass context is accessed as an
    interface to the available groups and matdb.yml database
    options.
    """
    def __init__(self):
        """Initialize global options and check for yaml file.
        vars:
            homedir(str): File path to the matdb/shell directory.
            database_log(file path): Load and save the current project state
                by storing the matdb.yml and log files to a project directory.
            multiline_commands(list): List of Cmd functions which receive
                multiple commandline entries.
            fname(str): The name of the project yaml file (always matdb.yml).
            config(yml object): Stores the current matdb.yml settings.
            ind(yml object): Stores the current matdb.yml indent.
            bsi(yml object): Stores matdb.yml Block Sequence Indent.
            contexts(list): List of available sub contexts.
        """
        self.database_log = './database/database.log'
        if (os.path.isfile(self.database_log)):
            pass
        else:
            os.system('touch {}'.format(self.database_log))
        super(ContextShell, self).__init__(use_ipython=True)
        readline.read_history_file(self.database_log)
        self.homedir = os.getcwd()
        self.db_groups = []
        self.setting_up = False
        self.sub_parameter = False
        self.current_dict = {}
        self.dat_index = 0
        self.group_index = 0

        try:
            self.fname = "matdb.yml"
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))
        except (IOError, OSError):
            self.fname = 'input.yml'
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))

    def preloop(self):
        """
        """
        try:
            self.config['databases']
        except KeyError:
            self.config.insert(len(self.config), 'databases', [])
            self.config['databases'].append({'steps': [{'type': ''}],
                                             'name': ''})
            ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                        indent=self.ind,
                                        block_seq_indent=self.bsi)
        self.poutput(self.colorize('TYPE: \'database\' to begin.', 'bold'))

    def do_database(self, name):
        """Choose a name for a new database instance or select an existing
        name to continue editing.
        """
        names = [('n', self.colorize('Add New Database', 'bold'))]
        # Check number of databases created.
        for i in range(len(self.config['databases'])):
            if(not self.config['databases'][i]['name'] == ''):
                names.append((i, self.config['databases'][i]['name']))
        names.append(('q', 'Quit back to MATDB context.'))

        selected = self.select(names, "Choose a database option: ")
        if(selected == 'q'):
            self.do_quit('\n')
            return -1
        elif(not selected == 'n'):
            self.dat_index = int(selected)
            self.poutput("Entering {} database".format(self.config[
                'databases'][self.dat_index]['name']))
            self.do_groups('\n')
        else:
            name = ''
            while(name == ''):
                name = input("Enter the new database name: ")
                if name == '':
                    self.poutput("Name cannot be an empty string.")
            if(self.config['databases'][0]['name'] == ''):
                self.dat_index = 0
                self.config['databases'][0]['name'] = '{}'.format(name)
            else:
                self.dat_index = len(self.config['databases'])
                self.config['databases'].append({'steps': [{'type': ''}],
                                                 'name': '{}'.format(name)})
            self.poutput('{} database created.'.format(name))
            ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                        indent=self.ind,
                                        block_seq_indent=self.bsi)
            self.do_groups('\n')

    @with_category('MATDB/Groups')
    def do_groups(self, args):
        """
        """
        tasks = []
        tasks.append(('n', self.colorize('Add New Group', 'bold')))
        for i in range(len(self.config['databases'][self.dat_index]['steps'])):
            group_type = self.config['databases'][self.dat_index]['steps'][i][
                'type']
            if(not group_type == ''):
                tasks.append((i, group_type))
        tasks.append(('q', 'Quit to database list.'))
        selected = self.select(tasks, "Choose a group option: ")
        if(selected == 'n'):
            step = self.config['databases'][self.dat_index]['steps']
            if(len(step) == 1, step[0]['type'] == ''):
                group_index = 0
            else:
                group_index = len(step)
            self._new_group(group_index)
        elif(selected == 'q'):
            self.poutput('Quitting back to the database selection...')
            self.poutput(self.colorize("TYPE \'database\'", 'bold'))
        else:
            group_index = selected
            group = self.config['databases'][self.dat_index]['steps'][
                group_index]['type']
            group_name = group[:-(len(group)+1)]
            tasks = [('e', 'Edit existing {} parameters'.format(group_name)),
                     ('c', 'Change this group type.')]
            selected = self.select(tasks, "Choose an action: ")
            if(selected == 'e'):
                self.group_index = group_index
                self.do_setup('\n')
            else:
                self._new_group(group_index)

    def _new_group(self, group_index):
        """
        """
        self.db_groups = []
        groups = os.listdir('./database')
        for each in groups:
            if(each.endswith('_group.py')):
                self.db_groups.append(each[:-9])
        tasks = []
        for group in self.db_groups:
            if(not self._check_group(group)):
                tasks.append((group, self.colorize('Do {}'.format(group),
                                                   'bold')))
            else:
                tasks.append((group, 'Add {} instance'.format(group)))
        tasks.append(('q', 'Quit to database list.'))
        selected = self.select(tasks, "Choose a group: ")
        if selected == 'q':
            self.poutput('Quitting back to database selection.')
            self.poutput(self.colorize("TYPE \'database\'", 'bold'))
        else:
            self.config['databases'][self.dat_index]['steps'][group_index][
                'type'] = '{}.{}'.format(selected, selected.capitalize())
            ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                        indent=self.ind,
                                        block_seq_indent=self.bsi)
            self.db_groups = selected
            self.poutput("Entering Group Setup.")
            self.do_setup('\n')

    def _check_group(self, name):
        """
        """
        if('databases' in self.config.keys()):
            database = self.config['databases'][self.dat_index]
            for types in database['steps']:
                if(name == types['type'][:-(len(name)+1)]):
                    return True
        else:
            return False

    def do_setup(self, args):
        """
        """
        db_group = self.db_groups
        group = importlib.import_module('{}_group'.format(db_group))
        parameters = group.group_options(name='setup')
        tasks = []
        for step in parameters.keys():
            # param = group.group_options(name=step)
            # self.pfeedback(param['{}_message'.format(param['type'])])
            if(self._check_parameters(step)):
                step_val = self.config['databases'][self.dat_index]['steps'][
                    self.group_index][step]
                tasks.append((step, 'Update {} = {}'.format(step, step_val)))
            else:
                tasks.append((step, self.colorize('Do {}'.format(step),
                                                  'bold')))
        tasks.append(('q', 'Quit to database list.'))
        selected = self.select(tasks, "Choose a task: ")
        if not selected == 'q':
            self.setting_up = True
            parameters = group.group_options(name=selected)
            self.step = selected
            if(parameters['type'] == 'dictionary'):
                self.do_dictionary('\n')
            elif(parameters['type'] == 'list'):
                self.do_list('\n')
            else:
                self.poutput(getattr(ContextShell, 'do_{}'.format(
                    parameters['type'])).__doc__)
                choice = input(
                    'type y to continue or any other key to return: ')
                if(not choice == 'y'):
                    self.do_setup('\n')
                else:
                    self.poutput(parameters['{}_message'.format(
                        parameters['type'])])
                    self.poutput(self.colorize('TYPE: {}'.format(
                        parameters['type']), 'bold'))
                    group = importlib.import_module('{}_group'.format(
                        self.db_groups))
                    self.parameters = group.group_options(self.step)
        else:
            self.setting_up = False
            self.poutput('Quitting back to database selection.')
            self.poutput(self.colorize("TYPE \'database\'", 'bold'))

    def _check_parameters(self, param):
        """
        """
        group_vals = self.config['databases'][self.dat_index]['steps'][
            self.group_index]
        for val in group_vals:
            if(val == param):
                return True
        return False

    @with_category('MATDB/Database')
    def do_choice(self, args):
        """
    The choose function is callod for any group option with a set of specifed
options. The options vary for each parameter. Type double tab to see the list
of options for a choice parameter.
        """
        parameters = self.parameters
        choice = Choice(options=parameters['options'],
                        exempt=parameters['exempt'],
                        exempt_message=parameters['exempt_message'],
                        choice_message=parameters['choice_message'])
        if(choice.validate(args)):
            try:
                if(float(args) == int(args)):
                    args = int(args)
                else:
                    args = float(args)
            except ValueError:
                args = str(args)
            if(not self.sub_parameter):
                self._write_to_file(args)
            elif(self.sub_parameter == 'dict'):
                self.do_dictionary('\n')
            elif(self.sub_parameter == 'list'):
                self.do_list('\n')
        if(self.setting_up):
            self.do_setup('\n')

    def complete_choice(self, text, line, start_index, end_index):
        """
        """
        if text:
            return [ele for ele in self.parameters['options'] if
                    ele.startswith(text)]
        else:
            return self.parameters['options']

    @with_category('MATDB/Database')
    def do_directory(self, args):
        """
    The directory function is callod for any group option which require a
relative path to a directory from a set location in the file structure.
        """
        parameters = self.parameters
        directory = Directory(root=parameters['root'],
                              exempt=parameters['exempt'],
                              exempt_message=parameters['exempt_message'],
                              directory_message=parameters[
                                  'directory_message'])
        if(directory.validate(args)):
            if(not self.sub_parameter):
                self._write_to_file(args)
            elif(self.sub_parameter == 'dict'):
                self.do_dictionary('\n')
            elif(self.sub_parameter == 'list'):
                self.do_list('\n')
        if(self.setting_up):
            self.do_setup('\n')

    complete_directory = functools.partialmethod(cmd2.Cmd.path_complete,
                                                 dir_only=True)

    @with_category('MATDB/Database')
    def do_file(self, args):
        """
    The directory function is callod for any group option which require a
relative path to a file of a certain type from a set location in the file
structure.
        """
        parameters = self.parameters
        group_file = File(root=parameters['root'],
                          file_type=parameters['file_type'],
                          type_exempt=parameters['type_exempt'],
                          type_exempt_message=parameters[
                              'type_exempt_message'],
                          file_exempt=parameters['file_exempt'],
                          file_exempt_message=parameters[
                                  'file_exempt_message'],
                          file_message=parameters['file_message'])
        if(group_file.validate(args)):
            if(not self.sub_parameter):
                self._write_to_file(args)
            elif(self.sub_parameter == 'dict'):
                self.do_dictionary('\n')
            elif(self.sub_parameter == 'list'):
                self.do_list('\n')
        if(self.setting_up):
            self.do_setup('\n')

    complete_file = cmd2.Cmd.path_complete

    @with_category('MATDB/Database')
    def do_dictionary(self, args):
        """
    The dictionary function is callod for any group option with a set of
specifed keys. The user will be prompted to submit values for each keys
or accept defaults where applicable.
        """
        # This is called after self.parameters has been set in setup
        # Includes the type, dict_message, key_names, defaults list.
        parameters = self.parameters
        # Print the dictionarf messagi iff quiet is set to False.
        self.pfeedback(parameters['dict_message'])
        # select list of keys to be edited and current values / state.
        options = []
        # For each parameter available in this dictionary.
        for key in parameters['key_names']:
            # Check if the value for each given key has been set.
            if(key not in self.current_dict.keys()):
                options.append((key, "Do {}".format(key)))
            else:
                # For values that have been set list the current value.
                options.append((
                    key, "Update {} = {}".format(key,
                                                 self.current_dict[
                                                     key])))
        # List dictionary options.
        options.append(('q', 'Save and exit dictionary.'))
        selected = self.select(options, "Choose an option: ")
        if(selected == 'q'):
            self._write_to_file(self.current_dict)
            self.current_dict = {}
        else:
            group = importlib.import_module(
                '{}_group'.format(self.db_groups))
            parameters = group.group_options(name=selected)
            self.step = selected
            self.poutput(getattr(ContextShell, 'do_{}'.format(
                parameters['type'])).__doc__)
            self.poutput(parameters['{}_message'.format(
                parameters['type'])])
            self.poutput(self.colorize('TYPE: {}'.format(
                parameters['type']), 'bold'))
            self.sub_parameter = 'dict'

    @with_category('MATDB/Database')
    def do_list(self, args):
        """
    The list function is called for any group option with a set of
specifed keys. The user will be prompted to submit values for each keys
or accept defaults where applicable.
        """
        if(self.parameters['type'] == 'list'):
            self.list_param = self.parameters
            options = []
        else:
            if(self.list_param['arg_names'][0] == 'length'):
                self.list_length = int(self.list_options)
                for i in ['arg_names', 'arg_format', 'arg_min', 'arg_max']:
                    self.list_param[i].pop(0)
            else:
                options.append(self.list_options)
        if(len(options) <= len(self.list_param['arg_names'])):
            arg = self.list_param['arg_names'][len(options)]
            group = importlib.import_module(
                '{}_group'.format(self.db_groups))
            parameters = group.group_options(name=arg)
            self.poutput(getattr(ContextShell, 'do_{}'.format(
                parameters['type'])).__doc__)
            self.poutput(parameters['{}_message'.format(
                parameters['type'])])
            self.poutput(self.colorize('TYPE: {}'.format(
                parameters['type']), 'bold'))
            self.sub_parameter = 'list'
            self.parameters = parameters
        else:
            index = 0
            for form in self.list_param['arg_format']:
                if(form is not ','):
                    joined = []
                    for i in range(len(options[index])):
                        joined.append(form.join(list(map(list, zip(
                            options[index], options[index + 1])))[i]))
                    options.pop(index)
                    options.pop(index + 1)
                    options.insert(0, joined)
                index += 1

    def do_range(self, args):
        """
    The range function is called for any group parameter with a specific range
of possible values. The options vary for each parameter. Type double tab to
see a list of some options within the range of this parameter.
        """
        db_group = self.db_groups
        group = importlib.import_module('{}_group'.format(db_group))
        parameters = group.group_options(self.step)
        self.parameters = parameters
        param_range = Range(bottom=parameters['bottom'], top=parameters['top'],
                            top_exempt=parameters['top_exempt'],
                            bottom_exempt=parameters['bottom_exempt'],
                            exempt_message=parameters['exempt_message'],
                            range_message=parameters['range_message'])
        try:
            args = int(args)
        except ValueError:
            args = float(args)
        if(param_range.validate(args)):
            if(self.sub_parameter is False):
                self._write_to_file(args)
            elif(self.sub_parameter == 'dict'):
                self.do_dictionary('\n')
            elif(self.sub_parameter == 'list'):
                self.do_list('\n')
        if(self.setting_up):
            self.do_setup('\n')

    def complete_range(self, text, line, start_index, end_index):
        """
        """
        values = np.arange(self.parameters['bottom'], self.parameters['top'],
                           self.parameters['step'])
        self.poutput(values)
        if text:
            return [val for val in values if
                    val.startswith(text)]
        else:
            return values

    def _write_to_file(self, args):
        """Write the change made to the matdb.yml file.
        """
        self.config['databases'][self.dat_index][
            'steps'][self.group_index].update({self.step: args})
        ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                    indent=self.ind,
                                    block_seq_indent=self.bsi)
        self.poutput("The {}: {} has been set.".format(self.step, args))

    @with_category('MATDB')
    def do_quit(self, args):
        """\nQuits the database shell.\n"""
        readline.write_history_file('database/database.log')
        readline.clear_history()
        self.poutput("\nQuitting Database Context to MATDB.\n")
        return -1


if __name__ == '__main__':
    prompt = ContextShell()
    message = ('\n Starting matdb database_builder...\n')
    prompt.prompt = '(matdb/database_builder) > '
    prompt.cmdloop(message)
