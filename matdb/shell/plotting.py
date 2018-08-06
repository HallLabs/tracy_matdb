#!/usr/bin/env python3.5
"""The Script for the plotting matdb_shell commands. This script is called to
edit the matdb.yml plotting context arguments for different projects and
interact with the global defaults"""
from shell import Choice, Range
from ruamel.yaml.util import load_yaml_guess_indent  # For parsing yaml args.
from cmd2 import Cmd, with_category  # Subclassed to build the shell.
import ruamel.yaml
import readline
import numpy as np
import os  # For saving and loading plottings.
import sys
import importlib
sys.path.insert(0, './plotting')


class ContextShell(Cmd):
    """An Interactive Shell Interface to matdb/plotting.
    The plotting subclass context is accessed as an interface to the
    available groups and matdb.yml plotting option.
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
        self.plotting_log = './plotting/plotting.log'
        if (os.path.isfile(self.plotting_log)):
            pass
        else:
            os.system('touch {}'.format(self.plotting_log))
        super(ContextShell, self).__init__(use_ipython=True)
        readline.read_history_file(self.plotting_log)
        self.homedir = os.getcwd()
        self.calc_groups = []
        self.setting_up = False
        self.calc_index = 0
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
            self.config['plotting']
        except KeyError:
            self.config.insert(len(self.config), 'plotting', [])

    def do_recompute(self, args):
        """
    Choose which, calculations to recompute for the scope of the
plotting being edited. Runs on setup if the plotting was previously
defined.
        """
        self.poutput("Editing an existing plotting. Choose a recompute" +
                     "option for the calculations that may be affected.\n")
        tasks = [self.colorize("Do not recompute.", "bold"),
                 "Recompute only failed calculations.",
                 "Recompute all calculations."]
        selected = self.choice(tasks, "Choose an Option: ")
        try:
            self.config['databases']
            scope = "global"
        except KeyError:
            try:
                self.config['type']
                scope = 'specific'
            except KeyError:
                try:
                    self.config['steps']
                    scope = 'local'
                except KeyError:
                    self.perror('Recompute called from an invalid location.' +
                                self.colorize("TYPE: \'quit\' for MATDB",
                                              "bold"))
                    return
        if(scope == 'global'):
            self.config['databases'].append({'rerun': selected})
        else:
            self.config.insert(0, 'rerun', selected)
        self.poutput(
            "The {} scope rerun parameter has been set.".format(scope))

    def do_plotting(self, name):
        """
        """
        tasks = []
        try:
            group_type = self.config['plotting']['name']
            tasks.append((group_type, self.colorize(
                "Edit {}".format(group_type), 'bold')))
            tasks.append(('n', self.colorize('Change plotting type',
                                             'bold')))
        except KeyError:
            tasks.append(('n', self.colorize('Add New Plotting',
                                             'bold')))
        tasks.append(('r', 'Delete the plotting and Exit.'))
        tasks.append(('q', 'Quit to database list.'))
        selected = self.select(tasks, "Choose a group option: ")
        if(selected == 'n'):
            self._new_group()
        elif(selected == 'r'):
            del self.config['plotting']
            ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                        indent=self.ind,
                                        block_seq_indent=self.bsi)
            self.do_quit('\n')
        elif(selected == 'q'):
            self.do_quit('\n')
        else:
            self.do_setup('\n')

    def _new_group(self):
        """
        """
        self.db_calcs = []
        calcs = os.listdir('./plotting')
        for each in calcs:
            if(each.endswith('_calc.py')):
                self.db_calcs.append(each[:-8])
        tasks = []
        for calc in self.db_calcs:
            tasks.append((calc, self.colorize('Do {}'.format(calc),
                                              'bold')))
        tasks.append(('q', 'Return to Plotting Start.'))
        selected = self.select(tasks, "Choose a group: ")
        if selected == 'q':
            self.poutput('Quitting back to plotting start.')
            self.poutput(self.colorize("TYPE \'plotting\'", 'bold'))
        else:
            self.config['plotting']['name'] = selected
            ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                        indent=self.ind,
                                        block_seq_indent=self.bsi)
            self.db_calc = selected
            self.poutput("Entering Plotting Setup.")
            self.do_setup('\n')

    def do_setup(self, args):
        """
        """
        db_group = self.db_group
        group = importlib.import_module('{}_group'.format(db_group))
        parameters = group.group_options(name='setup')
        tasks = []
        for step in parameters.keys():
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
            if('options' in parameters.keys()):
                self.poutput(getattr(ContextShell, 'do_choice').__doc__)
                self.poutput(parameters['choice_message'])
                self.poutput(self.colorize('TYPE: choice', 'bold'))
            else:
                self.poutput(getattr(ContextShell, 'do_range').__doc__)
                self.poutput(parameters['range_message'])
                self.poutput(self.colorize('TYPE: range', 'bold'))

            group = importlib.import_module('{}_group'.format(self.db_group))
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
            self._write_to_file(args)
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

    def do_range(self, args):
        """
    The range function is callod for any group parameter with a spocific range
of possible values. The options vary for each parameter. Type double tab to
see a list of some options within the range of this parameter.
        """
        db_group = self.db_group
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
            self._write_to_file(args)
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
        index = len(self.config['databases'][self.dat_index]['steps'][
                self.group_index])
        self.config['databases'][self.dat_index]['steps'][
            self.group_index].insert(index, self.step, args)
        ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                    indent=self.ind,
                                    block_seq_indent=self.bsi)
        self.poutput("The {}: {} has been set.".format(self.step, args))

    @with_category('MATDB')
    def do_quit(self, args):
        """Quits the plotting shell."""
        print("\nQuitting Plotting Context to MATDB.\n")
        self.loader.save()
        return -1


if __name__ == '__main__':
    prompt = ContextShell()
    message = ('\n Starting matdb plotting...\n')
    prompt.prompt = '(matdb/plotting) > '
    prompt.cmdloop(message)
