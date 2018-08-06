#!/usr/bin/env python3.5
"""The Script for the database matdb_shell commands. This script is called to
edit the matdb.yml database context arguments for the Distortion group
"""
from shell import Range, Choice
import importlib
import numpy as np
from ruamel.yaml.util import load_yaml_guess_indent  # For parsing yaml args.
from cmd2 import Cmd, with_category  # Subclassed to build the shell.
import ruamel.yaml  # Efficient reading of the input.yml
import os  # For saving and loading databases.
import sys

sys.path.insert(0, './database')


class GroupShell(Cmd):
    """An Interactive Shell Interface to the Distortion group
    Type help to see a list of available commands.
    """
    def __init__(self):
        """Initialize global options and check for yaml file.
        """
        super(GroupShell, self).__init__(
            use_ipython=True, persistent_history_file='database.log')
        self.homedir = os.getcwd()
        self.db_groups = []

        try:
            self.fname = "matdb.yml"
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))
        except (IOError, OSError):
            self.fname = 'input.yml'
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))

        self.db_group = 'distortion'
        self.step = 'rattle'

    def do_setup(self, args):
        """
        """
        db_group = self.db_group
        group = importlib.import_module('{}_group'.format(db_group))
        parameters, yaml = group.group_options(name=step)
        config = ruamel.yaml.load(yaml, Loader=ruamel.yaml.RoundTripLoader)
        for step in parameters.keys():
            

    @with_category('MATDB/Database')
    def do_choice(self, args, db_group=None, step=None):
        """
        """
        group = importlib.import_module('{}_group.py'.format(db_group))
        parameters, yaml = group(name=step)
        self.parameters = parameters
        choice = Choice(options=parameters['options'],
                        exempt=parameters['exempt'],
                        exempt_message=parameters['exempt_message'],
                        choice_message=parameters['choice_message'])
        if(choice.validate(args)):
            print('done')
            # self._write_to_file(args, parameter=step)

    def complete_choice(self, text, line, start_index, end_index):
        """
        """
        if text:
            return [ele for ele in self.parameters['options'] if
                    ele.startswith(text)]
        else:
            return [self.parameters['options']]

    def do_range(self, args, db_group=None, step=None):
        """
        """
        group = importlib.import_module('{}_group.py'.format(db_group))
        parameters, yaml = group(name=step)
        self.parameters = parameters
        param_range = Range(bottom=parameters['bottom'], top=parameters['top'],
                            top_exempt=parameters['top_exempt'],
                            bottom_exempt=parameters['bottom_exempt'],
                            exempt_message=parameters['exempt_message'],
                            range_message=parameters['range_message'])
        if(param_range.validate(args)):
            print('done')
            # self._write_to_file(args, parameter=step)

    def complete_range(self, text, line, start_index, end_index):
        """
        """
        values = np.arange(self.parameters['bottom'], self.parameters['top'],
                           self.parameters['step'])
        if text:
            return [val for val in values if val.startswith(text)]
        else:
            return [values]

    def _write_to_file(self, value, parameter=None):
        """Write the change made to the matdb.yml file.
        """
        self.config[parameter] = '{}'.format(value)
        ruamel.yaml.round_trip_dump(self.config, open('matdb.yml', 'w'),
                                    indent=self.ind,
                                    block_seq_indent=self.bsi)
        print("The {}: {} has been set.".format(parameter,
                                                self.config[parameter]))

    '''
    def do_distortion(self, args):
        import importlib
        try:
            group = importlib.import_module(args, 'GroupShell')
            prompt = group.GroupShell()
            message = ('\n Starting matdb {}...\n'.format(args))
            prompt.prompt = '(matdb/database/{}) > '.format(args)
            prompt.cmdloop(message)
            self.loader.load()
        except (IOError, OSError):
            print("{} is not an available matdb group. ".format(args) +
                  "Choose from the available")
            print("groups: {}".format(self.db_groups))
            print("or type > \"help groups\" to get more information on " +
                  "available options")

    def complete_distortion(self, text, line, start_index, end_index):
        """
        """
        groups = os.listdir('./database')
        for each in groups:
            if(each.endswith('_groups.py')):
                self.db_groups.append(each[:-9])
        if text:
            return [ele for ele in self.db_groups if
                    ele.startswith(text)]
        else:
            return [self.db_groups]
    '''

    def do_quit(self, args):
        """Quits the database shell."""
        print("\nQuitting Distortion group to Database Context.\n")
        self.loader.save()
        return -1


if __name__ == '__main__':
    prompt = GroupShell()
    message = ('\n Starting matdb distortion group...\n')
    prompt.prompt = '(matdb/database_builder/distortion_group/) > '
    prompt.cmdloop(message)
