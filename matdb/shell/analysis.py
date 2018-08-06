#!/usr/bin/env python3.5
"""The Script for the database matdb_shell commands. This script is called to
edit the matdb.yml database context arguments for different projects and
interact with the plotting global defaults"""
from shell import History_log
from ruamel.yaml.util import load_yaml_guess_indent  # For parsing yaml args.
from cmd2 import Cmd, with_category  # Subclassed to build the shell.
import os  # For saving and loading databases.
import sys

sys.path.insert(0, './database')


class ContextShell(Cmd):
    """An Interactive Shell Interface to matdb/database_builder.
    Type help to see a list of available commands.
    """
    def __init__(self):
        """Initialize global options and check for yaml file.
        """
        super(ContextShell, self).__init__(use_ipython=True)
        self.homedir = os.getcwd()
        self.loader = History_log()
        self.db_groups = []

        try:
            self.fname = "matdb.yml"
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))
        except (IOError, OSError):
            os.system("cp input.yml matdb.yml")
            self.fname = 'matdb.yml'
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))

    @with_category('MATDB/Groups')
    def do_groups(self, args):
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

    def complete_groups(self, text, line, start_index, end_index):
        """
        """
        groups = os.listdir('./database')
        for each in groups:
            if(each.endswith('_group.py')):
                self.db_groups.append(each[:-3])
        if text:
            return [ele for ele in self.db_groups if
                    ele.startswith(text)]
        else:
            return self.db_groups

    @with_category('MATDB')
    def do_quit(self, args):
        """Quits the database shell."""
        print("\nQuitting Database Context to MATDB.\n")
        self.loader.save()
        return -1


if __name__ == '__main__':
    prompt = ContextShell()
    message = ('\n Starting matdb database_builder...\n')
    prompt.prompt = '(matdb/database_builder) > '
    prompt.cmdloop(message)
