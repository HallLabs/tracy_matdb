#!/usr/bin/env python3.5
"""The Script for the fitting matdb_shell commands. This script is called to
edit the matdb.yml fitting context arguments for different projects and
interact with the plotting global defaults"""
from shell import History_log
from ruamel.yaml.util import load_yaml_guess_indent  # For parsing yaml args.
from cmd2 import Cmd, with_category  # Subclassed to build the shell.
import os  # For saving and loading fittings.
import sys

sys.path.insert(0, './fitting')


class ContextShell(Cmd):
    """An Interactive Shell Interface to matdb/fitting.
    Type help to see a list of available commands.
    """
    def __init__(self):
        """Initialize global options and check for yaml file.
        """
        super(ContextShell, self).__init__(use_ipython=True)
        self.homedir = os.getcwd()
        self.loader = History_log()
        self.fit_fitters = []

        try:
            self.fname = "matdb.yml"
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))
        except (IOError, OSError):
            os.system("cp input.yml matdb.yml")
            self.fname = 'matdb.yml'
            self.config, self.ind, self.bsi = load_yaml_guess_indent(
                open(self.fname))

    @with_category('MATDB/Fittings')
    def do_fittings(self, args):
        import importlib
        try:
            fitting = importlib.import_module(args, 'FittingShell')
            prompt = fitting.FittingShell()
            message = ('\n Starting matdb {}...\n'.format(args))
            prompt.prompt = '(matdb/fitting/{}) > '.format(args)
            prompt.cmdloop(message)
            self.loader.load()
        except (IOError, OSError):
            print("{} is not an available matdb fitting. ".format(args) +
                  "Choose from the available")
            print("fittings: {}".format(self.db_fittings))
            print("or type > \"help fittings\" to get more information on " +
                  "available options")

    def complete_fittings(self, text, line, start_index, end_index):
        """
        """
        fittings = os.listdir('./fitting')
        for each in fittings:
            if(each.endswith('_fitting.py')):
                self.fit_fittings.append(each[:-3])
        if text:
            return [ele for ele in self.db_fittings if
                    ele.startswith(text)]
        else:
            return self.fit_fittings

    @with_category('MATDB')
    def do_quit(self, args):
        """Quits the fitting shell."""
        print("\nQuitting Fitting Context to MATDB.\n")
        self.loader.save()
        return -1


if __name__ == '__main__':
    prompt = ContextShell()
    message = ('\n Starting matdb fitting...\n')
    prompt.prompt = '(matdb/fitting) > '
    prompt.cmdloop(message)
