#!/usr/bin/python

import tracy_wrapper
import datetime
import json

def examples():
    """Prints examples of using the script to the console using colored output.
    """
    from matdb import msg
    script = "MATDB Tracy Queue submission"
    explain = ("Once a database of configurations has been built and prepped "
               "for submission to the Tracy Queue this script will grab the data "
               " and submit it to the job script.")
    contents = [(("Submit the calculations to the queue."), 
                 "tracy_sub.py submission.json",
                 "This reads in the data and sends it to the endpoint.")]
    required = ("'submission.json' json file with the job specifications.")
    output = ("")
    details = ("")
    outputfmt = ("")

    msg.example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "subspec": {"help": "File containing the submission specifications."},
    "-test_run": {"action": "store_true",
                 "help": ("Allows for test runs that won't submit to the queue.")},
}
"""
dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse
    import sys
    from matdb import base
    pdescr = "MATDB Tracy Queue submission"
    parser = argparse.ArgumentParser(parents=[base.bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)
        
    args = base.exhandler(examples, parser)
    if args is None:
        return

    return args

def submit(subfil, test=False):
    """Submits the file to the queue if test is True.

    Args:
        subfil (str): the json file to be loaded.
        test (bool): True if a test run.
    """

    with open(subfil, "r") as f:
        payload = json.load(f)

    payload["input"] = str(payload["input"])
    payload["dateReady"] = datetime.datetime.now()

    tracy_client = tracy_wrapper.get_client()
    if test:
        tracy_client.api.Contract.add.add(results=payload)
    else: #pragma: no cover
        tracy_client.api.Contract.add.add(results=payload)

def run(args):
    """Runs the matdb setup and cleanup to produce database files.
    """
    if args is None:
        return
    
    #No matter what other options the user has chosen, we will have to create a
    #database controller for the specification they have given us.
    if args["test_run"]:
        submit(args["subspec"], test=True)
    else:
        submit(args["subspec"])
        
if __name__ == '__main__': # pragma: no cover
    run(_parser_options())

