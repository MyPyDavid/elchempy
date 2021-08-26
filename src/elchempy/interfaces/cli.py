"""
Created on Sat Aug 14 19:33:18 2021

@author: DW
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

# For CLI interface
import argparse

# File handling
from pathlib import Path
from os import PathLike

## Logging ###
import logging

logger = logging.getLogger(__name__)
# print('name', __name__)

logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler("elchempy.log")
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(lineno)s - %(message)s"
)
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


### Local imports
from .cli_helpers import print_CLI_args, cli_debug_settings

### Version
import importlib.metadata


try:
    _version = importlib.metadata.version("elchempy")
except Exception as e:
    _version = "version.not.found"

_version_text = f"\n=== CLI elchempy version: {_version} ===\n"
# print(_version_text)

# Args for CLI
RUN_MODES = ["normal", "testing", "make_index", "make_examples"]

EXP_TYPES = ["N2", "ORR", "HPRR", "EIS"]
SUPPORTED_TYPES = ["N2"]


def main():

    """
    The command line interface for elchempy
    """

    parser = argparse.ArgumentParser(
        description="Command-line interface for elchempy package main."
    )

    parser.add_argument(
        "-M",
        "--run-mode",
        type=str,
        choices=RUN_MODES,
        help="running mode of package, for testing",
        default="normal",
    )

    parser.add_argument(
        "-exp",
        "--exp-type",
        nargs="+",
        default=["any"],
        choices=EXP_TYPES,
        help=f"Selection of experimental types to perform, supported types: {', '.join(SUPPORTED_TYPES)}",
    )

    parser.add_argument(
        "-sIDs",
        "--sampleIDs",
        nargs="+",
        default=[],
        help="Selection of names of SampleIDs from index to run over.",
    )

    parser.add_argument(
        "-sGrps",
        "--samplegroups",
        nargs="+",
        default=[],
        help="Selection of names of sample groups from index to run over.",
    )

    parser.add_argument(
        "-idx",
        "--index",
        action="store_true",
        help="Performs the indexation of experimental files.",
    )
    # Enables debugging mode
    parser.add_argument(
        "-dbg", "--debug", action="store_true", help="enables debugging mode"
    )

    # Execute the parse_args() method
    parser.add_argument(
        "-vv", "--verbose", action="store_true", help="increase output verbosity"
    )

    parser.add_argument(
        "--version",
        # action=print(_version_text),
        action="version",
        version="%(prog)s {}".format(_version),
        # const=_version_text,
        help="Prints out the current version of the raman_fitting distribution, via importlib.metadata.version",
    )

    # Execute the parse_args() method
    args = parser.parse_args()

    # import the raman_fitting package
    import elchempy

    if args.run_mode == "testing":
        args = cli_debug_settings(args)
        # args.debug = True

    # set debugging mode
    if args.debug:
        args = cli_debug_settings(args)
        # args.verbose = True

    # increase verbosity
    if args.verbose:
        ch.setLevel(logging.DEBUG)
        print_CLI_args(args)

    print(f"CLI args: {args}")
    if args.run_mode == "normal":
        pass
        # _org_index = OrganizeRamanFiles()
        # RL = RamanLoop(_org_index, run_mode ='normal')
    elif args.run_mode.upper() == "DEBUG":
        args.run_mode = args.run_mode.upper()
        # TODO Add a FAST TRACK for DEBUG
    elif args.run_mode == "testing":
        from elchempy.experiments.N2.analyses import new_runner

        _result = new_runner()
        _ress = ",\n====\n".join(map(repr, _result))
        print(f"Finished:\n{_ress }")

    # _main_run = rf.MainDelegator(**vars(args))

    # return parser
