"""
Created on Tue Aug 24 12:00:20 2021

@author: DW
"""


import logging

logger = logging.getLogger(__name__)

#%% Logger CLI printing functions


def print_CLI_args(args):
    CLI_command_msg = " ".join(
        [f"--{arg} {value}" for arg, value in vars(args).items()]
    )
    _txt = f"CLI command: {CLI_command_msg}"
    logger.info(_txt)


def cli_debug_settings(args):
    logger.warning("=== Attention debugging mode is now enabled ===")
    args.debug, args.verbose = True, True
    args.exp_type = "N2"
    args.run_mode = "testing"

    return args
