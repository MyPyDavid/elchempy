#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 13:43:14 2020

@author: zmg
"""

# __import__("pkg_resources").declare_namespace(__name__)

import sys
import warnings

__author__ = "David Wallace"
__package_name__ = "elchempy"
__docformat__ = "markdown"
__status__ = "Development"


try:
    from ._version import version

    __version__ = version
except:
    __version__ = "__version__ = '0.0.1'"


from pathlib import Path
import logging

# from elchempy.config import config


class _mock_FindExpFolder:
    """Mock in place class to maintain FindExpFolder functionality"""

    def __init__(self, arg):
        if arg == "VERSASTAT":
            self.DestDir = Path.home().resolve().joinpath(f".{__package_name__}")
            self.PostDir = self.DestDir.joinpath("PostEC")
        else:
            raise ValueError


_format = "%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s"
try:
    import file_py_helper
    from file_py_helper.find_folders import FindExpFolder

    EXP_FOLDER = FindExpFolder("VERSASTAT")

except ImportError:
    FindExpFolder = _mock_FindExpFolder
    _format = _format + " : importerror : "
    print(
        f'File helper import error, please install...\nLocal logger file: {FindExpFolder("VERSASTAT").DestDir}'
    )
except ModuleNotFoundError:
    FindExpFolder = _mock_FindExpFolder
    _format = _format + " : modulenotfounderror: "
    print(
        f'File helper module not found, please install...\nLocal logger file:{FindExpFolder("VERSASTAT").DestDir}'
    )

logger = logging.getLogger(__package_name__)
# set log level
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
# define file handler and set formatter
fh = logging.FileHandler(EXP_FOLDER.DestDir.joinpath("PAR_DW_logfile.log"))
# fh = logging.FileHandler('spam.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter(_format)
# ('%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)
# create console handler and set level to debug
logger.info("=== Started logging {0}... ===".format(__package_name__))
logger.debug("=== Started logging {0}... ===".format(__package_name__))


# This code is written for Python 3.
if sys.version_info[0] != 3:
    logger.error("raman_fitting requires Python 3.")
    sys.exit(1)

# Let users know if they're missing any hard dependencies
hard_dependencies = ("numpy", "pandas", "scipy", "matplotlib", "lmfit")
soft_dependencies = {}
missing_dependencies = []

import importlib

for dependency in hard_dependencies:
    if not importlib.util.find_spec(dependency):
        missing_dependencies.append(dependency)

if missing_dependencies:
    raise ImportError(f"Missing required dependencies {missing_dependencies}")

for dependency in soft_dependencies:
    if not importlib.util.find_spec(dependency):
        warnings.warn(
            f"Missing important package {dependency}. {soft_dependencies[dependency]}"
        )

del hard_dependencies, soft_dependencies, dependency, missing_dependencies


import elchempy

# from elchempy.api import N2_test
## constants
# EvRHE = "E_vs_RHE"

# import experiments
# TODO main list
# TODO add logger.config file
# TODO remove file-py-helper dependency:
#   local Data Dir and destination paths handling
#   local version control and retrieval at runtime
#   local Electrode properties loading at __init__
