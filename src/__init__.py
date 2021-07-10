#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 13:43:14 2020

@author: zmg
"""
# import sys
from pathlib import Path
import logging

from config import config

__package_name__ = "ECpy"


class FindExpFolder:
    """Mock in place class to maintain logging functionality"""

    def __init__(self, arg):
        if arg == "VERSASTAT":
            self.DestDir = Path(__file__).parent.resolve()
        else:
            raise ValueError


_format = "%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s"
try:
    import file_py_helper
    from file_py_helper.find_folders import FindExpFolder

    EXP_FOLDER = FindExpFolder("VERSASTAT")

except ImportError:
    _format = _format + " : importerror : "
    print(
        f'File helper import error, please install...\nLocal logger file: {FindExpFolder("VERSASTAT").DestDir}'
    )
except ModuleNotFoundError:
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

# TODO main list
# TODO add logger.config file
