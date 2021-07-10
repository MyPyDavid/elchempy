# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:22:46 2020

@author: DWXMG
"""

print("THIS MODULE IS DEPRECATED", __name__, "\n", __file__)

import logging


# Gets or creates a logger
# logger = logging.getLogger(__name__)
#
## set log level
# logger.setLevel(logging.INFO)
#
## define file handler and set formatter
##    file_handler = logging.FileHandler(FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_DW_logfile.log'))
## create console handler and set level to debug
# ch = logging.StreamHandler()
# ch.setLevel(logging.WARNING)
# formatter    = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s')
# ch.setFormatter(formatter)
## add file handler to logger
# logger.addHandler(ch)
# logger.propagate = False
# logger.warning('=== Started logging {0}... ==='.format(__name__))


def start_logging(modname=__name__):
    # Gets or creates a logger
    logger = logging.getLogger(modname)

    # set log level
    logger.setLevel(logging.INFO)

    # define file handler and set formatter
    #    file_handler = logging.FileHandler(FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_DW_logfile.log'))
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    formatter = logging.Formatter(
        "%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s"
    )
    ch.setFormatter(formatter)
    # add file handler to logger
    logger.addHandler(ch)
    logger.info("=== Started logging {0}... ===".format(modname))
    return logger
