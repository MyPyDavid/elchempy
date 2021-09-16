"""
Created on Mon Aug 23 14:24:42 2021

@author: DW
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 15:05:57 2017
@author: David Wallace
This module reads PAR files in a certain folder and tries to analyze the experimental data
from N2, ORR, OER, HER, HPRR and EIS measurements.
"""
# import os
import sys

# sys.setrecursionlimit(3000)

from collections import namedtuple

import numpy as np
import pandas as pd
import re

# import multiprocessing
from pathlib import Path

import datetime as dt

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper import PostChar

if __name__ == "__main__":

    from indexer import prepare_input

    from experiments import run_EIS
    from experiments import run_N2
    from experiments import run_ORR
    from experiments import run_OER
    from experiments import run_HER
    from experiments import run_HPRR

    from experiments.EC_DataLoader.set_OCP_RHE import get_RHE_OCP

    # from PostEC.collect_load import get_EIS_pars
    # from PostEC.EIS_export import EIS_all_check_redchi, export_to_Origin

else:
    from .indexer import prepare_input

import logging

_logger = logging.getLogger(__name__)
# _logger.propagate = False
#%%
# Globally used Variable in every module
