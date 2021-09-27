"""
Created on Mon Sep 27 23:05:28 2021

@author: DW
"""
"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

## std lib
from typing import NamedTuple, Tuple, Dict
from collections import namedtuple
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local
import elchempy

# from elchempy.dataloaders.fetcher import ElChemData
from elchempy.experiments.N2.analysis import N2_analysis
from elchempy.experiments.N2.background_scan import (
    contains_background_scan,
    get_N2_background_data,
)
from elchempy.experiments.N2.calculations import N2_Cdl_calculations

from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

from elchempy.experiments._dev_datafiles._dev_fetcher import (
    get_files,
    _dev_test_read,
)

## 3rd party
import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

## constants
EvRHE = "E_vs_RHE"

#%%


class DevFuncs:
    """contains developer functions for testing"""

    def _test_data(exp_type: str):

        datacollection = _dev_test_read(get_files(exp_type))
        return datacollection

    def _test_runner():

        _results = []
        for ecdata in get_files("N2"):
            ecdata = N2_analysis(ecdata)
            _results.append(ecdata)
        return _results

    def test_single_file():

        filepath = get_files("N2_20cls_300_100_10_DW28_298")
        n2 = N2_analysis(filepath[0])
        return n2


def N2_runner():
    _result = []
    for fl in get_files("N2"):
        N2res = N2_analysis(fl)
        _result.append(N2res)
    return _result
