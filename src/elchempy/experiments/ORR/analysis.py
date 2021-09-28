"""
Created on Sat Aug 14 23:32:32 2021

@author: DW
"""

from typing import NamedTuple
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

import logging
logger = logging.getLogger(__name__)


import elchempy

from elchempy.experiments.dataloader.fetcher import ElChemData


def new_runner():
    from elchempy.experiments.dataloader._dev_fetcher import get_files, _dev_test_read
    _result = []
    for fl in get_files('O2'):
        N2res = N2_Data(fl)
        _result.append(N2res)
    return _result

# N2_results = namedtuple('N2', 'raw_data pars data N2_BG')

class O2_Results(NamedTuple):
    raw_data: pd.DataFrame
    pars: pd.DataFrame
    data: pd.DataFrame
    N2_BG: pd.DataFrame

EvRHE = "E_vs_RHE"

if 0:
    nn=N2_data('//mnt/DATA/APPS_SOFT/VENVS/repos/elchempy/data/raw/06.03.2018_DW28_HPRR_0.1MHClO4_RRDE22960/N2_20cls_300_100_10_DW28_298.par')

class O2_Data(ElChemData):

    def __init__(self, filepath: [Path, str], **kwargs):
        self.filepath = filepath
        self.kwargs = kwargs
        super().__post_init__()

        EC_data = self.select_data()

        O2_results = self.analyze(EC_data)
        self.add_analysis_method(O2_results)

    def select_data(self):
        # FIXME Select only CV types from Data segment
        # Select the data for N2 Cyclic Voltammograms
        N2_CVs = self.data.loc[self.data.ActionId == 38]
        N2_CVs = N2_CVs.dropna(subset=['scanrate']).loc[N2_CVs.scanrate_calc != 0]
        return N2_CVs

    def analyze(self, N2_CVs):
