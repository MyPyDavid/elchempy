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

from elchempy.dataloaders.fetcher import ElChemData
from elchempy.experiments.N2.background_scan import (
    contains_background_scan,
    get_N2_background_data,
)
from elchempy.experiments.N2.calculations import N2_Cdl_calculations

from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

## 3rd party
import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

## constants
EvRHE = "E_vs_RHE"

#%%


class DevFuncs:
    """contains developer functions for testing"""

    from elchempy.experiments._dev_datafiles._dev_fetcher import (
        get_files,
        _dev_test_read,
    )

    def _test_data(exp_type: str):
        datacollection = _dev_test_read(get_files(exp_type))
        return datacollection

    def _test_runner():

        _results = []
        for ecdata in _test_data("N2"):
            ecdata = N2_analysis(ecdata)
            _results.append(ecdata)
        return _results

    def test_single_file():

        filepath = get_files("N2_20cls_300_100_10_DW28_298")
        n2 = N2_analysis(filepath[0])
        return n2


def new_runner():
    from elchempy.experiments._dev_datafiles._dev_fetcher import (
        get_files,
        _dev_test_read,
    )

    _result = []
    for fl in get_files("N2"):
        N2res = N2_Data(fl)
        _result.append(N2res)
    return _result


class N2_analysis(ElChemData):
    """
    Performs the steps for the N2 analysis on data of a file
    """

    # PAR_exp = 'N2'

    def __init__(self, filepath: [Path, str], **kwargs):
        # self.filepath = Path(filepath, **kwargs)
        # self.kwargs = kwargs
        # self.data = None
        super().__init__(filepath, **kwargs)

        self.N2_CVs = N2_analysis.select_data(self.data)

        Cdl_pars, Cdl_data = N2_Cdl_calculations(self.N2_CVs, EvRHE=EvRHE)
        self.Cdl_pars, self.Cdl_data = Cdl_pars, Cdl_data

        self.N2_BG = get_N2_background_data(self.N2_CVs)
        # N2_results = N2_Results(N2_CVs, Cdl_pars, Cdl_data, N2_BG)
        # self.N2_results = N2_analysis.get_N2_analysis_results(self.N2_CVs)

    @staticmethod
    def select_data(data):
        # FIXME Select only CV types from Data segment
        # Select the data for N2 Cyclic Voltammograms
        try:
            N2_CVs = data.loc[data.ActionId == 38]
            N2_CVs = N2_CVs.dropna(subset=["scanrate"]).loc[N2_CVs.scanrate_calc != 0]
        except Exception as ex:
            logger.error(f"{self} select data error\m{ex}")
            N2_CVs = pd.DataFrame()
        else:
            if N2_CVs.empty:
                logger.warning("select_data is empty, file does not contain any N2 CVs")

        return N2_CVs

    def _test_plot_Cdl(self):

        if self.Cdl_pars.empty:
            logger.warning(f"N2_results is empty {self.filepath.name}")
            return

        self.Cdl_pars.groupby("SweepType").plot(
            x="E_vs_RHE", y="lin_slope", kind="scatter"
        )

    def _test_plot_scanrates(self):

        if self.N2_CVs.empty:
            logger.warning(f"N2_results is empty {self.filepath.name}")
            return

        N2_plot_raw_scans_scanrate(self.N2_CVs)
