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

from elchempy.experiments.N2.calculations import N2_Cdl_calculation
from elchempy.experiments.N2.background_scan import get_N2_background_data

from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

## 3rd party
import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

## constants
EvRHE = "E_vs_RHE"

#%%
def _dev_test_ORR_analysis():


class ORR_Analysis(ElChemData):
    """
    Inherits from ElChemData,
    performs the steps for the ORR analysis on data of a file

    Difficulties:
        Information from other files is required for processing,
        there should be a N2 background scan in another file in same folder

        there can be a RRDE measurement with Disk and Ring files, which
        need to be merged on Elapsed Time for the calculations

    Steps:
        N2 background file: look in folder or get file from argument (from index)

    """

    def __init__(
        self, filepath: [Path, str], ring_file=None, N2_background=None ** kwargs
    ):
        # self.filepath = Path(filepath, **kwargs)
        # self.kwargs = kwargs
        # self.data = None
        super().__init__(filepath, **kwargs)

        self.ORR = self.select_data(self.data)

        try:
            Cdl_pars, Cdl_data = N2_Cdl_calculation(self.N2_CVs, potential_key=EvRHE)
        except Exception as exc:
            logger.error(f"N2 Cdl calculations failed for {self.filepath}\n{exc}")
            # raise exc from exc
            Cdl_pars, Cdl_data = None, None
        self.Cdl_pars, self.Cdl_data = Cdl_pars, Cdl_data

        self.N2_BG = get_N2_background_data(self.N2_CVs)
        # N2_results = N2_Results(N2_CVs, Cdl_pars, Cdl_data, N2_BG)
        # self.N2_results = N2_analysis.get_N2_analysis_results(self.N2_CVs)

    def select_data(self, data):
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
