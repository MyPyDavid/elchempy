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
from elchempy.experiments.Logic.selection import DataSelection

from elchempy.experiments.N2.calculations import N2_Cdl_calculation
from elchempy.experiments.N2.background_scan import get_N2_background_data

from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

## 3rd party
import pandas as pd

## constants
EvRHE = "E_vs_RHE"

#%%


class N2_Analysis(ElChemData, DataSelection):
    """
    Inherits from ElChemData and DataSelection
    performs the steps for the N2 analysis on data of the given filepath

    """

    def __init__(self, filepath: [Path, str], **kwargs):
        # self.filepath = Path(filepath, **kwargs)
        # self.kwargs = kwargs
        # self.data = None
        super().__init__(filepath, **kwargs)

        # self.N2_CVs = self.select_data(self.data)
        try:
            Cdl_pars, Cdl_data = N2_Cdl_calculation(
                self.data_selection, potential_key=EvRHE
            )
        except Exception as exc:
            logger.warning(f"N2 Cdl calculations failed for {self.filepath}\n{exc}")
            # raise exc from exc
            Cdl_pars, Cdl_data = None, None
        self.Cdl_pars, self.Cdl_data = Cdl_pars, Cdl_data

        self.N2_BG = get_N2_background_data(self.data_selection)
        # N2_results = N2_Results(N2_CVs, Cdl_pars, Cdl_data, N2_BG)
        # self.N2_results = N2_analysis.get_N2_analysis_results(self.N2_CVs)

    def plot_Cdl(self):

        if not isinstance(self.Cdl_pars, pd.DataFrame):
            logger.warning(
                f"Unable to Plot, N2 Analysis did not calculate the Cdl pars for {self.filepath.name}"
            )
            return

        if self.Cdl_pars.empty:
            logger.warning(f"Unable to Plot, N2_results is empty {self.filepath.name}")
            return

        self.Cdl_pars.groupby("SweepType").plot(
            x="E_vs_RHE", y="lin_slope", kind="scatter"
        )

    def plot_scanrates(self):

        if self.N2_CVs.empty:
            logger.warning(f"N2_results is empty {self.filepath.name}")
            return

        N2_plot_raw_scans_scanrate(self.N2_CVs)
