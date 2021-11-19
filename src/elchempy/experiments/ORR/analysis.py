"""

this module analyzes the ORR experiments, including Disk, N2 background subtraction and possible Ring electrode selectivity data

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

from elchempy.experiments.ORR.RRDE.ring_helpers import get_file_from_secondary_electrode

from elchempy.experiments.N2.analysis import N2_Analysis

# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate


## 3rd party
import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

# constants
from elchempy.constants import EvRHE

#%%
def _dev_test_ORR_analysis():
    """function for testing"""


class ORR_Analysis(ElChemData, DataSelection):
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
        self,
        filepath: [Path, str],
        ring_file=None,
        N2_background_file=None,
        N2_background_data=pd.DataFrame(),
        **kwargs,
    ):
        # self.filepath = Path(filepath, **kwargs)
        # self.kwargs = kwargs
        # self.data = None
        super().__init__(filepath, **kwargs)

        # Load disk data from given filepath
        N2_BG = check_for_N2_BG()
        # self.disk_ecd = ElChemData(filepath, **kwargs)

        # Check for ring file in case of RRDE measurement
        self.ring_file = ring_file

        if not self.ring_file:
            self.ring_file = get_file_from_secondary_electrode(filepath)

        if self.ring_file:
            self.ring_ecd = ElChemData(ring_file, **kwargs)
        else:
            self.ring_ecd = None

        # self.ORR = self.select_data(self.data)
        if 0:
            try:
                pars, data = N2_Cdl_calculation(self.ORR, potential_key=EvRHE)
            except Exception as exc:
                logger.error(f"N2 Cdl calculations failed for {self.filepath}\n{exc}")
                # raise exc from exc
                Cdl_pars, Cdl_data = None, None
            self.Cdl_pars, self.Cdl_data = Cdl_pars, Cdl_data

            self.N2_BG = get_N2_background_data(self.N2_CVs)
        # N2_results = N2_Results(N2_CVs, Cdl_pars, Cdl_data, N2_BG)
        # self.N2_results = N2_analysis.get_N2_analysis_results(self.N2_CVs)

    def plot_ORR(self):
        """plotting fuctions"""

        if self.pars.empty:
            logger.warning(f"N2_results is empty {self.filepath.name}")
            return

        self.pars.groupby("SweepType").plot(x=EvRHE, y="lin_slope", kind="scatter")

    def _test_plot_scanrates(self):

        if self.data.empty:
            logger.warning(f"N2_results is empty {self.filepath.name}")
            return
