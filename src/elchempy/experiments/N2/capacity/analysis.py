"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

# std lib
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local

from elchempy.dataloaders.fetcher import ElChemData
from elchempy.experiments.Logic.selection import DataSelection

from elchempy.experiments.N2.capacity.calculations import N2_Cdl_calculation

# from elchempy.experiments.N2.background.background_scan_selectors import get_N2_background_data

from elchempy.experiments.N2.plotting import N2_Plotter

## 3rd party
import pandas as pd

## constants
from elchempy.constants import EvRHE

#%%


class N2_Analysis(ElChemData, DataSelection, N2_Plotter):
    """
    Inherits from ElChemData and DataSelection
    performs the steps for the N2 analysis on data of the given filepath

    """

    def __init__(self, filepath: [Path, str], calculate_cdl=True, **kwargs):

        super().__init__(filepath, **kwargs)

        self.calculate_cdl = calculate_cdl
        if self.calculate_cdl:
            self.Cdl_pars, self.Cdl_data = self.make_cdl_calculation()

    def make_cdl_calculation(self, potential_key=EvRHE):
        try:
            Cdl_pars, Cdl_data = N2_Cdl_calculation(
                self.data_selection, potential_key=EvRHE
            )
        except Exception as exc:
            # maybe raise exc from exc
            logger.warning(f"N2 Cdl calculations failed for {self.filepath}\n{exc}")

            Cdl_pars, Cdl_data = None, None
        return Cdl_pars, Cdl_data

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
