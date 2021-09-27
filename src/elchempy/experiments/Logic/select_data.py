"""
Created on Thu Sep 16 14:46:27 2021

@author: DW
"""
from dataclasses import dataclass

import logging

logger = logging.getLogger(__name__)

import pandas as pd


@dataclass
class ExperimentalDataSelector:
    """
    selects the relevant data for each type of experiment from the overall data,
    by slicing the data out of the pd.DataFrame
    """

    experiment_name: str
    data: pd.DataFrame
    selection = None

    def __post_init__(self):

        if not self.experiment_name:
            return

        if "N2" in self.experiment_name.upper():
            self.selection = self.N2(self.data)

    def N2(self, data):
        """returns Cyclic Voltammograms with scan rate above 0"""

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
