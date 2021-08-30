"""
Created on Tue Jul 20 13:18:07 2021

@author: DW
"""

import pandas as pd


import elchempy

from elchempy.experiments.dataloader.fetcher import ElChemData
from elchempy.experiments.N2.analyses import N2_analysis


class Analyze:  # TODO implement this class or move to BaseRunner class
    """
    Main class for N2 experiment analysis
    """

    def __init__(
        self,
        index_slice: pd.DataFrame = pd.DataFrame(),
        # multi_par_fit=False,
        **N2_kwargs,
    ):

        self.index_slice = index_slice
        # self.multi_par_fit

    def get_data_collection(self):
        results = []
        # for filepath in files:
        while True:
            try:
                filepath = next(filesgen)

                results.append(ElchemData(filepath))
            except StopIteration:
                print(f"data fetch finished len {len(results)}")
                break
        return results

    def run_parallel(self):
        pass

    def run_serial(self):
        pass

    def multi_par_fit(self):
        pass

    def _runner(datacollection):

        _results = []
        for ecdata in datacollection:
            ecdata = N2_analysis(ecdata)
            _results.append(ecdata)
