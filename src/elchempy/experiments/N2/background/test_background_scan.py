"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

## std lib
# from typing import NamedTuple, Tuple, Dict
# from collections import namedtuple
# from pathlib import Path

import logging

logger = logging.getLogger(__name__)

import unittest
import pytest

## local
# import elchempy

# from elchempy.dataloaders.fetcher import ElChemData
# from elchempy.experiments.N2.analysis import N2_Analysis
from elchempy.experiments.N2.background.background_scan import N2_Background

# from elchempy.experiments.N2.calculations import N2_Cdl_calculation

# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

# for developing
from elchempy.experiments._dev_datafiles._dev_fetcher import (
    get_files,
    _dev_test_read,
)

## 3rd party
# import numpy as np
import pandas as pd

## constants
EvRHE = "E_vs_RHE"

#%%


class Test_N2_Analysis(unittest.TestCase):
    """contains developer functions for testing"""

    def setUp(self):
        bad_BG_filepath = get_files("N2_20cls_300_100_10_DW28_298")[0]
        good_BG_filepath = get_files("N2_20cls_300_100_10_DW28_251")[0]
        self.filepath = good_BG_filepath
        self.N2_BG = N2_Background(good_BG_filepath)

        self.bad_N2_BG = N2_Background(bad_BG_filepath)
        self.bad_filepath = bad_BG_filepath

    def test_N2_setup(self):
        self.assertEqual(self.filepath, self.N2_BG.filepath)

    def test_BG_scan(self):

        self.assertTrue(hasattr(self.N2_BG, "N2_BG_scan"))

        if isinstance(self.N2_BG.N2_BG_scan, pd.DataFrame):
            self.assertFalse(self.N2_BG.N2_BG_scan.empty)
        self.assertEqual(28, self.N2_BG.N2_BG_scan["Segment #"].unique()[0])


# def skip():
#     def _not_test_data(exp_type: str):

#         datacollection = _dev_test_read(get_files(exp_type))
#         return datacollection

#     def _not_test_N2_batch(self):

#         _results = []
#         for filepath in get_files("N2"):
#             n2 = N2_Analysis(filepath)
#             _results.append(n2)
#         return _results


# def _plot():
#     [
#         i.N2_BG.plot(x=EvRHE, y="j_A_cm2")
#         for i in _results
#         if isinstance(i.N2_BG, pd.DataFrame)
#     ]


def _dev():
    self = Test_N2_Analysis()


# def N2_runner():
#     _result = []
#     for fl in get_files("N2"):
#         N2res = N2_Analysis(fl)
#         _result.append(N2res)
#     return _result


if __name__ == "__main__":
    unittest.main()
