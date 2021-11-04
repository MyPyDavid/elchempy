"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

## std lib
from typing import NamedTuple, Tuple, Dict
from collections import namedtuple
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

import unittest
import pytest

## local
import elchempy

# from elchempy.dataloaders.fetcher import ElChemData
from elchempy.experiments.N2.analysis import N2_Analysis
from elchempy.experiments.N2.background_scan import get_N2_background_data
from elchempy.experiments.N2.calculations import N2_Cdl_calculation

# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

# for developing
from elchempy.experiments._dev_datafiles._dev_fetcher import (
    get_files,
    _dev_test_read,
)

## 3rd party
# import numpy as np
# import pandas as pd

## constants
EvRHE = "E_vs_RHE"

#%%


class Test_N2_Analysis(unittest.TestCase):
    """contains developer functions for testing"""

    def setUp(self):
        filepath = get_files("N2_20cls_300_100_10_DW28_298")
        self.filepath = filepath[0]
        self.n2 = N2_Analysis(filepath[0])

    def test_N2_setup(self):
        self.assertEqual(self.filepath, self.n2.filepath)

    def test_Cdl_pars(self):
        Cdl_pars = self.n2.Cdl_pars

        cdl050 = Cdl_pars.loc[
            (0.49 < Cdl_pars[EvRHE])
            & (Cdl_pars[EvRHE] < 0.54)
            & (Cdl_pars["SweepType"] == "cathodic")
        ]
        cdl050_mean = cdl050["lin_slope"].mean()
        self.assertEqual(cdl050_mean, 0.0001221518868503618)

    def test_BG_scan(self):
        self.assertEqual(27, self.n2.N2_BG["Segment #"].unique()[0])

    def _test_data(exp_type: str):

        datacollection = _dev_test_read(get_files(exp_type))
        return datacollection

    def test_N2_batch(self):

        _results = []
        for filepath in get_files("N2"):
            n2 = N2_Analysis(filepath)
            _results.append(n2)
        return _results


def _plot():
    [
        i.N2_BG.plot(x=EvRHE, y="j_A_cm2")
        for i in _results
        if isinstance(i.N2_BG, pd.DataFrame)
    ]


def N2_runner():
    _result = []
    for fl in get_files("N2"):
        N2res = N2_Analysis(fl)
        _result.append(N2res)
    return _result


if __name__ == "__main__":
    unittest.main()
