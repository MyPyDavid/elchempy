"""
Tests for ElChemData fetcher. Fetches data from a file and constructs the electrochemical data columns
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

# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

from elchempy.experiments._dev_datafiles._dev_fetcher import (
    get_files,
    _dev_test_read,
)

## 3rd party
# import numpy as np
# import pandas as pd

## constants
EvRHE = "E_vs_RHE"


from collections import namedtuple
from dataclasses import dataclass, InitVar

from pathlib import Path
from typing import Union

import pandas as pd
import numpy as np

import elchempy
from elchempy.dataloaders.reader import DataReader
from elchempy.dataloaders.fetcher import ElChemData

# from .converters import get_current_density, get_potential_vs_RE, get_RPM_from_DAC_V
from elchempy.dataloaders.assigners import assign_all_EC_relevant_columns


#%%
def _dev():
    from elchempy.experiments._dev_datafiles._dev_fetcher import (
        get_files,
        _dev_test_read,
    )

    N2fls = get_files("N2")

    try:
        ecd = ElChemData(N2fls[0].with_suffix(".noname"))
        ba = ElChemData("a")
    except FileNotFoundError:
        pass
    except Exception as ex:
        print("unexpected error", ex)
        raise ex from ex

    ecd = ElChemData(N2fls[1])
    self = ecd


#%%


class Test_ElChemData(unittest.TestCase):
    """contains developer functions for testing"""

    def setUp(self):
        filepath = get_files("N2_20cls_300_100_10_DW28_298")
        self.filepath = filepath[0]
        self.ecd = ElChemData(filepath[0])

    def test_setup(self):
        self.assertEqual(self.filepath, self.ecd.filepath)

    def test_actions(self):
        self.assertEqual(len(self.ecd.actions), 10)

    def test_data(self):
        self.assertEqual(len(self.ecd.data), 49938)

    def test_metadata_only(self):
        self.ecd_meta = ElChemData(self.filepath, metadata_only=True)

        self.assertEqual(self.ecd_meta.metadata_only, True)
        self.assertEqual(self.ecd_meta.data.empty, True)

    def _test_data(exp_type: str):

        datacollection = _dev_test_read(get_files(exp_type))
        return datacollection

    def _test_batch(self):

        _results = []
        for filepath in get_files("N2"):
            n2 = ElChemData(filepath)
            _results.append(n2)
        return _results


def _plot():
    [
        i.N2_BG.plot(x=EvRHE, y="j_A_cm2")
        for i in _results
        if isinstance(i.N2_BG, pd.DataFrame)
    ]


def _dev():
    self = Test_ElChemData()
    self.setUp()


if __name__ == "__main__":
    unittest.main()
