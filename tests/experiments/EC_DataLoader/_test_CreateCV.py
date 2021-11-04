"""

tests for dataloaders

"""

## std lib
import sys

from typing import NamedTuple, Tuple, Dict
from collections import namedtuple
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

import unittest
import pytest

## local
import elchempy

from elchempy.dataloaders.fetcher import ElChemData

# from elchempy.experiments.N2.analysis import N2_Analysis
# from elchempy.experiments.N2.background_scan import get_N2_background_data
# from elchempy.experiments.N2.calculations import N2_Cdl_calculation

# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

# for developing and testing
from elchempy.experiments._dev_datafiles._dev_fetcher import (
    get_files,
    _dev_test_read,
)

## 3rd party
# import numpy as np
# import pandas as pd

## constants
EvRHE = "E_vs_RHE"


class TestCreateCV(unittest.TestCase):
    def setUp(self):

        N2_files = (i for i in get_files() if i.name.startswith("N2"))
        N2_file = next(N2_files)
        self.N2_file = N2_file
        self.ECD_empty = ElChemData("")

    def test_ElChemData(self):
        pass

    def test_N2(self):

        try:
            N2_data = ElChemData(self.N2_file)
        except Exception as exc:
            raise exc


def _dev():
    TestCreateCV()


if __name__ == "__main__":
    unittest.main()
