"""
Created on Thu Oct  7 15:20:24 2021

@author: DW
"""
import pytest
import unittest

import elchempy
from elchempy.api import N2_testrun


class TestN2Analysis(unittest.TestCase):
    def setUp(self):
        N2_scans = N2_testrun()

    def _test_plots(self):

        # inline plotting multiple figures
        [i._test_plot_scanrates() for i in N2_scans]
        [i._test_plot_Cdl() for i in N2_scans]
