"""
Created on Mon Jul 12 22:02:45 2021

@author: DW
"""

import sys
from pathlib import Path

import pytest
import unittest


import elchempy
from elchempy.experiments.EC_DataLoader.CreateCV import create_CVs




class TestCreateCV(unittest.TestCase):

    def setUp(self):

        N2_files = (i for i in test_data_files if i.name.startswith('N2'))
        N2_file  = next(N2_files)
        self.N2_file = N2_file


    def test_N2(self):

        try:
            N2_data = create_CVs(self.N2_file, kwargs)
        except Exception as e:
            pass




def _dev():
    dr= DataReader()


PAR_file_parser