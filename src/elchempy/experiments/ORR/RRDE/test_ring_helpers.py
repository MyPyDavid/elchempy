"""
testers for helpers of RRDE experiments.
file finder for ring
"""


# std lib

import logging

logger = logging.getLogger(__name__)

import unittest

## local
from elchempy.experiments.ORR.RRDE.ring_helpers import (
    find_file_from_secondary_electrode,
    count_matching_consecutive_strings,
)

## for developing and testing
from elchempy.experiments._dev_datafiles._dev_fetcher import get_files

## constants
# from elchempy.constants import EvRHE

### 3rd party
# import numpy as np
# import pandas as pd


#%%


class Test_ring_helpers(unittest.TestCase):
    """contains developer functions for testing"""

    def setUp(self):
        files = get_files("O2_ORR_DW28_0.5MH2SO4_PDX_Ch1_disk")
        self.disk_filepath = files[0]
        self.test_ring_filename = "O2_ORR_DW28_0.5MH2SO4_PDX_Ch2_ring.par"

    def test_find_file_from_secondary_electrode(self):
        """test the finder of secondary electrode file"""
        ring_file = find_file_from_secondary_electrode(self.disk_filepath)
        self.ring_file_found = ring_file
        self.assertEqual(self.ring_file_found.name, self.test_ring_filename)

    def test_count_matching_strings_filename(self):

        abc_count = count_matching_consecutive_strings(
            self.disk_filepath.name, self.test_ring_filename
        )
        self.assertEqual(abc_count, 28)

    def test_count_matching_strings_abc(self):

        abc_count = count_matching_consecutive_strings("abc1fg", "abc2de")
        self.assertEqual(abc_count, 3)


if __name__ == "__main__":
    unittest.main()
