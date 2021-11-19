"""

this module analyzes the ORR experiments, including Disk, N2 background subtraction and possible Ring electrode selectivity data

"""

## std lib
from typing import NamedTuple, Tuple, Dict, List
from collections import namedtuple
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local
import elchempy

from elchempy.dataloaders.fetcher import ElChemData
from elchempy.experiments.Logic.selection import DataSelection

from elchempy.experiments.ORR.RRDE.ring_helpers import (
    find_file_from_secondary_electrode,
)

from elchempy.experiments.N2.background.background_scan import N2_Background

# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate
from elchempy.experiments.ORR.RDE.plotting import RDE_Plotter


## 3rd party
import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

# constants
from elchempy.constants import EvRHE, Segment_lbl

#%%
def _dev_test_ORR_analysis():
    """function for testing"""
    from elchempy.experiments._dev_datafiles._dev_fetcher import get_files

    filepath = get_files("O2_ORR_DW28_0.5MH2SO4_PDX_Ch1_disk")[0]
    N2BG_file = get_files("N2_20cls_300_100_10_DW28_251")[0]
    ORR = ORR_Analysis(filepath, N2_background_scan_file=N2BG_file)


class ORR_Analysis(ElChemData, DataSelection, RDE_Plotter):
    """
    Inherits from ElChemData,
    performs the steps for the ORR analysis on data of a file

    Difficulties:
        Information from other files is required for processing,
        there should be a N2 background scan in another file in same folder

        there can be a RRDE measurement with Disk and Ring files, which
        need to be merged on Elapsed Time for the calculations

    Steps:
        N2 background file: look in folder or get file from argument (from index)

    """

    rpm_ref_list_mapper = {
        5: [0, 200, 400, 900, 1500],
        3: [0, 900, 1500],
        2: [1500, 1500],
        4: [0, 200, 400, 900],
        6: [0, 200, 400, 900, 1500, 2500],
    }

    def __init__(
        self,
        filepath: [Path, str],
        N2_background_scan_file=None,
        N2_background_scan_data=None,
        ORR_ring_file=None,
        auto_detect_ring_file=True,
        **kwargs
    ):
        # self.filepath = Path(filepath, **kwargs)
        # self.kwargs = kwargs
        # self.data = None
        super().__init__(filepath, **kwargs)

        self.N2_BG_file = N2_background_scan_file
        self.N2_BG_scan, self.N2_BG_message = self.get_N2_BG_scan(
            N2_background_scan_file=N2_background_scan_file,
            N2_background_scan_data=N2_background_scan_data,
        )

        # self.disk_ecd = ElChemData(filepath, **kwargs)
        if auto_detect_ring_file:
            self.ring_file = find_file_from_secondary_electrode(self.filepath)

    def get_rpm_list(
        self, disk_data, segment_label=Segment_lbl, rpm_mapper=None
    ) -> List:
        """guesses the RPM values from a reference mapper"""
        longsegs = [n for n, gr in disk_data.groupby(segment_label) if len(gr) > 1000]

        if not rpm_mapper:
            rpm_mapper = self.rpm_ref_list_mapper
        rpm_list = rpm_mapper.get(len(longsegs), [])
        logger.debug("Jkin calculation rpm list used: {0}".format(rpm_list))
        return rpm_list

    def get_N2_BG_scan(
        self, N2_background_scan_file=None, N2_background_scan_data=None
    ):
        message = ""
        if isinstance(N2_background_scan_data, pd.DataFrame):
            N2_BG_scan = N2_background_scan_data

            message = "side loaded"
            if N2_background_scan_data.empty:
                message += " frame empty"

            # N2_BG_file = N2_background_scan_file
            if N2_background_scan_file:
                message += ", file given"
            else:
                message += ", no file given"

        else:
            if N2_background_scan_file:
                # N2_BG_file = N2_background_scan_file
                # Load disk data from given filepath
                # TODO take class instance or only frame
                N2_BG = N2_Background(N2_background_scan_file)

                N2_BG_scan = N2_BG.N2_BG_scan
                message = "loaded from file"
                if N2_BG_scan.empty:
                    message += ", empty N2 BG scan"

            else:
                N2_BG_scan = pd.DataFrame()
                message = "no file and frame given"
        return N2_BG_scan, message
