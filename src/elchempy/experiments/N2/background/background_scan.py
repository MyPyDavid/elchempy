"""

this module handles the N2 background scan data for the ORR experiments

"""

## std lib
# from typing import NamedTuple, Tuple, Dict
# from collections import namedtuple
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local
# import elchempy

from elchempy.dataloaders.fetcher import ElChemData

# from elchempy.experiments.ORR.helpers import get_file_from_secondary_electrode

# from elchempy.experiments.N2.analysis import N2_Analysis
# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

from elchempy.experiments.Logic.selection import DataSelection
from elchempy.experiments.N2.background.helpers import (
    select_N2_background_scan_from_data,
)

from elchempy.experiments.N2.plotting import N2_Plotter

## 3rd party
# import numpy as np
import pandas as pd

# from scipy.stats import linregress, zscore

## constants
from elchempy.constants import EvRHE

#%%
class N2_Background(ElChemData, DataSelection, N2_Plotter):
    """loads N2 background data from file"""

    def __init__(self, filepath: [Path, str], check_for_excluded_files=True, **kwargs):

        # self.kwargs = kwargs
        # self.data = None
        self._excluded_filepath = None
        if check_for_excluded_files and filepath:
            filepath_from_exclustion_list = N2_BG_exclusion_list(filepath)

            if filepath_from_exclustion_list != filepath:
                self._excluded_filepath = filepath
                filepath = filepath_from_exclustion_list

        super().__init__(filepath, **kwargs)

        (
            self.N2_BG_scan,
            self._N2_BG_finder_message,
        ) = select_N2_background_scan_from_data(self.data_selection)
        if self.N2_BG_scan.empty:
            logger.error(
                f"N2_BG_scan empty for file {self}: {self._N2_BG_finder_message}"
            )
        # get_BG_scan_and_file(N2_background_file=self.filepath, N2_background_data = )


def N2_BG_side_loader(N2_background_file=None, N2_background_data=pd.DataFrame()):
    """performs the necesssary steps for loading the N2 background"""
    # Check for N2 background scan data and file
    N2_BG_file = N2_background_file
    N2_BG_data = N2_background_data

    if not N2_BG_data.empty and N2_BG_file:
        # both data and file are given
        # calculations can start everything is there
        # TODO start calculations
        N2_BG_scan = select_N2_background_scan_from_data(N2_BG_data)
    elif not N2_BG_data.empty and not N2_BG_file:
        # calculations can start but filename missing
        # TODO ask or search for filename and start calculations
        pass
    elif N2_BG_data.empty and N2_BG_file:
        # background data needs to be loaded from file
        # TODO load from file
        # N2_BG_ecd = N2_Analysis(N2_BG_file, calculate_cdl=False)
        # TODO check and select BG data
        pass
        # N2_BG_data = N2_BG_ecd.N2_BG

    elif N2_BG_data.empty and not N2_BG_file:
        # missing data and filename, skip calculations or try and search for file
        pass
    return N2_BG_scan


def N2_BG_exclusion_list(N2_BG_file):
    """takes a file and filters/replaces it depending on occurence in exclusion list mapper"""
    if not (isinstance(N2_BG_file, str) or isinstance(N2_BG_file, Path)):
        return None
    else:
        N2_BG_file = Path(N2_BG_file)
    # _topdir = "/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Organized_Data/VERSASTAT"
    N2_BG_exclusion_replace_mapper = {
        "25.01.2019_0.1MH2SO4_cell2/N2_20cls_300_100_10_JOS4_256.par": "24.01.2019_0.1MH2SO4_24490_cell3/N2_20cls_300_100_10_JOS4_270.par",
        "25.01.2019_0.1MH2SO4_cell2/N2_20cls_300_100_10_JOS4_256_fail.par": "24.01.2019_0.1MH2SO4_24490_cell3/N2_20cls_300_100_10_JOS4_270.par",
        # '2019-01-Jan/25.01.2019_0.1MH2SO4_cell3/N2_scans_v30/N2_20cls_300_100_10_JOS4_postAST-LC_256_BG.pkl',
    }
    N2_BG_file_key = f"{N2_BG_file.parent.name}/{N2_BG_file.name}"

    if N2_BG_file_key not in N2_BG_exclusion_replace_mapper.keys():
        return N2_BG_file

    N2_BG_file_mapper_val = N2_BG_exclusion_replace_mapper.get(N2_BG_file_key, None)
    if not N2_BG_file_mapper_val:
        return N2_BG_file
        # print("JOS4 look for alternative N2_BG files")
    p, name = N2_BG_file_mapper_val.split("/")
    N2_BG_file_replaced = N2_BG_file.parent.parent.joinpath(p) / name

    return N2_BG_file_replaced
