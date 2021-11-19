"""

this module prepares a background scan from the N2 data for an ORR measurement

"""

## std lib
from typing import List, Tuple

# from collections import namedtuple
# from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local
# import elchempy

# from elchempy.dataloaders.fetcher import ElChemData
# from elchempy.experiments.N2.background_scan import contains_background_scan, get_N2_background_data

# from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

## 3rd party
import numpy as np
import pandas as pd

# from scipy.stats import linregress, zscore

## constants
# from elchempy.constants import  EvRHE
#%%


def select_N2_background_scan_from_data(
    N2_CVs, segment_key="Segment #", verbose=False
) -> (pd.DataFrame, str):
    """finds the data segments for the N2 BG and selects the data"""
    if not isinstance(N2_CVs, pd.DataFrame):
        # optional warning for TypeError
        raise TypeError("N2_CVs data is not of type DataFrame")
        # return None,

    BG_segment_options, finder_message = find_valid_background_scan_segments(N2_CVs)
    if verbose:
        logger.error(finder_message)

    if not BG_segment_options:
        BG_scan = pd.DataFrame()
    else:
        BG_segment = BG_segment_options[-1]
        BG_scan = N2_CVs.loc[N2_CVs[segment_key] == BG_segment]

    return BG_scan, finder_message


def find_valid_background_scan_segments(
    N2_CVs, maximum_scanrate=0.011, scan_length=1950, segment_key="Segment #"
) -> List[int]:
    """checks if the data possibly contains a N2 background scan"""

    # finder_warning = None

    if not isinstance(N2_CVs, pd.DataFrame):
        # optional warning for TypeError
        return [], "wrong type"

    if N2_CVs.empty:
        # optional warning for empty data
        return [], "data is empty"

    sr_min = np.min([i for i in N2_CVs.scanrate.unique() if i > 0])

    if not (0 < sr_min <= maximum_scanrate):  # check presence of slow scanrates
        # logger.debug()
        # optional warning for missing data
        return (
            [],
            f"Minimum scanrate in data is larger than max {maximum_scanrate} or 0",
        )

    sr_grp_min = N2_CVs.loc[N2_CVs.scanrate == sr_min]
    # if len(sr_grp_min) < scan_length:
    #     logger.debug('Data selected at minimum scan rate is too short')
    #     return False
    BG_segment_options = [
        n
        for n, gr in sr_grp_min.groupby(segment_key)
        if (len(gr) < scan_length * 1.2 and len(gr) > scan_length * 0.8)
    ]
    if len(sr_grp_min) < scan_length:
        # logger.debug()
        return (
            [],
            f"Data selected at minimum scan rate {sr_min} is too short {len(sr_grp_min)} for scan {scan_length}",
        )

    if not BG_segment_options:
        # warning
        # logger.debug()
        # optional warning for having incomplete data
        return [], "No segments at minimum scan rate are valid"

    return BG_segment_options, 'valid segments {", ".join(map(str,BG_segment_options))}'


def validate_background_scan():
    pass


def _optional_manipulate_current_if_scan_is_missing(N2_CVs):
    """last resort option in order to retrieve a valid BG scan from data"""

    sr_min = np.min([i for i in N2_CVs.scanrate.unique() if i > 0])

    if sr_min > 0.01:
        N2_factor = sr_min / 0.01
        N2_CVs = N2_CVs.assign(
            **{"j_normalized_to_10mVs": N2_CVs.loc[:, "j A/cm2"] / N2_factor}
        )
        logger.warning(f"N2 scans minimun scanrate {sr_min} is larger than 10 mV/s")
    #                N2_scan.loc[:,'j A/cm2'] = N2_scan['j A/cm2']/N2_factor
    #                pd.DataFrame([ScanRates.min()])
    return N2_CVs
