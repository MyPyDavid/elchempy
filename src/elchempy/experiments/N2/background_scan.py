"""

this module prepares a background scan from the N2 data for an ORR measurement

"""

## std lib
from typing import NamedTuple, Tuple, Dict
from collections import namedtuple
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local
import elchempy

# from elchempy.dataloaders.fetcher import ElChemData
# from elchempy.experiments.N2.background_scan import contains_background_scan, get_N2_background_data

from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

## 3rd party
import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

## constants
EvRHE = "E_vs_RHE"
#%%


def select_background_scan_segments(
    N2_CVs, maximum_scanrate=0.011, scan_length=1950, segment_key="Segment #"
):
    """checks if the data possibly contains a N2 background scan"""

    if not isinstance(N2_CVs, pd.DataFrame):
        return False

    if N2_CVs.empty:
        return False

    sr_min = np.min([i for i in N2_CVs.scanrate.unique() if i > 0])

    if not (0 < sr_min <= maximum_scanrate):  # check presence of slow scanrates
        logger.debug("Minimum scanrate in data is larger than max or 0")
        return False

    sr_grp_min = N2_CVs.loc[N2_CVs.scanrate == sr_min]
    # if len(sr_grp_min) < scan_length:
    #     logger.debug('Data selected at minimum scan rate is too short')
    #     return False
    BG_segment_options = [
        n
        for n, gr in sr_grp_min.groupby(segment_key)
        if (len(gr) < scan_length * 1.2 and len(gr) > scan_length * 0.8)
    ]

    if not BG_segment_options:
        logger.debug("Data selected at minimum scan rate is too short")
        return False

    return BG_segment_options


def get_N2_background_data(N2_CVs, segment_key="Segment #"):

    BG_segment_options = select_background_scan_segments(N2_CVs)

    if not BG_segment_options:
        return None
    else:
        BG_segment = BG_segment_options[-1]
        BG_scan = N2_CVs.loc[N2_CVs[segment_key] == BG_segment]
        return BG_scan


def _optional_manipulate_current_if_scan_is_missing(N2_CVs):
    """last resort option in order to retrieve a valid BG scan from data"""

    sr_min = np.min([i for i in N2_CVs.scanrate.unique() if i > 0])

    if sr_min > 0.01:
        N2_factor = sr_min / 0.01
        N2_scan = preN2_scan.assign(
            **{"j_normalized_to_10mVs": preN2_scan.loc[:, "j A/cm2"] / N2_factor}
        )
        logger.warning(f"N2 scans minimun scanrate is larger than 10 mV/s")
    #                N2_scan.loc[:,'j A/cm2'] = N2_scan['j A/cm2']/N2_factor
    #                pd.DataFrame([ScanRates.min()])
