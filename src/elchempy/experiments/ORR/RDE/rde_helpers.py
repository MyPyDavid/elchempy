"""
RDE helper functions
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


def apply_BG_subtraction_to_data(disk_data=None, N2_BG_scan=None):
    """applies background subtraction on the data"""

    if not (
        isinstance(disk_data, pd.DataFrame) and isinstance(disk_data, pd.DataFrame)
    ):
        return

    disk_data
