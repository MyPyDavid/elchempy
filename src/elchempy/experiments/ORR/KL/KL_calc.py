import sys
from pathlib import Path
from collections import namedtuple
from datetime import datetime
import numpy as np

from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines

import os
import multiprocessing
from functools import partial
from itertools import repeat
import pandas as pd

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper.ExtraInfo import EC_Properties


if __name__ == "__main__":
    pass

import logging

logger = logging.getLogger(__name__)


class KL_operations:
    def __init__(self, KL_data, ORR_dest_dir_file, electrode_properties={}):
        self.KL_data = KL_data
        self.ORR_dest_dir_file = ORR_dest_dir_file
        self.electrode_properties = electrode_properties

    def prepare_data(self, KL_data):
        _KL_meta = (
            KL_data[[i for i in KL_data.columns if KL_data[i].nunique() == 1]]
            .iloc[0]
            .to_dict()
        )
        KL_coeff = KL_coefficients()

        KLcoeffpH = KL_coeff.loc[
            (KL_coeff.pH == _KL_meta.get("pH", 99))
            & (KL_coeff.Electrolyte == _KL_meta.get("Electrolyte", "Unknown")),
            :,
        ]
        F, D, nu, C = (
            96485,
            KLcoeffpH.D0.values[0],
            KLcoeffpH.kVis.values[0],
            KLcoeffpH.C0.values[0],
        )
        # FIXME
        Area = self.electrode_properties.get("area_cm2", 1)
        # TODO WE_SA_collection_eff('PINE')
