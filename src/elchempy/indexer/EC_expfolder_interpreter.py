""" Interpretes the files of index files of a folder"""

from pathlib import Path
from collections import Counter
from functools import wraps
from itertools import groupby

import datetime
from typing import Tuple, List, Dict, Union, Collection

import re
import copy

import logging

logger = logging.getLogger(__name__)

from elchempy.indexer.filepath_parser import FilePathParser
from elchempy.indexer.extra_EC_info import loading_ref, WE_surface_area_cm2

from elchempy.experiments.dataloaders.files_func_collector import run_func_on_files

from elchempy.indexer.EC_filepath_parser import ElChemPathParser
from elchempy.experiments.dataloaders.fetcher import ElChemData


### for Developing
from elchempy.config import LOCAL_FILES

### 3rd Party imports

import datefinder

import dateutil
from dateutil.parser import ParserError


#%%
def _dev():

    expm = ExpFolder(files)
    self = expm
    n2 = self.N2data[
        "/mnt/DATA/APPS_SOFT/REPOS/elchempy/src/elchempy/experiments/_dev_datafiles/06.03.2018_DW28_HPRR_0.1MHClO4_RRDE22960/N2_20cls_300_100_10_DW28_298.par"
    ]
    n2._test_plot_Cdl()
    [i._test_plot_scanrates() for i in self.N2data.values()]


class ExpFolder:
    """
    Takes in the files and groups by Experimental Folder and retrieves the relevant data.

    Typical types of electrochemical experiments:
        OCP: open circuit potential
            used for determining the RE_vs_RHE potential

        N2, N2_act: activation scans in N2 saturated electrolyte
            contains cycles of CVs at certain scanrates

        O2 ORR: oxygen reduction reaction scans in O2 saturated electrolyte
            contains cycles at slow scanrate (5, 10 mV/s) at several rotation speeds

        EIS: electrochemical impedance spectroscopy
            contains EIS spectra at certain potentials and/or rotation speeds


    Requirements for
        Any:
            a potential RE_vs_RHE value:
            - The potential RE_vs_RHE should be checked per experiment in the same cell.
            - Can be determined from the OCP measurement
            - If OCP is missing then there should be a value appended to the individual filename (eg. *_123.par)
            - Else the value can not be found and set to 0
            - Optionally, the value can be derived from other files in the same folder.

        O2 ORR:
            requires a N2 background scan at same scanrate (10 mV/s)

    """

    def __init__(self, files: Collection, index=None, multi_run=False):
        self._files = files
        self._index = index
        self._multi_run = multi_run

        self.ecpprs = self.run(ElChemPathParser)
        self.ecdata = self.run(ElChemData)

        # self.groupby = groupby(self.ecpprs.keys(), key=lambda x: Path(x).parent.name)
        self.groups, self.grpkeys = get_groups_keys(
            self.ecpprs.keys(), lambda x: Path(x).parent
        )

    def run(self, func):
        run_result = run_func_on_files(func, self._files, multi_run=self._multi_run)
        return run_result


def get_groups_keys(data, keyfunc):
    # self.groupby = groupby(self.ecpprs.keys(), key=lambda x: Path(x).parent.name)
    groups = []
    uniquekeys = []
    data = sorted(data, key=keyfunc)
    for k, g in groupby(data, keyfunc):
        groups.append(list(g))  # Store group iterator as a list
        uniquekeys.append(k)
    return groups, uniquekeys


#%%

if __name__ == "__main__":
    expfldr = ExpFolder(LOCAL_FILES, multi_run=True)
    self = expfldr
