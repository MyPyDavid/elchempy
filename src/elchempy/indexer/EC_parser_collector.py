""" Collects the files of index files of a folder"""

from pathlib import Path
from collections import Counter
from functools import wraps

import datetime
from typing import Tuple, List, Dict, Union, Collection

import re
import copy

import logging

logger = logging.getLogger(__name__)

from elchempy.indexer.filepath_parser import FilePathParser
from elchempy.indexer.extra_EC_info import loading_ref, WE_surface_area_cm2
from elchempy.indexer.helpers import (
    find_path_of_common_folder,
    relative_parent_paths_to_common_folder,
)

from elchempy.experiments.dataloaders.files_func_collector import run_func_on_files

from elchempy.indexer.EC_filepath_parser import ElChemPathParser
from elchempy.experiments.dataloaders.fetcher import ElChemData


### for Developing
from elchempy.config import LOCAL_FILES

### 3rd Party imports

import pandas as pd


class ElChemPathParserCollection:
    def __init__(self, files, multi_run=False):
        self._files = files
        self._multi_run = multi_run

        self.ecpps = run_func_on_files(
            ElChemPathParser, self._files, multi_run=self._multi_run
        )

        self.index = pd.DataFrame()
        self.index = pd.concat(
            [pd.DataFrame(val.EC_INDEX_entry).T for k, val in self.ecpps.items()],
            ignore_index=False,
        )
        self.common_folder = find_path_of_common_folder(self._files)
        self.rel_folders = relative_parent_paths_to_common_folder(
            self.index.index, self.common_folder
        )

        self.index = self.index.assign(
            **{
                "common_folder": self.common_folder,
                "relative_folders": self.rel_folders,
            }
        )

    def add_methods(self):
        """for persistence and loading of this 'collection' in a database of pkl file eg"""


#%%

if __name__ == "__main__":
    ecppcol = ElChemPathParserCollection(LOCAL_FILES, multi_run=True)
    self = ecppcol
