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
logger.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

from elchempy.indexer.filepath_parser import FilePathParser
from elchempy.indexer.extra_EC_info import loading_ref, WE_surface_area_cm2
from elchempy.indexer.helpers import (
    find_path_of_common_folder,
    relative_parent_paths_to_common_folder,
)

from elchempy.experiments.dataloaders.files_func_collector import run_func_on_files

# from elchempy.experiments.dataloaders.fetcher import ElChemData

from elchempy.indexer.EC_path_parser import ElChemPathParser
from elchempy.experiments.dataloaders.fetcher import ElChemData


### for Developing
from elchempy.config import LOCAL_FILES

### 3rd Party imports

import pandas as pd

#%%
class ElChemPathParserCollection:
    """
    Collects the parsers instances for a list of files.
    Can include metadata from file instrospection or only from parsing the filename

    """

    def __init__(self, files, include_metadata=False, multi_run=False):
        self._files = files
        self._multi_run = multi_run
        self._include_metadata = include_metadata

        self.ecpps = run_func_on_files(
            ElChemPathParser, self._files, multi_run=self._multi_run
        )

        self.ecds = None
        if self._include_metadata:
            self.ecds = run_func_on_files(
                ElChemData,
                self._files,
                multi_run=self._multi_run,
                metadata_only=self._include_metadata,
            )

        self.index = self.merge_and_make_index_from_collections()
        # self.index = pd.concat(
        #     [pd.DataFrame(val.EC_INDEX_entry).T for k, val in self.ecpps.items()],
        #     ignore_index=False,
        # )
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

    def merge_and_make_index_from_collections(self):
        entries = []
        for k, val in self.ecpps.items():
            ecpp_entry = pd.DataFrame(val.EC_INDEX_entry).T
            rename_columns = {
                i: f"exp_{i}"
                for i in ecpp_entry.columns
                if not (i.startswith("PAR") or i.startswith("File"))
            }
            ecpp_entry = ecpp_entry.rename(columns=rename_columns)
            if self.ecds:
                ecd = self.ecds.get(k)
                actions_names = ", ".join(ecd.actions.Name.unique())
                ecpp_entry = ecpp_entry.assign(**{"actions_names": actions_names})
                # val.EC_INDEX_entry.get(k).update()
                metadata = ecd.DR.metadata
                ecpp_entry = pd.merge(
                    ecpp_entry, metadata, left_index=True, right_index=True
                )
            entries.append(ecpp_entry)
        if entries:
            index = pd.concat(entries, ignore_index=False)
        else:
            index = pd.DataFrame()
        return index

    def add_methods(self):
        """for persistence and loading of this 'collection' in a database of pkl file eg"""

    def store_index(self):
        pass

    def load_index(self):
        pass


#%%

if __name__ == "__main__":
    ecppcol = ElChemPathParserCollection(
        LOCAL_FILES, multi_run=True, include_metadata=True
    )
    self = ecppcol
