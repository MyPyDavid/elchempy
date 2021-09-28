""" Collects the files of index files of a folder"""

from pathlib import Path
from collections import Counter
from functools import wraps

import datetime
from typing import Tuple, List, Dict, Union, Collection

# import re
# import copy

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

from elchempy.indexer.data import DATABASE

# from elchempy.dataloaders.files_func_collector import run_func_on_files
# from elchempy.experiments.dataloaders.fetcher import ElChemData

from elchempy.indexer.helpers import find_relevant_files_in_folder

from elchempy.indexer.creator import create_index


### for Developing
from elchempy.config import LOCAL_FILES, RAW_DATA_DIR, DEST_DATA_DIR

### 3rd Party imports

import pandas as pd
from pyarrow import ArrowInvalid

#%%
class ElChemIndex:
    """
    Collects the parsers instances for a list of files.
    Can include metadata from file instrospection or only from parsing the filename

    """

    supported_store_types = ["feather"]

    def __init__(
        self,
        folder,
        dest_dir=None,
        include_metadata=False,
        multi_run=False,
        store_type="feather",
        force_reload=False,
    ):
        self._folder = Path(folder)
        self._dest_dir = Path(dest_dir)
        self._multi_run = multi_run
        self._include_metadata = include_metadata
        self._store_type = store_type
        self._force_reload = force_reload

        self.files = find_relevant_files_in_folder(self._folder)
        # = files#[500::]
        self.store_file = self.get_store_file(
            store_type=self._store_type, dest_dir=self._dest_dir
        )

        loaded_index = self.load_index(self.store_file)
        if isinstance(loaded_index, pd.DataFrame) and not self._force_reload:
            index = loaded_index
            ecpps, ecds = None, None  # class instance objects are not reloaded
        else:
            index, ecpps, ecds = create_index(
                self.files,
                multi_run=self._multi_run,
                include_metadata=self._include_metadata,
            )
            self.store_index(index, self.store_file, overwrite=True)

        self.index = index
        self.ecpps = ecpps
        self.ecds = ecds

    def add_methods(self):
        """for persistence and loading of this 'collection' in a database of pkl file eg"""

    def get_store_file(self, store_type="", dest_dir=None, filename="index"):
        if not (store_type and dest_dir):
            return None

        if store_type not in self.supported_store_types:
            logger.warning("store type {store_type} is not supported")
            return None

        if dest_dir.is_file():
            dest_dir = dest_dir.parent

        daily_filepath = dest_dir.joinpath(f"{datetime.date.today()}_{filename}")

        if "feather" in store_type:
            store_file = daily_filepath.with_suffix(".feather")

        return store_file

    def store_index(self, index, store_file: Path = None, overwrite=False):
        if not store_file:
            logger.warning(f"No store file given: {store_file}")
            return None
        if not isinstance(index, pd.DataFrame):
            logger.warning(f"Index type is not pd.DataFrame: {type(index)}")
            return None
        if store_file.exists() and not overwrite:
            logger.warning(f"Index file exists and will not be overwritten.")
            return None

        index = index.reset_index()

        try:
            index.to_feather(store_file)
        except ArrowInvalid as exc:
            logger.error(f"error to_feather: {store_file}\n{exc}")

        logger.info(f"Index saved to : {store_file}")

    def load_index(self, store_file: Path = None):
        if not store_file.exists():
            logger.warning(f"Store file does not exist: {store_file}")
            return None
        try:
            index = pd.read_feather(store_file)
        except ArrowInvalid as exc:
            logger.error(f"error read_feather: {store_file}\n{exc}")
            index = None

        logger.info(f"Index loaded from : {store_file}")
        return index


#%%
def _dev_testing():
    files = self._files
    ElChemPathParser(files[159])


if __name__ == "__main__":

    folder = LOCAL_FILES[0].parent.parent
    # folder = '/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT'

    ecppcol = ElChemIndex(
        RAW_DATA_DIR,
        dest_dir=DEST_DATA_DIR,
        multi_run=True,
        include_metadata=True,
        force_reload=False,
    )
    self = ecppcol
    idx = self.index
    # print(self.index.ecpp_token_remainder.values)
    # print(self.index['ecpp_token_remainder'].dropna().unique())
