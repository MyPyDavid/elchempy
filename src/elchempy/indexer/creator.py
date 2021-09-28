""" Creates the index of files from the specific parser objects """

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


# from elchempy.indexer.filepath_parser import FilePathParser
# from elchempy.indexer.extra_EC_info import loading_ref, WE_surface_area_cm2
from elchempy.dataloaders.files_func_collector import run_func_on_files
from elchempy.dataloaders.fetcher import ElChemData
from elchempy.indexer.EC_path_parser import ElChemPathParser
from elchempy.indexer.helpers import (
    find_path_of_common_folder,
    relative_parent_paths_to_common_folder,
    find_relevant_files_in_folder,
)

### 3rd Party imports

import pandas as pd

#%%


def create_index(files, include_metadata=False, multi_run=False):
    ecpps = None
    ecpps = run_func_on_files(ElChemPathParser, files, multi_run=multi_run)

    ecds = None
    if include_metadata:
        ecds = run_func_on_files(
            ElChemData,
            files,
            multi_run=multi_run,
            metadata_only=include_metadata,
        )

    index = pd.DataFrame()
    if not ecpps:
        return index, ecpps, ecds

    try:
        index = merge_and_make_index_from_collections(ecpps, ecds)
    except Exception as exc:
        logger.error(f"create index error, index will be empty : {exc}")

    index = add_extra_folder_columns_to_index(files, index)

    return index, ecpps, ecds


def merge_and_make_index_from_collections(ecpps, ecds, prefix="ecpp"):
    entries = []
    for k, val in ecpps.items():
        if not val:
            continue
        ecpp_entry = pd.DataFrame(val.EC_INDEX_entry).T
        rename_columns = {
            i: f"{prefix}_{i}"
            for i in ecpp_entry.columns
            if not (i.startswith("PAR") or i.startswith("File"))
        }
        ecpp_entry = ecpp_entry.rename(columns=rename_columns)
        if ecds:
            ecd = ecds.get(k)
            actions_names = ", ".join(map(str, ecd.actions.Name.unique()))
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


def add_extra_folder_columns_to_index(files, index, to_string=True):
    common_folder = find_path_of_common_folder(files)
    rel_folders = relative_parent_paths_to_common_folder(index.index, common_folder)
    if to_string:
        common_folder = str(common_folder)
        rel_folders = list(map(str, rel_folders))

    index = index.assign(
        **{
            "common_folder": common_folder,
            "relative_folders": rel_folders,
        }
    )
    return index
