""" collects all function calls on list of files with options for multi or single core processing"""

import os
from pathlib import Path
from typing import List, Collection, Dict

# from functools import partial
from itertools import repeat

import concurrent.futures
from multiprocessing import Pool

import logging

logger = logging.getLogger(__name__)

# Local imports
# from elchempy.indexer.filename_parser import FilePathParser
# from elchempy.indexer.EC_filepath_parser import ElchemPathParser

# 3rd party
import pandas as pd

#%%


def wrapper_func(*args):
    # print(f'args: {arg}\nkwargs: {kwargs}')
    # **kwargs
    func, file, kwargs = args
    try:
        result = func(file, **kwargs)
        return result
    except Exception as exc:
        logger.info(f"Error in multiprocess wrapper {exc}")
        # result = None


def run_func_on_files(func, files, multi_run=False, **kwargs) -> Dict:
    collection = []

    if multi_run:
        collection = make_collection_multi(func, files, **kwargs)
    else:
        collection = make_collection_serial(func, files, **kwargs)
    # breakpoint()

    collect_dict = {str(i): i for i in collection}
    if not collect_dict:
        logger.warning(f"Collection len={len(collection)}, collect_dict is empty.")

    try:
        # _test = str(ecpp_collection[0])
        return collect_dict
    except TypeError as ex:
        raise ex from ex
    except Exception as ex:
        raise ex from ex


def make_collection_multi(func: callable, files: Collection, **kwargs) -> List:

    collection = []
    with Pool(os.cpu_count() - 2) as pool:
        try:
            #                    results = pool.map(EC_classifier_multi_core.EC_PAR_file_check, self.par_files_run)
            collection = pool.starmap(
                wrapper_func, zip(repeat(func), files, repeat(kwargs))
            )
        except Exception as ex:
            #                    print('FileHelper module not found:',e)
            logger.error(f"make_collection_multi multiprocessing error: {ex}")
            raise ex from ex
            # results = pool.map(PAR_file_parser, self.par_files_run)
    return collection


def make_collection_serial(func: callable, files: Collection, **kwargs) -> List:

    ecpp_collection = []
    for file in files:
        try:
            logger.debug(f"{__name__} calling {func} on\n{file}.")
            ecpp = func(file, **kwargs)

            ecpp_collection.append(ecpp)
        except Exception as ex:
            _err = {"PAR_file": file, "error": ex, "kwargs": kwargs}
            logger.warning(
                f"{__name__} make_collection unexpected error for calling {func} on\n{file}.\n{ex}"
            )
            raise ex from ex

    return ecpp_collection
