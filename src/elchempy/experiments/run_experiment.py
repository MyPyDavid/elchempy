"""
Created on Wed Jul  7 15:56:24 2021

@author: DW
"""

from typing import List, Tuple
from abc import ABC, abstractmethod
from collections.abc import (
    Collection,
    Mapping,
    Sequence,
    Generator,
    Iterable,
    MutableSequence,
    Hashable,
    Container,
    Callable,
)

from pathlib import Path

import pandas as pd

from elchempy.config import LOCAL_FILES
from elchempy.experiments.dataloaders.files_func_collector import run_func_on_files

from elchempy.indexer.EC_filepath_parser import ElchemPathParser
from elchempy.experiments.dataloaders.fetcher import ElChemData
from elchempy.experiments.N2.analyses import N2_Data

import logging

logger = logging.getLogger(__name__)


def _dev():

    expm = ExperimentManager(files)
    self = expm
    n2 = self.N2data[
        "/mnt/DATA/APPS_SOFT/REPOS/elchempy/src/elchempy/experiments/_dev_datafiles/06.03.2018_DW28_HPRR_0.1MHClO4_RRDE22960/N2_20cls_300_100_10_DW28_298.par"
    ]
    n2._test_plot_Cdl()
    [i._test_plot_scanrates() for i in self.N2data.values()]


class ExperimentManager:
    def __init__(self, files, multi_run=False):
        self._files = files
        self._multi_run = multi_run

        self.ecpfls = run_func_on_files(
            ElchemPathParser, self._files, multi_run=self._multi_run
        )

        # self.ecdata = run_func_on_files(ElChemData, self._files, multi_run= self._multi_run)

        self.N2data = run_func_on_files(N2_Data, self._files, multi_run=self._multi_run)

    def call_exp_file_interpreter(self):

        pass

    def call_exp_file_introspectors(self):

        ecdata = get_elchem_parsers_from_files(
            ElChemData, self._files, multi_run=self._multi_run
        )
        self.ecdata = ecdata

    def exp_file_explorers(self):
        pass

    def exp_developer(self):
        pass


class BaseRunner(ABC):
    def get_data_from_collection_or_other(self, arg):
        """
        Parameters
        ----------
        arg : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        if isinstance(arg, Callable):
            raise ValueError("arg can not be a Callable")

        if not isinstace(arg, Iterable) or isinstace(arg, str):
            return self.get_data_from_input(arg)

        if isinstace(arg, MutableSequence):
            self.start_loop_input_collection()
        # else:
        #     if isinstace(arg, MutableSequence):
        # self.start_loop_input_collection(arg)
        # for i in arg:

        else:
            self.get_data_from_input(arg)

    def start_loop_input_collection(arg: MutableSequence) -> pd.DataFrame:
        for i in arg:
            self.get_data_from_input(i)

    @staticmethod
    def check_data_input(self):
        pass

    # @abstractmethod
    def check_input_slice():
        pass

    # @abstractmethod
    def run_parallel(self):
        pass

    # @abstractmethod
    def run_serial(self):
        pass


if __name__ == "__main__":

    expm = ExperimentManager(LOCAL_FILES, multi_run=True)
    [i._test_plot_scanrates() for i in expm.N2data.values()]
