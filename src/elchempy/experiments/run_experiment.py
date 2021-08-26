"""
Created on Wed Jul  7 15:56:24 2021

@author: DW
"""

from typing import List, Tuple
from abc import ABC, abstractmethod
from collections.abc import Collection, Mapping, Sequence, Generator, Iterable, MutableSequence, Hashable , Container, Callable

from pathlib import Path

import pandas as pd

from elchempy.experiments.dataloader.fetcher import ElchemData

import logging
logger = logging.getLogger(__name__)

class ExperimentManager:

    def __init__(self,arg):
        pass

    def call_exp_file_interpreter(self):
        pass

    def exp_file_explorers(self):
        pass

    def exp_developer(self):



class BaseRunner(ABC):

    def get_data_from_collection_or_other(self, arg):
        '''
        Parameters
        ----------
        arg : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if isinstance(arg, Callable):
            raise ValueError('arg can not be a Callable')

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
    def get_data_from_input(arg):
        '''
        Parameters
        ----------
        arg : [Path, str, pd.DataFrame, ]
            DESCRIPTION.

        Returns
        -------
        data : pd.DataFrame
            DESCRIPTION.

        '''
        if not arg:
            return None

        data = None

        if isinstance(arg, Path) or isinstance(arg, str):
            try:
                data = ElchemData(filepath)
            except Exception as e:
                logger.warning(f"get_data_from_input failed for {arg}")

        elif isinstance(arg, pd.DataFrame):
            # TODO add possible double checks if hasattr(arg, column) ...
            self.check_data_input()

            data = arg

        return data

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
