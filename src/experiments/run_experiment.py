"""
Created on Wed Jul  7 15:56:24 2021

@author: DW
"""

from abc import ABC


class BaseRunner(ABC):
    @abstractmethod
    def check_input_slice():
        pass

    @abstractmethod
    def run_parallel(self):
        pass

    @abstractmethod
    def run_serial(self):
        pass
