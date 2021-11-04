"""
Fetches data from a file and constructs the electrochemical data columns
"""

from collections import namedtuple
from dataclasses import dataclass, InitVar

from pathlib import Path
from typing import Union

import pandas as pd
import numpy as np

import elchempy
from elchempy.dataloaders.reader import DataReader

# from .converters import get_current_density, get_potential_vs_RE, get_RPM_from_DAC_V
from elchempy.dataloaders.assigners import assign_all_EC_relevant_columns


#%%


class ElChemData:
    """

    Parses the data from a .PAR file into relevant electrochemical data.

    Calls the DataReader on a filepath and unpacks the relevant components.
    #TODO idea; turn this Class into subclass of a DataFrame and include functionalities
    This class contains all functions which add several collumns to an
    existing DataFrame with information that is important for the analsysis
    of electrochemimcal experiments.
    """

    args_keys = ["_filepath"]
    kwargs_keys = ["_metadata_only", "_kwargs"]
    called_class_keys = ["DR"]
    data_keys = ["metadata", "actions", "data", "start_time", "methods"]
    new_data_keys = ["new_EC_columns"]

    __slots__ = args_keys + kwargs_keys + called_class_keys + data_keys + new_data_keys

    def __init__(
        self, filepath: Union[Path, str], metadata_only: bool = False, **kwargs
    ):
        self._filepath = filepath
        self._metadata_only = metadata_only
        self._kwargs = kwargs
        self.DR = None

        if not self._filepath:
            return

        if not (isinstance(self._filepath, Path) or isinstance(self._filepath, str)):
            raise TypeError(
                f"{self.__class__.__qualname__} filepath arg is not Path nor string, {type(self.filepath)}"
            )

        # Calls the DataReader on the filepath
        self.DR = DataReader(self._filepath, metadata_only=self._metadata_only)

        # Unpacking of the DataReader instance in components
        self.metadata = self.DR.metadata.copy()
        # self.raw_actions = self.DR.actions
        self.actions = self.DR.actions.copy(deep=True)

        # self.raw_data = self.DR.data
        self.data = self.DR.data.copy(deep=True)

        # self.raw_metadata = self.DR.metadata
        self.start_time = self.DR.start_time
        # self.metadata = self.raw_metadata.copy(deep=True)

        if self.DR:
            if not self.DR.data.empty:
                self.data, self.new_EC_columns = assign_all_EC_relevant_columns(
                    self.data, self.actions, **kwargs
                )
        self.methods = []

        # in case of multiple inheritance subclass
        super().__init__()

    @property
    def filepath(self) -> Union[Path, str]:
        return self._filepath

    @property
    def metadata_only(self) -> bool:
        return self._metadata_only

    def add_analysis_method(self, results_from_method):
        methodname = results_from_method.__class__.__qualname__
        # if not hasattr(self, methodname):
        setattr(self, methodname, results_from_method)
        self.methods.append(methodname)
        # else:

    def __bool__(self):
        return bool(self.DR)

    def __len__(self):
        if self._metadata_only:
            return len(self.metadata)
        else:
            return len(self.data)

    def __repr__(self):
        if not self._filepath:
            return f"{self.__class__.__qualname__}()"
        _name = Path(self._filepath).name
        _txt = f"metadata={len(self.metadata)}, actions={len(self.actions)}, data={len(self.data)}"
        _methods = f'methods={", ".join(map(str, self.methods))}'
        return f"{self.__class__.__qualname__}: {_name}, {_txt}"
        # return f"{self.__class__.__qualname__}('{str(self.filepath)}')"

    def __str__(self):
        if not self._filepath:
            return f"{self.__class__.__qualname__}()"
        return str(Path(self.filepath))
