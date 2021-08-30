"""
Fetches data from a file and constructs the electrochemical data columns
"""

from collections import namedtuple
from dataclasses import dataclass

from pathlib import Path

import pandas as pd
import numpy as np

import elchempy
from elchempy.experiments.dataloaders.reader import DataReader

# from .converters import get_current_density, get_potential_vs_RE, get_RPM_from_DAC_V
from elchempy.experiments.dataloaders.assigners import assign_all_EC_relevant_columns


#%%
def _dev():
    from elchempy.experiments._dev_datafiles._dev_fetcher import (
        get_files,
        _dev_test_read,
    )

    N2fls = get_files("N2")

    try:
        ecd = ElChemData(N2fls[0].with_suffix(".noname"))
        ba = ElChemData("a")
    except FileNotFoundError:
        pass
    except Exception as ex:
        print("unexpected error", ex)
        raise ex from ex

    self = ecd


def reference_tables():
    # '.par' : {
    ActionId_to_Type_Reference = {
        "Cyclic Voltammetry (Multiple Cycles)": 38,
        "Chronoamperometry": 3,
        "Unknown": 0,
        "Potentiostatic EIS": 21,
        "Cyclic Voltammetry": 14,
    }
    return ActionId_to_Type_Reference


def _get_ECdata_from_file(arg):
    """
    Parameters
    ----------
    arg : [Path, str, pd.DataFrame, ]
        DESCRIPTION.

    Returns
    -------
    data : pd.DataFrame
        DESCRIPTION.

    """
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


@dataclass
class ElChemData:
    """
    Parses the data from a .PAR file into relevant electrochemical data.

    #TODO idea; turn this Class into subclass of a DataFrame and include functionalities
    This class contains all functions which add several collumns to an
    existing DataFrame with information that is important for the analsysis
    of electrochemimcal experiments.
    """

    filepath: [Path, str]

    def __post_init__(self, **kwargs):

        if not (isinstance(self.filepath, Path) or isinstance(self.filepath, str)):
            raise TypeError(
                f"{self.__class__.__qualname__} filepath arg is not Path nor string, {type(self.filepath)}"
            )

        # self.filepath = filepath
        self.DR = DataReader(self.filepath)
        self.raw_actions = self.DR.actions
        self.raw_data = self.DR.data
        self.data = self.raw_data.copy(deep=True)
        self.actions = self.raw_actions.copy(deep=True)

        if self.DR:
            self.data, self.new_EC_colums = assign_all_EC_relevant_columns(
                self.data, self.actions, **kwargs
            )
        self.methods = []

    def add_analysis_method(self, results_from_method):
        methodname = results_from_method.__class__.__qualname__
        # if not hasattr(self, methodname):
        setattr(self, methodname, results_from_method)
        self.methods.append(methodname)
        # else:

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        _name = Path(self.filepath).name
        _txt = f'actions={len(self.actions)}, data={len(self.data)}, methods={", ".join(map(str, self.methods))}'
        return f"{self.__class__.__qualname__}: {_name}, {_txt}"
        # return f"{self.__class__.__qualname__}('{str(self.filepath)}')"

    def __str__(self):
        return str(Path(self.filepath))
