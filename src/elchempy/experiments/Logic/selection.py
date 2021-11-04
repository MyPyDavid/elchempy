"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

# std lib

# from typing import NamedTuple, Tuple, Dict
# from collections import namedtuple
# from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local

### 3rd party
import pandas as pd

## constants


#%%
class DataSelection:
    """Mixin Class to provide selection methods on the data for each experiments"""

    def __init__(self):

        self._DS_expname = self.__class__.__qualname__
        self._DS_DR_parsername = self.DR.parser.__class__.__qualname__
        logger.info(
            f"Called DataSelection init from {self._DS_expname} with parser {self._DS_DR_parsername}"
        )

        # set empty DF as default value for attr
        self.data_selection = pd.DataFrame()

        # preconditional checks

        if not self._DS_DR_parsername:
            logger.warning(
                "Class does not have a DR.parser attribute, file contains data?"
            )
            return
        if not isinstance(self.data, pd.DataFrame):
            logger.warning(
                "Class instance does not have data,\n{self.filepath} file contains data?"
            )
            return
        else:
            if self.data.empty:
                logger.warning(
                    "Class instance data is empty,\n{self.filepath} file contains data?"
                )

        # Parent class calling
        self._qname = self.__class__.__qualname__

        # start selection based of self class qualname
        # TODO can develop into switch case statements

        if self._qname.startswith("N2"):
            data_selection = self.N2_selection(self.data, self._DS_DR_parsername)
        elif self._qname.startswith("ORR"):
            data_selection = pd.DataFrame()
            pass
        else:
            data_selection = pd.DataFrame()
            raise NotImplementedError(
                f"Data selection method for {self._qname} not implemented."
            )
        self.data_selection = data_selection

    def N2_selection(self, data: pd.DataFrame, parsername: str) -> pd.DataFrame:
        # FIXME Select only CV types from Data segment
        # Select the data for N2 Cyclic Voltammograms

        try:
            if "VersaStudio" in parsername:
                N2_CVs = data.loc[data.ActionId == 38]
                N2_CVs = N2_CVs.dropna(subset=["scanrate"]).loc[
                    N2_CVs.scanrate_calc != 0
                ]
            else:
                raise NotImplementedError(f"{parsername} ")

        except Exception as ex:
            logger.error(f"{self} select data error\m{ex}")
            N2_CVs = pd.DataFrame()
        else:
            if N2_CVs.empty:
                logger.warning("select_data is empty, file does not contain any N2 CVs")

        return N2_CVs
