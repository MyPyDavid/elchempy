"""

this module provides a mixin base class for data selection

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
            data_selection = self.ORR_disk_selection(self.data, self._DS_DR_parsername)
        else:
            data_selection = pd.DataFrame()
            raise NotImplementedError(
                f"Data selection method for {self._qname} not implemented."
            )
        self.data_selection = data_selection

    def N2_selection(
        self, data: pd.DataFrame, parsername: str, scanrate_min=0
    ) -> pd.DataFrame:
        # FIXME Select only CV types from Data segment
        # Select the data for N2 Cyclic Voltammograms

        try:
            if "VersaStudio" in parsername:
                N2_CVs = data.loc[data.ActionId == 38]

            else:
                raise NotImplementedError(
                    f"{parsername} selection for N2 not implemented "
                )

        except Exception as ex:
            logger.error(f"{self} select data error\m{ex}")
            N2_CVs = pd.DataFrame()
        finally:
            if N2_CVs.empty:
                logger.warning("select_data is empty, file does not contain any N2 CVs")
                return N2_CVs

        N2_CVs = N2_CVs.dropna(subset=["scanrate"]).loc[N2_CVs.scanrate != 0]

        return N2_CVs

    def ORR_disk_selection(
        self, data: pd.DataFrame, parsername: str, scanrate_min=0, scanrate_max=0.015
    ) -> pd.DataFrame:
        # FIXME Select only CV types from Data segment
        # Select the data for N2 Cyclic Voltammograms

        try:
            if "VersaStudio" in parsername:
                ORR_CVs = data.loc[data.ActionId == 38]
            else:
                raise NotImplementedError(
                    f"{parsername} selection for ORR disk not implemented "
                )

        except Exception as ex:
            logger.error(f"{self} select data error\m{ex}")
            ORR_CVs = pd.DataFrame()
        finally:
            if ORR_CVs.empty:
                logger.warning(
                    "select_data is empty, file does not contain any ORR_CVs"
                )
                return ORR_CVs

        ORR_CVs = ORR_CVs.dropna(subset=["scanrate"]).loc[
            (ORR_CVs.scanrate > scanrate_min) & (ORR_CVs.scanrate < scanrate_max)
        ]

        return ORR_CVs
