"""Class and methods for parsing electrochemically relevant information for a filepath"""
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import Counter
from functools import wraps

import datetime
from typing import Tuple, List, Dict, Union

import re
import copy

import logging

logger = logging.getLogger(__name__)

from elchempy.indexer.filepath_parser import FilePathParser
from elchempy.indexer.EC_tokenizer import tokenize_name_into_remainder

# from elchempy.indexer.extra_EC_info import loading_ref, WE_surface_area_cm2

### 3rd Party imports

import datefinder
import dateutil
from dateutil.parser import ParserError

#%%
def _dev():
    ### for Developing
    from elchempy.config import LOCAL_FILES
    import pandas as pd

    fname = LOCAL_FILES[0]
    ecpps = []
    for fname in LOCAL_FILES:
        sid = ElChemPathParser(fname)
        ecpps.append(sid)
    aa = pd.concat(
        [pd.DataFrame(i.EC_info_entry, index=[0]).T for i in ecpps], ignore_index=True
    )
    bb = pd.concat([pd.DataFrame(i.EC_info_undetermined).T for i in ecpps])
    cc = pd.concat([pd.DataFrame(i.EC_INDEX_entry, index=[0]) for i in ecpps])

    ecpp = ElChemPathParser(fname)

    return aa, ecpps


class ElChemPathParserError(Exception):
    """raises errors from ElChemPathParser"""


class ElChemPathParser(Path):
    """
    Parses all the relevant information contained in only the filename (no file introspection).
        - Filepath parts and file statistics
        - Relevant tokens for the experimental conditions

    Prepares a dictionary for entry into the Index or database
    """

    _flavour = type(Path())._flavour

    name_separators = ["_", "-"]

    def __init__(self, *args, **kwargs):

        # First call the FilePathParser to create all the info from the filepath
        self.fpp = FilePathParser(self)
        # super().__init__(*args, **kwargs)
        if any("metadata" in i for i in kwargs.keys()):
            _kws = ", ".join([f"{k}={val}" for k, val in kwargs.items()])
            logger.warning(
                f"{self.__class__.__qualname__} does not do introspection for metadata\n{_kws}"
            )

        parent_info = ElChemPathParser.call_tokenizer_on_part(self.parent.name)
        date_from_parent = None
        if parent_info.get("date_dt", None):
            date_from_parent = parent_info.get("date_dt", None)

        stem_info = ElChemPathParser.call_tokenizer_on_part(
            self.stem, date_from_parent=date_from_parent
        )

        self.EC_info_parts = {"parent": parent_info, "stem": stem_info}
        self.EC_info_merged = ElChemPathParser.merge_info_parts(self.EC_info_parts)

        self.EC_info = {}
        if self.EC_info_merged:
            self.EC_info = {
                k: val for k, val in self.EC_info_merged.items() if val["options"] == 1
            }
            self.EC_info_undetermined = {
                str(self): {
                    k: ",".join(str(val["value"]))
                    for k, val in self.EC_info_merged.items()
                    if val["options"] > 1
                }
            }
            self.EC_info_undetermined[str(self)].update(
                **{"parent": self.parent.name, "stem": self.stem}
            )
        self.EC_info_entry = {
            **{"PAR_file": str(self)},
            **{k: val["value"] for k, val in self.EC_info.items()},
        }
        if self.fpp.FP_info.get("PAR_file") == str(self):

            self.EC_INDEX_entry = {
                str(self): {**self.EC_info_entry, **self.fpp.FP_info}
            }

    @staticmethod
    def call_tokenizer_on_part(name, date_from_parent=None):
        try:
            info = tokenize_name_into_remainder(name, date_from_parent=date_from_parent)
        except ElChemPathParserError as ex:
            raise ElChemPathParserError(
                f"Error in tokenizer of name: {str(name)}"
            ) from ex
            info = {"tokenize_error": ex, "tokenize_name": name}
        return info

    @staticmethod
    def merge_info_parts(EC_info_parts):
        """
        merges relevant info from filepath parts: filename or parent
        depending on the info key decides to take the info from which part
        """
        # {str(self):
        keys = EC_info_parts.keys()
        subkeys = set([i for val in EC_info_parts.values() for i in val.keys()])
        merged_info = {}
        for skey in subkeys:
            key_skey_vals = {k: EC_info_parts[k].get(skey, None) for k in keys}
            true_skeys = {k: val for k, val in key_skey_vals.items() if val}
            settruekeys = set(true_skeys.values())

            mcol = {skey: {"value": None, "source": None, "options": 0}}
            if not true_skeys:
                pass
                # merge_col.append(mcol)
            if len(true_skeys) == 1 or len(settruekeys) == 1:
                mcol[skey]["value"] = list(true_skeys.values())[0]
                mcol[skey]["source"] = list(true_skeys.keys())[0]
                mcol[skey]["options"] = 1
            elif len(true_skeys) > 1 and not len(settruekeys) == 1:

                if "electrode" in skey:
                    for n, (elk, elval) in enumerate(true_skeys.items()):
                        mcol[f"{skey}_{n}"] = {
                            "value": elval,
                            "source": elk,
                            "options": 1,
                        }
                elif "fname" in skey:
                    mcol[skey]["source"] = f'{skey}_{"_".join(true_skeys.keys())}'
                    mcol[skey]["value"] = "/".join(true_skeys.values())
                    mcol[skey]["options"] = 1
                elif skey in ["Electrolyte", "PAR_exp"] or skey.startswith("Loading"):
                    # take value for merge from stem
                    stemk = true_skeys.get("stem", None)
                    if stemk:
                        mcol[skey]["value"] = stemk
                        mcol[skey]["source"] = "stem"
                        mcol[skey]["options"] = 1

                else:
                    mcol[skey]["value"] = tuple(true_skeys.values())
                    mcol[skey]["source"] = tuple(true_skeys.keys())
                    mcol[skey]["options"] = len(true_skeys)

            merged_info.update(**mcol)
        return merged_info


#%%
