#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from pathlib import Path
import datetime as dt
import hashlib

# import pandas as pd
# from bs4 import BeautifulSoup


# import lxml


import logging

logger = logging.getLogger(__name__)
# from functools import partial
# import multiprocessing
# from os import cpu_count
# import hashlib
# import lxml
# import logging
# import numpy as np
# import concurrent.futures
# from multiprocessing import Pool
from typing import Dict

### local imports
# from filename_helpers import sID_to_sgrpID, get_fstats

### 3rd party imports
# Special file-py-helper
from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations

# from file_py_helper.FindSampleID import GetSampleID

### for Developing
import elchempy
from elchempy.config import LOCAL_FILES

# tt = PAR_file_parser(cc.par_files_run[42])
# t1 = '/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT/2018-01-Jan/19.01.2018_DW18_bipot_0.5MH2SO4_RRDE22434/OCP_RHE.par'
# tt2 = PAR_file_parser(t1)
#%%


def _index_structure():
    dr = {
        "EXP_dir": exp_dir_path,
        "Dest_dir": Dest_dir,
        "EXP_date": exp_date,
        "PAR_date": PAR_date,
        "EXP_PAR_date_match": exp_DATE_match,
        "EXP_PAR_folder_match": sID_dict["sID_Match_file_dir"],
        "PAR_file": PARf,
        "PAR_hash": hash_f,
        "Gas": gas,
        "PAR_exp": tpnm,
        "postAST": postAST,
        "filesize": PARf.stat().st_size,
        "basename": basepf,
        "SampleID": sID_file,
        "SampleID_folder": sID_dict["sID_Folder"],
        "Creation_date": cr_date,
        "LastMod_date": mod_date,
        "Delta_PAR_LastMod": PAR_cr_Delta,
        "Electrode": electrode,
        "pH": pH["pH"][0],
        "Electrolyte": pH["Electrolyte"][0],
        "Comment": comment_act0_DF,
        "Loading_name": loading_name,
        "Loading_cm2": loading_cm2,
    }


index_file_primary_keys = {"fID": "string"}

index_file_path_keys = {"FileStem": "string", "FilePath": "Path", "PAR_file": "Path"}

index_folder_path_keys = {
    "DIR_name": "Path",
    "DIR_dest_name": "Path",
    "DIR_date": "datetime.date",
}

index_file_date_keys = {
    "PAR_date": "datetime.date",
    "PAR_introspec_data": "datetime.date",
}

index_file_sample_keys = {
    "SampleID": "string",
    "SampleGroup": "string",
}

index_file_read_text_keys = {"FileHash": "string", "FileText": "string"}

index_dtypes_collection = {
    **index_file_path_keys,
    **index_file_sample_keys,
    **index_file_read_text_keys,
}

# Extra name to sID mapper, if keys is in filename
# _extra_sID_name_mapper = {
#     "David": "DW",
#     "stephen": "SP",
#     "Alish": "AS",
#     "Aish": "AS"}
# Extra name to sID mapper, if key is in filepath parts
# _extra_sgrpID_name_mapper = {"Raman Data for fitting David": "SH"}


def _def_parse():
    for f in LOCAL_FILES:
        fp = FilePathParser(f)


class FilePathParser(Path):
    """parses the filename of an experimental file"""

    _flavour = type(Path())._flavour

    def __init__(self, *args, **kwargs):
        super().__init__()
        self._qcnm = self.__class__.__qualname__
        # self.stats_ = None
        self.data = None
        self.FP_info = {"FilePath": self, "file_exists": self.exists()}
        if self.exists() and self.is_file():
            self.FP_info.update(
                **{
                    "fID": self.get_rfID_from_path(self),
                    "FileStem": self.stem,
                    "PAR_file": str(self),
                    "FileParent": self.parent.name,
                }
            )
            self.filestats = self.get_dict_from_fstat(self)

            self.FP_info.update(**self.filestats)
            # self.parse_result = self.collect_parse_results(**kwargs)

    @staticmethod
    def get_rfID_from_path(path: Path) -> str:
        """
        Makes the ID from a filepath

        Parameters
        ----------
        path : Path
            DESCRIPTION.

        Returns
        -------
        str: which contains hash(parent+suffix)_stem of path
        """
        _parent_suffix_hash = hashlib.sha256(
            (str(path.parent) + path.suffix).encode("utf-8")
        ).hexdigest()
        _filestem = path.stem
        fnID = _parent_suffix_hash + "_" + _filestem
        return fnID

    def read_data(self, file):
        """optional read data function"""
        try:
            self.data = pd.DataFrame()
            # TODO optionally elchemdata
        except Exception as exc:
            logger.warning(f"{self._qcnm} {self} data reader failed.\n{exc}")

    @staticmethod
    def get_dict_from_fstat(filepath: Path) -> Dict:
        """converting creation time and last mod time to datetime object"""
        fstat = filepath.stat()
        index_file_stat_keys = {
            "FileCreationDate": "datetime64",
            "FileCreation": "float",
            "FileModDate": "datetime64",
            "FileMod": "float",
            "FileSize": "int64",
        }

        c_t = fstat.st_ctime
        m_t = fstat.st_mtime
        c_tdate, m_tdate = c_t, m_t
        try:
            c_t = dt.datetime.fromtimestamp(fstat.st_ctime)
            m_t = dt.datetime.fromtimestamp(fstat.st_mtime)
            c_tdate = c_t.date()
            m_tdate = m_t.date()
        except OverflowError as e:
            pass
        except OSError as e:
            pass
        stats = c_tdate, c_t, m_tdate, m_t, fstat.st_size
        return dict(zip(index_file_stat_keys.keys(), stats))


def _depr_collect_parse_results(
    self, read_data=False, store_data=False, **kwargs
) -> Dict:
    """performs all the steps for parsing the filepath"""
    parse_res_collect = {}

    self.stats_ = self.stat()

    _fnID = self.make_dict_from_keys(
        index_file_primary_keys, (self.get_rfID_from_path(self),)
    )
    _filepath = self.make_dict_from_keys(index_file_path_keys, (self.stem, self, self))
    # _sample = self.parse_sample_with_checks()
    _filestats = self.parse_filestats(self.stats_)
    if read_data == True:
        data = self.read_data(self)

        parse_res_collect = {**_fnID, **_filepath, **_filestats}
    else:
        logger.warning(f"{self._qcnm} {self} is not a file => skipped")
    # else:
    # logger.warning(f"{self._qcnm} {self} does not exist => skipped")
    return parse_res_collect


def _depr_parse_filestats(self, fstat) -> Dict:
    """get status metadata from a file"""

    filestats = get_fstats(fstat)
    return self.make_dict_from_keys(index_file_stat_keys, filestats)


def _depr_make_dict_from_keys(self, _keys_attr: Dict, _result: tuple) -> Dict:
    """returns dict from tuples of keys and results"""
    if not isinstance(_result, tuple):
        logger.warning(
            f"{self._qcnm} input value is not a tuple, {_result}. Try to cast into tuple"
        )
        _result = (_result,)

    _keys = _keys_attr.keys()

    if not len(_result) == len(_keys) and not isinstance(_keys, str):
        # if len not matches make stand in numbered keys
        _keys = [f"{_keys_attr}_{n}" for n, i in enumerate(_result)]
    return dict(zip(_keys, _result))


def _depr_parse_sample_with_checks(self):
    """parse the sID and sgrpID from stem"""

    _parse_res = filestem_to_sid_and_pos(self.stem)

    if len(_parse_res) == 2:
        sID, position = _parse_res

        try:
            sID = _extra_overwrite_sID_from_mapper(sID)
        except Exception as exc:
            logger.info(
                f"{self._qcnm} {self} _extra_overwrite_sID_from_mapper failed => skipped.\n{exc}"
            )

        sgrpID = sID_to_sgrpID(sID)

        try:
            sgrpID = _extra_overwrite_sgrpID_from_parts(self.parts, sgrpID)
        except Exception as exc:
            logger.info(
                f"{self._qcnm} {self} _extra_overwrite_sgrpID_from_parts failed => skipped.\n{exc}"
            )

        _parse_res = sID, position, sgrpID
    else:
        logger.warning(
            f"{self._qcnm} {self} failed to parse filename to sID and position."
        )
    return self.make_dict_from_keys(index_file_sample_keys, _parse_res)
