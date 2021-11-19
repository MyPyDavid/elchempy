#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path
import datetime as dt
import hashlib
from typing import Dict

import logging

logger = logging.getLogger(__name__)

### local imports
#%%


class FilePathParser(Path):
    """Subclass of path, creates a dictonary with filepath info (FP_info)"""

    _flavour = type(Path())._flavour

    def __init__(self, *args, **kwargs):
        super().__init__()
        self._qcnm = self.__class__.__qualname__
        # self.stats_ = None
        self.data = None
        self.FP_info = {"FilePath": str(self), "file_exists": self.exists()}
        if self.exists() and self.is_file():
            self.FP_info.update(
                **{
                    "FileID": self.get_rfID_from_path(self),
                    "FileStem": self.stem,
                    "PAR_file": str(self),
                    "FileParent": self.parent.name,
                }
            )
            self.filestats = self.get_dict_from_fstat(self)

            self.FP_info.update(**self.filestats)

    @staticmethod
    def get_fID_from_path(path: Path) -> str:
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
