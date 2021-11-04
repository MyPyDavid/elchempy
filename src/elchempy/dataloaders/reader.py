"""
empty placeholder
"""

from pathlib import Path

import pandas as pd
import logging

logger = logging.getLogger(__name__)

import elchempy
from elchempy.dataloaders.parse_filetypes.parser import parse_file
from elchempy.dataloaders.parse_filetypes.parser_helpers import (
    cast_parser_to_dict_of_dataframes,
    cast_metadata_to_flat_dict,
    get_starttime_from_parser,
)

#%%


class DataReader:
    """

    Wrapper for parsers in a class
    and handles the

    Class could be extended for reading out data
    from multiple source files.

    """

    supported_filetypes = [".par"]

    def __init__(self, filepath: Path, max_bytesize=1 * 10 ** 10, metadata_only=False):

        if not isinstance(filepath, Path):
            if isinstance(filepath, str):
                filepath = Path(filepath)
            else:
                raise TypeError("Argument given is not Path nor str")

        if not filepath.exists():
            raise FileNotFoundError(f"File does not exist:\n{filepath}")

        if filepath.suffix not in self.supported_filetypes:
            _warning = 'Filetype suffix is not supported, not in {", ".join(map(str,self.supported_filetypes))}'
            logger.warn(_warning)
            raise ValueError(_warning)

        filesize = filepath.stat().st_size
        if filesize > max_bytesize:
            _warning = f"File too large {filesize} > {max_bytesize}"
            logger.warn(_warning)
            raise ValueError(_warning)

        self.filepath = filepath
        self._metadata_only = metadata_only

        self.parser = None
        # Call the method read_file
        self.parser = parse_file(self.filepath, metadata_only=self._metadata_only)

        metadata = pd.DataFrame()
        actions = pd.DataFrame()
        data = pd.DataFrame()

        parser_dict = cast_parser_to_dict_of_dataframes(self.parser)
        self.parser_dict = parser_dict

        # metadata, actions, data = _frames
        self.actions = parser_dict["frames"]["actions"]
        self.data = parser_dict["frames"]["data"]

        flat_metadata = cast_metadata_to_flat_dict(parser_dict["dicts"]["metadata"])

        start_time, _start_time_source = get_starttime_from_parser(self.parser)
        self.start_time, self._start_time_source = start_time, _start_time_source
        flat_metadata.update({"meta_start_time": start_time})

        self.metadata = pd.DataFrame(flat_metadata, index=[str(self.filepath)])
        # self.double_check_data(self.data, expected_values=self.expected_keys_values)

    def __len__(self):
        if self._metadata_only:
            return len(self.metadata) + len(self.actions)
        else:
            return len(self.data)

    def __bool__(self):
        if self._metadata_only:
            return False if self.metadata.empty else True
        else:
            return False if self.data.empty else True

    def __repr__(self):
        _name = self.filepath.name
        _txt = f"actions = {len(self.actions)}, data = {len(self.data)}"
        return f"{self.__class__.__qualname__} on {_name}, {_txt}"
