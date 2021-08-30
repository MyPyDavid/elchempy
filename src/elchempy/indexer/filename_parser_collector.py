from pathlib import Path
from typing import List, Collection
import logging

logger = logging.getLogger(__name__)

from elchempy.indexer.filename_parser import FilePathParser
from elchempy.indexer.EC_filepath_parser import ElchemPathParser


def make_collection(files: Collection, **kwargs) -> List[ElchemPathParser]:

    ecpp_collection = []
    # _error_parse_filenames
    for file in files:
        try:
            ecpp = ElchemPathParser(file, **kwargs)
            ecpp_collection.append(ecpp)
        except Exception as e:
            _err = {"file": file, "error": e, "kwargs": kwargs}
            logger.warning(
                f"{__name__} make_collection unexpected error for calling PathParser on\n{file}.\n{e}"
            )
            # .append(file)
    ecpp_collection = sorted(ecpp_collection)
    return ecpp_collection
    # self.pp_collection = pp_collection
