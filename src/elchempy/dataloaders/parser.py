"""
Parsers for the experimental data files
"""

# import datetime
from pathlib import Path

from elchempy.dataloaders.parsers_filetypes.PAR_file import read_PAR_file

# 3rd party
# import datefinder
# import pandas as pd

#%%
class ParserError(ValueError):
    """unable to parse this file"""


def parse_file(filepath: Path, metadata_only=False):
    """
    Parameters
    ----------
    filepath : Path
        chooses the parsing method depending on the file suffix

    Returns
    -------
    parser_instance : object
        instance of parser

    """

    parser = None
    suffix = filepath.suffix

    if ".par" in suffix:
        parser = read_PAR_file(filepath, metadata_only=metadata_only)
        # this parser has dictionary attributes which need to be cast in to DataFrames

    elif ".other" in suffx:
        # add other possible filetypes here
        pass
    return parser
