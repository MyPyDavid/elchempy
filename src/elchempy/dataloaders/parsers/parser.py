"""
Parsers for the experimental data files
"""

# import datetime
from pathlib import Path

from elchempy.dataloaders.parsers.filetypes.PAR import read_PAR_file

#%%

# Parsers register
SUPPORTED_FILETYPES = {".par": read_PAR_file}

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

    suffix = filepath.suffix

    if suffix not in SUPPORTED_FILETYPES.keys():
        raise NotImplementedError(
            f'parse_file error, File with suffix {suffix} is not supported ({", ".join(SUPPORTED_FILETYPES .keys())}'
        )

    reader = SUPPORTED_FILETYPES.get(suffix, None)

    if not callable(reader):
        raise NotImplementedError(
            f"selected reader is not callable, File with suffix {suffix} is not supported."
        )

    try:
        # calling the selected read function
        parsed_file = reader(filepath, metadata_only=metadata_only)
    except Exception as exc:
        raise exc from exc
        parsed_file = None

    return parsed_file
