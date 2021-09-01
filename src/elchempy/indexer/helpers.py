"""
Helper functions related to file indexing.
"""

from pathlib import Path
from typing import Collection, Union, List

import logging

logger = logging.getLogger(__name__)


def find_path_of_common_folder(files) -> Union[None, Path]:
    """finds a common path in a list of files"""

    fparts = [Path(f).resolve().parts for f in files]

    maxlen = len(max(fparts, key=lambda x: len(x)))
    minlen = len(min(fparts, key=lambda x: len(x)))

    sets_per_part = [(i, list(set([fp[i] for fp in fparts]))) for i in range(0, minlen)]
    sets_len_is_1 = [i[1][0] for i in sets_per_part if len(i[1]) == 1]
    # [i for i in sets_len1 if i != '/']
    if not sets_len_is_1:
        return None

    if sets_len_is_1[0] == Path.home().root:
        # == '/' :
        # common_path = Path(Path.home().root).joinpath(Path("/".join(map(str,sets_len_is_1[1::] ))))
        common_path = Path(sets_len_is_1[0])
        for i in sets_len_is_1[1::]:
            common_path = common_path.joinpath(i)
        return common_path


def relative_parent_paths_to_common_folder(files, common_folder: Path) -> List[Path]:
    rel_paths = []
    for i in files:
        try:
            if isinstance(i, str):
                i = Path(i)
            relp = i.relative_to(common_folder).parent
        except ValueError as exc:
            relp = "destination_valuerror"
        except Exception as exc:
            logger.error(f"relative_parent_folders Unexpected errror {exc}")
            relp = "unkown_error"
        rel_paths.append(relp)
    return rel_paths
