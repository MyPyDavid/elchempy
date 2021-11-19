"""
helpers for ORR experiments.
Filefinder for ring
"""


def _dev_test():
    # for developing
    from elchempy.experiments._dev_datafiles._dev_fetcher import (
        get_files,
    )

    files = get_files()
    disk_filepath = files[8]
    ring_file = find_file_from_secondary_electrode(disk_filepath)


def find_file_from_secondary_electrode(disk_filepath, search_folder=None):
    """
    In case of ORR experiment.

    gets file from secondary electrode in same folder as filepath
    with a similar basename plus suffix

    """

    secondary_files_options = set()

    disk_fstem = disk_filepath.stem
    if not search_folder:
        search_folder = disk_filepath.parent

    files_with_same_suffix = list(search_folder.rglob(f"*{disk_filepath.suffix}"))

    files_in_folder = [i for i in files_with_same_suffix if i != disk_filepath]

    files_with_same_length = {
        f for f in files_in_folder if len(f.stem) == len(disk_fstem)
    }
    secondary_files_options.update(files_with_same_length)

    if not secondary_files_options:
        return None

    if len(secondary_files_options) == 1:
        return list(secondary_files_options)[0]

    match_count = [
        (i, count_matching_consecutive_strings(disk_fstem, i.stem))
        for i in secondary_files_options
    ]
    best_match = max(match_count, key=lambda x: x[1])
    return best_match


def count_matching_consecutive_strings(stem1, stem2) -> int:
    """zips two strings and returns the count of matching consecutive elements"""
    n = 0
    for n, (i1, i2) in enumerate(zip(stem1, stem2)):
        if i1 == i2:
            pass
        else:
            break
    return n
