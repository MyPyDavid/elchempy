"""
Created on Thu Jul 15 16:12:27 2021

@author: DW
"""
from pathlib import Path

from elchempy.dataloaders.reader import DataReader
from elchempy.dataloaders.fetcher import ElChemData
from elchempy.config import LOCAL_FILES


__all__ = ["_dev_get_files", "_dev_test_read", "_test_read"]


def get_files(name="", files=LOCAL_FILES):
    """returns all local data files for testing purposes"""

    _files = list([i for i in files if name in str(i)])
    if not _files:
        print("Warning, no files with name {name} found in:\n{datadir}")
    return _files


# from pathlib import Path
# _search = '*par'
# if name:
#     _search = f'**/**/*{name}*par'
# rel_data_folder = 'data/raw'
# CWD = Path(__file__)
# print(f'CURRENT WD: {CWD}')1
# if 'src' in CWD.parts:
#     _src_idx = [n for n,i in enumerate(CWD.parts) if i == 'src'][0]
# repodir = Path('/'.join(CWD.parts[0:_src_idx]))
# datadir= repodir.joinpath(rel_data_folder)
# print(datadir)


def _dev_test_read(files):
    # files = _dev()
    results = []
    # for filepath in files:
    for file in files:
        results.append(ElChemData(file))

    while False:
        try:
            filepath = next(filesgen)
            results.append(ElChemData(filepath))
        except StopIteration:
            print(f"data fetch finished len {len(results)}")
            break
    return results


def _false():
    if False:
        DR = DataReader(filepath)
        actions = DR.actions
        data = DR.data

        data = assign_electrochemical_data_columns(
            data, RHE_potential=2, geometric_SA=20, electrode_type=""
        )
        data = match_actions_data_segments(actions, data)

        results.append((DR, data, actions))


def _dev():
    _n2files = (
        Path.cwd().parent.parent.parent.parent.joinpath("data/raw").rglob("*N2*par")
    )
    return _n2files


def _test_read_files():
    """reads all test files and plots the data"""
    import matplotlib.pyplot as plt

    files = get_files()
    results = []
    for filepath in files:
        DR = DataReader(str(filepath))
        results.append(DR)
        actions = DR.actions
        data = DR.data
        # if True:
        try:
            fig, ax = plt.subplots()

            data.groupby("Segment #").plot(
                x="E(V)", y="I(A)", title=filepath.name, ax=ax
            )
            if any(i for i in actions.Name.unique() if "EIS" in i):
                data.groupby("Segment #").plot(
                    x="Z Real", y="Z Imag", kind="scatter", title=filepath.name, ax=ax
                )
            plt.show()
            plt.close()
        except Exception as exc:
            print(f"Plotting error for {filepath}:\n{exc}")
