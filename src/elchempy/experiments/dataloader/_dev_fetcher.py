"""
Created on Thu Jul 15 16:12:27 2021

@author: DW
"""
from pathlib import Path
from reader import DataReader
from fetcher import ElchemData

def _dev():
    from pathlib import Path
    _n2files = Path.cwd().parent.parent.parent.parent.joinpath('data/raw').rglob('*par')
    return _n2files

def _dev_test_read():
    files = _dev()
    results = []
    # for filepath in files:
    while True:
        filepath = next(files)

        results.append(ElchemData(filepath))

    if False:
        DR = DataReader(filepath)
        actions = DR.actions
        data = DR.data

        data = assign_electrochemical_data_columns(data,
                              RHE_potential = 2,
                              geometric_SA = 20,
                              electrode_type = ''
                              )
        data = match_actions_data_segments(actions, data)

        results.append((DR, data, actions))

def _dev():
    _n2files = Path.cwd().parent.parent.parent.parent.joinpath('data/raw').rglob('*N2*par')
    return _n2files

def _test_read():
    files = _dev()
    results = []
    for filepath in files:
        DR = DataReader(str(filepath))
        results.append(DR)
        actions = DR.actions
        data = DR.data
    # if True:
        data.plot(x='E(V)',y='I(A)', title=filepath.name)
        if any('EIS' in name for name in actions.Name.unique()):
            data.plot(x='Z Real',y='Z Imag', title=filepath.name)
