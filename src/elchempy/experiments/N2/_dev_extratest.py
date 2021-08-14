"""
Created on Sat Aug 14 11:35:59 2021

@author: DW
"""

def _dev_test_type(arg):

    from pathlib import Path
    import pandas as pd

    from collections.abc import Collection, Mapping, Sequence, Generator, Iterable, MutableSequence, Hashable , Container, Callable

    _tests = {'list' : [1,2,3],'dict' : {'1': 2}, 'str': '123','pd.DataFrame' : pd.DataFrame(),
     'gen' : (i for i in [1,2]), 'tuple' : (1,2,3), 'Path' : Path.home(), 'def' : lambda x : x**2}

    _results_lst = []
    for tname,tobj in _tests.items():


        _res = {'typename' : tname, 'typetype.name' : type(tobj).__name__}
        for ttype in [Container, Hashable, Iterable, Collection, Mapping, Sequence, Generator, MutableSequence, Callable]:

            _isinst = isinstance(tobj, ttype)
            _not = 'not' if not _isinst else ''

            print(f'{tname} ({type(tobj).__name__}) is {_not} an instance of {ttype.__name__}')
            _res.update({ttype.__name__ : _isinst})
        _results_lst.append(_res)
    typeres = pd.DataFrame(_results_lst)
    if 0:
        typeres.query('(Iterable == False) & (Callable == False)')
