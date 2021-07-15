"""
Created on Tue Jul 13 09:35:26 2021

@author: DW
"""
from pathlib import Path
import warnings

import pytest
import unittest

import elchempy


def get_test_data_folder():
    return Path(__file__).parent.parent.joinpath('data/raw')

def get_par_files(get_test_data_folder):
    return list(get_test_data_folder.rglob('*.par'))



@pytest.fixture(scope="session")
def test_data_files():
    files = get_par_files(get_test_data_folder())
    if not files:
        warnings.warn('No test data files were found, skipping these tests')
    return files

test_data_files = get_par_files(get_test_data_folder())