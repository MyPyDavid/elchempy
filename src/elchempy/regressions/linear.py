"""

holds functions for linear regressions

"""

## std lib
from typing import NamedTuple, Tuple, Dict
from collections import namedtuple

import logging

logger = logging.getLogger(__name__)

## local
# import elchempy
#
# from elchempy.dataloaders.fetcher import ElChemData
# from elchempy.experiments.N2.background_scan import contains_background_scan, get_N2_background_data

# from elchempy.regressions.linear import linregress_residual

## 3rd party
import numpy as np

# import pandas as pd
from scipy.stats import linregress, zscore

## constants


#%%


def linregress_residual(x_, y_) -> Tuple[Dict, Dict]:
    """
    linear regression of x and y values

    Parameters
    ----------
    x : array type
        DESCRIPTION.
    y : array type
        DESCRIPTION.

    Returns
    -------
    linfit_pars : dict
        contains the linear fit parameters.
    linfit_data : dict
        contains the linear fit data.
    """

    x = cast_in_array(x_)
    y = cast_in_array(y_)

    if not (x.any() and y.any()):
        return {}, {}

    # xcol="scanrate", ycol="j A/cm2"):
    _lindict = {}
    # _x, _y = gr[xcol], gr[ycol]
    _lin = linregress(x, y)
    _ymod = x * _lin.slope + _lin.intercept

    RMSD = (sum([(ym - y) ** 2 for y, ym in zip(_ymod, y)]) / len(y)) ** 0.5
    Z = zscore(abs(y))
    # !!! Optional, build in slice by Zscore before fitting...
    linfit_pars = {f"lin_{k}": getattr(_lin, k) for k in _lin._fields}
    linfit_data = {
        "data_x": x,
        "data_y": y,
        "lin_mod_y": _ymod,
        "lin_RMSD": RMSD,
        "lin_zscore_y": Z,
    }
    return linfit_pars, linfit_data


def cast_in_array(arg):
    """casts argument into np.array"""
    if isinstance(arg, np.ndarray):
        return arg

    msg = f"Can not cast type({type(arg)}), {arg} in array of floats.\n"
    try:
        arr = np.array(arg)
    except Exception as exc:
        raise exc(msg)

        arr = np.array(0)

    if not arr.dtype in (float, int):
        raise ValueError(msg + "Array is not of type float")

    if not arr.any():
        raise ValueError(msg + "Array is empty")

    return arr
