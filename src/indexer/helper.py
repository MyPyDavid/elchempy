# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:28:51 2020

@author: User
"""


import os
import glob
from functools import reduce
import datetime as dt

import numpy as np
import pandas as pd


def drop_cols(type_exp, DF):
    if type_exp == "ORR":
        N2_cols_drop_option = [
            "Segment #",
            "Point #",
            "Current Range",
            "Status",
            "Frequency(Hz)",
            "Z Real",
            "Z Imag",
            "ActionId",
            "AC Amplitude",
            "jcorr",
            "Gas",
            "EXP",
            "j_ring",
            "RPM",
            "Comment",
            "Measured_OCP",
            "Electrode",
        ]
    else:
        N2_cols_drop_option = [
            "Segment #",
            "Point #",
            "Current Range",
            "Status",
            "Frequency(Hz)",
            "Z Real",
            "Z Imag",
            "ActionId",
            "AC Amplitude",
            "jcorr",
            "Gas",
            "EXP",
            "j_ring",
            "RPM",
            "Comment",
            "Measured_OCP",
            "Electrode",
        ]
    drop_cols = [i for i in N2_cols_drop_option if i in DF.columns]
    DF = DF.drop(columns=drop_cols)
    return DF


#%%


def Plot_limits(DF, xMinMax, yMinMax):
    xlim = [DF[xMinMax].min(), DF[xMinMax].max() * 1.1]
    ylim = [np.min(DF[yMinMax].min() * 1.1), np.max(DF[yMinMax].max() * 1.1)]
    return xlim, ylim


def sigmoid(p, x):
    x0, y0, c, k = p
    y = c / (1 + np.exp(-k * (x - x0))) + y0
    return y


def residuals(p, x, y):
    return y - sigmoid(p, x)


def resize(arr, lower=0.0, upper=1.0):
    arr = arr.copy()
    if lower > upper:
        lower, upper = upper, lower
    arr -= arr.min()
    arr *= (upper - lower) / arr.max()
    arr += lower
    return arr
