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


def merge_SamplesOVV(top_data_folder, dest_plot_folder, O2, Cdl, HPRR):
    # O2, Cdl, HPRR = overall_rows, Cdl_fit_data, HPRR_pars
    #    O2,Cdl,HPRR = CALC_PARS_svd, Cdl_fit_svd, HPRR_svd
    DW_Samples = glob.glob(
        "%s\\Project*\\MG*\\DW_samples*.xlsx"
        % os.path.dirname(os.path.dirname(os.path.dirname(top_data_folder)))
    )[0]
    BET_Samples = glob.glob(
        "%s\\Organized_Data\\BET_Plots\\All_Samples_Overview.*"
        % os.path.dirname(os.path.dirname(top_data_folder))
    )[0]
    dfs = [i for i in [O2, Cdl, HPRR] if not i.empty]
    df_final = reduce(
        lambda left, right: pd.merge(left, right, on="SampleID", how="outer"), dfs
    )
    try:
        #        Samples_OVV = pd.read_excel(DW_Samples)
        #        BET_OVV = pd.read_excel(BET_Samples)
        #        OVV = pd.merge(Samples_OVV,BET_OVV,on='SampleID',how='outer')
        #        out = pd.merge(OVV,df_final,on='SampleID',how='outer')

        # make_sure_path_exists(dest_plot_folder + "\\EC_per_Sample_OVV")
        today_date = dt.datetime.today()
        for nm, gr in df_final.groupby(by="SampleID"):
            gr.to_csv(
                dest_plot_folder
                + "\\EC_per_Sample_OVV"
                + "\\EC_PARS_OVV_%s_%s.csv" % (today_date, nm)
            )
        print(
            "Merged to: %s" % (dest_plot_folder + "\\" + "PARS_OVV_%s.csv" % today_date)
        )
    except Exception as e:
        print("Problem merging Samples with BET. %s" % e)
    return df_final
