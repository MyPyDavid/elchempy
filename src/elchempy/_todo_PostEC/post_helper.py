# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:42:09 2020

@author: User
"""
import numpy as np
import pandas as pd
from file_py_helper.PostChar import (
    SampleSelection,
    Characterization_TypeSetting,
    SampleCodesChar,
)


def make_uniform_EvRHE(df, rounding_set=2):
    lst = []
    for E in df["E_AppV_RHE"].values:
        match = 0
        for i in np.arange(-0.10, 2, 0.05):
            if np.isclose(i, E, atol=0.025):
                match = 1
                lst.append((E, i))
        if match == 0:
            if E < 0 and E > -0.04:
                lst.append((E, i))
            else:
                lst.append((E, np.nan))
    #                print(E,i)

    if len(df["E_AppV_RHE"].values) == len(lst):
        df = df.assign(**{"E_RHE": [np.round(float(i[1]), rounding_set) for i in lst]})
        print(
            'Len({0}) matches, new column: "E_RHE"'.format(len(df["E_AppV_RHE"].values))
        )
    else:
        print(
            "make_uniform_EvRHE lengths do not match LenEIS : {0}, len(lst) : {1}".format(
                len(df["E_AppV_RHE"].values), len(lst)
            )
        )
    return df


def CheckCols(coll_lst, df):
    return [i for i in coll_lst if i in df.columns]


def serie_model(
    EIS_pars, sIDslice, ModEEC="Model(Singh2015_R3RQ)", neat_eis=True, RPM_lim=1000
):
    ECB = EIS_pars.loc[EIS_pars.SampleID.isin(sIDslice)]
    cols = EIS_pars.columns
    if "Model_EEC" in cols:
        ECB = ECB.query(f'Model_EEC == "{ModEEC}"')
    if "RPM_DAC" in cols and RPM_lim:
        ECB = ECB.query(f"RPM_DAC > {RPM_lim}")
    if "ECexp" in cols and "ECuniq" not in cols:
        ECB = ECB.dropna(subset=["ECexp"])
        ECuniq = ["_".join(i.split("_")[:-1]) for i in ECB.ECexp.values]
        ECB = ECB.assign(**{"ECuniq": ECuniq})
        if "E_RHE" in cols:
            ECuniqE = [f"{i[0]}_{i[1]:.2f}" for i in zip(ECuniq, ECB.E_RHE.values)]
            ECB = ECB.assign(**{"ECuniqE": ECuniqE})

    if neat_eis and all([i in EIS_pars.columns for i in ["Rs", "RedChisqr", "Rct"]]):
        prel = len(ECB)
        RedChiSq_limit = (
            ECB.query("Rs > 1").RedChisqr.mean()
            + 1 * ECB.query("Rs > 1").RedChisqr.std()
        )
        ECB = ECB.query("RedChisqr < @RedChiSq_limit & Rs > 2 & Rct < 9E05")
        print(f"Cleaned up {prel-len(ECB)} rows")
    return ECB


def loading_series(pars):
    Loading_pars = pd.concat(
        [
            i[1]
            for i in [
                (n, gr, gr.Loading_cm2.unique())
                for n, gr in pars.query("RPM_DAC > 1000").groupby(
                    ["SampleID"] + SampleSelection.EC_exp_cols[0:-2]
                )
                if gr.Loading_cm2.nunique() > 2
            ]
        ]
    )
    return Loading_pars
