"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

## std lib
from typing import NamedTuple, Tuple, Dict
from collections import namedtuple
from pathlib import Path

import logging

logger = logging.getLogger(__name__)

## local
import elchempy

# from elchempy.dataloaders.fetcher import ElChemData
# from elchempy.experiments.N2.background_scan import contains_background_scan, get_N2_background_data

from elchempy.experiments.N2.plotting import N2_plot_raw_scans_scanrate

## 3rd party
import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

## constants
EvRHE = "E_vs_RHE"

#%%


def N2_Cdl_calculations(
    N2_CVs: pd.DataFrame,
    EvRHE: str,
    current_density_key: str = "j_A_cm2",
    potential_key: str = EvRHE,
):
    """performs the calculations for the capacity (double layer) on the data"""

    if N2_CVs.empty:
        return None

    scanrates = N2_CVs.scanrate.unique()

    srgrpby = N2_CVs.groupby("scanrate")

    # _V08 = N2_CVs.loc[np.isclose(0.8, N2_CVs[EvRHE], atol=0.01)]
    grp = N2_CVs.groupby(["Segment #", "SweepType", EvRHE])

    _N2_scans_lst = []
    for sr in scanrates:

        sr_mV = int(sr * 1000)
        N2sr = srgrpby.get_group(sr)
        segments = N2sr["Segment #"].unique()

        N2sr_lastseg = check_if_last_segment_is_normal(N2sr)
        if len(segments) > 1:
            pass  # TODO do sth with 1st seg
        _N2_scans_lst.append(N2sr_lastseg)

    N2_scans = pd.concat(_N2_scans_lst)
    Cdl_dataprep = prepare_Cdl_frame(N2_scans)
    j_key = prepare_Cdl_frame.__defaults__[0]  # !!! check definition for key
    # Cdl_data.groupby("SweepType").plot(x='scanrate',y=j_key,kind='scatter')

    # Make a linear test and filter somehow???
    Cdl_pars, Cdl_data = make_cdl_pars_data_from_linregress(Cdl_dataprep)

    if False:
        for ycol in ["lin_slope", "lin_intercept"]:
            Cdl_pars.groupby("SweepType").plot(x=EvRHE, y=ycol)
    # Cdl_data.query("lin_zscore_y < 3")
    if False:
        Resmean = Cdl_data.groupby([EvRHE, "SweepType"]).mean()
        # Resmean = Resmean.query("lin_rmsd < 1E-3 & lin_zscore_y < -0.5")

    make_baseline_corr = True
    if make_baseline_corr:
        Cdl_pars = check_for_linear_baseline_correction_of_Cdl_values(Cdl_pars)

    Cdl_pars = get_Cdl_pars_per_E(Cdl_pars)
    # TODO SEND TO PLOTS AND EXPORT
    return Cdl_pars, Cdl_data


def get_Cdl_pars_per_E(
    Cdl_data: pd.DataFrame, E_linspace=np.linspace(0, 1, 51)
) -> pd.DataFrame:
    """selects certain potential E values from the data"""

    # _swplst = []
    Cdl_swp_E_selection_lst = []
    for swpnm, swpgrp in Cdl_data.groupby("SweepType"):
        Cdl_swp_E_selection = []
        for Et in E_linspace:
            try:
                Cdl_swp_E_selection_lst.append(
                    swpgrp.loc[
                        (np.isclose(swpgrp[EvRHE], Et, atol=50e-3) == True),
                        :,
                    ].head(1)
                )
            except Exception as e:
                pass
    Cdl_swp_E_selection = pd.concat(Cdl_swp_E_selection_lst)
    return Cdl_swp_E_selection


def _plot(Cdl_data):
    for ycol in ["lin_slope", "lin_slope_baseline_corr"]:
        Cdl_data.groupby("SweepType").plot(x=EvRHE, y=ycol, kind="scatter")


def make_cdl_pars_data_from_linregress(
    Cdl_dataprep: pd.DataFrame,
    xcol="scanrate",
    ycol="j_A_cm2",
    grpbykeys=[EvRHE, "SweepType"],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parameters
    ----------
    Cdl_dataprep : pd.DataFrame
        A frame that is prepared from the N2_CV raw data.
    xcol : str, optional
        DESCRIPTION. The default is "scanrate".
    ycol : str, optional
        DESCRIPTION. The default is "j_A_cm2".
    grpbykeys : list, optional
        List of keys to loop the groupby over. The default is [EvRHE, "SweepType"].

    Returns
    -------
    cdl_pars : pd.DataFrame
        contains the linear fit parameters(pars).
    cdl_data : pd.DataFrame
        contains the linear fit data columns.

    """
    cdl_pars, cdl_data = pd.DataFrame(), pd.DataFrame()
    # gr = gr.assign(**_lindict)
    if all(col in Cdl_dataprep.columns for col in [xcol, ycol]):

        cdl_pars_lst, cdl_data_lst = [], []
        for (Ev, nm), gr in Cdl_dataprep.groupby(by=grpbykeys):
            if not gr[[xcol, ycol]].dropna().empty:
                grpkeyval = dict(zip(grpbykeys, (Ev, nm)))

                pars, data = linregress_residual(gr[xcol], gr[ycol])
                # gr_pars = gr.drop(columns=[xcol, ycol]).assign(**pars)
                gr_pars = pd.DataFrame(pars, index=gr.index[0:1]).assign(**grpkeyval)
                cdl_pars_lst.append(gr_pars)

                gr_data = pd.DataFrame(data).assign(**grpkeyval)
                cdl_data_lst.append(gr_data)

        if cdl_data_lst:
            cdl_data = pd.concat(cdl_data_lst)
        else:
            logger.warning("make_linregress_Cdl_Ztest error, empty data")

        if cdl_pars_lst:
            cdl_pars = pd.concat(cdl_pars_lst)
        else:
            logger.warning("make_linregress_Cdl_Ztest error, empty pars ")

    else:
        logger.warning("make_linregress_Cdl_Ztest error, missing columns")

    return cdl_pars, cdl_data


def linregress_residual(x, y) -> Tuple[Dict, Dict]:
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
    # xcol="scanrate", ycol="j A/cm2"):
    _lindict = {}
    # _x, _y = gr[xcol], gr[ycol]
    _lin = linregress(x, y)
    _ymod = x * _lin.slope + _lin.intercept
    _rmsd = np.sqrt((y - _ymod) ** 2)
    Z = zscore(abs(y))
    # !!! Optional, build in slice by Zscore before fitting...
    linfit_pars = {f"lin_{k}": getattr(_lin, k) for k in _lin._fields}
    linfit_data = {
        "data_x": x,
        "data_y": y,
        "lin_mod_y": _ymod,
        "lin_rmsd": _rmsd,
        "lin_zscore_y": Z,
    }
    return linfit_pars, linfit_data


def check_for_linear_baseline_correction_of_Cdl_values(Cdl_pars) -> pd.DataFrame:
    """takes a linear fit over the Cdl values and substracts it as baseline correction"""
    _swplst = []
    for swpnm, swpgrp in Cdl_pars.groupby("SweepType"):

        swpgrp_r085 = swpgrp.query("lin_rvalue > 0.85")

        if swpgrp_r085.empty:
            continue

        lin_baseline = linregress(swpgrp_r085[EvRHE], swpgrp_r085["lin_slope"])

        lin_baseline_model = swpgrp[EvRHE] * lin_baseline.slope + lin_baseline.intercept

        lin_slope_baseline_corr = (
            swpgrp["lin_slope"] - lin_baseline_model + swpgrp["lin_slope"].mean()
        )

        swpgrp = swpgrp.assign(**{"lin_slope_baseline_corr": lin_slope_baseline_corr})
        _swplst.append(swpgrp)
    Cdl_data = pd.concat(_swplst)
    return Cdl_data


def prepare_Cdl_frame(
    N2_scans: pd.DataFrame,
    current_density_key: str = "j_A_cm2",
    potential_key: str = EvRHE,
    E_linspace=np.linspace(0.1, 1, 50),
) -> pd.DataFrame:
    """
    Prepares the raw data from CV at several scanrates to frame with
    a only relevant columns and potentials from a linspace.

    Parameters
    ----------
    N2_scans : pd.DataFrame
        DESCRIPTION.
    current_density_key : str, optional
        DESCRIPTION. The default is "j_A_cm2".
    potential_key : str, optional
        DESCRIPTION. The default is EvRHE.

    Returns
    -------
    Cdl_data : pd.DataFrame
        DESCRIPTION.

    """
    # Loop over Sweep Types
    _results = []
    for (sr, swpname), swpgrp in N2_scans.query("scanrate > 0").groupby(
        ["scanrate", "SweepType"]
    ):

        if swpname == ("chrono", "NA"):
            # skip these sweep types
            continue
        # Loop over potentials E This mean for example for E in [0.1, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.9]:
        for E in E_linspace:

            j_Cdl = np.abs(
                swpgrp.loc[(np.isclose(swpgrp[potential_key], E, atol=0.010))][
                    current_density_key
                ]
            ).mean()
            #                    Cdl_std = np.abs(Cdl_scan.loc[(np.isclose(Cdl_scan[EvRHE],E,atol=0.010)) & (Cdl_scan['Sweep_Type'] == sweep)]['j A/cm2']).std()/SR
            _results.append(
                {
                    "scanrate": sr,
                    f"{current_density_key}": j_Cdl,
                    "SweepType": swpname,
                    f"{potential_key}": E,
                    # "N2_CV_datafile": N2sr_fn_out,
                }
            )
    Cdl_data = pd.DataFrame(_results)
    return Cdl_data


def check_if_last_segment_is_normal(N2sr: pd.DataFrame) -> pd.DataFrame:
    """checks if the scan of the last segment is "normal" similar to the other scans at the same scan rate"""
    uniq_segs = N2sr["Segment #"].unique()
    nuniq_segs = N2sr["Segment #"].nunique()
    if nuniq_segs > 1:
        _N2sr_DF_out = pd.DataFrame()
        _idx = -1
        while _N2sr_DF_out.empty and abs(_idx) < nuniq_segs:
            _lastsegs = uniq_segs[-(nuniq_segs - int(nuniq_segs * 0.8)) :]
            _N2sr_segs = N2sr.loc[
                N2sr["Segment #"].isin(_lastsegs),
            ]
            _minmean = _N2sr_segs.groupby("Segment #").agg("j_mA_cm2").min().mean()

            _N2sr_test = N2sr.loc[
                N2sr["Segment #"] == uniq_segs[_idx],
            ]
            _last_mmn = _N2sr_test["j_mA_cm2"].min().mean()
            _dev_perc = 100 * (_minmean - _last_mmn) / _minmean
            #                            print(f'{_dev_perc}')
            if abs(_dev_perc) < 25:
                _N2sr_DF_out = _N2sr_test
            else:
                _idx -= 1
    else:
        _N2sr_DF_out = N2sr.loc[
            N2sr["Segment #"] == uniq_segs[0],
        ]
    return _N2sr_DF_out
