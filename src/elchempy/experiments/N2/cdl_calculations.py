"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

import elchempy

from elchempy.experiments.dataloader._dev_fetcher import get_files, _dev_test_read


#        grB = N2_CVs.groupby(by=['Gas','Type','EXP'])
#        for scan in grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')):
#            print(scan)
#                grA.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','0.1'))
#        grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')).to_csv(N2_dest_dir.joinpath('%s.csv' %N2_fn))
def _test_data():
    datacollection = _dev_test_read(get_files("N2"))


def runner(datacollection):

    for ecdata in datacollection:
        N2CV, Cdl_data = N2_selection(ecdata)


def N2_analyisis(ecdata):
    # FIXME Select only CV types from Data segment
    # Select the data for N2 Cyclic Voltammograms
    N2_CVs = ecdata.data.loc[ecdata.data.ActionId == 38]

    N2_CVs = N2_CVs.loc[N2_CVs.scanrate_calc != 0]

    # scanrates = N2_CVs.scanrate.unique()
    if N2_CVs.scanrate.nunique() > 2:
        # if multiple scanrates are present than calculations can start
        Cdl_data, Cdl_pars = CDL(N2_CVs)

    # Check if possible background (BG) scan is in the data
    BG_present = False
    if N2_CVs.scanrate.min() < 0.011:
        # check presence of slow scanrates in data
        BG_present = True
    if BG_present:
        N2_BG = get_N2_background_data(N2_CVs)

    return N2_CVs, Cdl_data


EvRHE = "E_AppV_RHE"


def CDL(N2_CVs):

    scanrates = N2_CVs.scanrate.unique()

    srgrpby = N2_CVs.groupby("scanrate")

    V08 = N2_CVs.loc[np.isclose(0.8, N2_CVs[EvRHE], atol=0.01)]
    grp = N2_CVs.groupby(["Segment #", "SweepType", "E_AppV_RHE"])

    _N2_scans_lst = []
    for sr in scanrates:
        sr_mV = int(sr * 1000)

        N2sr = srgrpby.get_group(sr)
        segments = N2sr["Segment #"].unique()

        N2sr_lastseg = check_if_last_seg_isnormal(N2sr)
        if len(segments) > 1:
            pass  # TODO do sth with 1st seg

        _N2_scans_lst.append(N2sr_lastseg)
    N2_scans = pd.concat(_N2_scans_lst)
    Cdl_data = prepare_Cdl_frame(N2_scans)
    j_key = prepare_Cdl_frame.__defaults__[0]  # !!! check definition for key
    # Cdl_data.groupby("SweepType").plot(x='scanrate',y=j_key,kind='scatter')

    # Make a linear test and filter somehow???
    Cdl_data = make_Cdl_Ztest(Cdl_data)

    if False:
        for ycol in ["lin_slope"]:
            Cdl_data.query("lin_zscore_y < 3").groupby("SweepType").plot(
                x=EvRHE, y="lin_intercept"
            )

    Cdl_data.query("lin_zscore_y < 3")
    if False:
        Resmean = Cdl_data.groupby([EvRHE, "SweepType", "scanrate"]).mean()
        Resmean = Resmean.query("lin_Res < 1E-3 & zscore_y < -0.5 & scanrate > 0.02")

    make_baseline_corr = True
    if make_baseline_corr:
        Cdl_data = check_for_linear_baseline_correction_of_Cdl_values(Cdl_data)

    Cdl_pars = get_Cdl_pars(Cdl_data)
    return Cdl_data, Cdl_pars
    # READY!!!
    # TODO SEND TO PLOTS AND EXPORT


def get_Cdl_pars(Cdl_data, E_linspace=np.linspace(0, 1, 21)):

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


def make_Cdl_Ztest(
    Cdl_data, xcol="scanrate", ycol="j A/cm2", grpbykeys=[EvRHE, "SweepType"]
):

    # gr = gr.assign(**_lindict)
    Cdl_test = pd.concat(
        [
            gr.assign(**linregress_residual(gr[xcol], gr[ycol]))
            for (Ev, nm), gr in Cdl_data.groupby(by=grpbykeys)
        ]
    )
    return Cdl_test


def linregress_residual(x, y):
    # xcol="scanrate", ycol="j A/cm2"):
    _lindict = {}
    # _x, _y = gr[xcol], gr[ycol]
    _lin = linregress(x, y)
    _ymod = x * _lin.slope + _lin.intercept
    _Res = abs(y) - abs(_ymod)
    Z = zscore(y)
    _lindict = {f"lin_{k}": getattr(_lin, k) for k in _lin._fields}
    _lindict.update({"lin_mod": _ymod, "lin_Res": _Res, "lin_zscore_y": Z})
    return _lindict

    # return gr


def check_for_linear_baseline_correction_of_Cdl_values(Cdl_data):
    _swplst = []
    for swpnm, swpgrp in Cdl_data.groupby("SweepType"):

        swpgrp_r085 = swpgrp.query("lin_rvalue > 0.85")

        lin_baseline = linregress(swpgrp_r085[EvRHE], swpgrp_r085["lin_slope"])

        lin_baseline_model = swpgrp[EvRHE] * lin_baseline.slope + lin_baseline.intercept

        lin_slope_baseline_corr = (
            swpgrp["lin_slope"] - lin_baseline_model + swpgrp["lin_slope"].mean()
        )

        swpgrp = swpgrp.assign(**{"lin_slope_baseline_corr": lin_slope_baseline_corr})
        _swplst.append(swpgrp)
    Cdl_data = pd.concat(_swplst)
    return Cdl_data


def prepare_Cdl_frame(N2_scans, current_density_key="j A/cm2", potential_key=EvRHE):
    # Loop over Sweep Types
    _results = []
    for (sr, swpname), swpgrp in N2_scans.query("scanrate > 0").groupby(
        ["scanrate", "SweepType"]
    ):
        if swpname == ("chrono", "NA"):
            continue
        # Loop over potentials E This mean for example for E in [0.1, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.9]:
        for E in np.linspace(0.1, 1, 50):

            j_Cdl = np.abs(
                swpgrp.loc[(np.isclose(swpgrp[EvRHE], E, atol=0.010))][
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


def check_if_last_seg_isnormal(N2sr):
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
            _minmean = _N2sr_segs.groupby("Segment #").agg("jmAcm-2").min().mean()

            _N2sr_test = N2sr.loc[
                N2sr["Segment #"] == uniq_segs[_idx],
            ]
            _last_mmn = _N2sr_test["jmAcm-2"].min().mean()
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
