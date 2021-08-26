"""

this module calculates the Cdl from Cyclic Voltammetries measured in N2 at several scanrates

"""

from typing import NamedTuple
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress, zscore

import logging

logger = logging.getLogger(__name__)


import elchempy

from elchempy.experiments.dataloader.fetcher import ElChemData
from elchempy.experiments.N2.background_scan import get_N2_background_data

#        grB = N2_CVs.groupby(by=['Gas','Type','EXP'])
#        for scan in grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')):
#            print(scan)
#                grA.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','0.1'))
#        grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')).to_csv(N2_dest_dir.joinpath('%s.csv' %N2_fn))


class _DevClass:
    from elchempy.experiments._dev_datafiles._dev_fetcher import (
        get_files,
        _dev_test_read,
    )

    def _test_data(exp_type: str):
        datacollection = _dev_test_read(get_files(exp_type))
        return datacollection

    def _test_runner():

        _results = []
        for ecdata in _test_data("N2"):
            ecdata = N2_analysis(ecdata)
            _results.append(ecdata)
        return _results


def new_runner():
    from elchempy.experiments.dataloader._dev_fetcher import get_files, _dev_test_read

    _result = []
    for fl in get_files("O2"):
        N2res = N2_Data(fl)
        _result.append(N2res)
    return _result


# N2_results = namedtuple('N2', 'raw_data pars data N2_BG')


class N2_Results(NamedTuple):
    raw_data: pd.DataFrame
    pars: pd.DataFrame
    data: pd.DataFrame
    N2_BG: pd.DataFrame


EvRHE = "E_vs_RHE"

if 0:
    nn = N2_data(
        "//mnt/DATA/APPS_SOFT/VENVS/repos/elchempy/data/raw/06.03.2018_DW28_HPRR_0.1MHClO4_RRDE22960/N2_20cls_300_100_10_DW28_298.par"
    )


class N2_Data(ElChemData):
    def __init__(self, filepath: [Path, str], **kwargs):
        self.filepath = filepath
        self.kwargs = kwargs
        super().__post_init__()

        N2_CVs = self.select_data()

        N2_results = self.analyze(N2_CVs)
        self.add_analysis_method(N2_results)

    def select_data(self):
        # FIXME Select only CV types from Data segment
        # Select the data for N2 Cyclic Voltammograms
        N2_CVs = self.data.loc[self.data.ActionId == 38]
        N2_CVs = N2_CVs.dropna(subset=["scanrate"]).loc[N2_CVs.scanrate_calc != 0]
        return N2_CVs

    def analyze(self, N2_CVs):
        """
        Performs the steps in the N2 analysis and add the
        results to the ElchemData instance.

        Parameters
        ----------
        ecdata : ElchemData
            contains the raw data.

        Returns
        -------
        ecdata : ElchemData
            contains the results in an added method as attribute "N2".

        """

        # ElChemData
        Cdl_pars, Cdl_data = pd.DataFrame(), pd.DataFrame()
        # scanrates = N2_CVs.scanrate.unique()
        if N2_CVs.scanrate.nunique() > 2:
            # if multiple scanrates are present than calculations can start
            Cdl_pars, Cdl_data = CDL(N2_CVs, EvRHE=EvRHE)
        # Check if possible background (BG) scan is in the data
        BG_present = False
        N2_BG = pd.DataFrame()
        if N2_CVs.scanrate.min() < 0.015:
            # check presence of slow scanrates in data
            BG_present = True
        if BG_present:
            N2_BG = get_N2_background_data()

        N2_results = N2_Results(N2_CVs, Cdl_pars, Cdl_data, N2_BG)

        return N2_results


def _old_N2_analysis(ecdata: ElChemData):
    """
    Performs the steps in the N2 analysis and add the
    results to the ElchemData instance.

    Parameters
    ----------
    ecdata : ElchemData
        contains the raw data.

    Returns
    -------
    ecdata : ElchemData
        contains the results in an added method as attribute "N2".
    """
    # ElChemData
    # FIXME Select only CV types from Data segment
    # Select the data for N2 Cyclic Voltammograms
    N2_CVs = ecdata.data.loc[ecdata.data.ActionId == 38]

    N2_CVs = N2_CVs.dropna(subset=["scanrate"]).loc[N2_CVs.scanrate_calc != 0]

    Cdl_pars, Cdl_data = pd.DataFrame(), pd.DataFrame()
    # scanrates = N2_CVs.scanrate.unique()
    if N2_CVs.scanrate.nunique() > 2:
        # if multiple scanrates are present than calculations can start
        Cdl_pars, Cdl_data = CDL(N2_CVs, EvRHE=EvRHE)

    # Check if possible background (BG) scan is in the data
    BG_present = False
    N2_BG = pd.DataFrame()
    if N2_CVs.scanrate.min() < 0.015:
        # check presence of slow scanrates in data
        BG_present = True
    if BG_present:
        N2_BG = get_N2_background_data()

    N2_ecdata = N2_method(N2_CVs, Cdl_pars, Cdl_data, N2_BG)

    ecdata.add_analysis_method(N2_ecdata)

    return ecdata


def CDL(N2_CVs: pd.DataFrame, EvRHE: str):

    scanrates = N2_CVs.scanrate.unique()

    srgrpby = N2_CVs.groupby("scanrate")

    # _V08 = N2_CVs.loc[np.isclose(0.8, N2_CVs[EvRHE], atol=0.01)]
    grp = N2_CVs.groupby(["Segment #", "SweepType", EvRHE])

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
    return Cdl_pars, Cdl_data

    # TODO SEND TO PLOTS AND EXPORT


def get_Cdl_pars_per_E(Cdl_data, E_linspace=np.linspace(0, 1, 21)):

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
    Cdl_dataprep, xcol="scanrate", ycol="j_A_cm2", grpbykeys=[EvRHE, "SweepType"]
):
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
        if cdl_pars_lst and cdl_data_lst:
            cdl_pars = pd.concat(cdl_pars_lst)
            cdl_data = pd.concat(cdl_data_lst)
        else:
            logger.warning("make_linregress_Cdl_Ztest error, empty pars and data")
    else:
        logger.warning("make_linregress_Cdl_Ztest error, missing columns")

    return cdl_pars, cdl_data


def linregress_residual(x, y):
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
    # return gr


def check_for_linear_baseline_correction_of_Cdl_values(Cdl_pars):
    """takes a linear fit over the Cdl values and substracts it as baseline correction"""
    _swplst = []
    for swpnm, swpgrp in Cdl_pars.groupby("SweepType"):

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


def prepare_Cdl_frame(
    N2_scans,
    current_density_key="j_A_cm2",
    potential_key=EvRHE,
    E_linspace=np.linspace(0.1, 1, 50),
):
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
    Cdl_data : TYPE
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


def check_if_last_seg_isnormal(N2sr):
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
