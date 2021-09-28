"""
Created on Fri Jul 16 10:20:30 2021

@author: DW
"""

import numpy as np
import pandas as pd

import elchempy

from elchempy.experiments.dataloader._dev_fetcher import get_files, _dev_test_read


#        grB = N2_CVs.groupby(by=['Gas','Type','EXP'])
#        for scan in grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')):
#            print(scan)
#                grA.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','0.1'))
#        grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')).to_csv(N2_dest_dir.joinpath('%s.csv' %N2_fn))
def _test_data():
    pass


def _test_N2_runner():
    datacollection = _dev_test_read(get_files("N2"))
    for ecdata in datacollection:
        N2CV = N2_selection(ecdata)


def N2_selection(ecdata):
    # FIXME Select only CV types from Data segment
    N2_CVs = ecdata.loc[ecdata.ActionId == 38]

    N2_CVs = N2_CVs.loc[N2_CVs.scanrate_calc != 0]

    scanrates = N2_CVs.scanrate.unique()

    # Check if possible background (BG) scan is in the data
    BG_present = False
    if any(sr < 0.02 for sr in scanrates):
        BG_present = True

    if len(scanrates) > 1:
        CDL(N2_CVs, scanrates)


def CDL(N2_CVs, scanrates):

    EvRHE = "E_AppV_RHE"

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


#%%
def calculate_Cdl():
    ###### ====== Analyze the Capacity (Cdl) of the N2 scan ========== #######
    # Cdl_out = pd.DataFrame({EvRHE: []})
    lst_dct, index_out = [], []
    #        N2_Cdl_ovv = pd.DataFrame(data=ScanRates, index=ScanRates, columns=['ScanRate'])
    #        index_info = {'PAR_file': N2_FileName
    # _Cdl_scans_plot = []
    if len(ScanRates) > 1:
        for SR in ScanRates:
            SR_mVs = int(SR * 1000)
            #               TODO export_cols = [EvRHE, 'j A/cm2','Scan Rate (V/s)']
            srA = grA.get_group(
                ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", SR)
            )
            nuniq_segs = srA["Segment #"].nunique()
            uniq_segs = srA["Segment #"].unique()
            N2sr_fn_out_base = N2_dest_dir.joinpath("CV_%s_%s.xlsx" % (N2_fn, SR_mVs))
            N2sr_DF_out = srA.loc[srA["Segment #"] == uniq_segs[-1], :]

            N2sr_DF_out = check_if_last_seg_isnormal(nuniq_segs, uniq_segs, srA)

            N2sr_DF_out = N2sr_DF_out.assign(
                **{"PAR_file": N2_FileName, "ScanRate_mVs": SR_mVs}
            )

            N2sr_fn_out = FileOperations.CompareHashDFexport(
                N2sr_DF_out, N2sr_fn_out_base, silent=True
            )
            N2sr_DF_out = N2sr_DF_out.assign(**{"N2_CV_datafile": N2sr_fn_out})

            #                , 'ScanRate': SR, 'DestFile': N2sr_fn_out, 'Segment': last_seg,
            #                              'Type_output': 'N2_CV'}
            #                index_out.append(index_info)
            #                    N2sr_DF_out.to_excel(N2sr_fn_out)

        def first_scan():
            if nuniq_segs > 1:
                #                           srA.loc[srA['Segment #'] == srA['Segment #'].unique()[0],[EvRHE,'j A/cm2']].to_excel(N2_dest_dir.joinpath('%s_%s_first.xlsx' %(N2_fn,int(SR*1000))))
                # Take first scan and do somethings
                N2_sr_DF = srA.loc[srA["Segment #"] == uniq_segs[0], :]
                N2_sr_fn_base = N2_dest_dir.joinpath(
                    "CV_%s_%s_first.xlsx" % (N2_fn, SR_mVs)
                )
                N2_sr_DF = N2_sr_DF.assign(
                    **{"PAR_file": N2_FileName, "ScanRate_mVs": SR_mVs}
                )
                N2_sr_fn = FileOperations.CompareHashDFexport(N2_sr_DF, N2_sr_fn_base)
            #                    index_info = {'PAR_file': N2_FileName, 'ScanRate': SR, 'DestFile': N2_sr_fn, 'Segment': first_seg,
            #                                  'Type_output': 'N2_CV'}
            #                    index_out.append(index_info)
            Cdl_scan = N2sr_DF_out
            # _Cdl_scans_plot.append(Cdl_scan)
            #                srA.loc[srA['Segment #'] == srA['Segment #'].unique()[-1]]

        #                    N2_Cdl_ovv.assign(**{'ScanRate': SR, f'j_{sweep}' : j_Cdl})
        # Cdl_scans = pd.concat(_Cdl_scans_plot)

    def exportingstuff():
        N2_plot_Cdl_scans_scanrate(
            Cdl_scans, N2_dest_dir.joinpath(f"N2_ScanRates_{N2_fn}.png")
        )

        Cdl_data = pd.DataFrame(lst_dct)
        Cdl_fn_out_base = N2_dest_dir.joinpath("Cdl_data_%s.xlsx" % N2_fn)
        Cdl_data = Cdl_data.assign(
            **{"PAR_file": N2_FileName, "DestFile_N2_BG": N2_act_BG}
        )
        Cdl_fn_out = FileOperations.CompareHashDFexport(Cdl_data, Cdl_fn_out_base)
        #            index_info = {'PAR_file': N2_FileName, 'DestFile': Cdl_fn_out, 'Type_output': 'N2_Cdl_data'}
        #            index_out.append(index_info)
        #                Cdl_pars.to_csv(Cdl_fn_out)

    def cdl_lin():
        N2_CV_datafilenames = ", ".join(
            [Path(i).name for i in Cdl_data["N2_CV_datafile"].unique()]
        )
        Cdl_cath = Cdl_data.loc[:, ["j_cathodic", "E_V", "ScanRate"]].dropna(
            axis=0, how="any"
        )
        Cdl_anod = Cdl_data.loc[:, ["j_anodic", "E_V", "ScanRate"]].dropna(
            axis=0, how="any"
        )
        fitl = []

        _meta_info = {
            "SampleID": SampleID,
            "Filename": N2_FileName,
            "Cdl_datafile": Cdl_fn_out,
            "N2_CV_datafilenames": N2_CV_datafilenames,
        }

        #           ====# UGLY STUFF FOR EMPTY Anodic Cdl scan ===
        if not Cdl_an_slice.empty:
            Cdl_an_corr = linregress(Cdl_an_slice[EvRHE], Cdl_an_slice["Cdl"])
            Cdl_an_slice = Cdl_an_slice.assign(
                **{
                    "Cdl_corr": Cdl_an_slice["Cdl"]
                    - (Cdl_an_slice[EvRHE].values * Cdl_an_corr[0] + Cdl_an_corr[1])
                    + Cdl_an_slice["Cdl"].mean()
                }
            )
        elif Cdl_an_slice.empty:
            Cdl_an_corr = linregress(Cdl_cath_slice[EvRHE], Cdl_cath_slice["Cdl"])
            Cdl_an_slice = Cdl_cath_slice
            Cdl_an_slice = Cdl_an_slice.assign(**{"Cdl_corr": 0})
        #                               columns=['SampleID','Sweep',EvRHE,'Cdl','Cdl_std','Cdl_R'])
        Cdlprs = {
            "Cdl_cath_max": Cdl_cath_slice.loc[Cdl_cath_slice["Cdl_corr"].idxmax(), :],
            "Cdl_cath_mean": Cdl_cath_slice.loc[:, ["Cdl_corr", "Cdl", "Cdl_R"]].mean(),
            "Cdl_an_max": Cdl_an_slice.loc[Cdl_an_slice["Cdl_corr"].idxmax(), :],
            "Cdl_an_mean": Cdl_an_slice.loc[:, ["Cdl_corr", "Cdl", "Cdl_R"]].mean(),
            "Analysis_date": datetime.now(),
            "SampleID": SampleID,
            "PAR_file": N2_FileName,
        }
        Cdl_Eout = []
        for Et in np.linspace(0, 1, 21):
            try:
                Cdl_Eout.append(
                    Cdl_cath_slice.loc[
                        (np.isclose(Cdl_cath_slice[EvRHE], Et, atol=50e-3) == True),
                        :,
                    ].head(1)
                )
                Cdl_Eout.append(
                    Cdl_an_slice.loc[
                        (np.isclose(Cdl_an_slice[EvRHE], Et, atol=50e-3) == True), :
                    ].head(1)
                )
            except Exception as e:
                logger.error("N2 Scan_Cdl Error at {0}: {1}".format(Et, e))
        Cdl_PARS = pd.concat([i for i in Cdl_Eout], axis=0)
        Cdl_PARS_ovv = pd.merge(Cdl_PARS, N2_ovv_file)

        Cdl_fn_pars_base = N2_dest_dir.joinpath("Cdl_pars_%s.xlsx" % N2_fn)
        Cdl_fn_pars_out = FileOperations.CompareHashDFexport(
            Cdl_PARS_ovv, Cdl_fn_pars_base
        )
        #            index_info = {'PAR_file': N2_FileName, 'DestFile': Cdl_fn_pars_out, 'Type_output': 'N2_Cdl_pars'}
        #            index_out.append(index_info)
        N2_plot_Cdl_sweeptype_scatter(
            SampleID,
            Cdl_cath_slice.ScanRates.unique()[0],
            Cdl_fit,
            Cdl_cath_slice,
            Cdl_an_slice,
            N2_dest_dir,
            N2_fn,
        )

    # elif len(ScanRates) <= 1:
    # logger.warning("N2_scans few ScanRates {0},{1}".format(ScanRates, N2_fn))

    # FIRST LAST N2 SCANS#
    preN2_scan = grA.get_group(
        ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.min())
    )
    cls_fastN2_scan = grA.get_group(
        ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.max())
    )["Segment #"].unique()
    fastN2_scan = grA.get_group(
        ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.max())
    )
    gr_cls_fastN2_scan = fastN2_scan.groupby(by="Segment #")
    first_fastN2_scan, last_fastN2_scan = gr_cls_fastN2_scan.get_group(
        cls_fastN2_scan[0]
    ), gr_cls_fastN2_scan.get_group(
        cls_fastN2_scan[int(len(cls_fastN2_scan) * 0.85) + 1]
    )


def GET_BACKGROUND_N2_FOR_ORR():
    #            N2_scan = grA.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','0.1'))
    ### === Prepare the background N2 scan (10 mV/s) for ORR with exactly 2000 rows === ###
    preN2_scan = preN2_scan.drop_duplicates()

    if float(ScanRates.min()) > 0.01:
        N2_factor = float(ScanRates.min()) / 0.01
        N2_scan = preN2_scan.assign(
            **{"j A/cm2": preN2_scan.loc[:, "j A/cm2"] / N2_factor}
        )
        logger.warning(
            f"!! N2 background wrong scan rate: {ScanRates.min()}, so j divided by {N2_factor} for {N2_fn}"
        )
    #                N2_scan.loc[:,'j A/cm2'] = N2_scan['j A/cm2']/N2_factor
    #                pd.DataFrame([ScanRates.min()])
    try:
        lenpreN2 = len(preN2_scan)
        N2_scan = preN2_scan
        if lenpreN2 > 2001:
            logger.warning(
                "!! N2 scan more than 2000 data points.. {0} !! len({1})".format(
                    N2_fn, lenpreN2
                )
            )
            #            N2_act_BG_over2000 = Path(N2_dest_dir.parent).joinpath(N2_fn+'_test2000.xlsx')# test
            #            preN2_scan.to_excel(N2_act_BG_over2000)
            #            logger.warning('!! N2 scan more than 2000 data points.. file saved {0} !! len({1})'.format(N2_fn,lenpreN2))
            preN2_scan = preN2_scan.loc[
                preN2_scan.PAR_file == preN2_scan["PAR_file"].unique()[-1], :
            ]
            if len(preN2_scan["Segment #"].unique()) > 1:
                for i in preN2_scan["Segment #"].unique():
                    if len(preN2_scan.loc[preN2_scan["Segment #"] == i]) == 2000:
                        N2_scan = preN2_scan.loc[preN2_scan["Segment #"] == i]
            elif len(preN2_scan["Segment #"].unique()) == 1:
                N2_scan = preN2_scan.drop_duplicates()
        elif len(N2_scan) != 2000:
            logger.warning(
                "!! N2 scan not 2000 data points.. {0} !! len({1})".format(
                    N2_fn, len(N2_scan)
                )
            )
            #            print('!! N2 scan less than 2000 data points.. !! %s' %lenpreN2)
            N2_scan = preN2_scan.loc[
                preN2_scan.PAR_file == preN2_scan["PAR_file"].unique()[0], :
            ]
            if len(preN2_scan["Segment #"].unique()) == 1:
                for i in preN2_scan["Segment #"].unique():
                    if len(preN2_scan.loc[preN2_scan["Segment #"] == i]) == 2000:
                        N2_scan = preN2_scan.loc[preN2_scan["Segment #"] == i]

        elif lenpreN2 < 2000:
            logger.warning(
                "!! N2 scan less than 2000 data points..  {0} !! len({1})".format(
                    N2_fn, lenpreN2
                )
            )
            N2_scan = preN2_scan.loc[
                preN2_scan.PAR_file == preN2_scan["PAR_file"].unique()[0], :
            ]
            if len(preN2_scan["Segment #"].unique()) == 1:
                for i in preN2_scan["Segment #"].unique():
                    if len(preN2_scan.loc[preN2_scan["Segment #"] == i]) == 2000:
                        N2_scan = preN2_scan.loc[preN2_scan["Segment #"] == i]
        #                            if len(preN2_scan.loc[preN2_scan['Segment #'] == i]) != 2000:
        #                            N2_scan = pd.DataFrame([])
        #                            print('!! Wrong size for N2 background!!')
        else:
            logger.info(
                "!! N2 scan has 2000 data points..success  {0} !! len({1})".format(
                    N2_fn, lenpreN2
                )
            )
            N2_scan = preN2_scan
    except Exception as e:
        logger.error("N2 scan length problems %s", e)
        N2_scan = preN2_scan
    if N2_scan.empty:
        logger.warning("!! N2 background is Empty!! {0}".format(N2_FileName))
    ### === Plot the N2 scans 1st, last,


def _get_N2_analysis_results(N2_CVs: pd.DataFrame) -> N2_Results:
    """
    Performs the steps in the N2 analysis and returns
    a namedtuple with results

    Parameters
    ----------
    N2_CVs : pd.DataFrame
        contains the raw data.

    Returns
    -------
    N2_results : namedtuple
        contains the results in an added method as attribute "N2".
    """

    class N2_Results(NamedTuple):
        raw_data: pd.DataFrame
        pars: pd.DataFrame
        data: pd.DataFrame
        N2_BG: pd.DataFrame

    if N2_CVs.empty:
        return None
    N2_Cdl_calculations
    # ElChemData
    Cdl_pars, Cdl_data = pd.DataFrame(), pd.DataFrame()
    # scanrates = N2_CVs.scanrate.unique()
    if N2_CVs.scanrate.nunique() > 2:
        # if multiple scanrates are present than calculations can start
        Cdl_pars, Cdl_data = N2_Cdl_calculations(N2_CVs, EvRHE=EvRHE)
    # Check if possible background (BG) scan is in the data
    # BG_present = False
    # N2_BG = pd.DataFrame()
    # if N2_CVs.scanrate.min() < 0.015:
    #     # check presence of slow scanrates in data
    #     BG_present = True
    if contains_background_scan(N2_CVs):
        N2_BG = get_N2_background_data(N2_CVs)

    N2_results = N2_Results(N2_CVs, Cdl_pars, Cdl_data, N2_BG)

    return N2_results
