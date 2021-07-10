import sys
from pathlib import Path
from collections import namedtuple
from datetime import datetime
import numpy as np

from scipy.stats import linregress

# from scipy.interpolate import UnivariateSpline
from scipy import constants
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import os
import multiprocessing
from functools import partial
from itertools import repeat
import pandas as pd

print("File", __file__, "\nName;", __name__)
if __name__ == "__main__":
    pass

import logging

logger = logging.getLogger(__name__)

EvRHE = "E_AppV_RHE"


def HER_scan(fit_run_arg, electrode_properties={}, **HER_kwargs):

    #    All_HER, HER_ovv_file, dest_dir)
    PAR_file, file_ovv, ovv = fit_run_arg
    file_ovv_row = file_ovv.iloc[0]
    EvRHE = "E_AppV_RHE"
    # Electrode_Type , CollEff, SA_disk, SA = electrode_properties.values
    # FIX ME WE_SA_collection_eff('PINE-ring').values()
    mA = 1000
    #    if All_HER.empty:
    #        print('!! Critical HER empty: %s!!' % dest_dir)
    #    if HER_ovv_file.empty:
    #        #           ovv[~ovv['SampleID'].str.contains('Pt_ring')].loc[:,['PAR_exp' == 'N2']].empty:
    #        logger.warning('!! Critical error HER empty: {0}!!'.format(dest_dir))
    HER_dest_dir = file_ovv_row[HER_kwargs.get("dest_dir_col") + "_ecexp"]
    #    file_ovv_row[HER_kwargs.get('dest_dir_col')]
    SampleID = file_ovv_row["SampleID"]

    HER_CVs, HER_action = create_CVs(file_ovv)
    HER_CV = HER_CVs.loc[HER_CVs.Type_action.str.contains("Cyclic Voltammetry")].query(
        'ScanRate_calc < 0.2 & SampleID != "Pt_ring"'
    )

    #    All_OER = All_OER.assign(**{'jmAcm-2' :  All_OER['j A/cm2']*1000, 'Abs_jmAcm-2' : np.abs(All_OER['j A/cm2']*1000),
    #                          'log_Abs_jmAcm-2' : np.log10(np.abs(All_OER['j A/cm2']*1000)),'RPM' : 1500 })
    #    make_sure_path_exists(HPRR_dest_dir)
    HER_CV = HER_CV.assign(
        **{
            "jmAcm-2": HER_CV["j A/cm2"] * 1000,
            "Abs_jmAcm-2": np.abs(HER_CV["j A/cm2"] * 1000),
            "log_Abs_jmAcm-2": np.log10(np.abs(HER_CV["j A/cm2"] * 1000)),
            "RPM_1500": 1500,
        }
    )

    #    HPRR_fn = Path(HPRR_ovv['File'].unique()[0]).stem
    HER_PAR_fn = PAR_file
    HER_fn = HER_PAR_fn.stem

    if HER_CV.empty:
        #           ovv[~ovv['SampleID'].str.contains('Pt_ring')].loc[:,['PAR_exp' == 'N2']].empty:
        logger.warning("!! Critical error HER CVs empty: {0}!!".format(HER_dest_dir))

    try:
        _scan_grp_cols = [
            "Gas",
            "Type_action",
            "EXP",
            "Scanrate",
            "RPM_DAC",
            "Segment #",
        ]
        grA = HER_CV.groupby(by=_scan_grp_cols)
        #        grB = HPRR_CV.groupby(by=['Gas','Type','EXP'])
        #        for scan in grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')):
        #            print(scan)
        #                grA.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','0.1'))
        #        grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','HPRR')).to_csv(HPRR_dest_dir.joinpath('%s.csv' %HPRR_fn))
        #        hp_data = grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','HPRR'))
        out, HER_pars_out, HER_out_lst = [], pd.DataFrame(), []
        for nm, gr in grA:
            for swnm, sweep in gr.groupby(by="Sweep_Type"):
                if swnm == "NA":
                    continue

                _sweep_meta_info = {
                    i: sweep[i].unique()[0]
                    for i in sweep.columns
                    if sweep[i].nunique() == 1
                }
                swp_target_file = HER_dest_dir.joinpath(
                    "{0}_{1}_{2}_{3}.pkl".format(swnm, nm[4], nm[5], HER_fn)
                )
                _sweep_meta_info.update({"HER_swp_targetfile": swp_target_file})
                #                try:
                #                    old_target = OER_dest_dir.joinpath('HPRR_Tafel_{0}_{1}.xlsx'.format(swnm,nm[4]))
                #                    if old_target.is_file():
                #                        old_target.unlink()
                #                        logger.warning('OER output deleted old target: {0}'.format(old_target))
                #                except:
                #                    logger.warning('HPRR output delete old target fail: {0}'.format(old_target))

                swp = sweep.loc[
                    :,
                    [
                        EvRHE,
                        "jmAcm-2",
                        "Abs_jmAcm-2",
                        "log_Abs_jmAcm-2",
                        "Sweep_Type",
                        "RPM_DAC",
                        "PAR_file",
                    ]
                    + _scan_grp_cols,
                ]
                j_use = "jmAcm-2_fltr"
                #                .rolling(window=5).mean()
                swp[j_use] = savgol_filter(swp["jmAcm-2"], 21, 3)
                swp = swp.assign(
                    **{
                        "Abs_jmAcm-2_fltr": np.abs(swp["jmAcm-2_fltr"]),
                        "log_Abs_jmAcm-2_fltr": np.log10(np.abs(swp["jmAcm-2_fltr"])),
                        "j/E": swp[j_use] / swp[EvRHE],
                        "dJ": swp[j_use].diff(),
                        "d/d2": swp[j_use].diff().diff(),
                        "dE": swp[EvRHE].diff(),
                        "E_overp": -1 * swp[EvRHE],
                    }
                )
                swp["dj/dE"] = swp["dJ"] / swp["dE"]
                #                swp.plot(x=EvRHE,y=j_use,logy=0)
                #                OER_slices.plot(x='E_overp',y=j_use,logy=0)

                ###### ======Analyzing HPRR CV and extracting kinetic parameters ========== #######
                #                HPOR = swp.loc[(np.isclose(swp[EvRHE],swp[EvRHE].max()-0.002,atol=0.010))][j_use].mean()
                HER_E_pars_dict = {}
                #                HER_slices = pd.concat([i for i in [swp.loc[np.isclose(swp[EvRHE],i,atol=0.005)].head(1) for i in np.arange(1.3,1.90,0.05)]])
                #                    OER_155 = swp.loc[[np.isclose(swp[EvRHE],i,atol=0.010).head(1) for i in np.arange(0.15,0.185,0.005)] ]

                J_options = [-1, -2, -3, -4, -5, -7, -10]
                TF_out = []
                for j_take_opt in J_options:
                    HER_E_pars_dict = {}
                    HER_onset = (
                        swp.loc[
                            (
                                np.isclose(swp[j_use], j_take_opt, rtol=0.01)
                                & (swp[EvRHE] < 0)
                            ),
                            :,
                        ]
                        .sort_values(by=EvRHE)
                        .head(1)
                    )
                    if HER_onset.empty:
                        rtol_set = 0.01
                        while HER_onset.empty:
                            HER_onset = (
                                swp.loc[
                                    (
                                        np.isclose(
                                            swp[j_use], -j_take_opt, rtol=rtol_set
                                        )
                                        & (swp[EvRHE] < 0)
                                    ),
                                    :,
                                ]
                                .sort_values(by=EvRHE)
                                .head(1)
                            )
                            rtol_set += 0.01

                    TF_onset_lowerE, TF_onset_upperE = (
                        HER_onset[EvRHE].values[0] + 0.0,
                        HER_onset[EvRHE].values[0] + 0.1,
                    )

                    swp_onset_TF = swp.loc[
                        (swp[EvRHE] <= TF_onset_upperE)
                        & (swp[EvRHE] >= TF_onset_lowerE),
                        :,
                    ]
                    #                swp_onset_TF.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                    #                swp.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                    TF_onset_fit = linregress(
                        swp_onset_TF["log_Abs_jmAcm-2_fltr"].values,
                        swp_onset_TF[EvRHE].values,
                    )
                    HER_E_pars_dict.update(
                        {
                            "HER_type": "j_slice_onset",
                            "HER_at_J_slice": j_take_opt,
                            "HER_E_onset": TF_onset_lowerE,
                            "HER_E_onset_upper": TF_onset_upperE,
                            "HER_Tafel_slope": TF_onset_fit[0] * -1000,
                        }
                    )
                    HER_E_pars_dict.update(
                        dict(
                            zip(
                                ["HER_" + i for i in TF_onset_fit._fields], TF_onset_fit
                            )
                        )
                    )

                    TF_onset_out = [
                        "E_onset",
                        TF_onset_lowerE,
                        TF_onset_upperE,
                        TF_onset_fit[0],
                        TF_onset_fit[1],
                        TF_onset_fit[2],
                        TF_onset_fit[0] * -1000,
                    ]
                    #                TF_onset_fit[0]*1000
                    TF_out.append(HER_E_pars_dict)

                #                Tafel_ovv = pd.DataFrame([])
                #                TFlst, TafelDir = [], HER_dest_dir.joinpath('TAFEL')
                #                TafelDir.mkdir(parents=True, exist_ok=True)
                E_lower, E_upper = swp[EvRHE].min(), swp[EvRHE].max()
                for Ev in np.arange(-0.8, 0.1, 0.04):
                    HER_E_pars_dict = {}
                    swp_TF_slice = swp.loc[
                        (swp[EvRHE] >= Ev) & (swp[EvRHE] <= Ev + 0.03), :
                    ]
                    j_lower, j_upper = (
                        swp_TF_slice.loc[swp_TF_slice[EvRHE].idxmin(), j_use],
                        swp_TF_slice.loc[swp_TF_slice[EvRHE].idxmax(), j_use],
                    )
                    #                    swp_TF_slice.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                    #                    OER_slices.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                    if not swp_TF_slice.empty:
                        TF_slice_fit = linregress(
                            swp_TF_slice["log_Abs_jmAcm-2_fltr"].values,
                            swp_TF_slice[EvRHE].values,
                        )
                        if np.abs(TF_slice_fit.rvalue) > 0.90:
                            HER_E_pars_dict.update(
                                {
                                    "HER_type": "E_slice",
                                    "HER_at_E_slice": Ev,
                                    "HER_J_lower": j_lower,
                                    "HER_J_upper": j_upper,
                                    "HER_Tafel_slope": TF_slice_fit[0] * -1000,
                                }
                            )
                            HER_E_pars_dict.update(
                                dict(
                                    zip(
                                        ["HER_" + i for i in TF_slice_fit._fields],
                                        TF_slice_fit,
                                    )
                                )
                            )
                            TF_out.append(HER_E_pars_dict)

                #                            TF_slc_out = ['E_slice', Ev, Ev + 0.03, TF_slice_fit[0], TF_slice_fit[1], TF_slice_fit[2],
                #                                          TF_slice_fit[0] * 1000, j_lower, j_upper]

                TF_pars_out = pd.DataFrame(TF_out)
                _HER_sweep_datafile = FileOperations.CompareHashDFexport(
                    swp, swp_target_file
                )
                TF_pars_out = TF_pars_out.assign(
                    **{**{"HER_DataFile": _HER_sweep_datafile}, **_sweep_meta_info}
                )

                try:
                    HER_plot_sweep(swp, j_use, swp_target_file, TF_pars_out)
                except Exception as e:
                    logger.error(f"HER plotting error: {e}\n for {swp_target_file}")

                #                TF_pars_out = pd.DataFrame(data=TF_out,
                #                                           columns=['E_type', EvRHE, EvRHE + '_upper', 'TF_a', 'TF_b', 'TF_fit_error',
                #                                                    'TafelSlope', 'j_lower', 'j_upper']).sort_values(EvRHE)

                #                TF_pars_out = TF_pars_out.assign(
                #                    **dict(zip(['Gas', 'Type_action', 'EXP', 'Scanrate', 'RPM_DAC', 'Segment #'], nm)))

                HER_out_lst.append(TF_pars_out)
        ###### ======Saving all HPRR CV kinetic parameters to file and index ========== #######
        HER_pars_out = pd.concat(
            [i for i in HER_out_lst], sort=False, ignore_index=True
        )
        HER_pars_base = HER_dest_dir.joinpath(HER_fn + "_pars.xlsx")
        HER_pars_target = FileOperations.CompareHashDFexport(
            HER_pars_out, HER_pars_base
        )
        logger.info(f"HER pars success!!: {HER_pars_target}")
    #        index_info_HER_TAFEL = pd.DataFrame(
    #            {'PAR_file': HER_PAR_fn, 'DestFile': HER_pars_target, 'Type_output': 'HER_Jkin_Tafel', 'Type_exp': 'HER'},
    #            index=[0])

    except Exception as e:
        #        print('No successfull HER: {0}'.format(e))
        logger.error("No successfull HER: {0}".format(e))
        #        index_info_HER_TAFEL = pd.DataFrame(
        #            {'PAR_file': HER_PAR_fn, 'DestFile': HER_ovv_file.index[0], 'Type_output': e, 'Type_exp': 'HER'}, index=[0])
        HER_pars_out = pd.DataFrame([])


#    return index_info_HER_TAFEL


def HER_plot_sweep(swp, j_use, swp_target_file, TF_pars_out):
    global EvRHE
    fig, ax = plt.subplots(figsize=(12, 8))
    #                ax = host.twinx()
    j_ax = ax.twinx()
    #                TF_pars_out.plot(x=EvRHE, y='TafelSlope', kind='scatter', c='k', s=30, title=swp_target_file.stem,
    #                                 ax=ax, label='TafelSlopes')
    swp.plot(x=EvRHE, y=j_use, kind="line", title=swp_target_file.stem, ax=j_ax, logy=0)
    TF_pars_out.query('HER_type == "j_slice_onset"').plot(
        x="HER_E_onset",
        y="HER_Tafel_slope",
        kind="scatter",
        s=80,
        c="red",
        title=swp_target_file.stem,
        ax=ax,
    )
    ax.legend(loc="upper left")
    j_ax.legend(loc="lower right")
    plt.savefig(swp_target_file.with_suffix(".png"), bbox_inches="tight", dpi=100)
    plt.close()


def HER_calc(HER_ovv, HER_out_fn, PathDB, plot_HER=True):
    HER_Pars = []
    fig, ax = plt.subplots(figsize=(10, 10))
    for a, gr in HER_ovv.groupby(by="PAR_file"):
        grCV, grInfo = create_CVs(gr)
        grCV = grCV.loc[grCV.Type == "Cyclic Voltammetry"]
        grOVV = HER_ovv.loc[HER_ovv["PAR_file"] == str(a)]
        if not grCV.empty or not grOVV.empty:
            SegUniq = grCV["Segment #"].unique()
            for Cycle, Seg in enumerate(SegUniq):
                Segr = grCV.loc[grCV["Segment #"] == Seg]

                E_HER_kin = -0.6
                HER_jkin = Segr.loc[
                    np.isclose(Segr[EvRHE], E_HER_kin, rtol=0.001), "jmAcm-2"
                ].mean()
                grOVV = grOVV.assign(
                    **{
                        "E_HER_kin": E_HER_kin,
                        "HER_jkin": HER_jkin,
                        "HER_cycle": Cycle + 1,
                    }
                )
                if Seg == SegUniq[-1]:
                    Segr.plot(
                        x=EvRHE,
                        y="jmAcm-2",
                        ax=ax,
                        label=" %s Cycle: %s (%.2f )"
                        % (grOVV.basename.values[0], Cycle + 1, HER_jkin),
                    )
                    HER_Pars.append([grOVV])
    HERpOut = pd.concat([i[0] for i in HER_Pars])
    FolderOps.FileOperations.CompareHashDFexport(HERpOut, HER_out_fn)
    if plot_HER:
        print("HER output to: %s" % HER_out_fn)
        plt.savefig(HER_out_fn.with_suffix(".png"))


#        HERpOut.to_excel(HER_out_fn)
