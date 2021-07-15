from pathlib import Path
import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import scipy
from scipy.stats import linregress
import matplotlib.pyplot as plt

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.ExtraInfo import EC_Properties
from .. import EvRHE


# All_OER = All_ovv.query('EXP == "OER"')
# All_OER, OER_ovv_file, dest_dir = (Samples_ovv_cv,OER_ovv_file, Path(OER_ovv_file.Dest_dir.iloc[0]))
#    OER_file_index =  OER_scan(Samples_ovv_cv,OER_ovv_file, Path(OER_ovv_file.Dest_dir.iloc[0]))
def OER_scan(All_OER, OER_ovv_file, dest_dir):
    CollEff, SA_disk, SA = EC_Properties.WE_SA_collection_eff("PINE-ring").values()
    mA = 1000
    if All_OER.empty:
        print("!! Critical OER empty: %s!!" % dest_dir)

    EvRHE = "E_AppV_RHE"
    OER_dest_dir = dest_dir.joinpath("OER_scans")
    OER_dest_dir.mkdir(parents=True, exist_ok=True)
    SampleID = All_OER["SampleID"].unique()[0]
    OER_CV = All_OER.loc[All_OER.Type_action.str.contains("Cyclic Voltammetry")].query(
        'EXP == "OER" & ScanRate_calc < 0.2 & SampleID != "Pt_ring"'
    )
    #    OER_CV = All_OER.query('EXP == "OER" & ScanRate_calc < 0.2 & SampleID != "Pt_ring" & Type_action == "Cyclic Voltammetry (Multiple Cycles)" ')
    #    All_OER = All_OER.assign(**{'jmAcm-2' :  All_OER['j A/cm2']*1000, 'Abs_jmAcm-2' : np.abs(All_OER['j A/cm2']*1000),
    #                          'log_Abs_jmAcm-2' : np.log10(np.abs(All_OER['j A/cm2']*1000)),'RPM' : 1500 })
    #    make_sure_path_exists(HPRR_dest_dir)
    OER_CV = OER_CV.assign(
        **{
            "jmAcm-2": OER_CV["j A/cm2"] * 1000,
            "Abs_jmAcm-2": np.abs(OER_CV["j A/cm2"] * 1000),
            "log_Abs_jmAcm-2": np.log10(np.abs(OER_CV["j A/cm2"] * 1000)),
            "RPM_1500": 1500,
        }
    )

    #    HPRR_fn = Path(HPRR_ovv['PAR_file'].unique()[0]).stem
    OER_PAR_fn = Path(OER_ovv_file.PAR_file.iloc[0])
    OER_fn = OER_PAR_fn.stem
    OER_out_lst = []
    if OER_ovv_file.empty:
        #           ovv[~ovv['SampleID'].str.contains('Pt_ring')].loc[:,['PAR_exp' == 'N2']].empty:
        logger.warning("!! Critical error OER empty: {0}!!".format(dest_dir))
    try:
        grA = OER_CV.groupby(
            by=["Gas", "Type_action", "EXP", "Scanrate", "RPM_DAC", "Segment #"]
        )
        #        grB = HPRR_CV.groupby(by=['Gas','Type','EXP'])
        #        for scan in grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')):
        #            print(scan)
        #                grA.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','0.1'))
        #        grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','HPRR')).to_csv(HPRR_dest_dir.joinpath('%s.csv' %HPRR_fn))
        #        hp_data = grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','HPRR'))
        out, OER_pars_out, OER_out_lst = [], pd.DataFrame(), []
        for nm, gr in grA:
            for swnm, sweep in gr.groupby(by="Sweep_Type"):
                if swnm == "NA":
                    continue

                swp_target_file = OER_dest_dir.joinpath(
                    "{0}_{1}_{2}_{3}.xlsx".format(swnm, nm[4], nm[5], OER_fn)
                )
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
                    ],
                ]
                j_use = "jmAcm-2_fltr"
                #                .rolling(window=5).mean()
                swp[j_use] = scipy.signal.savgol_filter(swp["jmAcm-2"], 21, 3)
                swp = swp.assign(
                    **{
                        "Abs_jmAcm-2_fltr": np.abs(swp["jmAcm-2_fltr"]),
                        "log_Abs_jmAcm-2_fltr": np.log10(np.abs(swp["jmAcm-2_fltr"])),
                        "j/E": swp[j_use] / swp[EvRHE],
                        "dJ": swp[j_use].diff(),
                        "d/d2": swp[j_use].diff().diff(),
                        "dE": swp[EvRHE].diff(),
                        "E_overp": swp[EvRHE] - 1.23,
                    }
                )
                swp["dj/dE"] = swp["dJ"] / swp["dE"]
                #                swp.plot(x=EvRHE,y=j_use,logy=0)
                #                OER_slices.plot(x='E_overp',y=j_use,logy=0)

                ###### ======Analyzing HPRR CV and extracting kinetic parameters ========== #######
                #                HPOR = swp.loc[(np.isclose(swp[EvRHE],swp[EvRHE].max()-0.002,atol=0.010))][j_use].mean()
                OER_E_pars_dict = {}
                #                OER_slices = pd.concat([i for i in [swp.loc[np.isclose(swp[EvRHE],i,atol=0.005)].head(1) for i in np.arange(1.3,1.90,0.05)]])
                #                    OER_155 = swp.loc[[np.isclose(swp[EvRHE],i,atol=0.010).head(1) for i in np.arange(0.15,0.185,0.005)] ]
                OER_onset = (
                    swp.loc[
                        (np.isclose(swp[j_use], 10, rtol=0.01) & (swp[EvRHE] > 1.23)), :
                    ]
                    .sort_values(by=EvRHE)
                    .head(1)
                )
                if OER_onset.empty:
                    rtol_set = 0.01
                    while OER_onset.empty or rtol_set > 10:
                        OER_onset = (
                            swp.loc[
                                (
                                    np.isclose(swp[j_use], 10, rtol=rtol_set)
                                    & (swp[EvRHE] > 1.23)
                                ),
                                :,
                            ]
                            .sort_values(by=EvRHE)
                            .head(1)
                        )
                        rtol_set += 0.01
                    if OER_onset.empty:
                        OER_onset = swp.loc[swp[j_use].idxmax()]
                        pass
                TF_out = []
                TF_onset_lowerE = OER_onset[EvRHE].values[0] + 0.0
                swp_onset_TF = swp.loc[
                    (swp[EvRHE] <= TF_onset_lowerE + 0.1)
                    & (swp[EvRHE] >= TF_onset_lowerE),
                    :,
                ]
                #                swp_onset_TF.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                #                OER_slices.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                TF_onset_fit = linregress(
                    swp_onset_TF["log_Abs_jmAcm-2_fltr"].values,
                    swp_onset_TF[EvRHE].values,
                )
                TF_onset_out = [
                    "E_onset",
                    TF_onset_lowerE,
                    TF_onset_lowerE + 0.1,
                    TF_onset_fit[0],
                    TF_onset_fit[1],
                    TF_onset_fit[2],
                    TF_onset_fit[0] * 1000,
                ]
                #                TF_onset_fit[0]*1000
                TF_out.append(TF_onset_out)
                Tafel_ovv = pd.DataFrame([])
                TFlst, TafelDir = [], OER_dest_dir.joinpath("TAFEL")
                TafelDir.mkdir(parents=True, exist_ok=True)

                for Ev in np.arange(1.3, 1.9, 0.03):
                    swp_TF_slice = swp.loc[
                        (swp[EvRHE] >= Ev) & (swp[EvRHE] <= Ev + 0.03), :
                    ]
                    #                    swp_TF_slice.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                    #                    OER_slices.plot(y=EvRHE,x='log_Abs_jmAcm-2_fltr',logy=0)
                    if not swp_TF_slice.empty:
                        j_lower, j_upper = (
                            swp_TF_slice.loc[swp_TF_slice[EvRHE].idxmin(), j_use],
                            swp_TF_slice.loc[swp_TF_slice[EvRHE].idxmax(), j_use],
                        )
                        TF_slice_fit = linregress(
                            swp_TF_slice["log_Abs_jmAcm-2_fltr"].values,
                            swp_TF_slice[EvRHE].values,
                        )
                        if np.abs(TF_slice_fit.rvalue) > 0.90:
                            TF_slc_out = [
                                "E_slice",
                                Ev,
                                Ev + 0.01,
                                TF_slice_fit[0],
                                TF_slice_fit[1],
                                TF_slice_fit[2],
                                TF_slice_fit[0] * 1000,
                                j_lower,
                                j_upper,
                            ]
                            TF_out.append(TF_slc_out)

                TF_pars_out = pd.DataFrame(
                    data=TF_out,
                    columns=[
                        "E_type",
                        EvRHE,
                        EvRHE + "_upper",
                        "TF_a",
                        "TF_b",
                        "TF_fit_error",
                        "TafelSlope",
                        "j_lower",
                        "j_upper",
                    ],
                ).sort_values(EvRHE)
                TF_pars_out = TF_pars_out.assign(
                    **{
                        "DataFile": swp_target_file,
                        "PAR_file": OER_PAR_fn,
                        "Sweep_Type": swnm,
                    }
                )
                TF_pars_out = TF_pars_out.assign(
                    **dict(
                        zip(
                            [
                                "Gas",
                                "Type_action",
                                "EXP",
                                "Scanrate",
                                "RPM_DAC",
                                "Segment #",
                            ],
                            nm,
                        )
                    )
                )

                fig, ax = plt.subplots()
                TF_pars_out.plot(
                    x=EvRHE,
                    y="TafelSlope",
                    kind="scatter",
                    title=swp_target_file.stem,
                    ax=ax,
                )
                j_ax = ax.twinx()
                swp.plot(
                    x=EvRHE,
                    y=j_use,
                    kind="line",
                    title=swp_target_file.stem,
                    ax=j_ax,
                    logy=0,
                )
                TF_pars_out.query('E_type == "E_onset"').plot(
                    x=EvRHE,
                    y="TafelSlope",
                    kind="scatter",
                    s=80,
                    c="red",
                    title=swp_target_file.stem,
                    ax=ax,
                    label="j = {0} mAcm^-2".format(),
                )
                ax.legend(loc="upper left")
                j_ax.legend(loc="lower right")
                plt.savefig(
                    TafelDir.joinpath("TAFEL_{0}.png".format(swp_target_file.stem)),
                    bbox_inches="tight",
                )
                plt.close()
            OER_out_lst.append(TF_pars_out)
        ###### ======Saving all HPRR CV kinetic parameters to file and index ========== #######
        OER_pars_out = pd.concat(
            [i for i in OER_out_lst], sort=False, ignore_index=True
        )
        OER_pars_base = OER_dest_dir.joinpath(OER_fn + "_pars.xlsx")
        # OER_pars_target = FolderOps.FileOperations.CompareHashDFexport(
        #     OER_pars_out, OER_pars_base
        # )
        # index_info_OER_TAFEL = pd.DataFrame(
        #     {
        #         "PAR_file": OER_PAR_fn,
        #         "DestFile": OER_pars_target,
        #         "Type_output": "OER_Jkin_Tafel",
        #         "Type_exp": "OER",
        #     },
        #     index=[0],
        # )

    except Exception as e:
        print("No successfull OER: {0}".format(e))
        logger.error("No successfull OER: {0}".format(e))
        index_info_OER_TAFEL = pd.DataFrame(
            {
                "PAR_file": OER_PAR_fn,
                "DestFile": OER_dest_dir.joinpath(OER_fn + "_fail.xlsx"),
                "Type_output": e,
                "Type_exp": "OER_failed",
            },
            index=[0],
        )
        OER_pars_out = pd.DataFrame([])
    return index_info_OER_TAFEL


def OER_calc(HER_ovv, HER_out_fn, PathDB, create_CVs, plot_OER=True):
    HER_Pars = []
    fig, ax = plt.subplots(figsize=(10, 10))
    for a, gr in HER_ovv.groupby(by="PAR_file"):
        grCV, grInfo = create_CVs(gr)
        if not grCV.empty or not gr.empty:
            grCV = grCV.loc[grCV.Type == "Cyclic Voltammetry"]
            grOVV = HER_ovv.loc[HER_ovv["PAR_file"] == str(a)]

            SegUniq = grCV["Segment #"].unique()
            for Cycle, Seg in enumerate(SegUniq):
                Segr = grCV.loc[grCV["Segment #"] == Seg]

                E_HER_kin = 1.6
                HER_jkin = Segr.loc[
                    np.isclose(Segr[EvRHE], E_HER_kin, rtol=0.001), "jmAcm-2"
                ].mean()
                grOVV = grOVV.assign(
                    **{
                        "E_OER_kin": E_HER_kin,
                        "OER_jkin": HER_jkin,
                        "OER_cycle": Cycle + 1,
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
    # FolderOps.FileOperations.CompareHashDFexport(HERpOut, HER_out_fn) # FIXME some EXPORT FUNC
    if plot_OER:
        print("OER output to: %s" % HER_out_fn)
        plt.savefig(HER_out_fn.with_suffix(".png"))


#        HERpOut.to_excel(HER_out_fn)
