#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:03:30 2020

@author: zmg
"""


class ExportECfromHER:
    """Exports from HER index"""

    EvRHE = "E_AppV_RHE"

    HER_cols = [
        "E_type",
        EvRHE,
        EvRHE + "_upper",
        "TF_a",
        "TF_b",
        "TF_fit_error",
        "TafelSlope",
        "j_lower",
        "j_upper",
    ]

    def __init__(self):
        pass

    def HER_plotting_parameters(series="*"):
        #%% ===== MAKE PLOTS of Parameters versus E v RHE ======
        OnlyRecentMissingOVV = run_PAR_DW.ECRunOVV(load=1).index
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        SeriesIDs = [
            SampleSelection.Series_CB_paper,
            SampleSelection.Series_Porhp_SiO2,
            SampleSelection.Series_Co_PANI,
        ]
        SeriesID_set = SeriesIDs[2]

        PostDestDirHERcom = PostDestDir.joinpath(
            "HER_Pars_Char_{0}".format(SeriesID_set["name"])
        )
        PPDHER = PostDestDirHERcom.joinpath("HER_slices")
        PPDHER.mkdir(parents=True, exist_ok=True)

        HER_pars = Load_from_Indexes.HER_pars_OVV(
            pd.DataFrame(), pd.DataFrame(), reload=False
        )  # EIS_Pars2
        #        HER_pars = HER_pars.drop_duplicates(subset=ExportECfromHER.HER_cols,keep='first')
        MatchPostAST = Load_from_Indexes.MatchPostASTs(postOVVout)
        MatchPostAST["postAST"] = [
            FileHelper.FindSampleID.determine_postAST_from_filename(i)
            for i in MatchPostAST.PAR_file.values
        ]
        HER_pars["RPM"] = 1500
        HER_pars = HER_pars.assign(
            **{
                "Loading_cm2": np.round(HER_pars["Loading_cm2"], 3),
                "postAST": [
                    FileHelper.FindSampleID.determine_postAST_from_filename(i)
                    for i in HER_pars.basename.values
                ],
            }
        )
        #        HER_pars['RPM'] = 1500
        HER_preAST = pd.merge(HER_pars, MatchPostAST, how="left", on="PAR_file")

        HER_serie = HER_pars.loc[HER_pars.SeriesID.isin(SeriesID_set["slice"])]
        OVV_CB = OnlyRecentMissingOVV.loc[
            OnlyRecentMissingOVV.PAR_exp.isin(HER_serie.Type_exp.unique())
            & OnlyRecentMissingOVV.SampleID.isin(HER_serie.SampleID.unique())
        ]

        #        eismiss_CB = OVV_CB.loc[OVV_CB.PAR_file.isin([i for i in OVV_CB.PAR_file.values if i not in v.PAR_file.values])]
        HER_E_slice = HER_serie.loc[
            (HER_pars["E_type"] == "E_slice") & (HER_serie["Segment #"] > 1)
        ]
        HER06_no = HER_E_slice.loc[
            np.isclose(HER_E_slice[EvRHE], -0.6, atol=0.015)
            & (HER_E_slice.postAST == "no")
        ].drop_duplicates(subset=ExportECfromHER.HER_cols)
        HER_E_slice.loc[np.isclose(HER_E_slice[EvRHE], -0.6, atol=0.015)].plot(
            y="j_lower", x="BET_cat_agg", kind="scatter"
        )
        HER06_no.plot(y="j_lower", x="BET_cat_agg", kind="scatter")
        HER_E_slice.loc[np.isclose(HER_E_slice[EvRHE], -0.6, atol=0.015)].plot(
            y="j_lower", x="N_content", kind="scatter"
        )
        HER06_no.plot(y="j_lower", x="N_content", kind="scatter")
        HER_E_slice.loc[np.isclose(HER_E_slice[EvRHE], -0.6, atol=0.015)].plot(
            y="TafelSlope", x="C/N_ratio", kind="scatter"
        )
        HER_E_slice.basename.unique()
        # chekcing...
        JOS3ast = HER_preAST.query('(SampleID == "JOS3")')  # & (postAST != "no")')
        plot_pd_SampleIDs(
            HER_E_slice, "j_lower", "BET_cat_agg", ddir=PPDHER, corr_val=0
        )
        ovv01 = OnlyRecentMissingOVV.query(
            "(EXP_date <= 20190126) & (EXP_date >= 20190124)"
        ).drop_duplicates()
        #        jandups = pd.concat([i[1] for i in [(n,gr[['SampleID','basename','PAR_hash','PAR_date','PAR_file','Creation_date']]) for n,gr in ovv01.groupby(['PAR_hash']) if len(gr) > 1 and gr.basename.nunique() > 1]])
        #        jos3ovv = OnlyRecentMissingOVV.loc[(OnlyRecentMissingOVV.PAR_date >= JOS3ast.PAR_date.unique().min()) & (OnlyRecentMissingOVV.PAR_date < JOS3ast.PAR_date.unique().max())]
        #        .query('(SampleID == "JOS3")')
        # /...
        EC_exp_uniq = [
            (i, HER_E_slice[i].unique()[0])
            for i in [
                c for c in SampleSelection.EC_exp_cols if c in HER_E_slice.columns
            ]
            if HER_E_slice[i].nunique() == 1
        ]
        EC_exp_non_uniq = [
            (i, HER_E_slice[i].unique())
            for i in [
                c for c in SampleSelection.EC_exp_cols if c in HER_E_slice.columns
            ]
            if HER_E_slice[i].nunique() != 1
        ]

        HERmc = HER_E_slice.set_index(
            [i[0] for i in EC_exp_uniq] + ["SampleID", "postAST"]
        )[[EvRHE, "j_upper"]]
        E_06_lst = []
        fig, ax = plt.subplots(figsize=(14, 14))
        for sID, sIDgrp in HER_E_slice.groupby(
            [i[0] for i in EC_exp_uniq] + ["SampleID"]
        ):
            for n, Hgr in sIDgrp.groupby("postAST"):
                if Hgr.preAST.nunique() == 1 and not "no-preAST" in Hgr.preAST.unique():
                    preAST_file = str(
                        FileHelper.FileOperations.find_CS_parts(
                            Hgr.preAST.unique()[0].strip("\(',)"),
                        )
                    )
                    preAST_pars_raw = HER_pars.query("(PAR_file == @preAST_file)")
                    preAST_pars = preAST_pars_raw.loc[preAST_pars_raw["Segment #"] > 1]
                    pre_mcols = [
                        i
                        for i in preAST_pars.columns
                        if Hgr[i].unique()[0] == preAST_pars[i].unique()[0]
                        and preAST_pars[i].nunique() == 1
                    ]
                    pre_mcols_nonuniq1 = [
                        i
                        for i in preAST_pars.columns
                        if Hgr[i].unique()[0] != preAST_pars[i].unique()[0]
                        and preAST_pars[i].nunique() == 1
                    ]
                    pre_mcols_nonuniq_more = [
                        i
                        for i in preAST_pars.columns
                        if Hgr[i].unique()[0] != preAST_pars[i].unique()[0]
                        and preAST_pars[i].nunique() != 1
                    ]
                    HERcols = pre_mcols_nonuniq_more[0:13]

                    #                     Eisclose = [[a for i in Hgr[EvRHE].values if np.isclose(i,a,atol=0.001)] for a in preAST_pars[EvRHE].values]
                    #                     [[i,min(Hgr[EvRHE].values, key=lambda x:abs(x-i))] for i in preAST_pars[EvRHE].values]
                    #                     list(zip(Hgr[EvRHE].values,preAST_pars[EvRHE].values))
                    #                     [i for i in Eisclose]
                    out_HER_post, out_HER_pre = Hgr[HERcols].rename(
                        columns=dict(zip(HERcols, ["post_" + i for i in HERcols]))
                    ).dropna(axis=1, how="all"), preAST_pars[HERcols].rename(
                        columns=dict(zip(preAST_pars, ["pre_" + i for i in HERcols]))
                    ).dropna(
                        axis=1, how="all"
                    )
                    pd.concat([out_HER_post, out_HER_pre], axis=1)
                    #                     post_preAST = pd.merge(Hgr[HERcols],preAST_pars[HERcols],how='left',on=EvRHE,suffixes=('_post','_pre'))
                    post_pre = pd.concat([Hgr, preAST_pars])
                    PPgr.loc[np.isclose(PPgr[EvRHE], -0.6, atol=0.015)].drop_duplicates(
                        subset=ExportECfromHER.HER_cols
                    )
                    for PPn, PPgr in post_pre.groupby("postAST"):
                        ls_set = "dashed" if PPn != "no" else "solid"
                        pAST_set = "" if PPn == "no" else PPn
                        ccode = PPgr.Colorcode.unique()
                        c_set = OriginColor.query("OriginCode == @ccode")
                        RGB = [
                            float(i) / 255
                            for i in (
                                OriginColor.query("OriginCode == @ccode")["RGB"].values[
                                    0
                                ]
                            ).split(",")
                        ]
                        label_set = "{0} ({1}) {2}".format(
                            PPgr.SampleLabel.unique()[0], sID[-1], pAST_set
                        )
                        for DT, DTgr in PPgr.groupby("PAR_date"):
                            PPgr.loc[
                                np.isclose(PPgr[EvRHE], -0.6, atol=0.015)
                            ].drop_duplicates(subset=ExportECfromHER.HER_cols)
                            E_06_lst.append(
                                [
                                    sID,
                                    PPn,
                                    PPgr.loc[
                                        np.isclose(PPgr[EvRHE], -0.6, atol=0.015)
                                    ].drop_duplicates(subset=ExportECfromHER.HER_cols),
                                ]
                            )
                            DTgr.iloc[5:-20].plot(
                                x=EvRHE,
                                y="j_lower",
                                c=RGB,
                                label=label_set,
                                ax=ax,
                                ls=ls_set,
                            )
                            xls_label = "{0}_{1}_{2}_{3}.xlsx".format(
                                DTgr.SampleLabel.unique()[0],
                                sID[-1],
                                pAST_set,
                                datetime.strftime(DT, format="%Y-%M-%d"),
                            )
                            DTgr[ExportECfromHER.HER_cols].to_excel(
                                PPDHER.joinpath(xls_label)
                            )
        ax.set_xlim((-0.8, 0))
        ax.set_ylabel(r"$\mathrm{j\/ / \/ mAcm^{-2}}$")
        ax.set_xlabel(r"$\mathrm{E\/ / \/ V_{RHE}}$")
        ax.set_title(", ".join(["{0} : {1}".format(i[0], i[1]) for i in EC_exp_uniq]))
        ax.legend(ncol=1, loc="lower right", fontsize=11)
        plt.savefig(
            PPDHER.joinpath("HER_CVs_{0}.png".format(HER_E_slice.PAR_date.nunique()))
        )
        #        bbox_to_anchor=(0.5,0)
        HER_06 = pd.concat([i[2] for i in E_06_lst]).set_index(
            ["SampleID", "SampleLabel"] + SampleSelection.EC_exp_cols
        )

        HER_06.plot().barh(x="SampleID", y="j_lower")
        #        for sID, gr HER_06.groupby('SampleID').groups
        corr_Cols = (
            SampleSelection.InterestingCols
            + SampleSelection.EC_EIS_par_cols
            + SampleSelection.RAMAN_cols_corr
            + SampleSelection.Inter_MS_cols
        )
        drop_corr_cols = [
            "Colorcode",
            "RedChisqr2",
            "RedChisqr1",
            "Area_in_cell",
            "BET_Area_RPT",
        ]
        corr_Cols_filtered = [i for i in corr_Cols if i not in drop_corr_cols]
        #        corr_method,corr_cutoff = 'spearman',0.1 # or default is pearson spearman
        #        rcorr = EIS_O2_065_acid_no[corr_Cols].corr(method=corr_method)
        out_topcorrs_lst = []
        for Gas_set in Gases:
            for pH_set in pHses:
                EIS_O2_no_query = EIS_CB_paper.query(
                    '(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
                ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
                #                ORREIS_O2_no_query = EIS_CB_paper.query('(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))').drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
                target_dir = PPDEIS.joinpath(
                    "EIS_corr_{0}_pH{1}".format(Gas_set, pH_set)
                )
                target_dir.mkdir(parents=True, exist_ok=True)
                for EvRHE, grE in EIS_O2_no_query.groupby("E_RHE"):
                    if len(grE) > 3:
                        EvRHE, grE
                        rcorr = grE[corr_Cols_filtered].corr(method=corr_method)
                        prom_corr = ExportECfromEIS.plot_triangle_hm(
                            grE,
                            rcorr,
                            np.round(EvRHE, 2),
                            target_dir,
                            corr_method_set=corr_method,
                            plot_option=False,
                        )
                        out_topcorrs_lst.append(
                            [Gas_set, pH_set, EvRHE, len(grE), prom_corr]
                        )

        #        pd.DataFrame(out_topcorrs_lst[0][-1],columns=[corr_method]).assign(**{'Gas' : out_topcorrs_lst[0][0], 'pH' : out_topcorrs_lst[0][1], 'E_RHE' : out_topcorrs_lst[0][2]})
        topcorrs = pd.concat(
            [
                pd.DataFrame(i[-1], columns=[corr_method]).assign(
                    **{"Gas": i[0], "pH": i[1], "E_RHE": i[2], "lenGr": i[3]}
                )
                for i in out_topcorrs_lst
            ]
        )
        topcorrs["score"] = np.abs(topcorrs[corr_method]) * topcorrs.lenGr
        topcorrs = topcorrs.sort_values("score", ascending=0)
        topcorr_best = topcorrs[
            (np.abs(topcorrs[corr_method]) > 0.7)
            & (topcorrs.lenGr >= 0.9 * topcorrs.lenGr.describe()["75%"])
        ].sort_values(by=corr_method)
        tpcb_top = pd.concat(
            [topcorr_best.head(50), topcorr_best.tail(50)], sort=corr_method
        )
