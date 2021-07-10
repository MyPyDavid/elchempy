#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 14:43:02 2020

@author: zmg
"""
#
# import FileHelper
# from FileHelper.PostChar import SampleSelection, Characterization_TypeSetting
# import pandas as pd
# import export


class MergeEISandCdl:
    """Made correlations by guessing for best,
    but should probably choose several E_RHE values eg. [0.2,0.4,0.6,0.8,1] to narrow it down."""

    EvRHE = "E_AppV_RHE"
    #    OriginColor = FileHelper.FindExpFolder().LoadOriginColor()
    """ Refit: pH1 : DW16,DW19,'DW17','DW28'
    """

    def __init__(self):
        pass

    @staticmethod
    def splitcol_Sweep_Cdl(Cdl_pars):
        d1 = Cdl_pars.iloc[0].to_dict()
        d1_set_cols = [
            key
            for key in d1
            if "set" in str(type(d1[key])) or "list" in str(type(d1[key]))
        ]
        for d1c in d1_set_cols:
            Cdl_pars = Cdl_pars.assign(
                **{d1c: [frozenset(i) for i in Cdl_pars[d1c].values]}
            )

        Cdl_cath, Cdl_an = Cdl_pars.query(
            'Sweep_Type_N2 == "cathodic"'
        ), Cdl_pars.query('Sweep_Type_N2 == "anodic"')
        merge_cols_catan = [i for i in Cdl_cath.columns if i in Cdl_an.columns]
        Cdl_catan = pd.merge(
            Cdl_cath,
            Cdl_an,
            on=[i for i in merge_cols_catan if i not in SampleSelection.EC_N2Cdl_cols],
            how="left",
            suffixes=["_cat", "_an"],
        )
        Cdl_catan["Cdl_sum"] = Cdl_catan["Cdl_an"] + Cdl_catan["Cdl_cat"]
        return Cdl_catan

    def merge_EIS_Cdl_pars(EIS_pars, Cdl_pars_catan):
        """Merge on EIS_pars with SampleCodes columns , basically adding Cdl_pars to it merging on ECexp and E_RHE"""
        EIScols = MakingCorrs.CheckCols(
            SampleSelection.EC_EIS_par_cols + ["ECexp", "E_RHE", "SampleID"], EIS_pars
        ) + MakingCorrs.CheckCols(SampleSelection.EC_exp_cols, Cdl_pars_catan)
        catancols = [
            i + a
            for a in ["_an", "_cat"]
            for i in SampleSelection.EC_N2Cdl_cols
            if i + a in Cdl_pars_catan.columns
        ] + ["Cdl_sum"]
        Cdlcatancols = Cdl_pars_catan[
            catancols
            + ["ECexp", "E_RHE"]
            + MakingCorrs.CheckCols(SampleSelection.EC_exp_cols, Cdl_pars_catan)
        ]

        EIS_Cdl_pars = pd.merge(
            EIS_pars,
            Cdlcatancols,
            on=["ECexp", "E_RHE"],
            how="left",
            suffixes=["_eis", "_cdl"],
        )
        eec = SampleSelection.EC_EIS_models[0]
        EC_uniq_data_set = EIS_Cdl_pars.loc[
            EIS_Cdl_pars.SampleID.isin(SampleSelection.Series_CB_paper["sIDs"])
        ].query(
            "E_RHE == 0.5 & pH_eis <= 2 & Model_EEC == @eec & ((Loading_cm2_eis < 0.6) & (Loading_cm2_eis > 0.3))"
        )
        EC_uniq_data_set_grpE = EC_uniq_data_set.groupby("E_RHE")
        EC_uniq_data_set.plot(
            x="Qad+Cdlp",
            y="Cdl_an",
            kind="scatter",
            c="BET_cat_agg",
            colormap="rainbow",
        )

        def correlate_Cdl_with_its_cath_an():
            EIS_Cdl_pars["Cdl_ratio_catan"] = (
                EIS_Cdl_pars["Cdl_cat"] / EIS_Cdl_pars["Cdl_an"]
            )
            EIS_Cdl_pars["Cdl_diff_catan"] = (
                EIS_Cdl_pars["Cdl_cat"] - EIS_Cdl_pars["Cdl_an"]
            )
            EC_uniq_data_set = EIS_Cdl_pars.loc[
                EIS_Cdl_pars.SampleID.isin(SampleSelection.Series_CB_paper["sIDs"])
            ].query(
                "pH_eis <= 2 & Model_EEC == @eec & ((Loading_cm2_eis < 0.6) & (Loading_cm2_eis > 0.3))"
            )
            EC_uniq_data_set.plot(
                x="E_RHE",
                y="Cdl_diff_catan",
                kind="scatter",
                c="BET_cat_agg",
                colormap="rainbow",
            )
            EC_uniq_data_set_grpE = EC_uniq_data_set.groupby("E_RHE")

            Cdl_itself_Ecorr_lst = []
            for E, grE_r in EC_uniq_data_set_grpE:
                grE = grE_r.dropna(subset=["Cdl_cat", "Cdl_an"])
                if not grE.empty:
                    corrE = grE[["Cdl_cat", "Cdl_an"]].corr().loc["Cdl_cat", "Cdl_an"]
                    linregress = scipy.stats.linregress(
                        grE["Cdl_cat"].values, grE["Cdl_an"].values
                    )
                    #                slope, intercept, r_value, p_value, std_err
                    Cdl_itself_Ecorr_lst.append(
                        {
                            "E_RHE": E,
                            "Cdl_catan_corr": corrE,
                            "slope": linregress.slope,
                            "intercept": linregress.intercept,
                            "r_val": linregress.rvalue,
                        }
                    )
                #                    grE.plot(x='Cdl_cat',y='Cdl_an',kind='scatter',c='BET_cat_agg',colormap='rainbow',title='E = {0} V'.format(E))
                else:
                    pass
            Cdl_itself_Ecorr = pd.DataFrame(Cdl_itself_Ecorr_lst).query("r_val > 0.5")
            Cdl_itself_Ecorr.plot(x="E_RHE", y="slope")
            Cdl_itself_Ecorr.plot(x="E_RHE", y="intercept")

    @staticmethod
    def EIS_Cdl(EIS_pars, Cdl_pars):
        EIS_pars
        Cdl_pars
        PDDmergeCdl = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
            "PostEC/MergeEIS_Cdl"
        )
        PDDmergeCdl.mkdir(parents=True, exist_ok=True)

        Cdl_Porph_SiO2, EIS_Porph_SiO2
        Cdl_Porph_SiO2_cath, Cdl_Porph_SiO2_an = Cdl_Porph_SiO2.query(
            'Sweep_Type_N2 == "cathodic"'
        ), Cdl_Porph_SiO2.query('Sweep_Type_N2 == "anodic"')
        merge_cols_catan = [
            i for i in Cdl_Porph_SiO2_cath.columns if i in Cdl_Porph_SiO2_an.columns
        ]
        Cdl_catan = pd.merge(
            Cdl_Porph_SiO2_cath,
            Cdl_Porph_SiO2_an,
            on=[i for i in merge_cols_catan if i not in SampleSelection.EC_N2Cdl_cols],
            how="left",
            suffixes=["_cat", "_an"],
        )
        Cdl_catan["Cdl_sum"] = Cdl_catan["Cdl_an"] + Cdl_catan["Cdl_cat"]

        EIS_pars["PAR_date_day"] = [
            pd.datetime.strftime(pd.to_datetime(i), format="%Y-%m-%d")
            for i in EIS_pars.PAR_date.fillna(0).values
        ]
        Cdl_pars["PAR_date_day"] = [
            pd.datetime.strftime(pd.to_datetime(i), format="%Y-%m-%d")
            for i in Cdl_pars.PAR_date.fillna(0).values
        ]
        merge_cols = [i for i in EIS_Porph_SiO2.columns if i in Cdl_Porph_SiO2.columns]

        ep2 = EIS_Porph_SiO2[
            SampleSelection.EC_EIS_par_cols
            + ["ECexp", "E_RHE"]
            + MakingCorrs.CheckCols(SampleSelection.EC_exp_cols, Cdl_pars)
        ]
        cp2cath = Cdl_Porph_SiO2_cath[
            SampleSelection.EC_N2Cdl_cols
            + ["ECexp", "E_RHE"]
            + MakingCorrs.CheckCols(SampleSelection.EC_exp_cols, Cdl_pars)
        ]
        cp2an = Cdl_Porph_SiO2_an[
            SampleSelection.EC_N2Cdl_cols
            + ["ECexp", "E_RHE"]
            + MakingCorrs.CheckCols(SampleSelection.EC_exp_cols, Cdl_pars)
        ]

        ec3 = pd.merge(
            ep2, cp2an, on=["ECexp", "E_RHE"], how="left", suffixes=["_eis", "_cdl"]
        )

        ec3.query("(pH_eis > 1)").plot(x="Qad+Cdlp", y="Cdl", kind="scatter")
        EIS_Cdl_pars_normal = EIS_Cdl_pars.query(
            '(pH_eis == 1) & (postAST_eis == "no") & ((Loading_cm2_eis < 0.5) & (Loading_cm2_eis > 0.3))'
        )
        EIS_Cdl_pars_normal.plot(x="Qad+Cdlp", y="Cdl_sum", kind="scatter")

        for yc in ["Cdlp", "Qad", "Qad+Cdlp", "nAd", "nDL"]:
            for Ev in [0, 0.2, 0.5, 0.7, 0.9]:
                fig, ax = plt.subplots()
                EIS_Cdl_pars_normal.loc[
                    EIS_Cdl_pars_normal.SampleID.isin(
                        SampleSelection.Series_Porhp_SiO2["sIDslice"]
                    )
                ].query("E_RHE == @Ev").plot(
                    x="Cdl_sum",
                    y=yc,
                    kind="scatter",
                    c="BET_cat_agg",
                    colormap="viridis",
                    ax=ax,
                )
                ax.set_xlim([0, 0.06])
                ax.set_ylim([0, 0.06])
                if yc in ["nAd", "nDL"]:
                    ax.set_ylim([0, 1])
                else:
                    ax.set_ylim([0, 0.03])
                #            ax.autoscale(True)
                plt.show()
                plt.close()  #

        for ph in [0.3, 1, 13]:
            fig, ax = plt.subplots(figsize=(14, 12))
            plt.suptitle("pH = {}".format(ph))
            Cdl_catan.query("(pH == @ph)").plot(
                x="Cdl_cat",
                y="Cdl_sum",
                kind="scatter",
                c="E_RHE",
                colormap="viridis",
                ax=ax,
                marker="^",
                label="PorphSiO2",
                s=60,
            )
            Cdl_CB_catan.query("(pH == @ph)").plot(
                x="Cdl_cat",
                y="Cdl_sum",
                kind="scatter",
                c="E_RHE",
                colormap="jet",
                ax=ax,
                marker="*",
                label="PorphCB",
                s=60,
            )
            #        Cdl_catan.query('(pH == 1)').plot(x='Cdl_cat',y='Cdl_sum',kind='scatter',c='E_RHE',colormap='viridis',ax=ax)
            ax.set_xlim([0, 0.06])
            ax.set_ylim([0, 0.06])
            ax.grid(True)
            plt.savefig(PDDmergeCdl.joinpath("Cdl_sum_an_Series_{0}.png".format(ph)))

        MergeCols = ["ECexp", "E_RHE"] + Load_from_Indexes.EC_label_cols[0:-1]
        EIS_merge_cols, Cdl_merge_cols = [
            i for i in EIS_pars.columns if i not in SampleCodes and i not in MergeCols
        ] + MergeCols, [
            i for i in Cdl_pars.columns if i not in SampleCodes and i not in MergeCols
        ] + MergeCols
        #          [i for i in EIS_merge_cols if i in Cdl_merge_cols]
        eis_cdl = pd.merge(
            EIS_Porph_SiO2[EIS_merge_cols],
            Cdl_Porph_SiO2[Cdl_merge_cols],
            on=MergeCols[1::],
            how="left",
        )

        eis_cdl.query("(pH == 1) ").plot(x="Qad+Cdlp", y="Cdl", kind="scatter")
        ecslice = eis_cdl.query("(pH == 1) & (E_RHE == 0.5)")
