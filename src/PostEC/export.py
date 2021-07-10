#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 14:44:15 2020

@author: zmg
"""


class MakingCorrs:

    #    ECnormal_query =

    def __init(self, Pars):
        self.Pars = Pars

    @staticmethod
    def CheckCols(coll_lst, df):
        return [i for i in coll_lst if i in df.columns]

    @staticmethod
    def make_correlations():
        """Making plots of correlations of EIS Parameters versus E v RHE"""
        Gases, pHses = ["N2", "O2"], [0.3, 1, 13]
        #        EIS_O2_no = EIS_CB_paper.query('(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))').drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        #                pH03, pH1 = EIS_O2_acid_no.query('(pH == 0.3)'), EIS_O2_acid_no.query('(pH == 1)')
        corr_Cols = (
            SampleSelection.InterestingCols
            + SampleSelection.EC_EIS_par_cols
            + SampleSelection.RAMAN_cols_corr
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

        PPDEISbest = PPDEIS.joinpath("top_corrs_{0}".format(corr_method))
        PPDEISbest.mkdir(parents=True, exist_ok=True)
        for tr in tpcb_top.iterrows():
            Gas_set, pH_set, Erhe_set, corr_val = (
                tr[1].Gas,
                tr[1].pH,
                tr[1].E_RHE,
                tr[1][corr_method],
            )
            grE = EIS_O2_no_normalload.query(
                "(Gas == @Gas_set) & (pH == @pH_set) & (E_RHE == @Erhe_set)"
            )
            xc, yc = tr[0][0], tr[0][1]
            plot_pd_SampleIDs(grE, xc, yc, corr_val, PPDEISbest)
            #            spectras = ExportECfromCV.make_uniform_EvRHE(pd.concat([pd.read_excel(i) for i in grE.SpectraFile.values],sort=False)).query('(E_RHE == @Erhe_set)')
            spectra_lst = []
            for n, spfrow in grE.iterrows():
                spf = spfrow.SpectraFile
                spdf = ExportECfromCV.make_uniform_EvRHE(pd.read_excel(spf))
                spdf_Ev = spdf.query("(E_RHE == @Erhe_set)")
                spdf_char = [
                    spdf_Ev.assign(
                        **{i: [spfrow[i]] * len(spdf_Ev.index) for i in spfrow.index}
                    )
                ][0]
                spectra_lst.append(spdf_char)
            spectras = pd.concat([i for i in spectra_lst], sort=False)
            plot_spectra_SampleIDs(xc, yc, spectras, PPDEISbest)

    def CharacterizationsOnly(SampleCodes):
        SeriesIDs = [
            SampleSelection.Series_CB_paper,
            SampleSelection.Series_Porhp_SiO2,
            {"name": "all_EIS"},
        ]
        SeriesID_set = SeriesIDs[0]

        if "all" in SeriesID_set["name"]:
            SampleCodesSerie = SampleCodes
        else:
            SampleCodesSerie = SampleCodes.loc[
                SampleCodes.SeriesID.isin(SeriesID_set["slice"])
            ]

        ddir = (
            FileHelper.FindExpFolder("VERSASTAT")
            .PostDir.joinpath("StructuralCorrs")
            .joinpath("{0}".format(SeriesID_set["name"]))
        )
        ddir.mkdir(parents=True, exist_ok=True)
        corr_Cols_unfiltered = (
            SampleSelection.InterestingCols
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
        corr_Cols = [
            i
            for i in corr_Cols_unfiltered
            if i in SampleCodes.columns and i not in drop_corr_cols
        ]
        corr_method, corr_cutoff = "spearman", 0.1  # or default is pearson spearman

        corr_stack = SampleCodesSerie[corr_Cols].corr(method=corr_method).stack()
        rcorr = corr_stack.sort_values()  # [(corr_stack < 1) & (corr_stack > -1)]
        bet1, bet2 = [
            1
            if i[0] in SampleSelection.Inter_BET_cols + SampleSelection.Sample_Lst_BET
            else 0
            for i in rcorr.index
        ], [
            1
            if i[1] in SampleSelection.Inter_BET_cols + SampleSelection.Sample_Lst_BET
            else 0
            for i in rcorr.index
        ]
        bets = [i[0] + i[1] for i in zip(bet1, bet2)]
        ea1, ea2 = [
            1 if i[0] in SampleSelection.Inter_EA_cols else 0 for i in rcorr.index
        ], [1 if i[1] in SampleSelection.Inter_EA_cols else 0 for i in rcorr.index]
        eas = [i[0] + i[1] for i in zip(ea1, ea2)]
        raman1, raman2 = [
            1
            if i[0]
            in SampleSelection.Inter_RAMAN_cols + SampleSelection.RAMAN_cols_corr
            else 0
            for i in rcorr.index
        ], [
            1
            if i[1]
            in SampleSelection.Inter_RAMAN_cols + SampleSelection.RAMAN_cols_corr
            else 0
            for i in rcorr.index
        ]
        ramans = [i[0] + i[1] for i in zip(raman1, raman2)]

        rcorr_fltr = pd.DataFrame(rcorr, columns=[corr_method]).assign(
            **{"BETs": bets, "EAs": eas, "RAMANs": ramans}
        )
        rcorr_fltr = rcorr_fltr.query(
            "BETs < 2 & EAs < 2 & RAMANs < 2 & {0} < 1 & {0} > -1".format(corr_method)
        ).drop_duplicates(keep="first")
        #        rcorr_fltr['sums'] = rcorr_fltr[['BETs','EAs','RAMANs']].sum(axis=1)
        rcorr_fltr.to_excel(
            ddir.joinpath("{0}_corr_fltrd.xlsx".format(SeriesID_set["name"]))
        )
        matcorrs_best = pd.concat([rcorr_fltr.head(20), rcorr_fltr.tail(20)])
        for tr in matcorrs_best.iterrows():
            x_col, y_col, corr_val = tr[0][0], tr[0][1], tr[1][corr_method]
            plot_pd_SampleIDs(SampleCodesSerie, x_col, y_col, corr_val, ddir)


class DataPerSample:
    def __init__(self):
        pass

    def PARS_data_per_Sample():
        #%% ===== Collect all Data and Parameters per Electrolyte for each Sample ======
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        postOVVout = LoadPostOVV()
        #    pd.read_excel(PostDestDir.joinpath('postEC_Organized.xlsx'),index_col=[0])
        #    postEIScom = postOVVout.loc[postOVVout['Type_Exp'] == 'EIS_Combined'].drop_duplicates()
        PostDestDirEIScom = PostDestDir.joinpath("EIS_PARS_ORR_N2")

        for Elec, ElecGr in postOVVout.groupby(by="Electrolyte"):
            PDDirEIScom = PostDestDirEIScom.joinpath(Elec)
            PDDirEIScom.mkdir(parents=True, exist_ok=True)
            for Sample, sGr in ElecGr.groupby(by="SampleID"):
                sPARs, sLst = pd.DataFrame(), []
                for sType, sTgr in sGr.groupby(by="Type_Exp"):
                    for sf, fGr in sTgr.groupby(by="SourceFilename"):
                        print("_".join([Elec, Sample, sType, sf]))
                        AllData_E_file = FileHelper.FileOperations.PDreadXLorCSV(sf)
                        AllData_E_file = AllData_E_file.assign(
                            **{
                                "Electrolyte": Elec,
                                " SampleID": Sample,
                                "Type_Exp": sType,
                                "SourceFilename": sf,
                                "Gas": fGr.Gas.values[0],
                                "Status": fGr.Status.values[0],
                                "EXP_date": fGr.EXP_date.values[0],
                            }
                        )
                        AllData_E_file.set_index(
                            [
                                "Electrolyte",
                                " SampleID",
                                "Type_Exp",
                                "Gas",
                                "Status",
                                "EXP_date",
                                "SourceFilename",
                            ]
                        )
                        AllData_E_file
                        sLst.append([AllData_E_file])
                #                    EISovv = fGr
                DestFile = PDDirEIScom.joinpath(
                    "_".join([Sample, Elec, str(len(sLst))])
                ).with_suffix(".xlsx")
                pd.concat([i[0] for i in sLst], sort=False).drop_duplicates().to_excel(
                    DestFile
                )

    #                PlotCombinedEIS(AllData_E_file,EISovv,DestFile)
    #                SwpFltr['Jcorr']= Jcorr_fltr
    #                SwpFltr['J_ring']= Jring_fltr
    #                SwpFltr['Frac_H2O2']= FracH2O2_fltr
    #                SwpFltr['J_N2_scan']= N2_fltr
    #                SwpFltr['n_ORR']= nORR_fltr
    #    def SelectSamples_DFG():
    #        DFG = PrepOVV.loc[PrepOVV.SeriesID.str.contains('MG') == True]

    def FilternaPlot(df, xl, yl, plotk, DestDir):
        fltrd = df.dropna(subset=[yl])
        fltrd = fltrd.sort_values(by="SeriesID", ascending=False)
        fig, ax = plt.subplots()
        fltrd.plot(x=xl, y=yl, kind=plotk, ax=ax, legend=False)
        #        plt.legend(False)
        fig.suptitle("%s" % yl)
        ax.set_xlabel("wt%")
        plt.savefig(
            DestDir.joinpath("%s" % yl.replace("/", "_")).with_suffix(".png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.show()
        return fltrd

    def Retrieve_All_ExpFiles(SerieOVV, DD_Series):
        Exp_data_dirs3 = [
            i
            for i in FileHelper.FindExpFolder("EA").DataDir.parent.glob("*")
            if i.is_dir()
        ]
        Exp_data_dirs = [
            i
            for i in Exp_data_dirs3
            if not i.name in ["VERSASTAT", "EA", "BET", "BAK", "BET_Plots", "PARSTAT"]
        ]
        extra_chars = []
        for sID in SerieOVV.SampleID.unique():
            for expd in Exp_data_dirs3:
                EA_Samples = []
                if expd.name == "EA":
                    EA_results = pd.read_excel(
                        list(FileHelper.FindExpFolder("EA").DestDir.rglob("*xlsx"))[0]
                    )
                    EA_Samples = EA_results.SampleID.unique()
                if list(expd.rglob("*%s*" % sID)) or sID in EA_Samples:
                    file_list_sID = list(expd.rglob("*%s*" % sID))
                    extra_chars.append(
                        [sID, expd.name, expd, file_list_sID, len(file_list_sID)]
                    )
        All_exp_files = pd.DataFrame(
            extra_chars,
            columns=["SampleID", "Exp_Type", "Exp_Dir", "Exp_files", "Exp_files_count"],
        )
        SerieOVV_exp_data = pd.merge(
            SerieOVV, All_exp_files.set_index("SampleID"), on="SampleID"
        )[
            [
                "SampleID",
                "SampleLabel",
                "SeriesID",
                "Exp_Type",
                "Exp_Dir",
                "Exp_files",
                "Exp_files_count",
            ]
        ]
        SerieOVV_exp_data.to_excel(
            DD_Series.parent.joinpath("{0}_ExpFiles.xlsx".format(DD_Series.name)),
            index=False,
        )
        #        SerieOVV_exp_data.to_excel(DD_Series.parent.joinpath('{0}_ExpFiles%'.format(DD_Series.name)),drop_index=True)
        return SerieOVV_exp_data
