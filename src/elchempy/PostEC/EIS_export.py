#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:03:37 2020

@author: zmg
"""
from pathlib import Path
import sys
from scipy.stats import linregress
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import shutil
import json

import matplotlib as mpl

mpl.style.use("seaborn")
mpl.rcParams["figure.dpi"] = 100

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper.PostChar import (
    SampleSelection,
    Characterization_TypeSetting,
    SampleCodesChar,
)


if __name__ == "__main__":
    print(f"Package: {__package__}, File: {__file__}")

    from ECpy.experiments.ORR import ORR_analyze_scans
    from ECpy.PostEC.plotting import eisplot
    from ECpy.PostEC.collect_load import Load_from_Indexes
    from ECpy.experiments.EIS.models import Model_Collection  # EEC_models_index

    #    logger = start_logger(Path.cwd())
    #    import plotting
    import post_helper

else:
    print("\n\n***** run_PAR_DW *****")
    print(__file__)

    from ECpy.PostEC.plotting import eisplot

    #    from ECpy.experiments.EIS import plotting
    from ECpy.PostEC import post_helper
    from ECpy.PostEC.collect_load import Load_from_Indexes
    from ECpy.experiments.ORR import ORR_analyze_scans
    from ECpy.experiments.EIS.models import EEC_models_index


OriginColor = FindExpFolder().LoadOriginColor()
EvRHE = "E_AppV_RHE"
SampleCodes = SampleCodesChar.load()


def check_file_and_backup(all_fits_df, filepath, destdir, _bakd={}, suffix=".pkl"):
    #    destdir=PDD_eischeck
    filepath = Path(destdir.joinpath(filepath).with_suffix(suffix))
    if filepath.is_file():
        date_prefix = datetime.fromtimestamp(filepath.stat().st_ctime).strftime(
            "%Y-%m-%d"
        )
        bak_stem = filepath.parent.joinpath(
            f"{date_prefix}_{filepath.stem}"
        ).with_suffix(suffix)
        _old_df = pd.read_pickle(filepath)
        _bakd.update({bak_stem: len(_old_df)})
        filepath.replace(bak_stem)
    all_fits_df.to_pickle(filepath)
    return _bakd


def unique_destdir(uniq_str):
    if uniq_str:
        PostDestDir = FindExpFolder("VERSASTAT").PostDir.joinpath(f"EIS/{uniq_str}")
        PostDestDir.mkdir(parents=True, exist_ok=True)
        print(f"Dest dir: {PostDestDir}")
        return PostDestDir


class EIS_selection:
    filter = "(lmfit_MSE < 65E4) & (Rct < 1E3) & (Rct > 2E-2) \
                                       & (Rs > 0.01)  & (Rs < 200) & (Cdlp < 0.075)\
                                       & (lmfit_redchi < 1E3)  & (Aw < 10E3)\
                                       & (Aw > 10E-2) & (Qad < 1) & (tau < 10)"

    _old_fast_checking_EEC_models = [
        "Model(R0-L0-p(R1,CPE1)-p(R2-W2,CPE2))",
        "Model(R0-L0-p(R1-Wo1,CPE1)-C2)",
        "Model(R0-L0-p(R1-Ws1,CPE1)-C2)",
        "Model(R0-L0-W0-p(R1,CPE1)-p(R2,CPE2))",
        "Model(R0-L0-p(R1-W1,CPE1)-C2)",
        "Model(R0-L0-W0-p(R1,CPE1)-p(R2,C2))",
        "Model(R0-L0-p(R1-W1,CPE1)-CPE2)",
        "Model(R0-L0-p(R1-Wo0,C1)-W0)",
        "Model(R0-L0-p(R1-W1,C1)-C2)",
        "Model(R0-L0-p(R1-Wo1,CPE1))",
        "Model(R0-L0-p(R1-Wo0,CPE1)-W0)",
        "Model(R0-L0-p(R1-W1,CPE1))",
        "Model(R0-L0-p(R1-Ws1,CPE1))",
    ]
    mod_select = {
        "N2": "Model(RL-TLM(Rct-Qad-W))",
        "O2": "Model(RL-TLM(Rct-p(Qad-W,Rorr)))",
    }

    loadgrp_cols = ["SampleID", "Electrolyte", "E_RHE", "Gas", "Model_EEC"]
    ECuniq_cols = ["SampleID", "Loading_cm2", "Electrolyte", "E_RHE", "Gas", "RPM_DAC"]

    def EIS_filter():
        _filter = "(EIS_pars.lmfit_MSE < 65E4) & (EIS_pars.Rct < 2E3) & (EIS_pars.Rct > 2E-4) \
                                     & (EIS_pars.Rs > 0.01)  & (EIS_pars.Rs < 200) & (EIS_pars.Cdlp < 0.075)\
                                     & (EIS_pars.lmfit_redchi < 1E3)  & (EIS_pars.Aw < 10E3) & (EIS_pars.Aw > 10E-2)\
                                     & (EIS_pars.Qad < 0.5) & (EIS_pars.tau < 3E2)"
        return _filter

    # def add_best_model_per_spectrum():
    #     _grpkeys = ['PAR_file', 'Segment #']

    #     _tt = EIS_pars.query('SampleID == "JOS4" & Gas == "O2"')
    #     for (pf,seg),PF_pars in _tt.groupby(_grpkeys):
    #         (pf,seg),PF_pars

    # def find_best_model_per_spectrum(PF_pars, var_lim = 1E5, var_err_lim = 3E3):

    #     PF_pars = PF_pars.loc[PF_pars.lmfit_message.str.contains('satisfied') == True]

    #     _lmfitcols = [i for i in PF_pars.columns if i.startswith('lmfit')]

    #     if PF_pars.Model_EEC.nunique() > 2 and PF_pars.Model_EEC.nunique() == len(PF_pars):
    #         _res = []
    #         for n,r in PF_pars.iterrows():
    #             _vars = r.lmfit_var_names.split(', ')
    #             _varserr = [i for i in [i+'_stderr' for i in _vars] if i in r.index]
    #             _vsum,_vmax = r[_vars].sum(),r[_vars].max()
    #             _bad_vars = set([i for i in _vars if r[i] > var_lim])
    #             _verrsum,_verrmax = 0, 0
    #             if _varserr:
    #                  _verrsum,_verrmax = r[_varserr].sum(),r[_varserr].max()
    #                  _bad_varserr = {i:{'val' : r[i],'rel_val' : r[i]/r[(i.split('_stderr')[0])], 'name' : (i.split('_stderr')[0])}  for i in _varserr}
    #                  _bad_varserr_val = [val['name'] for k,val in _bad_varserr.items() if val['val'] > var_err_lim]
    #                  _bad_varserr_perc = [val['name'] for k,val in _bad_varserr.items() if val['rel_val'] > var_err_lim]
    #                  _bad_vars_err = set(_bad_varserr_val + _bad_varserr_perc)
    #                  _bad_vars= _bad_vars.union(_bad_vars_err)
    #             _res.append([n,r.Model_EEC_name, len(_vars), _vsum, _vmax, ', '.join(_bad_vars), len(_bad_vars),_verrsum,_verrmax, r.lmfit_aic,r.lmfit_redchi, r.lmfit_chiqsr])
    #         var_res = pd.DataFrame(_res,columns=['pf_index','Model_EEC_name','len_vars','varsum','varmax','bad_vars','bad_vars_len','err_varsum','err_varmax',
    #                                    'lmfit_aic','lmfit_redchi','lmfit_chiqsr'])
    #         var_res = var_res.sort_values(by=['bad_vars_len','lmfit_aic','len_vars'])
    #         best_mod_row = var_res.head(1)
    #         _sorted_rank = ', '.join([str(i) for i in var_res.pf_index.values])
    #         _best_result = [best_mod_row.pf_index.iloc[0],best_mod_row.Model_EEC_name.iloc[0],_sorted_rank]
    #         var_res.bad_vars_len.unique()


def select_models_EISpars(EIS_pars):
    N2_mod = "Model(RL-TLM(Rct-Qad-W))"
    O2_mod = "Model(RL-TLM-Rct-W-p(Rorr,Qad))"
    N2_O2_models = EIS_pars.loc[
        ((EIS_pars.Model_EEC == N2_mod) & (EIS_pars.Gas == "N2"))
        | ((EIS_pars.Model_EEC == O2_mod) & (EIS_pars.Gas == "O2"))
    ]
    return N2_O2_models


def plot_pars_E(Ekin, E_lower=0.59, E_upper=0.76, Xcol="BET_cat_agg"):
    if Xcol in Ekin.columns:
        parcols = [
            i
            for i in Ekin.columns
            if i
            in [i for i in SampleSelection.EC_EIS_par_cols if not "lmfit_redchi" in i]
            + SampleSelection.EC_ORR_kin_par_cols
        ]
        for i in parcols:
            EKinsl = Ekin.query(
                f"E_RHE < {E_upper} & E_RHE > {E_lower} & Rct < 9E05 & Rorr < 1E09 & Qad < 35E-3 & Cdlp < 0.070"
            )
            for Elec, Elgr in EKinsl.groupby("Electrolyte"):
                fig, ax = plt.subplots()
                Elgr.plot(
                    x=Xcol, y=i, kind="scatter", c="E_RHE", colormap="rainbow_r", ax=ax
                )
                ax.set_xlim = (0.6, 0.9)
                ax.set_title(Elec)
                ps = plotting.eisplot(i)
                ax.set_ylim(ps.ylim)
                ax.set_yscale(ps.logyscale)
                plt.show()
                plt.close()


def EIS_all_check_redchi(EIS_pars, eischeck_plot=False):
    PostDestDir = FindExpFolder("VERSASTAT").PostDir
    PDD_eischeck = PostDestDir.joinpath("EIS/redchi_check")
    PDD_eischeck.mkdir(parents=True, exist_ok=True)

    SeriesIDs = [
        SampleSelection.Series_CB_paper,
        SampleSelection.Series_Porhp_SiO2,
        SampleSelection.Series_Co_PANI,
        SampleSelection.Series_ML_SiO2,
    ]
    all_sIDs = [a for i in [i.get("sIDs", []) for i in SeriesIDs] for a in i]
    EvRHE = "E_AppV_RHE"
    EIS_pars_all = EIS_pars.loc[
        EIS_pars.SampleID.isin(all_sIDs)
        & (EIS_pars.Rs > 1)
        & (EIS_pars.lmfit_message.str.contains("Fit succeeded"))
    ].drop_duplicates(subset=["PAR_file", EvRHE, "Model_EEC"])
    EIS_pars_all[EvRHE] = EIS_pars_all[EvRHE].round(3)
    all_good_lst, all_bad_lst, _bad_fit_suggestions = [], [], []
    mod_index = EEC_models_index()
    EIS_pars_all_grp_PF_mod = EIS_pars_all.groupby(["PAR_file", "Model_EEC"])
    for (pf, mod), Mgrp in EIS_pars_all_grp_PF_mod:
        #        pf,mod, Mgrp
        PDD_eischeck_figs = PDD_eischeck.joinpath(Path(pf).parent.name)
        PDD_eischeck_figs.mkdir(parents=True, exist_ok=True)
        #        egrp = Mgrp.groupby(EvRHE)
        if len(Mgrp) > 2:
            savestem = f'{Path(pf).stem}_{mod.split("Model")[-1]}'
            Mgrp2 = Mgrp
            var_names_uniq = Mgrp2.lmfit_var_names.unique()[0]
            if pd.isna(var_names_uniq):
                var_names = [i for i in mod_index if mod in i[1].name][0][1].param_names
            else:
                var_names = list(
                    map(
                        str.strip,
                        var_names_uniq.strip(")(").replace("'", "").split(","),
                    )
                )
            #            var_names = Mgrp2.lmfit_var_names.unique()[0]
            var_zscores = {
                var + "_zscore": (np.abs(stats.zscore(Mgrp[var])), Mgrp[var].index)
                for var in var_names
            }

            z_indexes = [
                i
                for key, val in var_zscores.items()
                for z, i in zip(val[0], val[1])
                if z < 1.45
            ]

            #            Mgrp[Mgrp[[var+'_zscore' for var in var_names]] < 1.8]
            #            rdchi_min_mean = (Mgrp2.lmfit_redchi.nsmallest(int(len(Mgrp2)/2))).mean()+2*Mgrp2.lmfit_redchi.nsmallest(int(np.round(0.7*len(Mgrp2)))).std()
            error_col = "lmfit_MSE"

            error_min_mean = (
                Mgrp2[error_col].nsmallest(int(len(Mgrp2) / 2))
            ).mean() + 2 * Mgrp2[error_col].nsmallest(
                int(np.round(0.7 * len(Mgrp2)))
            ).std()

            #           [Mgrp[var+'_zscore'] > 1.9 for var in var_names]1
            #            lin = linregress(x=Mgrp[EvRHE], y= Mgrp.lmfit_redchi)
            good_fits = Mgrp[
                ((stats.zscore(Mgrp[error_col])) < 1.4)
                & (Mgrp[error_col] < error_min_mean)
                & (Mgrp.Rs > 2)
                & (Mgrp.Rct < 20e3)
                & (Mgrp.index.isin(z_indexes))
            ]
            if good_fits.empty:
                s_ch = EIS_pars_all.loc[
                    (EIS_pars_all.SampleID == Mgrp.SampleID.unique()[0])
                    & (EIS_pars_all.PAR_file != pf)
                    & (EIS_pars_all.Model_EEC == mod)
                ]
                gs = s_ch[
                    ((stats.zscore(s_ch[error_col])) < 1.4)
                    & (s_ch[error_col] < error_min_mean)
                    & (s_ch.Rs > 2)
                    & (s_ch.Rct < 20e3)
                    & (s_ch.nAd > 0.1)
                ]
                if not gs.empty:
                    good_fits = gs

            var_min = {var + "_min": good_fits[var].min() for var in var_names}
            good_fits = good_fits.assign(**var_min)
            var_max = {var + "_max": good_fits[var].max() for var in var_names}
            good_fits = good_fits.assign(**var_max)

            bad_fits = Mgrp.loc[~Mgrp.index.isin(good_fits.index.values)]

            bad_fits = bad_fits.assign(**var_max)
            bad_fits = bad_fits.assign(**var_min)

            if not bad_fits.empty:
                for (Ev, RPM), Egrp in bad_fits.groupby([EvRHE, "RPM_DAC"]):

                    good_fits["E_diff"] = np.abs(good_fits["E_RHE"] - Ev)
                    top_Ev = int(len(good_fits) / 4 + 3)
                    gf_topE = good_fits.sort_values("E_diff").head(top_Ev)
                    Egrp = Egrp.assign(
                        **{var: gf_topE[var].mean() for var in var_names}
                    )
                    #                    bad_fit_guesses = good_fits.mean().drop(['E_AppV_RHE','E_RHE','RPM_DAC']).dropna().append(pd.Series({'PAR_file' : str(pf), 'Model_EEC' : mod, EvRHE : Ev, 'RPM_DAC' : RPM}))
                    _bad_fit_suggestions.append(Egrp)

            all_good_lst.append(good_fits)
            all_bad_lst.append(bad_fits)

            eischeck_plot = 0
            if eischeck_plot:
                for par in [error_col]:
                    #                    par = 'Rct'
                    fig, ax = plt.subplots()
                    good_fits.plot(x=EvRHE, y=par, c="g", ax=ax)
                    bad_fits.plot(
                        x=EvRHE, y=par, c="r", ax=ax, kind="scatter", marker="x", s=80
                    )
                    ax.set_title(savestem)
                    ax.set_yscale("log")
                    ax.autoscale(True)
                    ax.legend()
                    plt.savefig(
                        PDD_eischeck_figs.joinpath(f"{savestem}.png"),
                        bbox_inches="tight",
                    )
                    plt.close()

    all_good_fits = pd.concat(all_good_lst, sort=False).drop_duplicates()
    all_bad_fits = pd.concat(all_bad_lst, sort=False).drop_duplicates()
    all_bad_fits_suggestions = pd.concat(
        _bad_fit_suggestions, sort=False
    ).drop_duplicates()
    #    .dropna(axis=1, how='all')
    #    print(f'Good fits: {len(all_good_fits)}, bad fits: {len(all_bad_fits)}, diff {len(all_good_fits) - len(all_bad_fits)}')
    bak_dict = {}
    for (all_fits_df, filepath) in [
        (all_good_fits, "EIS_recheck_good_fits"),
        (all_bad_fits, "EIS_recheck_bad_fits"),
        (all_bad_fits_suggestions, "EIS_recheck_bad_fits_suggestions"),
    ]:
        bak_dict.update(check_file_and_backup(all_fits_df, filepath, PDD_eischeck))
    print(
        f"Good fits: {len(all_good_fits)}, bad fits: {len(all_bad_fits)}, diff {len(all_good_fits) - len(all_bad_fits)}"
    )
    print(f"{bak_dict}")


#    .to_excel(PDD_eischeck.joinpath('EIS_recheck_good_fits.xlsx'))
#    .to_excel(PDD_eischeck.joinpath())
#    all_bad_fits_suggestions.to_excel(PDD_eischeck.joinpath())
#    for (pf,mod),Mgrp in EIS_pars_all.groupby(['PAR_file','Model_EEC']):
#            Mgrp.lmfit_redchi.min(), Mgrp.lmfit_redchi.std()
def select_DRT_samples(ORR_acid_no_mod):

    # Reads DRT fitting results from EIS_pars, loops over SourceFilename and loads json
    _DRT_collect = {}
    for sf, sfgr in EIS_pars.groupby("SourceFilename"):
        #        sf,sfgr
        _DRT1 = sf.parent.joinpath(
            "GP_DRT", sf.stem.split("_pars")[0] + "_GP_DRT_res_params.json"
        )
        if _DRT1.is_file():
            _sfdrt = _DRT1
            _drt_json = json.loads(_sfdrt.read_text())
            _DRT_collect.update({sf: _drt_json})
        else:
            print(f"No file: {_DRT1}")

    DRT_pars = pd.DataFrame(data=_DRT_collect.values(), index=_DRT_collect.keys())
    DRT_pars = DRT_pars.rename(
        columns={
            i: "GP_DRT_" + i for i in DRT_pars.columns if not i.startswith("GP_DRT")
        }
    )

    DRT_pars_columns = DRT_pars.columns

    DRT_pars.index.rename("SourceFilename", inplace=True)
    DRT_pars = DRT_pars.reset_index()

    EIS_DRT_pars = pd.merge(EIS_pars, DRT_pars, on="SourceFilename", how="left")
    EIS_DRT_pars = EIS_DRT_pars.loc[
        (EIS_DRT_pars["GP_DRT_success"] == True)
        | (EIS_DRT_pars["GP_DRT_success"] == False)
    ]
    #

    for (msg, Elec), Elgr in EIS_DRT_pars.groupby(["GP_DRT_success", "Electrolyte"]):
        if not Elgr.empty:
            for i in ["sigma_n", "sigma_f", "ell"]:
                fig, ax = plt.subplots()
                #            Elgr.plot(x='BET_cat_agg', y= 'GP_DRT_'+i, kind='scatter', c='E_RHE',colormap='rainbow_r',ax=ax)
                Elgr.plot(
                    x="E_RHE",
                    y="GP_DRT_" + i,
                    kind="scatter",
                    c="BET_cat_agg",
                    colormap="rainbow_r",
                    ax=ax,
                )
                ax.set_xlim = (0.6, 0.9)
                ax.set_title(f"{Elec}, {msg}")
                #            ps = plotting.eisplot(i)
                #            ax.set_ylim(ps.ylim)
                #            ax.set_yscale(ps.logyscale)
                plt.show()
                plt.close()
    # Prepares Spectra from DRT pars and read pickles from Z_star filenames
    _lst_spectra = []
    ECexp_cols = [
        "PAR_file",
        "SampleID",
        "postAST",
        "Loading_cm2",
        "Electrolyte",
        "pH",
        "Gas",
        "RPM_DAC",
        "E_RHE",
    ]
    for grnm, sgrp in EIS_DRT_pars.drop_duplicates(
        subset=list(DRT_pars.columns) + ECexp_cols
    ).groupby(ECexp_cols):
        if len(sgrp) == 1:
            #            ['PAR_file','SampleID', 'E_RHE','Gas']
            _DRT_Zstar = pd.read_pickle(sgrp.GP_DRT_DF_Z_star.iloc[0])
            _DRT_meta = sgrp[DRT_pars_columns].to_dict(orient="records")[0]
            _DRT_Zstar = _DRT_Zstar.assign(
                **{**_DRT_meta, **dict(zip(ECexp_cols, grnm))}
            )
            _lst_spectra.append(_DRT_Zstar)
        else:
            print(grnm)
    DRT_spectra = pd.concat(_lst_spectra, sort=False)

    #            DRT_spectra.update({grnm : pd.read_pickle(sgrp.GP_DRT_DF_Z_star.iloc[0])})

    for sID, sgrp in DRT_spectra.groupby("SampleID"):
        sID, sgrp
        for Ev, Egrp in sgrp.groupby("E_RHE"):
            Ev, Egrp
            if len(Egrp.groupby("Gas")) > 1:
                fig, ax = plt.subplots()
                for n, gr in Egrp.groupby("Gas"):
                    gr.plot(
                        x="freq_vec_star", logx=True, y="gamma_vec_star", ax=ax, label=n
                    )
                ax.set_title(f"{sID} at {Ev}")
                plt.show()
                plt.close()
            else:
                print(sID, Ev)

    for ECexps, sgrp in DRT_spectra.loc[
        (DRT_spectra["GP_DRT_success"] == True)
    ].groupby(["SampleID", "Loading_cm2", "RPM", "Electrolyte", "Gas"]):
        sID, sgrp
        for Ev, Egrp in sgrp.groupby("E_RHE"):
            Ev, Egrp
            if len(Egrp.groupby("postAST")) > 1:
                fig, ax = plt.subplots()
                for n, gr in Egrp.groupby("postAST"):
                    gr.plot(
                        x="freq_vec_star", logx=True, y="gamma_vec_star", ax=ax, label=n
                    )
                ax.set_title(f"{ECexps} at {Ev}")
                plt.show()
                plt.close()
            else:
                print(sID, Ev)


#    lowR = ORR_acid_no_mod.query('Rct < 1.5E3 & Rorr < 1.5E3')
#    .plot(x='Rct',y='Rorr',kind='scatter')

# def export_origin():
def export_to_Origin(EIS_pars):

    SeriesIDs = [
        SampleSelection.Series_CB_paper,
        SampleSelection.Series_Porhp_SiO2,
        SampleSelection.Series_Co_PANI,
        SampleSelection.Series_ML_SiO2,
        {"name": "all_EIS"},
    ]
    SeriesID_set = SeriesIDs[1]
    PostDestDir = FindExpFolder("VERSASTAT").PostDir
    PPDN2 = PostDestDir.joinpath("EIS_{0}".format(SeriesID_set["name"]))
    PPDN2.mkdir(parents=True, exist_ok=True)
    #    all_sIDs = [a for i in [i.get('sIDs',[]) for i in SeriesIDs] for a in i]
    #    PPDCV = PPDN2.joinpath('O2_ORR')
    #    PPDCV.mkdir(parents=True, exist_ok=True)
    EIS_models_collection = Model_Collection()
    EIS_models = [i.model.name for i in EIS_models_collection.lmfit_models]
    EIS_parsSIO2 = EIS_pars.loc[
        (EIS_pars.SampleID.isin(SeriesID_set["sIDs"]))
        & (EIS_pars.Model_EEC.isin(EIS_models))
    ]
    #    EIS_parsSIO2 = EIS_parsSIO2.drop(columns=[i for i in EIS_pars.columns if i in SampleCodes.columns and i != 'SampleID'])
    #    N2_CV_index = postOVVout.groupby('Type_output').get_group('N2_CV')
    #    ORR_parsSIO2 = ORR_parsSIO2.loc[ORR_pars.SampleID.isin(SeriesID_set['sIDs'])]
    #    & ORR_pars.Source.str.contains('ExpDir')]
    Acid_slice = EIS_parsSIO2.query('pH < 7  & postAST == "no"').dropna(
        subset=["ECexp"]
    )
    jOS5 = Acid_slice.query('SampleID == "JOS5"')
    #    plotting.eisplot.fitpars
    ORR_acid_no = Acid_slice.loc[
        ((Acid_slice.PAR_date_day == "2019-05-06") & (Acid_slice.RPM_DAC > 1400))
    ]
    ORR_acid_no = ORR_acid_no.dropna(axis=1, how="all")

    ORR_acid_no

    # ORR_acid_no = check_best_model_per_segment(ORR_acid_no) OLD best model check

    #    'ORR_Jkin_calc_KL_data''ORR_Jkin_calc_RRDE''ORR_Jkin_calc_RRDE_Chrono'

    #    ['Model(Singh2015_RQRQ)', 'Model(Singh2015_RQRQR)', 'Model(Bandarenka_2011_RQRQR)',
    #                  'Model(Singh2015_RQRWR)', 'Model(Randles_RQRQ)', 'Model(Singh2015_R3RQ)']
    # model_select = EIS_models.lmfit_models[11].model.name
    #    ignore_cols = ['Unvalid_freqs','linKK_M', 'linKK_mu','ChiSqr','lmfit_redchi']
    #    selpost = postOVVout.loc[postOVVout.PAR_file.isin(list(ORR_acid_no.PAR_file.unique()))]
    PPD_bestmod = PPDN2.joinpath("best_models_select")
    PPD_bestmod.mkdir(parents=True, exist_ok=True)
    # ORR_acid_no_mod = ORR_acid_no.query('Model_EEC == @model_select')
    _ORR_acid_no_pars_best = ORR_acid_no.loc[
        ORR_acid_no.index.isin(ORR_acid_no.best_mod_index.unique())
    ]
    eisplot.plot_par_with_Ev(
        _ORR_acid_no_pars_best, SampleCodes, OriginColor, PPD_bestmod
    )
    #    spectras_comb = pd.read_excel(ORR_acid_no.SpectraFile.unique()[0],index_col=[0])

    for (bn, gas), bgrp in ORR_acid_no.groupby(["PAR_file", "Gas"]):
        newECexp = f"{bgrp.ECexp.unique()[0]}_{gas}.xlsx"
        _specfit_files = bgrp.File_SpecFit.unique()

        spectras_comb_all = pd.concat(
            [
                pd.read_excel(specfile, index_col=[0]).assign(
                    **{"File_ReadFile": specfile}
                )
                for specfile in _specfit_files
                if "_spectrumfit_" in specfile
            ],
            sort=False,
        )
        spectras_comb = spectras_comb_all.dropna(subset=["Model_EEC"])
        spectras_comb_false = spectras_comb_all.loc[spectras_comb_all.Model_EEC.isna()]

        _pars_best = bgrp.loc[bgrp.index.isin(bgrp.best_mod_index.unique())]

        _segbest = [
            (seg, seggrp.best_mod_name.unique()[0])
            for seg, seggrp in bgrp.groupby("Segment #")
        ]
        _segbest_qry = " | ".join(
            [f'(Segment == {i[0]} & Model_EEC_name == "{i[1]}")' for i in _segbest]
        )
        spectras_comb_bestmods = spectras_comb.query(_segbest_qry)
        _spec_bestmod_check = [
            (seg, seggrp.Model_EEC_name.unique())
            for seg, seggrp in spectras_comb_bestmods.groupby("Segment")
        ]

        DestFile = PPD_bestmod.joinpath(f"{Path(newECexp).stem}_bestmods.png")
        # TODO FIX FOR OR MAKE FIT USING PARS!!
        # EEC_models_index(select=mod)[0][1]
        plot_combined_model(spectras_comb_bestmods, bgrp, SampleCodes, DestFile)

        DestFilePars = PPD_bestmod.joinpath("Pars_Ev_bestmodels")
        # eisplot.plot_par_with_Ev(_pars_best,SampleCodes,OriginColor, PPD_bestmod)

    ORR_acid_no.Model_EEC.unique()

    def plot_best_select_N2_O2(ORR_acid_no):
        N2_mod = "Model(RL-TLM(Rct-Qad-W))"
        O2_mod = "Model(RL-TLM-Rct-W-p(Rorr,Qad))"
        PPD_bestmod = PPDN2.joinpath("best_models_O2_N2")
        PPD_bestmod.mkdir(parents=True, exist_ok=True)

        N2_O2_models = ORR_acid_no.loc[
            ((ORR_acid_no.Model_EEC == N2_mod) & (ORR_acid_no.Gas == "N2"))
            | ((ORR_acid_no.Model_EEC == O2_mod) & (ORR_acid_no.Gas == "O2"))
        ]
        eisplot.plot_par_with_Ev(N2_O2_models, SampleCodes, OriginColor, DestFilePars)
        # mgrp = N2_

    _meta = []
    for mod, mgrp in ORR_acid_no.groupby(["Model_EEC"]):
        #        mod,mgrp = model_select, ORR_acid_no_mod
        PPDmod = PPDN2.joinpath(f"{mod}")
        PPDmod.mkdir(parents=True, exist_ok=True)
        DestFilePars = PPDmod.joinpath(f'Pars_Ev_{mod.split("Model")[-1]}')

        if mod in EIS_models:
            _mod_select = [i for i in EIS_models if i == mod][0]

        eisplot.plot_par_with_Ev(mgrp, SampleCodes, OriginColor, DestFilePars)

        for (bn, gas), bgrp in mgrp.groupby(["PAR_file", "Gas"]):

            newECexp = f"{bgrp.ECexp.unique()[0]}_{gas}.xlsx"
            _specfit_files = bgrp.File_SpecFit.unique()

            spectras_comb_all = pd.concat(
                [
                    pd.read_excel(specfile, index_col=[0]).assign(
                        **{"File_ReadFile": specfile}
                    )
                    for specfile in _specfit_files
                    if "_spectrumfit_" in specfile
                ],
                sort=False,
            )
            spectras_comb = spectras_comb_all.dropna(subset=["Model_EEC"])
            spectras_comb_false = spectras_comb_all.loc[
                spectras_comb_all.Model_EEC.isna()
            ]

            #            for bng_file in bgrp.File_SpecFit.unique():
            #                spf = pd.read_excel(bng_file,index_col=[0])
            #                spf.columns
            ##            pd.read_excel(bgrp.File_SpecFit.unique()[0],index_col=[0])
            #            if not 'Model_EEC'  in spectras_comb.columns:
            #                spectras_comb = pd.concat([pd.read_excel(specfile,index_col=[0]) for specfile in bgrp.File_SpecRaw.unique()])

            #            [(n,gr.File_ReadFile.unique()) for n,gr in spec_comb_mod.groupby([EvRHE])]
            DestFile = PPDmod.joinpath(
                f'{Path(newECexp).stem}_{mod.split("Model")[-1]}.png'
            )
            # TODO FIX FOR OR MAKE FIT USING PARS!!
            # EEC_models_index(select=mod)[0][1]
            plot_combined_model(spec_comb_mod, bgrp, SampleCodes, DestFile)

            bgrp.to_excel(PPDmod.joinpath(Path(f"Pars_Ev_{newECexp}")))
            #            pars = plotting.eisplot.parlst + plotting.eisplot.extrapars + ignore_cols
            outlst = []
            uniqcols = {
                i: bgrp[i].unique()[0] for i in bgrp.columns if bgrp[i].nunique() == 1
            }
            parcols = [i for i in bgrp.columns if i not in uniqcols.keys()]
            for Ev, Egrp in bgrp.groupby("E_RHE"):

                parcols_Ev = [f"{i}_{Ev*1000:.0f}" for i in parcols if not "E_RHE" in i]
                Egrp_new = Egrp.rename(columns=dict(zip(parcols, parcols_Ev)))
                Egrp_new = Egrp_new.assign(
                    **{"basename": bn, "Gas": gas, "SampleID": uniqcols.get("SampleID")}
                )
                Egrp_new = Egrp_new.set_index(["basename", "Gas", "SampleID"])
                outlst.append(Egrp_new[parcols_Ev])
            basegrp_Ev_cols = pd.concat(outlst, ignore_index=False, axis=1)
            uniqcols.update({"OriginDestFile": PPDmod.joinpath(newECexp)})
            _meta.append(uniqcols)
            basegrp_Ev_cols.to_excel(PPDmod.joinpath(newECexp))
        meta = pd.concat([pd.DataFrame(i, index=[0]) for i in _meta])
        meta.to_excel(PPDmod.joinpath("meta_data_EIS_origin.xlsx"))


#    'EC_Porph_SiO2_0.1MH2SO4\EIS_Porph_SiO2'
def _flatten_columns(df):
    _c_lvl0, _c_lvl1 = df.columns.get_level_values(0), df.columns.get_level_values(1)
    df.columns = [
        f"{i[0]}_{i[1]}" if i[1] else f"{i[0]}" for i in zip(_c_lvl0, _c_lvl1)
    ]
    return df


def ML_SiO2_series_pars(EIS_pars):
    _sIDs = SampleSelection.Series.loc["ML_SiO2"].sIDs

    EIS_best_mods = EIS_pars.loc[
        EIS_pars.index.isin([i for i in EIS_pars.best_mod_n.unique() if not pd.isna(i)])
    ]

    PDD = unique_destdir("ML_SiO2")

    MLsio2 = EIS_pars_mod.loc[EIS_pars_mod.SampleID.isin(_sIDs)]
    MLgrpby = MLsio2.groupby(EIS_selection.ECuniq_cols[2:] + ["postAST"])
    [(n, gr.SampleID.unique()) for n, gr in MLgrpby if gr.SampleID.nunique() > 2]

    _grps = [
        ("0.1MH2SO4", 0.7, "O2", 1500.0, "no"),
        ("0.1MH2SO4", 0.5, "N2", 0.0, "no"),
    ]

    for _selgrp in _grps:

        MLgrp = MLgrpby.get_group(_selgrp)

        MLgrp.corr().stack()
        _corr = MLgrp.corr().loc["SizeNP_char"].dropna()
        _corr.name = "corr"
        _abs = _corr.abs()
        _abs.name = "corr_abs"
        _corr = pd.concat([_corr, _abs], axis=1).sort_values(
            "corr_abs", ascending=False
        )
        for n, r in _corr.loc[_corr["corr_abs"] > 0.5].iterrows():
            MLgrp.plot(y=f"{n}", x="SizeNP_char", kind="scatter")
        eisplot.plot_par_with_Ev(
            MLgrp, SampleCodes, OriginColor, PDD, xlabel="SizeNP_char", force=True
        )


def corr_lintangent(EIS_pars):
    PDD = unique_destdir("lintang")
    EIS_pars_mod = select_models_EISpars(EIS_pars)
    Epm_grpby = EIS_pars_mod.groupby(
        EIS_selection.ECuniq_cols[2:] + ["postAST", "Loading_cm2"]
    )
    _grps = [(n, len(gr)) for n, gr in Epm_grpby if gr.SampleID.nunique() > 2]
    _lint = [i for i in EIS_pars_mod.columns if "lintangent" in i]
    _flts = [i for i in EIS_pars_mod.columns if "float" in str(EIS_pars[i].dtype)]
    for ngrp in _grps:
        ngrp = ("0.1MH2SO4", 0.8, "N2", 1500.0, "no", 0.379)
        Egrp = Epm_grpby.get_group(ngrp)
        Egrp.corr().stack()
        _corr = (
            Egrp[_lint + _flts]
            .corr()
            .loc[_lint, [i for i in _flts if i not in _lint]]
            .dropna(axis=1, how="all")
            .drop_duplicates()
        )
        test = _corr.mask(_corr.abs() < 0.5)

        _grpcorr = Egrp[_lint + _flts].corr()
        grptest = _grpcorr.mask(_grpcorr.abs() < 0.5)
        _corr.name = "corr"
        _abs = _corr.abs()
        _abs.name = "corr_abs"
        _corr.loc["_lintangent_uf_freq_mean"]
        _x = [
            "_lintangent_uf_freq_mean",
            "_lintangent_lf_freq_mean",
            "_lintangent_hf_freq_mean",
        ]
        _ylbls = [
            "tBoer_extSA",
            "BET_cat_agg_char",
            "C_content",
            "N_content",
            "BET_RAW_calc_SA m2/g",
        ]
        _corr = pd.concat([_corr, _abs], axis=1).sort_values(
            "corr_abs", ascending=False
        )
        Egrp.plot(x=_x[1], y="N_content", kind="scatter")
        for xy in _ylbls:
            xyPD = unique_destdir(f"{PDD.name}/{xy}")
            eisplot.plot_par_with_Ev(
                Egrp, SampleCodes, OriginColor, xyPD, xlabel=xy, force=True, ylabels=[]
            )


def check_only_best_models(EIS_pars):
    EIS_best_mods = EIS_pars.loc[
        EIS_pars.index.isin([i for i in EIS_pars.best_mod_n.unique() if not pd.isna(i)])
    ]
    PPD_bestmod = PPDN2.joinpath("best_models_only")
    PPD_bestmod.mkdir(parents=True, exist_ok=True)
    DestFilePars = PPD_bestmod.joinpath("Pars_Ev_bestmodels")

    EIS_selection.ECuniq_cols[1:]
    ["SeriesID", "Loading_cm2", "Electrolyte", "RPM_DAC"]
    for (serID, _L, _E, _rpm), sergrp in EIS_best_mods.groupby(
        ["SeriesID", "Loading_cm2", "Electrolyte", "RPM_DAC"]
    ):
        _ndest = DestFilePars.joinpath(f"{serID}/{_E}/{_L}/{_rpm}")
        _ndest.mkdir(parents=True, exist_ok=True)
        eisplot.plot_par_with_Ev(sergrp, SampleCodes, OriginColor, _ndest)

    for (bn, gas), bgrp in EIS_best_mods.groupby(["PAR_file", "Gas"]):
        bgrp = bgrp.assign(**{"segm": bgrp["Segment #"].astype(int)})

        _segmod = bgrp[["segm", "Model_EEC"]].T.to_dict()
        _qry = [
            (f'( (segm == {int(val["segm"])}) & (Model_EEC == "{val["Model_EEC"]}")) ')
            for k, val in _segmod.items()
        ]
        _finalq = " | ".join(_qry)
        newECexp = f"{bgrp.ECexp.unique()[0]}_{gas}.xlsx"
        _dp = PPD_bestmod.joinpath(f"{Path(newECexp).stem}")
        _dp.mkdir(parents=True, exist_ok=True)
        DestFile = _dp.joinpath(f"{Path(bn).stem}.png")
        if not any(DestFile.stem in i.stem for i in DestFile.parent.rglob("*.png")):
            _specfit_files = bgrp.File_SpecFit.unique()
            spectras_comb_all = pd.concat(
                [
                    pd.read_excel(specfile, index_col=[0]).assign(
                        **{"File_ReadFile": specfile}
                    )
                    for specfile in _specfit_files
                    if "_spectrumfit_" in specfile
                ],
                sort=False,
            )
            spectras_comb = spectras_comb_all.dropna(subset=["Model_EEC"])
            spectras_comb = spectras_comb.assign(
                **{"segm": spectras_comb["Segment #"].astype(int)}
            )
            # spectras_comb_false = spectras_comb_all.loc[spectras_comb_all.Model_EEC.isna()]
            spectras_comb_qry = spectras_comb.query(_finalq)
            # TODO FIX FOR OR MAKE FIT USING PARS!!
            plot_combined_model(spectras_comb_qry, bgrp, SampleCodes, DestFile)

        eisplot.plot_par_with_Ev(bgrp, SampleCodes, OriginColor, DestFile.parent)
        plt.close()

    for ECn, ECgrp in EIS_best_mods.groupby(EIS_selection.ECuniq_cols[1:]):
        if ECgrp.SampleID.nunique() > 2:
            ECgrp.plot(
                x="BET_cat_agg",
                y="Cdlp",
                title=", ".join(str(i) for i in ECn),
                kind="scatter",
            )


def plot_only_best_models(ORR_acid_no):
    SeriesIDs = [
        SampleSelection.Series_CB_paper,
        SampleSelection.Series_Porhp_SiO2,
        SampleSelection.Series_Co_PANI,
        SampleSelection.Series_ML_SiO2,
        {"name": "all_EIS"},
    ]
    SeriesID_set = SeriesIDs[1]
    PostDestDir = FindExpFolder("VERSASTAT").PostDir
    PPDN2 = PostDestDir.joinpath("EIS_{0}".format(SeriesID_set["name"]))
    PPDN2.mkdir(parents=True, exist_ok=True)

    EIS_models = Model_Collection()
    EIS_models_lenpars = EIS_models.get_df_models_parameters()
    _ranks = []
    for (bn, gas), bgrp in ORR_acid_no.groupby(["PAR_file", "Gas"]):
        (bn, gas), bgrp
        print((bn, gas))
        # eisplot.plot_par_with_Ev(_mod_select, mgrp,SampleCodes,OriginColor, DestFilePars)
        # .loc[eval(_filter)]\
        bgrpmod = bgrp.groupby(["Model_EEC", "E_RHE"])
        errtype = "lmfit_MSE"
        best_models = (
            bgrpmod[errtype, "lmfit_errpars_sum"]
            .agg(["count", "mean", "std"])
            .sort_values(("lmfit_errpars_sum", "mean"), ascending=True)
        )
        PPDmod = PPDN2.joinpath(f"Best_mod_{gas}")
        PPDmod.mkdir(parents=True, exist_ok=True)

        _checkmods = best_models.loc[
            (best_models[("lmfit_errpars_sum", "mean")] > 0)
            & (best_models[("lmfit_errpars_sum", "mean")] < 1e3)
        ][[(errtype, "mean"), ("lmfit_errpars_sum", "mean")]]
        bgrp[["Model_EEC", "E_RHE", "lmfit_errpars_badpars"]]

        _checkmods = _flatten_columns(_checkmods)
        _checkmods.reset_index(inplace=True)

        _checkmods = pd.merge(
            _checkmods, bgrp[["Model_EEC", "E_RHE", "lmfit_errpars_badpars"]]
        )
        _mod_rank = [
            (n, len(gr), gr[errtype + "_mean"].mean())
            for n, gr in _checkmods.reset_index().groupby("Model_EEC")
        ]

        _mod_rank_best = sorted(_mod_rank, key=lambda x: x[-1])[0][0]
        best_models.loc[best_models[("lmfit_errpars_sum", "mean")] > 0].iloc[0].name

        _mod_select = _mod_rank_best
        _ranks.append([(Path(bn), gas, *i) for i in _mod_rank])

    Mod_Rank = pd.DataFrame(
        [a for i in _ranks for a in i],
        columns=["PAR_file", "Gas", "Model_EEC", "rank_len_gr", f"rank_{errtype}_mean"],
    )
    Best_Mod_Rank = Mod_Rank.loc[
        (Mod_Rank.rank_len_gr > 6) & (Mod_Rank[f"rank_{errtype}_mean"] < 30)
    ]

    Mod_Rank_grps = (
        Mod_Rank.groupby(["Model_EEC", "Gas"])["rank_len_gr", f"rank_{errtype}_mean"]
        .agg(["count", "mean", "sum", "std"])
        .sort_index(level=1)
    )

    _mrge = Mod_Rank_grps.reset_index()
    _mrge = _flatten_columns(_mrge)
    _mrge.columns
    Mod_Ranks = pd.merge(_mrge, EIS_models_lenpars, on="Model_EEC", how="left")

    N2 = Mod_Ranks.query('Gas == "N2" & rank_len_gr_count > 3')

    def plot_ranks():
        TT = bgrpmod.get_group(_mod_select)
        DestFilePars = PPDmod.joinpath(f'Pars_Ev_{_mod_select.split("Model")[-1]}')
        _mod_select_inst = [
            i for i in EIS_models.lmfit_models if i.model.name == _mod_select
        ][0]
        eisplot.plot_par_with_Ev(
            _mod_select_inst, bgrp, SampleCodes, OriginColor, DestFilePars
        )

        bgrp.groupby("Model_EEC").agg()

        newECexp = f"{bgrp.ECexp.unique()[0]}_{gas}.xlsx"
        _specfit_files = bgrp.File_SpecFit.unique()

        spectras_comb_all = pd.concat(
            [
                pd.read_excel(specfile, index_col=[0]).assign(
                    **{"File_ReadFile": specfile}
                )
                for specfile in _specfit_files
                if "_spectrumfit_" in specfile
            ],
            sort=False,
        )
        spectras_comb = spectras_comb_all.dropna(subset=["Model_EEC"])
        spectras_comb_false = spectras_comb_all.loc[spectras_comb_all.Model_EEC.isna()]

        #            for bng_file in bgrp.File_SpecFit.unique():
        #                spf = pd.read_excel(bng_file,index_col=[0])
        #                spf.columns
        ##            pd.read_excel(bgrp.File_SpecFit.unique()[0],index_col=[0])
        #            if not 'Model_EEC'  in spectras_comb.columns:
        #                spectras_comb = pd.concat([pd.read_excel(specfile,index_col=[0]) for specfile in bgrp.File_SpecRaw.unique()])
        #
        spec_comb_mod = (
            spectras_comb.groupby(["Model_EEC"])
            .get_group(_mod_select)
            .sort_values(by=[EvRHE, "Frequency(Hz)"])
        )
        #            [(n,gr.File_ReadFile.unique()) for n,gr in spec_comb_mod.groupby([EvRHE])]
        DestFile = PPDmod.joinpath(
            f'{Path(newECexp).stem}_{mod.split("Model")[-1]}.png'
        )
        # TODO FIX FOR OR MAKE FIT USING PARS!!
        # EEC_models_index(select=mod)[0][1]

        plot_combined_model(spec_comb_mod, bgrp, SampleCodes, DestFile)

        bgrp.to_excel(PPDmod.joinpath(Path(f"Pars_Ev_{newECexp}")))
        #            pars = plotting.eisplot.parlst + plotting.eisplot.extrapars + ignore_cols
        outlst = []
        uniqcols = {
            i: bgrp[i].unique()[0] for i in bgrp.columns if bgrp[i].nunique() == 1
        }
        parcols = [i for i in bgrp.columns if i not in uniqcols.keys()]
        for Ev, Egrp in bgrp.groupby("E_RHE"):

            parcols_Ev = [f"{i}_{Ev*1000:.0f}" for i in parcols if not "E_RHE" in i]
            Egrp_new = Egrp.rename(columns=dict(zip(parcols, parcols_Ev)))
            Egrp_new = Egrp_new.assign(
                **{"basename": bn, "Gas": gas, "SampleID": uniqcols.get("SampleID")}
            )
            Egrp_new = Egrp_new.set_index(["basename", "Gas", "SampleID"])
            outlst.append(Egrp_new[parcols_Ev])
        basegrp_Ev_cols = pd.concat(outlst, ignore_index=False, axis=1)
        uniqcols.update({"OriginDestFile": PPDmod.joinpath(newECexp)})
        _meta.append(uniqcols)
        basegrp_Ev_cols.to_excel(PPDmod.joinpath(newECexp))

    meta = pd.concat([pd.DataFrame(i, index=[0]) for i in _meta])
    meta.to_excel(PPDmod.joinpath("meta_data_EIS_origin.xlsx"))


# def check_best_model_per_segment(ORR_acid_no):
#     EIS_models = Model_Collection()
#     _errs = []
#     for (pf, nseg),grp in ORR_acid_no.groupby(['PAR_file', 'Segment #']):
#         (pf, nseg),grp
#         for modn, modgrp in grp.groupby('Model_EEC'):
#             modn,modgrp
#             modpars = EIS_models.modpars.get(modn)
#             if modpars:
#                 modpars_errs = [i+'_stderr' for i in modpars]
#                 moderr = modgrp[modpars_errs].iloc[0].T
#                 _badpars = ', '.join([key.split('_stderr')[0] for key,val in moderr.to_dict().items() if val > 1E3])
#                 _errs.append((pf,nseg, modn, moderr.sum(), moderr.mean(), moderr.std(), _badpars ))

#     _prfx = 'lmfit_errpars'
#     EIS_pars_errsum = pd.DataFrame(_errs,columns =['PAR_file', 'Segment #','Model_EEC',f'{_prfx}_sum', f'{_prfx}_mean', f'{_prfx}_std',f'{_prfx}_badpars'])
#     ORR_acid_no_err = pd.merge(ORR_acid_no,EIS_pars_errsum, how='left')
#     return ORR_acid_no_err


def plot_combined_model(spectras_comb, bgrp, SampleCodes, DestFile):

    for ZY in ["Z", "Y", "-Zangle"]:
        DestZY = DestFile.parent.joinpath(f"{DestFile.stem}_{ZY}.png")
        eisplot.PlotCombinedEIS(spectras_comb, bgrp, SampleCodes, DestZY, xEIS=ZY)


def plot_vars():
    for i in var_names:
        fig, ax = plt.subplots()
        Mgrp.plot(x="E_RHE", y=i, ax=ax)
        z_val = np.abs(stats.zscore(Mgrp[i]))
        for (k, v), z in zip(Mgrp[["E_RHE", i]].iterrows(), z_val):
            ax.annotate(
                np.round(z, 2),
                v,
                xytext=(10, -5),
                textcoords="offset points",
                family="sans-serif",
                fontsize=18,
                color="darkslategrey",
            )


# def plot_pars_with_Ev(bgrp,SampleCodes, DestFilePars):
#    for ZY in ['Z', 'Y', '-Zangle']:
#        DestZY = DestFile.parent.joinpath(f'{DestFile.stem}_{ZY}.png')
#        eisplot.plot_par_with_Ev(bgrp,SampleCodes,DestFilePars)
# plot_combined_model(spectras_comb, bgrp,SampleCodes, DestFile)


def MergeEISandORR():
    EODD = PostEC().DestDir.joinpath("EIS_ORR")
    EODD.mkdir(parents=True, exist_ok=True)
    # TODO 2020.01.16 trying to merge EIS and ORR pars

    EC_midx = ["SampleID", "pH", "Electrolyte", "Loading_cm2", "postAST"]
    EIScols = (
        post_helper.CheckCols(
            SampleSelection.EC_EIS_par_cols
            + ["ECexp", "E_RHE", "SampleID"]
            + ["lmfit_redchi", "Chisqr", "Model_EEC", "Model_index"],
            EIS_pars,
        )
        + post_helper.CheckCols(SampleSelection.EC_exp_cols, ORR_pars)
    )
    #    Cdlcatancols = Cdl_pars_catan[EC_midx+['ECexp','E_RHE']+MakingCorrs.CheckCols(SampleSelection.EC_exp_cols,Cdl_pars_catan)]

    fast_checking_EEC_models = [
        "Model(Singh2015_RQRQR)",
        "Model(Singh2015_RQRWR)",
        "Model(Singh2015_R3RQ)",
        "Model(Bandarenka_2011_RQRQR)",
    ]

    EIS_midx = EIS_pars.set_index(EC_midx).loc[:, EIScols].dropna(axis=1)
    sIDs_used = EIS_midx.index.get_level_values(level="SampleID").unique()
    ORR_midx = ORR_pars.set_index(EC_midx)
    ORR_mx_used = ORR_midx.loc[
        ORR_midx.index.get_level_values(level="SampleID").isin(sIDs_used)
    ]
    Cdl_catan_midx = Cdl_pars_catan.set_index(EC_midx)
    Cdl_mx_used = Cdl_catan_midx.loc[
        Cdl_catan_midx.index.get_level_values(level="SampleID").isin(sIDs_used)
    ]
    #    pd.merge(EIS_midx,ORR_midx,how='inner',left_index=True, right_index=True)
    EIS_midx.join(ORR_midx, lsuffix="_eis", rsuffix="_ORR")
    #    EIS_ORR_pars = pd.merge(EIS_pars, ORR_pars , on=['ECexp'],how='left',suffixes=['_eis','_ORR']) # 2013511
    #    EIS_ORR_pars.to_pickle(EODD.joinpath('EIS_ORR_pars_full.pkl.compress'),compression='xz')
    print(EIS_pars.Model_EEC.unique())
    EIS_ORR_fast = pd.merge(
        EIS_pars.query('Model_EEC == "Model(Singh2015_R3RQ)"'),
        ORR_pars,
        on=["ECexp"],
        how="left",
        suffixes=["_eis", "_ORR"],
    )  # 2013511
    print(EIS_ORR_fast.ECexp.unique())
    EIS_ORR_fast.to_pickle(
        EODD.joinpath("EIS_ORR_pars_mSingRQRQR.pkl.compress"), compression="xz"
    )
    faillst, ORReis_merge_lst = [], []
    merge_Model_EEC = [" Model(Singh2015_RQRQR)", "Model(Singh2015_R3RQ)"]
    merge_Model_EEC_set = merge_Model_EEC[1]

    add_ORR_cols = [
        i
        for i in ORR_pars.columns
        if i not in EIS_pars.columns and i.split("_")[-1] not in ["x", "y"]
    ]
    for exp, ECgr in EIS_pars.query(
        f'(Model_EEC == "{merge_Model_EEC_set}") & (RPM_DAC > 1000)'
    ).groupby("ECexp"):
        exp, ECgr
        ORR_ECexp_slice = ORR_pars.query("ECexp == @exp")

        if not ORR_ECexp_slice.empty:
            # TODO Merge on E_RHE with data from ORR scans and compare then...
            EIS_ECexp = ECgr
            nECexp = exp
            ORR_ECexp = ORR_ECexp_slice
            ORR_ECexp_PARS_1500 = ORR_ECexp.query("RPM > 1000")[add_ORR_cols]
            if len(ORR_ECexp_PARS_1500) == 1:
                for col in ORR_ECexp_PARS_1500.columns:
                    EIS_ECexp[col] = [ORR_ECexp_PARS_1500[col].values[0]] * len(
                        EIS_ECexp
                    )
                ORReis_merge_lst.append(EIS_ECexp)
            else:
                faillst.append([exp, ORR_ECexp_slice, "len_not_1"])
        else:
            faillst.append([exp, ORR_ECexp_slice, "ORR_empty"])
    failed_ECexps = pd.DataFrame(
        [(i[0], i[-1]) for i in faillst], columns=["ECexp", "MergeFailMsg"]
    )
    failed_ECexps_grps = pd.concat([i[1] for i in faillst], sort=False)
    ORReis_merge_raw = pd.concat(ORReis_merge_lst, sort=False)
    if not Cdl_pars_catan.empty:
        mcs1 = [i for i in Cdl_pars_catan.columns if i in ORReis_merge_raw.columns]
        olap = [i for i in Cdl_pars_catan.columns if i in ORReis_merge_raw.columns]
        [
            i
            for i in olap
            if Cdl_pars_catan[i].nunique() == 1 and ORReis_merge_raw[i].nunique() == 1
        ]

        mcs = [
            i
            for i in SampleSelection.EC_exp_cols
            if i in ORReis_merge_raw.columns and i in Cdl_pars_catan.columns
        ]

        EIS_ORR_Cdl_fast = pd.merge(
            ORReis_merge_raw.query(f'Model_EEC == "{merge_Model_EEC_set}"'),
            Cdl_pars_catan,
            on=["ECexp", "E_RHE"],
            how="left",
            suffixes=["", "_cdl"],
        )  # 2013511

    ORReis_merge_raw.to_pickle(EODD.joinpath("EIS_ORR_refit_pars.pkl.compress"))
    ORReis_mergeCB = ORReis_merge_raw.loc[
        ORReis_merge_raw.SampleID.isin(SampleSelection.Series_CB_paper["sIDs"])
    ]
    ORReis_merge = EIS_ORR_Cdl_fast
    RedChiSq_limit = (
        ORReis_merge.query("Rs > 1").lmfit_redchi.mean()
        + 1 * ORReis_merge.query("Rs > 1").lmfit_redchi.std()
    )

    ORReis_neat = ORReis_merge.query(
        "lmfit_redchi < @RedChiSq_limit & Rs > 2 & Rct < 9E05"
    )
    ORReis_neat.query("E_RHE < 0.8 & E_RHE > 0.74 & Rct < 9E05").plot(
        x="E_onset",
        y="Rct",
        kind="scatter",
        ylim=(0.1, 1e5),
        xlim=(0.5, 1),
        logy=True,
        logx=0,
        c="pH",
        colormap="viridis",
    )
    ORReis_neat.query("E_RHE < 0.8 & E_RHE > 0.74 & Rct < 9E05").plot(
        x="E_onset",
        y="Cdl_an",
        kind="scatter",
        xlim=(0.5, 1),
        logy=True,
        logx=0,
        c="pH_cdl",
        colormap="viridis",
    )
    ORReis_neat.SampleID.unique()
    EIS_pars.SampleID.unique()
    #    .loc[ORReis_neat.SampleID.isin(SampleSelection.Series_CB_paper['sIDs']+SampleSelection.Series__paper['sIDs']]
    sIDslice = (
        SampleSelection.Series_CB_paper["sIDs"]
        + SampleSelection.Series_Porhp_SiO2["sIDs"]
    )
    Ekin = (
        ORReis_neat.loc[ORReis_neat.SampleID.isin(sIDslice)]
        .query('E_RHE <= 0.9 & E_RHE >= 0.6 & Gas == "O2" & SampleID != "JOS5" ')
        .dropna(subset=["SeriesID"])
    )
    Ekin_grp = Ekin.groupby(["pH", "E_RHE"])

    ORReis_neat.query("E_RHE < 0.8 & E_RHE > 0.74 & Rct < 9E05").plot(
        x="Rs",
        y="Jkin_075",
        kind="scatter",
        ylim=(0.1, 100),
        xlim=(0.5, 70),
        logy=True,
        logx=0,
        c="BET_cat_agg",
        colormap="rainbow_r",
    )
    ORReis_neat.query("E_RHE < 0.8 & E_RHE > 0.74 & Rct < 9E05").plot(
        x="Cdl_cat",
        y="Rs",
        kind="scatter",
        logy=True,
        logx=0,
        c="BET_cat_agg",
        colormap="rainbow_r",
    )
    for i in [
        i for i in SampleSelection.EC_EIS_par_cols if not "lmfit_redchi" in i
    ] + SampleSelection.EC_ORR_kin_par_cols:
        EKinsl = Ekin.query(
            "E_RHE < 0.76 & E_RHE > 0.59 & Rct < 9E05 & Rorr < 1E09 & Qad < 35E-3 & Cdlp < 0.070"
        )
        for Elec, Elgr in EKinsl.groupby("Electrolyte"):
            fig, ax = plt.subplots()

            Elgr.plot(
                x="BET_cat_agg",
                y=i,
                kind="scatter",
                c="E_RHE",
                colormap="rainbow_r",
                ax=ax,
            )
            ax.set_xlim = (0.6, 0.9)
            ax.set_title(Elec)
            ps = plotting.eisplot(i)
            ax.set_ylim(ps.ylim)
            ax.set_yscale(ps.logyscale)
            plt.show()
            plt.close()

    def one_sample():
        DW25 = ORReis_neat.query('ECexp == "DW25_1.0_0.1MH2SO4_0.379_no_2019-01-28"')

    def count_corr_occurences(grEcorrstack, *args):
        lst = []
        for n, val in grEcorrstack.iteritems():
            t, c = 0, ""
            if any([(i in n) for i in eisplot.parlst]):
                t += 1
                c += "eis"
            if any([(i in n) for i in SampleSelection.EC_ORR_kin_par_cols]):
                t += 1
                c += "orr"
            if any([(i in n) for i in SampleSelection.Character_xcols]):
                c += "Xchar"

            odct = dict(
                zip(
                    ["corr_idx", "val", "ocurrc", "occur_s", "pH", "E_RHE"],
                    [n, val, t, c, args[0], args[1]],
                )
            )
            lst.append(odct)
        return lst

    Ekin_occur_lst = []
    corr_method, corr_cutoff = "spearman", 0.1  # or default is pearson spearman
    # TODO conti
    for E, grE in Ekin_grp:
        E, grE
        for corr_m in ["pearson", "spearman"]:
            grEcorrstack = grE.corr(method=corr_m).stack().sort_values()
            grEcorrstack = grEcorrstack[np.abs(grEcorrstack) < 0.99]
            occur_counts = count_corr_occurences(grEcorrstack, *E)
            occur_DF = pd.DataFrame(occur_counts)
            occur_DF["method"] = corr_m
            Ekin_occur_lst.append(occur_DF)
    Ekin_occur = (
        pd.concat(Ekin_occur_lst, sort=True)
        .drop_duplicates(subset=["val", "pH", "E_RHE"])
        .sort_values(by="val")
    )
    Ekin_occur.query("ocurrc == 2")
    grEtop = pd.concat(
        [
            Ekin_occur.query("ocurrc == 2").head(25),
            Ekin_occur.query("ocurrc == 2").tail(25),
        ],
        sort=False,
    )
    for pH, grppH in Ekin_occur.query("ocurrc == 2").groupby("pH"):
        grEtop = pd.concat(
            [grppH.query("ocurrc == 2").head(10), grppH.query("ocurrc == 2").tail(10)],
            sort=False,
        )
        PPDEISbest_overall_pH = EODD.joinpath(
            "MergeEIS_ORR_top_corrs_{0}_overall_{1}_{2}".format(
                "both", pH, merge_Model_EEC_set
            )
        )
        PPDEISbest_overall_pH.mkdir(parents=True, exist_ok=True)
        for n, row in grEtop.iterrows():
            #        n,row
            grE = Ekin_grp.get_group(tuple(row[["pH", "E_RHE"]]))
            plot_pd_SampleIDs(
                grE,
                row.corr_idx[0],
                row.corr_idx[1],
                row.val,
                PPDEISbest_overall_pH.joinpath(str(pH)),
                corr_method=row["method"],
            )
    """ From the correlations in pH 13 it seems that nDL strongly correlates with the TSb,j_diff lim, so seems relating
        to the exchange current density
        ORR_pars.plot(x='TSb_l',y='E_onset',kind='scatter',xlim=(0,1),c='pH',colormap='viridis')
        """

    PPDEISbest = EODD.joinpath("MergeEIS_ORR_top_corrs_{0}".format(corr_method))
    PPDEISbest.mkdir(parents=True, exist_ok=True)
    for E, grE in Ekin_grp:
        E, grE

        grEoccur = Ekin_occur.query(
            "ocurrc == 2 & pH == @E[0] & E_RHE == @E[1] "
        ).sort_values(by=["val"])
        rct_slice = [i for i in grEoccur.corr_idx.values if "Rct" in i]
        rct_occur = grEoccur.loc[grEoccur.corr_idx.isin(rct_slice)]
        grEtop = pd.concat([grEoccur.head(5), grEoccur.tail(5)], sort=False)
        for n, tr in grEtop.iterrows():
            #            Gas_set,pH_set,Erhe_set,corr_val = tr[1].Gas,tr[1].pH,tr[1].E_RHE,tr[1][corr_method]
            plot_pd_SampleIDs(
                grE,
                tr.corr_idx[0],
                tr.corr_idx[1],
                tr.val,
                PPDEISbest.joinpath(f"{E[0]}_{E[1]}"),
            )
            plt.close()
        for n, tr in rct_occur.iterrows():
            #            Gas_set,pH_set,Erhe_set,corr_val = tr[1].Gas,tr[1].pH,tr[1].E_RHE,tr[1][corr_method]
            plot_pd_SampleIDs(
                grE,
                tr.corr_idx[0],
                tr.corr_idx[1],
                tr.val,
                PPDEISbest.joinpath(f"Rct/{E[0]}_{E[1]}"),
            )
            plt.close()

    #    .loc[ORReis_merge_raw.SampleID.isin(SampleSelection.Series_CB_paper['sIDs'])]
    for E, grE in Ekin_grp:
        E, grE

        grEoccur = Ekin_occur.query("ocurrc == 2 & pH == @E[0] & E_RHE == @E[1] ")
        grEtop = pd.concat([grEoccur.head(5), grEoccur.tail(5)], sort=False)
        for n, tr in grEtop.iterrows():
            #            Gas_set,pH_set,Erhe_set,corr_val = tr[1].Gas,tr[1].pH,tr[1].E_RHE,tr[1][corr_method]
            plot_pd_SampleIDs(
                grE,
                tr.corr_idx[0],
                tr.corr_idx[1],
                tr.val,
                PPDEISbest.joinpath(f"{E[0]}_{E[1]}"),
            )
            plt.close()

    Ekincorr = Ekin.corr().stack()[SampleSelection.EC_ORR_kin_par_cols + eisplot.parlst]
    [eisplot.parlst]
    Ekincorr = Ekincorr.loc[(np.abs(Ekincorr) < 1)].sort_values()

    tpcb_top = pd.concat([Ekincorr.head(50), Ekincorr.tail(50)], sort=corr_method)

    for tr in tpcb_top.iteritems():
        #            Gas_set,pH_set,Erhe_set,corr_val = tr[1].Gas,tr[1].pH,tr[1].E_RHE,tr[1][corr_method]
        xc, yc = tr[0][0], tr[0][1]
        plot_pd_SampleIDs(Ekin, xc, yc, tr[1], PPDEISbest)

    for par in eisplot.parlst[::]:
        ylims, logys, lin = eisplot(par).pars_plot_lims()
        fig, ax = plt.subplots()
        Ekin.plot(
            x="E_RHE",
            y=par,
            kind="scatter",
            logy=logys,
            ylim=ylims,
            logx=0,
            c="Jkin_075",
            colormap="viridis",
            ax=ax,
        )
        ax.set_xlabel("E_RHE")


#        ORReis_neat.query('E_RHE < 0.81 & E_RHE > 0.59').plot(x='Jkin_075', y=par, kind='scatter',logy=0,ylim=(0,5000),logx=True,c='E_RHE',colormap='viridis')

#        ORReis_merge.groupby('E_RHE'):
# TODO


def EIS_CBpaper():
    EIS_exp_samples = [
        i
        for i in sorted(EIS_pars.SampleID.unique())
        if i in SampleSelection.Series_CB_paper["sIDs"]
    ]
    print(", ".join(EIS_exp_samples))
    print(EIS_pars.Model_EEC.unique())
    CBeis = post_helper.serie_model(EIS_pars, EIS_exp_samples, "Model(Singh2015_R3RQ)")
    CBcdl = post_helper.serie_model(Cdl_pars_catan, EIS_exp_samples)
    CBorr = post_helper.serie_model(ORR_pars, EIS_exp_samples)
    topcorrs, topcorr_best, topcorr_best_Ev = Correlations.MakeTopCorrs(
        CBeis, "pearson", PostEC().DestDir.joinpath("EIS_corr_CB")
    )
    topcorr_best.sort_values("SampleScore", ascending=False).head(30)


def EIS_loading():
    ELoadDD = FindExpFolder("VERSASTAT").PostDir.joinpath("Loading_series_EIS_ORR")
    ELoadDD.mkdir(parents=True, exist_ok=True)
    _filter = EIS_selection.filter
    _models = EIS_selection.fast_checking_EEC_models

    # EIS_pars = EIS_pars.assign(**{'R3_kin' : 1/EIS_pars.R3.values, 'Rorr_kin' : 1/EIS_pars.Rorr.values,  'Qad+Cdpl+Q3' : EIS_pars[['Qad','Cdlp','Q3']].sum(axis=1) })
    EISloading = post_helper.loading_series(
        EIS_pars.loc[EIS_pars.lmfit_message.str.contains("satisfied")].query(_filter)
    )
    filter_PARfls = "JOS15_0.3_0.5MH2SO4_0.379_no_2019-03-18"
    EISloading_fltr = EISloading.loc[
        EISloading.ECexp.str.contains(filter_PARfls) == False
    ].query(_filter)
    EIS_exp_samples = [
        i
        for i in sorted(EIS_pars.SampleID.unique())
        if i in SampleSelection.Series_CB_paper["sIDs"]
    ]
    EISloadmod = post_helper.serie_model(
        EISloading, EIS_exp_samples, "Model(R0-L0-p(R1-Ws1,CPE1)-C2)"
    ).query(_filter)

    loadgrp_cols = EIS_selection.loadgrp_cols
    llst = []
    for i in EISloadmod.lmfit_var_names.unique()[0].split(", "):
        EKinsl = EISloading.query(
            "Rct < 9E05 & Rorr < 1E09 & Qad < 35E-3 & Cdlp < 0.070"
        )
        Ekingrp = EKinsl.groupby(loadgrp_cols)
        for Elec, Elgr in Ekingrp:
            if Elgr.Loading_cm2.nunique() > 2:
                loading_corr = (
                    Elgr[["Loading_cm2", i]]
                    .corr(method="pearson")
                    .loc["Loading_cm2", i]
                )
                llst.append(
                    Elec
                    + (
                        "Loading_cm2",
                        i,
                        np.round(loading_corr, 3),
                        Elgr.lmfit_redchi.mean(),
                    )
                )
    loadcols = loadgrp_cols + ["corr_type", "eispar", "corr", "redchimean"]
    parsloadcorr = (
        pd.DataFrame(llst, columns=loadcols).dropna(subset=["corr"]).sort_values("corr")
    )

    _independent = ["Loading_cm2", "RPM_DAC"][1]
    loadgrp_cols = ["SampleID", "Electrolyte", "Gas", "Model_EEC"]
    _loading_corr = []
    for Elec, Elgr in EISloadmod.groupby(loadgrp_cols):
        for E, Egrp in [
            (n, gr) for n, gr in Elgr.groupby("E_RHE") if gr[_independent].nunique() > 2
        ]:
            _corr = Egrp.corr().stack()[_independent][
                EISloadmod.lmfit_var_names.unique()[0].split(", ")
            ]
            _corr.name = "corrval"
            _corrDF = pd.DataFrame(_corr).assign(
                **{**dict(zip(loadgrp_cols, Elec)), **{"E_RHE": E}}
            )
            _loading_corr.append(_corrDF)
            for par, corrval in _corr.iteritems():
                if np.abs(corrval) > 0.91:
                    Egrp.plot(
                        x=_independent,
                        y=par,
                        c=_independent,
                        kind="scatter",
                        cmap="rainbow",
                        title=f"{par},{E},{Elec},corr:{corrval:.2f}",
                    )
    _loading_corrs = pd.concat([i for i in _loading_corr])

    # fast_checking_EEC_models = ['Model(R0-L0-p(R1-W1,CPE1)-C2)']
    # 'Model(Singh2015_RQRQR)', 'Model(Singh2015_RQRWR)', 'Model(Singh2015_R3RQ)', 'Model(Bandarenka_2011_RQRQR)' ]
    EISloading_fltr_mod = post_helper.serie_model(
        EISloading.query(_filter), EIS_exp_samples, _models[2]
    )
    eflt_grp = EISloading_fltr_mod.groupby(["Gas", "SampleID"])

    for gasn, ggr in EISloading_fltr_mod.groupby(["Gas", "SampleID"]):
        ggr.plot(
            x="E_RHE",
            y="tau",
            kind="scatter",
            c="Loading_cm2",
            colormap="viridis",
            label=gasn,
        )

    eflt_grp_N2 = post_helper.serie_model(
        EISloading_fltr, ["JOS15"], _models[2]
    ).groupby(["Gas", "SampleID"])
    N2JOS15 = eflt_grp_N2.get_group(("N2", "JOS15"))
    N2lbl = N2JOS15[["Gas", "SampleCode"]].iloc[0].to_list()

    cmapset = "rainbow"
    fig, ax1 = plt.subplots()
    N2JOS15.plot(
        x="E_RHE",
        y="Rs",
        ylim=(0, 20),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        title=", ".join(N2lbl),
        ax=ax1,
    )
    plt.savefig(ELoadDD.joinpath("N2JOS15_Rs.png"), dpi=300, bbox_inches="tight")

    fig, ax1 = plt.subplots()
    N2JOS15.plot(
        x="E_RHE",
        y="Qad+Cdlp",
        ylim=(0, 0.035),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        title=", ".join(N2lbl),
        ax=ax1,
    )
    plt.savefig(ELoadDD.joinpath("N2JOS15_QadCdlp.png"), dpi=300, bbox_inches="tight")
    # parlst = ['Rct','Rs','Rorr','Rct_kin','Cdlp','Qad','nDL','nAd','Qad+Cdlp']
    for p in N2JOS15.lmfit_var_names.unique()[0].split(", "):
        fig, ax1 = plt.subplots(figsize=(8, 4))
        ylimset = eisplot(p, yvalues=N2JOS15[p].values).ylim
        N2JOS15.plot(
            x="E_RHE",
            y=p,
            ylim=ylimset,
            kind="scatter",
            c="Loading_cm2",
            colormap=cmapset,
            title=", ".join(N2lbl),
            ax=ax1,
        )
        plt.savefig(ELoadDD.joinpath(f"N2JOS15_{p}.png"), dpi=300, bbox_inches="tight")
    N2JOS15.plot(
        x="E_RHE",
        y="nDL",
        ylim=(0.5, 1),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        label=gasn,
    )
    N2JOS15.plot(
        x="E_RHE",
        y="nAd",
        ylim=(0.5, 1),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        label=gasn,
    )

    eflt_grp_O2 = post_helper.serie_model(
        EISloading_fltr, ["JOS15"], fast_checking_EEC_models[0]
    ).groupby(["Gas", "SampleID"])
    O2JOS15 = eflt_grp_O2.get_group(("O2", "JOS15"))
    O2lbl = O2JOS15[["Gas", "SampleCode"]].iloc[0].to_list()
    for p in O2JOS15.lmfit_var_names.unique()[0].split(", "):
        fig, ax1 = plt.subplots(figsize=(8, 4))
        eisplt = eisplot(p, yvalues=O2JOS15[p].values)
        ylimset = eisplt.ylim
        if "Rct_kin" in p:
            ylimset = (0, 2)
        O2JOS15.plot(
            x="E_RHE",
            y=p,
            ylim=ylimset,
            logy=eisplt.logy,
            kind="scatter",
            c="Loading_cm2",
            colormap=cmapset,
            title=", ".join(O2lbl),
            ax=ax1,
        )
        plt.savefig(ELoadDD.joinpath(f"O2JOS15_{p}.png"), dpi=300, bbox_inches="tight")

    O2JOS15.plot(
        x="E_RHE",
        y="Rs",
        ylim=(0, 20),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        label=", ".join(lbl),
    )
    O2JOS15.plot(
        x="E_RHE",
        y="Qad+Cdlp",
        ylim=(0, 0.05),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        label=gasn,
    )
    O2JOS15.plot(
        x="E_RHE",
        y="nDL",
        ylim=(0.5, 1),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        label=gasn,
    )
    O2JOS15.plot(
        x="E_RHE",
        y="nAd",
        ylim=(0.5, 1),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        label=gasn,
    )

    N2JOS15.EC_exp
    CdlJOS15 = Cdl_pars_catan.loc[Cdl_pars_catan.ECexp.isin(N2JOS15.ECexp.unique())]
    lbl = CdlJOS15[["Gas", "SampleCode"]].iloc[0].to_list()
    fig, ax = plt.subplots(figsize=(8, 4))
    CdlJOS15.plot(
        x="E_RHE",
        y="Cdl_cat",
        ylim=(0, 0.03),
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        title=", ".join(lbl),
        ax=ax,
    )
    plt.savefig(ELoadDD.joinpath("CdlJOS15_cat.png"), dpi=300, bbox_inches="tight")

    ORRJOS15 = ORR_pars.loc[ORR_pars.ECexp.isin(N2JOS15.ECexp.unique())].query(
        "RPM > 1000"
    )
    lbl = ORRJOS15[["Gas", "SampleCode"]].iloc[0].to_list()

    lst = []
    for n, r in ORRJOS15.iterrows():
        ls = r.Loading_cm2
        rrde = pd.read_excel(r.RRDE_DataFile)
        rrde["Loading_cm2"] = ls
        lst.append(rrde)
    ORRJOS15rrde = pd.concat(lst)

    fig, ax = plt.subplots(figsize=(8, 4))
    ORRJOS15.plot(
        x="Loading_cm2",
        y="Jkin_075",
        ylim=(0.1, 5),
        logy=True,
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        title=", ".join(lbl),
        ax=ax,
        s=100,
        marker="^",
    )
    plt.savefig(ELoadDD.joinpath("ORRJOS15_Jkin_L.png"), dpi=300, bbox_inches="tight")

    orreisJ15 = ORReis_merge_raw.loc[
        ORReis_merge_raw.ECexp.isin(N2JOS15.ECexp.unique())
    ].query(f'RPM > 1000 & Model_EEC == "{fast_checking_EEC_models[2]}"')
    fig, ax = plt.subplots(figsize=(8, 4))
    orreisJ15.plot(
        x="E_onset",
        y="Rs",
        ylim=(0.05, 7),
        logy=True,
        kind="scatter",
        c="Loading_cm2",
        colormap=cmapset,
        title=", ".join(lbl),
        ax=ax,
        s=100,
        marker="^",
    )
    plt.savefig(ELoadDD.joinpath("ORREISJOS15_Rs.png"), dpi=300, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(8, 4))
    #    for load,lgr in ORRJOS15rrde.groupby('Loading_cm2'):
    disk = ax.scatter(
        x=ORRJOS15rrde["E_AppV_RHE_disk"],
        y=ORRJOS15rrde["Jcorr"],
        c=ORRJOS15rrde["Loading_cm2"],
        cmap=cmapset,
    )
    ax.set_ylim(-6.5, 0.1)
    ax.set_xlabel("E / V_RHE")
    ax.set_ylabel("j / mAcm-2")
    cbaxes = fig.add_axes([1.05, 0.1, 0.03, 0.8])
    plt.colorbar(disk, cax=cbaxes)
    #    disk = ORRJOS15rrde.plot(x='E_AppV_RHE_disk',y='Jcorr',ylim=(-6.5,0.1),logy=0, kind='scatter',c='Loading_cm2',colormap=cmapset, title=', '.join(lbl),ax=ax)
    axr = ax.twinx()
    ring = axr.scatter(
        x=ORRJOS15rrde["E_AppV_RHE_disk"],
        y=ORRJOS15rrde["Frac_H2O2"],
        c=ORRJOS15rrde["Loading_cm2"],
        cmap=cmapset,
        alpha=0.5,
        marker="x",
        s=10,
    )
    axr.set_ylim(0, 40)
    axr.set_ylabel("% / H2O2")
    #    ORRJOS15rrde.plot(x='E_AppV_RHE_disk',y='Frac_H2O2',xlim=(0,0.9),ylim=(0,30),logy=0, kind='scatter',c='Loading_cm2',colormap=cmapset, title=', '.join(lbl),ax=axr,s=10)
    plt.savefig(ELoadDD.joinpath("ORRJOS15rrde.png"), dpi=300, bbox_inches="tight")

    for n, r in (
        parsloadcorr.query('Gas == "O2" & E_RHE < 0.72 &  E_RHE > 0.67')
        .head(20)
        .iterrows()
    ):
        #        n,r
        Elec = tuple(r.to_list()[0:5])
        N2_elec = Elec[0:3] + ("N2", Elec[-1])
        Ekingrpget = Ekingrp.get_group(Elec)
        fig, ax = plt.subplots()
        Ekingrpget.plot(
            x="Loading_cm2",
            y=r.eispar,
            kind="scatter",
            ax=ax,
            c="r",
            s=100,
            label=Elec[3],
        )
        try:
            Ekingrp.get_group(N2_elec).plot(
                x="Loading_cm2",
                y=r.eispar,
                kind="scatter",
                ax=ax,
                c="blue",
                s=100,
                label=N2_elec[3],
            )
        except:
            continue
        ax.set_xlim = (0.6, 0.9)
        ax.set_title(", ".join([str(i) for i in Elec[0:3] + (Elec[-1],)]))
        ps = plotting.eisplot(r.eispar, yvalues=Ekingrpget[r.eispar].values)
        ax.set_ylim(ps.ylim)
        ax.set_yscale(ps.logyscale)
        #        ax.legend(True)
        plt.show()
        plt.close()


class Correlations:
    def __init__(self):
        pass

    @staticmethod
    def plot_triangle_hm(
        grE,
        rcorr,
        Ev,
        target_dir,
        corr_method_set="pearson",
        corr_cutoff=0.5,
        plot_option=False,
    ):
        #            rcorr = dfcorr[corr_Cols].corr(method=corr_method)
        #        Ev = np.round(EvRHE,2)
        #        grE,rcorr,,target_dir,corr_method_set=corr_method,plot_option=False
        heatmap_title = "EIS_heatmap_{0}_{1}_{2}.png".format(
            corr_method_set, int(corr_cutoff * 100), Ev
        )
        print(heatmap_title)
        #        selcorr = rcorr[((rcorr > corr_cutoff) | (rcorr < corr_cutoff)) & (rcorr.abs() < 1)]
        corr_triu = rcorr.where(np.tril(np.ones(rcorr.shape)).astype(np.bool))
        scorr_triu = corr_triu.stack()
        #        filter_corr = scorr_triu[(scorr_triu.abs() > corr_cutoff) & (scorr_triu.abs() < 1)]
        #        fcorr = filter_corr.unstack()
        #        fcorr_sliced = fcorr.dropna(how='all')
        #        drop_cols = ['Colorcode','lmfit_redchi2','lmfit_redchi1']
        #        drop_cols_set = [i for i in drop_cols if i in fcorr_sliced.columns]
        #        unst_rcorr = rcorr.dropna(how='all')
        #        promising = fcorr_sliced.loc[fcorr_sliced.index.isin(SampleSelection.EC_EIS_par_cols[0:-2])].unstack().dropna().sort_values()
        promising = (
            rcorr.loc[rcorr.index.isin(SampleSelection.EC_EIS_par_cols[0:-2])]
            .unstack()
            .dropna()
            .sort_values()
        )
        promising = (
            promising.drop(
                index=[
                    ("Rct", "Rct_kin"),
                    ("Rct_kin", "Rct"),
                    ("Cdlp", "Qad+Cdlp"),
                    ("Qad+Cdlp", "Cdlp"),
                    ("Qad", "Qad+Cdlp"),
                    ("Qad+Cdlp", "Qad"),
                    ("Qad+Cdlp", "Qad+Cdlp"),
                ]
                + [i for i in promising.index if i[0] == i[1]],
                errors="ignore",
            )
            .drop_duplicates()
            .sort_values()
        )

        if not promising.unstack().empty and plot_option is "heatmap":
            plt.subplots(figsize=(40, 40))
            sns.heatmap(promising.unstack())
            #                heatmap_title = 'EIS_heatmap_{0}_{1}_{2}'.format(corr_method,int(corr_cutoff*100),EvRHE)
            plt.suptitle(heatmap_title)
            plt.savefig(
                target_dir.joinpath(heatmap_title), dpi=300, bbox_inches="tight"
            )
            plt.grid(True)
            plt.close()
        else:
            #            fcorr_sliced.loc[fcorr_sliced.index.isin(SampleSelection.EC_EIS_par_cols)]
            #                drop_cols_set = [i for i in drop_cols if i in fcorr_sliced.columns]
            #                if not promising.loc[[('Rct','Rct_kin'),('Rct_kin','Rct')]].dropna().empty:
            #            promising = promising.drop(index=[('Rct','Rct_kin'),('Rct_kin','Rct'),('Cdlp','Qad+Cdlp'), ('Cdlp','Qad+Cdlp'),('Qad','Qad+Cdlp'), ('Qad+Cdlp','Qad')],errors='ignore')
            toptail = pd.concat([promising.tail(20), promising.head(20)])
            PPDEIS_Ev = target_dir.joinpath(str(Ev))
            PPDEIS_Ev.mkdir(parents=True, exist_ok=True)
            if plot_option == "corr":
                for (xc, yc), corr_val in toptail.iteritems():
                    plot_pd_SampleIDs(grE, xc, yc, corr_val, PPDEIS_Ev)
        return promising

    # ===== Making plots of correlations of EIS Parameters versus E v RHE ======
    def MakeTopCorrs(EIS_CB_paper, corr_method, PPDEIS):
        Gases, pHses, Erheses = ["N2", "O2"], [0.3, 1, 13], [0, 0.4, 0.6, 0.8, 1]
        #        EIS_O2_no = EIS_CB_paper.query('(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))').drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        #                pH03, pH1 = EIS_O2_acid_no.query('(pH == 0.3)'), EIS_O2_acid_no.query('(pH == 1)')

        #        corr_method,corr_cutoff = 'spearman',0.1 # or default is pearson spearman
        #        rcorr = EIS_O2_065_acid_no[corr_Cols].corr(method=corr_method)
        corr_Cols_unfiltered = (
            SampleSelection.InterestingCols
            + SampleSelection.EC_EIS_par_cols
            + SampleSelection.RAMAN_cols_corr
            + SampleSelection.Inter_MS_cols
        )
        drop_corr_cols = [
            "Colorcode",
            "lmfit_redchi2",
            "lmfit_redchi1",
            "Area_in_cell",
            "BET_Area_RPT",
        ]
        corr_Cols = [i for i in corr_Cols_unfiltered if i not in drop_corr_cols]
        corr_Cols_filtered = [i for i in corr_Cols if i not in drop_corr_cols]

        out_topcorrs_lst = []
        for Gas_set in Gases:
            for pH_set in pHses:
                EIS_O2_no_query = EIS_CB_paper.query(
                    '(RPM_DAC > 700) & (Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
                )
                EIS_O2_no_query = EIS_O2_no_query.drop_duplicates(
                    subset=[
                        i
                        for i in SampleSelection.EC_EIS_par_cols
                        if i in EIS_O2_no_query
                    ]
                )
                #                ORREIS_O2_no_query = EIS_CB_paper.query('(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))').drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
                target_dir = PPDEIS.joinpath(
                    "EIS_corr_{0}_pH{1}".format(Gas_set, pH_set)
                )
                target_dir.mkdir(parents=True, exist_ok=True)
                for EvRHE, grE in EIS_O2_no_query.groupby("E_RHE"):
                    if len(grE) > 3:
                        EvRHE, grE
                        rcorr = (
                            grE[corr_Cols_filtered]
                            .corr(method=corr_method)
                            .drop_duplicates()
                        )
                        prom_corr = Correlations.plot_triangle_hm(
                            grE,
                            rcorr,
                            np.round(EvRHE, 2),
                            target_dir,
                            corr_method_set=corr_method,
                            plot_option=False,
                        )
                        out_topcorrs_lst.append(
                            [
                                Gas_set,
                                pH_set,
                                EvRHE,
                                len(grE),
                                grE.SampleID.nunique(),
                                prom_corr,
                            ]
                        )

        #        pd.DataFrame(out_topcorrs_lst[0][-1],columns=[corr_method]).assign(**{'Gas' : out_topcorrs_lst[0][0], 'pH' : out_topcorrs_lst[0][1], 'E_RHE' : out_topcorrs_lst[0][2]})
        topcorrs = pd.concat(
            [
                pd.DataFrame(i[-1], columns=[corr_method]).assign(
                    **{
                        "Gas": i[0],
                        "pH": i[1],
                        "E_RHE": i[2],
                        "lenGr": i[3],
                        "nSamples": i[4],
                    }
                )
                for i in out_topcorrs_lst
            ]
        )
        ECpar1, ECpar2 = [
            1 if i[0] in SampleSelection.EC_EIS_par_cols else 0 for i in topcorrs.index
        ], [1 if i[1] in SampleSelection.EC_EIS_par_cols else 0 for i in topcorrs.index]
        ECpars = [i[0] + i[1] for i in zip(ECpar1, ECpar2)]
        topcorrs = topcorrs.assign(
            **{
                "score": np.abs(topcorrs[corr_method]) * topcorrs.lenGr,
                "SampleScore": np.abs(topcorrs[corr_method]) * topcorrs.nSamples,
                "ECpars_duplicates": ECpars,
            }
        )

        topcorrs = topcorrs.sort_values("SampleScore", ascending=0)
        topcorr_best = topcorrs[
            (np.abs(topcorrs[corr_method]) > 0.3)
            & (topcorrs.nSamples > topcorrs.nSamples.describe()["25%"])
            & (topcorrs.lenGr >= topcorrs.lenGr.describe()["25%"])
        ].sort_values(by=corr_method)
        topcorr_best_Ev = topcorrs[
            topcorrs.E_RHE.isin(Erheses) & (topcorrs.ECpars_duplicates < 2)
        ].sort_values(by="SampleScore", ascending=False)
        return topcorrs, topcorr_best, topcorr_best_Ev


class ExportECfromEIS:
    """Made correlations by guessing for best,
    but should probably choose several E_RHE values eg. [0.2,0.4,0.6,0.8,1] to narrow it down."""

    EvRHE = "E_AppV_RHE"
    #    OriginColor = FileHelper.FindExpFolder().LoadOriginColor()
    """ Refit: pH1 : DW16,DW19,'DW17','DW28'
    """

    def __init__(self):
        pass

    def change_col():
        PostEISDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
            "EC_files_out"
        )
        EISl = list(PostEISDir.rglob("*EIS/*/EIS_*"))
        for fl in EISl:
            a = pd.read_excel(fl, index_col=0)
            a["Qad+Cdlp"] = a["Qad"] + a["Cdlp"]
            a.to_excel(fl)

    #                rcorr = gr_sID_desc_Codes[SampleSelection.InterestingCols+SampleSelection.RAMAN_cols_corr].corr(method=corr_method)

    def EIS_plotting_parameters():
        #%% ===== MAKE PLOTS of Parameters versus E v RHE ======
        # OnlyRecentMissingOVV = run_PAR_DW.ECRunOVV(load=1).index
        EIS_pars = Load_from_Indexes.EIS_pars_OVV(reload=False)  # EIS_Pars2
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")

        SeriesIDs = [
            SampleSelection.Series_CB_paper,
            SampleSelection.Series_Porhp_SiO2,
            SampleSelection.Series_Co_PANI,
            SampleSelection.Series_ML_SiO2,
            {"name": "all_EIS"},
        ]
        SeriesID_set = SeriesIDs[0]
        PostDestDirEIScom = PostDestDir.joinpath(
            "EIS_Pars_Char_{0}".format(SeriesID_set["name"])
        )
        PostDestDirEIScom.mkdir(parents=True, exist_ok=True)

        PPDEIS = PostDestDirEIScom.joinpath("EIS_slices")
        PPDEIS.mkdir(parents=True, exist_ok=True)

        EIS_CB_paper_raw = EIS_pars.loc[
            EIS_pars.SampleID.isin(SeriesID_set["sIDslice"])
        ]
        if "all" in SeriesID_set["name"]:
            EIS_CB_paper_raw = EIS_pars

        #        OVV_CB =  OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.PAR_exp.isin(EIS_CB_paper_raw.Type_exp.unique()) & OnlyRecentMissingOVV.SampleID.isin(EIS_CB_paper_raw.SampleID.unique())]
        #        eismiss_CB = OVV_CB.loc[OVV_CB.PAR_file.isin([i for i in OVV_CB.PAR_file.values if i not in EIS_CB_paper_raw.PAR_file.values])]
        #        EIS_CB_paper_raw,HPRR_CB_paper,ORR_CB_paper, Cdl_CB_paper, HER_CB_paper, OER_CB_paper = Load_from_Indexes.IndexPars_CB_paper()
        #        postOVVout  = PostEC.LoadPostOVV()

        def StandardEIS(EIS_pars, SeriesID_set):
            standard_EIS_pars_fn = (
                FileHelper.FindExpFolder("VERSASTAT")
                .DestDir.joinpath("PostEC")
                .joinpath("EIS_standard_pars_{0}.xlsx".format(SeriesID_set))
            )
            if standard_EIS_pars_fn.is_file():
                EIS_Standard_Pars = pd.read_excel(standard_EIS_pars_fn)
            else:
                EIS_std = []
                for nEC, ECgr in EIS_pars.groupby(
                    SampleSelection.EC_exp_cols[1:] + ["E_RHE"]
                ):
                    cond = ECgr[
                        [
                            i
                            for i in SampleSelection.EC_EIS_par_cols
                            if not "lmfit_redchi" in i
                        ]
                    ].query("(((Rs > 1) & (Rs < 150)) & (Rorr < 9E04)  & (Rct < 9E04))")
                    mean = cond.describe().loc["mean"]
                    means_out = dict(
                        zip(
                            list(mean.index)
                            + SampleSelection.EC_exp_cols[1:]
                            + ["E_RHE"],
                            list(mean.values) + list(nEC),
                        )
                    )
                    EIS_std.append(means_out)
                EIS_Standard_Pars = pd.DataFrame(EIS_std)
                EIS_Standard_Pars.to_excel()

        #    pd.read_excel(PostDestDir.joinpath('postEC_Organized.xlsx'),index_col=[0])
        #    EIS_parfiles = ['_'.join((Path(i).stem).split('_')[0:-1]) for i in postEIScom.SourceFilename.values]
        #        EIS_parfiles = [((Path(i).stem).split('_pars2')[0]) for i in postEIScom.SourceFilename.values]
        #        SampleCodes  = pd.read_excel(FileHelper.FindExpFolder('VERSASTAT').DestDir.parent.joinpath('SampleCodeLst.xlsx'))
        #    pd.read_excel(PostDestDir.joinpath('SampleCodeLst.xlsx'))
        #        EIS_pars = pd.concat([pd.read_excel(i,index_col=[0]) for i in EIS_Pars_files]).rename(columns={'File' : 'PAR_file'})
        OriginColor = FileHelper.FindExpFolder.LoadOriginColor()
        EIS_CB_paper = EIS_CB_paper_raw.query(
            "(((Rs > 4) & (Rs < 150))| (Rorr < 9E04))"
        ).dropna(subset=["SampleLabel"])
        EIS_CB_paper_wrongfits = EIS_CB_paper_raw.query(
            "~(((Rs > 4) & (Rs < 150))| (Rorr < 9E04))"
        )
        EIS_CB_paper = ExportECfromCV.make_uniform_EvRHE(EIS_CB_paper)
        EIS_CB_paper.Loading_cm2 = [
            np.round(i, 3) for i in EIS_CB_paper.Loading_cm2.values
        ]
        EIS_pars_set = EIS_CB_paper
        CBn, CBs = EIS_CB_paper.SampleID.nunique(), set(EIS_CB_paper.SampleID.unique())
        SampleSelection.EC_exp_cols
        EIS_pars_set.plot(y="nAd", x="Loading_cm2", kind="scatter")
        EIS_O2 = EIS_pars_set.query('Gas == "O2"')
        EIS_O2_065_acid_no = EIS_pars_set.query(
            '(Gas == "O2") & (E_RHE == 0.65) & (pH < 6) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        EIS_O2_065_acid1_no = EIS_O2_065_acid_no.query("(pH == 1)")
        EIS_O2_065_acid03_no = EIS_O2_065_acid_no.query("(pH == 0.3)")
        EIS_O2_045_acid_no = EIS_pars_set.query(
            '(Gas == "O2") & (E_RHE == 0.45) & (pH < 6) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        EIS_O2_acid_no = EIS_pars_set.query(
            '(Gas == "O2") & (pH < 6) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        EIS_O2_alk_no = EIS_pars_set.query(
            '(Gas == "O2") & (pH > 6) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        EIS_O2_no_normalload = EIS_CB_paper.query(
            '(postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        #        EIS_O2_065_acid_no.drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)

        corr_Cols_unfiltered = (
            SampleSelection.InterestingCols
            + SampleSelection.EC_EIS_par_cols
            + SampleSelection.RAMAN_cols_corr
            + SampleSelection.Inter_MS_cols
        )
        drop_corr_cols = ["Colorcode", "lmfit_redchi", "Area_in_cell", "BET_Area_RPT"]
        corr_Cols = [i for i in corr_Cols_unfiltered if i not in drop_corr_cols]
        corr_Cols_filtered = [i for i in corr_Cols if i not in drop_corr_cols]

        corr_method, corr_cutoff = "spearman", 0.1  # or default is pearson spearman
        rcorr = EIS_O2_065_acid_no[corr_Cols].corr()
        plt.subplots(figsize=(20, 15))
        selcorr = rcorr[(rcorr > corr_cutoff) | (rcorr < corr_cutoff)]
        sns.heatmap(selcorr)
        plt.savefig(
            PostDestDirEIScom.joinpath(
                "EIS_heatmap_{0}_{1}.png".format(corr_method, corr_cutoff)
            ),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
        #        melt= pd.melt(rcorr.reset_index(), id_vars='index')
        #        melt.columns
        ##            FileHelper.PlotAnalysis.corrplot(rcorr)
        #
        topcorrs, topcorr_best, topcorr_best_Ev = MakeTopCorrs(
            EIS_CB_paper, corr_method
        )
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
        #            spdf.update(spfrow)

        EIScorrsPPD = PPDEIS.joinpath("top_corrs_{0}_uniform".format(corr_method))
        EIScorrsPPD.mkdir(parents=True, exist_ok=True)
        EIS_loop_topcorrs(EIS_O2_no_normalload, topcorr_best_Ev, EIScorrsPPD)

        def EIS_loop_topcorrs(EIS_O2_no_normalload, topcorr_best_Ev, EIScorrsPPD):
            for pH_set in pHses:
                for Gas_set in Gases:
                    EIS_pH_gas = EIS_O2_no_normalload.query(
                        "(Gas == @Gas_set) & (pH == @pH_set) & (RPM > 100) & (Rct < 8000) & (Rct_kin < 20) & (Rs > 3)"
                    )
                    if not EIS_pH_gas.empty:
                        EISPPD_uniform_Ev = PostDestDirEIScom.joinpath(
                            "per_Erhe/pH{0}/{1}".format(pH_set, Gas_set)
                        )
                        EISPPD_uniform_Ev.mkdir(parents=True, exist_ok=True)
                        for ypar in eisplot.parlst:
                            plot_pd_SampleIDs(
                                EIS_pH_gas, "E_RHE", ypar, 0, EISPPD_uniform_Ev
                            )
                    for Erhe_set in Erheses:
                        #                        Gas_set,pH_set,Erhe_set,corr_val = tr[1].Gas,tr[1].pH,tr[1].E_RHE,tr[1][corr_method]
                        grE = EIS_O2_no_normalload.query(
                            "(Gas == @Gas_set) & (pH == @pH_set) & (E_RHE == @Erhe_set) & (RPM > 100)"
                        )
                        #                        grE = grE.query('(Rct < 250000)')

                        #                       Rct_slice = topcorr_best_Ev.loc[topcorr_best_Ev.index.isin([i for i in topcorr_best_Ev.index if 'Rct' in i])]
                        tpcb = topcorr_best_Ev.query(
                            "(Gas == @Gas_set) & (pH == @pH_set) & (E_RHE == @Erhe_set)"
                        )
                        #                        topcorr_best_Ev.index.isin(Rct_slice)
                        if not grE.empty and not tpcb.empty:
                            EISPPD_uniform = EIScorrsPPD.joinpath(
                                ("pH{0}/{1}").format(pH_set, Gas_set)
                            )
                            EISPPD_uniform.mkdir(parents=True, exist_ok=True)
                            xc, yc, corr_val = (
                                tpcb.index[0][0],
                                tpcb.index[0][1],
                                tpcb.iloc[0][corr_method],
                            )
                            plot_pd_SampleIDs(grE, xc, yc, corr_val, EISPPD_uniform)
                            #            spectras = ExportECfromCV.make_uniform_EvRHE(pd.concat([pd.read_excel(i) for i in grE.SpectraFile.values],sort=False)).query('(E_RHE == @Erhe_set)')
                            spectra_lst = []
                            for n, spfrow in grE.iterrows():
                                spf = spfrow.SpectraFile
                                spdf = ExportECfromCV.make_uniform_EvRHE(
                                    pd.read_excel(spf)
                                )
                                spdf_Ev = spdf.query("(E_RHE == @Erhe_set)")
                                spdf_char = [
                                    spdf_Ev.assign(
                                        **{
                                            i: [spfrow[i]] * len(spdf_Ev.index)
                                            for i in spfrow.index
                                        }
                                    )
                                ][0]
                                spectra_lst.append(spdf_char)
                            spectras = pd.concat([i for i in spectra_lst], sort=False)
                            plot_spectra_SampleIDs(xc, yc, spectras, EISPPD_uniform)

        def EIS_PostCorr_Refit():
            cond1 = {
                "pH": 0.3,
                "Gas": "O2",
                "Loading_cm2": 0.379,
                "postAST": "no",
                "E_RHE": [0.35, 0.45],
                "Samples": ["JOS12", "JOS13"],
            }
            cond3 = {
                "pH": 0.3,
                "Gas": "N2",
                "Loading_cm2": 0.379,
                "postAST": "no",
                "E_RHE": [0.25],
                "Samples": ["JOS15"],
            }
            cond4 = {
                "pH": 13,
                "Gas": "N2",
                "Loading_cm2": 0.379,
                "postAST": "no",
                "E_RHE": [0.75],
                "Samples": ["JOS12", "DW21"],
            }

            cond2 = {
                "pH": 1,
                "Gas": "O2",
                "Loading_cm2": 0.379,
                "postAST": "no",
                "E_RHE": [0.6],
                "Samples": ["DW19", "DW25", "DW24"],
            }

        #        def plot_Spectra_SampleIDs():
        topcorr_index = (
            topcorrs[(np.abs(topcorrs[corr_method]) > 0.7)]
            .reset_index()
            .set_index(["level_0", "level_1", "Gas", "pH", "E_RHE"])
            .groupby(level=[0, 1, 2, 3, 4])
            .sum()
        )
        Edepence_corr_lst = []
        for lbl, gr in topcorr_index.reset_index(level=[4]).groupby(level=[0, 1, 2, 3]):
            if len(gr) > 3:
                Ercorrcorr = (
                    gr[["E_RHE", corr_method]].corr(method=corr_method).iloc[0, 1]
                )
                sts = gr[corr_method].describe().to_list()
                #            if np.abs(Ercorrcorr) > 0.95:
                if sts[0] > 4:
                    gr.plot(x="E_RHE", y=corr_method, title=str(lbl), kind="scatter")
                Edepence_corr_lst.append(list(lbl) + [Ercorrcorr] + sts)
        corr_E_dependence = pd.DataFrame(
            Edepence_corr_lst,
            columns=[0, 1, "Gas", "pH", corr_method + "_EvRHE"]
            + ["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
        ).sort_values(by=corr_method + "_EvRHE")

        topcorrs.groupby(level=[0, 1]).size()
        topcorrs.groupby(level=[0, 1]).sum()
        topcorrs.loc[np.abs(topcorrs[corr_method] > 0.5)].groupby(
            level=[0, 1]
        ).sum().sort_values(by=corr_method).iloc[-10::].plot.barh(x=corr_method)
        topcorrs[corr_method].groupby(level=[0, 1]).sum().sort_values().iloc[
            -10::
        ].plot.barh(x=corr_method)
        top_score_sum = topcorrs[corr_method].groupby(level=[0, 1]).sum().sort_values()
        top_score_sum.plot.barh()
        top_score_sum[(12 < top_score_sum) | (-11 > top_score_sum)].plot.barh(
            figsize=(10, 16)
        )
        top_score_abssum = (
            np.abs(topcorrs[corr_method]).groupby(level=[0, 1]).sum().sort_values()
        )
        topabssum_best = top_score_abssum.loc[
            top_score_abssum > top_score_abssum.mean() + 1 * top_score_abssum.std()
        ]
        topabssum_best.plot.barh(figsize=(12, 12))

        topcorrs.query('(pH < 5) & (Gas == "O2")').index
        sumbest_acid_E = (
            topcorrs.query('(pH < 5) & (Gas == "O2")')[[corr_method, "E_RHE"]]
            .reset_index()
            .set_index(["level_0", "level_1", "E_RHE"])
            .groupby(level=[0, 1, 2])
            .sum()
        )
        sumbest_acid = (
            sumbest_acid_E[np.abs(sumbest_acid_E) > 0.7].dropna().reset_index("E_RHE")
        )
        sumbest_acid.plot(x="E_RHE", y=corr_method, kind="scatter")
        sumbest_acid.groupby(level=[0, 1])

        for nb, bgr in EIS_O2_no_normalload.groupby(by=["pH", "Gas", "E_RHE"]):
            best_cols = [
                i
                for i in bgr.columns
                if any([i in a for a in list(topabssum_best.index)])
            ]
            if len(bgr) > 2:
                best_target_dir = PPDEIS.joinpath(
                    "best_EIS_corr_{0}_pH{1}".format(nb[1], nb[0])
                )
                best_target_dir.mkdir(parents=True, exist_ok=True)
                best_rcorr = bgr[best_cols].corr(method=corr_method)
                best_corr2 = plot_triangle_hm(
                    bgr, best_rcorr, nb[2], best_target_dir, plot_option="corr"
                )

        #        for corrdx in topabssum_best.index:
        """ Alkaline
            Qad vs N_content /
            Qad vs H_content /
            Qad vs BET_Area_RPT \

            Cdlp vs BET_cat / low E

            nAd vs Yield_prec
            nAd vs ML \

            nDL vs N_content \

            Rs vs Metal_wt \
            Rs vs C/N_ratio /

            Rorr vs H_content /

            Rct vs BET_Area_RPT \

        """
        ylblE = ["C/N_ratio", "BET_cat_agg"]
        for E, gr in EIS_O2_065_acid_no.groupby("E_RHE"):
            titles = []
            if gr.Loading_cm2.nunique() < gr.SampleID.nunique():
                fig, ax = plt.subplots()
                for sID, Sgr in gr.groupby("SampleID"):
                    #                    if len(Sgr[ylblE].unique()) > 1:
                    titles.append(sID)
                    oc = (
                        OriginColor.loc[
                            OriginColor.OriginCode == Sgr.Colorcode.unique()[0], "color"
                        ]
                        .values[0]
                        .lower()
                    )
                    Sgr.plot(
                        y="Rct",
                        x=ylblE[0],
                        kind="scatter",
                        c=oc,
                        label="{} at {:.2f}".format(sID, E),
                        ylim=(0, 2e3),
                        ax=ax,
                    )
                plt.show()
                plt.close()
        PostDestDirEIScom = PostDestDir.joinpath("EIS_Pars_Char_CB_paper")
        CharLst = [
            i
            for i in EIS_pars.columns
            if i in SampleSelection.InterestingXcorrCols and i not in eisplot.parlst
        ]
        CharLst = [
            "BET_cat_agg",
            "N_content",
            "AD/AG",
            "SizeNP",
            "Metal_wt",
            "ML",
            "BET_cat_micro",
        ]
        for CY in CharLst:
            for yPar in eisplot.parlst:
                PDDirEIScomChar = PostDestDirEIScom.joinpath("{0}/{1}".format(CY, yPar))
                for Elec, ElecGr in EIS_pars_set.dropna(subset=[CY]).groupby(
                    by="Electrolyte"
                ):
                    PDDirEIScom = PDDirEIScomChar.joinpath(Elec)
                    PDDirEIScom.mkdir(exist_ok=True, parents=True)
                    #        ElecEIS = pd.concat([pd.read_excel(i) for i in ElecGr.SourceFilename.values])
                    sliceElecGr = ElecGr.loc[
                        ElecGr.E_RHE.isin(ElecGr.E_RHE.unique()[0::2])
                    ]
                    for St, stGr in sliceElecGr.groupby(by="postAST"):
                        for Ev, Egr in stGr.groupby(by="E_RHE"):
                            fig, ax = plt.subplots(1, 1)
                            ax.set_xlim(eisplot(yPar).lim)
                            for sID, sGr in Egr.groupby(by="SampleID"):
                                for Gas, gGr in sGr.groupby(by="Gas"):
                                    if not gGr.empty:
                                        ms = eisplot("Gas {0}".format(Gas)).ms
                                        ax.scatter(
                                            gGr[yPar],
                                            gGr[CY],
                                            label=str(
                                                "{0}_{1}_{2}".format(sID, Gas, St)
                                            ),
                                            s=40,
                                            alpha=ms[1],
                                            marker=ms[0],
                                        )
                            ax.legend(
                                bbox_to_anchor=(0.5, 1.8),
                                ncol=4,
                                loc="upper center",
                                fontsize=12,
                            )
                            ax.set_ylabel(CY)
                            ax.set_xlabel(yPar)
                            ax.set_title(
                                "{0} at {1} AST {2}".format(Elec, np.round(Ev, 2), St)
                            )
                            ax.grid(True)
                            DestFile = PDDirEIScom.joinpath(
                                "{0}_{1}_{2}.png".format(Elec, np.round(Ev, 2), St)
                            )
                            plt.savefig(DestFile, dpi=300, bbox_inches="tight")
                            plt.close()
        #                                    sGr.plot(x=yPar,y=CY,ax=ax,kind='scatter',label=sID)
        #                        plt.close()
        for Elec, ElecGr in EIS_pars.groupby(by="Electrolyte"):
            PDDirEIScom = PostDestDirEIScom.joinpath(Elec)
            for Sample, sGr in ElecGr.groupby(by="SampleID"):
                PDDirEIScom.joinpath("Pars_per_Sample").mkdir(
                    exist_ok=True, parents=True
                )
                PDDirEIScomSamples = PDDirEIScom.joinpath("Pars_per_Sample")
                for St, stGr in sGr.groupby(by="postAST"):
                    if not stGr.empty:
                        try:
                            ExportECfromEIS.EIS_ParsPlotting_Rs_per_Sample(
                                stGr, PDDirEIScomSamples
                            )
                        except Exception as e:
                            print(Elec, Sample, St, "\n{0}".format(e))

                    for sf, fGr in stGr.groupby(by="SourceFilename"):
                        Cdl = pd.DataFrame()
                        #                    EIS_parfile = '_'.join((Path(fGr.SourceFilename.unique()[0]).stem).split('_')[0:-1])
                        #                    ['_'.join((Path(i).stem).split('_')[0:-1]) for i in postEIScom.SourceFilename.values]
                        AllData_E_file = pd.read_excel(sf)
                        data_parts, sf_parts = FileHelper.FileOperations.find_CS_parts(
                            AllData_E_file.File.iloc[0]
                        ), FileHelper.FileOperations.find_CS_parts(sf)
                        AllData_E_file["EIS_parfile"] = [
                            (Path(i).stem) for i in AllData_E_file.File.values
                        ]
                        if data_parts[0] != sf_parts[0]:
                            update_file = [
                                sf_parts[0].joinpath(i[1])
                                for i in [
                                    FileHelper.FileOperations.find_CS_parts(i)
                                    for i in AllData_E_file.File.values
                                ]
                            ]
                            AllData_E_file.File = update_file
                        #                    fGr.EIS_parfile.unique()[0]
                        AllData_E_file = AllData_E_file.query(
                            "EIS_parfile == @fGr.EIS_parfile.unique()[0]"
                        )  # Need to select first right data from _pars2 Excel file!
                        #                    AllData_E_file.loc[(AllData_E_file['SampleID'] == Sample) & (AllData_E_file.File.str.contains(EIS_parfile) == True),:]
                        try:
                            AllData_E_file = pd.merge(
                                AllData_E_file, SampleCodes, on="SampleID"
                            )
                        except:
                            print(
                                "Merging Data and SampleCodes unsuccessfull, assigning NoName"
                            )
                            AllData_E_file["Sample"] = "NoName"
                        AllData_E_file["Status"] = fGr.Status.values[0]
                        Cdl_ovv = postOVVout.loc[
                            (postOVVout["Type_Exp"] == "N2_Cdl")
                            & (postOVVout["Electrolyte"] == Elec)
                            & (postOVVout["SampleID"] == Sample)
                            & (postOVVout.Status == St)
                        ].drop_duplicates()
                        Cdl = (
                            pd.concat([pd.read_csv(i) for i in Cdl_ovv.SourceFilename])
                            .drop_duplicates()
                            .query('(Sweep_Type_N2 == "anodic") & (Cdl_R > 0.85)')
                        )
                        if Cdl.empty:
                            print("Cdl is empty, skipped")
                        #                Cdl.plot(x=EvRHE,y='Cdl',kind='scatter')
                        DestFile = PDDirEIScom.joinpath(
                            "_".join(
                                [
                                    fGr.SampleID.values[0],
                                    fGr.Gas.values[0],
                                    fGr.Status.values[0],
                                    fGr.EXP_date.values[0],
                                ]
                            )
                        ).with_suffix(".png")
                        DestFile.parent.mkdir(exist_ok=True, parents=True)
                        print(DestFile)
                        #                AllData_E_file.to_excel(DestFile.with_suffix('.xlsx'))
                        #                    SampleCode = SampleCodes.loc[SampleCodes.SampleID == fGr['SampleID'].unique()[0],:]
                        #                    if SampleCode.empty:
                        #                        print('SampleCode empty')
                        #                        SampleCode = pd.DataFrame({'Sample' : 'NoName'},index=[0])
                        #                    SampleLabel = SampleCode.Sample.values[0]
                        if not AllData_E_file.empty:
                            EIS_ParsPlotting(AllData_E_file, Cdl, DestFile)
                            SampleData = pd.concat([SampleData, AllData_E_file])
                #            DestFile = PDDirEIScom.joinpath('_'.join([Sample).with_suffix('.png')
                #            SampleData.plot(x=EvRHE,y='Rct_kin',kind='scatter')
                if not SampleData.empty:
                    EIS_ParsPlotting_Rs_per_Sample(SampleData, PDDirEIScomSamples)

    def CB_RPM_ext():
        """Goal is to choose a proper model or depending on effect of rpm on parameters to determine wich ones are related
        to diffusion."""
        PostDestDir = FindExpFolder("VERSASTAT").PostDir

        EIS_pars_for_RPM_series = [
            (n, gr.RPM_DAC.unique())
            for n, gr in EIS_pars.groupby("SampleID")
            if gr.RPM_DAC.nunique() > 2
        ]
        print(EIS_pars_for_RPM_series)
        # EIS_pars_for_Loading_series = [(n, gr.Loading_cm2.unique()) for n,gr in EIS_pars.groupby('SampleID') if gr.Loading_cm2.nunique() > 2]
        # print(EIS_pars_for_Loading_series)
        CBsamples = SampleSelection.Series_CB_paper["sIDs"]
        EIS_CB_paper = EIS_pars.loc[EIS_pars.SampleID.isin(CBsamples)]
        print("Unique RPMs:", EIS_CB_paper.RPM_DAC.unique())
        #        EIS_CB_paper.groupby('E_RHE')
        #        CB_rpm_series =
        #        [(n,gr.RPM_DAC.unique(),gr) for n,gr in EIS_CB_paper.groupby(['E_RHE','SampleID']) if len(gr.groupby('RPM_DAC')) >= 3]
        # TODO:   2020-01-09 stopped at plotting settings for series, try make correlations evt....
        CB_EIS_rpm = EIS_pars.loc[
            EIS_pars.SampleID.isin([i[0] for i in EIS_pars_for_RPM_series])
        ]
        CB_EIS_rpm.to_pickle(PostDestDir.joinpath("EIS_RPM_series.pkl.compress"))

        CB_EIS_rpm = CB_EIS_rpm.assign(
            **{"Rorr_kin": 1 / CB_EIS_rpm.Rorr, "RPM_sqrt": np.sqrt(CB_EIS_rpm.RPM_DAC)}
        )
        PDD_rpm = PostDestDir.joinpath(
            "EIS_{0}/RPM_series".format(SampleSelection.Series_CB_paper["name"])
        )
        PDD_rpm.mkdir(parents=True, exist_ok=True)
        # rpm_serie_slice = eisplot.parlst+['RPM_DAC','R3_kin','Rorr_kin']

        _models = EIS_selection.mod_select
        rpm_series_test_corr_lst = []
        rpm_ser_idx_grp = ["SampleID", "Gas", "pH", "E_RHE", "Model_EEC"]
        CB_EIS_rpm_fltr = CB_EIS_rpm.query("RPM_DAC < 3000 & pH < 15")
        CB_EIS_rpm_grp = CB_EIS_rpm_fltr.groupby(rpm_ser_idx_grp)
        CB_EIS_rpm_grp_mods = CB_EIS_rpm.query("RPM_DAC < 3000 & pH < 15").groupby(
            rpm_ser_idx_grp[0:-1]
        )

        _cm = plt.get_cmap("tab10")
        for EV, EVgr in CB_EIS_rpm_grp:
            EVidx = {}
            if EVgr.RPM_DAC.nunique() > 3:
                # === PLOT RPM ===
                rpmgrp = EVgr
                _rawdata = pd.concat(
                    pd.read_excel(i) for i in rpmgrp.File_SpecFit.values
                )
                _rawdata_mod = _rawdata.loc[_rawdata.Model_EEC == EV[-1]]

                for pf, pfgrp in _rawdata_mod.groupby(["PAR_file"]):

                    _parpf = rpmgrp.query("PAR_file == @pf")
                    if _parpf.tau.sum() > 0 and len(_parpf) >= 2:
                        print(EV, Path(pf).stem)
                        fig, ax = plt.subplots(2, 2, figsize=(12, 12))
                        # _parpf.plot(x=['RPM_DAC']*2, y =['Aw','Rct'], kind='scatter', ax=ax[0][1])
                        _ps = [("red", "o"), ("green", "*")]
                        [
                            _parpf.plot(
                                x="RPM_sqrt",
                                y=_ypar,
                                label=_ypar,
                                kind="scatter",
                                ax=ax[0][1],
                                c=_ps[n][0],
                                marker=_ps[n][1],
                            )
                            for n, _ypar in enumerate(["Aw", "Rct"])
                        ]
                        # np.array(plt.get_cmap('Dark2')(n))
                        _lin = linregress(
                            _parpf.query("RPM_DAC > 0").RPM_sqrt,
                            _parpf.query("RPM_DAC > 0").tau,
                        )
                        _parpf["tau_lin"] = (
                            _parpf.RPM_sqrt * _lin.slope + _lin.intercept
                        )
                        _parpf.plot(
                            x="RPM_sqrt",
                            y="tau",
                            label="tau",
                            kind="scatter",
                            ax=ax[1][1],
                        )
                        _parpf.plot(
                            x="RPM_sqrt",
                            y="tau_lin",
                            label=f"lin: {_lin.slope:.3E} + {_lin.intercept:.3E}, r={_lin.rvalue:.3f}",
                            kind="line",
                            ax=ax[1][1],
                        )

                        for n, (rpm, rgrp) in enumerate(pfgrp.groupby(["RPM_DAC"])):
                            _c = np.array(_cm(n))

                            rgrp.sort_values("Frequency(Hz)", inplace=True)
                            # len(rgrp)*
                            print(n, rpm, _c)
                            ax[0][0].scatter(
                                rgrp["DATA_Zre"], rgrp["DATA_-Zim"], c=_c, label=rpm
                            )
                            ax[0][0].plot(
                                rgrp["FIT_Zre"], rgrp["FIT_-Zim"], c=_c, label=rpm
                            )
                            ax[1][0].scatter(
                                rgrp["DATA_Yre"], rgrp["DATA_Yim"], c=_c, label=rpm
                            )
                            ax[1][0].plot(
                                rgrp["FIT_Yre"], rgrp["FIT_Yim"], c=_c, label=rpm
                            )

                        ax[0][0].legend(loc=0)
                        plt.suptitle(
                            ", ".join(str(i) for i in EV) + " " + Path(pf).stem
                        )
                        plt.savefig(
                            PDD_rpm.joinpath(
                                "_".join(str(i) for i in EV) + Path(pf).stem + ".png"
                            ),
                            bbox_inches="tight",
                        )
                        plt.close()

                # EV_corr_dct = EVgr.corr().loc[rpm_serie_slice].stack()['RPM_DAC'][[i for i in rpm_serie_slice if i != 'RPM_DAC']].to_dict()
                # EVidx = dict(zip(rpm_ser_idx_grp,EV))
                # EVidx.update(EV_corr_dct)
                # rpm_series_test_corr_lst.append(EVidx)

        #                print(EV)
        #                print(EVgr)
        # ====== CHECK RPM versus tau ==========
        _j12grp = "JOS12", "O2", 1.0, 0.55, _models[2]
        rpmgrp = (
            CB_EIS_rpm.groupby(rpm_ser_idx_grp)
            .get_group(_j12grp)
            .dropna(axis=1, how="all")
        )
        for _j12grp, rpmgrp in CB_EIS_rpm.groupby(rpm_ser_idx_grp + ["PAR_file"]):
            _j12grp = (*_j12grp[0:-1], (Path(_j12grp[-1]).stem))
            # rpmgrp.query('Aw < 3E3 & Aw > 1E-02').plot(x='RPM_DAC', y = 'tau', title=', '.join(str(i) for i in _j12grp), kind='scatter')
            if not rpmgrp.query("RPM_DAC > 100 & RPM_DAC < 1200").empty:
                _lin = linregress(
                    rpmgrp.query("RPM_DAC > 100").RPM_sqrt,
                    rpmgrp.query("RPM_DAC > 100").tau,
                )
                rpmgrp["tau_lin"] = rpmgrp.RPM_sqrt * _lin.slope + _lin.intercept
                if _lin.slope and _lin.rvalue > 0.1:
                    # and _lin.rvalue > 0.5 and _lin.rvalue < 1.0:
                    fig, ax = plt.subplots(figsize=(8, 8))
                    rpmgrp.query("Aw < 3E3 & Aw > 1E-02").plot(
                        x="RPM_sqrt",
                        y="tau",
                        title=", ".join(str(i) for i in _j12grp),
                        kind="scatter",
                        ax=ax,
                    )
                    rpmgrp.query("Aw < 3E3 & Aw > 1E-02").plot(
                        x="RPM_sqrt",
                        y="tau_lin",
                        label=f"lin: {_lin.slope:.3E} + {_lin.intercept:.3E}, r={_lin.rvalue:.3f}",
                        kind="line",
                        ax=ax,
                    )
                    ax.legend(fontsize=16, fancybox=True)
                    plt.savefig(
                        PDD_rpm.joinpath(
                            "RPM_lin_tau_" + "_".join(str(i) for i in _j12grp) + ".png"
                        )
                    )

        EIS_rpm_corrs = pd.DataFrame(rpm_series_test_corr_lst)
        for par in eisplot.parlst + eisplot.extrapars:
            if par in EIS_rpm_corrs.columns:
                rpm_par_corrs = (
                    EIS_rpm_corrs[rpm_ser_idx_grp + [par]].sort_values(par).dropna()
                )
                top_corrs = pd.concat(
                    [rpm_par_corrs.head(7), rpm_par_corrs.tail(7)]
                ).set_index(rpm_ser_idx_grp[0:-1])

                PDD_par = PDD_rpm.joinpath("{0}".format(par))
                PDD_par.mkdir(parents=True, exist_ok=True)
                for n, row in top_corrs.iterrows():

                    fig, ax = plt.subplots(figsize=(8, 8))
                    prop_cycle = plt.rcParams["axes.prop_cycle"]
                    colors = prop_cycle.by_key()["color"]
                    ax.set_prop_cycle(prop_cycle)
                    for (modname, modindex), modgr in CB_EIS_rpm_grp_mods.get_group(
                        n
                    ).groupby(["Model_EEC", "Model_index"]):
                        if not modgr[par].dropna().empty:
                            corr_val = np.round(
                                EIS_rpm_corrs.set_index(rpm_ser_idx_grp).loc[
                                    n + (modname,), par
                                ],
                                2,
                            )
                            modgr.plot(
                                y=par,
                                x="RPM_DAC",
                                s=100,
                                alpha=0.6,
                                c=colors[modindex],
                                kind="scatter",
                                title="{0}".format(n),
                                label=modname + f" ({corr_val})",
                                ax=ax,
                            )
                            ylim, logytrue, logscale = eisplot(par).pars_plot_lims()
                            ax.set_ylim(ylim[0], ylim[1])
                            #                                     modgr[par].mean()+3*modgr[par].std()))
                            ax.set_yscale(logscale)
                    ax.legend(loc="best", fontsize=8)
                    plt.savefig(
                        PDD_par.joinpath(
                            "{0}.png".format(
                                "_".join([str(i) for i in [par] + list(n)])
                            )
                        ),
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close()
            else:
                print(f"passed {par}")

        #        for rpmn, rpmgr in CB_EIS_rpm.query('RPM_DAC < 3000 & pH < 15').groupby(['SampleID', 'Gas', 'pH']):
        #            for spf, spfgr in EVgr.groupby('SpectraFile'):
        #                        spectra_combined = pd.read_excel(spf,index_col=[0])
        #                        PAR_EIS_fit_V2.EIS_plotting_EvRHE(spectra_combined,EVgr, EVgr, PDD_rpm.joinpath('EIS_combined_{0}.jpg'.format('_'.join([str(i) for i in rpmn]))))
        #                PAR_EIS_fit_V2.EIS_plotting_EvRHE(AllData_E_file,EISovv, Pars, PDD_rpm.joinpath('EIS_combined.jpg'))

        for modEV in EVgr, groupb("Model_EEC"):
            modname = modEV
            #            for modname in fast_checking_EEC_models[:
            PDD_rpm_mod = PDD_rpm.joinpath("{0}".format(modname))
            PDD_rpm_mod.mkdir(parents=True, exist_ok=True)
            par_lst = [
                i
                for i in modEVgr.dropna(axis=1).columns
                if i in SampleSelection.EC_EIS_par_cols + ["Rorr_kin", "R3_kin"]
            ]
            [i for i in par_lst if ["Rorr"]]
            # TODO :Add spectra and fits to see if fits are oke
            #                    spfgr
            for par in par_lst:
                if len(EVgr[par].dropna()) > 2 and modEVgr.RPM_DAC.nunique() > 2:
                    fig, ax = plt.subplots(figsize=(8, 8))
                    modEVgr.plot(
                        y=par,
                        x="RPM_DAC",
                        c="E_RHE",
                        colormap="rainbow_r",
                        s=100,
                        alpha=0.7,
                        kind="scatter",
                        title="{0}, {1}".format(modname, EV),
                        ax=ax,
                    )
                    ylim, logytrue, logscale = eisplot(par).pars_plot_lims()

                    ax.set_ylim((ylim[0], EVgr[par].mean() + 3 * modEVgr[par].std()))
                    ax.set_yscale(logscale)
                    plt.savefig(
                        PDD_rpm_mod.joinpath(
                            "{0}_{1}.png".format(
                                "_".join([str(i) for i in EV[1:]]), par
                            )
                        ),
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close()
        #                EVgr.loc[EVgr['Model_EEC'] == modname].query('pH < 15').plot(y='Rs',x='RPM_DAC',kind='scatter',ylim=(0,100),title='{0}, {1}'.format(modname,EV))
        #                EVgr.loc[EVgr['Model_EEC'] == modname].query('pH < 15').plot(y='Rct',x='RPM_DAC',kind='scatter',ylim=(0,100),title=modname)
        CB_EIS_rpm.loc[CB_EIS_rpm["Model_EEC"] == modname].query("pH < 15").plot(
            y="nAd",
            x="E_RHE",
            c="BET_cat_agg",
            colormap="rainbow_r",
            kind="scatter",
            ylim=(0, 1),
            title=modname,
        )
        CB_EIS_rpm.loc[CB_EIS_rpm["Model_EEC"] == modname].query("pH < 15").plot(
            y="nDL",
            x="E_RHE",
            c="BET_cat_agg",
            colormap="rainbow_r",
            kind="scatter",
            ylim=(0, 1),
            title=modname,
        )
        CB_EIS_rpm.loc[CB_EIS_rpm["Model_EEC"] == modname].query("pH < 7").plot(
            y="Rct",
            x="E_RHE",
            c="BET_cat_agg",
            colormap="rainbow_r",
            kind="scatter",
            ylim=(1, 1e6),
            logy=True,
            title=modname,
        )
        CB_EIS_rpm.loc[CB_EIS_rpm["Model_EEC"] == modname].query("pH > 7").plot(
            y="Qad+Cdlp",
            x="E_RHE",
            c="BET_cat_agg",
            colormap="rainbow_r",
            kind="scatter",
            ylim=(0.1, 1e-4),
            logy=True,
            title=modname,
        )


def check_recent_pars(EIS_pars):
    EIS_pars.source_delta_mtime
    _fp = FindExpFolder("VERSASTAT").PostDir.joinpath("EIS_pars_nonrecent")

    plt.plot(EIS_pars.source_delta_mtime.unique())
    plt.xlim(2500, 3000)
    plt.ylim(0, 0.5e15)

    EIS_pars_nonrecent = EIS_pars[
        EIS_pars.source_delta_mtime > pd.Timedelta(2e14, unit="ns")
    ]
    EIS_pars_recent = EIS_pars[
        EIS_pars.source_delta_mtime < pd.Timedelta(2e14, unit="ns")
    ]
    EIS_pars_nonrecent_pf = EIS_pars.loc[
        EIS_pars.PAR_file.isin(
            [i for i in EIS_pars_nonrecent.PAR_file.unique() if not "fakeZmean" in i]
        )
    ]
    EIS_pars_nonrecent_pf.to_pickle(_fp)

    len([i for i in EIS_pars_recent.PAR_file.unique() if not "fakeZmean" in i])


def check_overall_features(EIS_pars):
    _fp = FindExpFolder("VERSASTAT").PostDir.joinpath("EIS_features")
    _fp.mkdir()

    EIS_pars.Model_EEC.unique()

    EISmod = EIS_pars.loc[eval(EIS_filter())].loc[
        (EIS_pars.lmfit_chiqsr < 5)
        & (EIS_pars.E_RHE > 0.5)
        & (EIS_pars.E_RHE < 0.8)
        & (EIS_pars.RPM_DAC > 1000)
        & (EIS_pars.Loading_name == "standard")
    ]

    ch_col = ["BET_cat_agg"]
    for n, gr in EISmod.groupby(["RPM_DAC", "Gas", "pH"]):
        n, gr
        gr_mod = gr.loc[gr.Model_EEC == mod_select[n[1]]]
        _txt = ", ".join([str(i) for i in n])
        mod_varnames = gr_mod.lmfit_var_names.unique()[0].split(", ")
        for ch in ch_col:
            for vn in mod_varnames:
                gr_mod = gr_mod.assign(
                    **{vn + "_stderr_rel": gr_mod[vn] / gr[vn + "_stderr"]}
                )
                gr_err = gr_mod.loc[
                    (gr[vn + "_stderr"] > 0) & (gr_mod[vn + "_stderr_rel"] < 1e2)
                ].dropna(subset=[ch])
                corr_vn = gr_err.corr()[vn]
                corr_vn.loc[np.abs(corr_vn) > 0.5]
                PFagg = gr_err.groupby("PAR_file")[[vn, ch]].agg("mean")

                if PFagg[ch].std() > 2:
                    PFagg.plot(x=ch, y=vn, kind="scatter", title=_txt)


def check_best_models():
    _err_type = "lmfit_MSE"
    _filter = EIS_filter()
    _filter += '& (EIS_pars.SampleID.str.contains("JOS1|JOS2|JOS3|JOS4|JOS5"))'
    _filter += "& (EIS_pars.EIS_fake == False)"
    _grps = ["Model_EEC", "Gas", "lmfit_var_names"][0:2]
    best_models = (
        EIS_pars.loc[eval(_filter)]
        .dropna(axis=0, subset=[_err_type])
        .groupby(_grps)[_err_type]
        .agg(["count", "mean", "std"])
        .sort_values("mean", ascending=True)
    )
    print(best_models)
    keep_models = (
        best_models.loc[(best_models["count"] > 5) & (best_models["std"] > 0)]
        .index.get_level_values(0)
        .unique()
    )
    EIS_pars = EIS_pars.loc[EIS_pars.Model_EEC.isin(keep_models)]
    best_models = (
        EIS_pars.loc[eval(_filter)]
        .dropna(axis=0, subset=[_err_type])
        .groupby(_grps)[_err_type]
        .agg(["count", "mean", "std"])
        .sort_values(["Gas", "mean"], ascending=True)
    )
    print(best_models)


class EIS_Warburg_lin:
    def __init__(self, EIS_pars, model_select="Model(R0-L0-p(R1-Wo1,CPE1)-C2)"):

        self.EIS_pars = EIS_pars
        self.filter_pars()
        self.model_select = model_select
        self.select_mod()

    def filter_pars(self):
        self.EIS_pars_raw = self.EIS_pars
        self.EIS_pars = self.EIS_pars.query(f"{EIS_selection.filter}")

    def select_mod(self):
        if self.model_select in self.EIS_pars.Model_EEC.unique():
            self.EIS_mod = self.EIS_pars.loc[
                self.EIS_pars.Model_EEC == self.model_select
            ]
            self.var_names = self.EIS_mod.lmfit_var_names.unique()[0].split(", ")

    def map_gas_values(self):
        _gas_map = {"O2": 1, "N2": 0}
        self.EIS_pars["Gas_n"] = self.EIS_pars["Gas"].map(_gas_map)

    def check_lin(self):

        c = "Gas"
        [i for i in self.EIS_mod.columns if "lin" in i and i.startswith("WB_")]
        _y = "WB_Zre_lin_slopes_slope"
        for par in self.var_names + [_y]:
            for (g, ser), ggrp in self.EIS_mod.query(
                "(WB_Zre_lin_slopes_rvalue > 0.8) & (RPM_DAC > 1000)"
            ).groupby(["Gas", "SeriesID"]):
                if not ggrp[_y].drop_duplicates().dropna(axis=0).empty:
                    ggrp.drop_duplicates(subset=_y).plot(
                        x="E_RHE",
                        y=par,
                        kind="scatter",
                        c="BET_cat_agg",
                        colormap="rainbow_r",
                        title=f"{ser}, {g}",
                    )


def _force_reload():
    EIS_pars = Load_from_Indexes.EIS_pars_OVV(
        reload=True, use_latest=True, use_daily=True
    )
    return EIS_pars


if __name__ == "__main__":

    mc = Model_Collection()
    _T = True
    _F = False

    try:
        self = EIS_Warburg_lin(EIS_pars)
    except Exception as e:
        print(f"EIS export: {e}")
        EIS_pars = Load_from_Indexes.EIS_pars_OVV(
            reload=True, use_latest=True, use_daily=True
        )
        # self = EIS_Warburg_lin(EIS_pars)
    # export_to_Origin(EIS_pars)
