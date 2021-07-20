#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:59:48 2021

@author: zmg
"""

import sys
from pathlib import Path
import functools

# import collections
from collections import Counter
import pickle

# import types
# import post_helper
# import plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.stats import linregress
import pandas as pd
import numpy as np
import datetime as dt


mpl.style.use("seaborn")
mpl.rcParams["figure.dpi"] = 100
# from sklearn.cluster import KMeans


# print ('Name prepare input:', __name__ )
if __name__ == "__main__":
    #    print(f'Package: {__package__}, File: {__file__}')
    #    FH_path = Path(__file__).parent.parent.parent.joinpath('FileHelper')
    #    sys.path.append(str(FH_path))
    #    sys.path.append(str(Path(__file__).parent.parent.joinpath('indexer')))
    sys.path.append(str(Path(__file__).parent.parent.parent))
    #    sys.path.append("..")
    #    print(sys.path)
    #    import FileHelper
    from FileHelper.PostChar import Characterization_TypeSetting, SampleCodesChar
    from FileHelper.PostPlotting import *
    from FileHelper.FindSampleID import GetSampleID
    from FileHelper.FindFolders import FindExpFolder

    #    from FileHelper.FileFunctions.FileOperations import PDreadXLorCSV
    from collect_load import Load_from_Indexes

    #    from FileHelper.FindExpFolder import FindExpFolder
    from plotting import eisplot
    import EIS_export
elif "prepare_input" in __name__:
    pass
else:
    from FileHelper.PostChar import Characterization_TypeSetting, SampleCodesChar
    from FileHelper.PostPlotting import *
    from FileHelper.FindSampleID import GetSampleID
    from FileHelper.FindFolders import FindExpFolder
#    import RunEC_classifier
#    from FileHelper.FindSampleID import FindSampleID

# from FileHelper.PostChar import SampleSelection, Characterization_TypeSetting
def mkfolder(folder):
    folder.mkdir(exist_ok=True, parents=True)
    return folder


def filter_cols(_df, n):
    if any(["startswith" in i for i in n]):
        _lst = [i for i in _df.columns if i.startswith(n[-1])]
    else:
        _lst = [i for i in _df.columns if n[-1] in i]
    return _lst


OriginColors = Characterization_TypeSetting.OriginColorList()
Pfolder = FindExpFolder().TopDir.joinpath(
    Path("Preparation-Thesis/SiO2_projects/SiO2_Me_ECdepth+LC")
)
plotsfolder = mkfolder(Pfolder.joinpath("correlation_plots/all_samples"))
EC_folder = Pfolder.joinpath("EC_data")
# EC_index, SampleCodes = Load_from_Indexes.get_EC_index()
print("finished")

SampleCodesChar = SampleCodesChar().load


def multiIndex_pivot(df, index=None, columns=None, values=None):
    # https://github.com/pandas-dev/pandas/issues/23955
    output_df = df.copy(deep=True)
    if index is None:
        names = list(output_df.index.names)
        output_df = output_df.reset_index()
    else:
        names = index
    output_df = output_df.assign(
        tuples_index=[tuple(i) for i in output_df[names].values]
    )
    if isinstance(columns, list):
        output_df = output_df.assign(
            tuples_columns=[tuple(i) for i in output_df[columns].values]
        )  # hashable
        output_df = output_df.pivot(
            index="tuples_index", columns="tuples_columns", values=values
        )
        output_df.columns = pd.MultiIndex.from_tuples(
            output_df.columns, names=columns
        )  # reduced
    else:
        output_df = output_df.pivot(
            index="tuples_index", columns=columns, values=values
        )
    output_df.index = pd.MultiIndex.from_tuples(output_df.index, names=names)
    return output_df


def get_float_cols(df):
    return [key for key, val in df.dtypes.to_dict().items() if "float64" in str(val)]


def cm2inch(value):
    return value / 2.54


# class PorphSamples():
#    def __init__(self):
#        self.template = PorphSamples.template()
def decorator(func):
    @functools.wraps(func)
    def wrapper_decorator(*args, **kwargs):
        # Do something before
        value = func(*args, **kwargs)
        # Do something after
        return value

    return wrapper_decorator


def read_load_pkl(_pklstem):
    _pklpath = EC_PorphSiO2.folder.joinpath(_pklstem).with_suffix(".pkl")
    if _pklpath.exists():
        try:
            print("pkl reloading:", _pklpath)
            DF_diff = pd.read_pickle(_pklpath)
            DF_diff.columns
            return DF_diff
        except Exception as e:
            print("reading error", e)
            return pd.DataFrame()
    else:
        print("read error not existing", _pklpath)
        return pd.DataFrame()


def save_DF_pkl(_pklstem, _DF):
    _pklpath = EC_PorphSiO2.folder.joinpath(_pklstem).with_suffix(".pkl")
    try:
        print("pkl saving to:", _pklpath)
        _DF.to_pickle(_pklpath)
    except Exception as e:
        print("pkl saving error", e, _pklpath)
    return _pklpath


def load_dict_pkl(_pklstem):
    _pklpath = EC_PorphSiO2.folder.joinpath(_pklstem).with_suffix(".pkl")
    if _pklpath.exists():
        try:
            print("pkl reloading:", _pklpath)
            with open(_pklpath, "rb") as file:
                _dict = pickle.load(file)
            return _dict
        except Exception as e:
            print("reading error", e)
            return {}
    else:
        print("read error not existing", _pklpath)
        return {}


def save_dict_pkl(_pklstem, _dict):
    _pklpath = EC_PorphSiO2.folder.joinpath(_pklstem).with_suffix(".pkl")
    try:
        print("pkl saving to:", _pklpath)
        with open(_pklpath, "wb") as file:
            pickle.dump(_dict, file)
    except Exception as e:
        print("pkl saving error", e, _pklpath)
    return _pklpath


def PorphSiO2_template():
    #        'SerieIDs' : ('Porph_SiO2')*5,
    Series_Porph_SiO2 = {
        "SampleID": ("JOS1", "JOS2", "JOS3", "JOS4", "JOS5"),
        "Metal": ("Fe", "Co", "MnTPP", "FeTPP", "H2"),
        "color": (2, 4, 6, 15, 3),
    }

    Porphyrins = {
        "TMPP": {"Formula": "C48H38N4O4", "MW": 734.8382},
        "TMPP-Fe(III)Cl": {"Formula": "C48H36ClFeN4O4", "MW": 824.1204},
        "TMPP-Co(II)": {"Formula": "C48H36CoN4O4", "MW": 791.7556},
        "TTP-Mn(III)Cl": {"Formula": "C44H28ClMnN4", "MW": 703.1098},
        "TPP-Fe(III)Cl": {"Formula": "C44H28ClFeN4", "MW": 704.0168},
        "TPP": {"Formula": "C44H30N4", "MW": 614.7346},
    }

    Porph_template = pd.DataFrame(Series_Porph_SiO2)
    return Porph_template


def SeriesID_mapper():
    _markers = {
        "CB3": "v",
        "CB7": "<",
        "CB4": "^",
        "CB6": ">",
        "CB_pyr": "s",
        "CB_bim": "d",
        "Porph_SiO2": "o",
        "Co_PANI": "P",
        "C3N4-Plu": "*",
        "AS27": "X",
    }
    return _markers


def SampleColor_mapper():
    _samplecolor = {
        i[0]: [int(f) / 255 for f in i[1].split(",")]
        for i in OriginColors[["OriginCode", "RGB"]].to_numpy()
    }
    return _samplecolor


# OriginColors[['OriginCode','RGB']].to_dict()


def _postchar_test():
    self = postChar()


class postChar:

    suffixes = ["R", "P", "EA", "B"]

    template = SampleCodesChar[list(SampleCodesChar.columns)[0:11]]

    def __init__(self, folder=""):

        self.folder = folder
        #        FindExpFolder().TopDir.joinpath(Path('Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc'))
        self.raman = self.postRaman
        self.bet = self.postBET
        self.ea = self.postEA
        self.prep = self.postPrep
        # self.prec_ea_ratios()
        self.merge_slice()
        # self.merged = self.merged()

    def merge_slice(self):
        # cols = ['SampleID']
        # list(PorphSiO2_template().columns)
        _rb = set(self.raman.columns).intersection(set(self.bet.columns))
        _rbea = set(_rb).intersection(self.ea.columns)
        _rbeap = set(_rbea).intersection(self.prep.columns)
        cols = list(_rbeap)
        ra_bet = pd.merge(self.raman, self.bet, on=cols, how="outer")
        ea_ra_bet = pd.merge(ra_bet, self.ea, on=cols, how="outer")
        p_ea_ra_bet = pd.merge(ea_ra_bet, self.prep, on=cols, how="outer")
        self.merged_raw = p_ea_ra_bet
        self.merged = self.merged_raw.loc[
            self.merged_raw.SampleID.isin(self.template.SampleID.unique())
        ]
        # return merged_out

    def decorator(func):
        @functools.wraps(func)
        def wrapper_decorator(*args, **kwargs):
            # Do something before
            value = postChar.slice_out(
                *args, slice_lst=postChar.template.SampleID, template=postChar.template
            )
            # Do something after
            return value

        return wrapper_decorator

    def slice_out(func, template=SampleCodesChar[list(SampleCodesChar.columns)[0:11]]):
        pars, suffx = func()
        try:
            slice_lst = template.SampleID
            pars_out = pars.loc[pars.SampleID.isin(slice_lst)]
            pars_out = pd.merge(template, pars_out, on="SampleID")
            if suffx != "":
                cols = [i for i in pars_out.columns if i not in template.columns]
                pars_out = pars_out.rename(
                    columns=dict(zip(cols, [f"{suffx}_" + i for i in cols]))
                )
        except:
            pars_out = pars
            print("No slice out for {func.__name__} ")
        return pars_out

    @slice_out
    def postRaman(peak_model="5peaks", plot_resid=True):
        ramandir = FindExpFolder("RAMAN").DestDir
        raman_files = ramandir.rglob("*xlsx")
        fitpars_fls = [
            i
            for i in raman_files
            if all([a in i.stem for a in ["FitParameters_Model", "peaks"]])
            and not "_Conflict" in i.stem
            and not "Porph_SiO2" in i.parts
        ]
        FitPars_raw = pd.concat(
            [
                pd.read_excel(i).assign(
                    **{"Model": i.stem.split("Model_")[-1], "SourceFile": i}
                )
                for i in fitpars_fls
            ],
            sort=False,
        )
        uniq_cols = ["D_center", "G_center", "G_amplitude", "D_amplitude"]
        FitPars = FitPars_raw.drop_duplicates(
            subset=["Model", "SampleID"] + uniq_cols
        ).sort_values(by=["Model", "SampleID"])
        # FitPars_raw.loc[FitPars_raw.SampleID.isin(PorphSiO2_template().SampleID.values)]
        if plot_resid:
            Model_redchi = pd.DataFrame(
                [
                    (
                        n,
                        np.abs(gr.redchi).mean(),
                        np.abs(gr.redchi).sum(),
                        gr.aic.mean(),
                    )
                    for n, gr in FitPars.groupby("Model")
                ],
                columns=["Model", "redchi_mean", "redchi_sum", "aic_mean"],
            ).set_index("Model")
            Model_chi = pd.DataFrame(
                [
                    (n, gr.chisqr.mean(), gr.chisqr.sum())
                    for n, gr in FitPars.groupby("Model")
                ],
                columns=["Model", "chi_mean", "chi_sum"],
            ).set_index("Model")
            Model_redchi.aic_mean.plot.bar()

            Model_chi.chi_mean.plot.bar()
        if peak_model:
            FPars_out_1st = FitPars.loc[FitPars.Model.isin([peak_model])].dropna(axis=1)
        else:
            FPars_out_1st = FitPars
        _error_cols = ["chisqr", "redchi", "bic", "aic", "nfev"]
        t2nd_mod = "2ndOrder_4peaks"
        if t2nd_mod in FitPars_raw.Model.unique():
            FPars_out_2nd = FitPars.loc[FitPars.Model == t2nd_mod].dropna(axis=1)
            flt2nd = [
                i
                for i in FPars_out_2nd.select_dtypes(include=["int64", "float"]).columns
                if i not in _error_cols
            ]
            get_float_cols(FPars_out_2nd)

            FPars_out_2nd = FPars_out_2nd.rename(
                columns=dict(zip(flt2nd, [f"2nd_{i}" for i in flt2nd]))
            )
            FPars_out = pd.merge(
                FPars_out_1st,
                FPars_out_2nd,
                on=["SampleID"],
                how="left",
                suffixes=["_1st", "_2nd"],
            )
        else:
            FPars_out = FPars_out_1st
        return FPars_out, "R"

    @slice_out
    def postBET():
        # betdir = FindExpFolder('PorphSiO2').folder.joinpath('SiO2_Me_EC+Struc/BET')
        betdir = FindExpFolder("BET").DestDir
        bet_files = betdir.rglob("BET_pars_index*pkl")
        BET_ovv = pd.concat([pd.read_pickle(i) for i in bet_files])
        BET_ovv_template = BET_ovv
        # .loc[BET_ovv.SampleID.isin(PorphSiO2_template().SampleID.unique())]
        BET_ovv_template = BET_ovv_template.loc[
            BET_ovv_template.B_RAW_ads_BET_C_constant >= 0
        ]
        BET_ovv_template = BET_ovv_template.loc[
            BET_ovv_template.fullPath.str.endswith("RAW")
            & ~BET_ovv_template.fullPath.str.contains("_FAIL")
        ]
        return BET_ovv_template, ""

    @slice_out
    def postEA():
        EAcols = [
            "SampleID",
            "C/N_ratio",
            "N_content",
            "C_content",
            "H_content",
            "100-CHN",
        ]
        EA_results = pd.read_excel(list(FindExpFolder("EA").DestDir.rglob("*xlsx"))[0])[
            EAcols
        ]
        EA_results["C/H_ratio"] = EA_results["C_content"] / EA_results["H_content"]
        return EA_results, "EA"

    @slice_out
    def postPrep():
        Prep_Porph_SiO2 = {
            "SampleID": ("JOS1", "JOS2", "JOS3", "JOS4", "JOS5"),
            "WL_precursor": (31.0, 56.0, 53.0, 38.0, 41.0),
            "Precursor_type": (
                "FeTMPPCl",
                "CoTMPPCl",
                "MnTPPCl",
                "FeTPPCl",
                "H2TMPPCl",
            ),
            "MW_g-mol": (824.12, 791.67, 703.11, 704.02, 734.84),
            "Metal_element": (26.0, 27.0, 25.0, 26.0, 1.0),
        }
        Prep_Porph = pd.DataFrame(Prep_Porph_SiO2)
        Prep_Porph["prec_N"] = 100 * (4 * 14) / Prep_Porph["MW_g-mol"]
        Prep_Porph["prec_C"] = (
            100 * (12 * np.array([48, 48, 44, 44, 48])) / Prep_Porph["MW_g-mol"]
        )
        Prep_Porph["prec_H"] = (
            100 * (1 * np.array([36, 36, 28, 28, 38])) / Prep_Porph["MW_g-mol"]
        )
        Prep_Porph["prec_100-CHN"] = 100 - (
            Prep_Porph["prec_H"] + Prep_Porph["prec_N"] + Prep_Porph["prec_C"]
        )
        Prep_Porph["prec_C/N_ratio"] = Prep_Porph["prec_C"] / Prep_Porph["prec_N"]

        return Prep_Porph, "P"

    def prec_ea_ratios(self):
        templ = PorphSiO2_template()
        ea = self.ea
        prep = self.postPrep
        # self.prep = self.prec_ea_ratios(, self.ea, postChar.postPrep)
        pcols = prep.columns
        _ms = []
        for eacol in [i for i in ea.columns if i not in templ.columns]:
            if eacol.endswith("_content"):
                _mcol = [
                    i for i in pcols if i.split("P_prec_")[-1] in eacol.split("EA_")[-1]
                ]
            else:
                _mcol = [
                    i for i in pcols if i.split("P_prec_")[-1] == eacol.split("EA_")[-1]
                ]

            if _mcol:
                _ms.append((eacol, _mcol[0]))
        _coldict = {}
        for mc in _ms:
            elem = mc[1].split("P_prec_")[-1]
            ratiocol = f"P_ratioEA_{elem}"
            prep = prep.assign(**{ratiocol: prep[mc[1]] / ea[mc[0]]})
            _coldict.update(
                {elem: {"ea_col": mc[0], "prep_col": mc[1], "ratio_col": ratiocol}}
            )
        print("Added Ratio cols to prep")
        self.prep = prep
        # return prep


def _testing():
    tt = test_corrs()
    self = tt
    self.prep_corrs()
    cc = self.corrs.stack().sort_values()


class test_corrs:
    def __init__(self):
        self.pc = postChar()
        self.pars = self.pc.merged
        self.charname = "test"
        self.prep_corrs()
        self.find_combos()

    def filter_char_pars(self):
        pars = self.pars
        _endsw = [
            "fit_R",
            "bic",
            "_aic",
            "_rmse",
            "_chisqr",
            "relPrange_min",
            "SAMPLE WEIGHT",
            "_stderr",
            "_pvalue",
            "_rvalue",
            "_redchi",
            "error(%)",
            "corr_coef",
            "_STATION",
            "nfev_2nd",
            "nfev_1st",
            "color",
            "_nmono",
            "_monolayer",
        ]
        _keepcols = [i for i in pars.columns if not any(i.endswith(a) for a in _endsw)]
        _contains = ["_bic", "_aic", "_rmse", "_chisqr", "_redchi", "_height", "_fwhm"]
        _extra = ["SizeNP", "_RPT_", "Colorcode"]
        _keepcols = [
            i for i in _keepcols if not any(a in i for a in _contains + _extra)
        ]

        _filteredcols = [i for i in pars.columns if i not in _keepcols]
        return pars[_keepcols]

    def prep_corrs(self, filter_min=0.0):
        pars = self.pars
        only_floats = [
            key for key, val in pars.dtypes.to_dict().items() if "float64" in str(val)
        ]
        parsfltr = self.filter_char_pars()
        _corrs = (
            parsfltr.dropna(thresh=3, axis=1)
            .corr(min_periods=6, method="pearson")
            .dropna(axis=0, how="all")
        )
        _crrstk = _corrs.stack()
        _crrstk = _crrstk.loc[_crrstk.abs() > filter_min]
        _filter = [i for i in _crrstk.index if str(i[0])[0] != str(i[1])[0]]
        _corrs = _crrstk.loc[_filter].unstack()
        self.corrs = _corrs
        # return _corrs1

    def plot_heatmap(self, charname):
        pars = self.pars

        _corrs = self.corrs

        fig, ax = plt.subplots(figsize=(24, 19))
        # plt.matshow(_corrs)
        # plt.matshow(_corrs, fignum=f.number)
        mask = np.triu(np.ones_like(_corrs, dtype=bool))
        # plt.xticks(range(pars.select_dtypes(['number']).shape[1]), pars.select_dtypes(['number']).columns, fontsize=14, rotation=45)
        # plt.yticks(range(pars.select_dtypes(['number']).shape[1]), pars.select_dtypes(['number']).columns, fontsize=14)

        # plt.yticks(np.arange(0.5, len(_corrs.index), 1), _corrs.index)
        # plt.xticks(np.arange(0.5, len(_corrs.columns), 1), _corrs.columns)
        cmap = sns.diverging_palette(230, 20, as_cmap=True)
        sns.heatmap(
            _corrs,
            mask=mask,
            cmap=cmap,
            vmax=1,
            center=0,
            square=True,
            linewidths=0.5,
            cbar_kws={"shrink": 0.5},
        )
        # ax = sns.heatmap(_corrs,  square=True, mask=np.zeros_like(_corrs, dtype=np.bool), ax=ax)
        plt.savefig(
            FindExpFolder("PorphSiO2").compare.joinpath(f"heatmap_{charname}.png"),
            bbox_inches="tight",
            dpi=300,
        )
        plt.close()

    def find_combos(self):
        _corrstk = self.corrs.stack()
        _highcorr = _corrstk.loc[_corrstk.abs() > 0.011]
        _highcorr.index
        len([set(i) for i in _highcorr.index])
        _setdex = set(
            [(", ".join(sorted(set(i))), _highcorr.loc[i]) for i in _highcorr.index]
        )
        [i for i in _setdex]
        _setdex = sorted(_setdex, key=lambda x: x[-1])
        self.setdex = _setdex

    def plot_setdex(self):

        _selx1 = [
            i for i in self.setdex if "B_RAW_ads_BET_SA" in i[0] and "EA_" in i[0]
        ]
        _selx2 = [
            i
            for i in self.setdex
            if "B_RAW_ads_BET_C_constant" in i[0] and "EA_" in i[0]
        ]

        _selx3 = [i for i in self.setdex if "C/N_ratio" in i[0]]

        _sel = self.setdex[:30] + self.setdex[-30:]
        self.setdex[-30:]

        SerIDmap = SeriesID_mapper()
        sColor = SampleColor_mapper()

        for si in _sel:
            xcol, ycol, corrval = (*si[0].split(", "), si[1])
            xstr = "".join([i for i in xcol if i != "/"])
            ystr = "".join([i for i in ycol if i != "/"])

            fig, ax = plt.subplots()
            _sprs = (
                self.pars[[xcol, ycol, "SampleID", "Colorcode", "SeriesID"]]
                .dropna(axis=0, how="any")
                .drop_duplicates()
            )
            for (ser, sID), sgrp in _sprs.groupby(["SeriesID", "SampleID"]):
                ax.scatter(
                    sgrp[xcol],
                    sgrp[ycol],
                    color=sColor.get(sgrp.Colorcode.unique()[0]),
                    marker=SerIDmap.get(ser),
                    label=f"{ser} {sID}",
                )
                # sgrp.plot(x=xcol,y=ycol, kind='scatter',ax=ax,title=corrval,marker=SerIDmap.get(ser),label=ser)
            ax.legend(
                ncol=4, fontsize=12, loc="upper left", bbox_to_anchor=(-0.05, 1.51)
            )
            ax.set_xlabel(xcol)
            ax.set_ylabel(ycol)
            plt.title(corrval)
            plt.savefig(
                plotsfolder.joinpath(f"{np.abs(corrval):.2f} {xstr}_{ystr}.png"),
                bbox_inches="tight",
            )

    def pair_plot(pars):
        pars = self.pars
        # _corrs = prep_corrs(pars)
        _corrs = self.corrs
        _corrsdex = _corrs.index
        ycol = [
            "B_RAW_ads_BET_C_constant",
            "B_RAW_ads_BET_SA",
            "R_ID3/(IG+ID2",
            "B_RAW_ads_B_RPT_iso_meta_Isotherm_P/P0_max",
            "R_D3_height",
        ][-1]
        xcol = ["EA_N_content", "SizeNP"][0]
        # _corrs[ycol]
        ycorrs = _corrs[ycol].dropna().sort_values()
        ycorrsdex = ycorrs.index
        # attend = sns.load_dataset("attention").query("subject <= 12")
        # fig,axes = plt.subplots(nrows=len(ycorrsdex), sharex=True)
        # for ax,y in zip(axes, ycorrsdex):
        #     sns.regplot(x=ycol, y=y, data=pars,ax=ax)

        xcols = list(ycorrsdex[0:5]) + list(ycorrsdex[-5:])

        fig, ax = plt.subplots(figsize=(24, 19))
        sns.pairplot(
            pars[xcols + [ycol, "SampleID"]],
            hue="SampleID",
            diag_kind="none",
            kind="reg",
        )
        plt.savefig(
            FindExpFolder("PorphSiO2").compare.joinpath(f"pairplot_{xcol}_{ycol}.png"),
            bbox_inches="tight",
            dpi=300,
        )
        plt.close()

        fig, ax = plt.subplots()
        _sprs = (
            pars[[xcol, ycol, "SampleID"]].dropna(axis=0, how="any").drop_duplicates()
        )
        _sprs.plot(x=xcol, y=ycol, kind="scatter")
        # sns.pairplot(_sprs, x_vars=[ycol],y_vars = [xcol]
        # ,hue='SampleID')
        plt.savefig(
            FindExpFolder("PorphSiO2").compare.joinpath(
                f"pairplot_xvar_{self.charname}_{ycol}.png"
            ),
            bbox_inches="tight",
            dpi=300,
        )
        plt.close()
