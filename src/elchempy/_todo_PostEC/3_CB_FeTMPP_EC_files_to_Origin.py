# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 11:33:55 2020

@author: User
"""
import sys
from pathlib import Path
import functools

# import collections
from collections import Counter

# import types
# import post_helper
# import plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.stats import linregress
import pandas as pd
import numpy as np
import datetime as dt


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
elif "prepare_input" in __name__:
    pass
#    import RunEC_classifier
#    from FileHelper.FindSampleID import FindSampleID

# from FileHelper.PostChar import SampleSelection, Characterization_TypeSetting
def mkfolder(folder):
    folder.mkdir(exist_ok=True, parents=True)
    return folder


OriginColors = Characterization_TypeSetting.OriginColorList()
Pfolder = FindExpFolder().TopDir.joinpath(Path("Preparation-Thesis/CB_projects"))
plotsfolder = mkfolder(Pfolder.joinpath("correlation_plots"))
EC_folder = Pfolder.joinpath("EC_data")
EC_index, SampleCodes = Load_from_Indexes.get_EC_index()
print("finished")

# SampleCodesChar().load


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


def PorphSiO2_template():
    #        'SerieIDs' : ('Porph_SiO2')*5,
    Series_Porph_SiO2 = {
        "SampleID": (
            "DW16",
            "DW17",
            "DW18",
            "DW19",
            "DW20",
            "DW21",
            "DW24",
            "DW25",
            "DW26",
            "DW27",
            "DW28",
            "DW29",
            "JOS12",
            "JOS13",
            "JOS14",
            "JOS15",
        )
    }
    #                         'Metal' : ('Fe','Co','MnTPP','FeTPP','H2'),
    #                         'color' : (2, 4, 6, 15, 3)}

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


def EC_types_grp():
    #    KL ['ORR_E_AppV_RHE', 'ORR_KL_E_AppV_RHE','Electrode']
    _basic_EC_cond = ["postAST_post", "Sweep_Type", "pH", "Loading_cm2"]
    _extra_EC_cond = {
        "N2CV": [],
        "ORR": ["RPM_DAC_uni"],
        "KL": ["Electrode", "ORR_E_AppV_RHE"],
        "EIS": ["E_RHE"],
        "HER": [],
        "OER": [],
    }
    _out = {key: _basic_EC_cond + val for key, val in _extra_EC_cond.items()}
    return _out


def save_EC_index_PorphSiO2(EC_index, EC_folder):
    _porph_index = EC_index.loc[EC_index.SampleID.isin(PorphSiO2_template().SampleID)]
    _porph_index.to_excel(EC_folder.joinpath("EC_index_PorphSiO2.xlsx"))
    # save_EC_index_PorphSiO2(EC_index, EC_folder)


class EC_PorphSiO2:
    folder = FindExpFolder("PorphSiO2").compare
    Porph_template = PorphSiO2_template()
    EIS_models = [
        "Model(EEC_Randles_RWpCPE)",
        "Model(EEC_2CPE)",
        "Model(EEC_2CPEpW)",
        "Model(EEC_RQ_RQ_RW)",
        "Model(EEC_RQ_RQ_RQ)",
        "Model(Randles_RQRQ)",
    ]
    ORR_pars_all = Load_from_Indexes.ORR_pars_OVV(reload=False, use_daily=True)
    mcols = [i for i in Load_from_Indexes.EC_label_cols if i not in ["PAR_file"]] + [
        "Sweep_Type"
    ]
    EC_idx_PorphSiO2 = pd.read_excel(list(EC_folder.rglob("*EC_index*"))[0])
    EC_idx_PorphSiO2 = EC_idx_PorphSiO2.assign(
        **{
            "PAR_date_day_dt": [
                dt.date.fromisoformat(np.datetime_as_string(np.datetime64(i, "D")))
                for i in EC_idx_PorphSiO2.PAR_date.to_numpy()
            ]
        }
    )

    #    ['Model(Singh2015_RQRQ)', 'Model(Singh2015_RQRQR)', 'Model(Bandarenka_2011_RQRQR)',
    #                  'Model(Singh2015_RQRWR)', 'Model(Randles_RQRQ)', 'Model(Singh2015_R3RQ)']
    #    model_select = EC_PorphSiO2.EIS_models[1]
    #    self = EC_PorphSiO2()
    def __init__(self):
        self.index, self.AST_days = EC_PorphSiO2.select_ECexps(EC_folder)

    #        = EC_PorphSiO2.get_AST_days()
    #        self.pars = EC_PorphSiO2.mergedEC()
    #        self.par_export = EC_OHC.to_excel(self.folder.joinpath('EC_ORR_HPRR.xlsx'))

    def get_AST_days():
        gr_idx = EC_PorphSiO2.EC_idx_PorphSiO2.groupby("PAR_date_day_dt")
        AST_days = []
        for n, gr in gr_idx:
            #            n,gr
            exps = gr.PAR_exp.unique()
            #                gr.PAR_date_day.unique()[0]
            if any(["AST" in i for i in exps]):
                #                print(n,exps)
                #                AST_days.append(n)
                if n + dt.timedelta(1) in gr_idx.groups.keys():
                    _post = gr_idx.get_group(n + dt.timedelta(1))
                    #                    print(n + dt.timedelta(1), gr_idx.get_group(n + dt.timedelta(1)))
                    AST_days.append((n, n + dt.timedelta(1)))
                else:
                    AST_days.append((n, n))
                    print(n + dt.timedelta(1), "grp missing")
        print(AST_days)
        return AST_days

    def select_ECexps(EC_folder):
        LC_idx_fp = list(EC_folder.rglob("*EC_index*"))[0].parent.joinpath(
            "LC_index.xlsx"
        )

        AST_days = EC_PorphSiO2.get_AST_days()

        if LC_idx_fp.exists():

            LC_fls = EC_PorphSiO2.EC_idx_PorphSiO2.loc[
                EC_PorphSiO2.EC_idx_PorphSiO2.PAR_date_day_dt.isin(
                    [i for a in AST_days for i in a]
                )
            ]
            LC_fls.to_excel(
                list(EC_folder.rglob("*EC_index*"))[0].parent.joinpath("LC_index.xlsx")
            )

        else:
            try:
                LC_fls = pd.read_excel(LC_idx_fp, index_col=[0])
            except Exception as e:
                print(f"Excel load fail: {e}\n,file: {LC_idx_fp}")
                LC_fls = pd.DataFrame()
        return LC_fls, AST_days

    def repr_index(self):
        PAR_exp_uniq = {grn: len(grp) for grn, grp in self.index.groupby("PAR_exp")}
        print(f"Len({len(self.index)},\n{PAR_exp_uniq}")

    def get_AST_matches(DF):
        LC_fls, AST_days = EC_PorphSiO2.select_ECexps(EC_folder)
        #        DF = ORR.drop_duplicates()
        #        DF = N2CV.drop_duplicates()
        #        DF = EIS.drop_duplicates()

        if "PAR_date_day_dt" not in DF.columns:
            DF = DF.assign(
                **{
                    "PAR_date_day_dt": [
                        dt.date.fromisoformat(
                            np.datetime_as_string(np.datetime64(i, "D"))
                        )
                        for i in DF.PAR_date.to_numpy()
                    ]
                }
            )
            DF.PAR_date_day_dt = pd.to_datetime(DF.PAR_date_day_dt)

        #        list((set(DF.columns).intersection(set(LC_fls.columns))).intersection(set(mcols) ))
        #        DF = pd.merge(DF,LC_fls,on=)
        _compare_cols = [
            i for i in ["SampleID", "pH", "Gas", "Loading_cm2"] if i in DF.columns
        ]
        _swp_rpm = [
            "Sweep_Type",
            "RPM_DAC_uni" if "RPM_DAC_uni" in DF.columns else "RPM_DAC",
        ]
        _coll = []
        _diffr = []
        for _pre, _post in [i for i in AST_days if len(i) == 2][::]:
            #            if len(_dates) == 2:
            #                _pre,_post = _dates
            #            elif (len_dates) == 1:
            _preslice = DF.loc[
                (DF.PAR_date_day == _pre.strftime("%Y-%m-%d")) & (DF.postAST == "no")
            ]
            pre = _preslice.groupby(_compare_cols)

            _postslice = DF.loc[
                (DF.PAR_date_day == _post.strftime("%Y-%m-%d")) & (DF.postAST != "no")
            ]
            post = _postslice.groupby(_compare_cols)

            _res = {"pre_PAR_date_day_dt": _pre, "post_PAR_date_day_dt": _post}
            print(_res, _postslice.postAST.unique())

            union = set(pre.groups.keys()).union(set(post.groups.keys()))

            matches = set(pre.groups.keys()).intersection(set(post.groups.keys()))

            _difference_pre = set(pre.groups.keys()).difference(set(post.groups.keys()))
            _difference_post = set(post.groups.keys()).difference(
                set(pre.groups.keys())
            )
            #            _diffr.append((_pre,_post,_difference_pre, _difference_post))

            #            if matches:
            for match in union:
                _res.update(dict(zip(_compare_cols, match)))
                _mgrpcols = ["PAR_file", "dupli_num", "postAST"]
                if match in matches:
                    _mpre = pre.get_group(match).groupby(_mgrpcols)
                    _mpost = post.get_group(match).groupby(_mgrpcols)
                elif match in _difference_pre:
                    _mpre = pre.get_group(match).groupby(_mgrpcols)
                    _mpost = pre.get_group(match).groupby(_mgrpcols)
                elif match in _difference_post:
                    _mpre = post.get_group(match).groupby(_mgrpcols)
                    _mpost = post.get_group(match).groupby(_mgrpcols)
                #                    print(_mpost.groups)

                for (_prePF, npr, _preAST), prgrp in _mpre:
                    _res.update(
                        {
                            "pre_dupli_num": npr,
                            "pre_PAR_file": _prePF,
                            "pre_postAST": _preAST,
                        }
                    )

                    for (_poPF, npo, _postAST), pogrp in _mpost:
                        _res.update(
                            {
                                "post_dupli_num": npo,
                                "post_PAR_file": _poPF,
                                "post_postAST": _postAST,
                                "dupli_num_combo": f"{npr}, {npo}",
                            }
                        )
                        if "postAST_sHA" in _postAST:
                            print(_res)

                        _pr1 = prgrp.groupby(_swp_rpm)
                        _po1 = pogrp.groupby(_swp_rpm)
                        _rpmswp_matches = set(_pr1.groups.keys()).intersection(
                            set(_po1.groups.keys())
                        )
                        for _m in _rpmswp_matches:
                            _res.update(dict(zip(_swp_rpm, _m)))
                            #                                print(_res)
                            _coll.append(_res.copy())

        AST_matches = pd.DataFrame(_coll)
        return AST_matches

    #                            prgrp,pogrp
    #                            prgrp.groupby(['Sweep_Type','RPM_DAC']).groups
    #                            prgrp['ORR_Jkin_min_700']-pogrp['ORR_Jkin_min_700']

    def mergedEC():

        mcols = [
            i for i in Load_from_Indexes.EC_label_cols if i not in ["PAR_file"]
        ] + ["Sweep_Type"]
        _mcols = [i for i in mcols if not i in ["Gas", "E_RHE"]]
        LC_fls, AST_days = EC_PorphSiO2.select_ECexps(EC_folder)

        EC_merged_dict = {}

        _reloadset = False

        template = PorphSiO2_template()
        HPRR = EC_PorphSiO2.HPRR()

        N2CV = EC_PorphSiO2().N2cv(reload=False, use_daily=True)
        N2_AST = EC_PorphSiO2.get_AST_matches(N2CV)
        N2_AST_diff = EC_PorphSiO2.compare_AST_pars(N2CV, N2_AST, reload=_reloadset)
        #        _DFtype = EC_PorphSiO2.sense_DF_type(N2CV)
        EC_merged_dict.update({"N2CV": N2_AST_diff})

        list(N2CV.columns)
        #        _renameN2 = {c : c.split('_')[-1] for c in [i for i in N2CV.columns if any([i.split('_')[-1] in mcols])]}
        #        N2CV = N2CV.rename(columns = _renameN2)

        ORR = EC_PorphSiO2().ORR_pars()
        ORR_AST = EC_PorphSiO2.get_AST_matches(ORR)

        #        ttpfs = ORR.loc[ORR.ORR_E_onset > 0.85].PAR_file.unique()

        ORR_AST_diff = EC_PorphSiO2.compare_AST_pars(ORR, ORR_AST, reload=_reloadset)
        EC_merged_dict.update({"ORR": ORR_AST_diff})
        #        _renameO2 = {c : c.split('_')[-1] for c in [i for i in ORR.columns if any([i.split('_')[-1] in mcols]) and not '_Ring' in i]}
        #        ORR = ORR.rename(columns = _renameO2)

        KL = EC_PorphSiO2().KL_pars()
        KL = KL.assign(**{"RPM_DAC": 0})
        KL_AST = EC_PorphSiO2.get_AST_matches(KL)
        KL_AST_diff = EC_PorphSiO2.compare_AST_pars(KL, KL_AST, reload=_reloadset)
        EC_merged_dict.update({"KL": KL_AST_diff})
        #        _renameKL = {c : c.split('_')[-1] for c in [i for i in KL.columns if any([i.split('_')[-1] in mcols]) and not '_Ring' in i]}
        #        KL = KL.rename(columns = _renameKL)

        EIS = EC_PorphSiO2.EIS_pars()
        EIS_AST = EC_PorphSiO2.get_AST_matches(EIS)
        EIS_AST_diff = EC_PorphSiO2.compare_AST_pars(EIS, EIS_AST, reload=_reloadset)
        EC_merged_dict.update({"EIS": EIS_AST_diff})
        #        _renameEIS = {c : c.split('_')[-1] for c in [i for i in EIS.columns if any([i.split('_')[-1] in mcols]) and not '_Ring' in i]}
        #        EIS = EIS.rename(columns = _renameEIS)

        return EC_merged_dict

    #        ECmerged = pd.merge(ORR,pd.merge(N2CV, EIS,on=_mcols),on=_mcols)
    #        EC_EIS = pd.merge(ECmerged,EIS,on=mcols)

    #        EC_OHN_merged = pd.merge(template, EC_EIS, on='SampleID')
    #        EC_PorphSiO2.export_to_xls(EC_OHN_merged)
    #        return EC_OHN_merged

    def _test_plots():
        # N2CV plot test
        _x = ["IndividualLabel", "SampleID"][0]

        _y = (["N2_Cdl_Fcm-2_600_diff_perc"], (-80, 100))
        _y = (["N2_Cdl_Fcm-2_600_pre", "N2_Cdl_Fcm-2_600_post"], (0, 50e-3))
        N2_pltqry = N2_AST_diff.groupby(
            ["postAST_post", "Sweep_Type", "pH", "Loading_cm2"]
        )
        #        .get_group(('postAST', 'cathodic',1.0))
        [
            gr.sort_values(by="IndividualLabel").plot.bar(
                x=_x, y=_y[0], ylim=_y[1], rot=60, title=", ".join([str(i) for i in n])
            )
            for n, gr in N2_pltqry
        ]  # .groupby(['postAST_post', 'Sweep_Type','pH'])

        N2_AST_diff.groupby(["postAST_post"]).plot(
            x="SampleID", y="N2_Cdl_Fcm-2_400_diff_abs", ylim=(-200, 200)
        )

        tt3 = EC_index.loc[
            (EC_index.PAR_date_day_dt.isin(set([a for i in AST_days for a in i])))
            & (EC_index.PAR_exp == "N2_act")
        ]

        # ORR
        list(ORR_AST_diff.columns)
        _y = "ORR_Frac_H2O2_500_diff_abs"
        _y = (["ORR_Frac_H2O2_500_pre", "ORR_Frac_H2O2_500_post"], (-0, 20))
        _y = (["ORR_Jkin_min_750_pre", "ORR_Jkin_min_750_post"], (-0, 3))
        _y = (["ORR_E_onset_pre", "ORR_E_onset_post"], (0.6, 0.95))

        ORR_pltqry = ORR_AST_diff.query("RPM_DAC_uni > 1000").groupby(
            ["postAST_post", "RPM_DAC_uni", "Sweep_Type", "pH"]
        )

        ORR_pltqry.plot(x=_x, y=_y[0], ylim=_y[1])
        [
            gr.sort_values(by="IndividualLabel").plot.bar(
                x=_x, y=_y[0], ylim=_y[1], title=", ".join([str(i) for i in n]), rot=60
            )
            for n, gr in ORR_pltqry
        ]

        # ORR KL
        _y = (["ORR_nElectrons_diff_perc"], (-0, 0.51))
        _y = (["ORR_nElectrons_pre", "ORR_nElectrons_post"], (0, 10))
        _grps = KL_AST_diff.groupby(
            ["postAST_post", "Electrode", "Sweep_Type", "pH", "ORR_E_AppV_RHE"]
        )

        _Etest = [
            (n, len(gr.dropna(subset=_y[0], how="any", axis=0)))
            for n, gr in _grps
            if n[1] == "KL_I_Disk" and len(gr) > 2
        ]
        _Etest = sorted(_Etest, key=lambda x: x[0][-1])
        E_KL = 0.675
        KL_qry = KL_AST_diff.query("ORR_E_AppV_RHE ==@E_KL").groupby(
            ["postAST_post", "Electrode", "Sweep_Type", "pH"]
        )

        [
            gr.sort_values(by="IndividualLabel").plot.bar(
                x=_x,
                y=_y[0],
                ylim=_y[1],
                title=", ".join([str(i) for i in n + (E_KL,)]),
                rot=60,
            )
            for n, gr in KL_qry
        ]

        # EIS
        _y = (["EIS_Rct_O2_pre", "EIS_Rct_O2_post"], (0, 600), (0.55, 0.65))
        EIS_qry = EIS_AST_diff.query(
            'E_RHE > @_y[2][0] & E_RHE < @_y[2][1] & Sweep_Type == "cathodic"'
        ).groupby(["postAST_post", "Sweep_Type", "pH", "E_RHE"])

        [
            gr.sort_values(by="IndividualLabel").plot.bar(
                x=_x, y=_y[0], ylim=_y[1], title=", ".join([str(i) for i in n]), rot=60
            )
            for n, gr in EIS_qry
        ]

        # _Dm.plot(x='E_RHE',y=['EIS_Rct_O2_pre', 'EIS_Rct_O2_post', 'EIS_Rct_O2_diff_abs']) # test plotting

    def EC_merge_postchar():
        EC_merged_dict = EC_PorphSiO2.mergedEC()
        postChars = postChar().merged
        _extracols = [i for i in SampleCodes.columns if not "SampleID" in i]

        _ECm_postchar = {}

        for _EC, _DF in EC_merged_dict.items():
            _initcols = _DF.columns
            _DF = _DF.dropna(axis=1, how="all")
            _DF = _DF.drop(columns=_DF.columns.intersection(_extracols))
            _DF = pd.merge(_DF, postChars, on="SampleID")
            _postcols = _DF.columns
            _ECm_postchar.update({_EC: _DF})
        return _ECm_postchar

    def compare_AST_pars(_DF, _AST_in, reload=False):
        #        _DF, _AST_in = EIS, EIS_AST
        #        _DF, _AST_in = N2CV, N2_AST
        #        _DF, _AST_in = ORR, ORR_AST
        #        _DF, _AST_in = KL, KL_AST
        #        _DF, _AST_in = args
        _DF = _DF.drop_duplicates()
        _DFtype = EC_PorphSiO2.sense_DF_type(_DF)
        _DFtype_prefix = _DFtype.split("_")[0]

        _pklpath = EC_PorphSiO2.folder.joinpath(f"AST_compared_pars_{_DFtype}.pkl")

        if _pklpath.exists() and not reload:
            try:
                print("AST compare reading:", _pklpath)
                DF_diff = pd.read_pickle(_pklpath)
                return DF_diff
            except Exception as e:
                return print("reading error", e)

        else:
            _prec = [i for i in _AST_in.columns if not i.startswith("post")]
            _precols = [
                i.split("pre_")[-1] if i.startswith("pre") else i for i in _prec
            ]

            _post = [i for i in _AST_in.columns if not i.startswith("pre")]
            _postcols = [
                i.split("post_")[-1] if i.startswith("post") else i for i in _post
            ]

            _dropnacols = set(_post + _prec)
            list(set(_prec).intersection(set(_post)))
            _AST = _AST_in.dropna(subset=_dropnacols, how="any")

            _DF_diff_out = []
            _errors = []
            for n, r in _AST.iterrows():

                #            EIS.query('  and  '.join([f'{k} == {repr(v)}' for k, v in _pred.items()]))
                _pred = dict(zip(_precols, r[_prec].to_dict().values()))
                _preQ = "  &  ".join(
                    [f"{k} == {repr(v)}" for k, v in _pred.items() if k in _DF.columns][
                        1:
                    ]
                )
                _Dpre = _DF.query(_preQ)

                _postd = dict(zip(_postcols, r[_post].to_dict().values()))
                _postQ = "  &  ".join(
                    [
                        f"{k} == {repr(v)}"
                        for k, v in _postd.items()
                        if k in _DF.columns
                    ][1:]
                )
                _Dpos = _DF.query(_postQ)

                #            pd.merge(_Dpre,_Dpos)
                _0 = [
                    (i, _Dpre[i].unique()[0])
                    for i in _Dpre.columns
                    if _Dpre[i].nunique() <= 1 and not i.startswith(_DFtype_prefix)
                ]
                _1 = [
                    (i, _Dpos[i].unique()[0])
                    for i in _Dpos.columns
                    if _Dpos[i].nunique() <= 1 and not i.startswith(_DFtype_prefix)
                ]

                _mcols = [
                    i[0]
                    for i in set(_0).intersection(set(_1))
                    if not i[0].startswith("dupli")
                ]
                _mcols.sort()

                _othercols = _DF.columns.difference(_mcols)
                t2 = _Dpre[_othercols]

                if "EIS" in _DFtype and all(
                    ["E_RHE" in i for i in [_Dpre.columns, _Dpos.columns]]
                ):
                    _mcols += ["E_RHE"]

                #                _Dm = pd.merge(_Dpre,_Dpos,on=_mcols + ['E_RHE'],suffixes=['_pre','_post'])
                elif "ORR" in _DFtype:
                    _KLcols = ["ORR_E_AppV_RHE", "ORR_KL_E_AppV_RHE", "Electrode"]
                    if all(i in _othercols for i in _KLcols):
                        _mcols += _KLcols
                #                _Dm = pd.merge(_Dpre, _Dpos, on = _mcols, suffixes = ['_pre','_post'])

                _Dm = pd.merge(_Dpre, _Dpos, on=_mcols, suffixes=["_pre", "_post"])

                _parcols = [
                    (i, i.replace("_pre", "_post"))
                    for i in _Dm.columns
                    if i.startswith(_DFtype_prefix)
                    and i.endswith("_pre")
                    and i.replace("_pre", "_post") in _Dm.columns
                ]

                for _c0, _c1 in _parcols:
                    try:
                        _diffabs = _Dm[_c0] - _Dm[_c1]

                        _diffperc = 100 * (_Dm[_c1] - _Dm[_c0]) / _Dm[_c0]

                        _Dm = _Dm.assign(
                            **{
                                _c0.split("_pre")[0] + "_diff_abs": _diffabs,
                                _c0.split("_pre")[0] + "_diff_perc": _diffperc,
                            }
                        )
                    except Exception as e:
                        #                   pass
                        _errors.append((_c0, _c1, e))
                _DF_diff_out.append(_Dm)
            #                    print(_c0, e)
            DF_diff = pd.concat(_DF_diff_out).drop_duplicates()
            DF_diff.to_pickle(_pklpath)
            print(f"AST compare len({len(DF_diff)}) saved to:{_pklpath}")
            return DF_diff

    # DF_diff.groupby(['postAST_post','SampleID']).plot(x='E_RHE', y='EIS_Rct_O2_diff_abs',ylim=(-200,200))

    def sense_DF_type(_DF):
        #        _c = [i[0] for i in Counter([i.split('_')[0] for i in _DF.columns]).most_common(5) if i[0] not in ['BET','tM']][0]
        _res = [
            i[0]
            for i in Counter(
                ["_".join(i.split("_")[0:2]) for i in _DF.columns]
            ).most_common(5)
            if not any([b in i[0] for b in ["BET", "tM"]])
        ]

        return _res[0]

    #        EC_all_merged_lst.append(EC_OHN_merged)
    #        EC_all_merged = pd.concat(EC_all_merged_lst)
    #        ORR_cath = EC_PorphSiO2.ORR_updated_pars(sweep_type_select='cathodic')
    #        ORR_an = EC_Pfrom collections import CounterorphSiO2.ORR_updated_pars(sweep_type_select='anodic')
    #        EC_OHN2 = pd.merge(template, pd.merge(ORR_an,pd.merge(HPRR, N2CV),on='SampleID'), on='SampleID')
    #        EC_OHN2_cath = pd.merge(template, pd.merge(ORR,pd.merge(HPRR, N2CV),on='SampleID'), on='SampleID')
    #        EC_OHN2.to_excel(FindExpFolder('PorphSiO2').compare.joinpath('EC_ORR_HPRR_N2.xlsx'))
    def export_to_xls(EC_OHN_merged):
        export_path = FindExpFolder("PorphSiO2").compare.joinpath(f"EC_pars_all.xlsx")
        if "Sweep_Type" in EC_OHN_merged.columns:
            with pd.ExcelWriter(export_path) as writer:
                for swp, swpgr in EC_OHN_merged.groupby("Sweep_Type"):
                    swpgr.to_excel(writer, sheet_name=swp)
                    swpgr.to_excel(export_path.with_name(f"EC_pars_{swp}.xlsx"))
        else:
            export_path = FindExpFolder("PorphSiO2").compare.joinpath(
                "EC_pars_no-sweep.xlsx"
            )
            EC_OHN_merged.to_excel(export_path)
        print(f"EC pars saved to:\n{export_path}")
        return export_path

    def edit_columns(func, template=pd.concat([PorphSiO2_template(), SampleCodes])):
        def wrapper(*args, **kwargs):
            if kwargs:
                pars_out, suffx = func(*args, **kwargs)
            else:
                pars_out, suffx = func(*args)

            _skipcols = set(
                EC_PorphSiO2.mcols
                + ["RPM_DAC_uni"]
                + list(PorphSiO2_template().columns)
                + list(EC_index.columns)
                + list(SampleCodes.columns)
            )
            cols = [
                i
                for i in pars_out.columns
                if i not in _skipcols and not i.startswith(f"{suffx}")
            ]
            pars_out = pars_out.rename(columns={i: f"{suffx}_" + i for i in cols})
            return pars_out

        return wrapper

    @edit_columns
    def HPRR(sweep_type_select=["anodic", "cathodic"]):
        hfs = []
        for swp in sweep_type_select:
            hprr_files = list(EC_PorphSiO2.folder.rglob(f"*{swp}*HPRR*disk*"))
            #            print(hprr_files)
            for hf in hprr_files:
                hprr_raw = pd.read_excel(hf)
                hprr_raw["file"] = hf.stem
                E_App_col = [i for i in hprr_raw.columns if "E_APP" in i.upper()][0]
                E_jmin = hprr_raw.iloc[np.abs(hprr_raw["jmAcm-2"]).idxmin()][E_App_col]
                sID = GetSampleID.try_find_sampleID(hf)[0]
                fit_lin_fit = linregress(hprr_raw[E_App_col], hprr_raw["HPRR_j0_Fit"])

                hfs.append(
                    {
                        "SampleID": sID,
                        "E_onset": E_jmin,
                        "dj/dE": fit_lin_fit[0],
                        "Sweep_Type": swp,
                    }
                )
        HPRR_pars_origin = pd.DataFrame(hfs)
        return HPRR_pars_origin, "HPRR"

    ##    @edit_columns
    #    def ORR_pars_all_old(self):
    #        LC_idx = self.index
    #        ORR_pars_all = Load_from_Indexes.ORR_pars_OVV(reload= False, extra_plotting=False, xls_out = False)
    #        ORR_needed = set(LC_idx.loc[LC_idx.PAR_exp == 'ORR'].PAR_file.unique())
    ##        ORR_pars_present = set(ORR_pars_all.loc[ORR_pars_all.PAR_file.isin(LC_idx.PAR_file.to_numpy())].PAR_file.unique())
    ##        ORR_pars_missing = ORR_needed.difference(ORR_pars_present)
    ##        orr_files,orrfs = list(EC_PorphSiO2.folder.rglob('*ORR_pars*')), []
    ##        for of in orr_files:
    ##            orr_raw = pd.read_excel(of,index_col=[0])
    ##            orr_raw.query('RPM > 1400')
    ##            orrfs.append(orr_raw.query('RPM > 1400'))
    ##        ORR_pars_origin = pd.concat(orrfs,ignore_index=True).reset_index(drop=True)
    #        return ORR_pars_all

    @edit_columns
    def ORR_pars(self):
        LC_idx = self.index
        ORR_pars = EC_PorphSiO2.ORR_pars_all.loc[
            (
                (EC_PorphSiO2.ORR_pars_all._type == "ORR_pars")
                & (
                    EC_PorphSiO2.ORR_pars_all.PAR_file.isin(
                        self.index.PAR_file.to_numpy()
                    )
                )
            )
        ]
        ORR_pars = ORR_pars.dropna(how="all", axis=1)

        return ORR_pars, "ORR"

    @edit_columns
    def KL_pars(self):
        LC_idx = self.index

        KL_pars = EC_PorphSiO2.ORR_pars_all.loc[
            (
                (EC_PorphSiO2.ORR_pars_all._type == "KL_pars")
                & (
                    EC_PorphSiO2.ORR_pars_all.PAR_file.isin(
                        self.index.PAR_file.to_numpy()
                    )
                )
            )
        ]
        KL_pars = KL_pars.dropna(how="all", axis=1)
        return KL_pars, "ORR"

    #    @edit_columns
    #    def ORR(self):
    #        LC_idx = self.index
    #        ORR_pars_all = Load_from_Indexes.ORR_pars_OVV(reload= False, extra_plotting=False, xls_out = False)
    #        ORR_needed = set(LC_idx.loc[LC_idx.PAR_exp == 'ORR'].PAR_file.unique())
    #        ORR_pars_present = set(ORR_pars_all.loc[ORR_pars_all.PAR_file.isin(LC_idx.PAR_file.to_numpy())].PAR_file.unique())
    #        ORR_pars_missing = ORR_needed.difference(ORR_pars_present)
    ##        orr_files,orrfs = list(EC_PorphSiO2.folder.rglob('*ORR_pars*')), []
    ##        for of in orr_files:
    ##            orr_raw = pd.read_excel(of,index_col=[0])
    ##            orr_raw.query('RPM > 1400')
    ##            orrfs.append(orr_raw.query('RPM > 1400'))
    ##        ORR_pars_origin = pd.concat(orrfs,ignore_index=True).reset_index(drop=True)
    #        return ORR_pars_all,'ORR'

    #    @edit_columns
    #    def ORR_updated_pars(sweep_type_select = ['anodic','cathodic']):
    #        orr_files,orrfs = list(EC_PorphSiO2.folder.rglob('PostEC*ORR*Jkin*pars*')), []
    #        for of in orr_files:
    #            orr_raw = pd.read_excel(of,index_col=[0])
    #            orr_raw.query('RPM > 1400')
    #            orrfs.append(orr_raw.query('RPM > 1400'))
    #        ORR_pars_origin = pd.concat(orrfs,ignore_index=True).reset_index(drop=True)
    #        sweep_col = [i for i in ORR_pars_origin.columns if 'SWEEP' in i.upper() and not 'RING'  in i.upper()][0]
    #        ORR_pars_origin_swp = ORR_pars_origin.loc[ORR_pars_origin[sweep_col].isin(sweep_type_select)]
    #        ORR_pars_origin_swp = ORR_pars_origin_swp.rename(columns={sweep_col : 'Sweep_Type'})
    ##        ORR_pars_swp  = {n : gr for n,gr in ORR_pars.groupby('Sweep_Type_disk')}
    ##        ORR_anod, ORR_cath = ORR_pars_swp.get('anodic'), ORR_pars_swp.get('cathodic')
    #        return ORR_pars_origin_swp,'ORR'

    @edit_columns
    def N2cv(
        self,
        sweep_type_select=["anodic", "cathodic"],
        unit="F",
        reload=False,
        use_daily=True,
        extra_plotting=False,
        xls_out=False,
    ):

        if not Pfolder.joinpath("N2_orig_data.pkl").exists() or reload == True:
            Cdl_pars_all = Load_from_Indexes.N2_pars_OVV(
                reload=reload,
                use_daily=use_daily,
                extra_plotting=extra_plotting,
                xls_out=xls_out,
            )

            Cdl_pars = Cdl_pars_all.loc[
                Cdl_pars_all.PAR_file.isin(self.index.PAR_file.to_numpy())
            ]
            #        IndexOVV_N2_pars_fn = FindExpFolder('VERSASTAT').PostDir.joinpath('N2Cdl_pars_IndexOVV_v{0}.pkl.compress'.format(FileOperations.version))
            Cdl_pars = Cdl_pars.assign(**{"E_RHE_mV": 1000 * Cdl_pars.E_RHE.to_numpy()})
            #            Cdl_pars.index = pd.MultiIndex.from_frame(Cdl_pars[['PAR_file','Sweep_Type_N2']])
            #        N2_files, N2fs = list(EC_PorphSiO2.folder.rglob('*CVs*xlsx')), []
            N2fs = []
            if unit == "mF":
                unit_factor = 1
            elif unit == "F":
                unit_factor = 1e-3
            else:
                unit_factor = 1
            for n2f, ngr in Cdl_pars.groupby("PAR_file"):
                idx_cols = [i for i in ngr.columns if ngr[i].nunique() == 1]
                _dc = [i for i in ngr.columns if ngr[i].nunique() > 1]
                #            sID = GetSampleID.try_find_sampleID(n2f)[0]
                ngr.index = pd.MultiIndex.from_frame(ngr[idx_cols])
                ngr.drop(columns=idx_cols, inplace=True)
                ngr = ngr.dropna(axis=1, how="all")
                for swp, swgrp in ngr.groupby("Sweep_Type_N2"):
                    if swp in sweep_type_select:
                        #                    anod = n2_raw.get(swp)
                        swgrp_Ev = swgrp.loc[
                            (swgrp.E_RHE_mV.isin(np.arange(0.0, 1000.0, 100)))
                            & (swgrp.Cdl_R > 0.8)
                        ]
                        _mgr = []
                        for n, gr in swgrp_Ev.groupby("E_RHE_mV"):
                            if len(gr) > 1:
                                _mean = pd.DataFrame(pd.DataFrame(gr.mean(axis=0)).T)
                                _mean.index = gr.take([0]).index
                                _mgr.append(_mean)
                            else:
                                _mgr.append(gr)
                        _swgr_Ev_mean = pd.concat(_mgr)
                        _pvt = _swgr_Ev_mean.pipe(
                            multiIndex_pivot,
                            index=None,
                            columns=["E_RHE_mV"],
                            values="Cdl",
                        )
                        _pvt = _pvt.assign(**{"Sweep_Type": swp})
                        N2fs.append(_pvt)
                    #                    evlst = {'SampleID' :  sID, 'Sweep_Type' : swp}
                    #                    for Ev in np.arange(0,1000,100):
                    #                        anod_05 = swgrp.loc[np.isclose(swgrp[EvRHE], Ev*1E-3,atol=0.001)]
                    #                        if not anod_05.empty:
                    #                            mean_anod05 = np.abs(anod_05.mean(axis=0)['Cdl'])
                    #                            if not np.isnan(mean_anod05):
                    #        #                        N2fs.append({'SampleID' :  sID, f'Cdl_mFcm-2_{Ev}' : mean_anod05, 'Sweep_Type' : swp})
                    #                                evlst.update({f'Cdl_{unit}cm-2_{Ev}' : mean_anod05*unit_factor})
                    #                    N2fs.append(pd.DataFrame(evlst,index=[(sID,swp)]))
                    else:
                        pass
            N2_orig = pd.concat([i.reset_index() for i in N2fs], ignore_index=True)
            N2_orig.columns = list(N2_orig.columns.get_level_values(0))
            #            N2_orig.index.names = N2fs[0].index.names
            N2_orig = N2_orig.rename(
                columns={
                    i: f"Cdl_{unit}cm-2_{int(i)}" for i in np.arange(0.0, 1000.0, 100)
                }
            )
            N2_orig = N2_orig.assign(**{"RPM_DAC": 0})
            N2_orig.to_pickle(Pfolder.joinpath("N2_orig_data.pkl"))
        else:
            N2_orig = pd.read_pickle(Pfolder.joinpath("N2_orig_data.pkl"))
        #        N2_orig = pd.DataFrame(N2fs) #.set_index('SampleID','Sweep_Type')
        return N2_orig, "N2"

    @edit_columns
    def EIS_pars(model_select=["Model(EEC_2CPE)"], _source="Load pars"):
        #        sweep_type_select = ['anodic','cathodic'], model_select = 'Model(Singh2015_RQRQR)'
        _source = "Load pars"
        if "files" in _source:
            eis_files, eisfs = (
                list(
                    EC_PorphSiO2.folder.parent.joinpath(
                        f"EIS_Porph_SiO2\{model_select}"
                    ).rglob("JOS*.xlsx")
                ),
                [],
            )
            if eis_files:
                for ef in eis_files:
                    eis_raw = pd.read_excel(ef, index_col=[0])
                    eisfs.append(eis_raw)
                EIS_pars_mod = pd.concat(eisfs, ignore_index=True).reset_index(
                    drop=True
                )
            else:
                print("EIS pars file list empty!!")
        else:
            EIS_pars_all = Load_from_Indexes.EIS_pars_OVV(
                reload=False, extra_plotting=False, xls_out=False
            )
            EIS_pars_mod = EIS_pars_all.loc[EIS_pars_all.Model_EEC.isin(model_select)]

        _sample_uniq_cols = set(
            [
                a
                for n, gr in EIS_pars_mod.groupby("SampleID")
                for a in [i for i in gr.columns if gr[i].nunique() == 1]
            ]
        )
        # Create EIS var columns with gas N2 or O2 as suffix names
        EPgrp = EIS_pars_mod.groupby("Gas")
        EP_N2, EP_O2 = EPgrp.get_group("N2").drop(columns="Gas"), EPgrp.get_group(
            "O2"
        ).drop(columns="Gas")

        EC_exp_index = [
            i for i in Load_from_Indexes.EC_label_cols if i not in ["PAR_file", "Gas"]
        ] + ["PAR_date_day"]

        _gasgrp = []
        for gas, grp in EPgrp:
            _varsgrp = [a for i in grp.lmfit_var_names.unique() for a in i.split(", ")]
            _varsgrp += ["Rct_kin" for i in _varsgrp if "Rct" in i] + [
                "Qad+Cdlp"
                for i in _varsgrp
                if all([i in _varsgrp for i in ["Qad", "Cdlp"]])
            ]
            #            grp.lmfit_var_names.unique()[0].split(', ')
            _grp = grp.rename(columns={i: i + f"_{gas}" for i in _varsgrp})
            _grp = _grp.drop(columns="Gas")
            #            _grp.set_index(EC_exp_index+[ i for i in list(_sample_uniq_cols) if i not in _varsgrp],inplace=True)
            #            _grp = _grp.drop(columns=EC_exp_index+[ i for i in list(_sample_uniq_cols) if i not in _varsgrp])
            #                    [i for i in Load_from_Indexes.EC_label_cols if i is not 'Gas']
            _gasgrp.append(_grp)
        #        _ggidx = [i.set_index(EC_exp_index) for i in _gasgrp]
        #        pd.concat(_ggidx,axis=0)
        #        _dups = [[(count,item) for item, count in collections.Counter(i.index.values).items() if count > 1] for i in _ggidx]
        #        _DUP_PARFILES = pd.concat(_ggidx).loc[[a[1] for i in _dups for a in i]].sort_values('PAR_date').PAR_file.unique()
        #        pd.merge(_gasgrp[0],_gasgrp[1], on =EC_exp_index+[ i for i in list(_sample_uniq_cols) if i not in _varsgrp])
        #        pd.merge(*_ggidx,left_index=True, right_index=True)
        EIS_N2O2 = pd.concat(_gasgrp, ignore_index=True)
        #        EIS_N2O2 = pd.merge(EP_N2,EP_O2, suffixes=['_N2','_O2'],on='SampleID')
        Rsis = [
            i
            for i in EIS_N2O2.columns
            if "Rs" in i and not any(c in i for c in ("stderr", "_kin_", "_setting"))
        ]
        Rct_cols = [
            i
            for i in EIS_N2O2.columns
            if "Rct" in i and not any(c in i for c in ("stderr", "_kin_"))
        ]
        #        EIS_pars_origin[Rsis] = EIS_pars_origin[Rsis].mask(EIS_pars_origin[Rsis] < 1)
        EIS_N2O2[Rsis] = EIS_N2O2[Rsis].mask(EIS_N2O2[Rsis] < 1)
        print("EIS Rs mask applied")
        EIS_N2O2[Rct_cols] = EIS_N2O2[Rct_cols].mask(EIS_N2O2[Rct_cols] > 1e5)
        print("EIS Rct mask applied")
        EIS_N2O2 = EIS_N2O2.dropna(axis=1, how="all")
        #        RedChiSq_limit = ORReis_merge.query('Rs > 1').RedChisqr.mean()+ 1*ORReis_merge.query('Rs > 1').RedChisqr.std()
        #        ORReis_neat = ORReis_merge.query('RedChisqr < @RedChiSq_limit & Rs > 2 & Rct < 9E05')
        EIS_N2O2_an, EIS_N2O2_cat = EIS_N2O2.copy(), EIS_N2O2.copy()
        EIS_N2O2_an["Sweep_Type"] = "anodic"
        EIS_N2O2_cat["Sweep_Type"] = "cathodic"
        EIS_N2O2_new = pd.concat([EIS_N2O2_an, EIS_N2O2_cat], axis=0)
        #        EIS_pars_orig_mod = EIS_pars_origin.query('Model_EEC == @model_select')
        return EIS_N2O2_new, "EIS"

    def EIS_spectra_origin_prep(model_select="Model(EEC_2CPE)"):
        eis_metaf, _specs = (
            list(
                EC_PorphSiO2.folder.parent.rglob(
                    "EIS_Porph_SiO2\meta_data*EIS*origin.xlsx"
                )
            ),
            [],
        )
        EISmeta = pd.read_excel(eis_metaf[0], index_col=[0])
        EISmeta.columns
        for (sID, gas), pgrp in EISmeta.groupby(["SampleID", "Gas"]):  # 'PAR_file'
            PF, pgrp
            EIScombined = pd.read_excel(pgrp.SpectraFile.iloc[0], index_col=[0])
            EISspectra_mod = EIScombined.query("Model_EEC == @model_select")
            EISspectra_mod = make_uniform_EvRHE(EISspectra_mod)
            for Ev, Egrp in EISspectra_mod.groupby("E_RHE"):
                Egrp = Egrp.assign(**{"SampleID": sID, "Gas": gas})
                _specs.append(Egrp)
        #                _specs.update({(sID,gas,Ev) : Egrp})
        spectra = pd.concat(_specs)
        spectra.to_excel(eis_metaf[0].with_name("clean_spectra.xlsx"))

    def EIS_spectra_origin(model_select="Model(EEC_2CPE)"):
        eis_metaf, _specs = (
            list(
                EC_PorphSiO2.folder.parent.rglob(
                    f"EIS_Porph_SiO2\{model_select}\meta_data*EIS*origin.xlsx"
                )
            ),
            [],
        )

        specdir = mkfolder(eis_metaf[0].parent.joinpath("spectra"))
        spectra = pd.read_excel(
            eis_metaf[0].with_name("clean_spectra.xlsx"), index_col=[0]
        )
        spectra.columns
        for ax_type in [("Zre", "-Zim"), ("Yre", "Yim")]:
            cols = [i + a for i in ["DATA_", "FIT_"] for a in ax_type]
            for gas, Ggrp in spectra.groupby("Gas"):
                for sID, sgrp in Ggrp.groupby("SampleID"):
                    with pd.ExcelWriter(
                        specdir.joinpath(f"{ax_type[0][0]}_{gas}_{sID}.xlsx")
                    ) as writer:
                        for Ev, Egrp in sgrp.groupby("E_RHE"):
                            # sheet_name = Ev
                            EmV = f"{1E3*Ev:.0f}"
                            Egrp[["Frequency(Hz)"] + cols].to_excel(
                                writer, sheet_name=EmV
                            )
                            # === plotting
                            fig, ax = plt.subplots()
                            Egrp.plot(
                                x=cols[0],
                                y=cols[1],
                                kind="scatter",
                                ax=ax,
                                label=cols[1],
                            )
                            Egrp.plot(x=cols[2], y=cols[3], c="r", ax=ax, label=cols[3])
                            plt.legend()
                            ax.set_xlabel(ax_type[0])
                            ax.set_ylabel(ax_type[1])
                            ax.set_title(f"{gas} {sID} {EmV}")
                            ax.grid(True)
                            plt.savefig(
                                specdir.joinpath(f"{ax_type[0][0]}_{gas}_{sID}_{EmV}"),
                                bbox_inches="tight",
                            )
                            plt.close()
                            # ===

    def save_load_AST_pars(func):
        #        argnames = func.func_code.co_varnames[:func.func_code.co_argcount]
        #        fname = func.func_name
        def wrapper(*args, **kwargs):
            func_args = inspect.signature(func).bind(*args, **kwargs).arguments
            func_args_str = ", ".join(
                "{} = {!r}".format(*item) for item in func_args.items()
            )
            print(f"{func.__module__}.{func.__qualname__} ( {func_args_str} )")
            #            args = list(args)
            #            print('must-have arguments are:')
            #            my_var_name = [ (k,v) for k,v in locals().items()]
            #            for i in my_var_name:
            #                print(f'{(i)}')
            ##            for i in args:
            ##                print(eval(i))
            #            print('optional arguments are:')
            #            for kw in kwargs.keys():
            #                print( kw+'='+str( kwargs[kw] ))
            return args

        return wrapper

    #    @save_load_AST_pars

    def corr_plots():
        EC_OHC.query('SampleID != "JOS5"').corr()
        corrstk = EC_OHC.query('SampleID != "JOS5"').corr().stack()

        EC_OHC.plot(x="E_onset", y="HPRR_E_onset", kind="scatter")
        EC_OHC.plot(x="FracH2O2_050", y="HPRR_E_onset", kind="scatter")
        EC_OHC.plot(x="N2_Cdl_mFcm-2_0.5", y="HPRR_dj/dE", kind="scatter")
        EC_OHC.plot(x="N2_Cdl_mFcm-2_0.5", y="E_half", kind="scatter")

        EC_OHC.corr(method="pearson")


# def EC_PorphSio():
##    folder = Path('F:\EKTS_CloudStation\CloudStation\Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc\EC_Porph_SiO2_0.1MH2SO4\Compare_parameters')
##    folder = Path('G:\CloudStation\Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc\EC_Porph_SiO2_0.1MH2SO4\Compare_parameters')
##    HPRR = pd.concat([pd.read_excel(i)['file'] for i in hprr_files])
#    EC_ORR_HPRR = pd.merge(ORR_pars_origin,HPRR_pars_origin)
#    HPRR_pars_origin.join(N2_orig, on='SampleID')
#    EC_OHC = pd.merge(ORR_pars_origin,pd.merge(HPRR_pars_origin, N2_orig),on='SampleID')
##        orr_raw.query('RPM > 1400')
##        orrfs.append(orr_raw.query('RPM > 1400'))
#    EC_OHC.to_excel(folder.joinpath('EC_ORR_HPRR.xlsx'))


class postChar:

    suffixes = ["R", "P", "EA", "B"]

    def __init__(self):
        self.folder = FindExpFolder("PorphSiO2").folder
        #        FindExpFolder().TopDir.joinpath(Path('Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc'))
        self.raman = postChar.postRaman
        self.bet = postChar.postBET
        self.ea = postChar.postEA
        self.prep = postChar.prec_ea_ratios(
            PorphSiO2_template(), self.ea, postChar.postPrep
        )
        self.merged = postChar.merged(self)

    def merged(self):
        cols = list(PorphSiO2_template().columns)
        merged_out = pd.merge(
            self.prep,
            pd.merge(self.ea, pd.merge(self.raman, self.bet, on=cols), on=cols),
            on=cols,
        )
        return merged_out

    def decorator(func):
        @functools.wraps(func)
        def wrapper_decorator(*args, **kwargs):
            # Do something before
            value = postChar.slice_out(*args, slice_lst=postChar.template().SampleID)
            # Do something after
            return value

        return wrapper_decorator

    def slice_out(func, template=PorphSiO2_template()):
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
    def postRaman(peak_model="5peaks", plot_resid=False):
        ramandir = FindExpFolder("PorphSiO2").folder / "SiO2_Me_EC+Struc" / "Raman"
        raman_files = ramandir.rglob("*xlsx")
        fitpars_fls = [
            i
            for i in raman_files
            if all([a in i.stem for a in ["FitParameters_Model", "peaks"]])
        ]
        FitPars_raw = pd.concat(
            [
                pd.read_excel(i).assign(**{"Model": i.stem.split("Model_")[-1]})
                for i in fitpars_fls
            ],
            sort=False,
        )
        FitPars = FitPars_raw.loc[
            FitPars_raw.SampleID.isin(PorphSiO2_template().SampleID.values)
        ]

        if plot_resid:
            Model_redchi = pd.DataFrame(
                [
                    (n, np.abs(gr.redchi).mean(), np.abs(gr.redchi).sum())
                    for n, gr in FitPars.groupby("Model")
                ],
                columns=["Model", "redchi_mean", "redchi_sum"],
            ).set_index("Model")
            Model_chi = pd.DataFrame(
                [
                    (n, gr.chisqr.mean(), gr.chisqr.sum())
                    for n, gr in FitPars.groupby("Model")
                ],
                columns=["Model", "redchi_mean", "redchi_sum"],
            ).set_index("Model")
            Model_redchi.plot.bar()
            Model_chi.plot.bar()
        if peak_model:
            FPars_out_1st = FitPars.loc[FitPars.Model.isin([peak_model])]
        else:
            FPars_out_1st = FitPars

        t2nd_mod = "2ndOrder_4peaks"
        if t2nd_mod in FitPars_raw.Model.unique():
            FPars_out_2nd = FitPars.loc[FitPars.Model == t2nd_mod].dropna(axis=1)
            flt2nd = get_float_cols(FPars_out_2nd)
            FPars_out_2nd = FPars_out_2nd.rename(
                columns=dict(zip(flt2nd, [f"2nd_{i}" for i in flt2nd]))
            )
            FPars_out = pd.merge(
                FPars_out_1st.dropna(axis=1),
                FPars_out_2nd,
                on="SampleID",
                how="left",
                suffixes=["_1st", "_2nd"],
            )

        else:
            FPars_out = FPars_out_1st
        return FPars_out, "R"

    @slice_out
    def postBET():
        betdir = FindExpFolder("PorphSiO2").folder.joinpath("SiO2_Me_EC+Struc/BET")
        bet_files = betdir.rglob("Overview_*xlsx")
        BET_ovv = pd.concat(
            [pd.read_excel(i, index_col=[0], sort=False) for i in bet_files]
        )
        return BET_ovv, "B"

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

    def prec_ea_ratios(templ, ea, prep):
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
        return prep


#        EA_results.loc[EA_results.SampleID.isin(isin_lst)]
class postCorrs:
    def __init__(self):
        pass

    def prep_ea(prep, ea, plots=False):
        prep_ea = pd.merge(prep, ea)
        corr_prep_ea = allchars_corr(prep_ea)
        if plots == True:
            testfolder = mkfolder(plotsfolder.joinpath("prep_ea_test"))
            for i in corr_prep_ea.index.values:
                plot(prep_ea, x=i[0], y=i[1], savefig=testfolder)
        return corr_prep_ea

    def get_EC_postchar():
        EC_postchars = EC_PorphSiO2.EC_merge_postchar()

        AllCorrs = {
            key: {"data": val, "corrs": postCorrs.allchars_corr(val, key)}
            for key, val in EC_postchars.items()
        }
        return AllCorrs

    #            Acorrs = postCorrs.allchars_corr(Allchars)

    def test_options():
        _postAST = ["postAST_LC", "postAST", "postAST_sHA"]
        _testopts = {}
        for _p in _postAST:
            _std = (f"{_p}", "cathodic", 1.0, 0.379)
            _tests = {
                "KL": {
                    "select": (*_std, "KL_I_Disk", 0.4),
                    "endswith": ["nElectrons_diff_perc"],
                },
                "ORR": {
                    "select": (*_std, 1500),
                    "endswith": ["Jkin_min_700_diff_perc"],
                },
                "N2CV": {"select": (_std), "endswith": ["Cdl_700_diff_perc"]},
                "EIS": {"select": (*_std, 0.6), "endswith": ["_diff_perc"]},
            }
            _testopts.update({_p: _tests})
        return _testopts

    def test_select_corrs(AllCorrs, swp="cathodic"):

        AllCorrs = postCorrs.get_EC_postchar()
        print("Experiment options:", ", ".join(AllCorrs.keys()))

        corvalnm = "corr_val"

        ECname = "EIS"
        ECgrp_keys = EC_types_grp().get(ECname)
        print("EC group format:", ECgrp_keys)
        ECgrp = ("postAST_sHA", "cathodic", 1.0, 0.379, 0.6)
        EC_select = dict(zip(ECgrp_keys, ECgrp))
        print(
            "Options:",
            "\n".join([str(i) for i in set(AllCorrs.get(ECname).get("corrs").keys())]),
        )
        chE = AllCorrs.get(ECname).get("corrs").get(ECgrp, "Empty dict")
        Allchars = AllCorrs.get(ECname).get("data")
        #        AllCan =  chE
        chE.chars.unique()
        charsECname = chE.loc[
            chE.chars.isin([i for i in chE.chars.unique() if ECname in i])
        ].sort_values(corvalnm)
        chars_struct = chE.loc[
            chE.chars.isin(
                [
                    i
                    for i in chE.chars.unique()
                    if any([c == ECname for c in i]) and not i == (ECname, ECname)
                ]
            )
        ].sort_values(corvalnm)
        #        chC = postCorrs.select(chE,postChar.suffixes,endswith_set = ['diff_abs'] ).sort_values(corvalnm)
        #        chC = postCorrs.select(chE,postChar.suffixes,endswith_set = ['C_content'] ).sort_values(corvalnm)
        chC = postCorrs.select(
            chE, postChar.suffixes + ["Rct"], endswith_set=[""]
        ).sort_values(corvalnm)
        set([i[0] for i in chC.index.unique()])
        postCorrs.combo_makers(
            chC,
            ECname,
            "",
            Allchars,
            EC_select,
            separate_dir="_".join([str(i) for i in (ECname, *ECgrp)]),
            corrval_zero=1,
        )
        comboCH = {}
        for Charnm in ["P", "EA", "B", "R"]:
            _out = postCorrs.combo_makers(
                chC,
                ECname,
                Charnm,
                Allchars,
                EC_select,
                include="Rct",
                separate_dir="_".join([str(i) for i in (ECname, *ECgrp)]),
                corrval_zero=1,
                force_replot=False,
                endswith="diff_perc",
            )
            _topidx = _out.query("corr_val > 0.75 or corr_val < 0.75")
            _counter = Counter([a for i in _topidx.index for a in i]).most_common(30)
            _ordering = [
                (
                    *_c,
                    _out.loc[[i for i in _topidx.index if any([c in _c[0] for c in i])]]
                    .corr_val.abs()
                    .sum()
                    .round(3),
                )
                for _c in _counter
            ]
            _ordering = sorted(_ordering, key=lambda x: x[-1], reverse=True)
            comboCH.update(
                {Charnm: {"combos": _out, "counter": _counter, "ordered": _ordering}}
            )

        [(key, val["ordered"][0]) for key, val in comboCH.items()]

        _topidx = comboCH["R"].query("corr_val > 0.8 or corr_val < 0.8").index

    #        chC.loc[chC.chars != ('EIS','EIS')]
    #        select(chE,endswith_set = 'diff_abs' )

    def test_combos():
        pass

    #    AllCan, main_ch,second_ch, Allchars,endswith, separate_dir, corrval_zero  = chE,'N2','',Allchars,['300_diff_abs'],'N2_test',1
    #        ['diff_perc', 'diff_abs']
    def combo_makers(
        AllCan,
        main_ch,
        second_ch,
        Allchars,
        EC_select,
        exclude_second="",
        include="",
        endswith="",
        separate_dir="",
        corrval_zero=1,
        cutoff_set=0,
        sep_topdir="",
        logy_set=False,
        line=False,
        force_replot=True,
        max_plots=100,
    ):
        chEA = postCorrs.uniq(AllCan, (main_ch, second_ch), ECgrp)
        chEA_excl = postCorrs.select(chEA, second_ch, endswith_set=endswith)

        if type(exclude_second) == type(""):
            exclude_second = [exclude_second]
        if type(exclude_second) == type([]):
            if exclude_second:
                for excl in exclude_second:
                    chEA_excl = postCorrs.select(chEA_excl, second_ch, exclude=excl)

        if type(include) == type(""):
            include = [include]
        if type(include) == type([]):
            for incl in include:
                chEA_excl = postCorrs.select(chEA_excl, incl)
        #                print(len(chEA_excl))
        #        chBET_SA = select(chEA, second_ch, exclude = exclude_second)
        #        chBET_SA = select(chBET_SA , '',exclude = 'Volume')
        cutoff = 0.4 if len(chEA_excl) > 300 else 0.1
        if cutoff:
            cutoff = cutoff_set
            corr_cutoff = chEA_excl[chEA_excl[corvalnm].abs() > cutoff]
            _ct = 0
            while len(corr_cutoff) < 30 and (cutoff + _ct) > 0:
                corr_cutoff = chEA_excl[chEA_excl[corvalnm].abs() > cutoff + _ct]
                _ct += -0.1
        if len(corr_cutoff) > max_plots:
            corr_cutoff = pd.concat(
                [chEA_excl.head(int(max_plots / 2)), chEA_excl.tail(int(max_plots / 2))]
            )
        print(
            f"Starting to plot ({main_ch}, {second_ch};+{include}+ ends({endswith}), -{exclude_second}): {len(corr_cutoff)}"
        )

        #        _ECgrps = EC_types_grp().get([i for i in EC_types_grp().keys() if main_ch in i][0])
        _Allchars_grp = Allchars.groupby(list(EC_select.keys())).get_group(
            tuple(EC_select.values())
        )

        for n, r in corr_cutoff.iterrows():
            combo_dest_dir = mkfolder(
                plotsfolder.joinpath(f"corrs_{r.chars[0]}_{r.chars[1]}")
            )
            if separate_dir:
                combo_dest_dir = mkfolder(combo_dest_dir.joinpath(separate_dir))
            if sep_topdir:
                combo_dest_dir = mkfolder(plotsfolder.joinpath(f"corrs_{sep_topdir}"))
            #            if corrval_zero == 0:
            _corr_val_set = 0 if corrval_zero != 1 else r[corvalnm]
            x, y = n[0], n[1]
            try:
                plot(
                    _Allchars_grp,
                    x=x,
                    y=y,
                    corr_val=_corr_val_set,
                    extra_title=EC_select,
                    savefig=combo_dest_dir,
                    logy=logy_set,
                    add_line=line,
                    force_replot=force_replot,
                )
                plt.clf()
                plt.close("all")
            except Exception as e:
                print(f"Plotting fail for x {x} and y {y},\n because: {e}")
                print(f"x: {_Allchars_grp[x].values},y: {_Allchars_grp[y].values}")
        return chEA_excl

    def all_combos_chars_EC(EC_name):
        _dict = {}
        for Charnm in ["P", "EA", "B", "R"]:
            N2P = combo_makers(AllCan, EC_name, Charnm, Allchars, exclude_second="")
            _dict.update({f"{EC_name}_{Charnm}": N2P})
        return _dict

    def allchars_corr(Allchars, ECname):
        AllCorrs = {}
        _ECgrps = EC_types_grp().get(ECname)
        #        corvalnm = f'corr_val_{swp}'
        corvalnm = "corr_val"

        if "Sweep_Type" not in Allchars.columns:
            Allchars["Sweep_Type"] = "none"

        for ngrp, swgrp in Allchars.groupby(_ECgrps):

            swp = [i for i in ngrp if i in ["cathodic", "anodic"]][0]

            only_floats = [
                key
                for key, val in swgrp.dtypes.to_dict().items()
                if "float64" in str(val)
            ]
            filter_cols_lst = [
                "error",
                "chi",
                "Coefficient",
                "rvalue",
                "pvalue",
                "Isotherm_P/P0",
                "relPrange",
                "wrong",
                "B_DFT",
                "AMBIENT ",
                "R_bic",
                "R_aic",
                "BET_P/Po_",
                "BET_1/(W((Po/P)-1))",
                "B_Area_in_cell_m2",
                "Delta_PAR_LastMod",
                "pH",
                "filesize",
                "Loading",
                "Date_PAR_EXP",
                "dupli_num",
                "index",
                "RHE_fn",
                "delta_mtime",
                "Instrument_SN",
                "Scanrate",
                "fitR",
            ]
            filter_cols_lst += [
                i for i in only_floats if i.startswith("B_") and "_RPT" in i
            ]
            filter_cols_lst += [
                i
                for i in only_floats
                if i.startswith("B_")
                and any(
                    i.endswith(c) for c in ["_microArea", "_Vmicro", "_r", "mass_RAW"]
                )
            ]
            filter_cols_lst += [
                i
                for i in only_floats
                if i.startswith("B_t") and any(i.endswith(c) for c in ["_slope"])
            ]
            filter_cols_lst += [
                i
                for i in only_floats
                if i.startswith("R_")
                and any(i.endswith(c) for c in ["_amplitude", "_sigma", "_bic", "_aic"])
            ]

            _EIS_filters_list = [
                "Chisqr",
                "stderr",
                "RedChisqr",
                "AppV",
                "aic",
                "bic",
                "trimming",
                "test_errcheck_msg",
                "lmfit_bic",
                "Scanrate",
                "deltaDT",
            ]
            filter_cols_lst += [
                i
                for i in only_floats
                if i.startswith("EIS_")
                and any(c in "_".join(i.split("EIS_")[1:]) for c in _EIS_filters_list)
            ]
            #            += [i for i in only_floats if i.endswith('diff_abs') and i.endswith('diff_perc')]

            #        filter_cols_lst += [i for i in only_floats if i.startswith('B_') and i.endswith()]
            #        filter_cols_lst += [i for i in only_floats if i.startswith('B_') and i.endswith()]
            filtred_cols = [
                i for i in only_floats if not any(t in i for t in filter_cols_lst)
            ]
            allcorrstk = swgrp[filtred_cols].corr(method="pearson").stack()
            allcorrs = allcorrstk[np.abs(allcorrstk) < 1]
            allcorrs.name = corvalnm
            allchrrs = pd.concat(
                [
                    allcorrs,
                    pd.Series(
                        [
                            (i[0].split("_")[0], i[1].split("_")[0])
                            for i in allcorrs.index
                        ],
                        name="chars",
                        index=allcorrs.index,
                    ),
                ],
                axis=1,
            )
            allchrrs = (
                allchrrs[np.abs(allchrrs[allcorrs.name]) < 0.991]
                .drop_duplicates(allcorrs.name)
                .sort_values(allcorrs.name)
            )

            allchflrtd = postCorrs.remove_corrdups(allchrrs, ECname)
            allchflrtd = allchflrtd.assign(**dict(zip(_ECgrps, ngrp)))
            AllCorrs.update({ngrp: allchflrtd})

        if "none" in Allchars["Sweep_Type"].unique()[0]:
            return AllCorrs.get("none").sort_values(by="corr_val_none")
        else:
            return AllCorrs

    def remove_corrdups(allchrrs, ECname):
        sl = allchrrs.loc[allchrrs.chars == (ECname, ECname)]
        sl_non_idx = allchrrs.loc[~allchrrs.index.isin(sl.index)]
        non_dupl = [
            i
            for i in sl.index
            if len(
                set(["_".join(i[0].split("_")[1:-1]), "_".join(i[1].split("_")[1:-1])])
            )
            > 1
        ]
        rings = ["J_ring", "Frac_H2O2", "n_ORR"]
        ring_dups = [
            i
            for i in non_dupl
            if not (any(r in i[0] for r in rings) and any(r in i[1] for r in rings))
        ]
        jkin_dups = [
            i for i in ring_dups if not (("Jkin_" in i[0] and "Jkin_" in i[1]))
        ]

        _incl_allchrrs = allchrrs.loc[list(sl_non_idx.index) + jkin_dups].sort_index()
        _dupl_diff_cols = [
            i
            for i in allchrrs.index
            if (
                (
                    i[0].endswith("_diff_abs")
                    and i[1].endswith("diff_perc")
                    or (i[1].endswith("_diff_abs") and i[0].endswith("diff_perc"))
                )
            )
        ]
        _excl_allchrrs = _incl_allchrrs.loc[~_incl_allchrrs.index.isin(_dupl_diff_cols)]
        return _excl_allchrrs

    def uniq(AllCan, combo_select, self_corr=False, Ev: int = 0):
        corvalnm = "corr_val"
        uniqchars = AllCan["chars"].unique()

        if any(i in combo_select for i in ["*", ""]):
            uniq_char = [i for i in uniqchars if any(c in i for c in combo_select)]
        else:
            uniq_char = [i for i in uniqchars if set(i) == set(combo_select)]
        if self_corr == True:
            uniq_char += [
                i for i in uniqchars if (combo_select[0], combo_select[0]) == i
            ]
        chC = AllCan.loc[AllCan["chars"].isin(uniq_char)].sort_values(corvalnm)
        if Ev > 50:
            AllC_Ev = get_corrs_E(chC, Ev).sort_values(corvalnm)
            return AllC_Ev
        else:
            return chC

    def select(chC, one_col: list = "", exclude="", endswith_set=""):
        idx_select = chC.index.values

        if one_col:
            if not any(a in i for i in ["", "*"] for a in one_col):
                idx_select = [
                    i for i in idx_select if any(col in a for a in i for col in one_col)
                ]
        if endswith_set:
            _endx = []
            if type(endswith_set) == type(""):
                _endx += [
                    i for i in idx_select if any(a.endswith(endswith_set) for a in i)
                ]
            elif type(endswith_set) == type([]):
                for _ends in endswith_set:
                    _endx += [
                        i for i in idx_select if any(a.endswith(_ends) for a in i)
                    ]
            idx_select = [i for i in idx_select if i in _endx]
        if exclude:
            idx_select = [i for i in idx_select if not any(exclude in a for a in i)]

        return chC.loc[idx_select].sort_values(
            by=[i for i in chC.columns if "corr_val" in i][0]
        )

    #            return chC

    def take_same_Ev(AllCan, allowed_diff=200):
        _take, _not = [], []
        for n, r in AllCan.iterrows():
            _num_idx0, _num_idx1 = -1, -1
            if n[0].startswith("EIS"):
                _num_idx0 = -2
            if n[1].startswith("EIS"):
                _num_idx1 = -2
            _n0, _n1 = n[0].split("_")[_num_idx0], n[1].split("_")[_num_idx1]

            if all([_n0.isnumeric(), _n1.isnumeric()]):
                if np.abs(float(_n0) - float(_n1)) < allowed_diff:
                    _take.append(n)
                else:
                    _not.append(n)
            else:
                _take.append(n)
        print(f"Rows left out: {len(AllCan) - len(_take)}")
        return AllCan.loc[_take]


def collect_pars():
    ECpars = EC_PorphSiO2().pars
    pst = postChar()
    raman, bet, ea, prep = pst.raman, pst.bet, pst.ea, pst.prep
    struc = pst.merged
    Allchars = pd.merge(ECpars, struc)

    char = {"EC": ECpars, "raman": raman, "bet": bet, "ea": ea, "prep": prep}
    return ECpars, raman, bet, ea, prep, struc, char, Allchars


def export_chars(char, suffix=""):
    for nm, df in char.items():
        swpcol = [i for i in df.columns if "Sweep" in i]
        if swpcol:
            for swp, swgr in df.groupby(swpcol[0]):
                swgr.to_excel(plotsfolder.joinpath(f"char_{nm}_{swp}{suffix}.xlsx"))
        else:
            df.to_excel(plotsfolder.joinpath(f"char_{nm}{suffix}.xlsx"))


def make_corrs(pars, charname):
    lst = []
    #    ECpars[[key for key,val in ECpars.dtypes.to_dict().items() if 'float64' in str(val)]]
    only_floats = [
        key for key, val in pars.dtypes.to_dict().items() if "float64" in str(val)
    ]
    for n, r in pars.iterrows():
        without1 = pars.drop(n)[only_floats]
        wcorrstk = without1.corr(method="pearson").stack()
        lst.append(wcorrstk.loc[wcorrstk < 1].rename(f"{charname}_{r.SampleID}"))
    corrs_outlier = pd.concat(lst, axis=1)
    abs_outliers = (
        np.abs(corrs_outlier).describe().sort_values("mean", axis=1, ascending=False)
    )
    abs_print = ",".join(
        [
            f"{key} : {val:.2f}"
            for key, val in (abs_outliers.loc["mean"]).to_dict().items()
        ]
    )
    print("Correlation without samples in collumns, indicates outliers...\n", abs_print)
    return {"corrs_outliers": corrs_outlier, "describe": abs_outliers}


def corrloop(char):
    out = []
    for charname, pars in char.items():
        out.append(make_corrs(pars, charname))
    return out


def corrplot():
    corr = PlotAnalysis.corrplot(ECpars.corr())
    corr = corr.loc[
        ~(corr.y.isin(template.columns) | corr.x.isin(template.columns))
        & (corr.value < 1)
    ]
    tops = corr.loc[corr.y.str.contains("HPRR")].sort_values(
        by="value", ascending=False
    )
    for n, r in tops.tail(10).iterrows():
        plot(ECpars, x=r.x, y=r.y)


def plot_heatmap():
    fig, ax = plt.subplots(figsize=(20, 20))
    ax = sns.heatmap(allchflrtd[allcorrs.name].unstack())
    plt.savefig(FindExpFolder("PorphSiO2").compare.joinpath("heatmap.png"))
    plt.close()


# ========= TDOD Main testing of correlations STARTS HERE --====---- ===========
def sort_corrs(AllCorrs):
    swp = "anodic"
    corvalnm = f"corr_val_{swp}"
    AllCan = AllCorrs.get(swp)
    AllChwoJOS5 = Allchars.query('SampleID != "JOS5"')
    AllCanwoJOS5 = allchars_corr(AllChwoJOS5).get(swp)
    Allcan_EV = take_same_Ev(AllCan, allowed_diff=200)
    AllCanwoJOS5_EV = take_same_Ev(AllCanwoJOS5, allowed_diff=200)

    Allcan_EV_dE50 = take_same_Ev(AllCan, allowed_diff=50)
    AllCanwoJOS5_EV_dE50 = take_same_Ev(AllCanwoJOS5, allowed_diff=50)

    chEA = uniq(AllCan, "EA", corvalnm)
    chBET_SA = select(AllCan, "P_", exclude="")
    chBET_SA = select(chEA, "P_", exclude="P_MW")
    char_cols = ["P", "EA", "B", "R"]
    #    for n,r in chBET_SA.head(10).iterrows():
    #        plot(Allchars,x=n[0],y=n[1])
    #    for elem in ['C','H','N']:
    #        plot(Allchars,x=f'P_prec_{elem}',y=f'EA_{elem}_content',equal=True,savefig=mkfolder(plotsfolder.joinpath('prep_ea_out')))
    #    plot(Allchars,x=f'P_prec_100-CNH',y=f'EA_100-CHN',equal=True,savefig=mkfolder(plotsfolder.joinpath('prep_ea_out')))
    #   main_ch,second_ch ='N2', 'B',
    #   main_ch,second_ch ='ORR', '*',
    #    AllCan, main_ch,second_ch, Allchars = Allcan_EV_dE50,'*', 'EIS', Allchars #, include = char_cols + eisplot.parlst,exclude_second = ['linKK'],sep_topdir=f'EIS_char_best',corrval_zero = 1,cutoff_set=0.6
    def combo_makers(
        AllCan,
        main_ch,
        second_ch,
        Allchars,
        exclude_second="",
        include="",
        endswith="",
        separate_dir="",
        corrval_zero=1,
        cutoff_set=0,
        sep_topdir="",
        logy_set=False,
        line=False,
        force_replot=True,
        max_plots=100,
    ):
        chEA = uniq(AllCan, (main_ch, second_ch), corvalnm)
        chEA_excl = select(chEA, second_ch, endswith_set=endswith)

        if type(exclude_second) == type(""):
            exclude_second = [exclude_second]
        if type(exclude_second) == type([]):
            for excl in exclude_second:
                chEA_excl = select(chEA_excl, second_ch, exclude=excl)

        if type(include) == type(""):
            include = [include]
        if type(include) == type([]):
            for incl in include:
                chEA_excl = select(chEA_excl, incl)
        #                print(len(chEA_excl))
        #        chBET_SA = select(chEA, second_ch, exclude = exclude_second)
        #        chBET_SA = select(chBET_SA , '',exclude = 'Volume')
        cutoff = 0.4 if len(chEA_excl) > 300 else 0.1
        if cutoff_set:
            cutoff = cutoff_set
        corr_cutoff = chEA_excl[chEA_excl[corvalnm].abs() > cutoff]
        if len(corr_cutoff) > max_plots:
            corr_cutoff = pd.concat(
                [chEA_excl.head(int(max_plots / 2)), chEA_excl.tail(int(max_plots / 2))]
            )
        print(
            f"Starting to plot ({main_ch}, {second_ch};+{include}+ ends({endswith}), -{exclude_second}): {len(corr_cutoff)}"
        )
        for n, r in corr_cutoff.iterrows():

            combo_dest_dir = mkfolder(
                plotsfolder.joinpath(f"corrs_{r.chars[0]}_{r.chars[1]}")
            )
            if separate_dir:
                combo_dest_dir = mkfolder(combo_dest_dir.joinpath(separate_dir))
            if sep_topdir:
                combo_dest_dir = mkfolder(plotsfolder.joinpath(f"corrs_{sep_topdir}"))
            #            if corrval_zero == 0:
            _corr_val_set = 0 if corrval_zero != 1 else r[corvalnm]
            plot(
                Allchars,
                x=n[0],
                y=n[1],
                corr_val=_corr_val_set,
                savefig=combo_dest_dir,
                logy=logy_set,
                add_line=line,
                force_replot=force_replot,
            )
            plt.clf()
            plt.close("all")
        return chEA_excl

    def all_combos_chars_EC(EC_name):
        _dict = {}
        for Charnm in ["P", "EA", "B", "R"]:
            N2P = combo_makers(AllCan, EC_name, Charnm, Allchars, exclude_second="")
            _dict.update({f"{EC_name}_{Charnm}": N2P})
        return _dict

    P_EA = combo_makers(AllCan, "P", "EA", Allchars, corrval_zero=0)
    P_P = combo_makers(AllCan, "P", "P", Allchars)

    R_EA = combo_makers(AllCan, "R", "EA", Allchars, exclude_second="2nd")
    R_P = combo_makers(
        AllCanwoJOS5,
        "P",
        "R",
        AllChwoJOS5,
        exclude_second="2nd",
        include="Metal",
        separate_dir="Metal2",
    )

    R_EA2 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="2nd",
        include="N_content",
        separate_dir="N_content",
    )

    R_EA4 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="D4D4",
        include="D4",
        separate_dir="D4_peak",
    )
    R_EA3 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="",
        include="D3",
        separate_dir="D3_peak",
    )
    R_EA3 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="",
        include="G",
        separate_dir="G_peak",
    )

    N2_combos = all_combos_chars_EC("N2")
    N2_B_Extra = combo_makers(
        AllCan,
        "N2",
        "B",
        Allchars,
        include=["N2_Cdl", "B_BET_RAW_calc_SA"],
        separate_dir="Cdl_extra",
        corrval_zero=0,
    )
    N2_B_Extra = combo_makers(
        AllCan,
        "N2",
        "B",
        Allchars,
        include=["N2_Cdl", "B_BET_RAW_calc_constant"],
        separate_dir="Cdl_extra_C",
        corrval_zero=0,
        logy_set=True,
    )

    ### EIS correlations
    #    EIS_combos = all_combos_chars_EC('EIS')
    N2_EIS_dE50 = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_nAd"],
        separate_dir="N2_eis_nAd_dE50",
        corrval_zero=0,
    )
    N2_EIS_dE50 = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_nDL"],
        separate_dir="N2_eis_nDL_dE50",
        corrval_zero=0,
    )

    N2_EIS_dE50 = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "Qad+Cdlp"],
        separate_dir="N2_eis_Qad+Cdlp_dE50",
        corrval_zero=0,
    )
    N2_EIS_dE50_Qad = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_Qad_"],
        separate_dir="N2_eis_Qad_dE50",
        corrval_zero=0,
    )
    N2_EIS_dE50_Qad = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_Cdlp_"],
        separate_dir="N2_eis_Cdl_dE50",
        corrval_zero=0,
    )

    best_EIS_char = combo_makers(
        Allcan_EV_dE50,
        "*",
        "EIS",
        Allchars,
        include="",
        exclude_second=["linKK"],
        sep_topdir=f"EIS_char_best",
        corrval_zero=1,
        cutoff_set=0.9,
    )

    def make_all_EIS_combos():
        for par in eisplot.parlst[:3]:
            not_par = [i for i in eisplot.parlst if i not in par]
            #        for C in ['Qad_','_Cdlp_','Qad+Cdlp']:
            plt.close("all")
            pdir = f"EIS_pars/{par}"
            pexcl = ["linKK"] + not_par

            for gas in ["O2", "N2"]:
                for orrC in [
                    "J_diff_lim",
                    "_Jkin_min_600",
                    "H2O2_500",
                    "n_ORR",
                    "J_ring",
                ][:]:
                    ORR_EIS_Extra = combo_makers(
                        AllCanwoJOS5_EV_dE50,
                        "ORR",
                        "EIS",
                        Allchars,
                        include=[orrC, par],
                        exclude_second=pexcl,
                        sep_topdir=f"{pdir}/{orrC}/{gas}",
                        corrval_zero=1,
                    )
                #            plt.close('all')
                #                N2P = combo_makers(Allcan_EV_dE50,'B', 'EIS', Allchars, include = ['B', par, 'B_BET_RAW_calc_SA m2/g'],exclude_second = pexcl, endswith = gas, sep_topdir=f'{pdir}/{Charnm}_BET/{gas}',corrval_zero = 1)
                #                N2P = combo_makers(Allcan_EV_dE50,'EA', 'EIS', Allchars, include = ['EA', par, 'EA_N_content'],exclude_second = pexcl, endswith = gas, sep_topdir=f'{pdir}/{Charnm}_EA_N_content/{gas}',corrval_zero = 1)
                #                N2P = combo_makers(Allcan_EV_dE50,'EA', 'EIS', Allchars, include = ['EA', par, 'EA_C_content'],exclude_second = pexcl, endswith = gas, sep_topdir=f'{pdir}/{Charnm}_EA_C_content/{gas}',corrval_zero = 1)
                N2P = combo_makers(
                    Allcan_EV_dE50,
                    "EA",
                    "EIS",
                    Allchars,
                    include=["EA", par, "EA_H_content"],
                    exclude_second=pexcl,
                    endswith=gas,
                    sep_topdir=f"{pdir}/{Charnm}_EA_H_content/{gas}",
                    corrval_zero=1,
                )
                for Charnm in char_cols:
                    N2P = combo_makers(
                        Allcan_EV_dE50,
                        Charnm,
                        "EIS",
                        Allchars,
                        include=[Charnm, par],
                        exclude_second=pexcl,
                        endswith=gas,
                        sep_topdir=f"{pdir}/{Charnm}/{gas}",
                        corrval_zero=1,
                    )
            plt.close("all")

            #                plt.close()
            #                    AllCan,EC_name, Charnm, Allchars, exclude_second = '')

            N2_EIS_dE50_Qad = combo_makers(
                Allcan_EV_dE50,
                "N2",
                "EIS",
                Allchars,
                include=["N2_Cdl", par],
                exclude_second=pexcl,
                sep_topdir=f"{pdir}",
                corrval_zero=0,
            )
            #            N2_EIS_dE50_Qad = combo_makers(Allcan_EV_dE50,'N2', 'EIS', Allchars, include = ['N2_Cdl','EIS_Qad_'],separate_dir='N2_eis_Qad_dE50',corrval_zero = 0)

            for HPRRc in ["HPRR_dj/dE", "HPRR_E_onset"]:
                HPRR_EIS_Extra = combo_makers(
                    AllCanwoJOS5_EV,
                    "HPRR",
                    "EIS",
                    Allchars,
                    include=[HPRRc, par],
                    exclude_second=pexcl,
                    sep_topdir=f'{pdir}/{HPRRc.replace("/","_")}',
                    corrval_zero=1,
                )
            #                plt.close()
            for orrC in ["J_diff_lim", "_Jkin_min_600", "H2O2_500", "n_ORR", "J_ring"][
                :
            ]:
                ORR_EIS_Extra = combo_makers(
                    AllCanwoJOS5_EV_dE50,
                    "ORR",
                    "EIS",
                    Allchars,
                    include=[orrC, par],
                    exclude_second=pexcl,
                    sep_topdir=f"{pdir}/{orrC}",
                    corrval_zero=1,
                )
            plt.close("all")

    make_all_EIS_combos()

    HPRR_Extra = combo_makers(
        AllCan,
        "HPRR",
        "ORR",
        Allchars,
        include=["HPRR_dj/dE"],
        separate_dir="djdE",
        corrval_zero=0,
    )
    HPRR_Extra = combo_makers(
        AllCan,
        "HPRR",
        "ORR",
        Allchars,
        include=["HPRR_dj/dE", "J_ring"],
        separate_dir="djdE",
        corrval_zero=0,
    )
    HPRR_orr = combo_makers(
        AllCan, "HPRR", "ORR", Allchars
    )  # , include = ['HPRR_dj/dE'],separate_dir='djdE',corrval_zero = 0)

    TS_corr = combo_makers(
        AllCan, "ORR", "EA", Allchars, include=["ORR_TSa_min"], sep_topdir="TSa"
    )
    TS_corr = combo_makers(
        AllCan, "ORR", "*", Allchars, include=["ORR_TSa_max"], sep_topdir="TSa_max"
    )

    TS_corr = combo_makers(
        AllCanwoJOS5, "ORR", "EA", Allchars, include=["ORR_TSa"], sep_topdir="TSa_EA"
    )

    HPRR_combos = all_combos_chars_EC("HPRR")

    HPRR_P = combo_makers(AllCan, "N2", "P", Allchars, exclude_second="")
    HPRR_EA = combo_makers(AllCan, "N2", "EA", Allchars, exclude_second="")
    HPRR_B = combo_makers(
        AllCan, "N2", "B_BET_RAW_calc_SA", Allchars, exclude_second=""
    )
    HPRR_R = combo_makers(AllCan, "N2", "R", Allchars, exclude_second="")

    HPPR_EIS = combo_makers(
        AllCanwoJOS5, "HPRR", "EIS", Allchars, exclude_second="linKK"
    )

    N2_EIS = combo_makers(Allcan_EV, "N2", "EIS", Allchars, exclude_second="linKK")

    ORR_EIS = combo_makers(
        AllCanwoJOS5_EV, "ORR", "EIS", Allchars, exclude_second="linKK"
    )

    plot(
        Allchars,
        x="EA_H_content",
        y="B_BET_RAW_calc_constant",
        equal=True,
        logy=True,
        savefig=mkfolder(plotsfolder.joinpath("prep_ea_out")),
    )

    def plot_only(
        Chars, Corrs, selection, number_plots=100, exclusion="", chars_excl=""
    ):
        if chars_excl:
            Corrs = Corrs.loc[Corrs.chars != ("R", "R")]
        chD2 = select(Corrs, selection, exclude=exclusion)
        for n, r in pd.concat([chD2.head(50), chD2.tail(50)]).iterrows():
            plot(
                Chars,
                x=n[0],
                y=n[1],
                corr_val=r[corvalnm],
                savefig=mkfolder(plotsfolder.joinpath(f"corrs_only_{selection}")),
            )

    plot_only(
        Allchars.query('SampleID != "JOS5"'),
        select(AllCanwoJOS5, "C/N_ratio"),
        "Metal_element",
    )
    plot_only(Allchars, select(AllCanwoJOS5, "C/N_ratio"), "C/N_ratio")

    plot_only(Allchars, AllCan, "2nd_", chars_excl=("R", "R"))

    def D2_only():
        chD2 = select(AllCan, "R_D2_height", exclude="_sigma")
        for n, r in pd.concat([chD2.head(50), chD2.tail(50)]).iterrows():
            plot(
                Allchars,
                x=n[0],
                y=n[1],
                corr_val=r[corvalnm],
                savefig=mkfolder(plotsfolder.joinpath(f"corrs_R_D2height_only")),
            )

    def HPRR_only():
        chD2 = select(AllCanwoJOS5, "HPRR_E_onset", exclude="_sigma")
        for n, r in chD2.iterrows():
            plot(
                Allchars,
                x=n[0],
                y=n[1],
                corr_val=r[corvalnm],
                savefig=mkfolder(plotsfolder.joinpath(f"corrs_HPRR_onset_only")),
            )

    for n, r in AllCan.head(10).iterrows():
        plot(swgrp, x=n[0], y=n[1])

    for ch in uniqchars:
        chc = AllCan.query("chars == @ch").sort_values(corvalnm)
        for n, r in chc.tail(1).iterrows():
            plot(swgrp, x=n[0], y=n[1])

    # == RAMAN corrs ==
    #    uniqR = [i for i in uniqchars if 'R' in i and i != ('R','R')]
    #    chR = AllCan.loc[AllCan['chars'].isin(uniqR)].sort_values(corvalnm)
    chR2 = uniq(AllCan, "R", corvalnm)
    for n, r in chR.tail(5).iterrows():
        plot(swgrp, x=n[0], y=n[1])
    for n, r in chR.head(5).iterrows():
        plot(swgrp, x=n[0], y=n[1])

    AllC_500 = get_corrs_E(AllCan, 600).sort_values(corvalnm)
    uniqEA = [i for i in uniqchars if "EA" in i]
    uniqEA = [("EA", "EA")]
    chEA = AllCan.loc[AllCan["chars"].isin(uniqEA)].sort_values(corvalnm)
    for n, r in chEA.tail(20).iterrows():
        plot(Allchars, x="WL_precursor", y=n[1])
    AllC_500

    AllC_500 = get_corrs_E(AllCan, 600).sort_values(corvalnm)

    uniqEA = [("B", "B")]
    chEA = AllCan.loc[AllCan["chars"].isin(uniqEA)].sort_values(corvalnm)
    for n, r in chEA.head(20).iterrows():
        plot(Allchars, x=n[0], y=n[1])


def get_corrs_E(AllCorrs_test, Ev):
    lst, sets = [], []
    for ix in AllCorrs_test.index:
        pt = []
        for n, pos in enumerate(ix):
            last = pos.split("_")[-1]
            match = re.match(r"^[0-9]{3}\b", last)
            if match:
                if float(match.group(0)) == Ev:
                    pt.append([n, pos])  # Not a Potential
                #                    lst += ix
                else:
                    pass  # Wrong Potential Ev
            else:
                pt.append([n, pos])
        if len(pt) == 2:
            idx = (pt[0][1], pt[1][1])
            _ridx = (pt[1][1], pt[0][1])
            #            sets.append(set(idx))
            if not _ridx in lst:
                lst.append(idx)
    return AllCorrs_test.loc[lst]


def make_uniform_EvRHE(df, rounding_set=2):
    lst = []
    for E in df["E_AppV_RHE"].values:
        match = 0
        for i in np.arange(-0.10, 1.2, 0.05):
            if np.isclose(i, E, atol=0.025):
                match = 1
                lst.append((E, i))
        if match == 0:
            if E < 0 and E > -0.04:
                lst.append((E, i))
            else:
                print(E, i)

    if len(df["E_AppV_RHE"].values) == len(lst):
        df = df.assign(**{"E_RHE": [np.round(float(i[1]), rounding_set) for i in lst]})
        print(
            'Len({0}) matches, new column: "E_RHE"'.format(len(df["E_AppV_RHE"].values))
        )
    else:
        print(
            "make_uniform_EvRHE lengths do not match LenEIS : {0}, len(lst) : {1}".format(
                len(df["E_AppV_RHE"].values), len(lst)
            )
        )
    return df


# def plot_chars_Ev(Allchars)
#         PlotAnalysis.corrplot(allchflrtd[allcorrs.name].unstack().dropna())


if __name__ == "__main__":

    EvRHE = "E_AppV_RHE"
    templ = PorphSiO2_template()


def old_test():
    ECpars, raman, bet, ea, prep, struc, char, Allchars = collect_pars()
    export_chars(char)
    export_chars({"all_chars": Allchars})

    list(struc.columns)
    list(ECpars.columns)
    corrs = corrloop(char)

    AllCorrs = allchars_corr(Allchars)
    corr_prep_ea = postCorrs.prep_ea(prep, ea)
    corr_prep_ea = postCorrs.prep_ea(ea, raman)

    makeplots = 0
    if not ECpars.empty and makeplots == True:

        for i in np.arange(100, 1000, 100):
            plot(ECpars, x=f"ORR_Frac_H2O2_{i}", y="HPRR_dj/dE", savefig=plotsfolder)

        plot(Allchars, x="EA_N_content", y="ORR_Frac_H2O2_400")
        plot(ECpars, x="ORR_Frac_H2O2_500", y="HPRR_dj/dE")

        plot(ECpars, x="ORR_Frac_H2O2_500", y="HPRR_E_onset")
        plot(ECpars, x="ORR_Frac_H2O2_500", y="ORR_Jkin_min_500")
        plot(ECpars, x="N2_Cdl_mFcm-2_700", y="ORR_Jkin_min_750")
        plot(ECpars, x="N2_Cdl_mFcm-2_700", y="ORR_Frac_H2O2_500")

        plot(ECpars, x="N2_Cdl_mFcm-2_500", y="ORR_n_ORR_500")
        #        === POWERPOINT slides
        ppslides = mkfolder(FindExpFolder("PorphSiO2").compare.joinpath("pp_slides"))
        plot(ECpars, x="HPRR_E_onset", y="HPRR_dj/dE", savefig=ppslides)
        plot(ECpars, x="ORR_E_onset", y="ORR_Jkin_min_750", logy=True, savefig=ppslides)

        plot(ECpars, x="ORR_E_onset", y="HPRR_E_onset", savefig=ppslides)
        plot(ECpars, x="HPRR_dj/dE", y="ORR_Jkin_min_750", savefig=ppslides)

        plot(ECpars, x="ORR_Frac_H2O2_500", y="ORR_Jkin_min_750", savefig=ppslides)
        plot(ECpars, x="ORR_Frac_H2O2_500", y="HPRR_E_onset", savefig=ppslides)

        plot(ECpars, x="N2_Cdl_mFcm-2_500", y="HPRR_dj/dE")  # , savefig=plotsfolder)
        ECpars["relHPRR"] = ECpars["HPRR_dj/dE"] / ECpars["N2_Cdl_mFcm-2_500"]
        plot(ECpars, x="ORR_Jkin_min_750", y="relHPRR")  # , savefig=plotsfolder)

        plot(ECpars, x="N2_Cdl_mFcm-2_500", y="ORR_J_diff_lim", savefig=plotsfolder)

        plot(
            Allchars,
            y="B_BET_RAW_calc_constant",
            x="P_MW_g-mol",
            savefig=plotsfolder,
            logy=True,
        )

        #        plot(Allchars,x='B_tBoer_microArea',y='P_MW_content',savefig=plotsfolder)
        plot(ECpars, x="HPRR_E_onset", y="HPRR_dj/dE", savefig=plotsfolder)
        plot(
            ECpars,
            x="ORR_E_onset",
            y="ORR_Jkin_min_750",
            logy=False,
            savefig=plotsfolder,
        )
        plot(
            Allchars,
            x="B_BET_RAW_calc_SA m2/g",
            y="HPRR_dj/dE",
            logy=False,
            savefig=plotsfolder,
        )
        plot(
            Allchars,
            x="B_BET_RAW_calc_SA m2/g",
            y="P_WL_precursor",
            logy=False,
            savefig=plotsfolder,
        )


#    ECpars = EC_PorphSiO2().pars
#    pst = postChar()
#    raman,bet,ea = pst.raman, pst.bet, pst.ea
#    struc = pst.merged


#
#    fig,ax= plt.subplots(figsize= (16,16))
#    scatter_matrix(ECpars[[key for key,val in ECpars.dtypes.to_dict().items() if 'float64' in str(val)]], alpha=0.5, ax=ax, diagonal='kde')
#    plt.savefig(FindExpFolder('PorphSiO2').compare.joinpath('EC_scmatrix.png'))
#    plt.close()
