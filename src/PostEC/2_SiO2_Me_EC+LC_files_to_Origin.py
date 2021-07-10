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
import pickle

# import types
# import post_helper
# import plotting
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.stats import linregress, zscore
import pandas as pd
import numpy as np
import datetime as dt
import pandas as pd

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
    from collect_load import Load_from_Indexes, CollectLoadPars

    #    from FileHelper.FindExpFolder import FindExpFolder
    from plotting import eisplot
    from prep_postchar import postChar
    import EIS_export

elif "prepare_input" in __name__:
    pass
#    import RunEC_classifier
#    from FileHelper.FindSampleID import FindSampleID

import logging

_logger = logging.getLogger(__name__)

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


def EC_types_grp():
    #    KL ['ORR_E_AppV_RHE', 'ORR_KL_E_AppV_RHE','Electrode']
    _basic_EC_cond = ["postAST_post", "Sweep_Type", "pH", "Loading_cm2"]
    _extra_EC_cond = {
        "N2CV": [],
        "N2": [],
        "ORR": ["RPM_DAC_uni"],
        "KL": ["Electrode", "ORR_E_AppV_RHE"],
        "EIS": ["E_RHE"],
        "HER": ["HER_RPM_post"],
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

    # globals EC_index
    #    ['Model(Singh2015_RQRQ)', 'Model(Singh2015_RQRQR)', 'Model(Bandarenka_2011_RQRQR)',
    #                  'Model(Singh2015_RQRWR)', 'Model(Randles_RQRQ)', 'Model(Singh2015_R3RQ)']
    #    model_select = EC_PorphSiO2.EIS_models[1]
    #    self = EC_PorphSiO2()
    def __init__(self):
        # self.index, self.AST_days = EC_PorphSiO2.select_ECexps(EC_folder)
        self.select_EC_ASTexps_from_ECindex()

    #        self.pars = EC_PorphSiO2.mergedEC()
    #        self.par_export = EC_OHC.to_excel(self.folder.joinpath('EC_ORR_HPRR.xlsx'))
    def select_EC_ASTexps_from_ECindex(self):
        EC_idx_PorphSiO2_samples = EC_index.loc[
            EC_index.SampleID.isin(self.Porph_template.SampleID.unique())
        ]
        # pd.read_excel(list(EC_folder.rglob('*EC_index*'))[0])

        EC_idx_PorphSiO2_samples = EC_idx_PorphSiO2_samples.assign(
            **{
                "PAR_date_day_dt": [
                    dt.date.fromisoformat(np.datetime_as_string(np.datetime64(i, "D")))
                    for i in EC_idx_PorphSiO2_samples.PAR_date.to_numpy()
                ]
            }
        )
        self.EC_idx_PorphSiO2_samples = EC_idx_PorphSiO2_samples
        self.get_AST_days()
        # LC_idx_fp = list(EC_folder.rglob('*EC_index*'))[0].parent.joinpath('LC_index.xlsx')
        EC_idx_PorphSiO2_AST = EC_idx_PorphSiO2_samples.loc[
            EC_idx_PorphSiO2_samples.PAR_date_day_dt.isin(
                [i for a in self.AST_days.to_numpy() for i in a]
            )
        ]
        # AST_days = EC_PorphSiO2.get_AST_days()
        # EC_idx_PorphSiO2_AST.to_excel(list(EC_folder.rglob('*EC_index*'))[0].parent.joinpath('LC_index.xlsx'))
        self.EC_idx_PorphSiO2 = EC_idx_PorphSiO2_AST
        # if LC_idx_fp.exists():
        # else:
        #     try:
        #         LC_fls = pd.read_excel(LC_idx_fp,index_col=[0])
        #     except Exception as e:
        #         print(f'Excel load fail: {e}\n,file: {LC_idx_fp}')
        #         LC_fls = pd.DataFrame()
        # return LC_fls, AST_days

    def get_AST_days(self):
        gr_idx = self.EC_idx_PorphSiO2_samples.groupby("PAR_date_day_dt")
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
        #        (AST_days[-1][0], AST_days[0][1])
        #        AST_days.append((dt.date(2019,5,6), dt.date(2019,1,25)))
        #        AST_days.append((dt.date(2019,5,6), dt.date(2019,1,26)))
        _extra_AST_days = [
            (dt.date(2019, 5, 6), dt.date(2019, 1, 25)),
            (dt.date(2019, 5, 6), dt.date(2019, 1, 26)),
        ]
        AST_days += _extra_AST_days
        AST_days = pd.DataFrame(
            AST_days, columns=["PAR_date_day_dt_pre", "PAR_date_day_dt_post"]
        )
        AST_days = AST_days.assign(
            **{
                "PAR_date_day_dt_diff": AST_days.PAR_date_day_dt_pre
                - AST_days.PAR_date_day_dt_post
            }
        )
        self.AST_days = AST_days

    # def select_ECexps(EC_folder):
    #     LC_idx_fp = list(EC_folder.rglob('*EC_index*'))[0].parent.joinpath('LC_index.xlsx')
    #     AST_days = EC_PorphSiO2.get_AST_days()
    #     if LC_idx_fp.exists():
    #         LC_fls = EC_PorphSiO2.EC_idx_PorphSiO2.loc[EC_PorphSiO2.EC_idx_PorphSiO2.PAR_date_day_dt.isin([i for a in AST_days.to_numpy() for i in a])]
    #         LC_fls.to_excel(list(EC_folder.rglob('*EC_index*'))[0].parent.joinpath('LC_index.xlsx'))
    #     else:
    #         try:
    #             LC_fls = pd.read_excel(LC_idx_fp,index_col=[0])
    #         except Exception as e:
    #             print(f'Excel load fail: {e}\n,file: {LC_idx_fp}')
    #             LC_fls = pd.DataFrame()
    #     return LC_fls, AST_days
    # def repr_index(self):
    #     PAR_exp_uniq = {grn : len(grp) for grn,grp in self.index.groupby("PAR_exp")}
    #     print(f'Len({len(self.index)},\n{PAR_exp_uniq}')


def _testing_():
    tt = EC_prepare_EC_merged(reload_AST=True, reload_merged=True, reload_pars=True)
    self = tt

    N2CV = self.N2cv(reload=False, use_daily=True)


#%% == EC_prepare_EC_merged == testing
class EC_prepare_EC_merged:
    EIS_models = EIS_export.EIS_selection.mod_select
    # ['Model(EEC_Randles_RWpCPE)', 'Model(EEC_2CPE)', 'Model(EEC_2CPEpW)',
    # 'Model(EEC_RQ_RQ_RW)', 'Model(EEC_RQ_RQ_RQ)', 'Model(Randles_RQRQ)']
    ORR_reload = dict(reload=True, use_daily=False)
    ORR_no_reload = dict(reload=False, use_daily=True)
    use_daily = True
    # global ParsColl
    # ParsColl = ParsColl
    mcols = [i for i in Load_from_Indexes.EC_label_cols if i not in ["PAR_file"]] + [
        "Sweep_Type"
    ]
    _pkl_EC_merged = "EC_merged_dict"

    def __init__(self, reload_AST=False, reload_merged=False, reload_pars=True):
        self.reload_AST = reload_AST
        self.reload_merged = reload_merged
        self.reload_pars = reload_pars

        self.set_pars_collection()

        self.reload_pars_kws = dict(reload=reload_pars, use_daily=self.use_daily)
        self.EC_merged_dict = {}
        self.load_EC_PorphSiO2()
        self.load_merged_EC()

    def set_pars_collection(self):
        if "ParsColl" in globals().keys():
            self.ParsColl = ParsColl
        else:

            Pars_Collection = CollectLoadPars(load_type="fast")
            # globals()['Pars_Collection'] = Pars_Collection
            ParsColl = Pars_Collection.pars_collection
            self.ParsColl = ParsColl

    def load_EC_PorphSiO2(self):
        self.EC_PorphSiO2 = EC_PorphSiO2()
        self.AST_days = self.EC_PorphSiO2.AST_days
        self.EC_idx_PorphSiO2 = self.EC_PorphSiO2.EC_idx_PorphSiO2

    def load_merged_EC(self):
        if self.reload_merged:
            self.reload_merged_EC()

        if not self.EC_merged_dict:
            _load_EC_merge = load_dict_pkl(self._pkl_EC_merged)
            if _load_EC_merge:
                self.EC_merged_dict = _load_EC_merge

    def reload_merged_EC(self):
        try:
            self.load_N2CV()
            self.load_ORR()
            self.load_KL()
            self.load_EIS()
            self.load_HER()
            self.add_filter_selection_of_EC_merged()
            save_dict_pkl(self._pkl_EC_merged, self.EC_merged_dict)
        except Exception as e:
            _logger.warning(f"EC_prepare_EC_merged, reload_merged_EC failure: {e}")

    def get_AST_matches(self, DF, _verbose=False):
        # LC_fls, AST_days = EC_PorphSiO2.select_ECexps(EC_folder)
        #        DF = ORR.drop_duplicates()
        #        DF = N2CV.drop_duplicates()
        #        DF = EIS.drop_duplicates()
        #        DF = HER.drop_duplicates()
        #        DF = ttpars
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
            DF.PAR_date_day_dt = pd.to_datetime(DF.PAR_date_day_dt, unit="D")
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
        #        AST_days_run_lst = [i for i in AST_days if len(i) == 2][-1:]
        for n, r in self.AST_days.iterrows():
            #            if len(_dates) == 2:
            #                _pre,_post = _dates
            #            elif (len_dates) == 1:
            _pre, _post = r.PAR_date_day_dt_pre, r.PAR_date_day_dt_post
            _preslice = DF.loc[
                (DF.PAR_date_day == _pre.strftime("%Y-%m-%d")) & (DF.postAST == "no")
            ]
            pre = _preslice.groupby(_compare_cols)

            _postslice = DF.loc[
                (DF.PAR_date_day == _post.strftime("%Y-%m-%d")) & (DF.postAST != "no")
            ]
            post = _postslice.groupby(_compare_cols)
            _res = {}
            _res = {
                "pre_PAR_date_day_dt": _pre,
                "post_PAR_date_day_dt": _post,
                "AST_days_n": n,
            }
            #            print(_res,[_preslice.postAST.unique()[0], _postslice.postAST.unique()[0]])
            union = set(pre.groups.keys()).union(set(post.groups.keys()))

            matches = set(pre.groups.keys()).intersection(set(post.groups.keys()))

            _difference_pre = set(pre.groups.keys()).difference(set(post.groups.keys()))
            _difference_post = set(post.groups.keys()).difference(
                set(pre.groups.keys())
            )
            #            _diffr.append((_pre,_post,_difference_pre, _difference_post))
            if not _preslice.empty and not _postslice.empty:
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
                            if _postAST in "postAST_sHA|postAST_LC" and _verbose:
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

    #                            prgrp.groupby(['Sweep_Type','RPM_DAC']).groups
    #                            prgrp['ORR_Jkin_min_700']-pogrp['ORR_Jkin_min_700']

    def load_N2CV(self):
        N2CV = self.edit_pars_N2cv(**self.reload_pars_kws)
        #        N2_pltqry = EC_merged_dict.get('N2CV')
        N2_AST = self.get_AST_matches(N2CV)
        N2_AST_diff = self.compare_AST_pars(N2CV, N2_AST, reload=self.reload_AST)
        #        _DFtype = EC_PorphSiO2.sense_DF_type(N2CV)
        #        EC_merged_dict.update({'N2CV' : N2_AST_diff})
        self.EC_merged_dict.update(
            {"N2CV": {"PARS": N2CV, "AST_matches": N2_AST, "AST_diff": N2_AST_diff}}
        )

    def load_ORR(self, _testing=False):
        ORR = self.edit_pars_ORR()
        ORR_AST = self.get_AST_matches(ORR)
        ORR_AST_diff = self.compare_AST_pars(ORR, ORR_AST, reload=self.reload_AST)
        if _testing:
            ttpars = ORR.query('RPM_DAC_uni > 1000 & Sweep_Type == "cathodic"')
            tt_AST = self.get_AST_matches(ttpars)
            tt = ORR_AST.query('RPM_DAC_uni > 1000 & Sweep_Type == "cathodic"')
            tt_diff = self.compare_AST_pars(ORR, tt, reload=reload_AST, save_pkl=False)
        #        ttpfs = ORR.loc[ORR.ORR_Jkin_max_700 > 0].PAR_file.unique()
        #        ttpfs = ORR.query('Sweep_Type == "mean"').loc[ORR.ORR_E_onset > 0.85].PAR_file.unique()
        #                ORR.loc[(ORR.ORR_E_onset > 0.85) & (ORR.Sweep_Type == 'cathodic')].PAR_file.unique()
        #        EC_merged_dict.update({'ORR' : ORR_AST_diff})
        self.EC_merged_dict.update(
            {"ORR": {"PARS": ORR, "AST_matches": ORR_AST, "AST_diff": ORR_AST_diff}}
        )

    def load_KL(self):
        KL = self.edit_pars_KL()
        KL = KL.assign(**{"RPM_DAC": 1500})
        KL_AST = self.get_AST_matches(KL)
        KL_AST_diff = self.compare_AST_pars(KL, KL_AST, reload=self.reload_AST)
        #        EC_merged_dict.update({'KL' : KL_AST_diff})
        self.EC_merged_dict.update(
            {"KL": {"PARS": KL, "AST_matches": KL_AST, "AST_diff": KL_AST_diff}}
        )

    def load_EIS(self):
        EIS = self.edit_pars_EIS()
        EIS_AST = self.get_AST_matches(EIS)
        EIS_AST_diff = self.compare_AST_pars(EIS, EIS_AST, reload=self.reload_AST)
        #        EC_merged_dict.update({'EIS' : EIS_AST_diff})
        self.EC_merged_dict.update(
            {"EIS": {"PARS": EIS, "AST_matches": EIS_AST, "AST_diff": EIS_AST_diff}}
        )

    def load_HER(self):
        HER = self.edit_pars_HER()
        HER_type_grp = HER.groupby("HER_type")
        HER.HER_at_E_slice = HER.HER_at_E_slice.round(3)
        HER_AST = self.get_AST_matches(HER)
        for Htype, Hgrp in HER_type_grp:
            #            Htype, Hgrp = 'E_slice', HER.loc[HER.groupby('HER_type').groups['E_slice']]
            HER_AST_diff = self.compare_AST_pars(
                Hgrp, HER_AST, reload=self.reload_AST, extra=Htype
            )
            try:
                if not HER_AST_diff.empty:
                    self.EC_merged_dict.update(
                        {
                            f"HER_{Htype}": {
                                "PARS": Hgrp,
                                "AST_matches": HER_AST,
                                "AST_diff": HER_AST_diff,
                            }
                        }
                    )
            except Exception as e:
                print(f"HER {Htype} fail, {e}")

    #            EC_merged_dict.update({f'HER_{Htype}' : HER_AST_diff})

    def finish_EC_merged(self):
        _pkl_EC_merged = "EC_merged_dict"
        EC_merged_dict = EC_PorphSiO2.add_filter_selection_of_EC_merged(EC_merged_dict)
        save_dict_pkl(_pkl_EC_merged, EC_merged_dict)
        # EC_merged_dict = load_dict_pkl(_pkl_EC_merged)

    def add_filter_selection_of_EC_merged(self):

        _drop_AST_row_pre = [
            "2019-01-25;N2_20cls_300_100_10_JOS5_256;no;0",
            "2019-01-25;N2_20cls_300_100_10_JOS4_256;no;0",
        ]
        _check_cols = [
            "SampleID",
            "AST_row",
            "PAR_date_day_dt_pre",
            "PAR_date_day_dt_post",
            "postAST_post",
        ]
        _srt2 = ["postAST_post", "SampleID"]
        _ch_match = [
            "SampleID",
            "pre_PAR_date_day_dt",
            "post_PAR_date_day_dt",
            "post_postAST",
            "pre_postAST",
        ]
        _sortcols = ["SampleID", "post_postAST"][::-1]
        pd.set_option("display.max_columns", 6)
        pd.set_option("display.width", 100)
        for _EC, _DF in self.EC_merged_dict.items():
            #             _EC, _DF =  'N2CV', EC_merged_dict['N2CV']
            #             _EC, _DF =  'ORR', EC_merged_dict['ORR']
            # print(_EC)
            if "AST_row_n" not in _DF["AST_diff"].columns:
                _DF["AST_diff"]["AST_row_n"] = [
                    int(i[-1]) for i in _DF["AST_diff"].AST_row.str.split("_").values
                ]

            AST_diff = _DF["AST_diff"].copy()
            AST_diff.loc[~AST_diff.AST_row_pre.isin(_drop_AST_row_pre)]
            AST_matches = (
                _DF["AST_matches"].copy().sort_values(by=["post_postAST", "SampleID"])
            )
            _rem1 = AST_matches.loc[
                (AST_matches.post_postAST == "postAST_LC")
                & (AST_matches.SampleID.isin(["JOS2", "JOS4", "JOS5"]))
                & (AST_matches.pre_PAR_date_day_dt == dt.date(2019, 1, 25))
            ].assign(**{"rem": 1})

            _rem2 = AST_matches.loc[
                (
                    (AST_matches.post_postAST == "postAST_LC")
                    & (AST_matches.pre_postAST == "no")
                    & (
                        AST_matches.SampleID.isin(
                            ["JOS1", "JOS2", "JOS3", "JOS4", "JOS5"]
                        )
                    )
                    & (AST_matches.pre_PAR_date_day_dt == dt.date(2019, 5, 6))
                    & (
                        AST_matches.post_PAR_date_day_dt.isin(
                            [dt.date(2019, 1, 25), dt.date(2019, 1, 26)]
                        )
                    )
                )
            ].assign(**{"rem": 2})
            #             _keep_clean.loc[2].to_dict()
            #            _jos3 = {'SampleID': 'JOS3', 'pre_PAR_date_day_dt': dt.date(2019, 1, 24), 'post_PAR_date_day_dt': dt.date(2019, 1, 25),
            #                     'post_postAST': 'postAST_LC', 'pre_postAST': 'no'}
            #            _jos3qry = ' & '.join([f'{k} == {repr(val)}' for k,val in _jos3.items()])
            #            AST_matches.query(_jos3qry)
            _rem3 = AST_matches.loc[
                (
                    (AST_matches.post_postAST == "postAST_LC")
                    & (AST_matches.pre_postAST == "no")
                    & (AST_matches.SampleID.isin(["JOS3"]))
                    & (AST_matches.pre_PAR_date_day_dt == dt.date(2019, 1, 24))
                    & (AST_matches.post_PAR_date_day_dt == dt.date(2019, 1, 25))
                )
            ].assign(**{"rem": 3})

            _rem4 = AST_matches.loc[(AST_matches.pre_postAST != "no")].assign(
                **{"rem": 4}
            )

            _edit = _rem1  # changed later 15.03

            _remove = pd.concat([_rem2, _rem4, _rem3])

            _keep = AST_matches.iloc[~AST_matches.index.isin(_remove.index.values)]
            AST_matches[_ch_match].drop_duplicates()
            _rem_clean = _remove[_ch_match + ["rem"]].sort_values(by=_sortcols)
            _keep_clean = _keep[_ch_match].sort_values(by=_sortcols)
            #             _remove[['SampleID','post_postAST']] # check
            #             _rem = _DF['AST_diff'].loc[_DF['AST_diff']['AST_row_n'].isin(_remove.index.values)]
            #             _rem[['SampleID','postAST_post','PAR_date_day_pre']] #check
            _filtered = AST_diff.loc[~AST_diff["AST_row_n"].isin(_remove.index.values)]
            #             DF['AST_diff'] = _filtered
            self.EC_merged_dict.update({_EC: {**_DF, **{"AST_diff_filter": _filtered}}})

        print(
            f'EC merged dict updated with dropped rows in "AST_diff_filter" for:\n {self.EC_merged_dict.keys()}'
        )
        # return EC_merged_dict

    #             _DF['AST_diff'].loc[_DF['AST_diff'].AST_row_n.isin(_rem]

    def EC_merge_postchar(_reloadset=False):

        _pkl_EC_postchar = "EC_merged_postchars"
        EC_postchars = load_dict_pkl(_pkl_EC_postchar)
        if not EC_postchars and _reloadset != True:
            EC_merged_dict = EC_PorphSiO2.mergedEC(_reloadset=True)
            #        EC_merged_dict_bak = EC_merged_dict.copy()
            EC_merged_dict = EC_PorphSiO2.add_filter_selection_of_EC_merged(
                EC_merged_dict
            )
            postChars = postChar().merged
            _extracols = [i for i in SampleCodes.columns if not "SampleID" in i]
            EC_postchars = {}
            for _EC, _DF_dict in EC_merged_dict.items():
                _DF = _DF_dict["AST_diff_filter"]
                _initcols = _DF.columns
                _DF = _DF.dropna(axis=1, how="all")
                _DF = _DF.drop(columns=_DF.columns.intersection(_extracols))
                _DF = pd.merge(_DF, postChars, on="SampleID")
                _postcols = _DF.columns
                EC_postchars.update({_EC: _DF})
            save_dict_pkl(_pkl_EC_postchar, EC_postchars)
        return EC_postchars

    def _fix_ORR_scans():
        EC_postchars = EC_PorphSiO2.EC_merge_postchar(_reloadset=True)
        _ORR = EC_postchars["ORR"]
        _J245 = _ORR.loc[
            _ORR.SampleID.isin(["JOS2,", "JOS4", "JOS5"])
            & (_ORR.postAST_post == "postAST_LC")
        ]
        _extracols = [i for i in SampleCodes.columns if not "SampleID" in i]

    def compare_AST_pars(self, _DF, _AST_in, reload=False, extra="", save_pkl=True):
        #        _DF, _AST_in = EIS, EIS_AST
        #        _DF, _AST_in = N2CV, N2_AST
        #        _DF, _AST_in = ORR, ORR_AST
        #        _DF, _AST_in = KL, KL_AST
        #        _DF, _AST_in = HER, HER_AST
        #        _DF, _AST_in = Hgrp, HER_AST
        #        _DF, _AST_in = ttN2CV, ttAST
        #        reload, extra = _reloadset, Htype
        _DF = _DF.drop_duplicates()
        _DFtype = self.sense_DF_type(_DF)
        _DFtype = "".join([i for i in _DFtype if str.isalpha(i)])

        _DFtype_prefix = _DFtype.split("_")[0]
        if extra:
            _pklpath = EC_PorphSiO2.folder.joinpath(
                f"AST_compared_pars_{_DFtype}_{extra}.pkl"
            )
        else:
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
            #            _AST = _AST_in.loc[(_AST_in.SampleID == "JOS4") ]
            # & (_AST_in.post_postAST.str.contains('LC'))]

            _DF_diff_out = []
            _errors = []
            _dms = []
            #            _AST.loc[_AST.SampleID == "JOS4"].tail(2)
            for n, r in _AST.iterrows():
                #                r[_dropnacols]
                _uniq_AST_row_pre = f"{r.pre_PAR_date_day_dt};{Path(r.pre_PAR_file).stem};{r.pre_postAST};{int(r.pre_dupli_num)}"
                _uniq_AST_row_post = f"{r.post_PAR_date_day_dt};{Path(r.post_PAR_file).stem};{r.post_postAST};{int(r.post_dupli_num)}"
                #            EIS.query('  and  '.join([f'{k} == {repr(v)}' for k, v in _pred.items()]))
                _pred = dict(zip(_precols, r[_prec].to_dict().values()))
                _preQ = "  &  ".join(
                    [f"{k} == {repr(v)}" for k, v in _pred.items() if k in _DF.columns][
                        1:
                    ]
                )
                _Dpre = _DF.query(_preQ).dropna(axis=1, how="all")

                _postd = dict(zip(_postcols, r[_post].to_dict().values()))
                _postQ = "  &  ".join(
                    [
                        f"{k} == {repr(v)}"
                        for k, v in _postd.items()
                        if k in _DF.columns
                    ][1:]
                )
                _Dpos = _DF.query(_postQ).dropna(axis=1, how="all")
                _dms.append((n, _pred, _postd))
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
                #                _dms.append((n, len(_Dm), _Dm ))

                _mcols = [
                    i[0]
                    for i in set(_0).intersection(set(_1))
                    if not i[0].startswith("dupli")
                ]
                _mcols = [
                    i
                    for i in _mcols
                    if i not in ["PAR_exp", "Dest_dir"] and not i.startswith("EXP_")
                ]
                _mcols.sort()

                _othercols = _Dpos.columns.difference(_mcols)
                t2 = _Dpos[_othercols]

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
                elif "HER" in _DFtype:
                    _addcols = [
                        i
                        for i in [
                            "HER_type",
                            "HER_at_J_slice",
                            "HER_at_E_slice",
                            "HER_Segnum",
                        ]
                        if i in set(_Dpre.columns).union(_Dpos.columns)
                    ]
                    _mcols += _addcols

                _Dm = pd.merge(_Dpre, _Dpos, on=_mcols, suffixes=["_pre", "_post"])
                _Dm = _Dm.assign(
                    **{
                        "AST_row": f"{_DFtype}_{n}",
                        "AST_row_n": int(n),
                        "AST_days_n": r.AST_days_n,
                        "AST_row_pre": _uniq_AST_row_pre,
                        "AST_row_post": _uniq_AST_row_post,
                    }
                )
                #                [(i, _Dpos[i].nunique(), _Dpos[i].unique()[0], _Dpre[i].nunique(), _Dpre[i].unique()[0], (_Dpos[i].unique(),_Dpre[i].unique()))
                #                    for i in _mcols if _Dpos[i].nunique() > 1]
                if _Dm.empty:
                    run_this
                #                    try:
                #                        _Dm = pd.merge_asof(_Dpre.sort_values(_mcols), _Dpos.sort_values(_mcols), on = _mcols, suffixes = ['_pre','_post'])
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
            if save_pkl == True:
                DF_diff.to_pickle(_pklpath)
            _logger.info(f"AST compare len({len(DF_diff)}) saved to:{_pklpath}")
            return DF_diff

    # DF_diff.groupby(['postAST_post','SampleID']).plot(x='E_RHE', y='EIS_Rct_O2_diff_abs',ylim=(-200,200))

    def sense_DF_type(self, _DF):
        #        _c = [i[0] for i in Counter([i.split('_')[0] for i in _DF.columns]).most_common(5) if i[0] not in ['BET','tM']][0]
        _excl = set(self.EC_idx_PorphSiO2.columns).union(SampleCodes.columns)
        _res = [
            i
            for i in Counter(
                ["_".join(i.split("_")[0:2]) for i in _DF.columns]
            ).most_common(20)
            if not any([i[0] in b for b in _excl]) and i[0][0].isalnum()
        ]

        _res2 = Counter(["_".join(i.split("_")[0:1]) for i, c in _res])
        _type = _res2.most_common(1)[0][0]
        _extraC = Counter(
            ["_".join(i.split("_")[1:2]) for i in _DF.columns if _type in i]
        ).most_common(1)
        if _extraC[0][1] > 4:
            _type = f"{_type}_{_extraC[0][0]}"
        #            if _res2.most_common(2)[1][1] > 3:
        #                _type = f'{_type}_{_res2.most_common(2)[1][0]}'
        return _type

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
                EC_prepare_EC_merged.mcols
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
    def edit_pars_HPRR(sweep_type_select=["anodic", "cathodic"]):
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

    def load_pars_HER(self):
        HER_pars_all = Load_from_Indexes.HER_pars_OVV(**self.reload_pars_kws)
        self.pars_HER = HER_pars_all

    @edit_columns
    def edit_pars_HER(self, sweep_type_select=["anodic", "cathodic"], unit="F"):
        # reload= False, use_daily = True,  extra_plotting=False, xls_out = False
        # LC_idx = self.index

        if (
            not Pfolder.joinpath("HER_orig_data.pkl").exists()
            or self.reload_pars == True
        ):
            self.load_pars_HER()
            HER_pars = self.pars_HER.loc[
                (
                    (self.pars_HER._type == "HER_pars")
                    & (self.pars_HER.PAR_file.isin(self.index.PAR_file.to_numpy()))
                )
            ]
            HER_pars.to_pickle(Pfolder.joinpath("HER_orig_data.pkl"))
        else:
            HER_pars = pd.read_pickle(Pfolder.joinpath("HER_orig_data.pkl"))

        HER_pars = HER_pars.dropna(how="all", axis=1)

        return HER_pars, "HER"

    def load_pars_ORR(self):
        ORR_pars_all = self.ParsColl["ORR_pars"]
        # Load_from_Indexes.ORR_pars_OVV(**self.reload_pars_kws)
        self.pars_ORR = ORR_pars_all

    @edit_columns
    def edit_pars_ORR(self):

        if not hasattr(self, "pars_ORR"):
            self.load_pars_ORR()

        ORR_pars = self.pars_ORR.loc[
            (
                (self.pars_ORR.source_type == "ORR_pars")
                & (
                    self.pars_ORR.PAR_file.isin(
                        self.EC_idx_PorphSiO2.PAR_file.to_numpy()
                    )
                )
            )
        ]
        ORR_pars = ORR_pars.dropna(how="all", axis=1)
        # Adding log cols to ORR pars
        ORR_pars = ORR_pars.assign(
            **{
                f'{"_".join(i.split("_")[0:-1])}_log_{i.split("_")[-1]}': np.log(
                    ORR_pars[i]
                )
                for i in ORR_pars.columns
                if "Jkin" in i
            }
        )
        return ORR_pars, "ORR"

    @edit_columns
    def edit_pars_KL(self):

        if not hasattr(self, "pars_ORR"):
            self.load_pars_ORR()

        KL_pars = self.pars_ORR.loc[
            (
                (self.pars_ORR.source_type == "KL_pars")
                & (
                    self.pars_ORR.PAR_file.isin(
                        self.EC_idx_PorphSiO2.PAR_file.to_numpy()
                    )
                )
            )
        ]
        KL_pars = KL_pars.dropna(how="all", axis=1)
        return KL_pars, "ORR"

    def load_pars_N2CV(self):
        # N2_loadpars = N2_LoadPars(reload = True, reload_raw = False )
        Cdl_pars_all = self.ParsColl["N2_pars"]
        # N2_loadpars.N2_pars
        # Load_from_Indexes.N2_pars_OVV(**self.reload_pars_kws)
        # (reload= self.reload_pars, use_daily = use_daily,  extra_plotting=extra_plotting, xls_out = xls_out)
        self.pars_N2CV = Cdl_pars_all

    @edit_columns
    def edit_pars_N2cv(
        self,
        sweep_type_select=["anodic", "cathodic"],
        unit="F",
        reload=False,
        use_daily=True,
        extra_plotting=False,
        xls_out=False,
    ):

        self.load_pars_N2CV()

        if not Pfolder.joinpath("N2_orig_data.pkl").exists() or reload == True:
            Cdl_pars_all = self.pars_N2CV
            Cdl_pars = Cdl_pars_all.loc[
                Cdl_pars_all.PAR_file.isin(self.EC_idx_PorphSiO2.PAR_file.to_numpy())
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

    def load_pars_EIS(self):
        _source = "Load pars"
        if "files" in _source:
            eis_files, eisfs = (
                list(
                    self.folder.parent.joinpath(f"EIS_Porph_SiO2\{model_select}").rglob(
                        "JOS*.xlsx"
                    )
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
            EIS_pars_mod = self.ParsColl["EIS_pars"]
            # Load_from_Indexes.EIS_pars_OVV(reload= False, extra_plotting=False, xls_out = False, use_daily = True, use_latest=True)
            # EIS_pars_mod = EIS_pars.loc[EIS_pars.Model_EEC.isin(self.EIS_models.values())]
        self.pars_EIS = EIS_pars_mod

    @edit_columns
    def edit_pars_EIS(self, _source="Load pars"):
        """Models used are selected in the EIS_export module
        via dict from EC_PorphSiO2.EIS_models"""

        self.load_pars_EIS()
        EIS_pars_mod = self.pars_EIS
        EIS_pars_mod = EIS_pars_mod.loc[
            EIS_pars_mod.index.isin(EIS_pars_mod.best_mod_index)
        ]
        _sample_uniq_cols1 = set(
            [
                a
                for n, gr in EIS_pars_mod.groupby("SampleID")
                for a in [i for i in gr.columns if gr[i].nunique() == 1]
            ]
        )
        _sample_uniq_cols2 = set(
            [
                a
                for n, gr in EIS_pars_mod.groupby("SampleID")
                for a in [i for i in gr.columns if gr[i].nunique() == 2]
            ]
        )
        _sample_uniq_cols2.difference(_sample_uniq_cols1)
        # Create EIS var columns with gas N2 or O2 as suffix names
        # EPgrp = EIS_pars_mod.groupby(['Gas','Model_EEC'])
        EPgrp_gas = EIS_pars_mod.groupby(["Gas"])

        # N2grp = ('N2',self.EIS_models.get('N2'))
        # O2grp = ('O2',self.EIS_models.get('O2'))
        # EP_N2,EP_O2 = EPgrp.get_group(N2grp).drop(columns='Gas'), EPgrp.get_group(O2grp).drop(columns='Gas')

        EC_exp_index = [
            i for i in Load_from_Indexes.EC_label_cols if i not in ["PAR_file", "Gas"]
        ] + ["PAR_date_day"]
        _gasgrp = []

        for gas in ["O2", "N2"]:
            # gasn_grp = (gas,self.EIS_models.get(gas))
            grp = EPgrp_gas.get_group(gas)
            _varsgrp = [a for i in grp.lmfit_var_names.unique() for a in i.split(", ")]
            _varsgrp += ["Rct_kin" for i in _varsgrp if "Rct" in i] + [
                "Qad+Cdlp"
                for i in _varsgrp
                if all([i in _varsgrp for i in ["Qad", "Cdlp"]])
            ]
            _sample_uniq_cols1 = set([i for i in grp.columns if grp[i].nunique() == 1])
            #            grp.lmfit_var_names.unique()[0].split(', ')
            _grp = grp.rename(columns={i: i + f"_{gas}" for i in set(_varsgrp)})
            #            _grp = _grp.drop(columns='Gas')
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

    def EIS_spectra_origin_prep(model_select=["Model(R0-L0-p(R1-Ws1,CPE1)-C2)"]):
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

    def EIS_spectra_origin(model_select=["Model(R0-L0-p(R1-Ws1,CPE1)-C2)"]):
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

    # #    @save_load_AST_pars
    #     def mergedEC( _reloadset = False):
    #         _pkl_EC_merged = 'EC_merged_dict'
    # #        EC_merged_dict = EC_PorphSiO2.mergedEC(_reloadset=True)
    #         if _reloadset == True:
    # #        EC_merged_dict_bak = EC_merged_dict.copy()
    # #            EC_merged_dict =  EC_PorphSiO2.take_selection_of_EC_merged(EC_merged_dict)
    #             mcols = [i for i in Load_from_Indexes.EC_label_cols if i not in ['PAR_file']]+['Sweep_Type']
    #             _mcols = [i for i in mcols if not i in ['Gas','E_RHE']]
    #             LC_fls, AST_days = EC_PorphSiO2.select_ECexps(EC_folder)

    #             EC_merged_dict = {}
    #     #        _reloadset = True
    #             template = PorphSiO2_template()
    #             HPRR = EC_PorphSiO2.HPRR()

    #             N2CV = EC_PorphSiO2().N2cv(reload= False, use_daily = True)
    #     #        N2_pltqry = EC_merged_dict.get('N2CV')
    #             N2_AST = EC_PorphSiO2.get_AST_matches(N2CV)
    #             N2_AST_diff = EC_PorphSiO2.compare_AST_pars(N2CV, N2_AST, reload = False)
    #     #        _DFtype = EC_PorphSiO2.sense_DF_type(N2CV)
    #     #        EC_merged_dict.update({'N2CV' : N2_AST_diff})
    #             EC_merged_dict.update({'N2CV' : {'PARS' : N2CV, 'AST_matches' : N2_AST, 'AST_diff' : N2_AST_diff}})
    #             # list(N2CV.columns)
    #     #        _renameN2 = {c : c.split('_')[-1] for c in [i for i in N2CV.columns if any([i.split('_')[-1] in mcols])]}
    #     #        N2CV = N2CV.rename(columns = _renameN2)
    #             ORR = EC_PorphSiO2().ORR_pars()
    #             ORR_AST = EC_PorphSiO2.get_AST_matches(ORR)
    #             ORR_AST_diff = EC_PorphSiO2.compare_AST_pars(ORR, ORR_AST, reload = _reloadset)
    #             ttpars = ORR.query('RPM_DAC_uni > 1000 & Sweep_Type == "cathodic"')
    #             tt_AST = EC_PorphSiO2.get_AST_matches(ttpars)
    #             tt = ORR_AST.query('RPM_DAC_uni > 1000 & Sweep_Type == "cathodic"')
    #             tt_diff = EC_PorphSiO2.compare_AST_pars(ORR, tt, reload = _reloadset, save_pkl = False)
    #     #        ttpfs = ORR.loc[ORR.ORR_Jkin_max_700 > 0].PAR_file.unique()
    #     #        ttpfs = ORR.query('Sweep_Type == "mean"').loc[ORR.ORR_E_onset > 0.85].PAR_file.unique()
    #     #                ORR.loc[(ORR.ORR_E_onset > 0.85) & (ORR.Sweep_Type == 'cathodic')].PAR_file.unique()
    #     #        EC_merged_dict.update({'ORR' : ORR_AST_diff})
    #             EC_merged_dict.update({'ORR' : {'PARS' : ORR, 'AST_matches' : ORR_AST, 'AST_diff' : ORR_AST_diff}})
    #     #        _renameO2 = {c : c.split('_')[-1] for c in [i for i in ORR.columns if any([i.split('_')[-1] in mcols]) and not '_Ring' in i]}
    #     #        ORR = ORR.rename(columns = _renameO2)
    #             KL = EC_PorphSiO2().KL_pars()
    #             KL = KL.assign(**{'RPM_DAC' : 0})
    #             KL_AST = EC_PorphSiO2.get_AST_matches(KL)
    #             KL_AST_diff = EC_PorphSiO2.compare_AST_pars(KL, KL_AST, reload = _reloadset)
    #     #        EC_merged_dict.update({'KL' : KL_AST_diff})
    #             EC_merged_dict.update({'KL' : {'PARS' : KL, 'AST_matches' : KL_AST, 'AST_diff' : KL_AST_diff}})
    #     #        _KLdatacols = ['ORR_KL_data_file_post','ORR_KL_data_x_post', 'ORR_KL_data_y_post', 'ORR_KL_fit_y_post', 'ORR_KL_fit_y_2e_post', 'ORR_KL_fit_y_4e_post']
    #     #        _renameKL = {c : c.split('_')[-1] for c in [i for i in KL.columns if any([i.split('_')[-1] in mcols]) and not '_Ring' in i]}
    #     #        KL = KL.rename(columns = _renameKL)
    #             EIS = EC_PorphSiO2.EIS_pars()
    #             EIS_AST = EC_PorphSiO2.get_AST_matches(EIS)
    #             EIS_AST_diff = EC_PorphSiO2.compare_AST_pars(EIS, EIS_AST, reload = _reloadset)
    #     #        EC_merged_dict.update({'EIS' : EIS_AST_diff})
    #             EC_merged_dict.update({'EIS' : {'PARS' : EIS, 'AST_matches' : EIS_AST, 'AST_diff' : EIS_AST_diff}})
    #     #        _renameEIS = {c : c.split('_')[-1] for c in [i for i in EIS.columns if any([i.split('_')[-1] in mcols]) and not '_Ring' in i]}
    #     #        EIS = EIS.rename(columns = _renameEIS)
    #             HER = EC_PorphSiO2().HER_pars(reload= False, use_daily = True)
    #             HER_type_grp = HER.groupby('HER_type')
    #             HER.HER_at_E_slice = HER.HER_at_E_slice.round(3)
    #             HER_AST = EC_PorphSiO2.get_AST_matches(HER)
    #             for Htype, Hgrp in HER_type_grp:
    #     #            Htype, Hgrp = 'E_slice', HER.loc[HER.groupby('HER_type').groups['E_slice']]
    #                 HER_AST_diff = EC_PorphSiO2.compare_AST_pars(Hgrp, HER_AST, reload = _reloadset,extra= Htype)
    #                 try:
    #                     if not HER_AST_diff.empty:
    #                         EC_merged_dict.update({f'HER_{Htype}' : {'PARS' : Hgrp, 'AST_matches' : HER_AST, 'AST_diff' : HER_AST_diff}})
    #                 except Exception as e:
    #                     print(f'HER {Htype} fail, {e}')
    #     #            EC_merged_dict.update({f'HER_{Htype}' : HER_AST_diff})
    #             EC_merged_dict =  EC_PorphSiO2.add_filter_selection_of_EC_merged(EC_merged_dict)
    #             save_dict_pkl(_pkl_EC_merged, EC_merged_dict)
    #         else:
    #             EC_merged_dict = load_dict_pkl(_pkl_EC_merged)
    #         return EC_merged_dict

    #        ECmerged = pd.merge(ORR,pd.merge(N2CV, EIS,on=_mcols),on=_mcols)
    #        EC_EIS = pd.merge(ECmerged,EIS,on=mcols)
    #        EC_OHN_merged = pd.merge(template, EC_EIS, on='SampleID')
    #        EC_PorphSiO2.export_to_xls(EC_OHN_merged)
    #        return EC_OHN_merged

    def corr_plots():
        EC_OHC.query('SampleID != "JOS5"').corr()
        corrstk = EC_OHC.query('SampleID != "JOS5"').corr().stack()
        EC_OHC.plot(x="E_onset", y="HPRR_E_onset", kind="scatter")
        EC_OHC.plot(x="FracH2O2_050", y="HPRR_E_onset", kind="scatter")
        EC_OHC.plot(x="N2_Cdl_mFcm-2_0.5", y="HPRR_dj/dE", kind="scatter")
        EC_OHC.plot(x="N2_Cdl_mFcm-2_0.5", y="E_half", kind="scatter")

        EC_OHC.corr(method="pearson")

    def _check_eis_plots():
        _par = ["Cdlp", "Rorr", "Rct", "Qad", "Aw"][-1]
        _checky = ["N_content", "BET_cat_agg"][0]
        for modn, mgrp in EIS_pars_all.loc[EIS_pars_all.pH < 3].groupby(
            ["pH", "postAST", "Model_EEC"]
        ):
            _ps = eisplot(_par)
            if len(mgrp[_par].dropna()) > 3:
                mgrp.plot(
                    x="E_RHE",
                    y=_par,
                    yerr=f"{_par}_stderr",
                    kind="scatter",
                    ylim=_ps.ylim,
                    logy=_ps.logy,
                    title=f"{modn}",
                    c=_checky,
                    cmap="rainbow",
                )


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
def _testing_():
    tt = EC_prepare_EC_merged()
    self = tt
    _pp = EC_post_plotting(tt)
    self = _pp
    N2CV = self.N2cv(reload=False, use_daily=True)


#%% == EC_post_plotting == testing
class EC_post_plotting:
    def __init__(self, _EC_prepare_EC_merged):

        self._EC_prepare_EC_merged = _EC_prepare_EC_merged
        self.add_attrs()

    def add_attrs(self):
        if hasattr(self._EC_prepare_EC_merged, "EC_merged_dict"):
            self.EC_merged = self._EC_prepare_EC_merged.EC_merged_dict
        else:
            self.EC_merged = {}  # self._EC_prepare_EC_merged

    def ORR_get_experiments(self):

        ORR_AST = self.EC_merged["ORR"]["AST_matches"]
        ORR_AST_mean1500 = ORR_AST.loc[
            (ORR_AST.Sweep_Type == "mean") & (ORR_AST.RPM_DAC_uni > 1000)
        ]
        ORR_AST_mean1500.to_excel(EC_folder.joinpath("ORR_AST_exp_overview.xlsx"))
        # N2_scan_index = EC_index.loc[(EC_index.SampleID.isin(_smpls)) & (EC_index.PAR_exp.str.contains('N2_act'))]
        # N2_scan_index.to_excel(EC_folder.joinpath('N2_scan_exp_overview.xlsx'))

    def N2_repr_Cdl(self):

        # ECname = 'N2'
        # Cdl_pars_all = Load_from_Indexes.N2_pars_OVV()
        _DF = self.EC_merged["N2CV"]["PARS"]

        # _DF = Cdl_pars_all

        ECname = "N2"
        _raw_data_folder = mkfolder(EC_folder.joinpath(f"{ECname}_reproducibility"))
        _grpcols = ["pH", "Loading_cm2", "SampleID"]
        _swpcol = [i for i in _DF.columns if "Sweep_Type" in i]
        _grpcols += _swpcol

        _sIDgrps = _DF.loc[
            _DF.SampleID.isin(PorphSiO2_template().SampleID.values) & (_DF.pH < 2)
        ]
        #            .query('postAST == "no"')
        _lst = []
        for sID, sgrp in _sIDgrps.groupby(_grpcols):
            sID,

            _sgpr_Cdl_mean = (
                sgrp.groupby("E_AppV_RHE").Cdl.mean().rename("Cdl_grp_mean")
            )
            _Cdl_cols = [i for i in sgrp.columns if i.startswith("N2_Cdl_F")]

            fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 10), sharex=True)
            _sgpr_Cdl_mean.plot(
                c="grey", alpha=0.5, ls="--", lw=5, label="mean Cdl", ax=ax1
            )
            pfgrps = sgrp.groupby("PAR_file")
            for pf, pfgrp in pfgrps:
                pfgrp = pd.merge(pfgrp, _sgpr_Cdl_mean, on=EvRHE)

                ls = "-" if "no" in pfgrp.postAST.unique() else "--"

                pfgrp = pfgrp.assign(
                    **{"Cdl_mean_diff": pfgrp.Cdl - pfgrp.Cdl_grp_mean}
                )
                _lst.append(pfgrp)

                pfgrp.plot(
                    x="E_AppV_RHE", y="Cdl_mean_diff", ax=ax2, legend=False, ls=ls
                )
                _dt = pfgrp.PAR_date_day.unique()[0]
                _lbl = f"{_dt}, {Path(pf).stem}"
                pfgrp.plot(x="E_AppV_RHE", y="Cdl", ax=ax1, label=_lbl, ls=ls)
            _title = ", ".join([f"{k} : {str(val)}" for k, val in (zip(_grpcols, sID))])
            _stem = "_".join([str(i) for i in sID]) + f"_{len(pfgrps)}"

            ax1.set_ylabel("Cdl")
            ax1.set_title(_title)
            ax1.legend(
                fontsize=15, bbox_to_anchor=(1.02, 1), loc="upper left", fancybox=True
            )
            ax2.set_ylabel("Cdl - Cdl_mean")
            #            ax2.legend(False)
            plt.savefig(
                _raw_data_folder.joinpath(_stem + ".png"), bbox_inches="tight", dpi=200
            )
            plt.close()
        N2_Cdl_pars_mean = pd.concat(_lst)

        def select_sID_N2(_sIDgrps):

            _grp_select = (1.0, 0.379, "JOS4", "cathodic")
            _jos4 = _sIDgrps.groupby(_grpcols).get_group(_grp_select)
            _raw_data_folder = mkfolder(
                EC_folder.joinpath(
                    f"{ECname}_reproducibility", "_".join([str(i) for i in _grp_select])
                )
            )

            _j4lc = _jos4.loc[_jos4.postAST == "postAST_LC"]
            j4post = pd.concat(
                [
                    pd.read_excel(_j4lc.sourceFilename.unique()[0].parent.joinpath(i))
                    for i in _j4lc.N2_CV_datafilenames.unique()[0].split(", ")
                ]
            )
            _j4no = _jos4.loc[
                (_jos4.PAR_date_day == "2019-05-06") & (_jos4.postAST == "no")
            ]

            j4no = pd.concat(
                [
                    pd.read_excel(_j4no.SourceFilename.unique()[0].parent.joinpath(i))
                    for i in _j4no.N2_CV_datafilenames.unique()[0].split(", ")
                ]
            )

            _j4no_pfgrps = _jos4.loc[(_jos4.postAST == "no")].groupby("PAR_file")
            for pf, pgr in _j4no_pfgrps:
                j4no_grps = pd.concat(
                    [
                        pd.read_excel(pgr.SourceFilename.unique()[0].parent.joinpath(i))
                        for i in pgr.N2_CV_datafilenames.unique()[0].split(", ")
                    ]
                ).groupby("ScanRate_mVs")
                for sr, sgrp in j4post.groupby("ScanRate_mVs"):
                    fig, ax = plt.subplots()
                    j4no_grps.get_group(sr).plot(
                        x="E_AppV_RHE",
                        y="jmAcm-2",
                        ax=ax,
                        label=f"pre,{pgr.PAR_date_day.unique()[0]} / {Path(pgr.PAR_file.unique()[0]).stem}",
                    )
                    sgrp.plot(
                        x="E_AppV_RHE", y="jmAcm-2", ax=ax, label="postAST_LC", title=sr
                    )
                    ax.legend(
                        fontsize=15,
                        bbox_to_anchor=(1.02, 1),
                        loc="upper left",
                        fancybox=True,
                    )
                    _stem = f"{sr}_{pgr.PAR_date_day.unique()[0]}_{Path(pgr.PAR_file.unique()[0]).stem}"
                    plt.savefig(
                        _raw_data_folder.joinpath(_stem + ".png"),
                        bbox_inches="tight",
                        dpi=200,
                    )

    def reproducibility_check_samples(_DF, ECname):

        ECname = "EIS"
        if ECname == "EIS":
            EIS = EC_PorphSiO2.EIS_pars()
            _DF = EIS

            _grpcols = ["pH", "Loading_cm2", "SampleID"]
            _eisqry = '(postAST == "no") &'
            _sIDgrps = _DF.loc[
                _DF.SampleID.isin(PorphSiO2_template().SampleID.values) & (_DF.pH < 2)
            ].query('(Sweep_Type == "cathodic")')

            _lst = []
            for sID, sgrp in _sIDgrps.groupby(_grpcols):
                sID
                pars, vars = eisplot.read_varnames(sgrp)
                for gas in ["O2", "N2"]:
                    #                    gas ='N2'
                    _vars_gas = [i for i in vars if i.endswith(gas)]
                    sgrp_gas = sgrp.copy()
                    sgrp_gas.dropna(subset=_vars_gas, axis=0, inplace=True)
                    sgrp_gas.dropna(axis=1, how="all", inplace=True)

                    sgrp_gas = sgrp_gas.loc[
                        sgrp_gas[[i for i in _vars_gas if f"Rct_{gas}" in i][0]] < 2000
                    ]

                    for var in _vars_gas:
                        _raw_data_folder = mkfolder(
                            EC_folder.joinpath(f"{ECname}_reproducibility/{var}")
                        )
                        #                    var = f'EIS_{_var}_{gas}'
                        _sgpr_var_mean = (
                            sgrp_gas.groupby("E_RHE")[var].mean().rename(f"{var}_mean")
                        )
                        fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 10), sharex=True)
                        _sgpr_var_mean.plot(
                            c="grey",
                            alpha=0.5,
                            ls="--",
                            lw=5,
                            label=f"{var}_mean",
                            ax=ax1,
                        )
                        pfgrps = sgrp_gas.groupby("PAR_file")
                        for pf, pfgrp in pfgrps:
                            pfgrp = pd.merge(pfgrp, _sgpr_var_mean, on="E_RHE")
                            pfgrp = pfgrp.assign(
                                **{
                                    f"{var}_mean_diff": pfgrp[var]
                                    - pfgrp[f"{var}_mean"]
                                }
                            )
                            _lst.append(pfgrp)
                            pfgrp.plot(
                                x="E_RHE", y=f"{var}_mean_diff", ax=ax2, legend=False
                            )
                            _dt = pfgrp.PAR_date_day.unique()[0]
                            _lbl = f"{_dt}, {Path(pf).stem}"
                            pfgrp.plot(
                                x="E_RHE", y=var, ax=ax1, label=_lbl
                            )  # TODO add yerr=
                        _title = (
                            ", ".join(
                                [f"{k} : {str(val)}" for k, val in (zip(_grpcols, sID))]
                            )
                            + f", par : {var}"
                        )
                        _stem = (
                            var
                            + "_"
                            + "_".join([str(i) for i in sID])
                            + f"_{len(pfgrps)}"
                        )
                        # TODO CHECK REPR
                        ax1.set_ylabel(var)
                        ax1.set_title(_title)
                        ax1.legend(
                            fontsize=15,
                            bbox_to_anchor=(1.02, 1),
                            loc="upper left",
                            fancybox=True,
                        )
                        ax2.set_ylabel(f"{var}_mean_diff")
                        #            ax2.legend(False)
                        plt.savefig(
                            _raw_data_folder.joinpath(_stem + ".png"),
                            bbox_inches="tight",
                            dpi=200,
                        )
                        plt.close()
            EIS_diff_means = pd.concat(_lst)

    def get_raw_data(EC_merged_dict, ECname):
        EC_merged_dict = EC_PorphSiO2.mergedEC(_reloadset=False)
        #        _CVdfls, _mcols, _setAST, pAST, ref_dict = _srcfls_fit, _mcols, _setAST, pAST, _uniq_id
        def _read_excel_df(_CVdfls, _mcols, _setAST, pAST, ref_dict={}):
            if _CVdfls:
                _stCVdata = pd.concat(
                    [pd.read_excel(i, index_col=[0]) for i in _CVdfls],
                    sort=False,
                    ignore_index=True,
                ).dropna(axis=1, how="all")
                _stCVdata = _stCVdata.rename(
                    columns={
                        i: f"{i}_{status}" for i in _stCVdata.columns if i not in _mcols
                    }
                )

                _good_cols = [
                    key
                    for key, val in ref_dict.items()
                    if key in _stCVdata.columns and _stCVdata[key].unique()[0] == val
                ]
                _bad_cols = [
                    key
                    for key, val in ref_dict.items()
                    if key in _stCVdata.columns and _stCVdata[key].unique()[0] != val
                ]
                _missing_cols = {
                    key: val
                    for key, val in ref_dict.items()
                    if key not in _stCVdata.columns
                }

                _stCVdata = _stCVdata.assign(
                    **{**{"postAST": _setAST, "postAST_post": pAST}, **_missing_cols}
                )

            else:
                _stCVdata = pd.DataFrame()
            return _stCVdata

        if "N2CV" in EC_merged_dict.keys():
            ECname = "N2"
            _raw_data_folder = mkfolder(EC_folder.joinpath(f"{ECname}_raw_data"))
            _loadDF = read_load_pkl(f"{ECname}_raw_data")

            if _loadDF.empty:
                AST_grp_cols = ["postAST_post", "Sweep_Type", "pH", "Loading_cm2"]
                N2_pltqry = EC_merged_dict.get("N2CV").groupby(AST_grp_cols)
                _mcols = [
                    EvRHE,
                    "ScanRate_mVs",
                    "Sweep_Type",
                    "SampleID",
                    "postAST",
                ] + ["postAST_post", "postAST_pre"]
                _grpcols = ("cathodic", 1.0, 0.379)

                _pAST_opts = ["postAST_LC", "postAST_sHA"]
                _st = []
                for pAST in _pAST_opts:
                    ASTgrp = N2_pltqry.get_group((pAST, *_grpcols))
                    _samplegrpcols = ["SampleID", "AST_row"]
                    for n, gr in ASTgrp.groupby(_samplegrpcols):
                        #                        n = list(ASTgrp.groupby(_samplegrpcols).groups)[-2]
                        #                        gr = ASTgrp.groupby(_samplegrpcols).get_group(n)
                        # TODO find missing JOS4 postAST_LC N2_act!
                        for status in ["pre", "post"]:
                            _sourcedfls = [
                                i
                                for i in gr[f"N2_SourceFilename_{status}"].unique()
                                if pd.notna(i)
                            ]
                            _setAST = "no" if "pre" in status else pAST
                            _CVdfls = [
                                f.parent.joinpath(e)
                                for f in _sourcedfls
                                for e in [
                                    a
                                    for i in gr[
                                        f"N2_CV_datafilenames_{status}"
                                    ].unique()
                                    for a in i.split(", ")
                                ]
                                if f.parent.joinpath(e).is_file()
                            ]
                            _uniq_id = dict(
                                zip(
                                    AST_grp_cols + _samplegrpcols,
                                    (_setAST, *_grpcols, *n),
                                )
                            )
                            if _CVdfls:
                                _stCVdata = _read_excel_df(
                                    _CVdfls, _mcols, _setAST, pAST, ref_dict=_uniq_id
                                )
                                #                                _stCVdata = pd.concat([pd.read_excel(i,  index_col=[0]) for i in _CVdfls],sort=False,ignore_index=True)
                                #                                _stCVdata = _stCVdata.assign(**{'postAST' :_setAST, 'postAST_post' : pAST })
                                #                                _stCVdata = _stCVdata.rename(columns = {i : f'{i}_{status}' for i in _stCVdata.columns if i not in _mcols})
                                _st.append(_stCVdata)
                            else:
                                print("Sources empty!!", n, status, pAST)

                N2_CVs = pd.concat([i for i in _st])
                save_DF_pkl(f"{ECname}_raw_data", N2_CVs)
            else:
                N2_CVs = _loadDF

            _select_cols = [
                c
                for c in N2_CVs.columns
                if any(
                    [
                        i in c
                        for i in [
                            EvRHE,
                            "jmAcm-2",
                            "postAST",
                            "Sweep_Type",
                            "pH",
                            "Loading_cm2",
                        ]
                    ]
                )
            ]

            _ScanRates_check_lim = [
                (
                    n,
                    gr[["jmAcm-2_pre", "jmAcm-2_post"]].max(),
                    gr[["jmAcm-2_pre", "jmAcm-2_post"]].min(),
                )
                for n, gr in N2_CVs.groupby(["ScanRate_mVs"])
            ]
            _ScanRates_ylim = {
                10: (-2.5, 2.5),
                100: (-11, 7),
                150: (-13, 10),
                200: (-16, 12),
                300: (-22, 16),
            }
            for n, gr in N2_CVs.groupby(["SampleID", "AST_row", "ScanRate_mVs"]):
                n, gr
                _name = "_".join([str(i) for i in (*_grpcols, *n)])
                fig, ax = plt.subplots(figsize=(6, 4))
                for _past, pgr in gr.groupby("postAST"):
                    _jcol = "pre" if _past == "no" else "post"
                    pgr.plot(
                        x=EvRHE,
                        y=f"jmAcm-2_{_jcol}",
                        title=f"{_name}",
                        label=f"{_jcol}, {_past}",
                        ax=ax,
                    )
                ax.set_ylim(_ScanRates_ylim.get(n[-1]))
                plt.savefig(
                    _raw_data_folder.joinpath(f"{_name}.png"),
                    dpi=100,
                    bbox_inches="tight",
                )
                plt.close()
                print(_name)
                gr[_select_cols].dropna(axis=1, how="all").to_excel(
                    _raw_data_folder.joinpath(f"{_name}.xlsx")
                )
            # TODO Check new plots with Update Pre scans for ORR and for FINAL PLOTS!!
        if "ORR" in EC_merged_dict.keys():
            ECname = "ORR"
            _raw_data_folder = mkfolder(EC_folder.joinpath(f"{ECname}_raw_data_3"))
            _loadDF = read_load_pkl(f"{ECname}_raw_data")
            _mcols = [EvRHE, "Sweep_Type", "SampleID", "postAST"] + [
                "postAST_post",
                "postAST_pre",
            ]
            _grpcols = (1500, "mean", 1)
            print(f"ORR selection: {_grpcols}")
            KL_AST_diff = EC_merged_dict.get("KL").get("AST_diff")
            if _loadDF.empty:
                ORR_pltqry = (
                    EC_merged_dict.get("ORR")
                    .get("AST_diff")
                    .query("RPM_DAC_uni > 1000")
                    .groupby(["postAST_post", "RPM_DAC_uni", "Sweep_Type", "pH"])
                )

                _st = []
                for pAST in ["postAST_LC", "postAST_sHA"]:
                    ASTgrp = ORR_pltqry.get_group((pAST, *_grpcols))
                    for n, gr in ASTgrp.groupby("SampleID"):
                        for status in ["pre", "post"]:
                            gr.ORR_RRDE_swp_data_file_pre
                            _sourcedfls = [
                                i
                                for i in gr[f"ORR_RRDE_swp_data_file_{status}"].unique()
                                if pd.notna(i)
                            ]
                            _setAST = "no" if "pre" in status else pAST
                            #                            _CVdfls = [f.parent.joinpath(e) for f in _sourcedfls
                            #                                        for e in [a for i in gr[f'ORR_datafilenames_{status}'].unique() for a in i.split(', ')]
                            #                                        if f.parent.joinpath(e).is_file() ]
                            if _sourcedfls:
                                _stCVdata = _read_excel_df(
                                    _sourcedfls, _mcols, _setAST, pAST
                                )
                                #                                _stCVdata = pd.concat([pd.read_excel(i,  index_col=[0]) for i in _sourcedfls],sort=False,ignore_index=True).dropna(axis=1,how='all')
                                #                                _stCVdata = _stCVdata.assign(**{'postAST' :_setAST, 'postAST_post' : pAST })
                                #                                _stCVdata = _stCVdata.rename(columns = {i : f'{i}_{status}' for i in _stCVdata.columns if i not in _mcols})
                                _st.append(_stCVdata)

                ORR_data = pd.concat([i for i in _st])
                save_DF_pkl(_raw_data_folder, ORR_data)
            else:
                ORR_data = _loadDF

        def ORR_prepost_AST_4plots():
            KL_grp = KL_AST_diff.query('Electrode == "KL_I_Disk"').groupby(
                ["postAST_post", "SampleID", "Sweep_Type"]
            )
            """ problem in data, there are double PAR_files in the _pre data
            01.03: fixed the plottings for pre, need to choose which pre version to take still...
            11.03: fixed and chose the ORR scans with N2 BGs, all use N2_jcorr, still check which ure used here...."""
            _take_sweep_type = "mean"
            # .query('postAST_post == "postAST_LC" ')
            for n, gr in ORR_data.groupby(["postAST_post", "SampleID"]):
                n, gr
                # gr = gr.sort_values(EvRHE,ascending=True)
                try:
                    _nKL = KL_grp.get_group((*n, gr.Sweep_Type.unique()[0]))
                    #                    _nKL = _nKL.assign(**{'ORR_KL_2e' : 2, 'ORR_KL_4e' : 4})
                    _nKL.to_excel(_raw_data_folder.joinpath(f"KL_{_name}.xlsx"))
                except:
                    print(f"no KL plot for: {n}")
                    _nKL = pd.DataFrame()

                _metal = PorphSiO2_template().query("SampleID == @n[-1]").Metal.iloc[0]
                # fig,axes = plt.subplots(2,2,figsize=(22/2.54,18/2.54))
                for prePF, pfgr in gr.groupby("PAR_file_disk_pre"):
                    prePF, pfgr

                    fig, axes = plt.subplots(2, 2, figsize=(22 / 2.54, 18 / 2.54))
                    # gr.sort_values(EvRHE,ascending=True).plot(x=EvRHE,y=[f'Jcorr_{_stat}'],ax=axes[1][0],ylim=(-6,0.5),xlim=(0,1))
                    _prepost_date = []
                    for _stat in ["pre", "post"]:
                        if _stat == "pre":

                            _date = [
                                np.datetime_as_string(i, unit="D")
                                for i in pfgr[f"EXP_date_day_dt_{_stat}"].unique()
                                if not pd.isna(i)
                            ][0]
                            _lbl = f"{_stat}_{_date}"
                            _ckws = {"c": "b", "label": _lbl}
                            pfgr.sort_values(EvRHE, ascending=True).plot(
                                x=EvRHE,
                                y="Jcorr_pre",
                                ax=axes[1][0],
                                ylim=(-6, 0.5),
                                xlim=(0, 1),
                                **_ckws,
                            )
                            # pfgr.plot(x=EvRHE,y='Jcorr_pre',ax=axes[1][0],ylim=(-6,0.5),xlim=(0,1),**_ckws)
                            pfgr.plot(
                                x=EvRHE,
                                y="Frac_H2O2_pre",
                                ylim=(0, 30),
                                ax=axes[0][0],
                                xlim=(0, 1),
                                **_ckws,
                            )
                            pfgr.plot(
                                x=EvRHE,
                                y=f"Jkin_min_{_stat}",
                                ylim=(0.01, 30),
                                logy=True,
                                xlim=(0.5, 0.9),
                                ax=axes[1][1],
                                **_ckws,
                            )
                            if not _nKL.empty:
                                _nKL.query("PAR_file_pre == @prePF").plot(
                                    x=f"ORR_{EvRHE}",
                                    y=f"ORR_nElectrons_{_stat}",
                                    ax=axes[0][1],
                                    ylim=(0, 7),
                                    xlim=(0, 1),
                                    **_ckws,
                                )

                            _prepost_date.append(_lbl)

                        elif _stat == "post":

                            _date = [
                                np.datetime_as_string(i, unit="D")
                                for i in gr[f"EXP_date_day_dt_{_stat}"].unique()
                                if not pd.isna(i)
                            ][0]
                            _lbl = f"{_stat}_{_date}"
                            _ckws = {"c": "r"}
                            gr.plot(
                                x=EvRHE,
                                y=f"Jcorr_{_stat}",
                                ax=axes[1][0],
                                ylim=(-6, 0.5),
                                xlim=(0, 1),
                                label=f"Jcorr_{_lbl}",
                                **_ckws,
                            )
                            gr.plot(
                                x=EvRHE,
                                y=f"Frac_H2O2_{_stat}",
                                ylim=(0, 30),
                                ax=axes[0][0],
                                xlim=(0, 1),
                                label=f"FracH2O2_{_lbl}",
                                **_ckws,
                            )
                            gr.plot(
                                x=EvRHE,
                                y=f"Jkin_min_{_stat}",
                                ylim=(0.01, 30),
                                logy=True,
                                xlim=(0.5, 0.9),
                                ax=axes[1][1],
                                label=f"Jkin_min_{_lbl}",
                                **_ckws,
                            )
                            if not _nKL.empty:
                                _nKL.sort_values(f"ORR_{EvRHE}", ascending=True).plot(
                                    x=f"ORR_{EvRHE}",
                                    y=f"ORR_nElectrons_{_stat}",
                                    ax=axes[0][1],
                                    ylim=(0, 7),
                                    xlim=(0, 1),
                                    **_ckws,
                                )

                            _prepost_date.append(_lbl)

                    axes[1][0].axhline(y=0, color="black", ls="--", alpha=0.2)
                    axes[0][1].axhline(y=2, color="grey", ls="-.", alpha=0.2)
                    axes[0][1].axhline(y=4, color="grey", ls="-.", alpha=0.4)

                    _name = "_".join([str(i) for i in (*_grpcols, *n, *_prepost_date)])
                    fig.suptitle(f"{_metal}, {_name}")
                    plt.savefig(
                        _raw_data_folder.joinpath(f"{_name}.png"),
                        dpi=100,
                        bbox_inches="tight",
                    )
                    plt.close()
                    gr.to_excel(_raw_data_folder.joinpath(f"{_name}.xlsx"))

        def old_ORR_plot():
            gr.plot(
                x=EvRHE,
                y=["Jcorr_pre", "Jcorr_post"],
                ax=axes[1][0],
                ylim=(-6, 0.5),
                xlim=(0, 1),
            )
            gr.plot(
                x=EvRHE,
                y=["Frac_H2O2_pre", "Frac_H2O2_post"],
                ylim=(0, 30),
                ax=axes[0][0],
                xlim=(0, 1),
            )
            gr.plot(
                x=EvRHE,
                y=["Jkin_min_pre", "Jkin_min_post"],
                ylim=(0.01, 30),
                logy=True,
                xlim=(0.5, 0.9),
                ax=axes[1][1],
            )

            _name = "_".join([str(i) for i in (*_grpcols, *n)])
            # TODO implement KL plotting raw data
            try:
                _nKL = KL_grp.get_group((*n, gr.Sweep_Type.unique()[0]))
                #                    _nKL = _nKL.assign(**{'ORR_KL_2e' : 2, 'ORR_KL_4e' : 4})
                _nKL.to_excel(_raw_data_folder.joinpath(f"KL_{_name}.xlsx"))
                _nKL.plot(
                    x=f"ORR_{EvRHE}",
                    y=["ORR_nElectrons_pre", "ORR_nElectrons_post"],
                    ax=axes[0][1],
                    ylim=(0, 7),
                    xlim=(0, 1),
                )
                #                    _nKL.plot(x=f'ORR_{EvRHE}',y=['ORR_KL_4e'],ax=axes[0][1],**{'ls' : '-.'},alpha=0.2, c='grey', label=None)
                #                    _nKL.plot(x=f'ORR_{EvRHE}',y=['ORR_KL_2e'],ax=axes[0][1],**{'ls' : '-.'},alpha=0.2, c='grey', label=None)
                axes[0][1].axhline(y=2, color="grey", ls="-.", alpha=0.2)
                axes[0][1].axhline(y=4, color="grey", ls="-.", alpha=0.4)
            except:
                print(f"no KL plot for: {n}")
            [ax.legend(fontsize=11) for ax1 in axes for ax in ax1]

            plt.savefig(
                _raw_data_folder.joinpath(f"{_name}.png"), dpi=100, bbox_inches="tight"
            )
            plt.close()
            gr.to_excel(_raw_data_folder.joinpath(f"{_name}.xlsx"))

        def ORR_plot_all():
            KL_grp = KL_AST_diff.query('Electrode == "KL_I_Disk"').groupby(
                ["postAST_post", "SampleID", "Sweep_Type"]
            )
            templ = PorphSiO2_template()
            fig, axes = plt.subplots(2, 2, figsize=(22 / 2.54, 18 / 2.54))
            for n, gr in ORR_data.groupby(["postAST_post", "SampleID"]):
                n, gr
                # gr = gr.sort_values(EvRHE,ascending=True)
                _metal = templ.query("SampleID == @n[-1]").Metal.iloc[0]
                _cn = templ.query("SampleID == @n[-1]").color.iloc[0]
                _RGB = OriginColors.iloc[_cn].RGB_255
                print(_metal, _cn, _RGB)
                _nKL = KL_grp.get_group((*n, gr.Sweep_Type.unique()[0]))
                _name = "_".join([str(i) for i in (*_grpcols, *n)])
                for _AST in ["pre", "post"]:
                    _ls = "-" if "pre" in _AST else ":"
                    kws = {"ls": _ls, "c": [*_RGB, 0.6]}
                    gr.plot(
                        x=EvRHE,
                        y=f"Jcorr_{_AST}",
                        ax=axes[1][0],
                        ylim=(-6, 0.5),
                        xlim=(0, 1),
                        **kws,
                    )
                    gr.plot(
                        x=EvRHE,
                        y=f"Frac_H2O2_{_AST}",
                        ylim=(0, 30),
                        ax=axes[0][0],
                        xlim=(0, 1),
                        **kws,
                    )
                    gr.plot(
                        x=EvRHE,
                        y=f"Jkin_min_{_AST}",
                        ylim=(0.01, 30),
                        logy=True,
                        xlim=(0.5, 0.9),
                        ax=axes[1][1],
                        **kws,
                    )
                    # TODO implement KL plotting raw data
                    try:
                        #                    _nKL = _nKL.assign(**{'ORR_KL_2e' : 2, 'ORR_KL_4e' : 4})
                        _nKL.plot(
                            x=f"ORR_{EvRHE}",
                            y=f"ORR_nElectrons_{_AST}",
                            ax=axes[0][1],
                            ylim=(0, 7),
                            xlim=(0, 1),
                        )
                        #                    _nKL.plot(x=f'ORR_{EvRHE}',y=['ORR_KL_4e'],ax=axes[0][1],**{'ls' : '-.'},alpha=0.2, c='grey', label=None)
                        #                    _nKL.plot(x=f'ORR_{EvRHE}',y=['ORR_KL_2e'],ax=axes[0][1],**{'ls' : '-.'},alpha=0.2, c='grey', label=None)
                        axes[0][1].axhline(y=2, color="grey", ls="-.", alpha=0.2)
                        axes[0][1].axhline(y=4, color="grey", ls="-.", alpha=0.4)
                    except:
                        print(f"no KL plot for: {n}")
            _lgnds = [ax.legend(fontsize=11) for ax1 in axes for ax in ax1]

            plt.legend(_lgnds, ["jA"])
            fig.suptitle(f'{_metal}, {", ".join(n)}')
            axes[1][0].axhline(y=0, color="black", ls="--", alpha=0.2)
            plt.savefig(
                _raw_data_folder.joinpath(f"{_name}.png"), dpi=100, bbox_inches="tight"
            )
            plt.close()
            # gr.to_excel(_raw_data_folder.joinpath(f'{_name}.xlsx'))

        if "EIS" in EC_merged_dict.keys():
            ECname = "EIS"

            EIS_merged = EC_merged_dict.get("EIS")

            model_select = EIS_merged["PARS"].EIS_Model_EEC.unique()[1]
            print(f"Chosen model: {model_select}")
            #            'Model(EEC_2CPE)'
            _raw_data_folder = mkfolder(
                EC_folder.joinpath(f"{ECname}_{model_select}_raw_data")
            )
            EIS_merg_matches_clean = EIS_merged["AST_matches"].drop_duplicates()
            EIS_merg_matches_clean = EIS_merg_matches_clean.loc[
                (EIS_merg_matches_clean.Sweep_Type == "cathodic")
                & (EIS_merg_matches_clean.RPM_DAC > 1000)
                & (EIS_merg_matches_clean.pH < 2)
            ]
            # EIS_merg_matches_clean.to_excel(_raw_data_folder.parent.joinpath('EIS_AST_matches_pH1_1500rpm.xlsx'))
            _loadDF = read_load_pkl(f"{ECname}_data_raw")
            _mcols = [EvRHE, "Sweep_Type", "SampleID", "postAST"] + [
                "postAST_post",
                "postAST_pre",
            ]
            _grpcols = ("cathodic", 1, 0.6)
            AST_grp_cols = ["postAST_post", "Sweep_Type", "pH", "E_RHE", "AST_row"]
            if _loadDF.empty:
                DF_pltqry = EIS_merged["AST_diff_filter"].groupby(AST_grp_cols)
                _grp_keys = [
                    (key, val)
                    for key, val in DF_pltqry.groups.items()
                    if "cathodic" in key and key[0] is not np.nan
                ]
                #                _st_raw, _st_fit, _st_source = [],[],[]
                _st_comb = []
                #                for pAST in ['postAST_LC','postAST_sHA']:
                for (pAST, Sweep_Type, pH, E_RHE, AST_n), _gr in _grp_keys:
                    _grpcols = (Sweep_Type, pH, E_RHE, AST_n)
                    #                    ASTgrp = DF_pltqry.get_group((pAST,*_grpcols))
                    ASTgrp = DF_pltqry.get_group((pAST, *_grpcols))
                    _samplegrpcols = ["SampleID", "Gas"]
                    for n, gr in ASTgrp.groupby(_samplegrpcols):

                        _prepost_cols = set(
                            [
                                "_".join(i.split("_")[0:-1])
                                for i in list(gr.columns)
                                if i.endswith("_pre") or i.endswith("_post")
                            ]
                        )

                        for status in ["pre", "post"]:
                            gr.EIS_File_SpecRaw_post
                            _setAST = "no" if "pre" in status else pAST
                            _uniq_id = dict(
                                zip(
                                    AST_grp_cols + _samplegrpcols,
                                    (_setAST, *_grpcols, *n),
                                )
                            )
                            _srcfls_raw = [
                                i
                                for i in gr[f"EIS_File_SpecRaw_{status}"].unique()
                                if pd.notna(i)
                            ]
                            _srcfls_fit = [
                                i
                                for i in gr[f"EIS_File_SpecFit_{status}"].unique()
                                if pd.notna(i)
                            ]
                            _srcfls_source = [
                                i
                                for i in gr[f"EIS_sourceFilename_{status}"].unique()
                                if pd.notna(i)
                            ]
                            #                            _srcfls_fit
                            #                            _test = list(_srcfls_fit[0].parent.rglob(f'*{_GPDRT_path.stem}*'))
                            #                            _GP_DRT_succes_rows = gr[f'EIS_GP_DRT_run_success_{status}'].dropna(axis=0)
                            if _srcfls_fit:
                                #                            not _GP_DRT_succes_rows.empty and all([i for i in _GP_DRT_succes_rows.unique()]):
                                #                                _GPDRT_path = Path(gr[f'EIS_GP_DRT_fit_{status}'].unique()[0])
                                def _GP_DRT_load(
                                    _srcfls_fit, DRTsuffix, _uniq_id, status
                                ):
                                    #                                    Path(_srcfls_fit[0]).parent.rglob(f'{Path(_srcfls_fit[0]).name.split("_spectrumfit_")[0]}')
                                    #                                    _test = list(_srcfls_fit[0].parent.rglob(f'*{_GPDRT_path.stem}*'))
                                    #                                    _test = list(_srcfls_raw[0].parent.rglob(f'*{_GPDRT_path.stem}*'))
                                    try:
                                        _DRT_path = list(
                                            Path(_srcfls_fit[0]).parent.rglob(
                                                f'*{Path(_srcfls_fit[0]).name.split("_spectrumfit_")[0]}*{DRTsuffix}*'
                                            )
                                        )
                                        #                                    _GPDRT_path.parent.joinpath(_GPDRT_path.stem+f'_GP_{DRTsuffix}.pkl')
                                        if _DRT_path:
                                            _DRT_DF = pd.read_pickle(_DRT_path[0])
                                            _DRT_DF = _DRT_DF.rename(
                                                columns={
                                                    i: f"{i}_{status}"
                                                    for i in _DRT_DF.columns
                                                }
                                            )
                                        else:
                                            _DRT_DF = pd.DataFrame()
                                    except:
                                        _DRT_DF = pd.DataFrame()
                                    # return _DRT_DF

                                    _DRT_DF = _DRT_DF.assign(**_uniq_id)
                                    return _DRT_DF

                                _DRT_star = _GP_DRT_load(
                                    _srcfls_fit, "DRT_Z_star", _uniq_id, status
                                )
                                _DRT_exp = _GP_DRT_load(
                                    _srcfls_fit, "Z_exp", _uniq_id, status
                                )
                                if not _DRT_star.empty and not _DRT_exp.empty:
                                    _DRT_star = pd.merge_asof(
                                        _DRT_star,
                                        _DRT_exp,
                                        left_on=f"freq_vec_star_{status}",
                                        right_on=f"freq_vec_{status}",
                                        suffixes=["", "_exp"],
                                    )
                                else:
                                    _DRT_star = pd.merge(_DRT_star, _DRT_exp)
                            else:
                                _DRT_star = pd.DataFrame()

                            if _srcfls_fit:
                                _mcols = set(_mcols).union(list(_uniq_id.keys()))
                                _src_data_raw = _read_excel_df(
                                    _srcfls_raw,
                                    _mcols,
                                    _setAST,
                                    pAST,
                                    ref_dict=_uniq_id,
                                )
                                _src_data_fit = _read_excel_df(
                                    _srcfls_fit,
                                    _mcols,
                                    _setAST,
                                    pAST,
                                    ref_dict=_uniq_id,
                                ).query(f'Model_EEC_{status} == "{model_select}"')
                                _src_data_source = _read_excel_df(
                                    _srcfls_source,
                                    _mcols,
                                    _setAST,
                                    pAST,
                                    ref_dict=_uniq_id,
                                ).query(f'Model_EEC_{status} == "{model_select}"')
                                _st_comb.append(
                                    (
                                        _src_data_raw,
                                        _src_data_fit,
                                        _src_data_source,
                                        _DRT_star,
                                    )
                                )
                #                                _test_plot_EIS_spectra(_src_data_fit,_DRT_star, status)
                #                                pd.concat([pd.read_excel(i,  index_col=[0]) for i in _sourcedfls],sort=False,ignore_index=True).dropna(axis=1,how='all')
                #                                _stCVdata = _stCVdata.assign(**{'postAST' :_setAST, 'postAST_post' : pAST })
                #                                _stCVdata = _stCVdata.rename(columns = {i : f'{i}_{status}' for i in _stCVdata.columns if i not in _mcols})
                #                                _st.append(_stCVdata)
                #                            if _srcfls_raw:
                EIS_data_raw = pd.concat([i[0] for i in _st_comb])
                save_DF_pkl(f"{ECname}_{model_select}_data_raw", EIS_data_raw)
                EIS_data_fit = pd.concat([i[1] for i in _st_comb])
                save_DF_pkl(f"{ECname}_{model_select}_data_fit", EIS_data_fit)
                EIS_data_source = pd.concat([i[2] for i in _st_comb])
                save_DF_pkl(f"{ECname}_{model_select}_data_source", EIS_data_source)
                EIS_data_DRT = pd.concat([i[3] for i in _st_comb])
                save_DF_pkl(f"{ECname}_{model_select}_data_DRT", EIS_data_DRT)

            else:
                EIS_data_raw = read_load_pkl(f"{ECname}_{model_select}_data_raw")
                EIS_data_fit = read_load_pkl(f"{ECname}_{model_select}_data_fit")
                EIS_data_source = read_load_pkl(f"{ECname}_{model_select}_data_source")
                EIS_data_DRT = read_load_pkl(f"{ECname}_{model_select}_data_DRT")

            _plt_grp_cols = ["SampleID", "AST_row", "Gas", "E_RHE"]
            EIS_data_DRT_grp = EIS_data_DRT.groupby(_plt_grp_cols)
            for n, gr in EIS_data_fit.groupby(_plt_grp_cols):
                n, gr
                "postAST_post"
                if n in EIS_data_DRT_grp.groups:
                    _DRT_grp = EIS_data_DRT_grp.get_group(n)
                else:
                    _DRT_grp = pd.DataFrame()
                _ndir = "_".join(
                    [
                        str(int(1000 * i)) if "float" in str(type(i)) else str(i)
                        for i in n[0:2]
                    ]
                )
                #                _ndir += f'_{gr.postAST_post.unique()[0]}'
                _png_stem = "_".join(
                    [
                        str(int(1000 * i)) if "float" in str(type(i)) else str(i)
                        for i in n
                    ]
                )
                _test_plot_EIS_spectra(
                    gr, _DRT_grp, save_path=_raw_data_folder.joinpath(_ndir, _png_stem)
                )
                # TODO check for missing status 'pre' DRT plots ...
                # TODO change plotting loop to include post and pre
            # TODO plot Pars vs E_RHE with pre and post combos besides raw data
            _plt_grp_cols_noE = ["SampleID", "AST_row"]
            for n, _pargrp in EIS_merged["AST_diff"].groupby(_plt_grp_cols_noE):
                n, _pargrp
                _ndir = "_".join(
                    [
                        str(int(1000 * i)) if "float" in str(type(i)) else str(i)
                        for i in n
                    ]
                )
                #                _ndir += f'_{_pargrp.postAST_post.unique()[0]}'
                _test_plot_EIS_Pars(_pargrp, save_path=_raw_data_folder.joinpath(_ndir))
        #                'postAST_post'
        #                if n in EIS_data_DRT_grp.groups:
        #                    _DRT_grp = EIS_data_DRT_grp.get_group(n)
        #                else:
        #                    _DRT_grp = pd.DataFrame()
        #                fig,axes = plt.subplots(2,2,figsize=(10,10))
        #                gr.plot(x=EvRHE,y=['Jcorr_pre','Jcorr_post'],ax=axes[1][0])
        #                gr.plot(x=EvRHE,y=['Frac_H2O2_pre','Frac_H2O2_post'],ylim=(0,50),ax=axes[0][0])
        #                gr.plot(x=EvRHE,y=['Jkin_min_pre','Jkin_min_post'],ylim=(0,3),xlim=(0.5,0.9),ax=axes[1][1])
        #                fig.suptitle(f'{n}')
        #                [ax.legend(fontsize=14) for ax1 in axes for ax in ax1]
        #                _name= '_'.join([str(i) for i in (*_grpcols,*n)])
        #                plt.savefig(_raw_data_folder.joinpath(f'{_name}.png'),dpi=100,bbox_inches='tight')
        #                plt.close()
        #                gr.to_excel(_raw_data_folder.joinpath(f'{_name}.xlsx'))
        if "HER" in EC_merged_dict.keys():
            print("export HER raw data")
            # TODO EHR
        if "OER" in EC_merged_dict.keys():
            print("export OER raw data")
            # TODO OER

    #    _src_data_fit, _DRT_star, status_opts = gr,EIS_data_DRT_grp.get_group(n), ['pre','post']
    def _test_plot_EIS_spectra(
        _src_data_fit, _DRT_star, status_opts=["pre", "post"], save_path=""
    ):
        fig, axes = plt.subplots(ncols=1, nrows=4, figsize=(8, 10))
        ax1, ax2, ax3, ax4 = axes
        _status = {}
        for status in status_opts:
            c_exp = "green" if "pre" in status else "orange"
            _src_data_fit.plot(
                x=f"FIT_Zre_{status}", y=f"FIT_-Zim_{status}", kind="line", ax=ax1
            )
            _src_data_fit.plot(
                x=f"DATA_Zre_{status}",
                y=f"DATA_-Zim_{status}",
                kind="scatter",
                c=c_exp,
                ax=ax1,
            )

            _src_data_fit.plot(
                x=f"FIT_Yre_{status}", y=f"FIT_Yim_{status}", kind="line", ax=ax2
            )
            _src_data_fit.plot(
                x=f"DATA_Yre_{status}",
                y=f"DATA_Yim_{status}",
                kind="scatter",
                c=c_exp,
                ax=ax2,
            )

            _status.update(
                {
                    status: _src_data_fit.dropna(
                        subset=[f"DATA_Zre_{status}"]
                    ).postAST.unique()[0]
                }
            )
            if not _DRT_star.empty:
                _DRT_star.plot(
                    x=f"freq_vec_star_{status}",
                    y=f"gamma_vec_star_{status}",
                    logx=True,
                    ax=ax4,
                )
                _DRT_star = _DRT_star.assign(
                    **{
                        f"-1_Z_exp_imag_{status}": -1
                        * _DRT_star[f"Z_exp_imag_{status}"],
                        f"-1_Z_im_vec_star_{status}": -1
                        * _DRT_star[f"Z_im_vec_star_{status}"],
                    }
                )
                _DRT_star.plot(
                    x=f"freq_vec_{status}",
                    y=f"-1_Z_exp_imag_{status}",
                    logx=True,
                    ax=ax3,
                    kind="scatter",
                    c=c_exp,
                )
                _DRT_star.plot(
                    x=f"freq_vec_star_{status}",
                    y=f"-1_Z_im_vec_star_{status}",
                    logx=True,
                    ax=ax3,
                )
                # TODO Add plot for -Zim vs Freq for DRT_fitting and Experimental data points...

        #        _src_data_fit.plot(x=f'FIT_Yre_{status}', y=f'FIT_Yim_{status}',kind='line',ax=axes[1][1])
        #        _src_data_fit.plot(x=f'DATA_Yre_{status}', y=f'DATA_Yim_{status}',kind='scatter', c='r',ax=axes[1][1])
        _titletxt = dict(
            [
                (i, _src_data_fit[i].unique()[0])
                for i in _src_data_fit.columns
                if _src_data_fit[i].nunique() == 1
                and not any([opt in i for opt in status_opts])
                and i in Load_from_Indexes.EC_label_cols[1::]
            ]
        )
        plt.suptitle(f'{", ".join([str(i) for i in _titletxt.values()])},\n {_status}')
        _newdir = f'{save_path.parent.name}_{_status["post"]}'
        mkfolder(save_path.parent.parent.joinpath(_newdir))
        if save_path:
            plt.savefig(
                save_path.parent.parent.joinpath(_newdir, save_path.stem).with_suffix(
                    ".png"
                ),
                dpi=200,
                bbox_inches="tight",
            )
        else:
            plt.show()
        plt.close()

    def _test_plot_EIS_Pars(_pargrp, status_opts=["pre", "post"], save_path=""):
        _pargrp.EIS_lmfit_var_names_pre
        #        ax1, ax2, ax3,ax4 = axes
        _varsgrp = [
            [
                a
                for i in _pargrp[f"EIS_lmfit_var_names_{_st}"].unique()
                for a in i.split(", ")
            ]
            for _st in status_opts
        ]
        _vars = set([a for i in _varsgrp for a in i])
        # mkfolder(save_path)
        for _par in _vars:
            _ps = eisplot(_par)
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 6))
            _status = {}
            _valsmax = []
            for status in status_opts:
                gas_pars = [
                    i
                    for i in _pargrp.columns
                    if i.startswith(f"EIS_{_par}")
                    and i.endswith(status)
                    and not "stderr" in i
                    and not _pargrp[i].isna().all()
                    and len(i.split("_")) == 4
                    and i.split("_")[1] == _par
                ]
                if "_" in _par:
                    gas_pars = [
                        i
                        for i in _pargrp.columns
                        if i.startswith(f"EIS_{_par}")
                        and i.endswith(status)
                        and not "stderr" in i
                        and not _pargrp[i].isna().all()
                        and len(i.split("_")) == 5
                        and "_".join(i.split("_")[1:3]) == _par
                    ]

                _status.update({status: _pargrp[f"postAST_{status}"].unique()[0]})

                for gp in gas_pars:
                    _stderr = [
                        i
                        for i in _pargrp.columns
                        if i.startswith(f'{"_".join(gp.split("_")[0:2])}')
                        and i.endswith(status)
                        and "stderr" in i
                    ]
                    c_exp = "green" if "pre" in status else "orange"
                    ms = "s" if not "guess" in gp else "*"
                    #                    _axnum = 0 if 'N2' in gp else 1
                    #                    print(gp,c_exp,_axnum)
                    _pargrp.plot(
                        x="E_RHE",
                        y=gp,
                        kind="scatter",
                        yerr=_stderr[0],
                        c=c_exp,
                        ax=ax,
                        logy=_ps.logy,
                        label=gp,
                        s=80,
                        marker=ms,
                    )
                    _clvals = [
                        i for i in _pargrp[gp].dropna().values if i > 0 and i < 1e6
                    ]
                    _vals_indx = [
                        n
                        for n, i in enumerate(zscore(_clvals))
                        if np.abs(i) < 3 and _clvals[n] > 0
                    ]
                    if _vals_indx:
                        _valsmax.append(np.array(_clvals)[_vals_indx].max())

            # _vals[_vals_indx].min()
            _ylim = _ps.ylim  # if not _par == 'Cdlp' else (0,0.002)
            if _valsmax:
                _ylim = (_ylim[0], np.min([_ylim[1], np.max(_valsmax).max() * 2]))
            ax.set_ylim(_ylim)
            ax.set_ylabel(_par)
            _titletxt = dict(
                [
                    (i, _pargrp[i].unique()[0])
                    for i in _pargrp.columns
                    if _pargrp[i].nunique() == 1
                    and not any([opt in i for opt in status_opts])
                    and i in Load_from_Indexes.EC_label_cols[1::]
                ]
            )
            plt.suptitle(
                f'{", ".join([str(i) for i in _titletxt.values()])},\n {_status}'
            )

            _newdir = save_path.parent.joinpath(f'{save_path.name}_{_status["post"]}')
            mkfolder(_newdir)
            _stem = f'{_par}_{_titletxt["Gas"]}_{_titletxt["SampleID"]}_{_titletxt["RPM_DAC"]}_{_pargrp.AST_row.unique()[0]}'
            if save_path:
                plt.savefig(
                    _newdir.joinpath(_stem).with_suffix(".png"),
                    dpi=200,
                    bbox_inches="tight",
                )
            else:
                plt.show()
            plt.close()

    def _test_plots():
        # N2CV plot test
        _x = [["IndividualLabel"], "SampleID"][0][0]

        N2_AST_diff = EC_merged_dict.get("N2CV").get("AST_diff_filter")
        ECname = "N2"
        _raw_data_folder = mkfolder(EC_folder.joinpath(f"{ECname}_raw_data"))
        _y = (["N2_Cdl_Fcm-2_600_diff_perc"], (-80, 100))
        _y = (["N2_Cdl_Fcm-2_600_pre", "N2_Cdl_Fcm-2_600_post"], (0, 50e-3))
        N2_pltqry = N2_AST_diff.groupby(
            ["postAST_post", "Sweep_Type", "pH", "Loading_cm2"]
        )
        n = ("postAST_LC", "cathodic", 1.0, 0.379)
        gr = N2_pltqry.get_group(n)
        for n, gr in N2_pltqry:  # .groupby(['postAST_post', 'Sweep_Type','pH']
            n, gr
            _title = "_".join([str(i) for i in n])
            fig, ax = plt.subplots()
            gr.set_index(_x).sort_values(by=_x).plot.bar(
                y=_y[0], ylim=_y[1], rot=0, title=_title, ax=ax
            )
            plt.savefig(
                _raw_data_folder.joinpath(f"bar_Cdl_{_title}.png"),
                dpi=200,
                bbox_inches="tight",
            )
            plt.close()

        N2_AST_diff.groupby(["postAST_post"]).plot(
            x="SampleID", y="N2_Cdl_Fcm-2_400_diff_abs", ylim=(-200, 200)
        )

        tt3 = EC_index.loc[
            (EC_index.PAR_date_day_dt.isin(set([a for i in AST_days for a in i])))
            & (EC_index.PAR_exp == "N2_act")
        ]

        gg = ("postAST_LC", "cathodic", 1.0, 0.379)
        ASTgrp = N2_pltqry.get_group(gg)

        # ORR
        _AST_diff = "AST_diff_filter"

        ORR_AST_diff = EC_merged_dict.get("ORR").get(_AST_diff).reset_index()
        ECname = "ORR"
        list(ORR_AST_diff.columns)
        _raw_data_folder = mkfolder(EC_folder.joinpath(f"{ECname}_raw_data"))
        _y = "ORR_Frac_H2O2_500_diff_abs"
        _y = (["ORR_Frac_H2O2_500_pre", "ORR_Frac_H2O2_500_post"], (0, 40), False)
        _y = (["ORR_Jkin_min_750_pre", "ORR_Jkin_min_750_post"], (0.01, 10), True)
        _y = (["ORR_Jkin_min_650_pre", "ORR_Jkin_min_650_post"], (0.1, 10), True)
        _y = (["ORR_E_half_pre", "ORR_E_half_post"], (0.35, 0.8), False)
        _y = (["ORR_E_onset_pre", "ORR_E_onset_post"], (0.65, 1), False)

        ORR_pltqry = ORR_AST_diff.query("RPM_DAC_uni > 1000").groupby(
            ["postAST_post", "RPM_DAC_uni", "Sweep_Type", "pH", "Loading_cm2"]
        )
        #        ORR_AST.query('RPM_DAC_uni > 1000').groupby(['post_postAST', 'RPM_DAC_uni','Sweep_Type','pH']).groups
        ORR_pltqry.groups
        #        oldORR.groupby(['postAST_post', 'RPM_DAC','Sweep_Type','pH'])

        [
            gr.plot(x=_x, y=_y[0], ylim=_y[1], title=", ".join(str(i) for i in n))
            for n, gr in ORR_pltqry
        ]

        for n, gr in ORR_pltqry:
            if n[1] > 1000 and "mean" in n[2]:
                fig, ax = plt.subplots()
                _title = "_".join([str(i) for i in n])
                gr.sort_values(by="IndividualLabel").plot.bar(
                    x=_x, y=_y[0], ylim=_y[1], logy=_y[2], title=_title, rot=60, ax=ax
                )
                _ylbl = _y[0][0].split("_pre")[0] if "pre" in _y[0][0] else _y[0][0]
                ax.set_ylabel(_ylbl)
                ax.legend(ncol=1, fontsize=10)
                plt.savefig(
                    _raw_data_folder.joinpath(f"bar_{ECname}_{_title}_{_ylbl}.png"),
                    dpi=200,
                    bbox_inches="tight",
                )

        # [gr.sort_values(by='IndividualLabel').plot.bar(x= _x , y=_y[0],ylim=_y[1],logy=_y[2], title=', '.join([str(i) for i in n]), rot=60)
        # for n,gr in ORR_pltqry if n[1] > 1000 and 'mean' in n[2]]

        [
            (
                n,
                "Pre",
                gr["PAR_file_pre"].unique(),
                "Post",
                gr["PAR_file_post"].unique(),
            )
            for n, gr in ORR_pltqry
        ]

        ORR_pltqry.get_group(("postAST_sHA", 1500.0, "cathodic", 1.0))[
            ["SampleID", "ORR_E_half_post"]
        ]
        ORR_pltqry.get_group((None, 1500.0, "cathodic", 1.0))[
            ["SampleID", *_y[0], "PAR_file_post"]
        ]
        ORR_pltqry.get_group(list(ORR_pltqry.groups.keys())[0])

        # ORR KL
        KL_AST_diff = EC_merged_dict.get("KL").get(_AST_diff).reset_index()
        ECname = "ORR"
        _raw_data_folder = mkfolder(EC_folder.joinpath(f"{ECname}_raw_data"))

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
        E_KL = 0.5
        KL_qry = KL_AST_diff.loc[
            np.isclose(KL_AST_diff.ORR_E_AppV_RHE, E_KL, atol=0.01)
        ].groupby(["postAST_post", "Electrode", "Sweep_Type", "pH"])
        KL_AST_diff.query("ORR_E_AppV_RHE ==@E_KL")

        for n, gr in KL_qry:
            if "Disk" in n[1] and "mean" in n[-2]:
                fig, ax = plt.subplots()
                _title = "_".join([str(i) for i in n + (f"E={E_KL}",)])
                _ylbl = _y[0][0].split("_pre")[0] if "pre" in _y[0][0] else _y[0][0]
                ax.set_ylabel(_ylbl)
                ax.legend(ncol=1, fontsize=10)
                #                _title = '_'.join([str(i) for i in n])
                gr.sort_values(by="IndividualLabel").plot.bar(
                    x=_x, y=_y[0], ylim=_y[1], title=_title, rot=60, ax=ax
                )
                plt.savefig(
                    _raw_data_folder.joinpath(f"bar_{ECname}_KL__{_title}_{_ylbl}.png"),
                    dpi=200,
                    bbox_inches="tight",
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
            if "Disk" in n[1]
        ]

        [
            gr.sort_values(by="IndividualLabel").plot.bar(
                x=_x,
                y=_y[0],
                ylim=_y[1],
                title=", ".join([str(i) for i in n + (E_KL,)]),
                rot=60,
            )
            for n, gr in KL_qry
            if "Disk" in n[1] and "mean" in n[2]
        ]
        # EIS
        _y = (["EIS_Rct_O2_pre", "EIS_Rct_O2_post"], (0, 600), (0.55, 0.65))
        EIS_qry = EIS_AST_diff.query(
            'E_RHE > @_y[2][0] & E_RHE < @_y[2][1] & Sweep_Type == "cathodic"'
        ).groupby(["postAST_post", "Sweep_Type", "pH", "E_RHE"])

        [
            gr.sort_values(by="IndividualLabel")
            .dropna(axis=1, how="all")
            .plot.bar(
                x=_x, y=_y[0], ylim=_y[1], title=", ".join([str(i) for i in n]), rot=60
            )
            for n, gr in EIS_qry
            if not gr.dropna(axis=1, how="all").empty
        ]
        # HER test plot
        HER_AST_diff_E = EC_merged_dict["HER_E_slice"]
        _y = (["HER_J_upper_pre", "HER_J_upper_post"], (-2, 2), (-0.49, -0.52))
        _y = (["HER_Tafel_slope_pre", "HER_Tafel_slope_post"], (0, 1e3), (-0.39, -0.42))
        E_350mV_slice = HER_AST_diff_E.loc[(HER_AST_diff_E["HER_Segnum"] > 1)].query(
            '(HER_type == "E_slice") & (HER_at_E_slice < @_y[2][0]) & (HER_at_E_slice > @_y[2][1])'
        )

        HER_pltqry = E_350mV_slice.groupby(
            ["postAST", "Sweep_Type", "pH", "Loading_cm2"]
        )
        HER_pltqry = E_350mV_slice.groupby(
            ["SampleID", "Sweep_Type", "pH", "Loading_cm2"]
        )
        [
            gr.sort_values(by="IndividualLabel").plot.bar(
                x=_x, y=_y[0], ylim=_y[1], title=", ".join([str(i) for i in n]), rot=60
            )
            for n, gr in HER_pltqry
        ]

        E_350mV_slice = HER_AST_diff.loc[(HER_AST_diff["HER_Segment #_pre"] > 1)].query(
            '(HER_type_pre == "E_slice") & (HER_at_E_slice_pre < -0.29) & (HER_at_E_slice_pre > -0.33)'
        )
        jmA2_slice = HER_AST_diff.loc[(HER_AST_diff["Segment #"] > 1)].query(
            '(HER_type == "j_slice_onset") & (HER_at_J_slice == -2)'
        )
        # _Dm.plot(x='E_RHE',y=['EIS_Rct_O2_pre', 'EIS_Rct_O2_post', 'EIS_Rct_O2_diff_abs']) # test plotting


#%%  == NOT_postChar ==


class NOT_postChar:

    suffixes = ["R", "P", "EA", "B"]

    def __init__(self):
        self.folder = FindExpFolder("PorphSiO2").folder
        #        FindExpFolder().TopDir.joinpath(Path('Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc'))
        self.raman = self.postRaman
        self.bet = self.postBET
        self.ea = self.postEA
        self.prec_ea_ratios()
        self.merged = self.merged()

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
        betdir = FindExpFolder("BET").DestDir
        bet_files = betdir.rglob("BET_pars_index*pkl")
        BET_ovv = pd.concat([pd.read_pickle(i) for i in bet_files])
        BET_ovv_template = BET_ovv
        # .loc[BET_ovv.SampleID.isin(PorphSiO2_template().SampleID.unique())]
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


#%%  === postCorrs ===
#        EA_results.loc[EA_results.SampleID.isin(isin_lst)]
def get_postChar_clean():
    chars = postChar()

    _ads_cols = [i for i in chars.bet.columns if "_ads" in i]
    _bet_cols = [
        "SampleID",
        "Metal",
        "color",
        "B_RAW_ads_BET_SA",
        "B_RAW_ads_BET_C_constant",
        "B_RAW_ads_t_Halsey_lin_0_Vmicro",
    ]
    _bet_rpt = [
        "SampleID",
        "B_RAW_ads_B_RPT_BET_Area",
        "B_RAW_ads_B_RPT_tmethod_corr_microporeA",
    ]
    chars.bet[_bet_cols]
    chars.bet[_bet_rpt]


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

    def filter_dups(DF, ECname):

        DF = DF.assign(**{"str_index": [";".join(sorted(i)) for i in DF.index]})
        DF = DF.drop_duplicates(subset=["str_index", "corr_val"])
        DF = DF.drop(columns="str_index")
        DF = postCorrs.exclude_ECname(DF, ECname)
        return DF

    def get_EC_postchar(_pklstem="EC_postchar_allcorrs", return_df=False):
        _pkl_dict = "EC_postchar_allcors_dict"
        if return_df:
            _read_pkl = read_load_pkl(_pklstem)
        AllCorrs = load_dict_pkl(_pkl_dict)

        if not AllCorrs:
            print("Allcorrs reloading...")
            EC_postchars = EC_PorphSiO2.EC_merge_postchar()
            AllCorrs = {
                ECname: {
                    "data": val,
                    "corrs": {
                        grp: postCorrs.filter_dups(dfcorr, ECname)
                        for grp, dfcorr in postCorrs.allchars_corr(val, ECname).items()
                    },
                }
                for ECname, val in EC_postchars.items()
            }
            # AllCorrs = {ECname: {'data' : val['data'], 'corrs' : {grp : postCorrs.filter_dups(dfcorr,ECname) for grp,dfcorr in val['corrs'].items()}} for ECname,val in AllCorrs.items() }
            save_dict_pkl(_pkl_dict, AllCorrs)
            if return_df:
                DF_AC = pd.concat(
                    [
                        v1.assign(**{"EC_type": [(i, *k1)] * len(v1)})
                        for i, val in AllCorrs.items()
                        for k1, v1 in val["corrs"].items()
                    ]
                )
                DF2 = DF_AC.loc[DF_AC.postAST_post != "postAST"]
                save_DF_pkl(_pklstem, DF_AC)
                return AllCorrs, DF_AC
            else:
                return AllCorrs
        else:
            return AllCorrs

    #            Acorrs = postCorrs.allchars_corr(Allchars)
    def test_options():
        _postAST = ["postAST_LC", "postAST", "postAST_sHA"]
        _testopts = {}
        for _p in _postAST:
            _std = (f"{_p}", "cathodic", 1.0, 0.379)
            _std_cols = ("postAST_post", "Sweep_Type", "pH", "Loading_cm2")
            _tests = {
                "KL": {
                    "select": (*_std, "KL_I_Disk", 0.4),
                    "columns": (*_std_cols, "Electrode", "ORR_E_AppV_RHE"),
                    "endswith": ["nElectrons_diff_perc"],
                },
                "ORR": {
                    "select": (*_std, 1500),
                    "columns": (*_std_cols, "RPM_DAC_uni"),
                    "endswith": ["Jkin_min_700_diff_perc"],
                },
                "N2CV": {
                    "select": _std,
                    "columns": _std_cols,
                    "endswith": ["Cdl_700_diff_perc"],
                },
                "EIS": {
                    "select": (*_std, 0.6),
                    "columns": (*_std_cols, "E_RHE"),
                    "endswith": ["_diff_perc"],
                },
                "HER": {
                    "select": (*_std, 0.6),
                    "columns": (*_std_cols, "E_RHE"),
                    "endswith": ["_diff_perc"],
                },
            }
            _testopts.update({_p: _tests})
        return _testopts

    def check_Allcorrs():
        DFAC = DF_AC.copy()
        DFAC = DFAC.reset_index()
        #        DFAC.drop_duplicates(subset=['level_0','level_1','corr_val','chars'])
        DFAC = DFAC.loc[DFAC.postAST_post != "postAST"].drop_duplicates(
            subset=["level_0", "level_1", "corr_val", "chars"]
        )

        for pAST, pgr in DFAC.loc[(DFAC.pH == 1) & (DFAC.Loading_cm2 == 0.379)].groupby(
            "postAST_post"
        ):
            pAST,
            pgr = pgr.sort_values(by="corr_val")
            _opts = postCorrs.test_options()[pAST][ECname]
            ECname = "ORR"
            ECgrp = (pAST,)
            Allchars = AllCorrs.get(ECname).get("data")

            _cORR = pgr.loc[
                (
                    pgr.level_0.str.startswith(ECname)
                    | pgr.level_1.str.startswith(ECname)
                )
                & pgr.Sweep_Type
                == "mean"
            ].set_index(["level_0", "level_1"])
            for (ERPM, charnm), Egr in _cORR.groupby([_opts["columns"][-1], "chars"]):
                ERPM, Egr
                ECname_corr = [
                    i
                    for i in set([a for i in Egr.chars.unique() for a in i])
                    if i in ECname
                ]

                postCorrs.combo_makers(
                    Erg,
                    charnm,
                    "",
                    Allchars,
                    EC_select,
                    separate_dir="_".join([str(i) for i in (ECname_corr, *ECgrp)]),
                    corrval_zero=1,
                )

    def _include_lst():
        _incl = {"N2": ["N2_Cdl_Fcm-2_600"]}
        return _incl

    def test_select_corrs(AllCorrs, swp="cathodic"):

        AllCorrs = postCorrs.get_EC_postchar()
        print("Experiment options:", ", ".join(AllCorrs.keys()))

        corvalnm = "corr_val"
        ECname = "EIS"
        # TODO for N2CV correlations postAST_LC for JOS4 and JOS5 take another pre scans!!
        ECgrp_keys = EC_types_grp().get(ECname)
        print("EC group format:", ECgrp_keys)
        print(
            "Options:",
            "\n".join(
                [str(i) for i in sorted(AllCorrs.get(ECname).get("corrs").keys())]
            ),
        )

        Allchars = AllCorrs.get(ECname).get("data")
        Allchars_grp = Allchars.groupby(ECgrp_keys)
        _AllC_grp_keys = AllCorrs.get(ECname).get("corrs").keys()
        print(f'{"____".join([str(i) for i in _AllC_grp_keys])}')
        ECgrp = ("postAST_sHA", "mean", 1.0, 0.379, 0.6)
        if not ECgrp in _AllC_grp_keys:
            print("ECgrp not in keys")
            ECgrp = ("postAST_sHA", "mean", 1.0, 0.379, 1500)
            ECgrp = ("postAST_LC", "mean", 1.0, 0.379, 1500)
            ECgrp = ("postAST_LC", "cathodic", 1.0, 0.379)
            ECgrp = ("postAST_sHA", "cathodic", 1.0, 0.379, 0.65)

        _AllC_grp_keys = [i for i in _AllC_grp_keys if i[0] != "postAST"]

        _ECgr_raw_data = Allchars_grp.get_group(_AllC_grp_keys[-5])
        if "ORR" in ECname:
            _AllC_grp_keys = [
                i for i in _AllC_grp_keys if (i[-1] > 1000) and i[0] != "postAST"
            ]
        if "EIS" in ECname:
            _AllC_grp_keys = [
                i for i in _AllC_grp_keys if i[1] == "cathodic" and i[0] != "postAST"
            ]
            _extra_dir = AllCorrs["EIS"]["data"].EIS_Model_EEC_pre.unique()[0]

        def get_ECname_corrs(df, ECname):
            ECname_corr_lst = [
                i
                for i in df.chars.unique()
                if (any([a in ECname for a in i]) and not "RPM" in i)
            ]
            ECname_corr = Counter([a for b in ECname_corr_lst for a in b]).most_common(
                1
            )[0][0]
            return ECname_corr

        def get_chars_ECname_Struct(chE, ECname):
            ECname_corr = get_ECname_corrs(chE, ECname)
            chars_EC_struct = chE.loc[
                chE.chars.isin(
                    [
                        i
                        for i in chE.chars.unique()
                        if any([c == ECname_corr for c in i])
                        and not i == (ECname_corr, ECname_corr)
                    ]
                )
            ].sort_values(corvalnm)
            return chars_EC_struct

        def get_allcorrs_include(
            AllCorrs,
            _AllC_grp_keys,
            ECname,
            ECgrp_keys,
            _extra_dir,
            include=[],
            pre_post=["pre", "post"],
            gases=["O2", "N2"],
        ):
            _collect = []

            _incl = [f"{i}_{g}_{p}" for i in include for p in pre_post for g in gases]
            for ECgrp in _AllC_grp_keys:

                chE_raw = AllCorrs.get(ECname).get("corrs").get(ECgrp, "Empty dict")
                ECname_corr = get_ECname_corrs(chE_raw, ECname)
                chars_EC_struct = get_chars_ECname_Struct(chE_raw, ECname_corr)
                chEA_excl = postCorrs.select(
                    chars_EC_struct, _incl, exclude=["_lintangent_", "_E_dc_RHE_mV"]
                )
                # chEA_excl = chEA_excl.loc[np.abs(chEA_excl.corr_val) < 0.99]
                _collect.append(chEA_excl)
            corrs_include = pd.concat(i for i in _collect)
            corrs_include[f"abs_{corvalnm}"] = np.abs(corrs_include[corvalnm])

            for ECgrp, _include_corrs in corrs_include.groupby(ECgrp_keys):
                _best_idx = (
                    _include_corrs.groupby(level=1)[f"abs_{corvalnm}"]
                    .agg(["mean", "std"])
                    .sort_values("mean", ascending=False)
                    .head(100)
                    .index
                )
                EC_select = dict(zip(ECgrp_keys, ECgrp))
                _sepdir = "_".join([str(i) for i in (ECname_corr, *ECgrp, _extra_dir)])
                postCorrs.combo_makers(
                    _include_corrs,
                    ECname_corr,
                    "",
                    Allchars,
                    EC_select,
                    include=_best_idx,
                    separate_dir=_sepdir,
                    corrval_zero=True,
                    sep_topdir="Allcorrs_best/",
                    force_replot=False,
                    corvalnm=corvalnm,
                    ECgrp=ECgrp,
                )

        get_allcorrs_include(
            AllCorrs,
            _AllC_grp_keys,
            ECname,
            ECgrp_keys,
            _extra_dir,
            include=[],
            pre_post=["pre", "post"],
            gases=["O2", "N2"],
        )
        for ECgrp in _AllC_grp_keys:
            EC_select = dict(zip(ECgrp_keys, ECgrp))
            print(EC_select)
            chE_raw = AllCorrs.get(ECname).get("corrs").get(ECgrp, "Empty dict")
            chE = postCorrs.exclude_ECname(chE_raw, ECname)
            #            Allchars = AllCorrs.get(ECname).get('data')
            chE.chars.unique()
            ECname_corr = get_ECname_corrs(chE, ECname)
            # charsECname = chE.loc[chE.chars.isin([i for i in chE.chars.unique() if any(c in ECname_corr for c in i)])].sort_values(corvalnm)
            chars_EC_struct = chE.loc[
                chE.chars.isin(
                    [
                        i
                        for i in chE.chars.unique()
                        if any([c == ECname_corr for c in i])
                        and not i == (ECname_corr, ECname_corr)
                    ]
                )
            ].sort_values(corvalnm)
            #        chC = postCorrs.select(chE,postChar.suffixes,endswith_set = ['diff_abs'] ).sort_values(corvalnm)
            chC = chars_EC_struct.drop_duplicates()
            #            postCorrs.select(chE,postChar.suffixes,endswith_set = [''] ).sort_values(corvalnm)
            set([i[0] for i in chC.index.unique()])
            postCorrs.combo_makers(
                chC,
                ECname_corr,
                "",
                Allchars,
                EC_select,
                include="R_ion",
                separate_dir="_".join(
                    [str(i) for i in (ECname_corr, *ECgrp, _extra_dir)]
                ),
                corrval_zero=True,
                sep_topdir="R_ion",
                force_replot=False,
                corvalnm=corvalnm,
                ECgrp="",
            )

            ECgrp_destdir = "_".join(
                [str(i) for i in (ECname_corr, *ECgrp, _extra_dir)]
            )
            comboCH = {}
            for Charnm in ["P", "EA", "B", "R"]:
                # EC_corr_char = [i for i in ECname_corr_lst if Charnm in i][0]
                _out = postCorrs.combo_makers(
                    chC,
                    ECname_corr,
                    Charnm,
                    Allchars,
                    EC_select,
                    include="",
                    separate_dir=ECgrp_destdir,
                    corrval_zero=1,
                    force_replot=False,
                    endswith="",
                    corvalnm=corvalnm,
                    ECgrp=ECgrp,
                )
                _out = _out.drop_duplicates()
                _topidx = _out.query("corr_val > 0.75 or corr_val < 0.75")
                _counter = Counter([a for i in _topidx.index for a in i]).most_common(
                    30
                )
                _ordering = [
                    (
                        *_c,
                        _out.loc[
                            [i for i in _topidx.index if any([c in _c[0] for c in i])]
                        ]
                        .corr_val.abs()
                        .sum()
                        .round(3),
                    )
                    for _c in _counter
                ]
                _ordering = sorted(_ordering, key=lambda x: x[-1], reverse=True)
                comboCH.update(
                    {
                        Charnm: {
                            "combos": _out,
                            "counter": _counter,
                            "ordered": _ordering,
                        }
                    }
                )

            Corrmeans = postCorrs.count_common_corrs(comboCH)
            postCorrs.plot_corrmeans(Corrmeans, _mincorr=0.5, dest_folder=ECgrp_destdir)

        #        [(key,val['ordered'][0]) for key,val in comboCH.items()]
        _topidx = (
            comboCH["R"]
            .get("combos", pd.DataFrame())
            .query("corr_val > 0.8 or corr_val < 0.8")
            .index
        )

    #        chC.loc[chC.chars != ('EIS','EIS')]
    #        select(chE,endswith_set = 'diff_abs' )

    def count_common_corrs(comboCH):
        #        chCr = comboCH['R'].get('combos',pd.DataFrame())
        _ends_ovv_lst = []

        for char, CHdict in comboCH.items():
            chCr = CHdict.get("combos", pd.DataFrame())
            _lvl = [n for n, i in enumerate(chCr.index[0]) if i.startswith(char)]

            for _ends in ["_perc", "_pre", "_post", "_abs"]:
                _end_fltr = [i for i in chCr.index if any(c.endswith(_ends) for c in i)]
                if _end_fltr:
                    chC_ends = chCr.loc[_end_fltr]

                    chC_ends = postCorrs.common_par_cols(chC_ends, _ends, char)
                    for parnm, pargrp in chC_ends.groupby("par_nm"):

                        _corrmeans = sorted(
                            [
                                (
                                    _ends,
                                    n,
                                    [
                                        round(
                                            gr.loc[
                                                [
                                                    i
                                                    for i in gr.index
                                                    if any(c.endswith(_ends) for c in i)
                                                ]
                                            ]
                                            .corr_val.abs()
                                            .mean(),
                                            4,
                                        ),
                                        round(
                                            gr.loc[
                                                [
                                                    i
                                                    for i in gr.index
                                                    if any(c.endswith(_ends) for c in i)
                                                ]
                                            ]
                                            .corr_val.abs()
                                            .std(),
                                            4,
                                        ),
                                        gr.loc[
                                            [
                                                i
                                                for i in gr.index
                                                if any(c.endswith(_ends) for c in i)
                                            ]
                                        ].index.unique(),
                                    ],
                                )
                                for n, gr in pargrp.groupby(level=_lvl)
                            ],
                            key=lambda x: x[1][0],
                            reverse=True,
                        )
                        _ovv_corrmeans = [
                            {
                                "charnm": char,
                                "parnm": parnm,
                                "ends": i[0],
                                "char": i[1],
                                "corr_mean": i[2][0],
                                "corr_mean_std": i[2][1],
                            }
                            for i in _corrmeans
                        ]
                        _ends_ovv_lst.append(_ovv_corrmeans)
        Corrmeans = pd.concat(
            [pd.DataFrame(i) for i in _ends_ovv_lst], ignore_index=True, sort=False
        )

        _fltr = ["R_nfev"]
        Corrmeans = Corrmeans.loc[
            Corrmeans.char.isin(
                [
                    i
                    for i in Corrmeans.char
                    if not any([i.startswith(fl) for fl in _fltr])
                ]
            )
        ]
        Corrmeans = Corrmeans.sort_values(by="corr_mean", ascending=False)
        return Corrmeans

    def plot_corrmeans(Corrmeans, _mincorr=0.5, dest_folder=""):
        if dest_folder:
            _dest = mkfolder(plotsfolder.joinpath("corr_means", dest_folder))
        else:
            dest = mkfolder(plotsfolder.joinpath("corr_means"))

        _cmap = dict(
            zip(["P", "EA", "B", "R"], ["blue", "orange", "steelblue", "green"])
        )
        for (n, npar), gr in Corrmeans.loc[Corrmeans.corr_mean > _mincorr].groupby(
            ["ends", "parnm"]
        ):
            try:
                fig, ax = plt.subplots(figsize=(20, 10))
                gr.plot.barh(
                    x="char",
                    y="corr_mean",
                    yerr="corr_mean_std",
                    color=gr.charnm.replace(_cmap),
                    rot=0,
                    title=n,
                    ax=ax,
                    fontsize=10,
                )

                plt.savefig(_dest.joinpath(f"{npar}{n}.png"), bbox_inches="tight")
            #            plt.show()
            except Exception as e:
                print(f"plot corrmeans fail {n}, {e}")
            finally:
                plt.close()

    def common_par_cols(chC_ends, _ends, char):
        _lst = []
        _p_lvl = [n for n, i in enumerate(chC_ends.index[0]) if not i.startswith(char)]
        for i in chC_ends.index:
            _par = i[_p_lvl[0]].strip(
                f"diff_{_ends}" if _ends in ["_abs", "_perc"] else _ends
            )
            try:
                _int = int(_par.split("_")[-1])
                _par = "_".join(_par.split("_")[:-1])

            except:
                pass
            _lst.append(_par)
        chC_ends = chC_ends.assign(**{"par_nm": _lst})
        return chC_ends

    def exclude_ECname(chE, ECname):
        _excl = {
            "ORR": [
                "ORR_lin_hor",
                "ORR_hor_overall",
                "ORR_Ring_E_programmed_modulation",
                "ORR_Ring_RPM",
            ],
            "Date": ["Date_PAR", "B_BET_index"],
            "N2CV": ["N2_delta_mtime", "N2_index", "N2_RHE_fn"],
            "EIS": [
                "EIS_index",
                "EIS_Segment #",
                "EIS_Current Range_",
                "EIS_delta_mtime",
                "EIS_I(A)",
            ],
        }
        _excl_lst = [i for val in _excl.values() for i in val]
        _selected = [
            i
            for i in chE.index
            if not any([c.startswith(_st) for _st in _excl_lst for c in i])
        ]
        #        selected = [i for i in _selected if not any([c.startswith(_st) for _st in _excl.get(Date,[]) for c in i])]
        return chE.loc[_selected]

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
        corvalnm="",
        ECgrp="",
    ):
        if ECgrp:
            chEA = postCorrs.uniq(AllCan, (main_ch, second_ch), ECgrp)
        else:
            chEA = postCorrs.uniq(AllCan, (main_ch, second_ch))
        chEA_excl = postCorrs.select(chEA, second_ch, endswith_set=endswith)

        if type(exclude_second) == type(""):
            if exclude_second:
                exclude_second = [exclude_second]
            else:
                exclude_second = []
        if type(exclude_second) == type([]):
            if exclude_second:
                for excl in exclude_second:
                    chEA_excl = postCorrs.select(chEA_excl, second_ch, exclude=excl)
        _incl_str = ""
        if type(include) == type(""):
            if include:
                include = [include]
            else:
                include = []
        if type(include) == type([]):
            if include:
                _incl_str = "+".join(str(i) for i in include)
                for incl in include:
                    chEA_excl = postCorrs.select(chEA_excl, [incl])
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
            if sep_topdir:
                combo_dest_dir = mkfolder(plotsfolder.joinpath(f"corrs_{sep_topdir}"))
            if separate_dir:
                combo_dest_dir = mkfolder(combo_dest_dir.joinpath(separate_dir))

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

    def chars_filter_cols(Allchars):

        only_floats = [key for key in Allchars.columns if Allchars[key].dtypes != "O"]
        #        [key for key,val in Allchars.dtypes.to_dict().items() if 'float64' in str(val)]

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
            "Delta_PAR_LastMod_pre",
            "filesize_post",
            "Delta_PAR_LastMod_post",
            "filesize_pre",
            "color",
        ]
        filter_cols_lst += [
            i for i in only_floats if i.startswith("B_") and "_RPT" in i
        ]
        _BET_endswith = [
            "_microArea",
            "_Vmicro",
            "_r",
            "mass_RAW",
            "_RPT_Correlation_Coefficient",
            "TEMPERATURE",
            "_error(%)",
            "_RAW_calc_relPrange_max",
            "_RAW_calc_relPrange_min",
            "Station#",
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("B_") and any(i.endswith(c) for c in _BET_endswith)
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("B_") and any(c in i for c in ["_Export_path_"])
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("B_t") and any(i.endswith(c) for c in ["_slope"])
        ]
        _RAMAN_endswith = ["_amplitude", "_sigma", "_bic", "_aic", "_chisqr", "_redchi"]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("R_") and any(i.endswith(c) for c in _RAMAN_endswith)
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
            "test_low",
            "test_sum",
            "test_high",
            "lmfit_",
            "Scanrate",
            "deltaDT",
            "index",
            "Current Range",
            "E_programmed_modulation",
            "RHE_fn",
            "ActionId",
            "Type_action",
            "Scanrate",
            "BRUTE",
            "FINAL",
            "level_0",
            "delta_mtime",
            "I(A)",
            "_lintangent_",
            "_E_dc_RHE_mV",
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("EIS_") and any([c in i[4:] for c in _EIS_filters_list])
        ]

        _ORR_excl = {
            "ORR": [
                "ORR_lin_hor",
                "ORR_hor_overall",
                "ORR_Ring_E_programmed_modulation",
                "ORR_Ring_RPM",
                "ORR_TF_TSerr",
                "ORR_Ring_Current Range",
                "ORR_Ring_Segment",
                "ORR_Scan Rate",
                "ORR_Ring_ActionId",
                "ORR_Ring_Scan",
                "ORR_Ring_E_programmed_modulation",
                "ORR_pH_index",
                "ORR_Current Range",
                "ORR_Instrument_SN",
                "ORR_AC Amplitude",
                "ORR_delta_mtime",
                "ORR_Segment #",
                "ORR_Ring_Status",
                "ORR_Ring_Instrument_SN",
                "ORR_ActionId",
            ],
            "N2": ["N2_delta_mtime"],
            "R_": ["R_nfev"],
        }
        filter_cols_lst += [
            i
            for i in only_floats
            if any([i.startswith(_st) for _st in _ORR_excl.get("ORR", [])])
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if any([i.startswith(_st) for _st in _ORR_excl.get("N2", [])])
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if any([i.startswith(_st) for _st in _ORR_excl.get("R_", [])])
        ]

        _strtswith = [
            "AST_",
            "Date_PAR",
            "_",
            "Delta_PAR_LastMod",
            "dupli_num",
            "filesize",
            "error",
            "Instrument",
            "relPrange",
            "rvalue",
            "wrong",
            "fitR",
        ]
        filter_cols_lst += [
            i for i in only_floats if any([i.startswith(_st) for _st in _strtswith])
        ]
        return filter_cols_lst

    #            += [i for i in only_floats if i.endswith('diff_abs') and i.endswith('diff_perc')]
    #        filter_cols_lst += [i for i in only_floats if i.startswith('B_') and i.endswith()]
    #        filter_cols_lst += [i for i in only_floats if i.startswith('B_') and i.endswith()]
    #        filtered_cols = [i for i in only_floats if i not any(t in i  for t in filter_cols_lst)]
    #        clean_cols = [i for i in only_floats if i not in filter_cols_lst and i in charDF.columns]
    #        return clean_cols

    def allchars_corr(
        Allchars,
        ECname,
        corr_methods=["pearson", "spearman"],
        max_corrval=1.1,
        min_corrval=0,
    ):
        AllCorrs = {}
        #        ECname = 'N2CV' #HER_E_slice'
        #        Allchars = EC_postchars.get(ECname)
        _ECgrps = EC_types_grp().get(
            [i for i in EC_types_grp().keys() if i in ECname][0]
        )
        #        corvalnm = f'corr_val_{swp}'
        corvalnm = "corr_val"

        print(f"allchars_corr start: {ECname}, groupby groups:\n {_ECgrps}")
        if "Sweep_Type" not in Allchars.columns:
            Allchars["Sweep_Type"] = "none"

        filter_cols_lst = postCorrs.chars_filter_cols(Allchars)

        for ngrp, swgrp in Allchars.groupby(_ECgrps):

            #            swp = [i for i in ngrp if i in ['cathodic','anodic','mean']][0]
            swgrp = swgrp.dropna(axis=1, how="all")
            if len(swgrp) > 2:
                clean_cols = [i for i in swgrp.columns if i not in filter_cols_lst]
                clean_cols = [i for i in clean_cols if swgrp[i].nunique() != 1]
                print(ngrp, f"len{len(swgrp)}")
                _allcorm = []
                for corr_method in corr_methods:
                    # allcorrstk = swgrp[clean_cols].corr(method=corr_method, min_periods=3).stack()
                    allcorr_m = swgrp[clean_cols].corr(
                        method=corr_method, min_periods=3
                    )
                    allcorr_m = allcorr_m.mask(
                        np.tril(np.ones(allcorr_m.shape)).astype(np.bool)
                    )
                    allcorrstk = allcorr_m.stack()
                    allcorrs = allcorrstk[
                        (min_corrval < np.abs(allcorrstk))
                        & (np.abs(allcorrstk) < max_corrval)
                    ]
                    if not allcorrs.empty:
                        allcorrs.name = corvalnm
                        # _tpls = [tuple(i) for i in allcorrs.index.values]
                        # _uniqtpls = [i for i in list(set([tuple(i) for i in [set(i) for i in _tpls]])) if len(i) == 2]
                        # allcorrs.iloc[pd.MultiIndex.from_tuples(_uniqtpls)]
                        allcorrs_DF = pd.concat(
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
                        allcorrs_DF["corr_method"] = corr_method
                        #                allcorrs['corr_method'] = corr_method
                        # _{corr_method}'
                        _allcorm.append(allcorrs_DF)
                #            allcorrstk = swgrp[clean_cols].corr(method='pearson', min_periods=3).stack()
                #            allcorrs = allcorrstk[np.abs(allcorrstk) < 1]
                #            allcorrs.name = corvalnm
                if _allcorm:
                    allchrrs = pd.concat(_allcorm)
                    allchrrs = allchrrs.assign(**dict(zip(_ECgrps, ngrp)))
                    AllCorrs.update({ngrp: allchrrs})
        #            allchrrs = allchrrs[np.abs(allchrrs[allcorrs.name]) < 0.991].drop_duplicates(allcorrs.name).sort_values(allcorrs.name)
        #            allchflrtd = postCorrs.remove_corrdups(allchrrs, ECname)
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
            if type(exclude) == str:
                _excl = list(exclude)
            elif type(exclude) == list:
                _excl = exclude
            else:
                _excl = list(exclude)
            idx_select = [
                i for i in idx_select if not any(ex in a for a in i for ex in _excl)
            ]
        #        .sort_values(by=[i for i in chC.columns if 'corr_val' in i][0])
        if len(idx_select) == len(chC):
            return chC
        else:
            return chC.loc[idx_select]

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


#%%  == testing_corr ==


class testing_corr:
    def __init__(DF):
        self.DF = DF

    def plot_bokeh(self):
        from bokeh.plotting import figure, show, output_notebook
        from bokeh.models import (
            ColumnDataSource,
            LinearColorMapper,
            ColorBar,
            BasicTicker,
        )
        from bokeh.transform import transform
        from bokeh.models.formatters import PrintfTickFormatter

        output_notebook()
        df = bet
        colors = ["#d7191c", "#fdae61", "#ffffbf", "#a6d96a", "#1a9641"]
        TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"
        data = df.corr().stack().rename("value").reset_index()
        p = figure(
            x_range=list(df.columns),
            y_range=list(df.SampleID),
            tools=TOOLS,
            toolbar_location="below",
            tooltips=[("Row, Column", "@level_0 x @level_1"), ("value", "@value")],
            height=500,
            width=500,
        )

        p.rect(
            x="level_1",
            y="level_0",
            width=1,
            height=1,
            source=data,
            fill_color={
                "field": "value",
                "transform": LinearColorMapper(
                    palette=colors, low=data.value.min(), high=data.value.max()
                ),
            },
            line_color=None,
        )
        color_bar = ColorBar(
            color_mapper=LinearColorMapper(
                palette=colors, low=data.value.min(), high=data.value.max()
            ),
            major_label_text_font_size="7px",
            ticker=BasicTicker(desired_num_ticks=len(colors)),
            formatter=PrintfTickFormatter(format="%f"),
            label_standoff=6,
            border_line_color=None,
            location=(0, 0),
        )
        p.add_layout(color_bar, "right")

        show(p)


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


def make_pars():
    pars = postChar().merged
    charname = "BET"


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


def filter_char_pars(pars):
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
    ]
    _keepcols = [i for i in pars.columns if not any(i.endswith(a) for a in _endsw)]
    _filteredcols = [i for i in pars.columns if i not in _keepcols]
    return pars[_keepcols]


def prep_corrs(pars, filter_min=0.5):
    only_floats = [
        key for key, val in pars.dtypes.to_dict().items() if "float64" in str(val)
    ]
    parsfltr = filter_char_pars(pars)
    _corrs = parsfltr.dropna(thresh=3, axis=1).corr().dropna(axis=0, how="all")
    _crrstk = _corrs.stack()
    _crrstk = _crrstk.loc[_crrstk.abs() > filter_min]
    _filter = [i for i in _crrstk.index if str(i[0])[0] != str(i[1])[0]]
    _corrs = _crrstk.loc[_filter].unstack()

    return _corrs


def plot_heatmap(pars, charname):

    _corrs = prep_corrs(pars)

    fig, ax = plt.subplots(figsize=(24, 19))

    # plt.matshow(_corrs)
    # plt.matshow(_corrs, fignum=f.number)
    mask = np.triu(np.ones_like(_corrs, dtype=bool))
    # plt.xticks(range(pars.select_dtypes(['number']).shape[1]), pars.select_dtypes(['number']).columns, fontsize=14, rotation=45)
    # plt.yticks(range(pars.select_dtypes(['number']).shape[1]), pars.select_dtypes(['number']).columns, fontsize=14)

    # plt.yticks(np.arange(0.5, len(_corrs.index), 1), _corrs.index)
    # plt.xticks(np.arange(0.5, len(_corrs.columns), 1), _corrs.columns)
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    #
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


def find_combos(_corrs):
    _corrstk = _corrs.stack()
    _highcorr = _corrstk.loc[_corrstk.abs() > 0.9]
    _highcorr.index
    len([set(i) for i in _highcorr.index])
    _setdex = set([(", ".join(set(i)), _highcorr.loc[i]) for i in _highcorr.index])
    [i for i in _setdex]


def pair_plot(pars):
    _corrs = prep_corrs(pars)
    _corrsdex = _corrs.index
    ycol = ["B_RAW_ads_BET_C_constant", "B_RAW_ads_BET_SA"][0]
    # _corrs[ycol]
    ycorrs = _corrs[ycol].dropna().sort_values()
    ycorrsdex = ycorrs.index
    # attend = sns.load_dataset("attention").query("subject <= 12")
    fig, axes = plt.subplots(nrows=len(ycorrsdex), sharex=True)
    for ax, y in zip(axes, ycorrsdex):
        sns.regplot(x=ycol, y=y, data=pars, ax=ax)

    xcols = list(ycorrsdex[0:5]) + list(ycorrsdex[-5:])

    fig, ax = plt.subplots(figsize=(24, 19))
    sns.pairplot(
        pars[xcols + [ycol, "SampleID"]], hue="SampleID", diag_kind="none", kind="reg"
    )
    plt.savefig(
        FindExpFolder("PorphSiO2").compare.joinpath(f"pairplot_{charname}_{ycol}.png"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.close()

    fig, ax = plt.subplots()
    sns.pairplot(pars, x_vars=[ycol], y_vars=xcols, hue="SampleID", kind="reg")
    plt.savefig(
        FindExpFolder("PorphSiO2").compare.joinpath(
            f"pairplot_xvar_{charname}_{ycol}.png"
        ),
        bbox_inches="tight",
        dpi=300,
    )
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


def _set_globals():
    try:
        ParsColl.keys()
        # global ParsColl
    except Exception as e:
        print(f"Load error: {e}")
        if not "Pars_Collection" in globals():
            Pars_Collection = CollectLoadPars(load_type="fast")
            # globals()['Pars_Collection'] = Pars_Collection
            ParsColl = Pars_Collection.pars_collection
            globals()["ParsColl"] = ParsColl


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
