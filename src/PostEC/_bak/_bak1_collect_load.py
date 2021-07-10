# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:30:31 2020

@author: User
"""
#%%

import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

# import os
from platform import system
import glob
import cycler
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from bs4 import BeautifulSoup
import re
from scipy.stats import linregress

# from sklearn import linear_model
import scipy.signal
import itertools
from itertools import chain, repeat
import logging
import datetime as dt


from pathlib import Path

# import h5py
from multiprocessing import Pool, cpu_count

# import timeit
# import time
import datetime as dt

matplotlib.rcParams.update({"font.size": 16})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Helvetica"
plt.rcParams["axes.edgecolor"] = "#333F4B"
plt.rcParams["xtick.color"] = "#333F4B"
plt.rcParams["ytick.color"] = "#333F4B"
try:
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    import seaborn as sns
except Exception as e:
    print("No modules: %s" % e)
##from FolderOrganizing import FileHelper
##from FolderOrganizing import RunEC_classifier
# from FolderOrganizing import PAR_EIS_fit_V2
# import run_PAR_DW
# PostOrganizeFolders(OnlyRecentMissingOVV)
# def start_logger(logdir = Path.cwd() , logger_option_set = True):
#    if logger_option_set == True:
#        # Gets or creates a logger
#        logger = logging.getLogger(__name__)
#        # set log level
#        logger.setLevel(logging.WARNING)
#        # define file handler and set formatter
#        file_handler = logging.FileHandler(logdir.joinpath('PostEC_logger.log'))
#        formatter    = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s')
#        file_handler.setFormatter(formatter)
#        # add file handler to logger
#        logger.addHandler(file_handler)
#        logger.warning('=========== Started logging PAR PostEC...  =============')
#    else:
#        print('=====x====== No logging PAR PostEC...  ======x=======')
#    return logger
print("Name prepare input:", __name__)
# if __name__ == '__main__':
##    from RunEC_classifier import EC_classifier_multi_core
#    sys.path.append(str(Path(__file__).parent.parent.parent))
#    import FileHelper
#    from FindFolders import FindExpFolder
#    from FileFunctions import FileOperations
#    from FindSampleID import GetSampleID
#    from runEC.EC_logging_config import start_logging
#    logger = start_logging(__name__)
#
# elif 'prepare_input' in __name__:
#    print(f'Name: {__name__} for file {__file__}')
#    sys.path.append(str(Path(__file__).parent))
#    from RunEC_classifier import EC_classifier_multi_core
##    import RunEC_classifier
#    sys.path.append(str(Path(__file__).parent.parent.parent))
##    import FileHelper
#    from FindFolders import FindExpFolder
#    from FindSampleID import GetSampleID
#    from FileFunctions import FileOperations
#
if __name__ == "__main__":
    print(f"Package: {__package__}, File: {__file__}")
    FH_path = Path(__file__).parent.parent.parent
    sys.path.append(str(FH_path))
    #    sys.path.append(str(Path(__file__).parent.parent.joinpath('indexer')))
    #    sys.path.append(str(Path(__file__).parent.parent.parent))
    #    sys.path.append("..")
    #    print(sys.path)
    from FileHelper.FindFolders import FindExpFolder
    from FileHelper.FileFunctions import FileOperations
    from FileHelper.FindSampleID import GetSampleID

    from PostChar import SampleSelection, Characterization_TypeSetting, SampleCodesChar
    from ECpy.main_run_PAR_DW import ECRunOVV

    #    from ECpy.experiments.EIS.helper import *
    import post_helper
    import merger

    #    import EC
    #    sys.path.append(list(FH_path.rglob('*.py')))
    #    import FH_path.joinpath('FindExpFolder.py')
    #    import FindExpFolder.py
    #    from FileHelper import FindExpFolder
    #    from FindExpFolder import *
    #    from .experiments import EIS
    #    from .runEC import run_PAR_DW
    from ECpy.runEC.EC_logging_config import start_logging

    logger = start_logging(__name__)

else:
    #    print('\n\n***** run_PAR_DW *****')
    print(f"File: {__file__}, Name:{__name__}, Package:{__package__}")
    #    FH_path = Path(__file__).parent.parent.parent
    #    sys.path.append(str(FH_path))
    #    import FileHelper
    from FileHelper.FindFolders import FindExpFolder
    from FileHelper.FileFunctions import FileOperations
    from FileHelper.FindSampleID import GetSampleID
    from FileHelper.PostChar import (
        SampleSelection,
        Characterization_TypeSetting,
        SampleCodesChar,
    )
    from ECpy.main_run_PAR_DW import ECRunOVV
    from ECpy.runEC.EC_logging_config import start_logging
    from ECpy.PostEC import post_helper, merger

    logger = start_logging(__name__)

#%%
EvRHE = "E_AppV_RHE"


class PostEC:
    AllColls = [
        "Unnamed: 0",
        "Segment #",
        "Point #",
        "E(V)",
        "I(A)",
        "Elapsed Time(s)",
        "Current Range",
        "Status",
        "E Applied(V)",
        "Frequency(Hz)",
        "Z Real",
        "Z Imag",
        "ActionId",
        "AC Amplitude",
        "RHE_OCP",
        "E_AppV_RHE",
        "E_Applied_VRHE",
        "j A/cm2",
        "jmAcm-2",
        "jcorr",
        "Gas",
        "EXP",
        "Electrode",
        "j_ring",
        "RPM",
        "Comment",
        "Measured_OCP",
        "pH",
        "Electrolyte",
        "ScanRate_calc",
        "SampleID",
        "File",
        "BaseName",
        "hash",
        "Instrument",
        "DATE",
        "EvRHE_diff",
        "DestFile",
        "Sweep_Type",
        "Type",
        "Cycle",
        "DAC_V",
        "Scanrate",
        "ORR_scan",
        "Jcorr",
        "J_N2_scan",
        "J_O2_diff",
        "J_O2_diff_diff",
        "Analysis_date",
        "J_2nd_diff",
        "Jkin_max",
        "Jkin_min",
        "E_onset",
        "Diff_lim",
        "E_half",
        "I(A)_ring",
        "I(A)_disk",
        "Frac_H2O2",
        "J_ring",
        "n_ORR",
    ]
    DropColls = [
        "Unnamed: 0",
        "Segment #",
        "Point #",
        "E(V)",
        "I(A)",
        "Elapsed Time(s)",
        "Current Range",
        "Status",
        "E Applied(V)",
        "Frequency(Hz)",
        "Z Real",
        "Z Imag",
        "ActionId",
        "AC Amplitude",
        "RHE_OCP",
        "E_AppV_RHE",
        "jmAcm-2",
        "jcorr",
        "Gas",
        "EXP",
        "Electrode",
        "j_ring",
        "RPM",
        "Comment",
        "Measured_OCP",
        "pH",
        "Electrolyte",
        "ScanRate_calc",
        "SampleID",
        "File",
        "BaseName",
        "hash",
        "Instrument",
        "DATE",
        "EvRHE_diff",
        "DestFile",
        "Sweep_Type",
        "Type",
        "Cycle",
        "DAC_V",
        "Scanrate",
        "ORR_scan",
        "Jcorr",
        "J_N2_scan",
        "J_O2_diff",
        "J_O2_diff_diff",
        "Analysis_date",
        "J_2nd_diff",
        "Jkin_max",
        "Jkin_min",
        "E_onset",
        "Diff_lim",
        "E_half",
        "I(A)_ring",
        "I(A)_disk",
        "Frac_H2O2",
        "J_ring",
        "n_ORR",
    ]
    KeepColls = [
        "E_AppV_RHE",
        "jmAcm-2",
        "Jcorr",
        "J_N2_scan",
        "Jkin_max",
        "Jkin_min",
        "Frac_H2O2",
        "J_ring",
        "n_ORR",
    ]
    #    SampleCodes = FindExpFolder.LoadSampleCode()
    #    FindExpFolder('VERSASTAT').SampleCodeLst
    #    PostDestDir.mkdir(parents=True,exist_ok=True)
    #    ExpPARovv = EC_loadOVV()
    #    OnlyRecentMissingOVV = runEC.MainPrepareList()
    #    ExpPARovv = ExpPARovv.iloc[100:120]
    OutParsID = pd.DataFrame()
    #    Go1, Go2, Go3 = True, False, False
    #    Go1, Go2, Go3 = False, True, False
    Go1, Go2, Go3 = False, False, True
    #    KL_coeff = KL_coefficients()
    EvRHE_List = [
        0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.45,
        0.5,
        0.55,
        0.6,
        0.65,
        0.7,
        0.75,
        0.8,
        0.9,
        1,
    ]

    def __init__(self):
        self.DestDir = FindExpFolder("VERSASTAT").PostDir

    @staticmethod
    def StartLogging(level_log="INFO"):
        #        level_log = kwargs['level']
        log_fn = FindExpFolder("VERSASTAT").PostDir.joinpath("PostEC_logger.log")
        logging.basicConfig(
            filename=log_fn,
            filemode="w",
            level=level_log,
            format="%(asctime)s %(levelname)s, %(lineno)d: %(message)s",
        )
        logging.warning("Started logging for PostEC script...")

    def applyParallel(dfGrouped, func):
        with Pool(cpu_count() - 1) as p:
            ret_list = p.map(func, [group for name, group in dfGrouped])
        return ret_list

    def check_status(file, verbose=False):
        """Check status will return (status,extra) of filename"""
        PAR_file_test = Path(str(file)).stem
        match = [
            re.search("(?<!VERS|Vers)(AST|postAST|pAST)", str(a))
            for a in PAR_file_test.split("_")
        ]
        if any(match):
            status = "EoL"
            extra = [
                a
                for a in PAR_file_test.split("_")
                if [i for i in match if i][0][0] in a
            ]
            if verbose:
                print(file, status, *extra)
            return status, extra[0]
        #            if any([re.search(r'', i) for i in str(Path(str(file)).stem.split('_'))]):
        else:
            return "BoL", 0

    #            status =
    #            extra =  [0]
    #        return status,extra
    def postEC_Status(files, verbose=False):
        #        files = ['N2_HER_1500rpm_JOS6_pAST-sHA_285_#3_Disc_Parstat']
        if len(files) > 1:
            status_lst, extra_lst = [], []
            for file in files:
                status, extra = PostEC.check_status(file)
                status_lst.append(status)
                extra_lst.append(extra)
            return status_lst, extra_lst
        else:
            return PostEC.check_status(files)

    def OLD_PostOrganizeFolders(TakeRecentList=True):
        postOVV = []
        PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        PAR_version = FileOperations.version

        RunOVV_fn_opts = list(
            FindExpFolder("VERSASTAT").DestDir.rglob(
                "RunOVV_v{0}.xlsx".format(PAR_version)
            )
        )
        RunOVV_fn = [i for i in RunOVV_fn_opts if not "_Conflict" in i.stem][0]
        if RunOVV_fn.is_file() and TakeRecentList == True:
            OvvFromFile = pd.read_excel(RunOVV_fn, index_col=[0])
            status, extra = PostEC.postEC_Status(OvvFromFile.PAR_file.values)
            OvvFromFile = OvvFromFile.assign(
                **{
                    "Date_PAR_EXP": OvvFromFile.PAR_date - OvvFromFile.EXP_date,
                    "Status": status,
                    "Extra": extra,
                }
            )
            OnlyRecentMissingOVV = OvvFromFile
            #            OvvFromFile['Date_PAR_EXP']  = OvvFromFile.PAR_date-OvvFromFile.EXP_date
            #            OvvFromFile['Status'] = OvvFromFile.PAR_file.values
            print("EC OVV loaded from file:{0}".format(RunOVV_fn))
        OnlyRecentMissingOVV = FileOperations.ChangeRoot_DF(
            OnlyRecentMissingOVV, ["Dest_dir", "EXP_dir", "PAR_file"]
        )
        #        CS_parts_PDD = FileOperations.find_CS_parts(PostDestDir)
        #        CS_parts_pOVV = FileOperations.find_CS_parts(OnlyRecentMissingOVV.Dest_dir.iloc[0])
        #        chLst =[]
        #        if CS_parts_PDD[0] != CS_parts_pOVV[0]:
        #            chLst = [CS_parts_PDD[0].joinpath(FileOperations.find_CS_parts(i)[1]) for i in  OnlyRecentMissingOVV.Dest_dir.values]
        #            OnlyRecentMissingOVV['Dest_dir'] = chLst
        #        else:
        #            pass
        postOVVlst, outLst = [], []
        postOVVcols = [
            "DestFilename",
            "SampleID",
            "Status",
            "Status_extra",
            "Electrolyte",
            "Gas",
            "RPM",
            "Scanrate",
            "EXP_date",
            "Type_Exp",
            "SourceFilename",
            "Exp_dir",
        ]
        #        postOVVout = PostEC.FromListgrp(group)
        #        postOVVlst = PostEC.applyParallel(OnlyRecentMissingOVV.groupby('Dest_dir'),PostEC.FromListgrp)
        #        postOVVlst = [outLst.append(PostEC.FromListgrp(i)) for i in OnlyRecentMissingOVV.groupby('Dest_dir')]
        #        for i in OnlyRecentMissingOVV.groupby('Dest_dir'):
        #            PostEC.FromListgrp(i)
        #        try:
        #            postOVVout = pd.DataFrame(postOVVlst,columns=)
        #        except Exception as e:
        #            postOVVout = pd.DataFrame(postOVVlst)
        #        for n,gr in OnlyRecentMissingOVV.groupby(by=['Dest_dir']):
        #            PostEC.FromListgrp(n,gr.EXP_dir.unique()[0])
        #            pass
        #        postOVVlst = [outLst.append(PostEC.FromListgrp(n,gr.EXP_dir.unique()[0])) for n,gr in OnlyRecentMissingOVV.groupby(by=['Dest_dir'])]
        postOVVout = pd.concat(
            [pd.DataFrame(i, columns=postOVVcols) for i in outLst],
            sort=False,
            ignore_index=True,
        )
        postOVVout.to_excel(PostDestDir.joinpath("postEC_Organized.xlsx"))
        return postOVVout


def get_daily_pickle(exp_type=""):
    exp_type = "N2_all"
    today = dt.datetime.now().date()
    _result = {"today": today}
    if exp_type:
        daily_pickle_path = FindExpFolder("VERSASTAT").PostDir.joinpath(
            f"{today:%Y-%m-%d}_{exp_type}{system()}.pkl.compress"
        )
        _daily_exists = daily_pickle_path.exists()
        daily_pickle_path_RAW = FindExpFolder("VERSASTAT").PostDir.joinpath(
            f"{today:%Y-%m-%d}_{exp_type}_{system()}_RAW.pkl.compress"
        )
        _daily_raw_exists = daily_pickle_path_RAW.exists()
        daily_pkl_options_RAW = list(
            FindExpFolder("VERSASTAT").PostDir.rglob(
                f"*_{exp_type}_{system()}_RAW.pkl.compress"
            )
        )
        _result.update(
            {
                "daily_path": daily_pickle_path,
                "_exists": _daily_exists,
                "daily_path_RAW": daily_pickle_path_RAW,
                "_raw_exists": _daily_raw_exists,
                "daily_options": daily_pkl_options_RAW,
            }
        )
    return _result


class CollectPostOVV:
    """Loops over all index files and merges them with the RunOVV"""

    def __init__():
        pass

    @staticmethod
    def LoadPostOVV(reload=False):
        PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        SampleCodes = FindExpFolder().LoadSampleCode()
        #        CS_parts_PDD = FileOperations.find_CS_parts(PostDestDir)
        if reload == True:
            postOVVout = CollectPostOVV.LoadIndexes(reload=True)
        else:
            try:
                postOVVout = CollectPostOVV.LoadIndexes(reload=False)
            except Exception as e:
                logging.warning(
                    "CollectPostOVV no Indexes available: {0}. Using postEC_Organized".format(
                        e
                    )
                )
                postOVVout = pd.read_excel(
                    PostDestDir.joinpath("postEC_Organized.xlsx"), index_col=[0]
                )
        #        pd.read_excel(PostDestDir.joinpath('SampleCodeLst.xlsx'))
        #        CS_parts_pOVV = FileOperations.find_CS_parts(postOVVout.Exp_dir.iloc[0])
        #        if CS_parts_PDD[0] != CS_parts_pOVV[0]:
        #            chLst = [CS_parts_PDD[0].joinpath(FileOperations.find_CS_parts(i)[1]) for i in  postOVVout.SourceFilename.values]
        #            postOVVout['SourceFilename'] = chLst
        #        else:
        #            pass
        postSample = pd.merge(postOVVout, SampleCodes, on="SampleID", how="left")
        print("Types:", " , ".join([str(i) for i in postSample.Type_output.unique()]))
        postSample.PAR_file = postSample.PAR_file.astype(str)
        postSample = FileOperations.ChangeRoot_DF(
            postSample,
            [
                "EXP_dir",
                "Dest_dir",
                "PAR_file",
                "PAR_file_Ring",
                "ORR_act_N2_bg",
                "DestFile",
                "SourceFilename",
            ],
        )
        return postSample

    #    def RunFolderCopy(serie):
    #        postOVVlst = [outLst.append(PostEC.FromListgrp(n,gr.EXP_dir.unique()[0])) for n,gr in serie.groupby(by=['Dest_dir'])]
    #        return postOVVlst
    @staticmethod
    def LoadIndexes(reload=False):
        IndexOVV_fn = FindExpFolder("VERSASTAT").DestDir.joinpath(
            "IndexOVV_v{0}.xlsx".format(FileOperations.version)
        )

        if IndexOVV_fn.exists() and not reload:
            Index_merged = pd.read_excel(IndexOVV_fn, index_col=[0])
            Index_merged = FileOperations.ChangeRoot_DF(
                Index_merged,
                [
                    "EXP_dir",
                    "Dest_dir",
                    "PAR_file",
                    "PAR_file_Ring",
                    "ORR_act_N2_bg",
                    "DestFile",
                    "SourceFilename",
                ],
            )
            logger.info("PostEC loaded IndexOVV from recent: {0}".format(IndexOVV_fn))
        else:
            logger.info(
                "PostEC reloading IndexOVV from Index files and Exp dir files!!"
            )
            OnlyRecentMissingOVV = ECRunOVV(load=1).index
            #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
            OnlyRecentMissingOVV = FileOperations.ChangeRoot_DF(
                OnlyRecentMissingOVV, []
            )
            OnlyRecentMissingOVV.PAR_file = OnlyRecentMissingOVV.PAR_file.astype(str)
            #        if index_source == 'ExpDirs':
            idx_files = [
                list(Path(i).rglob("**/*index*.xlsx"))
                for i in OnlyRecentMissingOVV.Dest_dir.unique()
                if list(Path(i).rglob("**/*index.xlsx"))
            ]
            #            for i in OnlyRecentMissingOVV.Dest_dir.unique():
            #        [idx_files.append([a for a in a if a]) for a in [(Path(i).rglob('index.xlsx')) for i in OnlyRecentMissingOVV.Dest_dir.unique()]]
            #        idx_dir = FindExpFolder('VERSASTAT').IndexDir
            #        idx_files = idx_dir.rglob('*.xlsx')
            #    subset=['PAR_file','DestFile','Type_output','Script_run_date']
            idx_lst = set([a for i in idx_files for a in i])
            idx_mtime = [
                (i, (dt.datetime.now() - dt.datetime.fromtimestamp(i.stat().st_mtime)))
                for i in idx_lst
            ]
            #            print(f'len {len(idx_lst)} and set {len(set(idx_lst))}')
            alst = (
                []
            )  #  Alternative = pd.concat([[pd.read_excel(c,index_col=[0]) for c in a ] for b in idx_files],sort=False,ignore_index=True)
            for idxfp in idx_lst:
                df = pd.read_excel(idxfp, index_col=[0])
                df["IndexSource"] = idxfp
                alst.append(df)
            Index_from_expdirs_all = pd.concat(
                [i for i in alst], sort=False, ignore_index=True
            )
            Index_from_expdirs_all.sort_values(
                "Script_run_date", ascending=False, inplace=True
            )

            Index_from_expdirs = Index_from_expdirs_all.drop_duplicates(keep="first")
            Index_from_expdirs = FileOperations.ChangeRoot_DF(Index_from_expdirs, [])

            idx_exp_tDelta = [
                (n, pd.to_datetime(dt.datetime.now()) - i["Script_run_date"])
                for n, i in Index_from_expdirs.iterrows()
            ]
            Index_from_expdirs = Index_from_expdirs.assign(
                **{
                    "Source": "ExpDirs",
                    "Time_since_run": [pd.to_timedelta(i[1]) for i in idx_exp_tDelta],
                }
            )
            #            Index_from_expdirs['Time_since_run'] = [pd.to_timedelta(pd.to_datetime(datetime.now())-i) for i in Index_from_expdirs['Script_run_date'].values]
            #            limit = pd.to_timedelta('7h')
            #            ['Time_since_run'] = [pd.to_timedelta(pd.to_datetime(datetime.now())-i) for i in Index['Script_run_date'].values]
            #        Index = Index.loc[Index['Time_since_run'] < limit]
            #            Index = Index.iloc[dups].loc[Index['Time_since_run'] < limit]
            #            else:
            #                dups.append(gr.Time_since_run.idxmin())
            #        1elif index_source == 'IndexDir':
            IndexDir_idxfiles = list(
                FindExpFolder("VERSASTAT").IndexDir.rglob("*.xlsx")
            )

            Index_from_idxdir_all = pd.concat(
                [
                    pd.read_excel(i, index_col=[0]).assign(IndexSource=i)
                    for i in IndexDir_idxfiles
                ],
                sort=False,
                ignore_index=True,
            )

            Index_from_idxdir_all.sort_values(
                "Script_run_date", ascending=False, inplace=True
            )
            Index_from_idxdir = Index_from_idxdir_all.drop_duplicates(keep="first")

            Index_from_idxdir = FileOperations.ChangeRoot_DF(Index_from_idxdir, [])
            Index_from_idxdir = Index_from_idxdir.assign(**{"Source": "IndexDir"})
            Index_from_idxdir["Time_since_run"] = [
                pd.to_timedelta(pd.to_datetime(dt.datetime.now()) - i)
                for i in Index_from_idxdir["Script_run_date"].values
            ]
            #        dup_idxdir = Index_from_idxdir.loc[Index_from_idxdir.DestFile.duplicated() == True]
            dups_date, singles, others, unused_dups = [], [], [], []
            for n, gr in Index_from_idxdir.groupby(
                ["PAR_file", "DestFile", "Type_output"]
            ):
                #        Indexes.groupby(['PAR_file','DestFile','Type_output','ScanRate','Segment']):
                if len(gr) > 1:
                    dgr = gr
                    #                print(n,gr.Time_since_run.unique())
                    dups_date.append(gr.Time_since_run.idxmin())
                    unused_dups.append(
                        list(set(gr.index) - {gr.Time_since_run.idxmin()})
                    )
                elif len(gr) == 1:
                    singles.append(gr.index[0])
                else:
                    others.append(gr.index)
            dup_fltr_idxdir = Index_from_idxdir.loc[singles + dups_date]
            #        Indexes = pd.merge(Index_from_expdirs,Index_from_idxdir, on=['PAR_file','DestFile','Type_output','ScanRate','Segment','Sweep_Type','Source'])
            Indexes = pd.concat([Index_from_expdirs, dup_fltr_idxdir], sort=False)
            #            Indexes['Time_since_run'] = [pd.to_timedelta(pd.to_datetime(datetime.now())-i) for i in Indexes['Script_run_date'].values]
            Indexes = Indexes.dropna(
                subset=["PAR_file", "DestFile", "Type_output"]
            ).reset_index()
            dups_date, singles, others = [], [], []

            Idxgr = Indexes.groupby(["PAR_file", "DestFile", "Type_output"])
            for n, gr in Idxgr:
                #        Indexes.groupby(['PAR_file','DestFile','Type_output','ScanRate','Segment']):
                if len(gr) > 1:
                    dgr = gr
                    idxmin = gr.Time_since_run.idxmin()
                    #                print(n,gr.Time_since_run.unique())
                    dups_date.append([idxmin, gr.loc[idxmin, "Source"]])
                elif len(gr) == 1:
                    singles.append(gr.index[0])
                else:
                    others.append(gr.index)
            #        for n2,gr2 in OnlyRecentMissingOVV.groupby('PAR_file'):
            #            if len(gr2) > 1:
            #                dgr2 = gr2
            #        Index = Index.iloc[dups].loc[Index['Time_since_run'] < limit]
            Index = Indexes.loc[singles + [i[0] for i in dups_date]].dropna(
                subset=["DestFile"]
            )
            #        for a in Index.DestFile.values:
            #            try: Path(a).is_file()
            #            except: print(a)
            #            if not any([Path(i).exists() for i in Index.DestFile.values]):
            #                Index = FileOperations.ChangeRoot_DF(Index,['PAR_file','DestFile']) 'EXP_dir','Dest_dir','PAR_file','PAR_file_Ring','ORR_act_N2_bg','DestFile','SourceFilename'
            Index = FileOperations.ChangeRoot_DF(Index, [])
            Index = Index.assign(
                **{
                    "Type_Exp": Index["Type_output"],
                    "SourceFilename": [Path(str(i)) for i in Index["DestFile"].values],
                }
            )
            #        Index['Type_Exp'] = Index['Type_output']
            #        Index['SourceFilename'] = [Path(str(i)) for i in Index['DestFile'].values]
            Index.PAR_file = Index.PAR_file.astype(str)
            Index_undup = Index.loc[
                (
                    Index.duplicated(
                        subset=[
                            "PAR_file",
                            "DestFile",
                            "Type_output",
                            "Time_since_run",
                            "Source",
                        ]
                    )
                    == False
                )
            ]
            idx_merge_cols = [
                i
                for i in Index_undup.columns
                if i in OnlyRecentMissingOVV.columns and not "Segment" in i
            ]
            Index_merged = pd.merge(
                Index_undup, OnlyRecentMissingOVV, on="PAR_file", how="left"
            )
            Index_merged.PAR_file = [
                Path(str(i)) for i in Index_merged["PAR_file"].values
            ]
            new_IndexOVV_target = FileOperations.CompareHashDFexport(
                Index_merged, IndexOVV_fn
            )
            try:
                logger.info(
                    "PostEC re-indexed and saved: {0}".format(new_IndexOVV_target)
                )
            except:
                print("no log")
        return Index_merged

    @staticmethod
    def MatchPostASTs(postOVVout):
        #        postOVVout.postAST.unique()
        #        [(n,len(gr)) for n,gr in postOVVout.groupby('postAST')]
        faillst, fail_index_gr = [], []
        matchAST_lst, non_uniq_lst = [], []
        for nAST, ASTgr in postOVVout.query(
            '(postAST != "no") & (postAST != "postORR")'
        ).groupby(["postAST", "PAR_date", "PAR_file"]):
            nAST, ASTgr
            #            for nDT,grDT in ASTgr.groupby(')
            if ASTgr.PAR_file.nunique() == 1 and ASTgr.Source.nunique() > 1:
                ASTgr_grSource = ASTgr.groupby("Source")
                ASTgr_info = [
                    (n, len(gr), gr.Time_since_run.mean()) for n, gr in ASTgr_grSource
                ]
                if len(set([i[1] for i in ASTgr_info])) == 1:
                    take_source = ASTgr_info[np.argmin([i[2] for i in ASTgr_info])][0]
                    ASTgr = ASTgr_grSource.get_group(take_source)
                    fail_index_source_gr = ASTgr_grSource.get_group(
                        ASTgr_info[np.argmax([i[2] for i in ASTgr_info])][0]
                    )
                    fail_index_gr.append(fail_index_source_gr)

            EC_exp_uniq = [
                (i, ASTgr[i].unique(), ASTgr[i].nunique())
                for i in [
                    c
                    for c in SampleSelection.EC_exp_cols
                    + ["SampleID", "Type_exp", "PAR_file"]
                    if c in ASTgr.columns
                ]
            ]
            EC_exp_non_uniq = [i for i in EC_exp_uniq if i[2] != 1]
            if EC_exp_non_uniq:
                print(
                    "Not unique PAR_date {0},multiple: {1}".format(
                        nAST[1], EC_exp_non_uniq
                    )
                )
                non_uniq_lst.append([nAST, EC_exp_non_uniq, EC_exp_uniq])
                faillst.append(ASTgr)
            EC_exp_query = " & ".join(
                [
                    '({0} == "{1}")'.format(i[0], i[1][0])
                    for i in EC_exp_uniq[1:-1] + [("postAST", ["no"])]
                    if not "Loading" in i[0]
                ]
            )
            past = nAST[1] - pd.to_timedelta(1, unit="D")
            past_slice = postOVVout.query("(PAR_date > @past) & (PAR_date < @nAST[1])")
            past_query = past_slice.query(EC_exp_query)
            if past_query.query(EC_exp_query).empty:
                # expand search to all OVV for similar conditions
                all_query = postOVVout.query(EC_exp_query)
                if not all_query.empty:
                    preAST = tuple(all_query.PAR_file.unique())
                else:
                    preAST = "no-preAST"
            else:
                # find previous preAST measurments
                preAST = tuple(past_query.PAR_file.unique())
            matchAST_lst.append(list(nAST) + [preAST])
        if fail_index_gr:
            fail_index_filter = pd.concat(fail_index_gr)
            postOVVout = postOVVout.loc[
                ~postOVVout.index.isin(fail_index_filter.index), :
            ]
        non_uniq = pd.DataFrame(non_uniq_lst)
        if faillst:
            fails = pd.concat(faillst)
        matchAST = pd.DataFrame(
            matchAST_lst, columns=["postAST", "PAR_date", "PAR_file", "preAST"]
        )
        postOVVout = pd.merge(
            postOVVout, matchAST[["PAR_file", "preAST"]], on="PAR_file", how="left"
        )
        return postOVVout

    #            ASTgr.SampleID.unique()

    @staticmethod
    def MatchECconditions(OVV_df):
        #        postOVVout.postAST.unique()
        #        [(n,len(gr)) for n,gr in postOVVout.groupby('postAST')]
        matchAST_lst = []
        #        'DW16_2018-03-06 00:00:00_no_0.1MHClO4+10mMH2O2_1.0_0.379'
        OVV_df["PAR_date_day"] = [
            dt.datetime.strftime(i, format="%Y-%m-%d")
            for i in OVV_df.PAR_date.fillna(dt.date(1970, 12, 12)).to_list()
        ]
        #        [pd.datetime.strftime(pd.to_datetime(i),format='%Y-%m-%d') for i in postOVVout.PAR_date.fillna(0).to_list()]
        EC_label_cols = [
            "SampleID",
            "pH",
            "Electrolyte",
            "Loading_cm2",
            "postAST",
            "PAR_date_day",
        ]
        post_prev_cols = OVV_df.columns
        #        +[i for i in SampleSelection.EC_exp_cols if i not in ['RPM','Gas']]
        for nAST, ASTgr in OVV_df.groupby(EC_label_cols):
            nAST, ASTgr
            #            for nDT,grDT in ASTgr.groupby(')
            minDT, maxDT = ASTgr.PAR_date.min(), ASTgr.PAR_date.max()
            deltaDT = maxDT - minDT
            #            par_Day = pd.datetime.strftime(nAST[-1],format='%Y-%m-%d')
            EC_exp_query = "_".join([str(i) for i in list(nAST)])
            EC_exp_nodate = "_".join([str(i) for i in list(nAST)[0:-1]])

            matchAST_lst.append(
                pd.DataFrame(
                    [
                        (i, EC_exp_query, EC_exp_nodate, deltaDT)
                        for i in ASTgr.PAR_file.unique()
                    ],
                    columns=["PAR_file", "ECexp", "ECuniq", "EC_deltaDT"],
                )
            )

        EC_exp_match = pd.concat(
            [i for i in matchAST_lst], ignore_index=True, sort=False
        )
        OVV_df = pd.merge(OVV_df, EC_exp_match, on=["PAR_file"], how="left")
        print(
            'Added columns: "{0}" to postOVV with len({1})'.format(
                ", ".join(list(set(post_prev_cols) - set(OVV_df.columns))), len(OVV_df)
            )
        )
        return OVV_df


#            ASTgr.SampleID.unique()
#        merge_cols = [i for i in Index.columns if i in OnlyRecentMissingOVV.columns and not 'Segment' in i]
#        p2,ovv2 = Index.set_index(merge_cols), OnlyRecentMissingOVV.set_index(merge_cols)
#        merge  = p2.update(ovv2)
#        merge  =  p2.combine_first(ovv2)
# else:
#        AllEIS_BoL = pd.concat([pd.read_excel(i) for i in list(PostDestDir.joinpath('EIS','0.1MH2SO4').rglob('*BoL*'))])
#        AllEIS_EoL = pd.concat([pd.read_excel(i) for i in list(PostDestDir.joinpath('EIS','0.1MH2SO4').rglob('*EoL*'))])
#        AllEIS_BoL = AllEIS_BoL.loc[(AllEIS_BoL['Unnamed: 0'] > 0.2901) & (AllEIS_BoL['Unnamed: 0'] < 0.301) & (AllEIS_BoL.SampleID != 'O2'),:]
#        AllEIS300_EoL = AllEIS_EoL.loc[(AllEIS_EoL['Unnamed: 0'] > 0.2901) & (AllEIS_EoL['Unnamed: 0'] < 0.301) & (AllEIS_EoL.SampleID != 'O2'),:]
#        .query('(EXP_date > 20181001)')
#        refl = []
#        for a in postOVVout.SampleID.values:
#            ScodeRef = SampleCodes.loc[SampleCodes.SampleID == a,:]
#            if ScodeRef.empty:
#                Scode = EISovv['SampleID'].unique()[0]
#            else:
#                Scode = ScodeRef.Sample.values[0]
#            refl.append(Scode)
#        postOVVout['SampleLabel'] = refl
#    return postOVVout
#    for a in postOVVout.SampleID.values:
#            ScodeRef = SampleCodes.loc[SampleCodes.SampleID == a,:]
#            if ScodeRef.empty:
#                Scode = EISovv['SampleID'].unique()[0]
#            else:
#                Scode = ScodeRef.Sample.values[0]
#            refl.append(Scode)
#        postOVVout['SampleLabel'] = refl
# postOVVout.loc[postOVVout.Type_Exp == 'EIS_Combined']
# def recently_modified(file,20):
#    file_mtime = pd.to_datetime(DestFile.stat().st_mtime,unit='s')
#%%
class Load_from_Indexes:
    """This class loads the parameters of Electrochemical Data files and merge it with the Overview"""

    SampleCodes = FindExpFolder().LoadSampleCode()
    EC_label_cols = [
        "SampleID",
        "pH",
        "Electrolyte",
        "Loading_cm2",
        "postAST",
        "PAR_date_day",
    ]
    PostDestDir = FindExpFolder("VERSASTAT").PostDir

    def __init__(self, **kwargs):
        if "reload" in kwargs:
            self.postOVVout = CollectPostOVV.LoadPostOVV(kwargs["reload"])
            print(
                "Exp types found in overview: {0}".format(
                    ", ".join([str(i) for i in self.postOVVout.Type_Exp.unique()])
                )
            )
            pass

    @staticmethod
    def PreparePostOVV(fastload=False):
        postOVV_pickle_path = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "PostOVVout_v20_{0}.pkl.compress".format(system())
        )
        if postOVV_pickle_path.is_file():
            tdelta = dt.datetime.now() - dt.datetime.fromtimestamp(
                postOVV_pickle_path.stat().st_mtime
            )
            if tdelta.seconds > 600:
                fastload = False
                print(f"Fastload overwrite to False, {tdelta}")

        if fastload == True:
            try:
                postOVVout = pd.read_pickle(postOVV_pickle_path, compression="xz")
                return postOVVout
            except Exception as e:
                print("Load postOVVout from pickle error: ", e)
                LoadOVV = Load_from_Indexes(reload=True)
        else:
            LoadOVV = Load_from_Indexes(reload=True)

        postOVVout = LoadOVV.postOVVout
        print("Types:", " , ".join([str(i) for i in postOVVout.Type_output.unique()]))
        postOVVout.Loading_cm2 = np.round(postOVVout.Loading_cm2, 3)
        postOVVout = CollectPostOVV.MatchPostASTs(postOVVout)
        postOVVout = CollectPostOVV.MatchECconditions(postOVVout)
        postOVVout.PAR_file = postOVVout.PAR_file.astype(str)

        postOVVout["PAR_date_day"] = [
            pd.datetime.strftime(pd.to_datetime(i), format="%Y-%m-%d")
            for i in postOVVout.PAR_date.fillna(0).values
        ]
        postOVVout = FileOperations.ChangeRoot_DF(postOVVout, [], coltype="string")

        postOVVout.to_pickle(postOVV_pickle_path, compression="xz")
        return postOVVout

    def CollectAllExpTypeOVV():
        PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        today = datetime.today()
        postOVVout = Load_from_Indexes.PreparePostOVV(fastload=False)  # len(22965)
        #        postOVVout.PAR_file = postOVVout.PAR_file.astype(str)
        # === Loading preparation overview of Samples and merging with the data from Characterization techniques === #
        SampleCodes = PostChar.SampleCodeChar()
        #
        Reload_set = True
        logger = start_logger()
        EIS_pars = Load_from_Indexes.EIS_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # EIS_Pars2 6745, 22813
        HPRR_pars = Load_from_Indexes.HPRR_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # HPRR 1668
        Cdl_pars = Load_from_Indexes.N2_pars_OVV(reload=Reload_set)  # Cdl runs 20322
        Cdl_pars_catan = MergeEISandCdl.splitcol_Sweep_Cdl(Cdl_pars)  # 10342
        HER_pars = Load_from_Indexes.HER_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # 2539
        OER_pars = Load_from_Indexes.OER_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # run 1347
        if list(
            PostDestDir.rglob(
                f"{today.year}-{today.month}-*_ORR_pars_{system()}.pkl.compress"
            )
        )[-1].is_file():
            ORR_pars = Load_from_Indexes.ORR_pars_OVV(
                postOVVout, SampleCodes, reload=Reload_set
            )  # ORR 1908

        ORR_pars.to_pickle(
            PostDestDir.joinpath(
                f"{today.year}-{today.month}-{today.day}_ORR_pars_{system()}.pkl.compress"
            )
        )
        EIS_pars.to_pickle(
            PostDestDir.joinpath(
                f"{today.year}-{today.month}-{today.day}_EIS_pars_{system()}.pkl.compress"
            )
        )

    # FindExpFolder().LoadSampleCode()
    #        SampleCodes =  ExportECfromCV.SampleCodes
    #        SampleSelect_all = SampleSelection('*','*')
    #        SampleCodesChar = SampleSelect_all.Prep_EA_BET
    #        SampleCodes = pd.merge(SampleCodes,SampleCodesChar,how='left',on='SampleID',suffixes=('','_char')).drop_duplicates(subset=['SampleID','N_content'])
    # === Start preparing pars OVV from index per Experimental type === #
    #        postOVVout,SampleCodes = pd.DataFrame(),pd.DataFrame()

    def extraPostOVV():
        OnlyRecentMissingOVV = run_PAR_DW.ECRunOVV(load=1).index
        # === Checking expirements from index to analyzed data=== #
        [
            (i)
            for i, gr in OnlyRecentMissingOVV.query('PAR_exp == "EIS"').groupby(
                "SampleID"
            )
            if gr.Loading_cm2.nunique() > 1
        ]
        [
            (i)
            for i, gr in postOVVout.query('PAR_exp == "EIS"').groupby("SampleID")
            if gr.Loading_cm2.nunique() > 1
        ]

        eismiss = OnlyRecentMissingOVV.loc[
            OnlyRecentMissingOVV.PAR_file.isin(
                [
                    i
                    for i in OnlyRecentMissingOVV.query(
                        'PAR_exp == "EIS"'
                    ).PAR_file.values
                    if i not in postOVVout.PAR_file.values
                ]
            )
        ].sort_values(
            by="PAR_date",
        )  # 40

        eismiss.to_excel(
            FindExpFolder("VERSASTAT").PostDir.joinpath("OVV_EIS_missing.xlsx")
        )
        orrmiss = OnlyRecentMissingOVV.loc[
            OnlyRecentMissingOVV.PAR_file.isin(
                [
                    i
                    for i in OnlyRecentMissingOVV.query(
                        'PAR_exp == "ORR" & Electrode != "Pt_ring"'
                    ).PAR_file.values
                    if i not in ORR_pars.PAR_file.values
                ]
            )
        ].sort_values(
            by="PAR_date",
        )  # 279
        #        orrmiss = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.PAR_file.isin([i for i in OnlyRecentMissingOVV.query('PAR_exp == "ORR"').PAR_file.values if i not in ORR_pars.PAR_file.values])].sort_values(by='PAR_date',)
        orrmiss.to_pickle(PostDestDir.joinpath("ORR_missing.pkl.compress"))
        SampleSelection.EC_exp_cols + "SampleID" + EvRHE

        for n, gr in Cdl_pars.groupby(
            [i for i in SampleSelection.EC_exp_cols if i in Cdl_pars.columns]
        ):
            fig, ax = plt.subplots()
            for sID, sgr in gr.groupby("SampleID"):
                sgr.plot(
                    y="Cdl",
                    x="Qad",
                    c="BET_cat_agg_x",
                    colormap="jet",
                    kind="scatter",
                    title="Cdl in acid",
                    ax=ax,
                )

        EIS_pars.query(SampleSelection.acid1500).query('Gas == "O2" & pH == 1 ').plot(
            x="BET_cat_agg", y="Rct", kind="scatter", c="N_content", colormap="viridis"
        )
        mcls = [i for i in EIS_pars.columns if i in Cdl_pars.dropna(axis=1).columns]
        mcls2 = [
            i
            for i in SampleSelection.EC_exp_cols + ["SampleID", "E_RHE"]
            if i in EIS_pars.columns and i in Cdl_pars.dropna(axis=1).columns
        ]
        mcls3 = [
            i
            for i in SampleSelection.EC_exp_cols + ["SampleID", "E_RHE"]
            if i in EIS_pars.columns
            and i in Cdl_pars.dropna(axis=1).columns
            and i in ORR_pars_char.columns
        ]
        [
            (i, EIS_pars[i].dtypes, Cdl_pars[i].dtypes)
            for i in mcls
            if EIS_pars[i].dtypes != Cdl_pars[i].dtypes
        ]
        EIS_Cdl = pd.merge(EIS_pars, Cdl_pars, on=mcls2, how="outer")
        EIS_Cdl_ORR = pd.merge(EIS_Cdl, ORR_pars_char, on=mcls3, how="outer")
        # [['E_RHE','Cdl','Cdlp']]
        ECdl = EIS_Cdl.dropna(how="any", axis=0, subset=["Cdl", "Cdlp"])
        ECdl_ORR = EIS_Cdl.dropna(how="any", axis=0, subset=["Cdl", "Cdlp"])
        test1_alk = ECdl.query(
            '(pH > 7) & (pH < 15) & (E_RHE > 0.494) & (E_RHE < 0.516) & (Sweep_Type_N2 == "cathodic")'
        )
        test1_acid = ECdl.query(
            '(pH < 7) & (E_RHE > 0.494) & (E_RHE < 0.516) & (Sweep_Type_N2 == "cathodic")'
        )
        test1_alk.plot(
            y="Cdl",
            x="Qad",
            c="BET_cat_agg_x",
            colormap="jet",
            kind="scatter",
            title="Cdl in alkaline",
        )
        test1_alk.plot(
            y="Cdl_corr",
            x="Rct",
            c="BET_cat_agg_x",
            colormap="jet",
            kind="scatter",
            title="Cdl in alkaline",
        )
        test1_acid.plot(
            y="Cdl",
            x="Qad",
            c="BET_cat_agg_x",
            colormap="jet",
            kind="scatter",
            title="Cdl in acid",
        )
        test1_acid.plot(
            y="Cdl",
            x="Rct_kin",
            c="BET_cat_agg_x",
            colormap="jet",
            kind="scatter",
            title="Cdl in acid",
        )
        #        HPRR_pars = pd.merge(HPRR_pars,postOVVout,on='PAR_file',how='left')
        #        print('Leftover SampleIDs: {0}'.format(set(HPRR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())))
        #        HPRR_pars = pd.merge(HPRR_pars,SampleCodes,on='SampleID',how='left')
        # @@ Check POST_AST status from OVV and PRM...
        print(
            "Leftover SampleIDs: {0}".format(
                set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())
            )
        )
        ORR_pars = pd.merge(ORR_pars, SampleCodes, on="SampleID", how="left")
        return HPRR_pars_ovv, EIS_pars_ovv

    def get_EC_index():
        OnlyRecentMissingOVV = ECRunOVV(load=1).index
        #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
        OnlyRecentMissingOVV = FileOperations.ChangeRoot_DF(OnlyRecentMissingOVV, [])
        OnlyRecentMissingOVV.PAR_file = OnlyRecentMissingOVV.PAR_file.astype(str)
        OnlyRecentMissingOVV["Loading_cm2"] = OnlyRecentMissingOVV["Loading_cm2"].round(
            3
        )
        SampleCodes = SampleCodesChar().load
        return OnlyRecentMissingOVV, SampleCodes

    @staticmethod
    def IndexPars_CB_paper():
        postOVVout, SampleCodes = pd.DataFrame(), pd.DataFrame()
        PostECddSeries = FindExpFolder("VERSASTAT").DestDir.joinpath(
            "PostEC/{0}".format(SampleSelection.Series_CB_paper["name"])
        )
        PostECddSeries.mkdir(exist_ok=True, parents=True)
        EIS_pars = Load_from_Indexes.EIS_pars_OVV(
            postOVVout, SampleCodes, reload=False
        )  # EIS_Pars2
        HPRR_pars = Load_from_Indexes.HPRR_pars_OVV(
            postOVVout, SampleCodes, reload=False
        )  # HPRR
        ORR_pars = Load_from_Indexes.ORR_pars_OVV(
            postOVVout, SampleCodes, reload=False
        )  # ORR
        Cdl_pars = Load_from_Indexes.N2_pars_OVV(reload=False)
        HER_pars = Load_from_Indexes.HER_pars_OVV(postOVVout, SampleCodes, reload=False)
        OER_pars = Load_from_Indexes.OER_pars_OVV(postOVVout, SampleCodes, reload=False)
        CBsamples = SampleSelection.Series_CB_paper["sIDs"]

        EIS_CB_paper = EIS_pars.loc[EIS_pars.SampleID.isin(CBsamples)]  # 7644
        HPRR_CB_paper = HPRR_pars.loc[HPRR_pars.SampleID.isin(CBsamples)]
        HPRR_CB_paper.to_excel(PostECddSeries.joinpath("HPRR_CB_paper.xlsx"))

        ORR_CB_paper = ORR_pars.loc[ORR_pars.SampleID.isin(CBsamples)]
        ORR_CB_paper.to_excel(PostECddSeries.joinpath("ORR_CB_paper.xlsx"))
        Cdl_CB_paper = Cdl_pars.loc[Cdl_pars.SampleID.isin(CBsamples)]
        Cdl_CB_paper.to_excel(PostECddSeries.joinpath("Cdl_CB_paper.xlsx"))

        HER_CB_paper = HER_pars.loc[HER_pars.SampleID.isin(CBsamples)]
        OER_CB_paper = OER_pars.loc[OER_pars.SampleID.isin(CBsamples)]

        Cdl_CB_cath, Cdl_CB_an = Cdl_CB_paper.query(
            'Sweep_Type_N2 == "cathodic"'
        ), Cdl_CB_paper.query('Sweep_Type_N2 == "anodic"')
        merge_cols_catan = [i for i in Cdl_CB_cath.columns if i in Cdl_CB_an.columns]
        Cdl_CB_catan = pd.merge(
            Cdl_CB_cath,
            Cdl_CB_an,
            on=[i for i in merge_cols_catan if i not in SampleSelection.EC_N2Cdl_cols],
            how="left",
            suffixes=["_cat", "_an"],
        )
        Cdl_CB_catan["Cdl_sum"] = Cdl_CB_catan["Cdl_an"] + Cdl_CB_catan["Cdl_cat"]

        return (
            EIS_CB_paper,
            HPRR_CB_paper,
            ORR_CB_paper,
            Cdl_CB_paper,
            HER_CB_paper,
            OER_CB_paper,
        )

    @staticmethod
    def IndexPars_Porph_SiO2():
        postOVVout, SampleCodes = pd.DataFrame(), pd.DataFrame()
        serie = SampleSelection.Series_Porhp_SiO2["sIDslice"]
        EIS_pars = Load_from_Indexes.EIS_pars_OVV(
            postOVVout, SampleCodes, reload=False
        )  # EIS_Pars2
        Cdl_pars = Load_from_Indexes.N2_pars_OVV(reload=False)
        EIS_Porph_SiO2 = EIS_pars.loc[EIS_pars.SampleID.isin(serie)]
        Cdl_Porph_SiO2 = Cdl_pars.loc[Cdl_pars.SampleID.isin(serie)]
        Cdl_Porph_SiO2_cath, Cdl_Porph_SiO2_an = Cdl_Porph_SiO2.query(
            'Sweep_Type_N2 == "cathodic"'
        ), Cdl_Porph_SiO2.query('Sweep_Type_N2 == "anodic"')
        HPRR_pars_char = Load_from_Indexes.HPRR_pars_OVV(
            postOVVout, SampleCodes, reload=False
        )  # HPRR
        ORR_pars_char = Load_from_Indexes.ORR_pars_OVV(
            postOVVout, SampleCodes, reload=False
        )  # ORR

        HER_pars = Load_from_Indexes.HER_pars_OVV(postOVVout, SampleCodes, reload=False)
        OER_pars = Load_from_Indexes.OER_pars_OVV(postOVVout, SampleCodes, reload=False)

        HPRR_Porph_SiO2 = HPRR_pars_char.loc[HPRR_pars_char.SampleID.isin(serie)]
        ORR_Porph_SiO2 = ORR_pars_char.loc[ORR_pars_char.SampleID.isin(serie)]

        HER_Porph_SiO2 = HER_pars.loc[Cdl_pars.SampleID.isin(serie)]
        OER_Porph_SiO2 = OER_pars.loc[Cdl_pars.SampleID.isin(serie)]

        return ORR_Porph_SiO2

    @staticmethod
    def EIS_pars_OVV(reload=False, extra_plotting=False, xls_out=False, **kwargs):

        #        IndexOVV_EISpars_fn_xls = PostDestDir.joinpath('EIS_pars_IndexOVV_v{0}.xlsx'.format(FileOperations.version))
        #        IndexOVV_EISpars_fn = PostDestDir.joinpath('EIS_pars_IndexOVV_v{0}.pkl.compress'.format(FileOperations.version))
        PostDestDir = Load_from_Indexes.PostDestDir
        #        FindExpFolder('VERSASTAT').PostDir
        eis_daily = get_daily_pickle(exp_type="EIS_pars")

        #        today = dt.datetime.now().date()
        #        eis_daily_pickle_path = PostDestDir.joinpath(f'{today.year}-{today.month}-{today.day}_EIS_pars_{system()}.pkl.compress')
        #        eis_daily_pickle_path_RAW = PostDestDir.joinpath(f'{today.year}-{today.month}-{today.day}_EIS_pars_{system()}_RAW.pkl.compress')

        if eis_daily.get("_exists", False) and not reload:
            EIS_pars = pd.read_pickle(eis_daily.get("daily_path"))
            EIS_pars = FileOperations.ChangeRoot_DF(EIS_pars, [], coltype="string")
            logger.info(
                f'Loaded EIS_pars OVV from  daily {today} pickle: {eis_daily.get("daily_path","")}'
            )
        else:
            # @@ Read EIS pars files and extend with columns from Samples
            # try other way::   idx_files_EIS = [list(Path(i).rglob('**/EIS/*pars_v20.xlsx')) for i in OnlyRecentMissingOVV.Dest_dir.unique() if list(Path(i).rglob('**/EIS/*pars_v20.xlsx'))]
            logger.info(f"START reloading EIS_pars OVV from  daily {today}")
            OnlyRecentMissingOVV = ECRunOVV(load=1).index
            #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
            OnlyRecentMissingOVV = FileOperations.ChangeRoot_DF(
                OnlyRecentMissingOVV, []
            )
            OnlyRecentMissingOVV.PAR_file = OnlyRecentMissingOVV.PAR_file.astype(str)
            OnlyRecentMissingOVV["Loading_cm2"] = OnlyRecentMissingOVV[
                "Loading_cm2"
            ].round(3)
            SampleCodes = SampleCodesChar().load

            def read_df(_par_fls):
                #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
                while True:
                    try:
                        i = next(_par_fls)
                        _pp = pd.read_excel(i, index_col=[0])
                        _pp = FileOperations.ChangeRoot_DF(_pp, [], coltype="string")
                        _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
                        _delta_mtime = dt.datetime.now() - _source_mtime
                        _pp = _pp.assign(
                            **{
                                "SourceFilename": i,
                                "source_mtime": _source_mtime,
                                "delta_mtime": _delta_mtime,
                                "basename": i.stem,
                            }
                        )
                        yield _pp
                    except StopIteration:
                        return "all done"
                        print("gen empty")

            #                    finally:
            #                        yield _pp
            #                    _pf = _pp.PAR_file.unique()[0]
            #                    _pfstem = Path(_pf).stem
            #                    _spectraf = list(Path(Path(i).parent).rglob(f'{_pfstem}_v{FileOperations.version}.xlsx' ))[0]
            #                    _spectradf = pd.read_excel(_spectraf )
            #                    yield _pp
            #            bn = 'O2_EIS-range_1500rpm_JOS1_285_5mV_1500rpm_pars_v20.xlsx'
            EIS_OVV = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.PAR_exp == "EIS"]

            col_names = ["File_SpecFit", "File_SpecRaw", "PAR_file"]
            #            +['PAR_file','Segment',EvRHE, 'RPM_DAC']
            #            [ Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' ) for d in EIS_OVV.Dest_dir.unique()]
            _par_files = [
                list(
                    Path(d.joinpath("EIS")).rglob(
                        f"*_pars_v{FileOperations.EIS_version}.xlsx"
                    )
                )
                for d in EIS_OVV.Dest_dir.unique()
            ]
            _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
            #            tt = (i for i in _par_fls if bn in i.name)

            #            __ttp = list(read_df(tt, col_names))
            if not eis_daily.get("_raw_exists", False):
                _pars_lst = list(read_df(_par_fls))
                EIS_pars_all = pd.concat(_pars_lst, sort=False)
                EIS_pars_all.sort_values("delta_mtime", inplace=True)
                EIS_pars_all = EIS_pars_all.reset_index()
                float_cols = set(
                    [
                        a
                        for i in EIS_pars_all.lmfit_var_names.unique()
                        if type(i) == str and not "(" in i
                        for a in i.split(", ")
                    ]
                )
                obj_flt_cols = [
                    i for i in float_cols if str(EIS_pars_all[i].dtype) == "object"
                ]
                EIS_pars_all = EIS_pars_all.dropna(subset=obj_flt_cols, axis=0)
                #                EIS_pars_all[obj_flt_cols] = EIS_pars_all[obj_flt_cols].fillna(value='0')
                wrong_fls = [
                    EIS_pars_all.loc[
                        EIS_pars_all[i].astype(str).str.contains("Parameter")
                    ]
                    for i in obj_flt_cols
                ]
                if wrong_fls:
                    wrong_objflt_df = pd.concat(wrong_fls)
                    fix_dct = {
                        i: [
                            float(v.split("value=")[-1].split(",")[0])
                            for v in wrong_objflt_df[i].values
                        ]
                        for i in obj_flt_cols
                    }
                    fixed_objflt_df = wrong_objflt_df.assign(**fix_dct)
                    EIS_pars_all = pd.concat(
                        [
                            EIS_pars_all.drop(index=wrong_objflt_df.index, axis=0),
                            fixed_objflt_df,
                        ],
                        axis=0,
                        sort=True,
                    )
                EIS_pars_all[list(float_cols)] = EIS_pars_all[list(float_cols)].astype(
                    float
                )
                EIS_pars_all.to_pickle(eis_daily.get("daily_path_RAW"))
            else:
                EIS_pars_all = pd.read_pickle(eis_daily.get("daily_path_RAW"))
            # === TAKING ONLY NEWEST FITTING PARS ===
            #
            #            for n ,gr in EIS_pars_all.groupby(by=col_names):
            #                n,gr
            E_dc_RHE_cols = [
                (np.round(i, 3), np.round(i, 3) * 1e3)
                for i in EIS_pars_all[EvRHE].values
            ]
            EIS_pars_all = EIS_pars_all.assign(
                **{
                    "E_dc_RHE": [i[0] for i in E_dc_RHE_cols],
                    "E_dc_RHE_mV": [i[1] for i in E_dc_RHE_cols],
                }
            )

            EIS_pars_recent = EIS_pars_all.loc[
                (EIS_pars_all.source_mtime > pd.Timestamp(dt.date(2020, 4, 20)))
                & (EIS_pars_all.PAR_file.str.contains("None") == False)
            ]

            EIS_pars_undup = EIS_pars_recent.dropna(subset=col_names).drop_duplicates(
                keep="first"
            )

            #           EIS_pars = EIS_pars.loc[EIS_pars.lmfit_var_names.str.contains('/(')]
            #                    set([a for i in EIS_pars_all.lmfit_var_names.unique() if not '(' in i for a in i.split(', ')])
            # === POST EDITING OF LOADED PARS ===
            EIS_pars_undup = EIS_pars_undup.assign(
                **{"Loading_cm2": EIS_pars_undup["Loading_cm2"].round(3)}
            )
            EIS_pars_undup = post_helper.make_uniform_EvRHE(EIS_pars_undup)
            EIS_pars_undup = CollectPostOVV.MatchECconditions(EIS_pars_undup)

            _oc_OVV = list(EIS_pars_undup.columns.intersection(EIS_OVV.columns))
            EIS_pars_undup = pd.merge(EIS_pars_undup, EIS_OVV, on=_oc_OVV, how="left")

            _oc_SC = list(EIS_pars_undup.columns.intersection(SampleCodes.columns))
            EIS_pars_undup = pd.merge(
                EIS_pars_undup, SampleCodes, how="left", on=_oc_SC
            )

            EIS_pars_BRUTE = EIS_pars_undup.loc[
                (EIS_pars_undup.BRUTE_FIT == 1) | (EIS_pars_undup.FINAL_FIT == 0)
            ]
            EIS_pars = EIS_pars_undup.loc[(EIS_pars_undup.FINAL_FIT == 1)]

            if extra_plotting == True:
                fast_checking_EEC_models = [
                    (1, "EEC_Randles_RWpCPE", 40),
                    (2, "EEC_2CPE", 50),
                    (3, "EEC_2CPEpW", 120),
                    (4, "EEC_RQ_RQ_RW", 100),
                    (5, "EEC_RQ_RQ_RQ", 100),
                    (6, "Randles_RQRQ", 60),
                ]
                #                ['Model(Singh2015_RQRQR)', 'Model(Singh2015_RQRWR)', 'Model(Singh2015_R3RQ)', 'Model(Bandarenka_2011_RQRQR)' ]
                for idx, _modname, n in [fast_checking_EEC_models[1]] + [
                    fast_checking_EEC_models[4]
                ]:
                    modname = f"Model({_modname})"
                    EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query(
                        "pH < 15"
                    ).plot(
                        y="Rs",
                        x="E_RHE",
                        c="pH",
                        colormap="rainbow_r",
                        kind="scatter",
                        ylim=(0, 100),
                        title=modname,
                    )
                    EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query(
                        "pH < 15"
                    ).plot(
                        y="Qad",
                        x="E_RHE",
                        c="pH",
                        colormap="rainbow_r",
                        kind="scatter",
                        ylim=(0, 0.03),
                        title=modname,
                    )
                    EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH < 7").plot(
                        y="nAd",
                        x="E_RHE",
                        c="BET_cat_agg",
                        colormap="rainbow_r",
                        kind="scatter",
                        ylim=(0, 1),
                        title=modname,
                    )
                    EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH < 7").plot(
                        y="nDL",
                        x="E_RHE",
                        c="BET_cat_agg",
                        colormap="rainbow_r",
                        kind="scatter",
                        ylim=(0, 1),
                        title=modname,
                    )
                    EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH < 7").plot(
                        y="Rct",
                        x="E_RHE",
                        c="BET_cat_agg",
                        colormap="rainbow_r",
                        kind="scatter",
                        ylim=(1, 1e6),
                        logy=True,
                        title=modname,
                    )
                    if (
                        not EIS_pars.loc[EIS_pars["Model_EEC"] == modname]
                        .query("pH > 7")
                        .empty
                    ):
                        EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query(
                            "pH > 7"
                        ).plot(
                            y="Qad+Cdlp",
                            x="E_RHE",
                            c="BET_cat_agg",
                            colormap="rainbow_r",
                            kind="scatter",
                            ylim=(0.1, 1e-4),
                            logy=True,
                            title=modname,
                        )
                plt.close()
                #                EIS_pars.query('pH < 17').groupby('Model_EEC').plot(y='RedChisqr',x='E_RHE',colormap='viridis',kind='scatter',ax=ax)
                _porph = EIS_pars.loc[EIS_pars.PAR_file.str.contains("06.05")]

                fig, ax = plt.subplots()
                for n, Hgr in _porph.query("pH < 7").groupby("postAST"):
                    c_set = "g" if n == "no" else "r"
                    Hgr.plot(
                        x="E_RHE",
                        y="Rct_kin",
                        s=50,
                        c=c_set,
                        kind="scatter",
                        label=n,
                        title="EIS, E vs Qad at",
                        ax=ax,
                        ylim=(1e-6, 1),
                        logy=True,
                    )
                plt.show()
                plt.close()

            EIS_pars.to_pickle(eis_daily_pickle_path)
            logger.info(f'EIS_pars OVV to daily pickle: {eis_daily.get("daily_path")}')
        return EIS_pars

    #             dest_files.append({'index' : n, 'PAR_file' : str(r.PAR_file),'EIS_dest_dir' : EIS_dest_dir,
    #         'EIS_dest_Pars' : EIS_dest_dir.joinpath( Path(r.PAR_file).stem + '_pars.xlsx'),
    #         'EIS_dest_spectra' :EIS_dest_dir.joinpath( Path(r.PAR_file).stem + '_Combined.xlsx')
    #         })
    #            EIS_pars_index_p1 = postOVVout.query('Type_output == "EIS_Pars1"')
    ##            EIS_pars_index_p2 = postOVVout.query('Type_output == "EIS_Pars2"')
    #            EIS_pars_indexes = postOVVout.query('Type_output == "EIS_Pars"')
    #            if 'source' in kwargs.keys():
    #                EIS_pars_indexes = EIS_pars_indexes.loc[EIS_pars_indexes.Source == kwargs.get('source','ExpDirs')]
    ##            pars_index_from_read = EIS_get_index_column_names()
    ##            EIS_pars_index = pd.concat([EIS_pars_index_p1,EIS_pars_index_p2])
    ##            EIS_pars_index = postOVVout.groupby('Type_output').get_group('EIS_Pars1')
    #            EIS_pars_spectra = postOVVout.groupby('Type_output').get_group('EIS_AllData_combined').drop_duplicates(subset=['PAR_file','DestFile','Time_since_run'])
    ##            EPtest = EIS_pars_indexes.loc[no_match] # a slice for testing purpose
    ##            test_load_nm = no_matches.loc[no_matches[2].str.contains('Columns not matching! "Loading_cm2" values:'),0].values
    ##            EPtest = EIS_pars_indexes.loc[EIS_pars_indexes.index.isin(test_load_nm)]
    #            EISlst,no_match,faillst = [],[],[]
    def old_EIS_pars_OVV():
        for n, r in EIS_pars_indexes.iterrows():
            try:
                PF_index, PF_index_stem = Path(r.PAR_file), Path(r.PAR_file).stem
                index_sourcefn = r.SourceFilename

                pars = pd.read_excel(index_sourcefn, index_col=[0])
                pars = FileOperations.ChangeRoot_DF(pars, [], coltype="string")
                if not pars.empty:
                    if "Loading_cm2" in pars.columns:
                        pars = pars.assign(
                            **{"Loading_cm2": np.round(pars["Loading_cm2"], 3)}
                        )

                PF_pars_nunq, PF_pars_unq = (
                    pars.PAR_file.nunique(),
                    pars.PAR_file.unique(),
                )
                #                    .rename(columns={'File' : 'PAR_file'})
                #                    pars_from_index = pd.DataFrame([i.split(', ') for i in pars.index],index=pars.index,columns=pars_index_from_read)
                #                    [i for i in pars.columns if i not in pars_from_index.columns]
                #                    [i for i in pars_from_index.columns if i not in pars.columns]
                #                    pd.merge(pars,pars_from_index,how='left',on=pars.index)
                #                    [i.strip('(') for i in pars[pars_index_from_read[0]].to_list()]
                #                    pars = pd.concat([pars,pars_from_index],axis=1)
                #                    [ i for i in pars.columns if pars[i].dtype not in ['float','int','datetime64[ns]'] and i not in force_skip_cols]
                if PF_pars_nunq > 1:
                    print(
                        "Multiple PAR files in read out pars file {}\n attempting to slice based on PAR_file column".format(
                            index_sourcefn
                        )
                    )
                    pars = pars.loc[
                        pars.PAR_file.str.contains("\\\\".join(PF_index.parts[-2::]))
                    ]

                elif PF_pars_nunq == 1:
                    if PF_index_stem != Path(PF_pars_unq[0]).stem:
                        print(
                            "!PAR_files not matching!\nIndex: {}, Pars: {}".format(
                                PF_index_stem,
                                "".join([Path(i).stem for i in PF_pars_unq]),
                            )
                        )
                        faillst.append([n, PF_index, PF_pars_unq[0]])
                    else:
                        pass
            except Exception as e:
                read_fail_msg = "EIS pars read fail: {}".format(e)
                print(read_fail_msg)
                if "No such file or directory" in read_fail_msg:
                    try:
                        print(f"Try to delete index file {e}")
                        FileOperations.unlink_missing(Path(r.IndexSource))
                    except Exception as e2:
                        print(f"Error to try to delete index file {e2}")

                pars = pd.DataFrame()
                faillst.append([n, r.PAR_file, read_fail_msg])

            spectra_files = EIS_pars_spectra.loc[
                EIS_pars_spectra.PAR_file == r.PAR_file
            ].DestFile.unique()
            if len(spectra_files) == 1:  # pd.read_excel(spectra,index_col=[0])
                spf = spectra_files[0]
            elif len(spectra_files) > 1:  # pd.read_excel(spectra,index_col=[0])
                spf = spectra_files[0]
                print(
                    "EIS prep Took 1st spectra file: {0} of {1}".format(
                        spectra_files[0], spectra_files
                    )
                )
            elif len(spectra_files) == 0:
                print("EIS prep Missing spectra file: {0}".format(spectra_files))
                faillst.append([n, r.PAR_file, "missing spectra"])
                spf = None
            #                pars.assign(**{'SpectraFile' : spf})
            if not pars.empty:
                lenp = len(pars)
                if lenp > 0:
                    overlap_cols = [i for i in r.index if i in pars.columns]
                    mismatch_cols = [
                        i for i in overlap_cols if (r[i] != pars[i].unique()[0])
                    ]
                    PF_index_str = r.PAR_file
                    if any([pd.isna(r[i]) for i in mismatch_cols]):
                        index_isna_cols = [i for i in mismatch_cols if pd.isna(r[i])]
                        for col in index_isna_cols:
                            r[col] = pars[col].unique()[0]
                    mismatch_cols = [
                        i for i in overlap_cols if (r[i] != pars[i].unique()[0])
                    ]
                    if any(
                        c in mismatch_cols
                        for c in ["pH", "Electrolyte", "Loading_cm2", "Loading_name"]
                    ):
                        for col in [
                            i
                            for i in mismatch_cols
                            if i in ["pH", "Electrolyte", "Loading_cm2", "Loading_name"]
                        ]:
                            print(
                                "changed for {0} from pars {1} to {2} from index for {3}".format(
                                    col, pars[col].unique()[0], r[col], PF_index_stem
                                )
                            )
                            pars[col] = r[col]
                    mismatch_cols = [
                        i for i in overlap_cols if (r[i] != pars[i].unique()[0])
                    ]
                    [
                        (r[i], pars[i].unique()[0])
                        for i in overlap_cols
                        if (r[i] != pars[i].unique()[0])
                    ]
                else:
                    overlap_cols, mismatch_cols = [], [1, 2]
            else:
                overlap_cols, mismatch_cols = [], "pars empty"

            if len(mismatch_cols) > 0:
                mismatch_values = [
                    (r[i], pars[i].unique()[0])
                    for i in overlap_cols
                    if (r[i] != pars[i].unique()[0])
                ]
                mismatch_msg = (
                    'Columns not matching! "{1}" values: {0} ,\n Skipped: {2}'.format(
                        *mismatch_values, *mismatch_cols, r.SourceFilename
                    )
                )
                print(mismatch_msg)
                no_match.append([n, r.PAR_file, mismatch_msg])
            else:
                #                print('Columns matching ok!'.format(mismatch_cols))
                not_overlap_cols = list(set(r.index) - set(overlap_cols))
                for i in not_overlap_cols:
                    pars = pars.assign(**{i: [r[i]] * lenp})
                pars = pars.assign(**{"SpectraFile": spf})
                EISlst.append(pars)
        FAILS, no_matches = pd.DataFrame(faillst), pd.DataFrame(
            no_match
        )  # for testing purpose
        EIS_pars = pd.concat(
            [i for i in EISlst if not i.empty], sort=False, ignore_index=True
        )
        EIS_pars = post_helper.make_uniform_EvRHE(EIS_pars)
        EIS_pars_char_mcols = [i for i in EIS_pars.columns if i in SampleCodes.columns]
        nonmatching_dtypes = [
            (i, EIS_pars[i].dtype, SampleCodes[i].dtype)
            for i in EIS_pars_char_mcols
            if EIS_pars[i].dtype != SampleCodes[i].dtype
        ]
        nonmt_cls = [i[0] for i in nonmatching_dtypes]
        #            for a,d1,d2 in nonmatching_dtypes:
        #                try:#    SampleCodes[a] = SampleCodes[a].astype(d1)
        #                except:#    SampleCodes[a].fillna(value=0).str.replace(',','.').astype(d1)
        skip_merging_chars = True
        if skip_merging_chars == True:
            EIS_pars_char = EIS_pars
            print("skipped merging chars with EIS pars? {skip_merging_chars}")
        else:
            EIS_pars_char = pd.merge(
                EIS_pars.drop(columns=nonmt_cls),
                SampleCodes.drop(columns=nonmt_cls),
                on=[i for i in EIS_pars_char_mcols if i not in nonmt_cls],
                how="left",
            )
        EIS_pars_char = EIS_pars_char.loc[EIS_pars.Model_EEC != "Model(Singh2015_RQR)"]
        EIS_pars_char.to_pickle(IndexOVV_EISpars_fn)

    #            new_IndexOVV_EISpars_target = FileOperations.CompareHashDFexport(EIS_pars_char,IndexOVV_EISpars_fn)
    #            try:
    #                logger.info('PostEC EIS re-indexed and saved: {0}'.format(new_IndexOVV_EISpars_target))
    #            except:
    #                pass
    #            EIS_pars_char.query('pH < 17').groupby('Model_EEC').plot(y='RedChisqr',x='E_RHE',colormap='viridis',kind='scatter',yscale='log')

    @staticmethod
    def HPRR_pars_OVV(
        postOVVout, SampleCodes, reload=False, extra_plotting=False, xls_out=False
    ):
        #        exp_type = 'H
        IndexOVV_HPRRpars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_HPRR_v{0}.xlsx".format(FileOperations.version)
        )
        if IndexOVV_HPRRpars_fn.exists() and reload != True:
            HPRR_pars_char = pd.read_excel(IndexOVV_HPRRpars_fn, index_col=[0])
            HPRR_pars_char = FileOperations.ChangeRoot_DF(
                HPRR_pars_char, [], coltype="string"
            )
        else:
            # === Making destination directories === #
            PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
            PPD_HPRR = PostDestDir.joinpath("HPRR")
            PPD_HPRR.mkdir(parents=True, exist_ok=True)
            PPD_HPRR_data = PPD_HPRR.joinpath("DataFiles")
            PPD_HPRR_data.mkdir(parents=True, exist_ok=True)

            #        # === Loading Index files for HPRR and reading the Parameters files into one DataFrame === #
            HPRR_pars_index = postOVVout.groupby("Type_output").get_group("HPRR")
            HP_Pars_files = [
                Path(i)
                for i in HPRR_pars_index["SourceFilename"].unique()
                if "_Pars" in Path(i).stem
            ]
            HPRR_pars_raw = pd.concat(
                [pd.read_excel(i, index_col=[0]) for i in HP_Pars_files], sort=False
            )
            HPRR_pars_raw = FileOperations.ChangeRoot_DF(
                HPRR_pars_raw, [], coltype="string"
            )
            HPRR_merge_cols = [
                i
                for i in HPRR_pars_raw.columns
                if i in HPRR_pars_index.columns and not "Segment" in i
            ]
            HPRR_p2, HPRR_ovv2 = HPRR_pars_raw.set_index(
                HPRR_merge_cols
            ), HPRR_pars_index.set_index(HPRR_merge_cols)
            HPRR_pars_ovv = HPRR_p2.join(HPRR_ovv2, rsuffix="_ovv").reset_index()

            HPRR_pars_merge_cols = [
                i
                for i in HPRR_pars_ovv.columns
                if i in postOVVout.columns and not "Segment" in i and not "Unnamed" in i
            ]
            HPRR_pars = pd.merge(
                HPRR_pars_ovv, postOVVout, on=HPRR_pars_merge_cols, how="left"
            )
            #        HPRR_pars = pd.merge(HPRR_pars_ovv,postOVVout,on='PAR_file',how='left')
            print(
                "Leftover SampleIDs: {0}".format(
                    set(HPRR_pars.SampleID.unique())
                    - set(SampleCodes.SampleID.unique())
                )
            )
            HPRR_char_merge_cols = [
                i
                for i in HPRR_pars_ovv.columns
                if i in SampleCodes.columns
                if not "Unnamed" in i
            ]
            HPRR_pars_char = pd.merge(
                HPRR_pars_ovv, SampleCodes, on=HPRR_char_merge_cols, how="left"
            )
            HPRR_pars_char = HPRR_pars_char.drop(
                columns=[i for i in HPRR_pars_char.columns if "Unnamed" in i]
            )
            new_IndexOVV_HPRRpars_target = FileOperations.CompareHashDFexport(
                HPRR_pars_char, IndexOVV_HPRRpars_fn
            )
            logger.info(
                "PostEC HPRR re-indexed and saved: {0}".format(
                    new_IndexOVV_HPRRpars_target
                )
            )
        if extra_plotting:
            try:
                HPRR_pars_char.query(
                    '(RPM_HPRR > 700) & (Loading_cm2 > 0.1) &  (E_name == "E_j0")'
                ).plot(x="AD/AG", y="fit_slope_HPRR", kind="scatter")
            except Exception as e:
                print("HPRR plot fail:", e)
            try:
                HPRR_pars_char.query(
                    '(RPM_HPRR > 700) & (Loading_cm2 > 0.1) &  (E_name == "E_j0")'
                ).plot(x="N_content", y="fit_slope_HPRR", kind="scatter")
            except Exception as e:
                print("HPRR plot fail:", e)
        return HPRR_pars_char

    @staticmethod
    def HER_pars_OVV(
        postOVVout, SampleCodes, reload=False, extra_plotting=False, xls_out=False
    ):
        #        exp_type = 'H
        IndexOVV_HER_pars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_HER_v{0}.xlsx".format(FileOperations.version)
        )
        if postOVVout.empty or SampleCodes.empty:
            reload = False

        if IndexOVV_HER_pars_fn.exists() and reload != True:
            HER_pars_char = pd.read_excel(IndexOVV_HER_pars_fn, index_col=[0])
            HER_pars_char = FileOperations.ChangeRoot_DF(
                HER_pars_char, [], coltype="string"
            )
        else:
            # === Making destination directories === #
            PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
            PPD_HER_OER = PostDestDir.joinpath("HER_OER")
            PPD_HER_OER.mkdir(parents=True, exist_ok=True)
            PPD_HER_OER_data = PPD_HER_OER.joinpath("DataFiles")
            PPD_HER_OER_data.mkdir(parents=True, exist_ok=True)
            #        # === Loading Index files for HPRR and reading the Parameters files into one DataFrame === #
            HER_pars_index = postOVVout.groupby("Type_output").get_group(
                "HER_Jkin_Tafel"
            )
            #            HP_Pars_files = [i for i in HER_pars_index['SourceFilename'].unique() if '_pars' in i.stem]
            HER_pars_raw = pd.concat(
                [
                    pd.read_excel(i, index_col=[0])
                    for i in HER_pars_index["SourceFilename"].unique()
                ]
            )
            HER_pars_raw = FileOperations.ChangeRoot_DF(
                HER_pars_raw,
                [i for i in HER_pars_raw.columns if re.search("([F-f]ile)", i)],
                coltype="string",
            )
            HER_merge_cols = [
                i
                for i in HER_pars_raw.columns
                if i in HER_pars_index.columns
                and not "Segment" in i
                and not "Sweep_Type" in i
            ]
            HER_p2, HER_ovv2 = HER_pars_raw.set_index(
                HER_merge_cols
            ), HER_pars_index.set_index(HER_merge_cols)
            HER_pars_ovv = HER_p2.join(HER_ovv2, rsuffix="_ovv").reset_index()
            #            HER_pars = pd.merge(HER_pars_ovv,postOVVout,on=HEpars_merge_cols,how='left')
            #            OER_pars = pd.merge(HPRR_pars_ovv,postOVVout,on=HPRR_pars_merge_cols,how='left')
            #    #        HPRR_pars = pd.merge(HPRR_pars_ovv,postOVVout,on='PAR_file',how='left')
            #            print('Leftover SampleIDs: {0}'.format(set(HER_.SampleID.unique()) - set(SampleCodes.SampleID.unique())))
            HER_char_merge_cols = [
                i for i in HER_pars_ovv.columns if i in SampleCodes.columns
            ]
            HER_pars_char = pd.merge(
                HER_pars_ovv, SampleCodes, on=HER_char_merge_cols, how="left"
            )
            new_IndexOVV_HERpars_target = FileOperations.CompareHashDFexport(
                HER_pars_char, IndexOVV_HER_pars_fn
            )
            logger.info(
                "PostEC HPRR re-indexed and saved: {0}".format(
                    new_IndexOVV_HERpars_target
                )
            )
        HER_pars_char.loc[(HER_pars_char["Segment #"] > 1)].query(
            '(E_type == "E_onset")'
        ).plot(x="AD/AG", y="TafelSlope", kind="scatter")
        HER_pars_char.loc[(HER_pars_char["Segment #"] > 1)].query(
            '(E_type == "E_onset")'
        ).plot(x="N_content", y="TafelSlope", s=50, c="g", kind="scatter")
        HER_atE = HER_pars_char.loc[
            (HER_pars_char["Segment #"] > 1)
            & np.isclose(HER_pars_char[EvRHE + "_upper"], -0.3, atol=0.02)
        ].query('(E_type == "E_slice")')
        if extra_plotting:
            fig, ax = plt.subplots()
            for n, Hgr in HER_atE.groupby("postAST"):
                c_set = "g" if n == "no" else "r"
                Hgr.plot(
                    x="N_content",
                    y="j_upper",
                    s=50,
                    c=c_set,
                    kind="scatter",
                    label=n,
                    title="HER at -0.3 Vrhe, j vs N_content",
                    ax=ax,
                )
            HER_atE.plot(
                x="N_content", y="j_upper", kind="bar", title="HER, j vs N_content at"
            )
            HER_atE.plot(
                x="AD/AG",
                y="j_upper",
                s=50,
                c="g",
                kind="scatter",
                title="HER, j vs N_content at",
            )
        return HER_pars_char

    @staticmethod
    def OER_pars_OVV(
        postOVVout, SampleCodes, reload=False, extra_plotting=False, xls_out=False
    ):
        #        exp_type = 'H
        IndexOVV_OERpars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_OER_v{0}.xlsx".format(FileOperations.version)
        )
        if IndexOVV_OERpars_fn.exists() and reload != True:
            OER_pars_char = pd.read_excel(IndexOVV_OERpars_fn, index_col=[0])
            OER_pars_char = FileOperations.ChangeRoot_DF(
                OER_pars_char, [], coltype="string"
            )

        else:
            # === Making destination directories === #
            PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
            PPD_HER_OER = PostDestDir.joinpath("HER_OER")
            PPD_HER_OER.mkdir(parents=True, exist_ok=True)
            PPD_HER_OER_data = PPD_HER_OER.joinpath("DataFiles")
            PPD_HER_OER_data.mkdir(parents=True, exist_ok=True)
            #        # === Loading Index files for HPRR and reading the Parameters files into one DataFrame === #
            OER_pars_index = postOVVout.groupby("Type_output").get_group(
                "OER_Jkin_Tafel"
            )
            OER_pars_raw = pd.concat(
                [
                    pd.read_excel(i, index_col=[0])
                    for i in OER_pars_index["SourceFilename"].unique()
                ]
            )
            OER_pars_raw = FileOperations.ChangeRoot_DF(
                OER_pars_raw,
                [i for i in OER_pars_raw.columns if re.search("([F-f]ile)", i)],
                coltype="string",
            )
            OER_merge_cols = [
                i
                for i in OER_pars_raw.columns
                if i in OER_pars_index.columns
                and not "Segment" in i
                and not "Sweep_Type" in i
            ]
            OER_p2, OER_ovv2 = OER_pars_raw.set_index(
                OER_merge_cols
            ), OER_pars_index.set_index(OER_merge_cols)
            OER_pars_ovv = OER_p2.join(OER_ovv2, rsuffix="_ovv").reset_index()

            #            HER_pars = pd.merge(HER_pars_ovv,postOVVout,on=HEpars_merge_cols,how='left')
            #            OER_pars = pd.merge(HPRR_pars_ovv,postOVVout,on=HPRR_pars_merge_cols,how='left')
            #    #        HPRR_pars = pd.merge(HPRR_pars_ovv,postOVVout,on='PAR_file',how='left')
            #            print('Leftover SampleIDs: {0}'.format(set(HER_.SampleID.unique()) - set(SampleCodes.SampleID.unique())))
            OER_char_merge_cols = [
                i
                for i in OER_pars_ovv.columns
                if i in SampleCodes.columns and not "Unnamed" in i
            ]
            OER_pars_char = pd.merge(
                OER_pars_ovv, SampleCodes, on=OER_char_merge_cols, how="left"
            )
            new_IndexOVV_OERpars_target = FileOperations.CompareHashDFexport(
                OER_pars_char, IndexOVV_OERpars_fn
            )
            logger.info(
                "PostEC OER re-indexed and saved: {0}".format(
                    new_IndexOVV_OERpars_target
                )
            )
        OER_pars_char.loc[(OER_pars_char["Segment #"] > 1)].query(
            '(E_type == "E_onset")'
        ).plot(x="AD/AG", y="TafelSlope", kind="scatter")
        OER_pars_char.loc[(OER_pars_char["Segment #"] > 1)].query(
            '(E_type == "E_onset")'
        ).plot(x="N_content", y="TafelSlope", s=50, c="g", kind="scatter")
        if extra_plotting:
            OER_atE = OER_pars_char.loc[
                (OER_pars_char["Segment #"] > 1)
                & np.isclose(OER_pars_char[EvRHE + "_upper"], 1.7, atol=0.02)
            ].query('(E_type == "E_slice")')
            fig, ax = plt.subplots()
            for n, Ogr in OER_atE.groupby("postAST"):
                c_set = "g" if n == "no" else "r"
                Ogr.plot(
                    x="N_content",
                    y="j_upper",
                    s=50,
                    c=c_set,
                    kind="scatter",
                    label=n,
                    title="OER, j vs N_content at",
                    ax=ax,
                )
        return OER_pars_char

    @staticmethod
    def ORR_pars_OVV(reload=False, extra_plotting=False, xls_out=False):
        #        exp_type = 'H
        IndexOVV_ORRpars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_ORR_v{0}.pkl.compress".format(FileOperations.version)
        )
        PostDestDir = Load_from_Indexes.PostDestDir

        orr_daily = get_daily_pickle(exp_type="ORR_pars")

        if IndexOVV_ORRpars_fn.exists() and reload is not True:
            #            ORR_pars_char = pd.read_excel(IndexOVV_ORRpars_fn,index_col=[0])
            ORR_pars_char = pd.read_pickle(IndexOVV_ORRpars_fn)
            ORR_pars_char = FileOperations.ChangeRoot_DF(
                ORR_pars_char, [], coltype="string"
            )
            ORR_pars_char = ORR_pars_char.drop_duplicates(
                subset=ORR_pars_char.columns[0:19]
            )
        elif reload == "pickle":
            IndexOVV_ORRpars_fn_pkl = list(
                PostDestDir.rglob(
                    f"{today.year}-{today.month}-*_ORR_pars_{system()}.pkl.compress"
                )
            )[-1]
            ORR_pars_char = pd.read_pickle(IndexOVV_ORRpars_fn_pkl)
        else:
            # @@ Check POST_AST status from OVV and PRM
            ORR_pars_index = postOVVout.groupby("Type_output").get_group(
                "ORR_Jkin_calc_Pars"
            )
            ORR_pars_index_RRDE = postOVVout.groupby("Type_output").get_group(
                "ORR_Jkin_calc_RRDE"
            )
            ORR_pars_index_RRDE_Chrono = (
                postOVVout.groupby("Type_output")
                .get_group("ORR_Jkin_calc_RRDE_Chrono")
                .drop_duplicates(subset=["PAR_file", "DestFile", "Time_since_run"])
            )  # cathodic
            ORR_Pars_files = [
                i
                for i in ORR_pars_index["SourceFilename"].unique()
                if re.search("(?i)(_pars|_v20)", Path(i).stem) and Path(i).exists()
            ]
            ORR_pars_raw = pd.concat(
                [pd.read_excel(i, index_col=[0]) for i in ORR_Pars_files], sort=False
            )
            ORR_pars_raw.PAR_file.fillna(value=ORR_pars_raw.File, inplace=True)
            ORR_pars = ORR_pars_raw.drop(columns=["File"], axis=1)
            #            .rename(columns={'File' : 'PAR_file'})
            ORR_pars = FileOperations.ChangeRoot_DF(
                ORR_pars,
                [i for i in ORR_pars.columns if re.search("([F-f]ile)", i)],
                coltype="string",
            )
            ORR_pars.PAR_file = ORR_pars.PAR_file.astype(str)
            ORR_pars_index.PAR_file = ORR_pars_index.PAR_file.astype(str)

            OnlyRecentMissingOVV, SampleCodes = Load_from_Indexes.get_EC_index()

            def read_df(_par_fls, read_types=["ORR_pars"]):
                #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
                while True:
                    try:
                        i = next(_par_fls)
                        _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
                        _delta_mtime = dt.datetime.now() - _source_mtime
                        _i_stem = i.stem
                        _pparts = i.parent.parts

                        if "KL" == _pparts[-1]:
                            if _i_stem.startswith("KL_"):
                                _type = "KL_data"
                            else:
                                _type = "KL_unknown"
                        elif "RingDisk" == _pparts[-1]:
                            _type = "ORR_ringdisk"
                        elif "TAFEL" == _pparts[-1]:
                            _type = "Tafel"
                        else:
                            if _i_stem.startswith("ORR_pars"):
                                _type = "ORR_pars"
                            elif _i_stem.startswith("KL_pars"):
                                _type = "KL_pars"
                            elif _i_stem.startswith("O2_ORR") and _i_stem.endswith(
                                f"_RRDE_v{FileOperations.version}"
                            ):
                                _type = "ORR_RRDE"
                            else:
                                _type = "O2_ORR_unknown"

                        _meta = {
                            "SourceFilename": i,
                            "source_mtime": _source_mtime,
                            "delta_mtime": _delta_mtime,
                            "basename": _i_stem,
                            "_type": _type,
                        }
                        if _type in read_types:
                            _pp = pd.read_excel(i, index_col=[0])
                            _pp = FileOperations.ChangeRoot_DF(
                                _pp, [], coltype="string"
                            )
                            _pp = _pp.assign(**_meta)
                        else:
                            _pp = pd.DataFrame(_meta, index=[0])
                        yield _pp
                    except StopIteration:
                        return "all done"
                        print("gen empty")

            if orr_daily.get("_raw_exists", False):
                ORR_pars_all = pd.read_pickle(orr_daily.get("daily_path_RAW"))
            elif ORR_pars_all.get("daily_options", False):
                ORR_pars_all = pd.read_pickle(orr_daily.get("daily_options")[-1])
            else:  # Construct new N2 pars ovv from reading in files
                ORR_OVV = OnlyRecentMissingOVV.loc[
                    OnlyRecentMissingOVV.PAR_exp == "ORR"
                ]
                _par_files = [
                    list(
                        Path(d.joinpath(f"ORR_v{FileOperations.version}")).rglob(
                            "*xlsx"
                        )
                    )
                    for d in ORR_OVV.Dest_dir.unique()
                ]
                _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
                _par_reads = read_df(_par_fls, read_types=["ORR_pars"])
                ORR_pars_all = pd.concat(_par_reads, sort=False, ignore_index=True)
                for n, gr in ORR_pars_all.groupby("PAR_file"):
                    print(
                        n,
                        f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}',
                        ",".join(gr._type.unique()),
                    )
                ORR_pars_all.to_pickle(orr_daily["daily_path_RAW"])

            rrde_fls, emptylst = [], []
            for fl in ORR_pars_index.PAR_file.values:
                rrde_df_slice = ORR_pars_index_RRDE_Chrono.loc[
                    (ORR_pars_index_RRDE_Chrono.PAR_file == fl)
                ]
                if not rrde_df_slice.empty:
                    rrde_df = rrde_df_slice.loc[(rrde_df_slice.Time_since_run.idxmin())]
                    if "Series" in str(type(rrde_df)):
                        spf = rrde_df.DestFile
                    else:
                        if len(rrde_df) == 1:  # pd.read_excel(spectra,index_col=[0])
                            spf = rrde_df.DestFile.unique()[0]
                        elif len(rrde_df) > 1:  # pd.read_excel(spectra,index_col=[0])
                            spf = rrde_df.DestFile.unique()[0]
                            two = rrde_Df
                        #                   print('ORR prep Took 1st spectra file: {0}'.format(rrde_df))
                        else:
                            print("ORR prep Missing spectra file: {0}".format(rrde_df))
                            miss = rrde_df
                            spf = None
                    rrde_fls.append(spf)
                else:
                    emptylst.append(fl)
                    rrde_fls.append(None)
            if len(ORR_pars_index.PAR_file.values) != len(rrde_fls):
                print("ORR mismatch length for adding data")
            else:
                print("ORR pars length matches RRDE Datafiles... adding column")
                ORR_pars_index = ORR_pars_index.assign(**{"RRDE_DataFile": rrde_fls})

            ORR_merge_cols = [
                i
                for i in ORR_pars.columns
                if i in ORR_pars_index.columns and not "Segment" in i
            ]
            p2, ovv2 = ORR_pars.dropna(subset=ORR_merge_cols).set_index(
                ORR_merge_cols
            ), ORR_pars_index.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols)
            ORR_pars_ovv = p2.join(ovv2, rsuffix="_ovv").reset_index()
            #            ORR_pars_ovv.query('(pH < 7)').plot(y='E_onset',x='Loading_cm2',kind='scatter',logy=False)
            #        ORR_pars_ovv = pd.merge(ORR_pars,ORR_pars_index,on=ORR_merge_cols,suffixes=('','_ovv'),how='left')
            #        ORR_pars = pd.merge(ORR_pars,postOVVout,on=['PAR_file','SampleID','Electrolyte','pH','postAST'],how='left',suffixes=('','_ovv'))
            print(
                "Leftover SampleIDs: {0}".format(
                    set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())
                )
            )
            ORR_pars = pd.merge(ORR_pars_ovv, SampleCodes, on="SampleID", how="left")
            # TODO : taking out duplicates based on time_since_run....
            Load_na = ORR_pars.loc[ORR_pars.Loading_cm2.isna()]
            Load_na_missingvalues = [
                (n, *GetSampleID.ink_loading_from_filename(i.PAR_file))
                for n, i in Load_na.iterrows()
            ]
            Load_na_vals = (
                pd.DataFrame(Load_na_missingvalues)
                .rename(columns={1: "Loading_name", 2: "Loading_cm2"})
                .set_index([0])
            )
            ORR_pars.Loading_cm2.fillna(value=Load_na_vals.Loading_cm2, inplace=True)

            ORR_char_merge_cols = [
                i for i in ORR_pars.columns if i in SampleCodes.columns
            ]
            ORR_pars_char = pd.merge(
                ORR_pars, SampleCodes, on=ORR_char_merge_cols, how="left"
            )
            ORR_pars_char = ORR_pars_char.drop(
                columns=[i for i in ORR_pars_char.columns if "Unnamed" in i]
            )
            if ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna()].empty == False:
                ORR_pars_char.Loading_cm2 = ORR_pars_char.Loading_cm2.fillna(
                    value=0.379
                )  # fillna for Loading_cm2
            #            ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna() == True]
            if xls_out:
                IndexOVV_ORRpars_fn = FileOperations.CompareHashDFexport(
                    ORR_pars_char, IndexOVV_ORRpars_fn
                )
            ORR_pars_char.to_pickle(IndexOVV_ORRpars_fn)
            logger.info(
                "PostEC ORR re-indexed and saved: {0}".format(IndexOVV_ORRpars_fn)
            )
        if extra_plotting:
            ORR_pars_char.query("(pH < 14) & (RPM > 900)").plot(
                y="Jkin_075",
                x="E_onset",
                c="pH",
                kind="scatter",
                logy=True,
                colormap="rainbow_r",
                xlim=(0.5, 1),
            )
        return ORR_pars_char

    @staticmethod
    def ORR_KL_pars_OVV(reload=False, extra_plotting=False, xls_out=False):
        #        exp_type = 'H
        PostDestDir = Load_from_Indexes.PostDestDir
        IndexOVV_ORR_KLpars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_ORR-KL_v{0}.pkl.compress".format(FileOperations.version)
        )

        if IndexOVV_ORRpars_fn.exists() and reload is not True:
            ORR_KL_pars = pd.read_excel(IndexOVV_ORR_KLpars_fn, index_col=[0])
            ORR_pars_char = FileOperations.ChangeRoot_DF(
                ORR_pars_char, [], coltype="string"
            )
            ORR_pars_char = ORR_pars_char.drop_duplicates(
                subset=ORR_pars_char.columns[0:19]
            )
        elif reload == "pickle":
            IndexOVV_ORRpars_fn_pkl = list(
                PostDestDir.rglob(
                    f"{today.year}-{today.month}-*_ORR_pars_{system()}.pkl.compress"
                )
            )[-1]
            ORR_pars_char = pd.read_pickle(IndexOVV_ORRpars_fn_pkl)

        else:
            # @@ Check POST_AST status from OVV and PRM
            #            ORR_index_KL_pars = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_KL_pars')
            #            ORR_index_KL_data = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_KL_data')
            #            ORR_pars_index_RRDE_Chrono = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_RRDE_Chrono').drop_duplicates(subset=['PAR_file','DestFile','Time_since_run']) # cathodic
            #            ORR_Pars_files = [i for i in ORR_pars_index['SourceFilename'].unique() if re.search('(?i)(_pars|_v20)', Path(i).stem) and Path(i).exists()]
            #            ORR_pars_raw = pd.concat([pd.read_excel(i,index_col=[0]) for i in ORR_Pars_files],sort=False)

            if orr_daily_pickle_path_RAW.exists():
                N2_pars_all = pd.read_pickle(orr_daily_pickle_path_RAW)
            elif orr_daily_pickle_path_RAW:
                if orr_daily_pickle_path_RAW[-1].exists():
                    N2_pars_all = pd.read_pickle(orr_daily_pickle_path_RAW[-1])
            else:  # Construct new N2 pars ovv from reading in files
                N2_OVV = OnlyRecentMissingOVV.loc[
                    OnlyRecentMissingOVV.PAR_exp == "N2_act"
                ]
                _par_files = [
                    list(Path(d.joinpath("N2_scans_v30")).rglob(f"*.xlsx"))
                    for d in N2_OVV.Dest_dir.unique()
                ]
                _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
                _par_reads = read_df(_par_fls)
                N2_pars_all = pd.concat(_par_reads, sort=False)

                for n, gr in N2_pars_all.groupby("PAR_file"):
                    print(
                        n,
                        f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}',
                        ",".join(gr.N2_type.unique()),
                    )
                N2_pars_all.to_pickle(n2_daily_pickle_path_RAW)

            ORR_pars_raw.PAR_file.fillna(value=ORR_pars_raw.File, inplace=True)
            ORR_pars = ORR_pars_raw.drop(columns=["File"], axis=1)
            #            .rename(columns={'File' : 'PAR_file'})
            ORR_pars = FileOperations.ChangeRoot_DF(
                ORR_pars,
                [i for i in ORR_pars.columns if re.search("([F-f]ile)", i)],
                coltype="string",
            )
            ORR_pars.PAR_file = ORR_pars.PAR_file.astype(str)
            ORR_pars_index.PAR_file = ORR_pars_index.PAR_file.astype(str)

            rrde_fls, emptylst = [], []
            for fl in ORR_pars_index.PAR_file.values:
                rrde_df_slice = ORR_pars_index_RRDE_Chrono.loc[
                    (ORR_pars_index_RRDE_Chrono.PAR_file == fl)
                ]
                if not rrde_df_slice.empty:
                    rrde_df = rrde_df_slice.loc[(rrde_df_slice.Time_since_run.idxmin())]
                    if "Series" in str(type(rrde_df)):
                        spf = rrde_df.DestFile
                    else:
                        if len(rrde_df) == 1:  # pd.read_excel(spectra,index_col=[0])
                            spf = rrde_df.DestFile.unique()[0]
                        elif len(rrde_df) > 1:  # pd.read_excel(spectra,index_col=[0])
                            spf = rrde_df.DestFile.unique()[0]
                            two = rrde_Df
                        #                   print('ORR prep Took 1st spectra file: {0}'.format(rrde_df))
                        else:
                            print("ORR prep Missing spectra file: {0}".format(rrde_df))
                            miss = rrde_df
                            spf = None
                    rrde_fls.append(spf)
                else:
                    emptylst.append(fl)
                    rrde_fls.append(None)
            if len(ORR_pars_index.PAR_file.values) != len(rrde_fls):
                print("ORR mismatch length for adding data")
            else:
                print("ORR pars length matches RRDE Datafiles... adding column")
                ORR_pars_index = ORR_pars_index.assign(**{"RRDE_DataFile": rrde_fls})

            ORR_merge_cols = [
                i
                for i in ORR_pars.columns
                if i in ORR_pars_index.columns and not "Segment" in i
            ]
            p2, ovv2 = ORR_pars.dropna(subset=ORR_merge_cols).set_index(
                ORR_merge_cols
            ), ORR_pars_index.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols)
            ORR_pars_ovv = p2.join(ovv2, rsuffix="_ovv").reset_index()
            #            ORR_pars_ovv.query('(pH < 7)').plot(y='E_onset',x='Loading_cm2',kind='scatter',logy=False)
            #        ORR_pars_ovv = pd.merge(ORR_pars,ORR_pars_index,on=ORR_merge_cols,suffixes=('','_ovv'),how='left')
            #        ORR_pars = pd.merge(ORR_pars,postOVVout,on=['PAR_file','SampleID','Electrolyte','pH','postAST'],how='left',suffixes=('','_ovv'))
            print(
                "Leftover SampleIDs: {0}".format(
                    set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())
                )
            )
            ORR_pars = pd.merge(ORR_pars_ovv, SampleCodes, on="SampleID", how="left")
            # TODO : taking out duplicates based on time_since_run....
            Load_na = ORR_pars.loc[ORR_pars.Loading_cm2.isna()]
            Load_na_missingvalues = [
                (n, *GetSampleID.ink_loading_from_filename(i.PAR_file))
                for n, i in Load_na.iterrows()
            ]
            Load_na_vals = (
                pd.DataFrame(Load_na_missingvalues)
                .rename(columns={1: "Loading_name", 2: "Loading_cm2"})
                .set_index([0])
            )
            ORR_pars.Loading_cm2.fillna(value=Load_na_vals.Loading_cm2, inplace=True)

            ORR_char_merge_cols = [
                i for i in ORR_pars.columns if i in SampleCodes.columns
            ]
            ORR_pars_char = pd.merge(
                ORR_pars, SampleCodes, on=ORR_char_merge_cols, how="left"
            )
            ORR_pars_char = ORR_pars_char.drop(
                columns=[i for i in ORR_pars_char.columns if "Unnamed" in i]
            )
            if ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna()].empty == False:
                ORR_pars_char.Loading_cm2 = ORR_pars_char.Loading_cm2.fillna(
                    value=0.379
                )  # fillna for Loading_cm2
            #            ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna() == True]
            if xls_out:
                IndexOVV_ORRpars_fn = FileOperations.CompareHashDFexport(
                    ORR_pars_char, IndexOVV_ORRpars_fn
                )
            ORR_pars_char.to_pickle(IndexOVV_ORRpars_fn)
            logger.info(
                "PostEC ORR re-indexed and saved: {0}".format(IndexOVV_ORRpars_fn)
            )
        #        ORR_pars_char.query('(pH < 7) & (RPM > 900)').plot(y='Jkin_075',x='AD/AG',kind='scatter',logy=False)
        return ORR_pars_char

    @staticmethod
    def N2_pars_OVV(reload=False, extra_plotting=False, xls_out=False):
        #        exp_type = 'H
        PostDestDir = Load_from_Indexes.PostDestDir
        #        IndexOVV_N2_pars_fn_xls = FindExpFolder('VERSASTAT').PostDir.joinpath('Pars_IndexOVV_CdlN2_v{0}.xlsx'.format(FileOperations.version))
        #        IndexOVV_N2_pars_fn = FindExpFolder('VERSASTAT').PostDir.joinpath('N2Cdl_pars_IndexOVV_v{0}.pkl.compress'.format(FileOperations.version))
        n2_daily = get_daily_pickle(exp_type="N2_all")

        if _n2_daily.get("_exists", False) and reload != True:
            #            Cdl_pars_char = pd.read_excel(IndexOVV_N2_pars_fn,index_col=[0])
            Cdl_pars_char = pd.read_pickle(n2_daily.get("daily_path"))
            Cdl_pars_char = FileOperations.ChangeRoot_DF(
                Cdl_pars_char, [], coltype="string"
            )
        else:
            # @@ Check POST_AST status from OVV and PRM
            logger.info(f"START reloading N2_pars OVV from daily {today:%Y-%m-%d}")
            #            OnlyRecentMissingOVV = ECRunOVV(load=1).index
            #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
            #            OnlyRecentMissingOVV = FileOperations.ChangeRoot_DF(OnlyRecentMissingOVV,[])
            #            OnlyRecentMissingOVV.PAR_file = OnlyRecentMissingOVV.PAR_file.astype(str)
            #            OnlyRecentMissingOVV['Loading_cm2'] = OnlyRecentMissingOVV['Loading_cm2'].round(3)
            #            SampleCodes = SampleCodesChar().load
            OnlyRecentMissingOVV, SampleCodes = Load_from_Indexes.get_EC_index()

            def read_df(_par_fls):
                #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
                while True:
                    try:
                        i = next(_par_fls)
                        _pp = pd.read_excel(i, index_col=[0])
                        _pp = FileOperations.ChangeRoot_DF(_pp, [], coltype="string")
                        _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
                        _delta_mtime = dt.datetime.now() - _source_mtime
                        _i_stem = i.stem
                        if _i_stem.endswith("_BG"):
                            _N2_type = "BG"
                        else:
                            if _i_stem.startswith("CV_"):
                                _N2_type = "CV"
                                if _i_stem.endswith(
                                    f"_first_v{FileOperations.version}"
                                ):
                                    _N2_type = "CV_first"
                            #                                if not 'Scan Rate' in _pp.columns:
                            #                                    'N2_CV_raw =  N2_CV_raw.assign(**{'ScanRate' : [i.split(f'_v{FileOperations.version}')[0].split('_')[-1] for i in N2_CV_raw.basename.to_numpy()]})

                            elif _i_stem.startswith("Cdl_data_"):
                                _N2_type = "Cdl_data"
                            elif _i_stem.startswith("Cdl_pars"):
                                _N2_type = "Cdl_pars"
                            else:
                                _N2_type = "N2_unknown"

                        _pp = _pp.assign(
                            **{
                                "SourceFilename": i,
                                "source_mtime": _source_mtime,
                                "delta_mtime": _delta_mtime,
                                "basename": _i_stem,
                                "N2_type": _N2_type,
                            }
                        )
                        yield _pp
                    except StopIteration:
                        return "all done"
                        print("gen empty")

            if n2_daily.get("_raw_exists", False):
                N2_pars_all = pd.read_pickle(n2_daily.get("daily_path_RAW"))
            elif n2_daily.get("daily_options", False):
                if n2_daily.get("daily_options")[-1]:
                    N2_pars_all = pd.read_pickle(n2_daily.get("daily_options")[-1])
            else:  # Construct new N2 pars ovv from reading in files
                N2_OVV = OnlyRecentMissingOVV.loc[
                    OnlyRecentMissingOVV.PAR_exp == "N2_act"
                ]
                _par_files = [
                    list(Path(d.joinpath("N2_scans_v30")).rglob(f"*.xlsx"))
                    for d in N2_OVV.Dest_dir.unique()
                ]
                _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
                _par_reads = read_df(_par_fls)
                N2_pars_all = pd.concat(_par_reads, sort=False)

                for n, gr in N2_pars_all.groupby("PAR_file"):
                    print(
                        n,
                        f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}',
                        ",".join(gr.N2_type.unique()),
                    )
                N2_pars_all.to_pickle(n2_daily_pickle_path_RAW)

            #
            N2_type_grps = N2_pars_all.groupby("N2_type")

            # N2 CVs TODO add Scan Rate column
            N2_CV_raw = N2_type_grps.get_group("CV").dropna(axis=1, how="all")
            #            N2_CV_raw.plot(x=EvRHE,y='jmAcm-2')

            N2_CV_pivot_SR_lst = []
            for PF, PFgr in N2_CV_raw.groupby("PAR_file"):
                #                PF ,PFgr
                for swp, swgrp in PFgr.groupby("Sweep_Type"):
                    #                    swp, swgrp
                    #                    swgrp.plot(x=EvRHE,y='jmAcm-2')
                    #                    E_T_idx = pd.MultiIndex.from_tuples(zip(swgrp['Elapsed Time(s)'].to_numpy(),swgrp[EvRHE].to_numpy()),names=['Elapsed_Time_s',EvRHE])
                    #                    swgrp.index = E_T_idx
                    #                    {n : len(gr) for n,gr in swgrp.groupby('Segment #')}
                    pvt = swgrp.pivot(
                        index="Elapsed Time(s)",
                        columns="ScanRate_mVs",
                        values=[EvRHE, "jmAcm-2", "Segment #"],
                    )
                    #                    pvt = swgrp.pivot(index=EvRHE,columns='ScanRate_mVs',values='jmAcm-2')
                    pvt.columns = pd.MultiIndex.from_tuples(
                        [(f"{i[0]}_{int(i[1])}", i[1]) for i in pvt.columns]
                    )
                    #                    pvt.rename(columns=pd.MultiIndex.from_tuples([(f'{i[0]}_{int(i[1])}', i[1]) for i in pvt.columns],names=['data','ScanRate_mVs']),inplace=True)
                    indx = pd.MultiIndex.from_tuples(
                        zip(repeat(PF), repeat(swp), pvt.index),
                        names=["PAR_file", "Sweep_Type", EvRHE],
                    )
                    pvt.index = indx
                    N2_CV_pivot_SR_lst.append(pvt)
            #                    for sr, srgrp in PFgr.groupby('ScanRate_mVs'):
            #                        SR = int(sr)
            N2_CV_pivot_SR = pd.concat(N2_CV_pivot_SR_lst, sort=False)

            #            N2Cdl_pars_index = N2_grps.groupby('N2_type').get_group('Cdl_pars')
            #            N2Cdl_pars_files = [Path(i) for i in N2Cdl_pars_index['SourceFilename'].unique() if re.search('(?i)(_pars|_v20)',Path(i).stem) and Path(i).exists()]
            #            cdl = pd.read_excel(N2Cdl_pars_files[0],index_col=[0])
            #                                                                     [i for i in N2Cdl_pars_raw.columns if re.search('([F-f]ile|ORR_act)',i)]
            #            N2Cdl_pars.rename(columns={'Filename' : 'PAR_file'})
            #            EPtest = N2Cdl_pars_index.loc[no_match] # a slice for testing purpose
            #            pd.merge(N2Cdl_pars_raw,N2_CV_index[['PAR_file','DestFile']],on='PAR_file',how='left')
            #            N2Cdl_pars_raw =  N2_type_grps.get_group('Cdl_pars').dropna(axis=1,how='all')
            #            N2Cdl_data_index = postOVVout.groupby('Type_output').get_group('N2_Cdl_data')
            #            N2_CV_index = postOVVout.groupby('Type_output').get_group('N2_CV')
            #            lst, no_match, non_exist = [],[],[]
            #            for n,r in N2Cdl_pars_raw.iterrows():
            #                Cdl_data_file = N2Cdl_data_index.loc[N2Cdl_data_index.PAR_file == r.PAR_file].DestFile.unique()
            #                CV_files = N2_CV_index.loc[N2_CV_index.PAR_file == r.PAR_file].DestFile.unique()
            #                lst.append([set(Cdl_data_file),set(CV_files)])
            #            if len(N2Cdl_pars_raw) == len(lst):
            #                N2Cdl_pars_raw = N2Cdl_pars_raw.assign(**{'Cdl_data_file' : [i[0] for i in lst], 'Cdl_CV_data_files' : [i[1] for i in lst]})
            #            Cdl_pars = pd.concat([i for i in lst],sort=False,ignore_index=True)
            N2Cdl_pars_raw = N2_type_grps.get_group("Cdl_pars").dropna(
                axis=1, how="all"
            )
            N2Cdl_pars_raw.drop_duplicates(
                subset=N2Cdl_pars_raw.columns[0:19], keep="first", inplace=True
            )
            N2Cdl_pars_raw = FileOperations.ChangeRoot_DF(
                N2Cdl_pars_raw, [], coltype="string"
            )
            Cdl_pars = post_helper.make_uniform_EvRHE(N2Cdl_pars_raw)
            Cdl_pars.drop_duplicates(subset=Cdl_pars.columns[0:19], inplace=True)
            Cdl_pars_merge_cols = [
                i
                for i in Cdl_pars.columns
                if i in SampleCodes.columns and not "Unnamed" in i
            ]
            Cdl_pars_char = pd.merge(
                Cdl_pars, SampleCodes, on=Cdl_pars_merge_cols, how="left"
            )
            Cdl_pars_char.drop_duplicates(
                subset=Cdl_pars_char.columns[0:19], inplace=True
            )
            if xls_out:
                new_N2_pars_char_target = FileOperations.CompareHashDFexport(
                    Cdl_pars_char, IndexOVV_N2_pars_fn
                )
                logger.info(
                    "PostEC Cdl N2 CVs re-indexed and saved: {0}".format(
                        new_N2_pars_char_target
                    )
                )
            Cdl_pars_char.to_pickle(IndexOVV_N2_pars_fn)

        Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)').plot(
            y="Cdl",
            x="E_RHE",
            c="BET_cat_agg",
            colormap="viridis",
            kind="scatter",
            ylim=(0, 0.08),
            title="checking plot: Cdl in acid",
        )
        #        Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)').groupby('BET_cat_agg').plot(y='Cdl',x='E_RHE',colormap='viridis',kind='scatter',ylim=(0,0.08),title='Cdl in acid')
        if extra_plotting:
            Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH > 7)').plot(
                y="Cdl",
                x="E_RHE",
                c="BET_cat_agg",
                colormap="viridis",
                kind="scatter",
                ylim=(0, 0.03),
                title="Cdl in alkaline",
            )
            alkCdl = Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH > 7)')
            acidCdl = Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)')
            # 3d plotting
            #        fig = plt.figure()
            #        ax = fig.add_subplot(111, projection='3d')
            #        ax.plot_trisurf(alkCdl.E_RHE,alkCdl.Cdl,alkCdl.BET_cat_agg,cmap=cm.viridis)
            Cdl_atE = Cdl_pars_char.loc[
                (Cdl_pars_char.Sweep_Type_N2 == "cathodic")
                & (np.isclose(Cdl_pars_char["E_RHE"], 0.5, atol=0.02))
            ]
            fig, ax = plt.subplots()
            for n, Ogr in Cdl_atE.query(
                '(Sweep_Type_N2 == "cathodic") & (pH < 7)'
            ).groupby("postAST"):
                c_set = "g" if n == "no" else "r"
                Ogr.plot(
                    x="BET_cat_agg",
                    y="Cdl",
                    s=50,
                    c=c_set,
                    kind="scatter",
                    label=n,
                    title="N2 Cdl vs BET in acid",
                    ax=ax,
                    ylim=(0, 50e-3),
                )
            fig, ax = plt.subplots()
            for n, Ogr in Cdl_atE.query(
                '(Sweep_Type_N2 == "cathodic") & (pH > 7)'
            ).groupby("postAST"):
                c_set = "g" if n == "no" else "r"
                Ogr.plot(
                    x="BET_cat_agg",
                    y="Cdl",
                    s=50,
                    c=c_set,
                    kind="scatter",
                    label=n,
                    title="N2 Cdl vs BET in alk",
                    ax=ax,
                    ylim=(0, 50e-3),
                )
            # ==== #
        return Cdl_pars_char


#        .plot(y='Cdl',x='E_RHE',c='BET_cat_agg',colormap='viridis',kind='scatter',ylim=(0,0.03),title='Cdl in alkaline')


def get_EIS_pars(reload=False):
    EIS_pars = Load_from_Indexes.EIS_pars_OVV(
        reload=reload, extra_plotting=True
    )  # EIS_Pars2 6745, 17994
    return EIS_pars


#
#    reloadOVV = 1
#    reloading, Reload_set = True, True


def ORR_stat():
    for i in SampleSelection.EC_ORR_kin_par_cols:
        sc = ORR_pars[i]
        print(
            "{} in {:.3f} min {:.3f} mean {:.3f} std {:.3f}".format(
                i, sc.max(), sc.min(), sc.mean(), sc.std()
            )
        )


#%%
if __name__ == "__main__":
    reloadOVV = False
    if reloadOVV:
        PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        today = dt.today()
        postOVVout = Load_from_Indexes.PreparePostOVV(fastload=False)  # len(22965)
        SampleCodes = FindExpFolder().LoadSampleCode()

    #        postOVVout.PAR_file = postOVVout.PAR_file.astype(str)
    # === Loading preparation overview of Samples and merging with the data from Characterization techniques === #
    reloading, Reload_set = True, False
    if 0:
        #        logger = start_logger()
        EIS_pars = Load_from_Indexes.EIS_pars_OVV(
            reload=Reload_set, source="ExpDirs"
        )  # EIS_Pars2 6745, 22813
        HPRR_pars = Load_from_Indexes.HPRR_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # HPRR 1668
        Cdl_pars = Load_from_Indexes.N2_pars_OVV(
            reload=Reload_set, extra_plotting=False, xls_out=False
        )  # Cdl runs 20322
        Cdl_pars_catan = merger.MergeEISandCdl.splitcol_Sweep_Cdl(Cdl_pars)  # 10342
        HER_pars = Load_from_Indexes.HER_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # 2539
        #    OER_pars = Load_from_Indexes.OER_pars_OVV(postOVVout,SampleCodes,reload= Reload_set) # run 1347
        #    if list(PostDestDir.rglob(f'{today.year}-{today.month}-*_ORR_pars_{system()}.pkl.compress'))[-1].is_file():
        ORR_pars = Load_from_Indexes.ORR_pars_OVV(
            postOVVout, SampleCodes, reload=False
        )  # ORR 1908
        ORR_pars.RPM_DAC = ORR_pars.RPM
        HPRR_pars.to_pickle(
            PostDestDir.joinpath(
                f"{today.year}-{today.month}-{today.day}_HPRR_pars_{system()}.pkl.compress"
            )
        )
        Cdl_pars.to_pickle(
            PostDestDir.joinpath(
                f"{today.year}-{today.month}-{today.day}_Cdl_pars_{system()}.pkl.compress"
            )
        )
        ORR_pars.to_pickle(
            PostDestDir.joinpath(
                f"{today.year}-{today.month}-{today.day}_ORR_pars_{system()}.pkl.compress"
            )
        )
        EIS_pars.to_pickle(
            PostDestDir.joinpath(
                f"{today.year}-{today.month}-{today.day}_EIS_pars_{system()}.pkl.compress"
            )
        )
