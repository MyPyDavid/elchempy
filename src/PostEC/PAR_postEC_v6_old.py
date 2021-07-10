# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 21:01:16 2019

@author: DWXMG
@version : 6
"""

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
from itertools import chain
import logging

from pathlib import Path

# import h5py
from multiprocessing import Pool, cpu_count

# import timeit
# import time
from datetime import datetime

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
print("DONT EDIT THIS ONE!!")
##from FolderOrganizing import FileHelper
##from FolderOrganizing import RunEC_classifier
# from FolderOrganizing import PAR_EIS_fit_V2
# import run_PAR_DW
# PostOrganizeFolders(OnlyRecentMissingOVV)
if __name__ == "__main__":
    print(f"Package: {__package__}, File: {__file__}")
    FH_path = Path(__file__).parent.parent.parent.joinpath("FileHelper")
    sys.path.append(str(FH_path))
    sys.path.append(str(Path(__file__).parent.parent.joinpath("indexer")))
    sys.path.append(str(Path(__file__).parent.parent.parent))
    sys.path.append("..")
    #    print(sys.path)
    import FileHelper
    import prepare_input
    import runEC

    #    import EC
    #    sys.path.append(list(FH_path.rglob('*.py')))
    #    import FH_path.joinpath('FindExpFolder.py')
    #    import FindExpFolder.py
    #    from FileHelper import FindExpFolder
    #    from FileHelper.FindExpFolder import *
    #    from .experiments import EIS
    #    from .runEC import run_PAR_DW
    OriginColor = FileHelper.FindExpFolder.LoadOriginColor()
    logger = start_logger()
#    print('\n\n***** run_PAR_DW *****')
#    for k in __file__.__dict__.keys():
#        print(k)
else:
    logger = start_logger()


def start_logger(logger_option_set=True):
    if logger_option_set == True:
        # Gets or creates a logger
        logger = logging.getLogger(__name__)
        # set log level
        logger.setLevel(logging.WARNING)
        # define file handler and set formatter
        file_handler = logging.FileHandler(
            FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath("PostEC_logger.log")
        )
        formatter = logging.Formatter(
            "%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s"
        )
        file_handler.setFormatter(formatter)
        # add file handler to logger
        logger.addHandler(file_handler)
        logger.warning("=========== Started logging PAR PostEC...  =============")
    else:
        print("=====x====== No logging PAR PostEC...  ======x=======")
    return logger


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
    #    SampleCodes = FileHelper.FindExpFolder.LoadSampleCode()
    #    FileHelper.FindExpFolder('VERSASTAT').SampleCodeLst

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
        self.DestDir = FileHelper.FindExpFolder("VERSASTAT").PostDir

    @staticmethod
    def StartLogging(level_log="INFO"):
        #        level_log = kwargs['level']
        log_fn = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "PostEC_logger.log"
        )
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
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        PAR_version = FileHelper.FileOperations.version

        RunOVV_fn_opts = list(
            FileHelper.FindExpFolder("VERSASTAT").DestDir.rglob(
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
        OnlyRecentMissingOVV = FileHelper.FileOperations.ChangeRoot_DF(
            OnlyRecentMissingOVV, ["Dest_dir", "EXP_dir", "PAR_file"]
        )
        #        CS_parts_PDD = FileHelper.FileOperations.find_CS_parts(PostDestDir)
        #        CS_parts_pOVV = FileHelper.FileOperations.find_CS_parts(OnlyRecentMissingOVV.Dest_dir.iloc[0])
        #        chLst =[]
        #        if CS_parts_PDD[0] != CS_parts_pOVV[0]:
        #            chLst = [CS_parts_PDD[0].joinpath(FileHelper.FileOperations.find_CS_parts(i)[1]) for i in  OnlyRecentMissingOVV.Dest_dir.values]
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


class CollectPostOVV:
    """Loops over all index files and merges them with the RunOVV"""

    def __init__():
        pass

    @staticmethod
    def LoadPostOVV(reload=False):
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        SampleCodes = FileHelper.FindExpFolder().LoadSampleCode()
        #        CS_parts_PDD = FileHelper.FileOperations.find_CS_parts(PostDestDir)
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
        #        CS_parts_pOVV = FileHelper.FileOperations.find_CS_parts(postOVVout.Exp_dir.iloc[0])
        #        chLst =[]
        #        if CS_parts_PDD[0] != CS_parts_pOVV[0]:
        #            chLst = [CS_parts_PDD[0].joinpath(FileHelper.FileOperations.find_CS_parts(i)[1]) for i in  postOVVout.SourceFilename.values]
        #            postOVVout['SourceFilename'] = chLst
        #        else:
        #            pass
        postSample = pd.merge(postOVVout, SampleCodes, on="SampleID", how="left")
        print("Types:", " , ".join([str(i) for i in postSample.Type_output.unique()]))
        postSample.PAR_file = postSample.PAR_file.astype(str)
        postSample = FileHelper.FileOperations.ChangeRoot_DF(
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
        IndexOVV_fn = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
            "IndexOVV_v{0}.xlsx".format(FileHelper.FileOperations.version)
        )

        if IndexOVV_fn.exists() and not reload:
            Index_merged = pd.read_excel(IndexOVV_fn, index_col=[0])
            Index_merged = FileHelper.FileOperations.ChangeRoot_DF(
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
            OnlyRecentMissingOVV = runEC.ECRunOVV(load=1).index
            #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
            OnlyRecentMissingOVV = FileHelper.FileOperations.ChangeRoot_DF(
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
            #        idx_dir = FileHelper.FindExpFolder('VERSASTAT').IndexDir
            #        idx_files = idx_dir.rglob('*.xlsx')
            #    subset=['PAR_file','DestFile','Type_output','Script_run_date']
            alst = (
                []
            )  #  Alternative = pd.concat([[pd.read_excel(c,index_col=[0]) for c in a ] for b in idx_files],sort=False,ignore_index=True)
            for a in idx_files:
                for i in a:
                    df = pd.read_excel(i, index_col=[0])
                    alst.append(df)
            Index_from_expdirs_all = pd.concat(
                [i for i in alst], sort=False, ignore_index=True
            )
            Index_from_expdirs_all.sort_values(
                "Script_run_date", ascending=False, inplace=True
            )

            Index_from_expdirs = Index_from_expdirs_all.drop_duplicates(keep="first")
            Index_from_expdirs = FileHelper.FileOperations.ChangeRoot_DF(
                Index_from_expdirs, []
            )

            idx_exp_tDelta = [
                (n, pd.to_datetime(datetime.now()) - i["Script_run_date"])
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
                FileHelper.FindExpFolder("VERSASTAT").IndexDir.rglob("*.xlsx")
            )
            Index_from_idxdir_all = pd.concat(
                [pd.read_excel(i, index_col=[0]) for i in IndexDir_idxfiles],
                sort=False,
                ignore_index=True,
            )

            Index_from_idxdir_all.sort_values(
                "Script_run_date", ascending=False, inplace=True
            )
            Index_from_idxdir = Index_from_idxdir_all.drop_duplicates(keep="first")

            Index_from_idxdir = FileHelper.FileOperations.ChangeRoot_DF(
                Index_from_idxdir, []
            )
            Index_from_idxdir = Index_from_idxdir.assign(**{"Source": "IndexDir"})
            Index_from_idxdir["Time_since_run"] = [
                pd.to_timedelta(pd.to_datetime(datetime.now()) - i)
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
            #                Index = FileHelper.FileOperations.ChangeRoot_DF(Index,['PAR_file','DestFile']) 'EXP_dir','Dest_dir','PAR_file','PAR_file_Ring','ORR_act_N2_bg','DestFile','SourceFilename'
            Index = FileHelper.FileOperations.ChangeRoot_DF(Index, [])
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
            new_IndexOVV_target = FileHelper.FileOperations.CompareHashDFexport(
                Index_merged, IndexOVV_fn
            )
            try:
                logger.info(
                    "PostEC re-indexed and saved: {0}".format(new_IndexOVV_target)
                )
            except:
                print("no log")
        return Index_merged


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
class Load_from_Indexes:
    """This class loads the parameters of Electrochemical Data files and merge it with the Overview"""

    SampleCodes = FileHelper.FindExpFolder().LoadSampleCode()
    EC_label_cols = [
        "SampleID",
        "pH",
        "Electrolyte",
        "Loading_cm2",
        "postAST",
        "PAR_date_day",
    ]

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

        postOVV_pickle_path = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "PostOVVout_v20_{0}.pkl.compress".format(system())
        )
        #        if postOVV_pickle_path.exists():
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
        postOVVout = Load_from_Indexes.MatchPostASTs(postOVVout)
        postOVVout = Load_from_Indexes.MatchECconditions(postOVVout)
        postOVVout.PAR_file = postOVVout.PAR_file.astype(str)

        postOVVout["PAR_date_day"] = [
            pd.datetime.strftime(pd.to_datetime(i), format="%Y-%m-%d")
            for i in postOVVout.PAR_date.fillna(0).values
        ]
        postOVVout = FileHelper.FileOperations.ChangeRoot_DF(
            postOVVout, [], coltype="string"
        )

        postOVVout.to_pickle(postOVV_pickle_path, compression="xz")
        return postOVVout

    def CollectAllExpTypeOVV():
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        today = datetime.today()
        postOVVout = Load_from_Indexes.PreparePostOVV(fastload=False)  # len(22965)
        #        postOVVout.PAR_file = postOVVout.PAR_file.astype(str)

        # === Loading preparation overview of Samples and merging with the data from Characterization techniques === #
        #        SampleCodes = FileHelper.FindExpFolder().LoadSampleCode()
        #        SampleSelect_all = SampleSelection('*','*')
        #        SampleCodesChar = SampleSelect_all.Prep_EA_BET
        SampleCodes = FileHelper.PostChar.SampleCodesChar()
        # pd.merge(SampleCodes,SampleCodesChar,how='left',on='SampleID',suffixes=('','_char')).drop_duplicates(subset=['SampleID','N_content'])
        # === Start preparing pars OVV from index per Experimental type === #
        #        postOVVout,SampleCodes = pd.DataFrame(),pd.DataFrame()
        Reload_set = True
        logger = start_logger()
        EIS_pars = Load_from_Indexes.EIS_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # EIS_Pars2 6745, 22813
        HPRR_pars = Load_from_Indexes.HPRR_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # HPRR 1668
        Cdl_pars = Load_from_Indexes.N2_pars_OVV(
            postOVVout, SampleCodes, reload=Reload_set
        )  # Cdl runs 20322
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
            FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
                "OVV_EIS_missing.xlsx"
            )
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

    @staticmethod
    def IndexPars_CB_paper():
        postOVVout, SampleCodes = pd.DataFrame(), pd.DataFrame()
        PostECddSeries = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
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
        Cdl_pars = Load_from_Indexes.N2_pars_OVV(postOVVout, SampleCodes, reload=False)
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
        Cdl_pars = Load_from_Indexes.N2_pars_OVV(postOVVout, SampleCodes, reload=False)
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
        fail_index_filter = pd.concat(fail_index_gr)
        postOVVout = postOVVout.loc[~postOVVout.index.isin(fail_index_filter.index), :]
        non_uniq = pd.DataFrame(non_uniq_lst)
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
    def MatchECconditions(postOVVout):
        #        postOVVout.postAST.unique()
        #        [(n,len(gr)) for n,gr in postOVVout.groupby('postAST')]
        matchAST_lst = []
        #        'DW16_2018-03-06 00:00:00_no_0.1MHClO4+10mMH2O2_1.0_0.379'
        postOVVout["PAR_date_day"] = [
            pd.datetime.strftime(pd.to_datetime(i), format="%Y-%m-%d")
            for i in postOVVout.PAR_date.fillna(0).values
        ]
        EC_label_cols = [
            "SampleID",
            "pH",
            "Electrolyte",
            "Loading_cm2",
            "postAST",
            "PAR_date_day",
        ]
        post_prev_cols = postOVVout.columns
        #        +[i for i in SampleSelection.EC_exp_cols if i not in ['RPM','Gas']]
        for nAST, ASTgr in postOVVout.groupby(EC_label_cols):
            nAST, ASTgr
            #            for nDT,grDT in ASTgr.groupby(')
            minDT, maxDT = ASTgr.PAR_date.min(), ASTgr.PAR_date.max()
            deltaDT = maxDT - minDT
            #            par_Day = pd.datetime.strftime(nAST[-1],format='%Y-%m-%d')
            EC_exp_query = "_".join([str(i) for i in list(nAST)])
            matchAST_lst.append(
                pd.DataFrame(
                    [(i, EC_exp_query, deltaDT) for i in ASTgr.PAR_file.unique()],
                    columns=["PAR_file", "ECexp", "EC_deltaDT"],
                )
            )

        EC_exp_match = pd.concat(
            [i for i in matchAST_lst], ignore_index=True, sort=False
        )
        postOVVout = pd.merge(postOVVout, EC_exp_match, on=["PAR_file"], how="left")
        print(
            'Added columns: "{0}" to postOVV with len({1})'.format(
                ", ".join(list(set(post_prev_cols) - set(postOVVout.columns))),
                len(postOVVout),
            )
        )
        return postOVVout

    #            ASTgr.SampleID.unique()

    @staticmethod
    def EIS_pars_OVV(postOVVout, SampleCodes, reload=False, extra_plotting=False):
        IndexOVV_EISpars_fn_xls = FileHelper.FindExpFolder(
            "VERSASTAT"
        ).PostDir.joinpath(
            "EIS_pars_IndexOVV_v{0}.xlsx".format(FileHelper.FileOperations.version)
        )
        IndexOVV_EISpars_fn = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "EIS_pars_IndexOVV_v{0}.pkl.compress".format(
                FileHelper.FileOperations.version
            )
        )
        #        EIS_pars_spectra
        if postOVVout.empty or SampleCodes.empty:
            reload = True
        if IndexOVV_EISpars_fn.exists() and reload is not True:
            EIS_pars_char = pd.read_excel(IndexOVV_EISpars_fn, index_col=[0])
            EIS_pars_char = FileHelper.FileOperations.ChangeRoot_DF(
                EIS_pars_char, [], coltype="string"
            )
        else:
            # @@ Read EIS pars files and extend with columns from Samples
            # try other way::   idx_files_EIS = [list(Path(i).rglob('**/EIS/*pars_v20.xlsx')) for i in OnlyRecentMissingOVV.Dest_dir.unique() if list(Path(i).rglob('**/EIS/*pars_v20.xlsx'))]
            columns_containing_EIS = postOVVout.loc[
                postOVVout.Type_output.str.contains("EIS"), "Type_output"
            ].unique()
            #            EIS_pars_index_p1 = postOVVout.query('Type_output == "EIS_Pars1"')
            #            EIS_pars_index_p2 = postOVVout.query('Type_output == "EIS_Pars2"')
            EIS_pars_indexes = postOVVout.query('Type_output == "EIS_Pars"')
            pars_index_from_read = PAR_EIS_fit_V2.EIS_set_index_columns()
            #            EIS_pars_index = pd.concat([EIS_pars_index_p1,EIS_pars_index_p2])
            #            EIS_pars_index = postOVVout.groupby('Type_output').get_group('EIS_Pars1')
            #            EIS_pars_index
            EIS_pars_spectra = (
                postOVVout.groupby("Type_output")
                .get_group("EIS_AllData_combined")
                .drop_duplicates(subset=["PAR_file", "DestFile", "Time_since_run"])
            )
            #            EPtest = EIS_pars_indexes.loc[no_match] # a slice for testing purpose
            #            test_load_nm = no_matches.loc[no_matches[2].str.contains('Columns not matching! "Loading_cm2" values:'),0].values
            #            EPtest = EIS_pars_indexes.loc[EIS_pars_indexes.index.isin(test_load_nm)]
            EISlst, no_match, faillst = [], [], []
            for n, r in EIS_pars_indexes.iterrows():
                try:
                    PF_index, PF_index_stem = Path(r.PAR_file), Path(r.PAR_file).stem
                    index_sourcefn = r.SourceFilename

                    pars = pd.read_excel(index_sourcefn, index_col=[0])
                    pars = FileHelper.FileOperations.ChangeRoot_DF(
                        pars, [], coltype="string"
                    )
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
                            pars.PAR_file.str.contains(
                                "\\\\".join(PF_index.parts[-2::])
                            )
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
                            index_isna_cols = [
                                i for i in mismatch_cols if pd.isna(r[i])
                            ]
                            for col in index_isna_cols:
                                r[col] = pars[col].unique()[0]
                        mismatch_cols = [
                            i for i in overlap_cols if (r[i] != pars[i].unique()[0])
                        ]
                        if any(
                            c in mismatch_cols
                            for c in [
                                "pH",
                                "Electrolyte",
                                "Loading_cm2",
                                "Loading_name",
                            ]
                        ):
                            for col in [
                                i
                                for i in mismatch_cols
                                if i
                                in ["pH", "Electrolyte", "Loading_cm2", "Loading_name"]
                            ]:
                                print(
                                    "changed for {0} from pars {1} to {2} from index for {3}".format(
                                        col,
                                        pars[col].unique()[0],
                                        r[col],
                                        PF_index_stem,
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
                    mismatch_msg = 'Columns not matching! "{1}" values: {0} ,\n Skipped: {2}'.format(
                        *mismatch_values, *mismatch_cols, r.SourceFilename
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
            EIS_pars = ExportECfromCV.make_uniform_EvRHE(EIS_pars)
            EIS_pars_char_mcols = [
                i for i in EIS_pars.columns if i in SampleCodes.columns
            ]
            nonmatching_dtypes = [
                (i, EIS_pars[i].dtype, SampleCodes[i].dtype)
                for i in EIS_pars_char_mcols
                if EIS_pars[i].dtype != SampleCodes[i].dtype
            ]
            nonmt_cls = [i[0] for i in nonmatching_dtypes]
            #            for a,d1,d2 in nonmatching_dtypes:
            #                try:#    SampleCodes[a] = SampleCodes[a].astype(d1)
            #                except:#    SampleCodes[a].fillna(value=0).str.replace(',','.').astype(d1)

            EIS_pars_char = pd.merge(
                EIS_pars.drop(columns=nonmt_cls),
                SampleCodes.drop(columns=nonmt_cls),
                on=[i for i in EIS_pars_char_mcols if i not in nonmt_cls],
                how="left",
            )
            EIS_pars_char = EIS_pars_char.loc[
                EIS_pars.Model_EEC != "Model(Singh2015_RQR)"
            ]
            EIS_pars_char.to_pickle(IndexOVV_EISpars_fn)
            #            new_IndexOVV_EISpars_target = FileHelper.FileOperations.CompareHashDFexport(EIS_pars_char,IndexOVV_EISpars_fn)
            try:
                logger.info(
                    "PostEC EIS re-indexed and saved: {0}".format(
                        new_IndexOVV_EISpars_target
                    )
                )
            except:
                pass

        #            EIS_pars_char.query('pH < 17').groupby('Model_EEC').plot(y='RedChisqr',x='E_RHE',colormap='viridis',kind='scatter',yscale='log')
        if extra_plotting == True:
            fast_checking_EEC_models = [
                "Model(Singh2015_RQRQR)",
                "Model(Singh2015_RQRWR)",
                "Model(Singh2015_R3RQ)",
                "Model(Bandarenka_2011_RQRQR)",
            ]
            for modname in fast_checking_EEC_models[0:1]:
                EIS_pars_char.loc[EIS_pars_char["Model_EEC"] == modname].query(
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
                EIS_pars_char.loc[EIS_pars_char["Model_EEC"] == modname].query(
                    "pH < 7"
                ).plot(
                    y="nAd",
                    x="E_RHE",
                    c="BET_cat_agg",
                    colormap="rainbow_r",
                    kind="scatter",
                    ylim=(0, 1),
                    title=modname,
                )
                EIS_pars_char.loc[EIS_pars_char["Model_EEC"] == modname].query(
                    "pH < 7"
                ).plot(
                    y="nDL",
                    x="E_RHE",
                    c="BET_cat_agg",
                    colormap="rainbow_r",
                    kind="scatter",
                    ylim=(0, 1),
                    title=modname,
                )
                EIS_pars_char.loc[EIS_pars_char["Model_EEC"] == modname].query(
                    "pH < 7"
                ).plot(
                    y="Rct",
                    x="E_RHE",
                    c="BET_cat_agg",
                    colormap="rainbow_r",
                    kind="scatter",
                    ylim=(1, 1e6),
                    logy=True,
                    title=modname,
                )
                EIS_pars_char.loc[EIS_pars_char["Model_EEC"] == modname].query(
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
            fig, ax = plt.subplots()
            EIS_pars_char.query("pH < 17").groupby("Model_EEC").plot(
                y="RedChisqr", x="E_RHE", colormap="viridis", kind="scatter", ax=ax
            )
            for n, Hgr in EIS_pars_char.query("pH < 7").groupby("postAST"):
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
            plt.close()
        today = datetime.now().date()
        eis_daily_pickle_path = PostDestDir.joinpath(
            f"{today.year}-{today.month}-{today.day}_EIS_pars_{system()}.pkl.compress"
        )
        EIS_pars.to_pickle(eis_daily_pickle_path)
        print(f"EIS to daily pickle: {eis_daily_pickle_path}")
        return EIS_pars_char

    @staticmethod
    def HPRR_pars_OVV(postOVVout, SampleCodes, reload=False):
        #        exp_type = 'H
        IndexOVV_HPRRpars_fn = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_HPRR_v{0}.xlsx".format(FileHelper.FileOperations.version)
        )
        if IndexOVV_HPRRpars_fn.exists() and reload != True:
            HPRR_pars_char = pd.read_excel(IndexOVV_HPRRpars_fn, index_col=[0])
            HPRR_pars_char = FileHelper.FileOperations.ChangeRoot_DF(
                HPRR_pars_char, [], coltype="string"
            )
        else:
            # === Making destination directories === #
            PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
                "PostEC"
            )
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
            HPRR_pars_raw = FileHelper.FileOperations.ChangeRoot_DF(
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
            new_IndexOVV_HPRRpars_target = (
                FileHelper.FileOperations.CompareHashDFexport(
                    HPRR_pars_char, IndexOVV_HPRRpars_fn
                )
            )
            logger.info(
                "PostEC HPRR re-indexed and saved: {0}".format(
                    new_IndexOVV_HPRRpars_target
                )
            )
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
    def HER_pars_OVV(postOVVout, SampleCodes, reload=False):
        #        exp_type = 'H
        IndexOVV_HER_pars_fn = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_HER_v{0}.xlsx".format(FileHelper.FileOperations.version)
        )
        if postOVVout.empty or SampleCodes.empty:
            reload = False

        if IndexOVV_HER_pars_fn.exists() and reload != True:
            HER_pars_char = pd.read_excel(IndexOVV_HER_pars_fn, index_col=[0])
            HER_pars_char = FileHelper.FileOperations.ChangeRoot_DF(
                HER_pars_char, [], coltype="string"
            )

        else:
            # === Making destination directories === #
            PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
                "PostEC"
            )
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
            HER_pars_raw = FileHelper.FileOperations.ChangeRoot_DF(
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
            new_IndexOVV_HERpars_target = FileHelper.FileOperations.CompareHashDFexport(
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
    def OER_pars_OVV(postOVVout, SampleCodes, reload=False):
        #        exp_type = 'H
        IndexOVV_OERpars_fn = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_OER_v{0}.xlsx".format(FileHelper.FileOperations.version)
        )
        if IndexOVV_OERpars_fn.exists() and reload != True:
            OER_pars_char = pd.read_excel(IndexOVV_OERpars_fn, index_col=[0])
            OER_pars_char = FileHelper.FileOperations.ChangeRoot_DF(
                OER_pars_char, [], coltype="string"
            )

        else:
            # === Making destination directories === #
            PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
                "PostEC"
            )
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
            OER_pars_raw = FileHelper.FileOperations.ChangeRoot_DF(
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
            new_IndexOVV_OERpars_target = FileHelper.FileOperations.CompareHashDFexport(
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
    def ORR_pars_OVV(postOVVout, SampleCodes, reload=False):
        #        exp_type = 'H
        IndexOVV_ORRpars_fn = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_ORR_v{0}.pkl.compress".format(
                FileHelper.FileOperations.version
            )
        )

        if IndexOVV_ORRpars_fn.exists() and reload is not True:
            ORR_pars_char = pd.read_excel(IndexOVV_ORRpars_fn, index_col=[0])
            ORR_pars_char = FileHelper.FileOperations.ChangeRoot_DF(
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
            ORR_pars = FileHelper.FileOperations.ChangeRoot_DF(
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
                (n, *FileHelper.FindSampleID.ink_loading_from_filename(i.PAR_file))
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
            #            IndexOVV_ORRpars_fn = FileHelper.FileOperations.CompareHashDFexport(ORR_pars_char,IndexOVV_ORRpars_fn)
            ORR_pars_char.to_pickle(IndexOVV_ORRpars_fn)
            logger.info(
                "PostEC ORR re-indexed and saved: {0}".format(IndexOVV_ORRpars_fn)
            )
        #        ORR_pars_char.query('(pH < 7) & (RPM > 900)').plot(y='Jkin_075',x='AD/AG',kind='scatter',logy=False)
        fig, ax = plt.subplots()
        ORR_pars_char.query("(pH < 15) & (RPM > 900)").plot(
            y="J_diff_lim",
            x="Loading_cm2",
            kind="scatter",
            logy=0,
            c="BET_cat_agg",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars_char.query("(pH < 14) & (RPM > 900)").plot(
            x="TSb_l",
            y="E_onset",
            kind="scatter",
            xlim=(0.5, 1),
            ylim=(0.5, 1),
            c="pH",
            colormap="rainbow_r",
            ax=ax,
        )
        #        plot(y='J_diff_lim',x='Loading_cm2',kind='scatter',logy=0,c='BET_cat_agg',colormap='viridis',ax=ax)
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars_char.query("(pH < 14) & (RPM > 900)").plot(
            x="TSa_l",
            y="E_onset",
            kind="scatter",
            xlim=(10, 200),
            ylim=(0.5, 1),
            c="pH",
            colormap="rainbow_r",
            ax=ax,
        )
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars_char.query("(pH < 7) & (RPM > 900)").plot(
            y="J_diff_lim",
            x="Loading_cm2",
            kind="scatter",
            logy=0,
            c="BET_cat_agg",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars_char.query(
            '(pH < 15) & (RPM > 1100) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        ).plot(
            y="E_onset",
            x="N_content",
            kind="scatter",
            logy=0,
            c="pH",
            colormap="viridis",
            ax=ax,
            xlim=(0, 10),
        )
        plt.show()
        plt.close()
        fig, ax = plt.subplots()
        ORR_pars_char.query("(pH < 15) & (RPM > 900)").plot(
            y="E_half",
            x="AD/AG",
            kind="scatter",
            logy=0,
            c="pH",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()
        fig, ax = plt.subplots()
        ORR_pars_char.query("(pH < 7) & (RPM > 900)").plot(
            y="Jkin_075",
            x="BET_cat_agg",
            kind="scatter",
            logy=0,
            c="N_content",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()
        fig, ax = plt.subplots()
        ORR_pars_char.query("(pH < 7) & (RPM > 1200)").plot(
            y="Jkin_075",
            x="D1_pop_Fe_wt",
            kind="scatter",
            logy=0,
            c="N_content",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()
        return ORR_pars_char

    @staticmethod
    def N2_pars_OVV(postOVVout, SampleCodes, reload=False):
        #        exp_type = 'H
        IndexOVV_N2_pars_fn = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_CdlN2_v{0}.xlsx".format(FileHelper.FileOperations.version)
        )
        if IndexOVV_N2_pars_fn.exists() and reload != True:
            Cdl_pars_char = pd.read_excel(IndexOVV_N2_pars_fn, index_col=[0])
            Cdl_pars_char = FileHelper.FileOperations.ChangeRoot_DF(
                Cdl_pars_char, [], coltype="string"
            )
        else:
            # @@ Check POST_AST status from OVV and PRM

            N2Cdl_pars_index = postOVVout.groupby("Type_output").get_group(
                "N2_Cdl_pars"
            )
            N2Cdl_pars_files = [
                Path(i)
                for i in N2Cdl_pars_index["SourceFilename"].unique()
                if re.search("(?i)(_pars|_v20)", Path(i).stem) and Path(i).exists()
            ]
            #            cdl = pd.read_excel(N2Cdl_pars_files[0],index_col=[0])
            N2Cdl_pars_raw = pd.concat(
                [
                    pd.read_excel(i, index_col=[0]).drop_duplicates()
                    for i in N2Cdl_pars_files
                ],
                sort=False,
            )
            N2Cdl_pars_raw.drop_duplicates(
                subset=N2Cdl_pars_raw.columns[0:11], inplace=True
            )
            N2Cdl_pars_raw = FileHelper.FileOperations.ChangeRoot_DF(
                N2Cdl_pars_raw, [], coltype="string"
            )
            #                                                                     [i for i in N2Cdl_pars_raw.columns if re.search('([F-f]ile|ORR_act)',i)]
            #            N2Cdl_pars.rename(columns={'Filename' : 'PAR_file'})
            #            EPtest = N2Cdl_pars_index.loc[no_match] # a slice for testing purpose
            #            pd.merge(N2Cdl_pars_raw,N2_CV_index[['PAR_file','DestFile']],on='PAR_file',how='left')
            N2Cdl_data_index = postOVVout.groupby("Type_output").get_group(
                "N2_Cdl_data"
            )
            N2_CV_index = postOVVout.groupby("Type_output").get_group("N2_CV")
            lst, no_match, non_exist = [], [], []
            for n, r in N2Cdl_pars_raw.iterrows():
                Cdl_data_file = N2Cdl_data_index.loc[
                    N2Cdl_data_index.PAR_file == r.PAR_file
                ].DestFile.unique()
                CV_files = N2_CV_index.loc[
                    N2_CV_index.PAR_file == r.PAR_file
                ].DestFile.unique()
                lst.append([set(Cdl_data_file), set(CV_files)])

            if len(N2Cdl_pars_raw) == len(lst):
                N2Cdl_pars_raw = N2Cdl_pars_raw.assign(
                    **{
                        "Cdl_data_file": [i[0] for i in lst],
                        "Cdl_CV_data_files": [i[1] for i in lst],
                    }
                )
            #            Cdl_pars = pd.concat([i for i in lst],sort=False,ignore_index=True)
            Cdl_pars = ExportECfromCV.make_uniform_EvRHE(N2Cdl_pars_raw)
            Cdl_pars = pd.merge(
                Cdl_pars,
                N2Cdl_pars_index[["PAR_file", "ECexp", "preAST", "EC_deltaDT"]],
                on="PAR_file",
                how="left",
            )
            Cdl_pars.drop_duplicates(subset=Cdl_pars.columns[0:11], inplace=True)
            Cdl_pars_merge_cols = [
                i
                for i in Cdl_pars.columns
                if i in SampleCodes.columns and not "Unnamed" in i
            ]
            Cdl_pars_char = pd.merge(
                Cdl_pars, SampleCodes, on=Cdl_pars_merge_cols, how="left"
            )
            Cdl_pars_char.drop_duplicates(
                subset=Cdl_pars_char.columns[0:11], inplace=True
            )
            new_N2_pars_char_target = FileHelper.FileOperations.CompareHashDFexport(
                Cdl_pars_char, IndexOVV_N2_pars_fn
            )
            logger.info(
                "PostEC Cdl N2 CVs re-indexed and saved: {0}".format(
                    new_N2_pars_char_target
                )
            )
        Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)').plot(
            y="Cdl",
            x="E_RHE",
            c="BET_cat_agg",
            colormap="viridis",
            kind="scatter",
            ylim=(0, 0.08),
            title="Cdl in acid",
        )
        #        Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)').groupby('BET_cat_agg').plot(y='Cdl',x='E_RHE',colormap='viridis',kind='scatter',ylim=(0,0.08),title='Cdl in acid')
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
        for n, Ogr in Cdl_atE.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)').groupby(
            "postAST"
        ):
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
        for n, Ogr in Cdl_atE.query('(Sweep_Type_N2 == "cathodic") & (pH > 7)').groupby(
            "postAST"
        ):
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

    @staticmethod
    def HPRR_scans():
        # === Loading Index overview files from folders === #
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        LoadOVV = ExportECfromCV(reload=True)
        postOVVout = LoadOVV.postOVVout
        #        print('Types:',' , '.join([str(i) for i in postOVVout.Type_output.unique()]))
        postOVVout.PAR_file = postOVVout.PAR_file.astype(str)

        N2_CV = postOVVout.groupby("Type_output").get_group("N2_CV")
        N2_CV.query(
            '(RPM_HPRR > 700) & (Loading_cm2 > 0.1) &  (E_name == "E_j0")'
        ).plot(x="BET_cat_agg", y="fit_slope_HPRR", kind="scatter")

        ORR_pars_index = postOVVout.groupby("Type_output").get_group(
            "ORR_Jkin_calc_Pars"
        )
        ORR_index_data = postOVVout.groupby("Type_output").get_group(
            "ORR_Jkin_calc_RRDE"
        )
        # @@ Check POST_AST status from OVV and PRM
        ORR_Pars_files = [
            i
            for i in ORR_pars_index["SourceFilename"].unique()
            if re.search("(?i)(_pars|_v20)", i.stem) and i.is_file()
        ]
        ORR_pars = pd.concat(
            [pd.read_excel(i, index_col=[0]) for i in ORR_Pars_files]
        ).rename(columns={"File": "PAR_file"})
        ORR_pars.PAR_file = ORR_pars.PAR_file.astype(str)
        ORR_pars_index.PAR_file = ORR_pars_index.PAR_file.astype(str)
        ORR_merge_cols = [
            i
            for i in ORR_pars.columns
            if i in ORR_pars_index.columns and not "Segment" in i
        ]
        p2, ovv2 = ORR_pars.set_index(ORR_merge_cols), ORR_pars_index.set_index(
            ORR_merge_cols
        )
        ORR_pars_ovv = p2.join(ovv2, rsuffix="_ovv").drop_duplicates()
        #        ORR_pars_ovv.query('(RPM > 1400) & ( RPM < 1600) & (Loading_cm2 > 0.1)').plot(y='J_diff_lim',x='Loading_cm2',kind='scatter')
        ORR_pars_ovv.query("(RPM > 100)").plot(
            y="E_onset", x="Loading_cm2", kind="scatter"
        )
        acid_BoL = ORR_pars_ovv.query(
            '(RPM > 1400) & ( RPM < 1600) & (Loading_cm2 > 0.1) & (pH < 7) & (postAST == "no")'
        )

        Jdiff_lim_6 = (
            ORR_pars_ovv.query("(J_diff_lim < -6)")
            .sort_values(by="J_diff_lim")
            .reset_index()
        )

        acid_BoL.plot(y="FracH2O2_050", x="Loading_cm2", kind="scatter")
        DW21 = ORR_pars_ovv.query(
            '(RPM > 900) & ( RPM < 1600) & (Loading_cm2 > 0.1) & (pH < 7) & (SampleID == "DW21")'
        )

        #            ORR_pars_ovv.query('(RPM > 900) & ( RPM < 1600) & (Loading_cm2 > 0.1) & (pH < 7)').plot(y=col,x='Loading_cm2',kind='scatter',label=col,s=100,c='g')
        #            ORR_pars_ovv.query('(RPM > 900) & ( RPM < 1600) & (Loading_cm2 > 0.1) & (SampleID == "DW21") & (pH < 7)').plot(y=col,x='Loading_cm2',kind='scatter',label=col,s=100,c='g')

        ORR_pars_ovv.query(
            '(RPM > 900) & ( RPM < 1600) & (Loading_cm2 > 0.1) & (SampleID == "DW21")'
        ).plot(y="J_diff_lim", x="Loading_cm2", kind="scatter")
        DW21.plot(y="Jkin_075", x="Loading_cm2", kind="scatter")
        DW21.plot(y="Jkin_080", x="Loading_cm2", kind="scatter")
        DW21.plot(y="FracH2O2_050", x="Loading_cm2", kind="scatter")
        #        ORR_pars_ovv.query('(RPM_HPRR > 700) & (Loading_cm2 > 0.1) &  (E_name == "E_j0")').plot(y='J_diff_lim',x='fit_slope_HPRR',kind='scatter')
        #        ORR_pars_ovv = pd.merge(ORR_pars,ORR_pars_index,on=ORR_merge_cols,suffixes=('','_ovv'),how='left')
        #        ORR_pars = pd.merge(ORR_pars,postOVVout,on=['PAR_file','SampleID','Electrolyte','pH','postAST'],how='left',suffixes=('','_ovv'))
        print(
            "Leftover SampleIDs: {0}".format(
                set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())
            )
        )
        ORR_pars = pd.merge(ORR_pars, SampleCodes, on="SampleID", how="left")
        HP_ORR_pars = pd.merge(HPRR_pars, ORR_pars, on="SampleID", how="left")

        #            Multidx = pd.MultiIndex.from_tuples([r.to_list()]*len(pars) ,names=r.index)
        #            parsIDX = pars.set_index(Multidx)

        #        EIS_Pars = pd.merge(EIS_Pars,postOVVout,on=['PAR_file','SampleID','Electrolyte','pH','Gas'],how='left',suffixes=('','_ovv'))
        #        EIS_Pars = pd.merge(EIS_Pars,postOVVout,on=['PAR_file','SampleID'],how='left',suffixes=('','_ovv'))
        #        ORR_data_files = [i for i in ORR_index_data['SourceFilename'].unique()]
        print(
            "Leftover SampleIDs: {0}".format(
                set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())
            )
        )
        ORR_pars = pd.merge(ORR_pars, SampleCodes, on="SampleID", how="left")

    @staticmethod
    def plot_ECkin_Loading(ORR_pars_ovv, X_col="Loading_cm2"):

        X_col = "pH"
        if "Loading" in X_col:
            acid_BoL = ORR_pars_ovv.query(
                '(RPM > 1400) & ( RPM < 1600) & (Loading_cm2 > 0.1) & (pH < 1) & (PAR_date > 20180101) & (postAST == "no")'
            )
        elif "pH" in X_col:
            acid_BoL = ORR_pars_ovv.query(
                '(RPM > 900) & ( RPM < 1600) & (Loading_cm2 > 0.3) & (pH <= 14)  & (Loading_cm2 < 0.45) & (PAR_date > 20180101) & (postAST == "no")'
            ).reset_index()
        acid_BoL = acid_BoL.drop_duplicates()
        DW28_normal = acid_BoL.query(
            '(SampleID == "DW28") & (Loading_cm2 < 0.4) & (Loading_cm2 > 0.3)'
        )
        JOS15_normal = acid_BoL.query(
            '(SampleID == "JOS15") & (Loading_cm2 < 0.4) & (Loading_cm2 > 0.3)'
        )
        PostEC_plotdir = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "ECkin_corr_plots"
        )
        PostEC_plotdir.mkdir(parents=True, exist_ok=True)
        acid_BoL.to_excel(PostEC_plotdir.joinpath("corr_data_{0}.xlsx".format(X_col)))
        for col in SampleSelection.EC_ORR_kin_par_cols:
            fig, ax = plt.subplots(figsize=(8, 8))
            plot_path = PostEC_plotdir.joinpath("corr_{0}_{1}".format(X_col, col))

            for sID, sgr in acid_BoL.groupby("SampleID"):
                if sgr[X_col].nunique() > 1 and len(sgr[col].dropna()) > 1:
                    x, y = sgr[X_col].values, sgr[col].values
                    res = smf.ols(
                        "sgr.{0} ~ sgr.{1}".format(col, X_col), data=sgr
                    ).fit()
                    poly_res = np.polyfit(x, y, 2, full=True)
                    intercept, slope, rsquared = (
                        res.params[0],
                        res.params[1],
                        res.rsquared,
                    )
                    x_lin = np.arange(sgr[X_col].min(), sgr[X_col].max() + 0.1, 0.05)
                    y_lin, y_poly = x_lin * slope + intercept, np.poly1d(poly_res[0])
                    #                    sgr.plot(y=col,x='Loading_cm2',kind='scatter',label=sID,ax=ax)
                    ax.scatter(
                        sgr[X_col].values,
                        sgr[col].values,
                        label="{0}_{1:.2f}".format(sID, rsquared),
                    )
                    #                    ax.plot(x_lin,y_lin,lw=2, ms=10,alpha=rsquared+0.1)
                    print(poly_res)
                    ax.plot(
                        x_lin, np.poly1d(poly_res[0])(x_lin), lw=2, ms=10, alpha=0.5
                    )
            ax.set_xlabel(X_col)
            ax.set_ylabel(col)
            ax.legend(loc="best", ncol=2, fontsize=10)
            plt.tight_layout()
            plt.savefig(plot_path.with_suffix(".png"), dpi=300)
            plt.close()

    @staticmethod
    def plot_ECkin_AST(ORR_pars_ovv, X_col="postAST"):
        #        idx = pd.IndexSlice
        #        ORR_pars_ovv.loc[idx[:,:,:,:,:]]
        #        AST_Samples = ORR_pars_ovv.query('postAST != "no"').reset_index().SampleID.unique()
        AST_samples = ORR_pars_ovv.loc[
            (
                ORR_pars_ovv.query('postAST != "no"')
                .index.get_level_values("SampleID")
                .unique()
            ),
            :,
        ]
        acid_AST = AST_samples.query(
            "(pH < 7) & (RPM > 900) & ( RPM < 1600) & (Loading_cm2 > 0.3) & (Loading_cm2 < 0.45) & (PAR_date > 20180101)"
        )
        alk_AST = AST_samples.query(
            "(pH > 7) & (RPM > 900) & ( RPM < 1600) & (Loading_cm2 > 0.3) & (Loading_cm2 < 0.45) & (PAR_date > 20180101)"
        )

        PostEC_plotdir = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "ECkin_corr_plots/AST"
        )
        PostEC_plotdir.mkdir(parents=True, exist_ok=True)
        for col in SampleSelection.EC_ORR_kin_par_cols:

            plot_path = PostEC_plotdir.joinpath("corr_{0}_{1}".format(X_col, col))
            for sID, sgr in acid_AST.groupby("SampleID"):
                fig, ax = plt.subplots(figsize=(8, 8))
                if sgr.index.get_level_values(X_col).nunique() > 1:
                    sgr.reset_index().plot.bar(x="postAST", y=col, title=sID, ax=ax)
                    x, y = sgr[X_col].values, sgr[col].values
                    res = smf.ols(
                        "sgr.{0} ~ sgr.{1}".format(col, X_col), data=sgr
                    ).fit()
                    poly_res = np.polyfit(x, y, 2, full=True)
                    intercept, slope, rsquared = (
                        res.params[0],
                        res.params[1],
                        res.rsquared,
                    )
                    x_lin = np.arange(sgr[X_col].min(), sgr[X_col].max() + 0.1, 0.05)
                    y_lin, y_poly = x_lin * slope + intercept, np.poly1d(poly_res[0])
                    #                    sgr.plot(y=col,x='Loading_cm2',kind='scatter',label=sID,ax=ax)
                    ax.scatter(
                        sgr[X_col].values,
                        sgr[col].values,
                        label="{0}_{1:.2f}".format(sID, rsquared),
                    )
                    #                    ax.plot(x_lin,y_lin,lw=2, ms=10,alpha=rsquared+0.1)
                    print(poly_res)
                    ax.plot(
                        x_lin, np.poly1d(poly_res[0])(x_lin), lw=2, ms=10, alpha=0.5
                    )
            ax.set_xlabel(X_col)
            ax.set_ylabel(col)
            ax.legend(loc="best", ncol=2, fontsize=10)
            plt.tight_layout()
            plt.savefig(plot_path.with_suffix(".png"), dpi=300)
            plt.close()

        if len(list(PPD_HPRR_data.rglob("*"))) < len(HPRR_pars.DataFile.unique()):
            for sIDcode, sCgrp in HPRR_pars.groupby(["SampleID"]):
                for DTfile, DTgr in sCgrp.groupby("DataFile"):
                    DTfilePD = pd.read_excel(DTfile, index_col=[0])
                    rpm, date = DTgr.RPM_HPRR.iloc[
                        0
                    ], FileHelper.FindSampleID.Determine_date_from_filename(
                        DTfile
                    ).strftime(
                        "%Y-%m-%d"
                    )
                    swptp = DTgr.Sweep_Type_HPRR.unique()[0]
                    DT_dest = PPD_HPRR_data.joinpath(
                        "{0}_{1}_{2}_{3}_{4}.xlsx".format(
                            swptp,
                            rpm,
                            sIDcode,
                            date,
                            Path(DTfile).stem.replace(
                                "{0}_{1}_".format(swptp, rpm), ""
                            ),
                        )
                    )
                    DTfilePD.to_excel(DT_dest)

        HPRR_pars_1500 = HPRR_pars.query("RPM_HPRR > 1400")
        ORR_pars_1500 = ORR_pars.query("(RPM > 1400) & ( RPM < 1600)")
        #        & SeriesID == "Porph_SiO2"
        HPRR_pars_1500_grp = HPRR_pars_1500.groupby(["Sweep_Type_HPRR", "E_name"])
        HPRR_pars_1500 = HPRR_pars_1500.assign(**{"all": "all"})

        #%%
        for nm, gr1500 in HPRR_pars_1500.groupby(["Sweep_Type_HPRR", "E_name"]):
            for SerID, SerGrp in gr1500.groupby(["all"]):
                Y_cols = ["fit_slope_HPRR", "E_fit_HPRR"]
                Y_lbls = ["slope HPRR / dJ/dE", "E / V v RHE"]
                Y_units = dict(zip(Y_cols, Y_lbls))
                smpllst = [
                    "SampleCode",
                    "SampleID",
                    "MetalPrec",
                    "SeriesID",
                    "ReferenceID",
                    "SizeNP",
                    "ML",
                    "C_N",
                    "Colorcode",
                    "Metal_wt",
                    "IndividualLabel",
                ]
                PPD_HPRR_data_Serie = PPD_HPRR_data.joinpath(SerID)
                PPD_HPRR_data_Serie.mkdir(parents=True, exist_ok=True)

                if len(SerGrp) < 2:
                    continue
                gr_sID = SerGrp.groupby(["SampleCode"])

                gr_sID_desc = gr_sID.describe(include="all")
                slice_cols = [
                    a
                    for a in gr_sID_desc.columns
                    if "mean" in a or "top" in a or "std" in a
                ] + [
                    i
                    for i in gr_sID_desc.columns
                    if ("count" in i) and gr_sID_desc[i].sum() < 1e-10
                ]
                slice_cols.sort()
                gr_sID_desc_Codes = gr_sID_desc[slice_cols].dropna(axis=1, how="all")
                #                    gr_sID_desc_Codes = pd.merge(gr_sID_desc,SampleCodes.set_index('SampleCode'),left_index=True, right_index=True,how='inner')
                idx_len = [(i, len(i)) for i in gr_sID_desc_Codes.index]
                max_len = np.max([i[1] for i in idx_len])
                newlabels = [i.ljust(max_len, " ") for i in gr_sID_desc_Codes.index]
                gr_sID_desc_Codes.index = newlabels
                gr_sID_desc_Codes.to_excel(
                    PPD_HPRR_data_Serie.joinpath(
                        "{0}_{1}_{2}.xlsx".format(SerID, nm[0], nm[1])
                    )
                )

                corr_method, corr_cutoff = (
                    "spearman",
                    0.5,
                )  # or default is pearson spearman
                #                rcorr = gr_sID_desc_Codes[SampleSelection.InterestingCols+SampleSelection.RAMAN_cols_corr].corr(method=corr_method)
                gr_Corr = gr_sID_desc_Codes[
                    [i for i in gr_sID_desc_Codes.columns if "mean" in i]
                ]
                rcorr = gr_Corr.corr(method=corr_method)
                corr_triu = rcorr.where(np.tril(np.ones(rcorr.shape)).astype(np.bool))
                scorr_triu = corr_triu.stack()
                filter_corr = scorr_triu[
                    (scorr_triu.abs() > corr_cutoff) & (scorr_triu.abs() < 0.97)
                ]
                #                for n,i in filter_corr.iteritems():
                #                    print(n,i)
                #                    try:
                #                        plot_pd_SampleIDs(gr_sID_desc_Codes,n[0],n[1],PPD_HPRR_data_Serie)
                #                    except:
                #                        pass
                fcorr = filter_corr.unstack()
                fcorr_sliced = fcorr[Y_cols].dropna(how="all")
                if not fcorr_sliced.empty:
                    plt.subplots(figsize=(10, 20))
                    sns.heatmap(fcorr_sliced)
                    heatmap_title = "{0}_{1}_heatmap_{2}_{3}".format(
                        corr_method, nm[0], nm[1], int(corr_cutoff * 100)
                    )
                    plt.suptitle(heatmap_title)
                    plt.savefig(
                        PPD_HPRR_data_Serie.joinpath(heatmap_title).with_suffix(".png"),
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close()
                # === Plot for different metal wt% or ML coverage ===
                PPD_HPRR_data_ser_scatter = PPD_HPRR_data_Serie.joinpath(
                    "scatter_{0}".format(corr_method)
                )
                PPD_HPRR_data_ser_scatter.mkdir(parents=True, exist_ok=True)
                for Y_type in Y_cols:
                    y_lbl = Y_units[Y_type]
                    if "Tafel" in nm[1] and "fit_slope_HPRR" in Y_type:
                        y_lbl = "Tafel slope HPRR / mV dec^-1"

                    corr_head = (
                        fcorr_sliced.sort_values(by=[(Y_type, "mean") for i in Y_cols])
                        .head(7)
                        .index
                    )
                    corr_tail = (
                        fcorr_sliced.sort_values(by=[(Y_type, "mean") for i in Y_cols])
                        .tail(7)
                        .index
                    )

                    X_cols = (
                        ["ML", "Metal_wt", "BET_cat_agg", "N_content"]
                        + [i[0] for i in corr_head]
                        + [i[0] for i in corr_tail]
                    )
                    for X_type in X_cols:

                        fig, ax = plt.subplots(1, figsize=(12, 12))
                        ax.set_ylabel(y_lbl)
                        ax.set_xlabel(X_type)

                        for sID, sgrX in gr_sID_desc_Codes.iterrows():
                            #                    marker_style = dict(color='tab:blue', linestyle=':', marker='o',
                            #                    markersize=15, markerfacecoloralt='tab:red')
                            marker_style = (
                                Characterization_TypeSetting.Scatter_Marker_Style(
                                    sgrX[("MetalPrec", "top")], serie=SerID
                                )
                            )
                            ax.scatter(
                                sgrX[(X_type, "mean")],
                                sgrX[(Y_type, "mean")],
                                **marker_style,
                                label=sgrX.SeriesID,
                            )
                            ax.set_yscale("linear")
                        #                    ax.legend(True)
                        X_filename = "".join([i for i in X_type if i.isalnum()])
                        plt.savefig(
                            PPD_HPRR_data_ser_scatter.joinpath(
                                "{0}_{1}_{2}_{3}_{4}.png".format(
                                    SerID, nm[0], nm[1], Y_type, X_filename
                                )
                            ),
                            bbox_inches="tight",
                        )
                        plt.close()
                    #            nm,gr
                    #            sortgr.sort_values(by='fit_slope_HPRR',inplace=True)
                    #                for i in idx_len:
                    #                    diff = max_len-i[1]
                    #                    if diff > 0:
                    #                        newlabels.append(str(i[0]+[' ']*diff))
                    #                    else:
                    #                        newlabels.append(i)
                    #            fig,ax = plt.subplots(1,2,figsize=(12,12))
                    PPD_HPRR_data_ser_bar = PPD_HPRR_data_Serie.joinpath("bar_plots")
                    PPD_HPRR_data_ser_bar.mkdir(parents=True, exist_ok=True)
                    gr_sID_desc_Codes.sort_values(by=[(Y_type, "mean")], inplace=True)
                    fig, ax = plt.subplots(1, figsize=(12, 12))

                    for sID, sgr in gr_sID_desc_Codes.iterrows():
                        #                    sID_Char = SampleCodes.query('SampleID == @sID')
                        #                for sData,sDgr in gr.groupby('DataFile'):
                        #                    data = pd.read_excel(sData,index_col=[0])
                        #                    data = data.query('Sweep_Type == @nm[0]')
                        #                           desc = gr.groupby('SampleID').describe()
                        marker_style = (
                            Characterization_TypeSetting.Scatter_Marker_Style(
                                sgr[("MetalPrec", "top")], serie=SerID
                            )
                        )
                        ax.bar(
                            sID,
                            sgr[(Y_type, "mean")],
                            yerr=sgr[(Y_type, "std")],
                            color=marker_style["color"],
                        )
                    #                    y_col = [i for i in data.columns if nm[1] in i]
                    #                    if not y_col:
                    #                        y_col = [i for i in data.columns if nm[1].split('_')[-1] in i]
                    #                    data.plot(x='E_Applied_VRHE',y=y_col,ls='--',ax=ax[0],label='')
                    #                    if 'Tafel' in y_col[0]:
                    #                        exp_data = 'log_Abs_jmAcm-2'
                    #                    else:
                    #                        exp_data = 'jmAcm-2'
                    #                    data.plot(x='E_Applied_VRHE',y=exp_data,ax=ax[0],label=sID)
                    plt.title("{0} sweep at {1} for {2} ".format(nm[0], nm[1], Y_type))
                    plt.ylabel(y_lbl)
                    plt.xticks(rotation=90, ha="center")
                    plt.savefig(
                        PPD_HPRR_data_ser_bar.joinpath(
                            "{0}_{1}_{2}_{3}.png".format(SerID, nm[0], nm[1], Y_type)
                        ),
                        bbox_inches="tight",
                    )
                    plt.close()
        #%%
        #                desc.plot(y=[('fit_slope_HPRR','mean')],kind='bar',title='{0} sweep for {1}'.format(nm[0],nm[1]))
        gr = HPRR_pars_1500_grp.get_group(("cathodic", "E_j0"))
        desc = gr.groupby("SampleID").describe()

    # list(Path('F:\EKTS_CloudStation\CloudStation\Experimental data\Organized_Data\VERSASTAT\PostEC\EIS\0.1MH2SO4').rglob('*'))
    def ORR_filter_Sweep(plot_Filter=False, **kwargs):
        EvRHE = "E_AppV_RHE"
        postOVVout = PostEC.LoadPostOVV()
        if "take_tail" in kwargs:
            postOVVout = postOVVout.sort_values(by="EXP_date").tail(300)
        if "take_series" in kwargs:
            postOVVout = postOVVout.loc[
                postOVVout["SeriesID"] == kwargs["take_series"], :
            ]

        for n, r in postOVVout.query('Type_Exp == "ORR_RRDE"').iterrows():
            ORR_PDD = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
                "EC_files_out", "ORR_RRDE", r.Electrolyte
            )
            #            ORR_PDD = FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PostEC'
            ORR_PDD.mkdir(parents=True, exist_ok=True)
            TF = pd.read_excel(r.SourceFilename)
            print(n, r.DestFilename)
            SwpCats = CategoricalDtype(["anodic", "cathodic"], ordered=False)
            TF[[EvRHE, "J_N2_scan"]].to_excel(
                ORR_PDD.joinpath("N2_" + Path(r.DestFilename).stem).with_suffix(".xlsx")
            )
            TF = TF.dropna(subset=["J_ring", "Jcorr"])
            if not TF.empty:
                TF["EvRHE_diff"] = TF[EvRHE].diff()
                TF = TF.assign(
                    **{
                        "Sweep_Type": np.where(
                            TF.EvRHE_diff > 0,
                            "anodic",
                            np.where(TF.EvRHE_diff < 0, "cathodic", np.nan),
                        )
                    }
                )
                TF.Sweep_Type = TF.Sweep_Type.astype(SwpCats).fillna(method="bfill")
                TF = TF.iloc[5:-5]
                Fltr_Win = 71
                for Swp, swGr in TF.groupby(by="Sweep_Type"):
                    swGr = swGr.iloc[5:-5]
                    Jcorr_fltr = scipy.signal.savgol_filter(
                        swGr.Jcorr.values, Fltr_Win, 3
                    )
                    Jring_fltr = scipy.signal.savgol_filter(
                        swGr.J_ring.values, Fltr_Win, 3
                    )
                    FracH2O2_fltr = scipy.signal.savgol_filter(
                        swGr.Frac_H2O2.values, Fltr_Win, 3
                    )
                    N2_fltr = scipy.signal.savgol_filter(
                        swGr.J_N2_scan.values, Fltr_Win, 3
                    )
                    nORR_fltr = scipy.signal.savgol_filter(
                        swGr.n_ORR.values, Fltr_Win, 3
                    )
                    SwpFltr = pd.DataFrame(swGr[EvRHE])
                    SwpFltr = SwpFltr.assign(
                        **{
                            "Jcorr": Jcorr_fltr,
                            "J_ring": Jring_fltr,
                            "Frac_H2O2": FracH2O2_fltr,
                            "J_N2_scan": N2_fltr,
                            "n_ORR": nORR_fltr,
                        }
                    )
                    SwpFltr.to_excel(
                        ORR_PDD.joinpath(
                            "%s_" % Swp + Path(r.DestFilename).stem
                        ).with_suffix(".xlsx")
                    )
                    if plot_Filter == True:
                        fig, ax1 = plt.subplots()
                        ax2 = ax1.twinx()
                        #                ax.plot(swGr[EvRHE].values,swGr.J_ring.values)
                        #                ax.plot(SwpFltr[EvRHE].values,SwpFltr.J_ring.values)
                        #                        ax.plot(swGr[EvRHE].values,swGr.Frac_H2O2.values)
                        ax2.plot(
                            SwpFltr[EvRHE].values, SwpFltr.Frac_H2O2.values, c="orange"
                        )
                        ax1.plot(SwpFltr[EvRHE].values, SwpFltr.Jcorr.values)
                        ax1.plot(SwpFltr[EvRHE].values, SwpFltr.J_N2_scan.values)
                        ax2.set_ylim(0, 50)
                        ax1.set_ylim(-6.2, 0.1)
                        ax1.set_xlim(0, 1.1)
                        plt.suptitle(Path(r.DestFilename).stem)
                        plt.savefig(
                            ORR_PDD.joinpath(
                                "%s_" % Swp + Path(r.DestFilename).stem
                            ).with_suffix(".png"),
                            dpi=50,
                            bbox_inches="tight",
                        )
                        plt.show()
                        plt.close()


class SD_CNT_Samples:
    def __init__(self):
        pass


#    SDdates = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.SampleID.str.contains('SD')].EXP_date.unique()
#    def StructuralCharacterizations(group_selection = '',type_of_selection = ''):
# a =SampleSelection('Porhp_SiO2','SeriesID',['EA','BET','EC_select2','EC_select_table2'])
#    stgr, x_col, y_col,corr_val,ddir,kwargs = grE,tr.corr_idx[0],tr.corr_idx[1],tr.val,PPDEISbest.joinpath(f'{E[0]}_{E[1]}'), {} for testing

# for a in ww.SeriesOVV.SeriesID.unique()[-6::]:
#    bb = SampleSelection(a,'SeriesID',['EA','BET','EC_select','EC_select_table'])
def CharResTest():
    #    a= SampleSelection(group_selection='CB_paper',type_of_selection='SeriesID',extra_selection_outputs=['EA','BET','RAMAN','MS' ])
    a = SampleSelection(
        group_selection="Porph_SiO2",
        type_of_selection="SeriesID",
        extra_selection_outputs=["RAMAN"],
    )
    a = SampleSelection(
        group_selection="CB_paper",
        type_of_selection="SeriesID",
        extra_selection_outputs=["RAMAN"],
    )
    selEC4 = a.RAMAN_selection_ovv()
    ddir = a.DD_Series.joinpath("corr")
    ddir.mkdir(parents=True, exist_ok=True)
    corr_method, corr_cutoff = "pearson", 0.5  # or default is pearson spearman
    rcorr = selEC4[
        SampleSelection.InterestingCols + SampleSelection.RAMAN_cols_corr
    ].corr(method=corr_method)
    corr_triu = rcorr.where(np.tril(np.ones(rcorr.shape)).astype(np.bool))
    scorr_triu = corr_triu.stack()
    filter_corr = scorr_triu[(scorr_triu.abs() > corr_cutoff) & (scorr_triu.abs() < 1)]
    for n, i in filter_corr.iteritems():
        print(n, i)
        plot_pd_SampleIDs(a.Prep_EA_BET, n[0], n[1], ddir)

    fcorr = filter_corr.unstack()
    plt.subplots(figsize=(20, 15))
    sns.heatmap(fcorr)

    plt.savefig(
        ddir.joinpath(
            "heatmap_{0}_{1}.png".format(corr_method, int(corr_cutoff * 100))
        ),
        dpi=300,
    )
    plt.close()
    plot_pd_SampleIDs(a.Prep_EA_BET, "BET_cat_agg", "AD/AG", ddir)
    ORR = a.EC_retrieve_ORR_pars()
    return a, selEC4, ORR


class SampleSelection:
    """This class initializes and merges the Characterization results for each Sample"""

    SampleLabel_cols = ["SampleLabel", "SampleID", "SeriesID"] + [
        "Colorcode",
        "IndividualLabel",
    ]
    Synthesis_cols = [
        "Overall Weight loss, post-Pyr, post-AL (%)",
        "Yield_precursors(%)",
        "postHT_final_weight_mg",
        "ML",
        "Metal_wt",
        "SizeNP",
        "C_N",
        "MetalPrec",
    ]
    #'''These are used as fillna in prep_ovv : 'Precursor WL overall %', 'Final Mass 1st Batch (mg)'''
    Inter_EA_cols = ["C/N_ratio", "N_content", "C_content", "H_content", "100-CHN"]
    Inter_BET_cols = [
        "BET_SA m2/g_RAW",
        "BET_constant_RAW",
        "BET_Area_RPT",
        "tBoer_Vmicro",
        "tBoer_extSA",
        "t-Method micropore area (m2/g)",
        "Area_in_cell",
    ]
    Sample_Lst_BET = ["BET_CB_meso", "BET_cat_agg", "BET_cat_micro", "BET_cat_meso"]

    Inter_EC_cols = [
        "postAST",
        "Electrolyte",
        "pH",
        "Gas",
        "RPM",
        "Scanrate",
        "EXP_date",
        "Type_Exp",
        "SourceFilename",
        "Loading_cm2",
    ]
    EC_exp_cols = ["postAST", "Electrolyte", "pH", "Gas", "RPM", "Loading_cm2"]
    Inter_RAMAN_cols = [
        "D5_center",
        "D4_center",
        "D3_center",
        "D2_center",
        "D_center",
        "G_center",
        "D_fwhm",
        "D4_fwhm",
        "D3_fwhm",
        "G_fwhm",
        "ID/IG",
        "ID3/IG",
        "ChiSqr",
        "RedChiSqr",
        "ResMean",
        "I_D",
        "I_D3",
        "I_D4",
        "I_D5",
        "I_G",
        "FitModel",
        "I_D2",
        "L_a",
        "L_eq",
        "GD1_sigma",
        "GD1_center",
        "GD1_amplitude",
        "D4D4_sigma",
        "D4D4_center",
        "D4D4_amplitude",
        "D2D2_sigma",
        "D2D2_center",
        "D2D2_amplitude",
        "D1D1_sigma",
        "D1D1_center",
        "D1D1_amplitude",
        "D1D1_fwhm",
        "D2D2_fwhm",
        "D4D4_fwhm",
        "GD1_fwhm",
        "I2D/I2G",
        "ChiSqr.1",
        "RedChiSqr.1",
        "ResMean.1",
        "I_D1D1",
        "I_D2D2",
        "I_D4D4",
        "I_GD1",
    ]
    RAMAN_cols_corr = [
        "D5_center",
        "D4_center",
        "D3_center",
        "D2_center",
        "D_center",
        "G_center",
        "D5_fwhm",
        "D_fwhm",
        "D4_fwhm",
        "D3_fwhm",
        "G_fwhm",
        "I_D",
        "I_D2",
        "I_D3",
        "I_D4",
        "I_D5",
        "I_G",
        "ID/IG",
        "ID3/IG",
        "AD/AG",
        "AD2/AG",
        "AD3/AG",
        "AD4/AG",
        "L_a",
        "L_eq",
        "D4D4_center",
        "D2D2_center",
        "D1D1_center",
        "D1D1_fwhm",
        "D2D2_fwhm",
        "D4D4_fwhm",
        "GD1_fwhm",
        "I2D/I2G",
    ]
    Inter_MS_cols = [
        "D1_pop",
        "D2_pop",
        "D3_pop",
        "D4_pop",
        "sextet1_pop",
        "sextet2_pop",
    ]

    Series_CB_paper = {
        "name": "CB_paper",
        "SerieIDs": ["CB3", "CB4", "CB6", "CB_pyr"],
        "sIDs": [
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
            "DW38A",
            "DW38B",
            "DW38C",
            "DW38D",
            "DW38E",
            "DW38F",
        ],
    }  # removed 'CB7'== JOS11
    Series_Porhp_SiO2 = {
        "name": "Porph_SiO2",
        "SerieIDs": ["Porph_SiO2"],
        "sIDs": ["JOS1", "JOS4", "JOS2", "JOS3", "JOS5"],
    }
    Series_Co_PANI = {
        "name": "Co_PANI",
        "SerieIDs": ["Co_PANI"],
        "sIDs": ["JOS6", "JOS7", "JOS8"],
    }
    Series_ML_SiO2 = {
        "name": "ML_SiO2",
        "SerieIDs": ["ML_SiO2"],
        "sIDs": [
            "JOS1",
            "JOS4",
            "PTA7",
            "PTA200",
            "PTA20",
            "PTA50",
            "PTA100",
            "PTB 0,5",
            "PTB 1",
            "PTB 2",
            "PTB 5",
            "PTB 20",
        ],
    }

    Series = pd.DataFrame(
        [Series_CB_paper, Series_Co_PANI, Series_ML_SiO2, Series_Porhp_SiO2]
    ).set_index("name")

    EC_EIS_par_cols = [
        "Qad",
        "nAd",
        "Cdlp",
        "nDL",
        "Rct",
        "Rct_kin",
        "Rs",
        "Rorr",
        "Qad+Cdlp",
        "RedChisqr2",
        "RedChisqr1",
        "R3",
        "Q3",
        "n3",
    ]
    EC_EIS_models = [
        "Model(Singh2015_RQRQR)",
        "Model(Singh2015_RQRWR)",
        "Model(Singh2015_R3RQ)",
        "Model(Bandarenka_2011_RQRQR)",
    ]
    EC_ORR_kin_par_cols = [
        "E_onset",
        "E_half",
        "J_diff_lim",
        "Jkin_075",
        "Jkin_080",
        "TSa_l",
        "TSb_l",
        "TSa_h",
        "TSb_h",
        "Jring_050",
        "FracH2O2_050",
    ]

    EC_HER_par_cols = [
        "E_AppV_RHE",
        "E_AppV_RHE_upper",
        "TF_a",
        "TF_b",
        "TF_fit_error",
        "TafelSlope",
        "j_lower",
        "j_upper",
    ]

    EC_N2Cdl_cols = [
        "Sweep_Type_N2",
        "ScanRates",
        "j_Acm2_fit",
        "E_AppV_RHE",
        "Cdl",
        "Cdl_R",
        "Cdl_fit",
        "EASA_m2",
        "Cdl_corr",
    ]

    acid1500 = '(pH < 6) & (RPM_DAC > 1100) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
    alk1500 = '(pH > 6) & (RPM_DAC > 1100) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'

    EC_out_cols = SampleLabel_cols + Inter_EC_cols
    Inter_Prep_cols = SampleLabel_cols + Synthesis_cols + Sample_Lst_BET
    Character_xcols = Inter_BET_cols[:-1] + Inter_EA_cols

    InterestingCols = (
        SampleLabel_cols
        + Synthesis_cols
        + Inter_EA_cols
        + Inter_BET_cols
        + Sample_Lst_BET
    )

    InterestingXcorrCols = InterestingCols + RAMAN_cols_corr + EC_EIS_par_cols

    def __init__(
        self, group_selection, type_of_selection, extra_selection_outputs=None
    ):
        if "CB_paper" in group_selection:
            self.group_selection = SampleSelection.Series_CB_paper["slice"]
            self.group_name = SampleSelection.Series_CB_paper["name"]
        else:
            self.group_selection = group_selection
            self.group_name = group_selection
        self.type_of_selection = type_of_selection
        self.extra_selection_outputs = extra_selection_outputs
        self.DD_Series = SampleSelection.DestinationDir_Selection(self)
        self.Prep_EA_BET = SampleSelection.Preparation_ovv(self, write_prep=False)
        #        if group_selection == '*' or 'all' in group_selection:
        #            SampleSelection.Ovv_Selection(self)
        #
        if extra_selection_outputs:
            #            self.Prep_EA_BET = SampleSelection.Preparation_ovv(self,write_prep=False)
            if "write_prep" in extra_selection_outputs:
                self.SeriesOVV = SampleSelection.Preparation_ovv(self, write_prep=True)
            elif extra_selection_outputs[0] in ["", "*", "all"]:
                self.selEC = pd.DataFrame()
            #                self.SeriesOVV = SampleSelection.Preparation_ovv(self,write_prep=False)
            else:
                print(type_of_selection)
                #                if extra_selection_outputs:
                for i in extra_selection_outputs:
                    if i == "EA":
                        self.EA_DF = SampleSelection.EA_selection_ovv(self)
                    elif i == "BET":
                        self.BET_DF = SampleSelection.BET_selection_ovv(self, "plot")
                    elif i == "RAMAN":
                        self.RAMAN_sel = SampleSelection.RAMAN_selection_ovv(self)
                    elif i == "MS":
                        self.MS_sel = SampleSelection.Mossbauer_selection_ovv(self)
                    #                        print(a)
                    elif "EC" in i:
                        self.DD_Series_EC = self.DD_Series.joinpath("EC")
                        self.DD_Series_EC.mkdir(parents=True, exist_ok=True)
                        if "reload" in i:
                            print("Reconstructing EC selection preparation ...")
                            self.selEC = SampleSelection.EC_selection_ovv_prep(
                                self, force_reload=True
                            )
                        else:
                            print("Loading EC selection preparation ...")
                            self.selEC = SampleSelection.EC_selection_ovv_prep(self)
                            if any([a for a in ["pars", "kinetic"] if a in i]):
                                self.orr_kin_pars = self.EC_retrieve_ORR_pars()
                                self.eis_pars = self.EC_retrieve_EIS_pars()

                        if i == "EC_select_table":
                            #                                self.EC_sel = SampleSelection.EC_selection_ovv_prep(self)
                            self.SerieOVV_EC_all = (
                                SampleSelection.Retrieve_EC_Table_from_OVVs(self)
                            )
                            self.SerieOVV_exp_table = (
                                SampleSelection.Retrieve_EC_Table_ExpFiles(self)
                            )

        else:
            pass

    #            self.Prep_EA_BET = SampleSelection.Preparation_ovv(self,write_prep=False)

    def DestinationDir_Selection(self):
        if self.group_selection in ["", "*", "all"]:
            DD_Series = FileHelper.FindExpFolder("VERSASTAT").DestDir.parent.joinpath(
                "Series", "all"
            )
        elif self.group_selection == SampleSelection.Series_CB_paper:
            DD_Series = FileHelper.FindExpFolder("VERSASTAT").DestDir.parent.joinpath(
                "Series", "CB_paper"
            )
        else:
            DD_Series = FileHelper.FindExpFolder("VERSASTAT").DestDir.parent.joinpath(
                "Series", "%s" % self.group_selection
            )
        DD_Series.mkdir(parents=True, exist_ok=True)
        return DD_Series

    #        self.SerieOVV,self.DD_Series = SampleSelection().OvvSelection(group_selection,type_of_selection)
    def Ovv_Selection(self):
        #%% ===== loading of all preparation overview and experimental folders/files ======
        #        group_selection = 'BP2000'
        #        SeriesID = 'CB_bim'
        #        SeriesID = 'Porph_SiO2'
        #        type_of_selection = 'C_N'
        #        group_selection,type_of_selection = 'Porph_SiO2','SeriesID'
        group_selection, type_of_selection = (
            self.group_selection,
            self.type_of_selection,
        )
        #        group_selection,type_of_selection = a.group_selection,a.type_of_selection
        if type(group_selection) == str:
            group_selection = [group_selection]

        try:
            PrepOVV = FileHelper.FindExpFolder().LoadPrepOVV()
            SampleCodeLst = FileHelper.FindExpFolder().LoadSampleCode()
        except:
            print("No preparation overviews!")
        PrepOVV_ext = pd.merge(
            PrepOVV, SampleCodeLst, on=["SampleID", "SeriesID"], how="outer"
        )
        PrepOVV_ext = PrepOVV_ext.assign(
            **{
                "postHT_final_weight_mg": PrepOVV_ext[
                    "EXP post-pyr-post-AL FINAL weight (g)"
                ]
                * 1000
            }
        )
        PrepOVV_ext.SampleLabel = PrepOVV_ext.SampleLabel.fillna(
            value=PrepOVV_ext.SampleCode
        )
        #        print(PrepOVV_ext.SeriesID.unique())
        PrepOVV_ext["postHT_final_weight_mg"] = PrepOVV_ext[
            "postHT_final_weight_mg"
        ].fillna(PrepOVV_ext["Final Mass 1st Batch (mg)"])
        PrepOVV_ext["Overall Weight loss, post-Pyr, post-AL (%)"] = PrepOVV_ext[
            "Overall Weight loss, post-Pyr, post-AL (%)"
        ].fillna(PrepOVV_ext["Precursor WL overall %"])
        PrepOVV_ext["Yield_precursors(%)"] = (
            100 - PrepOVV_ext["Overall Weight loss, post-Pyr, post-AL (%)"]
        )
        PrepOVV_ext["SizeNP"] = PrepOVV_ext["SizeNP"].astype(float)

        if group_selection in ["", "*", "all"] or type_of_selection in ["", "*", "all"]:
            print(
                "No selection of Series/Group or Columns, output all from Preparations in one DF"
            )
            SerieOVV = PrepOVV_ext

        elif type_of_selection in PrepOVV_ext.columns:
            SerieOVV = PrepOVV_ext.loc[
                PrepOVV_ext[type_of_selection].str.contains("|".join(group_selection))
                == True
            ]
            print(
                'Selecting on column {0} containing "{3}"\nlen before: {1} and after {2}'.format(
                    type_of_selection, len(PrepOVV_ext), len(SerieOVV), *group_selection
                )
            )
        elif type_of_selection not in PrepOVV_ext.columns:
            print(
                "Please choose another (not %s) column: %s"
                % (type_of_selection, PrepOVV_ext.columns)
            )
            SerieOVV = pd.DataFrame()
        else:
            print(
                "No selection of Series/Group or Columns, output all from Preparations in one DF"
            )
            SerieOVV = pd.DataFrame()
        #            SerieOVV = PrepOVV_ext
        #        prep_ovv_selection = SerieOVV
        try:
            SerieOVV = SerieOVV.sort_values(by=["SeriesID", "SampleID"])
        except:
            pass
        return SerieOVV[SampleSelection.Inter_Prep_cols]

    #        elif group_selection != '':
    #            SerieOVV = PrepOVV_ext.loc[PrepOVV_ext[type_of_selection].str.contains('{0}'.format(group_selection)) == True]
    #        else:
    #            print('No selection of Series/Group or Columns, output all from Preparations in one DF')
    #            SerieOVV = PrepOVV_ext
    # IMPORTANT Adaptations to PostEC overview.... implement later to Constructor PostEC method!!
    #%%
    def Preparation_ovv(self, write_prep=True):
        #        [i for i in FileHelper.FindExpFolder('EA').DestDir.parent.glob('*') if i.is_dir()]
        """Read Characeterization experimental results from Organized Data folder and Merge with the
        Preparation overview on Samples in the SeriesID"""

        group_selection = self.group_selection
        group_name = self.group_name
        prep_ovv_selection = SampleSelection.Ovv_Selection(self)
        extra_selection_outputs = self.extra_selection_outputs
        DD_Series = SampleSelection.DestinationDir_Selection(self)

        EA_results_imp = pd.read_excel(
            list(FileHelper.FindExpFolder("EA").DestDir.rglob("*xlsx"))[0]
        )
        EA_results = EA_results_imp[["SampleID"] + SampleSelection.Inter_EA_cols]
        #        'EA_mass', 'C/N_ratio', 'N_content', 'C_content','H_content', '100-CHN']]

        BET_results_lst = [
            i
            for i in list(FileHelper.FindExpFolder("BET").DestDir.rglob("BET_OVV*xlsx"))
            if not "_Conflict" in i.stem
        ]
        BET_results = pd.concat(
            [pd.read_excel(i, index_col=0) for i in BET_results_lst], sort=False
        ).drop_duplicates()
        BET_SampleIDs = [
            FileHelper.FindSampleID.try_find_sampleID((Path(str(i)).stem))[0]
            for i in BET_results.BET_filename
        ]
        BET_results = BET_results.assign(**{"SampleID": BET_SampleIDs})

        RAMAN_results_files = [
            i
            for i in list(
                FileHelper.FindExpFolder("RAMAN").DestDir.rglob("AllPars*xlsx")
            )
            if not "_Conflict" in i.stem
        ]
        RAMAN_results = pd.concat(
            [pd.read_excel(i, index_col=0) for i in RAMAN_results_files], sort=False
        )
        #        pd.read_excel(list(FileHelper.FindExpFolder('RAMAN').DestDir.rglob('RAMAN_OVV*'))[0])
        RAMAN_results["SampleID"] = [
            i[0] for i in RAMAN_results.SampleID.str.split("_mean")
        ]
        #        '*Wallace*/*overview*xlsx
        Moss_data_file = list(
            FileHelper.FindExpFolder("Mossbauer").DestDir.rglob("2019-05*OVV*.xlsx")
        )[0]
        MossB_results = pd.read_excel(Moss_data_file).dropna(axis=1, how="all")
        print("Mossbauer file used: {0}".format(Moss_data_file))
        #        MossB_results['MS_par_file'] = Moss_data_file

        Prep_EA = pd.merge(
            prep_ovv_selection,
            EA_results,
            on=["SampleID"],
            how="left",
            suffixes=["", "_EA"],
        )
        Prep_EA_BET = pd.merge(
            Prep_EA,
            BET_results,
            left_on="SampleID",
            right_on="Sample type",
            how="left",
            suffixes=["", "_BET"],
        )
        Prep_EA_BET_RAMAN = pd.merge(
            Prep_EA_BET,
            RAMAN_results,
            on=["SampleID"],
            how="left",
            suffixes=["", "_RAMAN"],
        )
        Prep_EA_BET_RAMAN_MS = pd.merge(
            Prep_EA_BET_RAMAN,
            MossB_results,
            on=["SampleID"],
            how="left",
            suffixes=["", "_MS"],
        )
        Prep_EA_BET_RAMAN_MS = Prep_EA_BET_RAMAN_MS.dropna(subset=["SampleID"])
        MS_Fewt_add = pd.concat(
            [
                pd.DataFrame(i)
                for i in [
                    {
                        i
                        + "_Fe_wt": 0.01
                        * Prep_EA_BET_RAMAN_MS[i].fillna(0).astype(float)
                        * Prep_EA_BET_RAMAN_MS.Metal_wt
                    }
                    for i in SampleSelection.Inter_MS_cols
                ]
            ],
            axis=1,
        )
        Prep_EA_BET_RAMAN_MS = pd.concat([Prep_EA_BET_RAMAN_MS, MS_Fewt_add], axis=1)
        #        pd.merge(SerieOVV,postOVVout,on='SampleID',how='left',suffixes=['','_B'])
        #        Prep_EA_BET[InterestingCols]
        #        selEA = Prep_EA_BET.dropna(subset=['N_content']), selEA.plot(x='SampleLabel',y='H_content',kind='barh')
        #        selEA.plot(x='SampleID',y='H_content',kind='barh')
        #        Prep_EA_BET.plot(x='SampleLabel',y='Overall Weight loss, post-Pyr, post-AL (%)',kind='barh')
        #        Prep_EA_BET[Inter_Prep_cols]
        serOVV_out_lst = []
        prep_out_path = DD_Series.joinpath("Prep_%s.xlsx" % group_name)
        #        prep_out_path_prec = self.DD_Series.joinpath('Prep_%s_Prec.xlsx' %group_selection)
        #        Prep_EA_BET = a.Prep_EA_BET
        #        Prep_EA_BET_RAMAN.loc[Prep_EA_BET_RAMAN.SeriesID.str.contains('|'.join(group_selection)),:]
        #        prep_out_path = a.DD_Series.joinpath('Prep_%s.xlsx' %group_selection)
        #        prep_out_path.parent.mkdir(parents=True,exist_ok=True)
        #        print(write_prep,extra_selection_outputs[0] in ['','*','all'] )
        if extra_selection_outputs:
            if (
                write_prep == True
                or prep_out_path.is_file() == False
                and extra_selection_outputs[0] not in ["", "*", "all"]
            ):

                for Sn, Sgr in Prep_EA_BET_RAMAN_MS.groupby(by="SeriesID"):
                    if Sgr.empty:
                        print('Group "%s" is empty , no output' % Sn)
                        continue
                    #                    with pd.ExcelWriter(prep_out_path) as writer:
                    slbl = [
                        i for i in Sgr.SampleLabel.unique() if len(i.split("_")) > 1
                    ][0]
                    PrecLabel = "_".join(str(slbl).split("_")[:-1])
                    serOVV_out_lst.append(
                        (Sn, PrecLabel, ", ".join(Sgr.SampleID.values))
                    )
                    print(Sn)
                    #                        Sgr[SampleSelection.Inter_Prep_cols+SampleSelection.Inter_EA_cols[1::]+['BET_Area_RPT','t-Method micropore area (m2/g)']].set_index('SampleID').to_excel(writer,
                    #                           sheet_name='%s'%Sn, float_format='%.2f' )
                    Sgr[
                        SampleSelection.Inter_Prep_cols
                        + SampleSelection.Inter_EA_cols[1::]
                        + ["BET_Area_RPT", "t-Method micropore area (m2/g)"]
                    ].set_index("SampleID").to_excel(
                        self.DD_Series.joinpath("Prep_%s.xlsx" % group_selection),
                        float_format="%.2f",
                    )
                    PrecSerie = pd.DataFrame(
                        serOVV_out_lst,
                        columns=["PrecursorID", "PrecursorLabel", "Samples_DW"],
                    )
                    PrecSerie.to_excel(
                        DD_Series.joinpath("Prep_%s_Prec.xlsx" % group_name)
                    )
                print("Preparation saved to:%s" % prep_out_path)
            elif write_prep == False and extra_selection_outputs[0] in ["", "*", "all"]:

                Prep_EA_BET_RAMAN_MS.sort_values(by="SampleID", inplace=True)
                Prep_EA_BET_RAMAN_MS.to_excel(DD_Series.joinpath("Prep_ovv_all.xlsx"))
        else:
            print("No extra selection outputs")
        return Prep_EA_BET_RAMAN_MS

    #        Prep_EA_BETPrep_EA_BET['EXP post-pyr-post-AL FINAL weight (g)'
    #         EA_cols = ['C/N_ratio', 'N_content', 'C_content','H_content', '100-CHN']
    #        for EAcol in Inter_EA_cols:aa
    #            PlotDF = DataPerSample.FilternaPlot(Prep_EA_BET,'SampleLabel','%s'%EAcol,'barh',DD_Series)
    #    Extra (stupid) color stuff
    #        norm = matplotlib.colors.Normalize(vmin=0,vmax=int(SerieOVV.Colorcode.max()))
    #        clst = cm.get_cmap(plt.get_cmap('Set1'))(range(0,int(SerieOVV.Colorcode.max())))
    #        ====

    def EA_selection_ovv(self):
        #%% ===== Plot and Output Elementary Analysis ======
        Prep_EA_BET = self.Prep_EA_BET
        DD_Series = self.DD_Series
        group_selection = self.group_selection
        group_name = self.group_name
        EA_out_cols = set(
            SampleSelection.SampleLabel_cols
            + SampleSelection.Inter_EA_cols
            + SampleSelection.InterestingCols
        )
        EA_DF = Prep_EA_BET[EA_out_cols]
        EA_DF = (
            EA_DF.dropna(subset=SampleSelection.Inter_EA_cols)
            .sort_values(by="SampleID")
            .drop_duplicates()
        )
        if not EA_DF.empty:
            Lbls = list(
                zip(
                    EA_DF.SampleID.values,
                    EA_DF.SampleLabel.values,
                    EA_DF.SeriesID.values,
                )
            )
            if len(EA_DF.SeriesID.unique()) > 1:
                sample_labels = ["%s, %s : %s" % (i[0], i[2], i[1]) for i in Lbls]
            else:
                sample_labels = ["%s : %s" % (i[0], i[1]) for i in Lbls]

            fig, ax = plt.subplots(figsize=(12, 10))
            #        sample_labels = ['%s : %s'%(i[0],i[1]) for i in  list(zip(EA_DF.SampleID.values,EA_DF.SampleLabel.values, EA_DF.SeriesID.values))]
            y_pos = np.arange(len(sample_labels))
            Nv, Cv, Hv, Xv = (
                EA_DF["N_content"].values,
                EA_DF["C_content"].values,
                EA_DF["H_content"].values,
                EA_DF["100-CHN"].values,
            )
            #        PlotDF.set_index('SampleLabel')[['N_content', 'C_content','H_content', '100-CHN']].plot.barh(stacked=True,legend=True)
            Cc = plt.barh(y_pos, Cv, label="[C]", color="dimgrey")
            Nc = plt.barh(y_pos, Nv, left=Cv, label="[N]", color="limegreen")
            Xc = plt.barh(
                y_pos, Xv, left=(Cv + Nv), label="[O,Fe,..]", color="firebrick"
            )
            Hc = plt.barh(
                y_pos, Hv, left=(Cv + Nv + Xv), label="[H]", color="royalblue"
            )
            ax.legend(
                bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
                loc="lower left",
                ncol=2,
                mode="expand",
                borderaxespad=0.0,
            )
            ax.set_yticks(y_pos)
            ax.set_yticklabels(sample_labels)
            ax.set_xlabel("wt%")
            #        ax.set_title('Elementary Analysis C,H,N',bbox_to_anchor=(0., 1.02, 1., .102))
            #        ax.text(0,1.2,'Elementary Analysis C,H,N')
            fig.suptitle(
                "Elementary Analysis C, H, N for serie %s" % group_name, y=1.15
            )
            #        , 'C_content','H_content', '100-CHN']].plot.barh(stacked=True,legend=True)
            EA_path = DD_Series.joinpath("EA_series_%s" % group_name)
            plt.savefig(EA_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
            EA_DF[EA_out_cols].set_index("SampleID").to_excel(
                EA_path.with_suffix(".xlsx")
            )
        else:
            print("EA selection empty")
            EA_DF = pd.DataFrame()
        EA_sel = EA_DF
        return EA_sel

    #        ====== ====
    #%% ===== Plot and Output BET ======
    def BET_selection_ovv(self, *args):
        #        BET_ylbl = Inter_BET_cols[2]
        Prep_EA_BET = self.Prep_EA_BET
        DD_Series = self.DD_Series
        group_selection = self.group_selection
        group_name = self.group_name
        BET_out_cols = (
            SampleSelection.Inter_Prep_cols
            + SampleSelection.Inter_BET_cols
            + SampleSelection.Sample_Lst_BET
        )
        BET_path = DD_Series.parent.joinpath(
            "%s_series_%s" % (DD_Series.name, group_name)
        )
        #        + ['Colorcode']
        #        BET_ylbl = 'BET_Area_RPT'
        #        Prep_EA_BET = a.Prep_EA_BET
        #        DD_Series = a.DD_Series
        #        group_selection = a.group_selection
        #        group_name = a.group_name
        #
        selBET = Prep_EA_BET[BET_out_cols].dropna(subset=SampleSelection.Inter_BET_cols)
        selBET = selBET.dropna(subset=BET_out_cols, how="all").sort_values(
            by="SampleID"
        )
        selBET = selBET.drop_duplicates(subset=BET_out_cols)

        selBET[BET_out_cols].set_index("SampleID").to_excel(
            BET_path.with_suffix(".xlsx")
        )

        if not selBET.empty and "plot" in args:

            Lbls = list(
                zip(
                    selBET.SampleID.values,
                    selBET.SampleLabel.values,
                    selBET.SeriesID.values,
                )
            )
            OriginColors = Characterization_TypeSetting.OriginColorList()
            if len(selBET.SeriesID.unique()) > 1:
                sample_labels = ["%s, %s : %s" % (i[0], i[2], i[1]) for i in Lbls]
            else:
                sample_labels = ["%s : %s" % (i[0], i[1]) for i in Lbls]
            selBET.Colorcode = selBET.Colorcode.fillna(1)
            for BET_ylbl in SampleSelection.Inter_BET_cols:
                if selBET[BET_ylbl].empty:
                    continue
                y_pos = np.arange(len(sample_labels) * 2)[0::2]
                #                np.arange(len(sample_labels))
                x_vals = selBET[BET_ylbl].values
                #                if 'BET_RPT' in BET_ylbl:
                #                    selBET[]
                #        Inter_BET_cols
                fig, ax = plt.subplots(figsize=(12, 8))
                prop_cycle = plt.rcParams["axes.prop_cycle"]
                colors = prop_cycle.by_key()["color"]
                # (UN)NECESSARILY COMPLICATED ORIGIN COLOR CODE REFERENCE !!!!
                OrCls = [
                    [float(b) / 255 for b in i.split(",")]
                    for i in [
                        (
                            OriginColors.loc[
                                OriginColors.OriginCode == int(i), "RGB"
                            ].values[0]
                        )
                        for i in selBET.Colorcode.values
                    ]
                ]
                #        cl = [colors[int(i)] for i in selBET.Colorcode.values]
                cmap = plt.cm.Set2
                clrs = [cmap(i) for i in selBET.Colorcode.values]
                BET = plt.barh(y_pos, x_vals, 1.5, color=OrCls)
                #        selBET.plot(x='SampleLabel',y='BET_Area_RPT',kind='barh',ax=ax,legend=False,title='BET')
                #        selBET.plot(x='SampleLabel',y='BET_SA m2/g_RAW',kind='barh',ax=ax,legend=False,title='BET')
                #        selBET.plot(x='SampleLabel',y=BET_ylbl, kind='barh',ax=ax,legend=False,title='%s'%BET_ylbl,color=OrCls)
                for i, v in enumerate(selBET[BET_ylbl].values):
                    ax.text(
                        v + v * 0.050,
                        y_pos[i] - 0.1,
                        "%.0f " % v,
                        color="k",
                        fontweight="bold",
                    )
                #        selBET.plot(x='SampleLabel',y='BET_SA m2/g_RAW',kind='barh',ax=ax,legend=False,title='BET')
                #        ax.legend(bbox_to_anchor=(0., 1.2, 1., .102), loc='lower left',
                #           ncol=2, mode="expand", borderaxespad=0.)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(sample_labels, fontsize=14)
                fig.suptitle("%s" % BET_ylbl)
                if selBET[BET_ylbl].dropna().min() < 0:
                    xmin = selBET[BET_ylbl].dropna().min() * 1.2
                else:
                    xmin = 0
                ax.set_xlim([xmin, selBET[BET_ylbl].dropna().max() * 1.2])
                ax.set_xlabel("BET surface area / $m^{2}g^{-1}$")
                BET_path_ylbl = DD_Series.joinpath(
                    "%s_series_%s"
                    % (re.compile("[\W_]+", re.UNICODE).sub("", BET_ylbl), group_name)
                )
                plt.tight_layout()
                try:
                    plt.savefig(
                        BET_path_ylbl.with_suffix(".png"), dpi=100, bbox_inches="tight"
                    )
                except:
                    pass
                plt.show()
                plt.close()
        elif selBET.empty:
            print("BET selection empty")
            selBET = pd.DataFrame()
        else:
            pass
        return selBET

    #%% ===== Plot and Output RAMAN ======
    def RAMAN_selection_ovv(self, *args):
        #        BET_ylbl = Inter_BET_cols[2]
        Prep_EA_BET = self.Prep_EA_BET
        Type_Exp = "RAMAN"
        DD_Series = self.DD_Series.joinpath(Type_Exp)
        DD_Series.mkdir(parents=True, exist_ok=True)
        group_selection = self.group_selection
        group_name = self.group_name
        RAMAN_cols = SampleSelection.Inter_RAMAN_cols
        RAMAN_path = DD_Series.parent.joinpath(
            "%s_series_%s" % (DD_Series.name, group_name)
        )
        #        Prep_EA_BET = a.Prep_EA_BET
        #        DD_Series = a.DD_Series.joinpath('RAMAN')
        #        DD_Series.mkdir(parents=True,exist_ok=True)
        #        group_selection = a.group_selection
        #        group_name = a.group_name
        #        RAMAN_cols = a.Inter_RAMAN_cols
        try:
            selRAMAN = Prep_EA_BET.dropna(subset=RAMAN_cols, how="all").sort_values(
                by="SampleID"
            )
            selRAMAN = Prep_EA_BET.drop_duplicates(subset=RAMAN_cols)
        except:
            print("Raman selection problems")
            selRAMAN = pd.DateFrame()
        if "raw" in args:
            DDir = FileHelper.FindExpFolder(Type_Exp).DestDir
            DD_Series_raw = DD_Series.joinpath("RAW_data")
            DDir_list = [
                list(DDir.rglob("*/{0}_*xlsx".format(i)))
                for i in selRAMAN.SampleID.values[:]
            ]
            DD_Series_raw.mkdir(parents=True, exist_ok=True)
            for a in DDir_list:
                for src in a:
                    dst = DD_Series_raw.joinpath(src.name)
                    shutil.copy2(src, dst)
            print("Raman Raw files copied to: {0}".format(DD_Series_raw))

        if not selRAMAN.empty and "plot" in args:
            selRAMAN[RAMAN_cols + SampleSelection.InterestingCols].set_index(
                "SampleID"
            ).to_excel(RAMAN_path.with_suffix(".xlsx"))
            Lbls = list(
                zip(
                    selRAMAN.SampleID.values,
                    selRAMAN.SampleLabel.values,
                    selRAMAN.SeriesID.values,
                )
            )
            OriginColors = Characterization_TypeSetting.OriginColorList()
            if len(selRAMAN.SeriesID.unique()) > 1:
                sample_labels = ["%s, %s : %s" % (i[0], i[2], i[1]) for i in Lbls]
            else:
                sample_labels = ["%s : %s" % (i[0], i[1]) for i in Lbls]
            selRAMAN.Colorcode = selRAMAN.Colorcode.fillna(1)
            for RAMAN_ylbl in RAMAN_cols:
                if selRAMAN[RAMAN_ylbl].empty:
                    continue
                #                y_pos = np.arange(len(sample_labels))
                y_pos = np.arange(len(sample_labels) * 5)[0::5]
                try:
                    x_vals = [float(i) for i in selRAMAN[RAMAN_ylbl].values]
                except Exception as e:
                    print(
                        "No Raman plot of column: {0} because \, {1}".format(
                            RAMAN_ylbl, e
                        )
                    )
                    continue

                #                if 'BET_RPT' in BET_ylbl:
                #                    selBET[]
                #        Inter_BET_cols
                fig, ax = plt.subplots(figsize=(12, 10))
                prop_cycle = plt.rcParams["axes.prop_cycle"]
                colors = prop_cycle.by_key()["color"]
                # (UN)NECESSARILY COMPLICATED ORIGIN COLOR CODE REFERENCE !!!!
                OrCls = [
                    [float(b) / 255 for b in i.split(",")]
                    for i in [
                        (
                            OriginColors.loc[
                                OriginColors.OriginCode == int(i), "RGB"
                            ].values[0]
                        )
                        for i in selRAMAN.Colorcode.values
                    ]
                ]
                #        cl = [colors[int(i)] for i in selBET.Colorcode.values]
                cmap = plt.cm.Set2
                clrs = [cmap(i) for i in selRAMAN.Colorcode.values]
                BET = plt.barh(y_pos, x_vals, 4, color=OrCls)
                #        selBET.plot(x='SampleLabel',y='BET_Area_RPT',kind='barh',ax=ax,legend=False,title='BET')
                #        selBET.plot(x='SampleLabel',y='BET_SA m2/g_RAW',kind='barh',ax=ax,legend=False,title='BET')
                #        selBET.plot(x='SampleLabel',y=BET_ylbl, kind='barh',ax=ax,legend=False,title='%s'%BET_ylbl,color=OrCls)
                for i, v in enumerate(selRAMAN[RAMAN_ylbl].values):
                    try:
                        if v <= 10:
                            ylbl_text = "%.2f " % v
                        elif v > 10 and v <= 100:
                            ylbl_text = "%.1f " % v
                        elif v > 100:
                            ylbl_text = "%.0f " % v
                        else:
                            ylbl_text = "%.2f " % v
                        ax.text(
                            v + v * 0.050,
                            y_pos[i] - 0.1,
                            ylbl_text,
                            color="k",
                            fontweight="bold",
                        )
                    except Exception as e:
                        print("No text plot: {0}".format(e))
                #        selBET.plot(x='SampleLabel',y='BET_SA m2/g_RAW',kind='barh',ax=ax,legend=False,title='BET')
                #        ax.legend(bbox_to_anchor=(0., 1.2, 1., .102), loc='lower left',
                #           ncol=2, mode="expand", borderaxespad=0.)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(sample_labels)
                fig.suptitle("%s" % RAMAN_ylbl)
                if selRAMAN[RAMAN_ylbl].dropna().min() < 0:
                    xmin = selRAMAN[RAMAN_ylbl].dropna().min() * 1.2
                else:
                    xmin = 0
                ax.set_xlim([xmin, selRAMAN[RAMAN_ylbl].dropna().max() * 1.2])
                ax.set_xlabel("{0}".format(RAMAN_ylbl))
                RAMAN_path_ylbl = DD_Series.joinpath(
                    "%s_series_%s"
                    % (re.compile("[\W_]+", re.UNICODE).sub("", RAMAN_ylbl), group_name)
                )
                plt.tight_layout()
                try:
                    plt.savefig(
                        RAMAN_path_ylbl.with_suffix(".png"),
                        dpi=100,
                        bbox_inches="tight",
                    )
                except:
                    pass
                plt.show()
                plt.close()
        else:
            print("RAMAN selection empty")
        #            selRAMAN = pd.DataFrame()
        return selRAMAN

    #%% ===== Plot and Output Mößbauer ======
    def Mossbauer_selection_ovv(self, *args):
        #        BET_ylbl = Inter_BET_cols[2]
        Prep_EA_BET = self.Prep_EA_BET
        DD_Series = self.DD_Series.joinpath("Mossbauer")
        DD_Series.mkdir(parents=True, exist_ok=True)
        group_selection = self.group_selection
        group_name = self.group_name
        #        MS_cols = SampleSelection.Inter_MS_cols
        #        Prep_EA_BET = a.Prep_EA_BET
        #        DD_Series = a.DD_Series.joinpath('RAMAN')
        #        DD_Series.mkdir(parents=True,exist_ok=True)
        #        group_selection = a.group_selection
        #        RAMAN_cols = a.Inter_RAMAN_cols
        MS_pop_cols = [
            i
            for i in Prep_EA_BET.columns
            if i.split("_")[-1] in ["pop", "CS", "d", "H"]
        ]
        MS_colors = {
            "D1": "green",
            "D2": "blue",
            "D3": "violet",
            "sextet1": "yellow",
            "sextet2": "yellow",
            "singlet": "orange",
        }
        #        [i for i in a.Prep_EA_BET.columns if 'pop' in i.split('_')[-1]]
        selMS = (
            Prep_EA_BET.dropna(subset=MS_pop_cols, how="all")
            .drop_duplicates(subset=MS_pop_cols)
            .sort_values(by="SampleID")
        )
        selMS["ML"] = selMS["ML"].astype(float)
        for pop_col in [i for i in MS_pop_cols if "pop" in i]:
            selMS[pop_col].fillna(0, inplace=True)
            pop_Fe_wt_col = "{0}_Fe_wt".format(pop_col)
            selMS = selMS.assign(
                **{pop_Fe_wt_col: 1e-02 * selMS[pop_col] * selMS["Metal_wt"]}
            )
            MS_pop_cols.append(pop_Fe_wt_col)
        selmS_output_cols_all = (
            SampleSelection.SampleLabel_cols
            + SampleSelection.Synthesis_cols
            + SampleSelection.Inter_EA_cols
            + MS_pop_cols
            + SampleSelection.Sample_Lst_BET
        )
        selmS_output_cols = [i for i in selmS_output_cols_all if i in selMS.columns]
        EA_MS = selMS[selmS_output_cols].drop_duplicates()
        EA_MS.to_excel(DD_Series.parent.joinpath("MS_EA_selection.xlsx"))
        if not selMS.empty and "plot" in args:

            for xscatter in SampleSelection.InterestingCols:
                #            ['ML','N_content']
                figp, axp = plt.subplots()
                for pop in [i for i in MS_pop_cols if "pop" in i and "_Fe_wt" in i]:
                    if not selMS[xscatter].dropna(how="all").empty:
                        try:

                            selMS.loc[selMS[pop] > 0].plot(
                                x=xscatter,
                                y=pop,
                                kind="scatter",
                                color=MS_colors[pop.split("_")[0]],
                                ax=axp,
                                logx=False,
                                logy=0,
                                s=80,
                                label=pop,
                            )
                        except:
                            plt.close()
                            continue
                            pass
                    else:
                        plt.close()

                        continue
                axp.legend(
                    bbox_to_anchor=(0.0, 1.2, 1.0, 0.102),
                    loc="lower left",
                    ncol=2,
                    mode="expand",
                    borderaxespad=0.0,
                )
                #                axp.legend(ncol=2)
                plt.show()
                MS_path_sct = DD_Series.joinpath(
                    "%s_scatter_%s"
                    % (re.compile("[\W_]+", re.UNICODE).sub("", xscatter), pop)
                )
                figp.savefig(MS_path_sct.with_suffix(".png"))
                print("MS Scatter fig saved to: {0}".format(MS_path_sct))
                plt.close()

            #                selMS.plot(x='ML',y=pop,kind='scatter',ax=axp)
            #            axp.legend(False)
            Lbls = list(
                zip(
                    selMS.SampleID.values,
                    selMS.SampleLabel.values,
                    selMS.SeriesID.values,
                )
            )
            OriginColors = Characterization_TypeSetting.OriginColorList()
            if len(selMS.SeriesID.unique()) > 1:
                sample_labels = ["%s, %s : %s" % (i[0], i[2], i[1]) for i in Lbls]
            else:
                sample_labels = ["%s : %s" % (i[0], i[1]) for i in Lbls]
            selMS.Colorcode = selMS.Colorcode.fillna(1)
            for MS_ylbl in MS_pop_cols:
                if (
                    selMS[MS_ylbl].empty
                    or selMS[MS_ylbl].dropna(how="all").empty
                    or selMS[MS_ylbl].astype(float).sum() == 0
                ):
                    continue
                y_pos = np.arange(len(sample_labels))
                try:
                    x_vals = [float(i) for i in selMS[MS_ylbl].values]
                except Exception as e:
                    print("No MS plot of column: {0} because \, {1}".format(MS_ylbl, e))
                    continue
                #                if 'BET_RPT' in BET_ylbl:
                #                    selBET[]
                #        Inter_BET_cols
                fig, ax = plt.subplots()
                prop_cycle = plt.rcParams["axes.prop_cycle"]
                colors = prop_cycle.by_key()["color"]
                # (UN)NECESSARILY COMPLICATED ORIGIN COLOR CODE REFERENCE !!!!
                OrCls = [
                    [float(b) / 255 for b in i.split(",")]
                    for i in [
                        (
                            OriginColors.loc[
                                OriginColors.OriginCode == int(i), "RGB"
                            ].values[0]
                        )
                        for i in selMS.Colorcode.values
                    ]
                ]
                #        cl = [colors[int(i)] for i in selBET.Colorcode.values]
                cmap = plt.cm.Set2
                clrs = [cmap(i) for i in selMS.Colorcode.values]
                BET = plt.barh(y_pos, x_vals, color=OrCls)
                #        selBET.plot(x='SampleLabel',y='BET_Area_RPT',kind='barh',ax=ax,legend=False,title='BET')
                #        selBET.plot(x='SampleLabel',y='BET_SA m2/g_RAW',kind='barh',ax=ax,legend=False,title='BET')
                #        selBET.plot(x='SampleLabel',y=BET_ylbl, kind='barh',ax=ax,legend=False,title='%s'%BET_ylbl,color=OrCls)
                for i, v in enumerate(selMS[MS_ylbl].values):
                    try:
                        if v <= 10:
                            ylbl_text = "%.2f " % v
                        elif v > 10 and v <= 100:
                            ylbl_text = "%.1f " % v
                        elif v > 100:
                            ylbl_text = "%.0f " % v
                        else:
                            ylbl_text = "%.2f " % v
                        ax.text(
                            v + v * 0.050,
                            i - 0.1,
                            ylbl_text,
                            color="k",
                            fontweight="bold",
                        )
                    except Exception as e:
                        print("No text plot: {0}".format(e))
                #        selBET.plot(x='SampleLabel',y='BET_SA m2/g_RAW',kind='barh',ax=ax,legend=False,title='BET')
                #        ax.legend(bbox_to_anchor=(0., 1.2, 1., .102), loc='lower left',
                #           ncol=2, mode="expand", borderaxespad=0.)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(sample_labels)
                fig.suptitle("%s" % MS_ylbl)
                if selMS[MS_ylbl].dropna().min() < 0:
                    xmin = selMS[MS_ylbl].dropna().min() * 1.2
                else:
                    xmin = 0
                ax.set_xlim([xmin, selMS[MS_ylbl].dropna().max() * 1.2])
                ax.set_xlabel("{0}".format(MS_ylbl))
                MS_path = DD_Series.joinpath(
                    "%s_series_%s"
                    % (re.compile("[\W_]+", re.UNICODE).sub("", MS_ylbl), group_name)
                )

                try:
                    plt.savefig(
                        MS_path.with_suffix(".png"), dpi=300, bbox_inches="tight"
                    )
                except:
                    pass
                selMS[MS_pop_cols + SampleSelection.InterestingCols].set_index(
                    "SampleID"
                ).to_excel(MS_path.with_suffix(".xlsx"), float_format="%.0f")
                print("MS excel and fig saved to: {0}".format(MS_path))
                plt.show()
                plt.close()
        else:
            print("MS selection empty or not plot")
        #            selMS = pd.DataFrame()
        return selMS

    #        ==== ====
    #%% ===== Plot and Output EC selection OVV prep ======
    def EC_selection_ovv_prep(self, force_reload=False, force_postOVV=False):
        Prep_EA_BET = self.Prep_EA_BET
        DD_Series_EC = self.DD_Series_EC
        file_EC = DD_Series_EC.joinpath("EC_selection_OVV.xlsx")

        if file_EC.is_file() == False or force_reload == True:
            postOVVout_Load = PostEC.LoadPostOVV(force_postOVV)
            postOVVout = postOVVout_Load.loc[
                postOVVout_Load.SampleID.str.contains(
                    "|".join(Prep_EA_BET.SampleID.values)
                )
            ]
            postOVVout.loc[
                (postOVVout.SourceFilename.str.contains("_300_first") == True)
                & (postOVVout.Type_Exp == "N2_CV"),
                "Scanrate",
            ] = 299
            postOVVout.SourceFilename = postOVVout.SourceFilename.fillna(
                "no_name"
            ).astype(str)
            postOVVout = postOVVout.loc[
                (postOVVout.SourceFilename.str.contains("_Conflict") == False) == True
            ]
            #        postOVVout.loc[postOVVout.SourceFilename.str.contains('_Conflict') == False]
            postOVVout["SourceStem"] = [
                "_".join(Path(i).stem.split("300_first")[0].split("_")[0:-1])
                for i in postOVVout.SourceFilename.values
            ]
            postOVVout = postOVVout[
                [
                    "DestFilename",
                    "SampleID",
                    "Status",
                    "Electrolyte",
                    "Gas",
                    "RPM",
                    "Scanrate",
                    "EXP_date",
                    "Type_Exp",
                    "SourceFilename",
                    "SourceStem",
                    "Exp_dir",
                ]
            ]
            postOVVout["expDT"] = pd.to_datetime(postOVVout["EXP_date"])
            Bol_file = []
            BoL_found_match = []
            for n, r in postOVVout.iterrows():
                if r.Status == "EoL":
                    EoL_date = pd.to_datetime(r.expDT)
                    p_date_slice = postOVVout.loc[
                        (postOVVout.expDT <= EoL_date)
                        & (postOVVout.expDT >= EoL_date - pd.to_timedelta(1, unit="D")),
                        :,
                    ]
                    BoL_match = p_date_slice.query(
                        '(SampleID == @r.SampleID) & (Electrolyte == @r.Electrolyte) & (Status == "BoL") & (Type_Exp == @r.Type_Exp)'
                    )
                    if r.Type_Exp == "N2_CV":
                        #                    BoL_match = BoL_match.loc[BoL_match.Scanrate == r.Scanrate,:]
                        BoL_match = BoL_match.query("Scanrate == @r.Scanrate")
                    #                print(p_date_slice)
                    if not BoL_match.empty:
                        BoL_date_match = [
                            i for i in BoL_match.EXP_date.unique() if i != r.EXP_date
                        ]
                        #                    print('row:' ,Path(r.DestFilename).name,'\n match :',BoL_date_match)
                        if BoL_date_match:
                            BoL_found_match.append(BoL_date_match[0])
                        else:
                            BoL_date_match_inkl_EoL = [
                                i for i in BoL_match.EXP_date.unique()
                            ]
                            if BoL_date_match_inkl_EoL:
                                BoL_found_match.append(BoL_date_match_inkl_EoL[0])
                            else:
                                BoL_found_match.append("None")
                    elif BoL_match.empty:
                        BoL_found_match.append("None")
                    else:
                        BoL_found_match.append("None")
                #                BoL_file = 1
                elif r.Status == "BoL":
                    BoL_found_match.append(r.EXP_date)
                else:
                    BoL_found_match.append(r.EXP_date)
            print(len(postOVVout), len(BoL_found_match))
            postOVVout["BoL_date"] = BoL_found_match
            pH_lst = [
                FileHelper.FindSampleID.determine_pH_from_filename(i)["pH"][0]
                for i in postOVVout.Electrolyte.values
            ]
            postOVVout["pH"] = pH_lst
            #        SerieOVV[SampleLabel_cols],
            selEC = pd.merge(
                Prep_EA_BET[SampleSelection.InterestingCols],
                postOVVout,
                how="left",
                on=["SampleID"],
            )
            #            selEC = selEC.T.drop_duplicates().T
            selEC.to_excel(file_EC)
        else:
            selEC = pd.read_excel(file_EC, index_col=[0])
            #            selEC = selEC.T.drop_duplicates().T
            if selEC.empty:
                print("Loaded empty file... restart with force reload...")
                selEC = self.EC_selection_ovv_prep(force_reload=True)
            else:
                print(
                    "Loaded EC prep OVV file...{0}\n with length {1}".format(
                        file_EC, len(selEC)
                    )
                )
        #        selEC = pd.merge(,postOVVout,on=['SampleID','SeriesID','Colorcode'],how='left')
        return selEC
        # IMPORTANT Adaptations to PostEC overview.... implement later to Constructor PostEC method!!

    #        Prep_EA_EC  = pd.merge(EA_results,postOVVout,on='SampleID',how='outer')
    #        selEC = pd.merge(SerieOVV,postOVVout,on=['SampleID','SeriesID','Colorcode'],how='left')
    #%% ===== Run PAR_DW script over selection of Samples ======
    #    def EC_run_PAR_DW_selection(self):
    #        selEC = self.selEC
    #        DD_Series_EC = self.DD_Series_EC
    #        prep = self.Prep_EA_BET
    #        a= SampleSelection(group_selection='Porph_SiO2',type_of_selection='SeriesID',extra_selection_outputs=['EC'])
    #        selEC,DD_Series_EC,prep  =a.selEC, a.DD_Series_EC,a.Prep_EA_BET
    ##        file_EC = DD_Series_EC.joinpath('EC_ORR_kin_pars.xlsx')
    #
    #        runOVV = SampleSelection.EC_load_RunOVV().sort_values(by='expDT',ascending=False)
    #        runSeries = pd.merge(prep[SampleSelection.SampleLabel_cols],runOVV, on=['SampleID'], how='left')
    #        runspecial = runSeries.loc[runSeries.SampleID.str.contains('PT*|JOS1|JOS4')]
    ##        runspecial = runSeries.loc[runSeries.SampleID.str.contains('PT*')]
    ##        runspecial = runSeries.loc[runSeries.SampleID.str.contains('JOS1|JOS4')]
    #        for dt, gr in runspecial.groupby(by=['EXP_date', 'SampleID']):
    #            print(dt[0],dt[1])
    #            pp = runOVV.query('EXP_date == @dt[0] & ((SampleID == @dt[1])| (PAR_exp == "RHE"))')
    ##            pp = runOVV.query('EXP_date == 20190215 & ((SampleID == "JOS8")| (PAR_exp == "RHE"))')
    #            PAR_DW_v18.runEC.EC_Analysis_Run(pp,runOVV)
    #        print('Ran over all PAR files from Series')

    #%% ===== Plot and Output EC ORR kin pars ======
    def EC_retrieve_ORR_pars(self):
        selEC = self.selEC
        DD_Series_EC = self.DD_Series_EC
        prep = self.Prep_EA_BET
        #        selEC,DD_Series_EC,prep  =a.selEC, a.DD_Series_EC,a.Prep_EA_BET
        file_EC = DD_Series_EC.joinpath("EC_ORR_kin_pars.xlsx")
        #        selEC.query('Type_Exp == "ORR_RRDE"').drop_duplicates()
        #        runOVV = pd.read_excel(FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('RunOVV.xlsx'),index_col=[0])
        #        runOVV['expDT'] = pd.to_datetime(runOVV['EXP_date'])
        #        runOVV.rename(columns={'EXP_dir' : 'Exp_dir'},inplace=True)
        ORR_exp = runOVV.query('(SampleID == "JOS2") & (PAR_exp == "ORR")')

        statii, ast_type, ovv = [], [], []
        ORR_pars = (
            selEC.query('Type_Exp == "ORR_pars"')
            .drop_duplicates(subset="SourceFilename")
            .sort_values(by="expDT", ascending=False)
        )
        ORR_pars["Gas"] = "O2"
        merge_cols = ["SampleID", "Electrolyte", "Gas", "Exp_dir", "expDT", "pH"]
        #        [i for i in ORR_pars.columns if i in (runOVV.columns) and not 'EXP_date' in i]
        t1, t2 = ORR_pars[merge_cols].head(5), runOVV.loc[
            (runOVV.expDT == ORR_pars[merge_cols].head(1).expDT.values[0])
            & (runOVV["Gas"] == "O2")
            & (runOVV["PAR_exp"] == "ORR")
        ].head(5)
        #        ORR_pars_OVV = pd.merge(ORR_pars,runOVV,how='left',on=merge_cols).fillna(method='ffill')
        """Merge the selection of samples in series with the actual OVV of all experimental PAR files """
        ORR_pars_OVV = pd.merge(
            ORR_pars,
            runOVV,
            how="left",
            on=merge_cols,
            sort=False,
            suffixes=["_Xdrop", ""],
        )
        ORR_pars_OVV = ORR_pars_OVV.drop(
            columns=[i for i in ORR_pars_OVV.columns if "_Xdrop" in i]
        )
        #        for i,n in ORR_pars.iterrows():
        #            SF = n.SourceFilename
        #            pSF = Path(SF)
        #            sID = FileHelper.FindSampleID.try_find_sampleID(pSF.stem)
        #            lst = [i.stem for i in list(pSF.parent.rglob('*{0}*xlsx'.format(sID[0])))]
        #            ASTs_matches = ([re.match('(.*|_?)(AST)',i) for i in lst])
        #            if any(ASTs_matches):
        #                statii = 'EoL'
        #                ovv = list(pSF.parent.parent.rglob('*OVV*xlsx'.format(sID[0])))
        #                if ovv:
        #                    ovvf = pd.read_excel(ovv[0],index_col=[0])
        #        ORR_pars['Gas'] = 'O2'
        out = []
        for nm, gr in ORR_pars_OVV.groupby(by=["SampleID", "SourceFilename"]):
            #            print(nm[1])
            orr_par_file = pd.read_excel(nm[1], index_col=[0])
            orr_par_file["SourceFilename"] = nm[1]
            for fl, fgr in orr_par_file.groupby(by="File"):
                st, xtr = PostEC.check_status(fl)
                fgr["Status"] = st
                fgr["Status_extra"] = xtr
                out.append([fgr])
        #            if len(orr_par_file.File.unique()) > 1:
        #                print(nm,orr_par_file.File.unique())
        #            for fl,fgr in orr_par_file.groupby(by='File'):
        #                ORR_pars_OVV.loc[(ORR_pars_OVV['PAR_file'] == fl)]
        ##            orr_par_file['Status'] = nm[1]
        #                out.append([orr_par_file])
        #            pars_1500 = par_file.query('RPM == 1500')
        #            if pars_1500.empty:
        #                out.update({'SampleID' : nm[0], 'SourceFilename' : nm[1], '})
        #            else:
        orr_kin_pars_raw = pd.concat([i[0] for i in out], axis=0, sort=False)
        orr_kin_pars_raw.rename(columns={"File": "PAR_file"}, inplace=True)
        orr_kin_pars_raw["Gas"] = "O2"
        orr_kin_pars = pd.merge(orr_kin_pars_raw, prep, how="left", on="SampleID")
        #        ORR = pd.merge(ORR_pars,orr_kin_pars,on=[['SourceFilename','SampleID']+SampleSelection.EC_exp_cols])
        pars_merge_cols = [
            "PAR_file",
            "SourceFilename",
            "SampleID",
            "Electrolyte",
            "RPM",
            "pH",
            "Gas",
        ]
        #        ORR = pd.merge(ORR_pars_OVV,orr_kin_pars,on=pars_merge_cols,how='outer')
        ORR = pd.merge(
            ORR_pars_OVV,
            orr_kin_pars,
            on=pars_merge_cols,
            how="outer",
            suffixes=["_ovv", ""],
        )
        #        for n,ph in orr:
        #        orr_kin_1500 = orr_kin_pars.query('RPM == 1500')
        #        orr_kin_1500.plot(x='E_onset',y='Jkin_075',kind='scatter',s=80)
        ORR1500 = ORR.query("RPM > 1000")
        fig, ax = plt.subplots()
        for nm, gr in ORR1500.groupby(by=["SampleID", "pH"]):
            gr.plot(
                x="E_onset",
                y="Jkin_075",
                kind="scatter",
                s=80,
                c="red",
                label="{0}_{1}".format(*nm),
                ax=ax,
            )
        #        ORR.query('pH <= 1').plot(x='E_onset',y='Jkin_075',kind='scatter',s=80)
        fig, ax = plt.subplots()
        acid_BoL = ORR1500.query('(pH <= 1) & (Status == "BoL")')
        ORR1500.plot(
            x="E_onset",
            y="Jkin_075",
            kind="scatter",
            s=80,
            c="red",
            ax=ax,
            label="acid",
        )
        acid_BoL.to_excel(file_EC.parent.joinpath(file_EC.stem + "_acid_BoL" + ".xlsx"))
        for nm, gr in acid_BoL.groupby(by=["SampleID", "Status"]):
            gr.plot(
                x="E_onset",
                y="Jkin_075",
                kind="scatter",
                s=80,
                c="red",
                label="{0}_{1}".format(*nm),
            )
        alk = ORR1500.query("pH >= 12")
        alk.plot(
            x="E_onset",
            y="Jkin_075",
            kind="scatter",
            s=80,
            c="b",
            ax=ax,
            label="alkaline",
        )

        ORR.to_excel(file_EC)
        return ORR

    #%% ===== Plot and Output EC EIS pars ======
    def EC_retrieve_EIS_pars(self):
        selEC = self.selEC
        DD_Series_EC = self.DD_Series_EC
        prep = self.Prep_EA_BET
        #        selEC,DD_Series_EC,prep  =a.selEC, a.DD_Series_EC,a.Prep_EA_BET
        file_EC = DD_Series_EC.joinpath("EC_EIS_pars.xlsx")

        runOVV = SampleSelection.EC_load_RunOVV()
        EIS_pars = selEC.query('Type_Exp == "EIS"').drop_duplicates(
            subset="SourceFilename"
        )
        merge_cols = ["SampleID", "Electrolyte", "Gas", "Exp_dir", "expDT", "pH"]
        #        [i for i in ORR_pars.columns if i in (runOVV.columns) and not 'EXP_date' in i]
        #        t1,t2 = ORR_pars[merge_cols].head(5),runOVV.loc[(runOVV.expDT == ORR_pars[merge_cols].head(1).expDT.values[0]) & (runOVV['Gas'] == 'O2') & (runOVV['PAR_exp'] == 'ORR') ].head(5)
        #        merge_cols =[i for i in EIS_pars.columns if i in (runOVV.columns) and not 'EXP_date' in i]
        EIS_pars_OVV = pd.merge(
            EIS_pars, runOVV, how="left", on=merge_cols, suffixes=["_Xdrop", ""]
        )

        out = []
        for nm, gr in EIS_pars_OVV.drop_duplicates(subset="SourceFilename").groupby(
            by=["SampleID", "SourceFilename"]
        ):
            #            print(nm[1])
            par_file = pd.read_excel(nm[1], index_col=[0])
            par_file["SourceFilename"] = nm[1]
            for fl, fgr in par_file.groupby(by="File"):
                st, xtr = PostEC.check_status(fl)
                fgr["Status"] = st
                fgr["Status_extra"] = xtr
                print(
                    nm[0],
                    Path(nm[1]).name,
                    fgr.SampleID.unique(),
                    Path(fl).name,
                    st,
                    xtr,
                )
                out.append([fgr])
        #            out.append([par_file])
        #            pars_1500 = par_file.query('RPM == 1500')
        #            if pars_1500.empty:
        #                out.update({'SampleID' : nm[0], 'SourceFilename' : nm[1], '})
        #            else:
        pars_raw = pd.concat([i[0] for i in out], axis=0, sort=False)
        pars_raw = pars_raw_dup.drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        pars_raw["Qad+Cdlp"] = pars_raw["Qad"] + pars_raw["Cdlp"]
        pars_raw.rename(columns={"File": "PAR_file"}, inplace=True)
        pars = pd.merge(pars_raw, prep, how="left", on="SampleID")
        pars_merge_cols = [
            "PAR_file",
            "SourceFilename",
            "SampleID",
            "Electrolyte",
            "RPM",
            "pH",
            "Gas",
            "Status",
            "Status_extra",
        ]
        #        === Take care of the OVV and pars merger !!! ===
        #        EIS_raw = pd.merge(EIS_pars_OVV,pars_raw,on=pars_merge_cols,how='outer')
        #        EIS_plus = pd.merge(EIS_raw,prep,how='left',on='SampleID')
        #        EIS = pd.merge(pars,EIS_pars_OVV,on=pars_merge_cols,how='left')
        #        === Take care of the OVV and pars merger !!! ===
        pars_acid_BoL = pars.query(
            '(pH <= 1) & (Status == "BoL") & (Gas == "O2") & (RPM > 1000)  & (RedChisqr2 < 0.1E-03)'
        )
        #        .dropna(axis=1,how='all')
        try:
            pars_acid_BoL.query("(0.58 < E_AppV_RHE < 0.61)").plot(
                x="BET_cat_agg", y="Qad+Cdlp", kind="scatter"
            )
            for st, gr in EIS.groupby(by="Status"):
                gr.plot(x="E_AppV_RHE", y="Qad", kind="scatter")
        except Exception as e:
            print("EIS plot fail {0}".format(e))

        rcorr = pars_acid_BoL[
            SampleSelection.InterestingCols
            + SampleSelection.EC_EIS_par_cols
            + SampleSelection.RAMAN_cols_corr
            + ["Qad+Cdlp"]
        ].corr()
        plt.subplots(figsize=(20, 15))
        sns.heatmap(rcorr)
        plt.savefig(
            DD_Series_EC.joinpath("EIS_heatmap.png"), dpi=300, bbox_inches="tight"
        )
        plt.close()

        #        EIS = EIS.T.drop_duplicates().T
        #        orr_kin_1500 = orr_kin_pars.query('RPM == 1500')
        #        orr_kin_1500.plot(x='E_onset',y='Jkin_075',kind='scatter',s=80)
        EIS.to_excel(file_EC)
        EIS_acid_BoL = EIS.query('(pH <= 1) & (Status == "BoL") & (Gas == "O2")')
        EIS_acid_BoL.plot(x="E_AppV_RHE", y="Qad", kind="scatter", ylim=[0, 1e-04])
        smpl = ["SampleID"]
        #        pd.merge(acid_BoL,EIS_acid_BoL,how='outer',on=smpl)
        acid_merge = pd.merge(
            acid_BoL, EIS_acid_BoL, how="left", on=SampleSelection.Inter_Prep_cols
        )

        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        for nm, gr in acid_merge.groupby(by="SampleID"):
            gr.plot(x="ColorCode", y="Jkin_075", ax=ax)

        ovm = [i for i in ORR.columns if i in (EIS.columns) and not "EXP_date" in i]
        ovm = ["SampleID", "Electrolyte", "pH", "Gas", "Status", "Status_extra"]
        EIS[ovm]
        ORR[ovm]
        EIS_ORR = pd.merge(EIS, ORR, on=ovm, how="left")
        ORR_EIS = pd.merge(EIS, ORR, on=ovm, how="right")

        #        for a,gr

        GR = EIS_ORR.groupby(
            by=(SampleSelection.InterestingCols + SampleSelection.EC_exp_cols)
        )
        for n in EIS_ORR.groupby(
            by=(SampleSelection.InterestingCols + SampleSelection.EC_exp_cols)
        ):
            print(n)

        #%%
        #        EIS.query('(0.65 < E_AppV_RHE < 0.78) & (Gas == "O2")').Rct.values, ORR.query('RPM == 1500').Jkin_075.values
        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        EIS_pars.query('(0.65 < E_AppV_RHE < 0.78) & (Gas == "O2")').plot(
            x="Colorcode", y="Qad", kind="scatter", ax=ax, c="r", s=80, ylim=(0, 1e-3)
        )
        ORR.query("RPM > 1000").dropna(subset=["Jkin_075"]).plot(
            x="Colorcode", y="Jkin_075", kind="scatter", ax=ax2, c="b", s=80
        )
        #        plt.scatter()
        #%%

        return EIS

    #%% ===== Plot and Output EC N2 CV and Cdl pars ======
    def EC_retrieve_N2_pars(self):
        selEC = self.SerieOVV_EC

        #        selEC.query('Type_Exp == "N2_CV"').drop_duplicates()

        N2_CV_pars = selEC.query('Type_Exp == "N2_CV"').drop_duplicates(
            subset="SourceFilename"
        )
        N2_Cdl_pars = selEC.query('Type_Exp == "N2_Cdl"').drop_duplicates(
            subset="SourceFilename"
        )
        out = []
        for nm, gr in N2_CV_pars.groupby(by=["SampleID", "SourceFilename"]):
            print(nm[1])
        for nm, gr in N2_Cdl_pars.groupby(by=["SampleID", "SourceFilename"]):
            print(nm[1])
            N2_Cdl_file = pd.read_csv(nm[1], index_col=[0])
        #            out.append([orr_par_file])
        #            pars_1500 = par_file.query('RPM == 1500')
        #            if pars_1500.empty:
        #                out.update({'SampleID' : nm[0], 'SourceFilename' : nm[1], '})
        #            else:
        orr_kin_pars = pd.concat([i[0] for i in out], axis=0)
        orr_kin_1500 = orr_kin_pars.query("RPM == 1500")
        orr_kin_1500.plot(x="E_onset", y="Jkin_075", kind="scatter", s=80)
        return N2_kin_pars

    #%% ===== Plot and Output EC Tables ======

    def Retrieve_EC_Table_from_OVVs(self):
        #        Porph = SampleSelection('Porph_SiO2','SeriesID',[]).EC_selection_ovv_prep()
        #        PorphSiO2 = SampleSelection('Porph_SiO2','SeriesID',['EC_select'])
        #        DD_Series = PorphSiO2.DD_Series
        #        SerieOVV_EC = self.EC_sel
        SerieOVV_EC = self.Prep_EA_BET
        lst = []
        EC_ovvs_fp = list(
            FileHelper.FindExpFolder("VERSASTAT").DestDir.rglob("OVV_*xlsx")
        )
        for i in EC_ovvs_fp:
            try:
                lst.append(pd.read_excel(i))
            except:
                pass
        EC_ovv_all = pd.concat(lst, sort=False)
        #        EC_ovvs_fp = list(FileHelper.FindExpFolder('VERSASTAT').DestDir.rglob('OVV_*xlsx'))
        SerieOVV_EC_all = pd.merge(
            SerieOVV_EC[SampleSelection.InterestingCols],
            EC_ovv_all,
            how="left",
            on="SampleID",
            suffixes=["", ""],
        )
        status, extra = PostEC.postEC_Status(SerieOVV_EC_all.basename.values)
        SerieOVV_EC_all = SerieOVV_EC_all.assign(
            **{"Status": status, "Status_extra": extra}
        )
        #        SerieOVV_EC_all['Status_extra'] = extra
        return SerieOVV_EC_all

    #    @staticmethod
    def Retrieve_EC_Table_ExpFiles(self):
        SerieOVV_EC = self.EC_sel
        DD_Series = self.DD_Series
        #        SerieOVV.SampleID.unique()
        #        Serie_EC = postOVVout.loc[postOVVout.SampleID.str.contains('|'.join(SerieOVV.SampleID.unique())),:]
        SerieOVV_EC_all = self.SerieOVV_EC_all
        #        SerieOVV_EC_all = pd.merge(SerieOVV_EC[SampleSelection.InterestingCols],EC_ovv_all,how='left',on='SampleID',suffixes=['',''])
        #        SerieOVV_EC_all['Status'] = PostEC.postEC_Status(SerieOVV_EC_all.basename.values)
        EC_lst = []
        print("Starting to create table from Electrochemical measurements...")
        #        SerieOVV_EC_all
        for n, gr in SerieOVV_EC_all.loc[
            SerieOVV_EC_all.Electrolyte != "False"
        ].groupby(by="SampleID"):
            smlbl = gr.SampleLabel.unique()[0]
            serID = gr.SeriesID.unique()[0]
            if gr.empty:
                continue
            elif gr.dropna(subset=["Electrolyte"]).empty:
                continue
            for En, Egr in gr.groupby(by="Electrolyte"):
                for STn, STgr in Egr.groupby(by="Status"):
                    for Gn, Ggr in STgr.groupby(by="Gas"):
                        #                        Un_Exp_pre = Ggr.loc[Ggr.Type_Exp.str.contains('N2|EIS_Combined|ORR_RRDE|HER|OER'),:].drop_duplicates(subset='SourceStem' )
                        Un_Exp_pre = Ggr.loc[
                            (Ggr.PAR_exp.str.contains("N2|EIS|ORR|HER|OER|AST"))
                            & (Ggr.Electrode != "Pt_ring"),
                            :,
                        ].drop_duplicates(subset="basename")
                        #                        Unique_Exps =Un_Exp_pre['Type_Exp'].value_counts()
                        ovvs_sel = Un_Exp_pre.loc[
                            (Un_Exp_pre.Electrolyte == En)
                            & (Un_Exp_pre.SampleID == n)
                            & (Un_Exp_pre.Gas == Gn),
                            :,
                        ]
                        Unique_Exps = ovvs_sel["PAR_exp"].value_counts()
                        Unique_Dates = ", ".join(ovvs_sel.EXP_date.astype(str).unique())
                        Uniq_ExpDir = ovvs_sel.EXP_dir.unique()
                        if STn == "EoL" or STn == "BoL":
                            EoL_dirs = Uniq_ExpDir
                            #                            ovv_files =[list(FileHelper.FindExpFolder('VERSASTAT').DestDir.rglob('*/*%s*/*OVV_%s*xlsx'%(Path(i).name,Path(i).name))) for i in EoL_dirs]
                            #                            if ovv_files:
                            #                                ovvs = pd.concat([pd.read_excel(i[0]) for i in ovv_files if i])
                            #                            ovvs = Ggr

                            #                                ovvs_sel = ovvs.loc[(ovvs.Electrolyte == En) & (ovvs.SampleID == n),:]
                            ASTs_matches = [
                                re.match("(?<!(post))AST", i)
                                for i in ovvs_sel.basename.values
                            ]
                            LC_matches = [
                                re.match("(_LC_)", i) for i in ovvs_sel.basename.values
                            ]
                            #                                ASTs_matches
                            if any(ASTs_matches + LC_matches):
                                #                                    matches = [i for i in ASTs_matches+LC_matches if i != None]
                                if Gn == "N2":
                                    extraST = "_SST"
                                elif Gn == "O2":
                                    extraST = "_LC"
                            else:
                                extraST = ""
                            ASTs = ovvs_sel[ovvs_sel.basename.str.contains("_AST|_LC_")]
                            ORR_files = Un_Exp_pre.PAR_exp.str.contains("ORR")
                            #                                ovvs_sel.loc[ovvs_sel.PAR_exp.str.contains('AST')]
                            #                                Un_AST = ASTs.PAR_exp.unique()
                            STn = STn + extraST
                            #                                (','.join(Un_AST))
                            Unique_AST_names = ", ".join(ASTs.basename.unique())
                        else:
                            Unique_AST_names = ""
                        #                            [list(i.glob('*AST')) for i in Eol_dirs]
                        #                            list(Path(EoL_dirs[0]).rglob(r'(?<!(post|VERS))AST'))
                        if not Unique_Exps.empty and En != "FALSE":
                            EC_row = {
                                "SeriesID": serID,
                                "SampleLabel": smlbl,
                                "SampleID": n,
                                "Electrolyte": En,
                                "Status": STn,
                                "Gas": Gn,
                                "Dates": Unique_Dates,
                                "Type_AST": Unique_AST_names,
                            }
                            EC_row.update(Unique_Exps.to_dict())
                            #                            EC_row = [n,En,STn,Gn,'%s'%Unique_Exps.to_dict()]
                            #                            print(EC_row)
                            EC_lst.append(EC_row)
        EC_Table = pd.DataFrame(EC_lst).fillna(0)
        if not EC_Table.empty:
            EC_Table_indx = EC_Table.set_index(
                ["SeriesID", "SampleLabel", "SampleID", "Electrolyte", "Status", "Gas"]
            )
            EC_Table_indx.to_excel(
                DD_Series.parent.joinpath("EC_Table_%s" % DD_Series.name).with_suffix(
                    ".xlsx"
                )
            )
        else:
            print("EC table is empty")
            EC_Table_indx = pd.DataFrame()
        #        EC_Table_indx.to_excel(DD_Series.parent.joinpath('EC_Table_All').with_suffix('.xlsx'))
        return EC_Table_indx

    #%% ===== EC Load run ovv ======
    def EC_load_RunOVV(PAR_version=18):
        RunOVV_fn = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
            "RunOVV_v{0}.xlsx".format(PAR_version)
        )
        runOVV = pd.read_excel(RunOVV_fn, index_col=[0])
        runOVV["expDT"] = pd.to_datetime(runOVV["EXP_date"])
        runOVV.Date_PAR_EXP = pd.to_timedelta(runOVV.Date_PAR_EXP, unit="h")
        if FileHelper.FindExpFolder("BET").DataDir.drive != (
            Path(runOVV.EXP_dir.iloc[0]).drive and Path(runOVV.PAR_file.iloc[0]).drive
        ):
            CS_parts_PDD = FileHelper.FileOperations.find_CS_parts(
                FileHelper.FindExpFolder("BET").DataDir
            )
            chLst_ExpDir = [
                CS_parts_PDD[0].joinpath(FileHelper.FileOperations.find_CS_parts(i)[1])
                for i in runOVV.EXP_dir.values
            ]
            chLst_PARfile = [
                CS_parts_PDD[0].joinpath(FileHelper.FileOperations.find_CS_parts(i)[1])
                for i in runOVV.PAR_file.values
            ]
            runOVV = runOVV.assign(
                **{
                    "EXP_dir": chLst_ExpDir,
                    "Exp_dir": chLst_ExpDir,
                    "PAR_file": chLst_PARfile,
                }
            )

        #        runOVV.rename(columns={'EXP_dir' : 'Exp_dir'},inplace=True)
        elec = [
            FileHelper.FindSampleID.determine_pH_from_filename(i)["Electrolyte"][0]
            for i in runOVV.EXP_dir.values
        ]
        status, extra = PostEC.postEC_Status(
            [Path(i).stem for i in runOVV.PAR_file.values]
        )
        runOVV = runOVV.assign(
            **{
                "Electrolyte": elec,
                "PAR_file_stem": [Path(i).stem for i in runOVV.PAR_file.values],
                "Status": status,
                "Status_extra": extra,
            }
        )
        #        runOVV['Electrolyte'] = elec

        #        FileHelper.FindExpFolder('BET').DataDir.drive,xyPath.drive

        return runOVV

    #%% ===== EC EIS series ======
    def SampleSelection_EC_EIS_serie(selEC, group_selection):
        #        group_selection = 'Porph_SiO2'
        selEC1 = SampleSelection("Porph_SiO2", "SeriesID", []).EC_selection_ovv_prep()
        selEC2 = SampleSelection("Co_PANI", "SeriesID", []).EC_selection_ovv_prep()
        selEC3 = SampleSelection("CB_bim", "SeriesID", []).EC_selection_ovv_prep()
        OriginColors = SampleSelection.OriginColorList()
        EIS_DD_Series = FileHelper.FindExpFolder("VERSASTAT").DestDir.parent.joinpath(
            "Series", "%s" % group_selection, "EIS"
        )
        EIS_DD_Series.mkdir(parents=True, exist_ok=True)
        selEC.dropna(subset=["SourceFilename"], inplace=True)
        EC_grp_cols = [
            "SampleID",
            "Electrolyte",
            "Gas",
            "RPM",
            "Scanrate",
            "Status",
            "SourceFilename",
        ]
        EC_EIS_grp_cols = [
            "SampleID",
            "Electrolyte",
            "Gas",
            "RPM",
            "Status",
            "SourceFilename",
        ]
        #        DD_Series = FileHelper.FindExpFolder('VERSASTAT').DestDir.parent.joinpath('Series', '%s'%SeriesID)
        selEC.query('Type_Exp == "EIS"').columns
        acid1500o2 = selEC.query(
            'RPM == 1500 & Gas == "O2" & Electrolyte == "0.1MH2SO4"'
        )
        EC_EIS_pars2 = selEC.loc[
            (selEC["Type_Exp"] != "EIS_Combined")
            & (selEC.SourceFilename.str.contains("EIS|_pars2") == True),
            :,
        ]
        EC_EIS_pars2 = selEC.loc[(selEC["Type_Exp"] == "EIS_Combined"), :]
        selEC.loc[(selEC.SourceFilename.str.contains("_pars2") == True), :]
        EC_EIS_pars2 = selEC.query(
            'RPM == 1500 & Gas == "O2" & Electrolyte == "0.1MH2SO4"'
        )
        for El, Egr in selEC.query('Type_Exp == "EIS"').groupby(by="Electrolyte"):
            for gas, Gas_gr in Egr.groupby(by="Gas"):
                EC_out_df = []
                for rpm, rpm_gr in Egr.groupby(by="RPM"):
                    for sID, sID_gr in Egr.groupby(by="SampleID"):
                        EIS_pars = pd.read_excel(fn, drop_index=True)
                for fn, fGr in Rgr.groupby(by="SourceFilename"):
                    EIS_pars = pd.read_excel(fn, drop_index=True)
        EC_EIS_grp_cols = [
            "SampleID",
            "Electrolyte",
            "Gas",
            "RPM",
            "Status",
            "SourceFilename",
        ]
        #%%
        ## MISSING PARS2 files from measuremetns!! Were they created/ where did they go??

    #        for y_col in ['Qad','Cdlp','Rct','Rct_kin','Rs','Rorr']:
    def Read_EIS_ParFiles_to_OVV(selEC, group_selection):
        OriginColors = SampleSelection.OriginColorList()
        EIS_DD_Series = FileHelper.FindExpFolder("VERSASTAT").DestDir.parent.joinpath(
            "Series", "%s" % group_selection, "EIS"
        )
        EIS_DD_Series.mkdir(parents=True, exist_ok=True)
        EC_EIS_pars1 = selEC1.loc[(selEC1["Type_Exp"] == "EIS_Combined"), :]
        EC_EIS_pars2 = selEC2.loc[(selEC2["Type_Exp"] == "EIS_Combined"), :]
        #        .EC_selection_ovv_prep()
        EC_EIS_pars_all = selEC4.loc[
            (selEC4["Type_Exp"] == "EIS_Combined") & (selEC4["SampleID"] == "DW25"), :
        ]
        EC_EIS_pars_all = selEC_all.loc[(selEC_all["Type_Exp"] == "EIS_Combined"), :]
        EC_EIS_pars = pd.concat([EC_EIS_pars1, EC_EIS_pars2])
        all_EC_pars, case_EC_pars = [], []
        for y_col in SampleSelection.EC_EIS_par_cols[:-2]:
            #        for y_col in ['Qad']:
            for nms, gr in EC_EIS_pars_all.groupby(
                by=["Electrolyte", "Gas", "RPM", "Status"]
            ):
                #            ['Electrolyte','Gas','RPM','Status','SampleID' ,'SourceFilename']):
                print(", ".join([str(i) for i in list(nms)]))
                case_EC_pars, y_maxes = [], []
                #                print(gr[['SampleLabel','SourceFilename']].values)
                #%%
                Cdl_plotted = []
                fig, ax = plt.subplots()
                fig.suptitle(
                    "%s " % y_col
                    + "{0} with {1} at {2}, {3}".format(
                        nms[0], nms[1], str(int(nms[2])), nms[3]
                    )
                )
                nms_path_png = EIS_DD_Series.joinpath(
                    "{0}/{1}_{2}_{3}_{4}.png".format(
                        y_col, nms[0], nms[1], str(int(nms[2])), nms[3]
                    )
                )
                for IDnm, IDgr in gr.groupby(by=["SampleID", "SourceFilename"]):
                    #                        selEC_slice = selEC.loc[((selEC.SampleID == IDnm[0]) & (selEC.Electrolyte == nms[0]) & (selEC.Gas == nms[1]) & (selEC.Status == nms[-1])),:]
                    pars2 = pd.DataFrame()
                    ID_colorcode = [
                        i for i in IDgr.Colorcode.unique() if str(i) != "nan"
                    ][0]

                    RGB = [
                        float(i) / 255
                        for i in (
                            OriginColors.loc[
                                OriginColors.OriginCode == ID_colorcode, "RGB"
                            ].values[0]
                        ).split(",")
                    ]
                    #                RGB = OriginColors.loc[OriginColors.OriginCode == IDgr.Colorcode.unique()[0],'color'].values[0]
                    pars2_SourceFile = list(
                        Path(IDnm[1]).parent.rglob(
                            "%s_pars2*" % Path(IDnm[1]).stem.split("_Combined")[0]
                        )
                    )
                    if len(pars2_SourceFile) > 1:
                        print(
                            "Multiple source files pars2: ",
                            ", ".join([Path(i).stem for i in pars2_SourceFile]),
                        )
                    elif not pars2_SourceFile:
                        print(IDnm, pars2_SourceFile)
                        continue
                    pars2 = pd.read_excel(pars2_SourceFile[0])
                    try:
                        pars2.RPM.astype(float, inplace=True)
                    except:
                        continue
                    if y_col not in pars2.columns:
                        continue
                    #            pars2.query('SampleID == @nms[0] & Electrolyte == @nms[1] & Gas == @nms[2] & RPM == @nms[3]')
                    pars2_grp = pars2.loc[
                        (pars2.SampleID == IDnm[0])
                        & (pars2.Electrolyte == nms[0])
                        & (pars2.Gas == nms[1])
                        & (pars2.RPM == nms[2]),
                        :,
                    ]
                    pars2_grp["SourceStem"] = [Path(i).stem for i in pars2_grp.File]
                    if y_col in "Rorr":
                        pars2_grp = pars2_grp.query("Rorr < 1E5")
                    if pars2_grp.empty:
                        print(
                            "EIS pars empty:",
                            "%s " % y_col + "{0} with {1} at {2}, {3}".format(*nms),
                        )
                        print("%s " % y_col + "{0} in {1}".format(*IDnm))
                        continue

                    if "BoL" in nms[-1]:
                        pars2_grp = pars2_grp.loc[
                            pars2_grp.SourceStem.str.contains("AST") == False, :
                        ]
                    elif "EoL" in nms[-1]:
                        pars2_grp = pars2_grp.loc[
                            pars2_grp.SourceStem.str.contains("AST") == True, :
                        ]

                    if pars2_grp.empty:
                        continue

                    Erhe = pars2_grp["E_AppV_RHE"].values
                    y = pars2_grp[y_col].values
                    #                    pars2_add_cols = pars2_add_cols = dict(zip(['Electrolyte_add','Gas_add','RPM_add','Status_add']+['SampleID_add','SourceFilename_add','Pars2_SourceFile_add'],list(nms)+[*IDnm,pars2_SourceFile]))
                    #                    dict(zip(['Electrolyte_add','Gas_add','RPM_add','Status_add']+['SampleID_add','SourceFilename_add','Pars2_SourceFile_add']
                    #                    ,[[i]*len(pars2_grp) for i in list(nms)]+[[i]*len(pars2_grp) for i in IDnm],[[i]*len(pars2_grp) for i in pars2_SourceFile]))
                    #                    pars2_grp.assign(**pars2_add_cols)
                    EC_pars = pd.merge(
                        pars2_grp,
                        gr,
                        on=["Electrolyte", "Gas", "RPM", "SampleID", "SourceStem"],
                        how="left",
                    )

                    all_EC_pars.append(EC_pars)
                    case_EC_pars.append(EC_pars)
                    #                pars2_grp.plot(x='E_AppV_RHE',y='Qad',label='{0}'.format(IDnm[0]),kind='scatter',ax=ax,c=RGB)
                    ax.scatter(
                        Erhe,
                        y,
                        label="{0},{1}".format(IDnm[0], Path(pars2_SourceFile[0]).stem),
                        color=RGB,
                    )
                    ym = y.max() * 1.3 if y.max() < 2 * y.mean() else y.mean() * 1.5
                    y_maxes.append(ym)
                    #                    ax.set_ylim(0,ym)
                    ax.set_xlim(0, 1.2)
                    ax.legend(
                        bbox_to_anchor=(-0.3, 1.15, 1.8, 0.102),
                        loc="lower left",
                        ncol=1,
                        mode="expand",
                        borderaxespad=0.0,
                    )
                    nms_path = EIS_DD_Series.joinpath(
                        "{0}/{1}/{2}/{3}".format(
                            y_col, nms[0], nms[1], str(int(nms[2]))
                        )
                    )
                    nms_path.mkdir(parents=True, exist_ok=True)
                    IDnm_path = "{0}_{1}_{2}.xlsx".format(
                        IDnm[0], nms[3], "_".join(Path(IDnm[1]).stem.split("_")[-3::])
                    )
                    EC_pars.to_excel(nms_path.joinpath(IDnm_path))

                    if "Qad111" in y_col or "N2222" in nms:
                        selEC_slice = selEC.loc[
                            (
                                (selEC.SampleID == IDnm[0])
                                & (selEC.Electrolyte == nms[0])
                                & (selEC.Gas == "N2")
                                & (selEC.Status == nms[-1])
                                & (selEC.Type_Exp == "N2_Cdl")
                            ),
                            :,
                        ]

                        if not selEC_slice.empty:
                            for Cdl_nm, Cdl_gr in selEC_slice.groupby(
                                by="SourceFilename"
                            ):
                                print(Cdl_nm)
                                Cdl_path = nms_path_png.parent.joinpath(
                                    "Cdl_"
                                    + nms_path_png.stem
                                    + "_"
                                    + "_".join(Path(Cdl_nm).stem.split("_")[-3:])
                                    + ".xlsx"
                                )
                                try:
                                    Cdl_sample = pd.read_csv(Cdl_nm).query(
                                        'Sweep_Type_N2 == "cathodic"'
                                    )
                                    Cdl_sample.to_excel(Cdl_path)
                                except:
                                    Cdl_sample = pd.DataFrame()
                                    continue
                                #                                Cdl_sample.plot(x='E_AppV_RHE',y='Cdl',ax=ax)
                                if not Cdl_nm in Cdl_plotted and not Cdl_sample.empty:
                                    Cdl_y = Cdl_sample["Cdl"].values
                                    ax.scatter(
                                        Cdl_sample["E_AppV_RHE"].values,
                                        Cdl_y,
                                        label="{0}".format(Path(Cdl_nm).stem),
                                        marker="^",
                                        color=RGB,
                                        alpha=0.2,
                                    )
                                    Cdl_plotted.append(Cdl_nm)
                                    Cdl_ym = (
                                        Cdl_y.max() * 1.3
                                        if Cdl_y.max() < 2 * Cdl_y.mean()
                                        else Cdl_y.mean() * 1.5
                                    )

                                    y_maxes.append(Cdl_ym)
                #                                    if Cdl_ym > ym:
                #                                        ax.set_ylim(0,Cdl_ym)

                #                nms_path = EIS_DD_Series.joinpath('_'.join([y_col]+[str(i) for i in list(nms)]))
                #                nms_path_png = EIS_DD_Series.joinpath('{0}/{1}/{2}/{3}/{4}.png'.format(y_col,nms[0],nms[1],str(int(nms[2])),nms[3]))
                if y_maxes:
                    ax.set_ylim(0, np.max(y_maxes))
                else:
                    ax.set_ylim(0, ym)
                #                nms_path_xlsx = EIS_DD_Series.joinpath('{0}/{1}_{2}_{3}_{4}.xlsx'.format(y_col,nms[0],nms[1],str(int(nms[2])),nms[3]))
                plt.show()
                #%%
                plt.savefig(nms_path_png, dpi=300, bbox_inches="tight")
        #                plt.close()
        #                pd.concat([i for i in case_EC_pars]).drop_duplicates().to_excel(nms_path_xlsx)
        EIS_pars_char = pd.concat([i for i in all_EC_pars]).drop_duplicates(
            subset=SampleSelection.EC_EIS_par_cols
        )
        EIS_pars_char["E_Vrhe"] = EIS_pars_char["E_AppV_RHE"].round(2)
        EIS_pars_char.to_excel(EIS_DD_Series.joinpath("EIS_pars_char_new_all.xlsx"))
        #        for E,grE in EIS_pars_char.groupby(by='E_Vrhe'):
        #%%
        EIS_pars_char = pd.read_excel(
            EIS_DD_Series.parent.parent.joinpath("EIS_pars_char_new_all.xlsx")
        )
        lin_regs, log_regs = [], []
        #
        #        x_char_cols = ['BET_Area_RPT','N_content','100-CHN','t-Method micropore area (m2/g)']
        x_char_cols = ["BET_Area_RPT", "N_content", "100-CHN"]
        Series_pars_char = EIS_pars_char.query('SeriesID == "CB4"')
        #        for x_col in SampleSelection.Character_xcols[2:]:
        for x_col in x_char_cols:
            for y_col in SampleSelection.EC_EIS_par_cols[:-2]:
                for E in np.arange(0, 1.55, 0.3):
                    grE = Series_pars_char.query(
                        "(E_Vrhe < @E+0.03) & (E_Vrhe > @E-0.03) & (Rorr < 1E06)"
                    )
                    print(grE)
                    grE = grE.loc[grE[y_col] < 1e6]
                    #            EIS_pars_char.loc[((EIS_pars_char.E_Vrhe < E+0.003) & (EIS_pars_char.E_Vrhe > E-0.003)),:]
                    if (
                        len(grE) >= 3
                        and grE[y_col].mean() != 0
                        and grE[x_col].mean() != 0
                    ):
                        #                        'Electrolyte','Gas','RPM','Status'
                        for nms, stgr in grE.groupby(
                            by=["Electrolyte", "Gas", "RPM", "Status"]
                        ):
                            #                        for nms, stgr in grE.groupby(by=['Electrolyte','Gas','Status']):
                            #                            nms = (nms[0],nms[1],'1500',nms[-1])
                            stgr = stgr.dropna(subset=[x_col])
                            if len(stgr) >= 3 and stgr[x_col].std() > 0.1:
                                #                                reg = linear_model.LinearRegression()
                                #                                reg.fit(stgr[x_col].values.reshape(1,-1),stgr[y_col].values.reshape(1,-1))
                                reg = sm.OLS(
                                    stgr[y_col], sm.add_constant(stgr[x_col])
                                ).fit()
                                log_reg = sm.OLS(
                                    np.log(stgr[y_col]), sm.add_constant(stgr[x_col])
                                ).fit()
                                lin_regs.append(
                                    [
                                        (x_col, y_col),
                                        reg.rsquared,
                                        reg.params[x_col],
                                        len(stgr),
                                        nms,
                                        stgr.SampleID.unique(),
                                    ]
                                )
                                log_regs.append(
                                    [
                                        (x_col, y_col),
                                        log_reg.rsquared,
                                        log_reg.params[x_col],
                                        len(stgr),
                                        nms,
                                        stgr.SampleID.unique(),
                                    ]
                                )
                                if reg.rsquared > 0.01 or log_reg.rsquared > 0.08:
                                    stem_path = "{0}_{1}_{2}_{3}_{4}_{5}_{6}".format(
                                        nms[0],
                                        y_col,
                                        x_col,
                                        int(E * 1000),
                                        nms[1],
                                        str(int(nms[2])),
                                        nms[3],
                                    )
                                    png_path = stem_path.replace("/", "-") + ".png"
                                    xls_path = stem_path.replace("/", "-") + ".xlsx"
                                    #                                png_fullpath =
                                    #                                    folder_path = EIS_DD_Series.joinpath('E_Vrhe4')
                                    folder_path = Path.home().joinpath("E_Vrhe4")
                                    folder_path.mkdir(parents=True, exist_ok=True)
                                    #                                    Path.home().joinpath('E_Vrhe2').mkdir(parents=True,exist_ok=True)

                                    fig, ax = plt.subplots()
                                    fig.suptitle(
                                        "%s " % y_col
                                        + "{0} with {1} at {2}, {3}".format(
                                            nms[0], nms[1], str(int(nms[2])), nms[3]
                                        )
                                        + " at %.2f V_RHE " % E
                                    )

                                    #                                    ym = stgr[y_col].max()*1.3 if stgr[y_col].max() < 2*stgr[y_col].mean() else stgr[y_col].mean()*1.5
                                    ym = stgr[y_col].max() * 1.3
                                    EC_st_out = []
                                    for IDnm, IDgr in stgr.groupby(
                                        by=["SampleID", "SourceFilename"]
                                    ):
                                        ID_cl1 = [
                                            i
                                            for i in IDgr.Colorcode.unique()
                                            if str(i) != "nan"
                                        ]
                                        if ID_cl1:
                                            ID_colorcode = ID_cl1[0]
                                        else:
                                            ID_colorcode = 1

                                        RGB = [
                                            float(i) / 255
                                            for i in (
                                                OriginColors.loc[
                                                    OriginColors.OriginCode
                                                    == ID_colorcode,
                                                    "RGB",
                                                ].values[0]
                                            ).split(",")
                                        ]
                                        x = IDgr[x_col].values
                                        y = IDgr[y_col].values
                                        ax.scatter(
                                            x, y, label="{0}".format(IDnm[0]), color=RGB
                                        )
                                        EC_st_out.append(IDgr)
                                    pd.concat([i for i in EC_st_out]).to_excel(
                                        folder_path.joinpath(xls_path)
                                    )
                                    ax.set_ylim(0, ym)
                                    #                                    if log_reg.rsquared > reg.rsquared:
                                    #                                        ax.autoscale(True)
                                    #                                        ax.set_yscale('log')

                                    ax.set_xlabel(x_col)
                                    ax.set_ylabel(y_col)
                                    #                                    ax.autoscale(True)
                                    ax.legend(
                                        bbox_to_anchor=(-0.3, 1.15, 1.8, 0.102),
                                        loc="lower left",
                                        ncol=2,
                                        mode="expand",
                                        borderaxespad=0.0,
                                    )

                                    plt.savefig(
                                        folder_path.joinpath(png_path),
                                        dpi=300,
                                        bbox_inches="tight",
                                    )
                                    plt.show()
                                    plt.close()

        #                                 grE.plot(x='BET_Area_RPT',y='Qad',label=E,kind='scatter')
        LinResults = pd.DataFrame(
            lin_regs,
            columns=["cols", "Rsq", "Xcol", "len_gr", "conditions", "SampleIDs"],
        )
        LinResults095 = LinResults.query("Rsq > 0.95").sort_values(
            by=["Rsq", "len_gr"], ascending=False
        )
        #%%
        #        df = pd.DataFrame(np.random.randn(50, 7), columns=list('ABCDEFG'))
        df = EIS_pars_char[
            ["BET_Area_RPT", "N_content"] + SampleSelection.EC_EIS_par_cols[:-2]
        ]
        #        df = EIS_pars_char[SampleSelection.Character_xcols+['Qad']]

        corr2 = df.corr(method="pearson")
        sns.heatmap(corr2)
        # initiate empty dataframe
        corr = pd.DataFrame()
        for a in SampleSelection.Character_xcols + SampleSelection.EC_EIS_par_cols:
            for b in list(df.columns.values):
                corr.loc[a, b] = df.corr().loc[a, b]

            #%%

    #            aa =pars2.loc[(pars2.SampleID == nms[0]) & (pars2.Electrolyte == nms[1]) & (pars2.Gas == nms[2])]
    #            & (pars2.RPM == float(nms[3]))),:]
    def SampleSelection_EC_EIS_serie(selEC, group_selection):
        for Sn, Sgr in EC3.query('Type_Exp == "EIS_Combined"').groupby(by="SeriesID"):
            print(Sgr["SourceFilename"].unique())
            for El, Egr in Sgr.groupby(by="Electrolyte"):
                for ScR, Rgr in Egr.groupby(by="Scanrate"):
                    EC_out_df = []
                    fig, ax = plt.subplots()
                    try:
                        for fn, fGr in Rgr.groupby(by="SourceFilename"):
                            EIS = pd.read_excel(fn, drop_index=True)
                            for Ev, EISgr in EIS.groupby("E_AppV_RHE"):

                                ECout = pd.DataFrame()
                                #                                CV['j_mAcm2'] = CV['j A/cm2']*1E3
                                outlbl = "%s (%s)" % (
                                    fGr.SampleLabel.values[0],
                                    fGr.SampleID.values[0],
                                )
                                outlbl2 = Path(fn).stem
                                fig, ax = plt.subplots()
                                EISgr.plot(
                                    x="DATA_Yre",
                                    y="DATA_Yim",
                                    label="data_%s" % outlbl,
                                    ax=ax,
                                    kind="scatter",
                                )
                                EISgr.plot(
                                    x="FIT1_Yre", y="FIT1_Yim", ax=ax, label="fit"
                                )
                                EISgr.plot(
                                    x="FIT2_Yre", y="FIT2_Yim", ax=ax, label="fit"
                                )
                                ax.axis("equal")
                                ax
                    #                                ECout = CV[['E_AppV_RHE','j_mAcm2']].rename({'E_AppV_RHE' : 'EV_%s'%Path(fn).stem,'j_mAcm2' : 'j_mA_%s'%Path(fn).stem}, axis=1).reset_index(drop=True)
                    #                                EC_out_df.append(ECout)
                    except:
                        pass
                    #                        pd.concat([i for i in EC_out_df],sort=False,axis=1).to_excel(DD_Series.joinpath('EC_%s_%s_%.0f.xlsx'%(El,Sn,int(ScR)))) # export to XLS
                    #                    ax.set_ylim(-1E-4,1E-4)
                    ax.autoscale(True)
                    ax.set_xlim(0, 1.2)
                    ax.legend(
                        bbox_to_anchor=(-0.3, 1.15, 1.8, 0.102),
                        loc="lower left",
                        ncol=1,
                        mode="expand",
                        borderaxespad=0.0,
                    )
                    ax.set_xlabel("$\mathrm{E \//\/ V_{RHE}}$")
                    ax.set_ylabel("$\mathrm{j \//\/ mAcm^{-2}}$")
                    fig.suptitle("CV in N2 at %.2f mV/s" % ScR, y=0.95)
                    plt.savefig(
                        DD_Series.joinpath("EC_%s_%s_%.0f.png" % (El, Sn, int(ScR))),
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close()

    def SampleSelection_EC_N2_CV(selEC, group_selection):
        DD_Series = FileHelper.FindExpFolder("VERSASTAT").DestDir.parent.joinpath(
            "Series", "%s" % group_selection
        )
        DD_Series.mkdir(parents=True, exist_ok=True)
        EC3 = selEC.dropna(subset=["SourceFilename"])
        #        DD_Series = FileHelper.FindExpFolder('VERSASTAT').DestDir.parent.joinpath('Series', '%s'%SeriesID)
        for Sn, Sgr in EC3.query('Type_Exp == "N2_CV"').groupby(by="SeriesID"):
            #                print(Sn,Sgr)
            print(Sgr["SourceFilename"].unique())
            for El, Egr in Sgr.groupby(by="Electrolyte"):
                for ScR, Rgr in Egr.groupby(by="Scanrate"):
                    EC_out_df = []
                    fig, ax = plt.subplots()
                    try:
                        for fn, fGr in Rgr.groupby(by="SourceFilename"):
                            CV = pd.read_excel(fn)
                            ECout = pd.DataFrame()
                            CV["j_mAcm2"] = CV["j A/cm2"] * 1e3
                            outlbl = "%s (%s)" % (
                                fGr.SampleLabel.values[0],
                                fGr.SampleID.values[0],
                            )
                            outlbl2 = Path(fn).stem
                            CV.plot(x="E_AppV_RHE", y="j_mAcm2", ax=ax, label=outlbl)
                            ECout = (
                                CV[["E_AppV_RHE", "j_mAcm2"]]
                                .rename(
                                    {
                                        "E_AppV_RHE": "EV_%s" % Path(fn).stem,
                                        "j_mAcm2": "j_mA_%s" % Path(fn).stem,
                                    },
                                    axis=1,
                                )
                                .reset_index(drop=True)
                            )
                            EC_out_df.append(ECout)
                    except:
                        pass
                    #                        pd.concat([i for i in EC_out_df],sort=False,axis=1).to_excel(DD_Series.joinpath('EC_%s_%s_%.0f.xlsx'%(El,Sn,int(ScR)))) # export to XLS

                    #                    ax.set_ylim(-1E-4,1E-4)
                    ax.autoscale(True)
                    ax.set_xlim(0, 1.2)
                    ax.legend(
                        bbox_to_anchor=(-0.3, 1.15, 1.8, 0.102),
                        loc="lower left",
                        ncol=1,
                        mode="expand",
                        borderaxespad=0.0,
                    )
                    ax.set_xlabel("$\mathrm{E \//\/ V_{RHE}}$")
                    ax.set_ylabel("$\mathrm{j \//\/ mAcm^{-2}}$")
                    fig.suptitle("CV in N2 at %.2f mV/s" % ScR, y=0.95)
                    plt.savefig(
                        DD_Series.joinpath("EC_%s_%s_%.0f.png" % (El, Sn, int(ScR))),
                        dpi=300,
                        bbox_inches="tight",
                    )
                    plt.close()

    def EC_EIS_ORR_N2():
        a = SampleSelection("CB_paper", "SeriesID", [])
        EC_ovv = a.Retrieve_EC_Table_from_OVVs()
        special = EC_ovv.loc[EC_ovv.SampleLabel.str.contains("-200") == True, :]
        for date, dtgr in special.groupby(by="EXP_date"):
            par_Exps = dtgr.PAR_exp.unique()
            if "EIS" and "ORR" in par_Exps:
                print(par_Exps, dtgr)
        #%%


#        ['BET_SA m2/g_RAW','BET_Area_RPT','Area_in_cell']
##         if 'BET' in yl:
##            for n,r in fltrd.iterrows():
##                r['Isotherm_P/P0'].val
#        selBET = Prep_EA_BET.dropna(subset=['PyRunHash'])
#        selBET.plot(x='SampleLabel',y='BET_Area_RPT',kind='barh')
#        selBET.plot(x='SampleID',y='BET_SA m2/g_RAW',kind='barh')

#        Prep_EA_BET.plot(x='N_content',y='Overall Weight loss, post-Pyr, post-AL (%)',kind='scatter')


def postEC_Reproducibility():
    PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
    postOVVout = LoadPostOVV()
    #    pd.read_excel(PostDestDir.joinpath('postEC_Organized.xlsx'),index_col=[0])
    SampleCodes = pd.read_excel(
        FileHelper.FindExpFolder("VERSASTAT").DestDir.parent.joinpath(
            "SampleCodeLst.xlsx"
        )
    )
    #    SampleCodes  = pd.read_excel(PostDestDir.joinpath('SampleCodeLst.xlsx'))
    postRRDE = postOVVout.loc[postOVVout["Type_Exp"] == "ORR_RRDE"].drop_duplicates()

    PostDestDirRRDE = PostDestDir.joinpath("ORR_RRDE")
    for Elec, ElecGr in postRRDE.groupby(by="Electrolyte"):
        for Sample, sGr in ElecGr.groupby(by="SampleID"):
            splOut = []
            for Fn, Fgr in sGr.groupby(by="SourceFilename"):
                suffix = [
                    Sample,
                    Fgr.Status.values[0],
                    Fgr.Electrolyte.values[0],
                    Fgr.Gas.values[0],
                    Fgr.EXP_date.values[0],
                ]
                SamFn = pd.read_excel(Fn, index_col=[0]).add_suffix(
                    "_%s" % ("_".join(suffix))
                )
                splOut.append([SamFn])
            SampleCols = pd.concat([i[0] for i in splOut], axis=1)
            SampleCols.to_excel(
                PostDestDirRRDE.joinpath("%s_%s_%s.xlsx" % (Sample, Elec, len(splOut)))
            )

    EISpars = pd.concat(
        [
            pd.read_excel(i, index_col=EvRHE)
            for i in [
                i for i in postOVVout.query('Type_Exp == "EIS"').SourceFilename.values
            ]
        ]
    )
    EISpars.columns
    status, extra = postEC_Status(EISpars.File.values)
    EISpars = EISpars.assign(
        **{
            "DATE": pd.to_datetime(
                ["%s-%s-%s" % (i.year, i.month, i.day) for i in EISpars.PAR_date]
            ),
            "Status": status,
            "Status_extra": extra,
        }
    )
    EIS_500mV = EISpars.loc[((EISpars.index < 0.62) & (EISpars.index > 0.58)), :]
    EIS_200mV = EISpars.loc[((EISpars.index < 0.22) & (EISpars.index > 0.18)), :]
    #    EISpars.PAR_date
    #    EISpars = EISpars.set_index(['Electrolyte','Gas','RPM','SampleID'])
    Cdlpars = pd.concat(
        [
            pd.read_excel(i)
            for i in [i for i in list(PostDestDir.rglob("*/*Cdl_FIT_N2*"))]
        ]
    )
    Cdlpars["Gas"] = "N2"
    if "Electrolyte" not in Cdlpars.columns:
        Cdlpars = Cdlpars.assign(
            **{
                "Electrolyte": [
                    FileHelper.FindSampleID.determine_pH_from_filename(Path(i))[
                        "Electrolyte"
                    ][0]
                    for i in Cdlpars.Filename.values
                ],
                "RPM": len(Cdlpars) * [0],
                "Status": postEC_Status(Cdlpars.Filename.values),
            }
        )
    Cdlpars = Cdlpars.set_index(["Electrolyte", "Gas", "RPM", "SampleID"])
    Cdlpars.columns
    ParsMerge = pd.concat([EISpars, Cdlpars], sort=False)
    list(PostDestDir.rglob("*/*ORR_pars*"))
    postOVVout.query('Type_Exp == "EIS"').SourceFilename.values
    ORRpars = pd.concat(
        [
            pd.read_excel(i)
            for i in [
                i
                for i in postOVVout.query(
                    'Type_Exp == "ORR_pars"'
                ).SourceFilename.values
            ]
        ]
    )
    ORRpars["Gas"] = "O2"
    Ostatus, Oextra = postEC_Status(EISpars.File.values)
    ORRpars = ORRpars.assign(**{"Status": Ostatus, "Status_extra": Oextra})
    ORRpars.columns
    #    ORRpars = ORRpars.set_index(['Electrolyte','Gas','RPM','SampleID'])

    ORR_EIS_500 = pd.merge(
        ORRpars,
        EIS_500mV,
        on=["Electrolyte", "Gas", "RPM", "SampleID", "DATE", "Status"],
    )
    ORR_EIS_200 = pd.merge(
        ORRpars,
        EIS_200mV,
        on=["Electrolyte", "Gas", "RPM", "SampleID", "DATE", "Status"],
    )
    refl, cref = [], []
    for a in ORR_EIS_200.SampleID.values:
        ScodeRef = SampleCodes.loc[SampleCodes.SampleID == a, :]
        if ScodeRef.empty:
            Scode = EISovv["SampleID"].unique()[0]
            Ccode = 1
        else:
            Scode = ScodeRef.Sample.values[0]
            Ccode = ScodeRef.Colorcode.values[0]
        refl.append(Scode)
        cref.append(Ccode)
    ORR_EIS_200["SampleLabel"] = refl
    ORR_EIS_200["ColorLabel"] = cref

    for Stat, grSt in ORR_EIS_200.groupby(by=["Status"]):
        for Elec, gr in grSt.groupby(by=["Electrolyte"]):
            for rpm, Rgr in grSt.groupby(by=["RPM"]):
                DestFile = PostDestDirRRDE.joinpath(
                    "PARS_" + "_".join([Stat, Elec, str(rpm), "200mV"]) + ".xlsx"
                )
                Rgr.to_excel(DestFile)

    ParsMerge = pd.concat([ParsMerge, ORRpars], sort=False)

    ParsMerge = ParsMerge.reset_index()
    #%%
    for Stat, grSt in ORR_EIS_200.groupby(by=["Status"]):
        for Elec, gr in grSt.groupby(by=["Electrolyte"]):
            Elec, gr
            if len(gr) > 1:
                fig, ax = plt.subplots(
                    nrows=3, ncols=3, constrained_layout=True, figsize=(8, 8)
                )
                plt.suptitle("%s, %s" % (Elec, Stat))
                ax[0, 0].scatter(gr["Jkin_075"], gr["Rct_kin"])
                ax[0, 1].scatter(gr["E_half"], gr["Rct_kin"])
                ax[1, 0].scatter(gr["E_half"], gr["Qad"], c="red")
                ax[1, 0].set_ylim = (0, 0.01)
                ax[1, 2].scatter(gr["FracH2O2_050"], gr["Qad"], c="red")
                ax[1, 2].set_ylim = (0, 0.01)
                ax[1, 1].scatter(gr["J_diff_lim"], gr["Qad"], c="red")
                ax[1, 1].set_ylim = (0, 0.01)
                ax[2, 2].scatter(gr["E_half"], gr["Cdlp"], c="red")
                ax[2, 0].scatter(gr["TSa_l"], gr["Rct"], c="red")
                ax[2, 1].scatter(gr["TSb_l"], gr["Rct"], c="red")
                ax[0, 2].scatter(gr["FracH2O2_050"], gr["Rct"], c="red")
    #%%
    for Elec, OrrElec in ORR_EIS.groupby(by=["Electrolyte"]):
        [(i, OrrElec.loc[OrrElec.SampleID == i]) for i in OrrElec.SampleID.unique()]
        for sID, grsID in OrrElec.groupby(by=["SampleID"]):
            grsIDstat = grsID.groupby(by=["Status"])
            if len(grsIDstat.groups) > 1:
                fig, ax = plt.subplots(
                    nrows=3, ncols=3, constrained_layout=True, figsize=(8, 8)
                )
                plt.suptitle("%s, %s" % (Elec, sID))
                for Stat, gr in grsIDstat:
                    if Stat == "EoL":
                        color = "blue"
                    elif Stat == "BoL":
                        color = "red"
                    else:
                        color = "k"
                    ax[0, 0].scatter(gr["Jkin_075"], gr["Rct_kin"], c=color)
                    ax[0, 1].scatter(gr["E_half"], gr["Rct_kin"], c=color)
                    ax[1, 0].scatter(gr["E_half"], gr["Qad"], c=color)
                    #                ax[1,0].set_ylim=(0,0.01)
                    ax[1, 2].scatter(gr["FracH2O2_050"], gr["Qad"], c=color)
                    #                ax[1,2].set_ylim=(0,0.01)
                    ax[1, 1].scatter(gr["J_diff_lim"], gr["Qad"], c=color)
                    #                ax[1,1].set_ylim=(0,0.01)
                    ax[2, 2].scatter(gr["E_half"], gr["Cdlp"], c=color)
                    ax[2, 0].scatter(gr["TSa_l"], gr["Rct"], c=color)
                    ax[2, 1].scatter(gr["TSb_l"], gr["Rct"], c=color)
                    ax[0, 2].scatter(gr["FracH2O2_050"], gr["Rct"], c=color)

    #%%
    #        RPM1500 = gr.query('RPM == 1500')
    #    EISpars.join(Cdlpars,on='SampleID')
    #    pd.concat
    [i for i in list(PostDestDir.rglob("*/*Cdl_FIT_N2*"))]
    test2 = postOVVout.loc[postOVVout.SampleID.str.contains("JOS3") == True]
    Pars500mV = pd.DataFrame()
    for nm, gr in test2.groupby(["Electrolyte", "SampleID", "EXP_date"]):
        for nm2, gr2 in gr.groupby(["Status"]):
            print(nm, nm2)
            #            for fn in EISnm2.Filename:
            EISnm2 = gr2.loc[
                (gr2.Type_Exp.str.contains("EIS") == True)
                & (gr2.Filename.str.contains("O2"))
            ]
            if not EISnm2.empty:
                EISnm2pars = pd.read_excel(EISnm2.iloc[0].Filename)
                EISnm2pars = EISnm2pars.loc[
                    EISnm2pars.File.str.contains(
                        Path(EISnm2.iloc[0].Filename).stem.split("_pars2")[0]
                    )
                    == True
                ]
                EISat500mV = EISnm2pars.loc[
                    ((EISnm2pars[EvRHE] < 0.52) & (EISnm2pars[EvRHE] > 0.48)), :
                ].iloc[0]
                if Pars500mV.empty:
                    Pars500mV = EISat500mV
                else:
                    Pars500mV = pd.concat([EISat500mV, Pars500mV])
            ORRnm2 = gr2.loc[(gr2.Type_Exp.str.contains("ORR") == True)]
            if not ORRnm2.empty:
                ORRnm2pars = pd.read_excel(ORRnm2.iloc[0].Filename)
                ORRnm2parsRPM = ORRnm2pars.iloc[ORRnm2pars.RPM.idxmax()]
                if Pars500mV.empty:
                    Pars500mV = ORRnm2parsRPM
                else:
                    Pars500mV = pd.concat([ORRnm2parsRPM, Pars500mV])
            Cdlnm2 = gr2.loc[
                (gr2.Type_Exp.str.contains("N2_Cdl") == True)
                & (gr2.Filename.str.contains("Cdl_FIT"))
            ]
            if not Cdlnm2.empty:
                #                for dt, grCdl in Cdlnm2
                CdlFN = Cdlnm2.iloc[0].Filename
                if Path(CdlFN).suffix == ".csv":
                    Cdlnm2pars = pd.read_csv(Cdlnm2.iloc[0].Filename)
                elif Path(CdlFN).suffix == ".xlsx":
                    Cdlnm2pars = pd.read_excel(Cdlnm2.iloc[0].Filename)
                Cdlat500mV = Cdlnm2pars.loc[
                    ((Cdlnm2pars[EvRHE] < 0.52) & (Cdlnm2pars[EvRHE] > 0.48)), :
                ].iloc[0]
                if Pars500mV.empty:
                    Pars500mV = Cdlat500mV
                else:
                    Pars500mV = pd.concat([Cdlat500mV, Pars500mV])

            #                pd.concat([Cdlat500mV,ORRnm2parsRPM])

            #                Pars500mV = Cdlat500mV.join(Pars500mV,on='SampleID')
            print(nm, nm2)
        #            EISat500mV.join(ORRnm2parsRPM,on='SampleID').join(Cdlat500mV,on='SampleID')

        # a.to_excel(PostDestDir.joinpath('postEC_Organized.xlsx'))
        #%%
        try:
            pH = gr.pH.unique()[0]
            Elec = gr.Electrolyte.unique()[0]
            sID = FindSampleID.match_SampleID(d)
            sDATE = gr.EXP_date.dt.strftime("%Y-%m-%d").unique()[0]
            DestPth = PostDestDir.joinpath(Elec, sID, sDATE)
            DestPth.mkdir(parents=True, exist_ok=True)
            if Go1:
                hits1 = list(Path(d).rglob("*TAFEL*xlsx"))
                #            if len(hits) > 1:
                for t1 in hits1:
                    #                continue
                    xt1 = pd.read_excel(t1)
                    xt1.to_excel(DestPth.joinpath(t1.name))
            #            print(hits)
            if Go2:
                hits2 = list(Path(d).rglob("*O2*ORR*RING*csv"))
                for t2 in hits2:
                    #                print(t2.name)
                    #                xt2 = pd.read_excel(t2)
                    xt2 = pd.read_csv(t2)
                    #                print(xt2.columns)
                    xt2[KeepColls].to_excel(
                        DestPth.joinpath(t2.with_suffix(".xlsx").name)
                    )
            if Go3:
                hits3 = list(Path(d).rglob("*KL_PARs*xlsx"))
                for t3 in hits3:
                    #                print(t2.name)
                    #                xt2 = pd.read_excel(t2)
                    xt3 = pd.read_excel(t3)

                    T3 = xt3.loc[
                        xt3[EvRHE].isin(
                            [
                                i
                                for i in xt3[EvRHE].values
                                if [
                                    a
                                    for a in np.linspace(0, 1, 30)
                                    if np.isclose(a, i, atol=0.020)
                                ]
                            ]
                        )
                        & (xt3["Electrode"] == "I_Disk")
                        & (xt3["Sweep_Type"] == "cathodic"),
                        :,
                    ]
                    #                    for Elnm, Elgr in T3.groupby(by='Electrolyte'):
                    #                    n_Elec = 1/((KLfitting[0])*(0.62*1**(2/3)*nu**(-1/6)*F*C*Area))
                    KLp = KL_coeff.loc[KL_coeff["Electrolyte"] == Elec, :]
                    #                    print(KLp,Elec,pH)
                    T3 = T3.assign(
                        **{
                            "SampleID": sID,
                            "pH": pH,
                            "Electrolyte": Elec,
                            "Date": sDATE,
                            "nElectrons2": 1
                            / (
                                (T3["KL_slope"])
                                * (
                                    0.2
                                    * KLp["C0"].values[0]
                                    * KLp["D0"].values[0] ** (2 / 3)
                                    * KLp["kVis"].values[0] ** (-1 / 6)
                                    * 96485
                                    * WE_SA_collection_eff("PINE")["Disk_cm2"]
                                )
                            ),
                        }
                    )

                    #                T3.plot(x=EvRHE,y='nElectrons',kind='scatter')
                    #                xt3.loc[xt3[EvRHE] == np.linspace(0.1,1,10),:]
                    #                xt3.loc[xt3[EvRHE] == [i for i in xt3[EvRHE].values if [a for a in np.linspace(0,1,10) if np.isclose(a,i,atol=0.010)]],:]
                    OutParsID = pd.concat([OutParsID, T3])
        except Exception as e:
            print("=== Skipped:%s,\n%s ====" % (d, e))
    KLaggOut = pd.DataFrame({EvRHE: []})
    for pHn, pHgr in OutParsID.groupby(by="pH"):
        pHgr.to_excel(PostDestDir.joinpath("KL_Pars_OVV_%s.xlsx" % pHn))
        for KLnm, KLgr in pHgr.groupby(by="SampleID"):
            #            KLagg = KLgr.groupby(by=EvRHE).agg({'nElectrons' : ['mean','count','std']})
            KLagg = (
                KLgr.groupby(by=EvRHE)
                .agg({"nElectrons": ["mean", "std", "count"]})
                .rename(columns={"nElectrons": "%s_%.1f" % (KLnm, pHn)})
                .reset_index()
            )
            KLaggOut = pd.merge(KLagg, KLaggOut, on=EvRHE, how="outer")
    #            KLgr.groupby(EvRHE).agg([np.mean,np.std])
    #            print(pHn,KLnm,KLagg)
    KLaggOutB = KLaggOut.set_index(EvRHE).T
    for pHo in [13, 0.3]:
        pHKLagg = KLaggOutB.loc[
            KLaggOutB.index.get_level_values(0).str.contains(str(pHo)), :
        ]
        pHKLagg.sort_index().to_excel(PostDestDir.joinpath("KL_Pars_AGG_%s.xlsx" % pHo))

    pHKLagg03 = pd.read_excel(
        PostDestDir.joinpath("KL_Pars_AGG_0.3.xlsx"), index_col=[0, 1]
    )
    pHKLagg13 = pd.read_excel(
        PostDestDir.joinpath("KL_Pars_AGG_13.xlsx"), index_col=[0, 1]
    )
    pHKLagg03sl = pHKLagg03.loc[
        :,
        [
            i
            for i in pHKLagg03.columns
            for a in EvRHE_List
            if np.isclose(i, a, atol=0.02)
        ],
    ].unstack(level=-1)
    pHKLagg13sl = pHKLagg13.loc[
        :,
        [
            i
            for i in pHKLagg13.columns
            for a in EvRHE_List
            if np.isclose(i, a, atol=0.02)
        ],
    ].unstack(level=-1)

    sA, sB, CB = (
        ["DW16", "DW17", "DW18", "DW19", "DW20", "DW21"],
        ["DW24", "DW25", "DW26", "DW27", "DW28", "DW29"],
        ["DW38A|DW38B|DW38C|DW38D|DW38E|DW38F"],
    )
    sABstr = "DW16|DW17|DW18|DW19|DW20|DW21|DW24|DW25|DW26|DW27|DW28|DW29"

    pHKLagg03_Paper2018 = pHKLagg03sl.loc[
        pHKLagg03sl.index.get_level_values(0).str.contains(sABstr), :
    ]
    pHKLagg13_Paper2018 = pHKLagg13sl.loc[
        pHKLagg13sl.index.get_level_values(0).str.contains(sABstr), :
    ]
    pHKLagg03_Paper2018.to_excel(PostDestDir.joinpath("KL_Pars_AGG_Origin_0.3.xlsx"))
    pHKLagg13_Paper2018.to_excel(PostDestDir.joinpath("KL_Pars_AGG_Origin_13.xlsx"))

    pHKLagg03_Paper2018_mean = pHKLagg03_Paper2018.loc[
        :, [i for i in pHKLagg03_Paper2018.columns if "mean" in i]
    ]
    pHKLagg13_Paper2018_mean = pHKLagg13_Paper2018.loc[
        :, [i for i in pHKLagg13_Paper2018.columns if "mean" in i]
    ]
    #%%
    # CHEKC ALKALINE PARAMETERS!!!! #
    for i in pHKLagg03_Paper2018_mean.columns:
        fig, ax = plt.subplots(1, 1)

        ax.scatter(
            pHKLagg03_Paper2018_mean.index.values,
            pHKLagg03_Paper2018_mean[i].values,
            s=100,
            c="orange",
        )
        ax.scatter(
            pHKLagg03_Paper2018_mean.index.values,
            pHKLagg13_Paper2018_mean[i].values,
            s=100,
            c="blue",
        )
        ax.set_ylim(0, 6)
        fig.suptitle("%.2f V_RHE" % i[0])
        plt.xticks(
            pHKLagg03_Paper2018_mean.index.values,
            pHKLagg03_Paper2018_mean.index.values,
            rotation=-60,
        )
        plt.savefig((PostDestDir.joinpath("KL_Pars_AGG_%.2f V_RHE.png" % i[0])))
    #        Ekl =
    return print("Post EC succes: %s" % PostDestDir)


#        fig,ax = plt.subplots()
#        ax.barh(pHKLagg.index,pHKLagg.iloc[:,Ekl])
#        pHKLagg.plot(x=pHKLagg.index,y=pHKLagg.iloc[:,5],kind='barh')
#    fig = plt.figure()
#    ax =fig.add_subplot(111, projection='3d')
#    Xa = pHKLagg13_Paper2018_mean.unstack().index.get_level_values(2).values
#    X = [int(x.split('DW')[1].split('_')[0]) for x in Xa]
#    Y = pHKLagg13_Paper2018_mean.unstack().index.get_level_values(0)
#    Z = pHKLagg13_Paper2018_mean.unstack().values
#
#    ax.set_zlim(0, 6)
#    for angle in range(0, 360):
#
#        ax.view_init(30, angle)
#        surf = ax.scatter(X, Y, Z, antialiased=True)
#        plt.draw()
#        plt.pause(.001)
#        plt.show()
#        plt.close()
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#%%
class PrepareHypotheses:
    def __init__():
        pass

    def Cdl_CV_origin():
        xls_Cdl = Path(
            "G:\\CloudStation\\Presentations\\2019-06-14_Fe-SiO2_ORR+EIS\\data\\N2\\0.1MH2SO4_JOS4_standard_03-25\\Cdl_data_N2_20cls_300_100_10_JOS4_283_v20.xlsx"
        )
        cdl = pd.read_excel(xls_Cdl, index_col=[0])
        cdl.E_V = cdl.E_V.apply(lambda x: np.round(x, 3))
        [(n, gr) for n, gr in cdl.groupby("ScanRate")]
        cdl_cols = pd.concat(
            [
                i[1]
                for i in [
                    (
                        n,
                        gr.rename(
                            columns=dict(
                                zip(
                                    gr.columns,
                                    [i + "_{0}".format(n) for i in gr.columns],
                                )
                            )
                        ).set_index("E_V_{0}".format(n)),
                    )
                    for n, gr in cdl.groupby("ScanRate")
                ]
            ],
            axis=1,
        )
        cdl_cols.iloc[0::11].to_excel(xls_Cdl.parent.joinpath("Cdl_fits_origing.xlsx"))

        cath = (
            cdl.loc[cdl.E_V.isin(cdl.dropna(subset=["j_cathodic"]).E_V.unique()[0::7])]
            .dropna(subset=["j_cathodic"])
            .drop(columns="j_anodic")
            .reset_index()
        )
        anod = (
            cdl.loc[cdl.E_V.isin(cdl.dropna(subset=["j_anodic"]).E_V.unique()[0::7])]
            .dropna(subset=["j_anodic"])
            .drop(columns="j_cathodic")
            .reset_index()
        )
        cath.drop_duplicates(["ScanRate", "E_V", "j_cathodic"])
        cath_cdl = cath.drop_duplicates(["ScanRate", "E_V"]).pivot(
            index="ScanRate", columns="E_V", values="j_cathodic"
        )
        anod_cdl = anod.drop_duplicates(["ScanRate", "E_V"]).pivot(
            index="ScanRate", columns="E_V", values="j_anodic"
        )

        linfits_results_cath = [
            [
                i,
                linregress(cath_cdl.index, cath_cdl[i].values),
                linregress(np.sqrt(cath_cdl.index), cath_cdl[i].values),
            ]
            for i in cath_cdl.columns
        ]
        linfits_results_anod = [
            [
                i,
                linregress(anod_cdl.index, anod_cdl[i].values),
                linregress(np.sqrt(anod_cdl.index), anod_cdl[i].values),
            ]
            for i in anod_cdl.columns
        ]

        linfits_cath = pd.concat(
            [
                i
                for i in [
                    pd.DataFrame(i, index=cath_cdl.index)
                    for i in [
                        {i[0]: cath_cdl.index.values * i[1][0] + i[1][1]}
                        for i in linfits_results_cath
                    ]
                ]
            ],
            axis=1,
        )
        linfits_anod = pd.concat(
            [
                i
                for i in [
                    pd.DataFrame(i, index=anod_cdl.index)
                    for i in [
                        {i[0]: anod_cdl.index.values * i[1][0] + i[1][1]}
                        for i in linfits_results_anod
                    ]
                ]
            ],
            axis=1,
        )
        fig, ax = plt.subplots()
        linfits_cath.plot(ax=ax)
        for i in cath_cdl.columns:
            ax.scatter(cath_cdl.index.values, cath_cdl[i])

        cath_cdl.to_excel(xls_Cdl.parent.joinpath("Cdl_data_cath.xlsx"))
        anod_cdl.to_excel(xls_Cdl.parent.joinpath("Cdl_data_anod.xlsx"))

        linfits_cath.to_excel(xls_Cdl.parent.joinpath("Cdl_cath_fits.xlsx"))
        linfits_anod.to_excel(xls_Cdl.parent.joinpath("Cdl_anod_fits.xlsx"))

    def EIS_example_fit_data():
        DestDirTop = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath(
            "EIS_example"
        )
        DestDirTop.mkdir(parents=True, exist_ok=True)
        #    default_pars = {'Rs' : 60, 'Cdlp' : 2E-04,'nDL' : 1, 'Rct' : 100,'Qad' : 1E-03,'nAd' : 0.5,'Rorr' : 7E5}
        default_pars = {
            "Rs": 60,
            "Cdlp": 2e-04,
            "nDL": 0.5,
            "Rct": 0.1,
            "Qad": 1e-03,
            "nAd": 0.99,
            "Rorr": 7e5,
        }
        ranges = {
            "Rct": [10, 20, 50, 100, 500, 1e3, 5e3],
            "Rs": [10, 20, 40, 80, 200],
            "nDL": [0.5, 0.7, 0.8, 0.9, 1],
            "nAd": [0.5, 0.7, 0.8, 0.9, 1],
            "Qad": [1e-05, 5e-05, 1e-04, 5e-04, 1e-03],
            "Cdlp": [1e-06, 5e-06, 1e-05, 5e-05],
            "Rorr": [0, 1e0, 1e1, 1e2, 1e3, 1e4],
        }

        DestTopDir = DestDirTop.joinpath("test_example_hypotheses")
        DestTopDir.mkdir(parents=True, exist_ok=True)

        EIS_pars.query('SampleID == "JOS4" & pH == 1')
        jos4 = EIS_pars.query(
            'SampleID == "JOS4" & pH == 1 & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        )
        #        jos4_20190506 = jos4.query('(PAR_date_day == "2019-03-25")')
        jos4_20190325 = jos4.query('(PAR_date_day == "2019-03-25")')
        JOS4_n2 = jos4_20190325.query('Gas == "N2"').loc[
            jos4_20190325.basename.str.contains("post") == False
        ]
        JOS4_O2 = jos4_20190325.query('Gas == "O2"').loc[
            jos4_20190325.basename.str.contains("post") == False
        ]
        Elst = [0.9, 0.65, 0.5, 0.1]
        for par in eisplot.parlst:
            jos4_20190325
            fig, ax = plt.subplots()

            JOS4_n2.plot(x="E_RHE", y=par, kind="scatter", ax=ax, c="b")
            JOS4_O2.plot(x="E_RHE", y=par, kind="scatter", ax=ax, c="r")
            ax.set_ylim(eisplot(par).pars_plot_lims()[0])
        jos4_20190325.loc[jos4_20190325.E_RHE.isin(Elst)]

        Y, Z = eisplot.AmdImp()["Y"], eisplot.AmdImp()["Z"]
        fit_pars_lst = []
        jos4_data = pd.concat(
            [
                ExportECfromCV.make_uniform_EvRHE(pd.read_excel(n)).assign(
                    **{"SpectraFile": n, "Gas": gr.Gas.unique()[0]}
                )
                for n, gr in jos4_20190325.groupby("SpectraFile")
            ]
        )
        #        jos4_data = jos4_data.assign(**{ 'DATA_Z_complex' : [np.complex(i)  for i in jos4_data.DATA_Z.values]})
        jos4_data = jos4_data.assign(
            **{
                "DATA_Zabs": np.abs(
                    [
                        np.complex(a, b)
                        for a, b in zip(
                            jos4_data.DATA_Zre.values, jos4_data.DATA_Zim.values
                        )
                    ]
                ),
                "DATA_Zangle": np.angle(
                    [
                        np.complex(a, b)
                        for a, b in zip(
                            jos4_data.DATA_Zre.values, jos4_data.DATA_Zim.values
                        )
                    ],
                    deg=True,
                ),
            }
        )

        ranges_jos4 = {
            "Rct": [1, 25, 50, 100, 500],
            "Rs": [5, 20, 40, 80],
            "nDL": [0.5, 0.7, 0.8, 0.9, 1],
            "nAd": [0.5, 0.7, 0.8, 0.9, 1],
            "Qad": [20e-04, 32e-04, 44e-04, 56e-04, 68e-04],
            "Cdlp": [4e-4, 14e-04, 24e-04, 34e-04],
            "Rorr": [0, 1e2, 1e3, 2e3, 1e4],
        }

        for Ev in Elst:
            fig, ax = plt.subplots()
            DestDir = DestTopDir.joinpath(str(Ev))
            DestDir.mkdir(parents=True, exist_ok=True)
            spectra = jos4_data.query("E_RHE == @Ev")
            for gas, gasgrp in jos4_data.groupby("SpectraFile"):
                #                spf = gasgrp.SpectraFile.unique()[0]
                #                spdf = ExportECfromCV.make_uniform_EvRHE(pd.read_excel(spf))
                fit_pars_item = jos4_20190325.query(
                    "E_RHE == @Ev & SpectraFile == @gas"
                ).drop_duplicates(subset=eisplot.parlst)[
                    eisplot.parlst + ["E_RHE", "Gas", "SpectraFile"]
                ]
                fit_pars_lst.append(fit_pars_item)

                gasgrp.query("E_RHE == @Ev").plot(
                    x=Z["x_data"],
                    y=Z["y_data"],
                    ax=ax,
                    kind="scatter",
                    title="E / V v RHE = {0}".format(Ev),
                    label="{0}".format(gasgrp.Gas.unique()[0]),
                )
                gasgrp.query("E_RHE == @Ev").iloc[::].plot(
                    x=Z["x_fit"], y=Z["y_fit"], ax=ax, kind="line"
                )
                jos4_20190325.query("E_RHE == @Ev & SpectraFile == @gas")
                jos4_20190325.query(
                    "E_RHE == @Ev & SpectraFile == @gas"
                ).drop_duplicates(subset=eisplot.parlst)
                fit_pars_item = jos4_20190325.query(
                    'E_RHE == @Ev & Gas == "O2"'
                ).drop_duplicates(subset=eisplot.parlst)[
                    eisplot.parlst + ["E_RHE", "Gas", "SpectraFile"]
                ]
                fit_pars_lst.append(fit_pars_item)
                #                spectra = gasgrp.query('E_RHE == @Ev')
                PrepareHypotheses.jos4_examples(
                    spectra, fit_pars_item, ranges_jos4, DestDir
                )

            plt.show()
            plt.close()
        fit_pars = pd.concat(fit_pars_lst)
        default_pars = {
            "Rs": 60,
            "Cdlp": 2e-04,
            "nDL": 0.5,
            "Rct": 0.1,
            "Qad": 1e-03,
            "nAd": 0.99,
            "Rorr": 7e5,
        }

    def jos4_examples(spectra, fit_pars_item, ranges, Destdir):
        default_pars = fit_pars_item[eisplot.fitparlst].squeeze().to_dict()

        Y, Z = eisplot.AmdImp()["Y"], eisplot.AmdImp()["Z"]

        for par in ranges_jos4.keys():
            outpar, outdata = [], []

            color = plt.cm.viridis(np.linspace(0, 1, len(ranges_jos4[par])))
            matplotlib.rcParams["axes.prop_cycle"] = cycler.cycler("color", color)
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12))
            ax, ax2 = axes[0, 0], axes[0, 1]
            ax3, ax4 = axes[1, 0], axes[1, 1]
            EC_exp_nuniq = [
                i
                for i in [
                    i
                    for i in spectra.columns
                    if i in SampleSelection.EC_exp_cols + ["E_RHE"]
                ]
                if spectra[i].nunique() == 1
            ]
            EC_exp_cond1 = [
                "{0} : {1}".format(i, spectra[i].unique()[0]) for i in EC_exp_nuniq
            ]
            par_title = ", ".join(
                "{:s} = {:2.2g}".format(k, v)
                for k, v in default_pars.items()
                if k not in par
            )
            EC_title = ", ".join(EC_exp_cond1) + "\n" + par_title
            plt.suptitle(EC_title)
            #        matplotlib.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
            #        hexcolor = map(lambda rgb:'#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255),
            #               tuple(color[:,0:-1]))
            for spn, spgr in spectra.groupby(["SpectraFile", "Gas"]):
                color = "darkred" if "O2" in spn else "blue"
                ax.scatter(
                    spgr[Z["x_data"]].values,
                    spgr[Z["y_data"]].values,
                    c=color,
                    s=60,
                    alpha=0.5,
                )
                ax2.scatter(
                    spgr[Y["x_data"]].values,
                    spgr[Y["y_data"]].values,
                    c=color,
                    s=60,
                    alpha=0.5,
                    label=spn[1],
                )
                ax3.scatter(
                    spgr["Frequency(Hz)"].values,
                    spgr["DATA_Zabs"].values,
                    c=color,
                    s=60,
                    alpha=0.5,
                )
                ax4.scatter(
                    spgr["Frequency(Hz)"].values,
                    spgr["DATA_Zangle"].values,
                    c=color,
                    s=60,
                    alpha=0.5,
                )

            for val, cl in zip(ranges_jos4[par], color):
                dfile = DestDir.joinpath("{0}_{1}_data.xlsx".format(par, val))
                #            print(par,val)
                od, oprs = PAR_EIS_fit.make_example_plot(
                    par=par, val=val, default=default_pars, single_plot=0
                )
                od["color"] = val
                od.to_excel(dfile)
                #            matplotlib.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
                outpar.append(oprs)
                outdata.append(od)
                #            od.plot(x='INIT2_Zre',y='INIT2_-Zim',kind='scatter',s=80, c=[cl]*len(od),ax=ax,label='{0} = {1}'.format(par,val))
                ax.plot(
                    od["INIT2_Zre"].values,
                    od["INIT2_-Zim"].values,
                    lw=3,
                    label="{:s} = {:2.2g}".format(par, val),
                )

                ax2.plot(od["INIT2_Yre"].values, od["INIT2_Yim"].values, lw=3)
                #                ax2.scatter(spectra['INIT2_Yre'].values, spectra['INIT2_Yim'].values)
                ax3.plot(
                    od["Frequency(Hz)"].values,
                    od["INIT2_Zabs"].values,
                    lw=3,
                    label="{:s} = {:2.2g}".format(par, val),
                )

                ax4.plot(
                    od["Frequency(Hz)"].values,
                    od["INIT2_Zangle"].values,
                    lw=3,
                    label="{:s} = {:2.2g}".format(par, val),
                )
            #            od.plot(x='INIT2_Yre',y='INIT2_Yim',kind='scatter',s=80, c=[cl]*len(od),ax=ax2,alpha=0.5)

            all_oprs = pd.concat(
                [pd.DataFrame(i, index=[0]) for i in outpar],
                ignore_index=True,
                sort=False,
            )
            all_data = pd.concat([i for i in outdata], ignore_index=True, sort=False)
            """ax2: Admittance Plots"""
            ax.set_xlim(0, all_data["INIT2_Zre"].max())
            ax.set_xlabel("Z real")
            ax.set_ylim(0, all_data["INIT2_Zre"].max())
            ax.set_ylabel("Z imag")
            ax.legend(loc="upper left")
            ax.grid(True), ax2.grid(True)
            ax2.set_xlim(0, all_data["INIT2_Yre"].max())
            ax2.set_xlabel("Y real")
            ax2.set_ylabel("Y imag")
            ax2.set_ylim(0, all_data["INIT2_Yre"].max())
            ax2.legend(loc="upper right")
            ax3.set_xscale("log")
            ax3.set_yscale("log")
            ax3.set_ylabel("log(|Z| / Ohm)")
            ax3.set_xlabel("log(Frequency / Hz)")

            ax4.set_xscale("log")
            ax4.set_ylabel("phase agle / degree")
            ax4.set_xlabel("log(Frequency / Hz)")
            ax4.grid(True), ax3.grid(True)
            plt.savefig(DestDir.joinpath("{0}.png".format(par)), dpi=300)
            plt.show()
            print("Figure saved to: {0}".format(DestDir))
            plt.close()
            #        all_data.plot(x='INIT2_Zre',y='INIT2_-Zim',kind='scatter',c='color',lw=2.5,label='{0} = {1}'.format(par,val))
            #        all_data.plot(x='INIT2_Yre',y='INIT2_Yim',kind='scatter',c='color',lw=2.5,label='{0} = {1}'.format(par,val))
            all_oprs.to_excel(DestDir.joinpath("{0}.xlsx".format(par)))


class ReadPostEC:
    """reading files only"""

    def __init__():
        pass

    def EIS_read():
        EIS_file = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "EIS_pars_IndexOVV_v20.xlsx"
        )
        ORR_file = FileHelper.FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_ORR_v20.xlsx"
        )
        if EIS_file.exists():
            EIS_pars = pd.read_excel(EIS_file, index_col=[0])
        if ORR_file.exists():
            ORR_pars = pd.read_excel(ORR_file, index_col=[0])

        CB = EIS_pars.query(
            '(SeriesID == "CB3") | (SeriesID == "CB4") | (SeriesID == "CB5") | (SeriesID == "CB6")'
        )

        EIS_KOH_O2 = EIS_pars.loc[
            (EIS_pars.Electrolyte.str.contains("KOH")) & (EIS_pars.Gas == "O2")
        ]
        EIS_KOH_O2.SeriesID.unique()
        EIS_KOH_O2_CB = EIS_KOH_O2.query('SeriesID = "CB3"')
        loc[(EIS_KOH_O2.SeriesID.str.contains("CB"))]
        ORR_KOH_O2 = ORR_pars.loc[
            (ORR_pars.Electrolyte.str.contains("KOH")) & (ORR_pars.Gas == "O2")
        ]

        fig, ax = plt.subplots()
        for nm, gr in EIS_pars.groupby(by="SampleID"):
            gr.plot(x="E_RHE", y="Rs", kind="scatter", c="pH", ax=ax, label=nm)

        CB.query("pH < 7").plot(
            x="E_RHE",
            y="Cdlp",
            kind="scatter",
            c="Loading_cm2",
            colormap="viridis",
            ylim=(0, 10e-03),
        )
        ORR_KOH_O2.plot(
            x="BET_cat_agg",
            y="Jkin_075",
            kind="scatter",
            c="ML",
            colormap="viridis",
            label=nm,
            logy=1,
        )

    # from FolderOrganizing import FileHelper


# from FolderOrganizing import RunEC_classifier
# from FolderOrganizing import PAR_EIS_fit_V2
# import run_PAR_DW
#### INITIALIZE EASY ####
#        pass
#        try:
#            SampleSelect_all = SampleSelection('*','*')
#            SampleCodesChar = SampleSelect_all.Prep_EA_BET
#        except Exception as e:
#            print('ERRROR',e)
#            SampleCodesChar = pd.DataFrame()
#        pass
#    PostEC.StartLogging()
