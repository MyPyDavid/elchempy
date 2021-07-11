#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 18:48:27 2020

@author: zmg
"""


import datetime as dt
import re
from pathlib import Path

# import hashlib

import pandas as pd
import numpy as np

# import sys
import subprocess
import glob

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations

# from file_py_helper.FindSampleID import GetSampleID
#  form .. import EXP_FOLDER # FIXME
EXP_FOLDER = FindExpFolder("VERSASTAT").DataDir

from .RunEC_classifier import EC_classifier_multi_core

import logging

logger = logging.getLogger(__name__)
#%%


def RunEC(ExpDirs, EXP_FOLDER, mode=""):
    for month, explst in ExpDirs:
        alldirs = explst[::]
        print(alldirs)
        # === MODE DEPENDENT; Load previous made files ===
        if mode == "append" or mode == 11:
            print("#====== Running with previously saved data =====#")
            dest_plot_folder = EXP_FOLDER.DestDir
            try:
                saved = pd.read_csv(dest_plot_folder.joinpath("STATUS_OVV.csv"))
                status = pd.merge(
                    alldirs, saved, on="dir", how="outer", indicator="indicator_column"
                )
                CALC_PARS_svd = pd.read_csv(
                    glob.glob(dest_plot_folder.joinpath("CALC_PARS_ORR_OVERVIEW*"))[-1]
                )
                ORR_JKIN_svd = pd.read_csv(
                    glob.glob(dest_plot_folder.joinpath("ORR_JKIN_calc_OVV_*"))[-1]
                )
                Cdl_fit_svd = pd.read_csv(
                    glob.glob(dest_plot_folder.joinpath("N2_Cdl_calc_OVV_*"))[-1]
                )
                HPRR_svd = pd.read_csv(
                    glob.glob(dest_plot_folder.joinpath("HPRR_PARS_calc_OVV_*.csv"))[-1]
                )
            except Exception as e:
                print(e)
        elif mode == "new" or mode == 0:
            print("#====== Starting from zero =====#")
            saved = pd.DataFrame(
                [[i, 0, 0, 0, dt.datetime.now()] for i in alldirs],
                columns=["dir", "N2", "ORR", "HPRR", "Last_run"],
            )
            status = pd.merge(alldirs, saved, on="dir", indicator="indicator_column")

            CALC_PARS_svd, ORR_JKIN_svd, Cdl_fit_svd, HPRR_svd = (
                pd.DataFrame([]),
                pd.DataFrame([]),
                pd.DataFrame([]),
                pd.DataFrame([]),
            )


# === ---------------- ====


def EC_loadOVV(ForceRecalculation="no", exp_folder_search_set=None):
    # === OVV filter ====
    ExpParOVVraw = pd.DataFrame()
    index_fls = EXP_FOLDER.DestDir.glob("EXP_PAR_class_20*.xlsx")
    ExpLsOVV = sorted([(i, i.stat().st_ctime) for i in index_fls], key=lambda x: x[1])[
        -1
    ][0]
    try:
        #        sorted(FindExpFolder('VERSASTAT').DestDir.glob('EXP_PAR_class_20*.xlsx'))[-1]
        ExpParOVVraw = pd.read_excel(ExpLsOVV, index_col=[0])
    except Exception as e:
        print(e)

    #    logger.info('EC_loadOVV used: {0}'.format(ExpLsOVV))
    #    or np.abs(len(ExpParOVVraw)-len(FileHelper.FindExpFolder('VERSASTAT').EC_find_folders()[1])
    td = dt.date.today()
    out_path = FindExpFolder("VERSASTAT").DestDir.joinpath(
        f"EXP_PAR_class_{td.year}-{td.month}-{td.day}.xlsx"
    )
    if (
        pd.datetime.now() - pd.to_datetime(ExpLsOVV.stat().st_ctime, unit="s")
    ).days > 20 or ForceRecalculation != "no":
        if re.search("index|reindex", ForceRecalculation, re.IGNORECASE):
            print("Replace and Re-index PAR files from existing")
            logger.info(
                "Replace and Re-index PAR files, forced? {0}. Start".format(
                    ForceRecalculation
                )
            )

        #            EC_classifier_multi_core.main(exp_folder_search= exp_folder_search_set, pre_load_RunOVV = ExpParOVVraw)
        elif "force" in ForceRecalculation:
            print("Comple Replace and Re-index of all PAR files")
            logger.info(
                "Replace and Re-index PAR files, forced? {0}. Start".format(
                    ForceRecalculation
                )
            )
            EC_class = EC_classifier_multi_core(
                exp_folder_search=exp_folder_search_set, main_run=True
            )
        #            EC_class.main_run()
        #            EC_classifier_multi_core.main(exp_folder_search= exp_folder_search_set)
        else:
            logger.warning(
                "Replace and Re-index PAR files, forced? {0}. Finish".format(
                    ForceRecalculation
                )
            )
        index_fls = FindExpFolder("VERSASTAT").DestDir.glob("EXP_PAR_class_20*.xlsx")
        ExpLsOVV = sorted(
            [(i, i.stat().st_ctime) for i in index_fls], key=lambda x: x[1]
        )[-1][0]
        #        ExpLsOVV = sorted(FindExpFolder('VERSASTAT').DestDir.glob('EXP_PAR_class_20*.xlsx'))[-1]
        ExpParOVVraw = pd.read_excel(ExpLsOVV, index_col=[0])
    else:
        logger.info("Keep present file:%s" % str(ExpLsOVV))

    ExpParOVV = ExpParOVVraw.loc[
        ~ExpParOVVraw.EXP_dir.str.contains(
            "RRDE_Protocols|data_extract|old|settings|_xx_"
        ),
        :,
    ]
    #    dt =[]
    #    for i,a in ExpParOVV.iterrows():
    #        try:
    #            d = pd.to_datetime(np.datetime64(a.EXP_date))
    #        except:
    #            d = GetSampleID.Determine_date_from_filename(str(a.EXP_dir))
    #        dt.append(d)
    #    ExpParOVV.loc[:,'EXP_date'] = dt
    ExpParOVV = ExpParOVV.loc[ExpParOVV.EXP_dir.str.contains("old") == False, :]
    ExpParOVV = ExpParOVV.assign(
        **{"Date_PAR_EXP": ExpParOVV.PAR_date - ExpParOVV.EXP_date}
    )
    ExpParOVV = FileOperations.ChangeRoot_DF(
        ExpParOVV, ["Dest_dir", "EXP_dir", "PAR_file"]
    )
    ExpParOVV_target_base = FindExpFolder("VERSASTAT").DestDir.joinpath("RunOVV.xlsx")
    ExpParOVV_target = FileOperations.CompareHashDFexport(
        ExpParOVV, ExpParOVV_target_base
    )
    logger.info(
        "EC_loadOVV used: {0}\n to make: {1}".format(ExpLsOVV, ExpParOVV_target)
    )
    return ExpParOVV


class MainEC:
    """This is the main class for running the EC analysis script in a big loop over the experimental PAR files.
    First it indexes the files and creates an Overview in a DataFrame, or takes an existing index (OVV, EC_index).
    Then it loops over a groupby of the Experimental Folders (Exp_Dir). It also starts the logging for debugging and info on the script.
    """

    #%%  MainEC
    def __init__(self):
        pass

    #        log_fn = FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('EC_logger.log')
    #        logging.basicConfig(filename=log_fn, filemode='w', level=level_log,format='%(asctime)s %(levelname)s, %(lineno)d: %(message)s')
    def MainPrepareList(ForceRecalc="no"):
        """Loads the Overview (OVV) depending on if it exist.
        BYPASSED by using EC_LoadOVV directly..."""
        mode, OneTestDir = "new", False
        #        [(i,a) for i,a in enumerate(pathlib.Path.cwd().parts) if a == 'CloudStation']
        #        ExpPAR_DB =  pd.HDFStore(PathDB)
        ExpParOVV = EC_loadOVV(ForceRecalculation=ForceRecalc)
        EC_index = ExpParOVV
        #        .query('(EXP_date > 20181001)')
        return EC_index

    #        ExpPARs = pd.read_hdf(FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_Files_class.hdf5'))
    #        ExpParGrp = aa.groupby(by='EXP_dir')
    #        for exp_dir,gr in EC_index.loc[EC_index.EXP_dir.str.contains('JOS1')].groupby(by='EXP_dir'):
    #%%
    #        test1 = EC_index.query('(EXP_date > 20190214) & (EXP_date < 20190301)')
    #       EC_index.loc[EC_index.EXP_dir.str.contains('25.01.2019_0.1MH2SO4_cell2*') == True]
    #        for exp_dir,gr in EC_index.groupby(by='EXP_dir'):
    def EC_Analysis_Input(TakeRecentList=True, exp_folder=None, force_recalc="index"):
        """Loads the acutal Overview (OVV) depending on if it exists or not and if it need to be rec.
        BYPASSED by using EC_LoadOVV directly ..."""
        version = FileOperations.version
        _RunOVV_fn_format = f"RunOVV_v{version}.pkl"
        RunOVV_fn_lst = list(FindExpFolder("VERSASTAT").DestDir.glob(_RunOVV_fn_format))
        if not RunOVV_fn_lst:
            EC_index = EC_loadOVV(ForceRecalculation=force_recalc)
            RunOVV_fn = FindExpFolder("VERSASTAT").DestDir.joinpath(_RunOVV_fn_format)
        #            RunOVV_fn = list(FileHelper.FindExpFolder('VERSASTAT').DestDir.rglob('RunOVV*xlsx'))[0]
        else:
            RunOVV_fn = [i for i in RunOVV_fn_lst if not "_Conflict" in i.stem][0]
            try:
                OvvFromFile = FileOperations.PDreadXLorCSV(RunOVV_fn_lst[0])
            #                pd.read_excel(RunOVV_fn,index_col=[0])
            except Exception as e:
                OvvFromFile = pd.DataFrame()
                logger.error("Index, error reading RunOVV: {e}\nFile: {RunOVV_fn}")
            if RunOVV_fn.exists() and not OvvFromFile.empty:
                logger.info("Index, took recent file: {0}".format(RunOVV_fn))
                OvvFromFile["Date_PAR_EXP"] = (
                    OvvFromFile.PAR_date - OvvFromFile.EXP_date
                )
                if TakeRecentList == True:
                    EC_index = OvvFromFile
                else:
                    EC_index = EC_loadOVV(
                        ForceRecalculation=force_recalc,
                        exp_folder_search_set=exp_folder,
                    )
            else:
                #            hashTarget = pd.util.hash_pandas_object(EC_index,index=False).sum()
                EC_index = EC_loadOVV(
                    ForceRecalculation="force", exp_folder_search_set=exp_folder
                )
        #        EC_index.Date_PAR_EXP = pd.to_timedelta(EC_index.Date_PAR_EXP,unit='ns')

        _EC_date_min = [
            (
                f"{pd.to_datetime(i):%Y-%m-%d}",
                f"{pd.to_datetime(i):%Y-%m-%dT%H:%M}",
                dt.date.fromisoformat(np.datetime_as_string(np.datetime64(i, "D"))),
            )
            for i in EC_index.PAR_date.to_numpy()
        ]

        EC_index = EC_index.assign(
            **{
                "PAR_date_day": [i[0] for i in _EC_date_min],
                "PAR_date_min": [i[1] for i in _EC_date_min],
                "PAR_date_day_dt": [i[2] for i in _EC_date_min],
            }
        )
        EC_index["Loading_cm2"] = EC_index["Loading_cm2"].round(3)

        EC_index = add_detect_duplicate_measurements(EC_index)
        #            EC_index = MainEC.MainPrepareList()
        #        hashSource = pd.util.hash_pandas_object(OvvFromFile,index=False).sum()
        #        test1 = EC_index
        #        test2 = EC_index.query('(EXP_date == 20190326) | (EXP_date == 20190212)')
        #        test3 = EC_index.query('(EXP_date == 20190301) & (EXP_date < 20190401)')
        #        test3 = EC_index.query('(EXP_date >= 20190301)  & ((SampleID >= "DW01") | (SampleID == "JOS12") | (SampleID == "JOS13") | (SampleID == "JOS14") | (SampleID == "JOS15")) & (PAR_exp == "EIS")')
        #        test3.SampleID.values
        #        filepath_cols = ['Dest_dir','EXP_dir','PAR_file']
        #        for i in filepath_cols:
        #            if  Path(EC_index[i].iloc[0]).is_file() == False:
        #                CS_parts_PDD = FileHelper.FileOperations.find_CS_parts(FileHelper.FindExpFolder('BET').DataDir)
        #                chLstr = [CS_parts_PDD[0].joinpath(FileHelper.FileOperations.find_CS_parts(str(i))[1]) for i in  EC_index[i].values]
        #                EC_index[i] = chLstr
        #                logger.info('Index file paths readjusted for local drive:{0}'.format(i))
        EC_index = EC_index.query('SampleID != "TEMPLATE"')
        EC_index.PAR_file = [Path(i) for i in EC_index.PAR_file.to_numpy()]
        EC_index = FileOperations.ChangeRoot_DF(
            EC_index, ["Dest_dir", "EXP_dir", "PAR_file"]
        )
        FindExpFolder("VERSASTAT").IndexDir.mkdir(parents=True, exist_ok=True)
        new_RunOVV_target = FileOperations.CompareHashDFexport(EC_index, RunOVV_fn)
        logger.info("Index, took recent file and saved: {0}".format(new_RunOVV_target))
        return EC_index


#        return test2
def add_detect_duplicate_measurements(EC_index):
    if not ("dupli_name" and "dupli_num") in EC_index.columns:
        _compare_cols = [
            "PAR_date_day",
            "PAR_exp",
            "SampleID",
            "pH",
            "Gas",
            "Loading_cm2",
            "postAST",
        ]
        _compgrps = EC_index.groupby(_compare_cols)

        _col = []
        for ngr, gr in _compgrps:
            _name = "_".join([str(i) for i in ngr])
            for n, r in enumerate(gr.itertuples()):
                _col.append({"index": r.Index, "dupli_name": _name, "dupli_num": n})

        _dupliDF = pd.DataFrame(_col).set_index("index").sort_index()
        return pd.concat([EC_index, _dupliDF], axis=1)
    else:
        return EC_index

    #        (EXP_date > 20180101)
    #        test2 = EC_index.loc[EC_index.EXP_dir.str.contains('26.10.2018_PTA200_bipot_0.1MH2SO4_RRDE22775')]
    #        EC_index.query('(EXP_date == 26.10.2018_PTA200_bipot_0.1MH2SO4_RRDE22775)')
    #                print('Fail: %s\n because %s' %(exp_dir,e))
    #            num_dir = exp_dir.Index
    #%%
    #        return print('Main EC loop Success!')

    def MergeIndexOvv(results, ovv):
        IndexDir = FindExpFolder("VERSASTAT").IndexDir
        index = pd.read_excel(results, index_col=[0])
        ovv_index = pd.merge(
            ovv, index, how="left", on="PAR_file", left_index=False, right_index=False
        )
        for dt, gr in ovv_index.groupby(by=["EXP_date", "EXP_dir"]):
            exp_date, exp_dir = dt
            date = pd.to_datetime(exp_date).strftime("%Y_%m_%d")
            xl_out = IndexDir.joinpath(
                "{0}_{1}_{2}.xlsx".format(date, Path(exp_dir).name, len(gr))
            )
            xl_out.parent.mkdir(parents=True, exist_ok=True)
            gr.to_excel(xl_out)
            logger.info(
                "Merged Index with OVV from groupby and exported to: {0}".format(xl_out)
            )


def CleanUpCrew(
    EC_index=pd.DataFrame(),
    exp_folder="",
    filename_part="",
    list_of_files=[],
    delete=False,
    **kwargs,
):

    _removed = []

    if not isinstance(list_of_files, list):
        list_of_files = list(list_of_files)
    #    filename_part
    def fln_del_log(dd, fln, delete, _removed):
        if delete:
            fln.unlink()
        _removed.append([dd, fln, dt.datetime.now(), "unlink", delete])
        return _removed

    if not EC_index.empty:
        for dd in EC_index.Dest_dir.unique():
            if (exp_folder or filename_part) and not list_of_files:
                old_N2_scans = dd.joinpath(exp_folder)
                if old_N2_scans.is_dir():
                    if filename_part and exp_folder:
                        for fln in list(old_N2_scans.rglob(f"*{filename_part}*")):
                            _removed = fln_del_log(dd, fln, delete, _removed)
                    elif filename_part and not exp_folder:
                        for fln in list(old_N2_scans.rglob(f"*{filename_part}*")):
                            #                    fln in old_N2_scans.iterdir():
                            if fln.is_file():
                                if filename_part in fln.stem:
                                    _removed = fln_del_log(dd, fln, delete, _removed)
                    else:
                        if delete:
                            _del = subprocess.run(["rm", "-rf", str(old_N2_scans)])
                        _removed.append(
                            [dd, old_N2_scans, dt.datetime.now(), _del, delete]
                        )
            elif list_of_files and not (exp_folder or filename_part):
                for fln in list_of_files:
                    _removed = fln_del_log(dd, fln, delete, _removed)
            else:
                pass

    elif EC_index.empty and list_of_files:
        for fln in list_of_files:
            _removed = fln_del_log(Path(fln).parent, fln, delete, _removed)
    else:
        pass

    removed_log = pd.DataFrame(
        _removed,
        columns=[
            "Dest_dir_index",
            "target_Dest_dir",
            "date_now",
            "subprocess",
            "deleted",
        ],
    )
    if not removed_log.empty:
        removed_log.to_excel(
            FindExpFolder("VERSASTAT").IndexDir.with_name(
                f"{dt.date.today():%Y-%m-%d}_CleanUp_{exp_folder}-{filename_part}_{len(list_of_files)}_{FileOperations.version}.xlsx"
            )
        )
    else:
        print("Nothing deleted, log empty")


#    FolderOps.FindExpFolder('VERSASTAT').IndexDir.mkdir(parents=True,exist_ok=True)
#    self.index = EC_index
#    EC_index = ECRunOVV(load=1).index.sort_values('PAR_date',ascending=False)
#    #    run = input('Want to start the fitting run?')
#    run = 'ytest ORR N2_act'
# run = 'ytest short EIS'
#    if 'y' in run:
#        PAR_init_text()
#        EC_index = MainEC.EC_Analysis_Input(True)
#        PostDestDir = FolderOps.FindExpFolder('VERSASTAT').DestDir.joinpath('PostEC')
#%%
#            RHE_OCP_data = RHE_index[1], RHE_index[2]
# Forget about HDF files for now 21.01.19 ====
#            try:
#                RHE_ovv.to_hdf(PathDB,str(dest_dir.joinpath('RHE_OVV').as_posix()),table=True,mode='w')
#            except Exception == 'HDF5ExtError':
#                os.remove(PathDB)
#                RHE_ovv.to_hdf(PathDB,str(dest_dir.joinpath('RHE_OVV').as_posix()),table=True,mode='w')
#             === ====#  , 'dest_dir' : dest_dir})
#        Experiments, Gases  = ovv.PAR_exp.unique(), ovv.Gas.unique()
#        dest_dir.mkdir(parents=True,exist_ok=True)
#        outPF_ovv = dest_dir.joinpath('OVV_'+dest_dir.name+'.xlsx')
#        FileHelper.FileOperations.CompareHashDFexport(ovv,outPF_ovv)
#                ovv.to_excel(outPF_ovv)
