# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 15:41:04 2018

@author: User
"""
# from FindExpFolder import *

from functools import partial
import multiprocessing
from datetime import datetime
import os
from os import cpu_count
import hashlib

from pathlib import Path
import re

# from FindExpFolder import *
# FindExpFolder('VERSASTAT').PARfiles
import concurrent.futures
from multiprocessing import Pool

# import lxml
import sys

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper.FindSampleID import GetSampleID

# from file_py_helper import PostChar

import logging

logger = logging.getLogger(__name__)

from .file_inspector import PAR_file_parser

#%%
def getting_date_from_PARf(PARf, exp_dir):
    exp_date = ""
    up_dir = -2
    while exp_date == "":
        #                print(up_dir)
        try:
            exp_date = pd.to_datetime(
                exp_dir.split("_")[0], format="%d.%m.%Y", errors="raise"
            )
        #            exp_date = pd.to_datetime(exp_dir.split('_')[0],errors='ignore')
        except Exception as e1:
            #                faillst1.append((PARf,exp_dir,exp_date))
            try:
                exp_date = pd.to_datetime(
                    "_".join(exp_dir.split("_")[0:3]), format="%Y_%m_%d", errors="raise"
                )
            except Exception as e2:
                #                exp_date = pd.to_datetime(exp_dir.split('_')[0],format='%d.%m.%Y',errors='raise')
                exp_dir = PARf.parts[-up_dir]
                up_dir -= 1
                logger.warning(
                    "EXP date classifier error: {0} for {1}".format(
                        str(e1) + str(e2), exp_dir
                    )
                )
        #                        exp_date = pd.to_datetime(PARf.stat().st_ctime,unit='s')
        #                    faillst2.append((PARf,exp_dir,exp_date))
        if up_dir == -(len(PARf.parts) - 1):
            exp_date = pd.to_datetime(PARf.stat().st_ctime, unit="s")
            logger.warning(
                "EXP date classifier error setting stat ctime: for {0}".format(exp_dir)
            )
    #                 faillst3.append((PARf,exp_dir,exp_date))
    return exp_date


class EC_classifier_multi_core:
    """Class that initiates the loop over all experimental .PAR files and indexes each of them.
    Tries to run a multi core process for speed.
    attr: par_files_run
    func: main_run
    """

    def __init__(
        self, exp_folder_search=None, pre_load_RunOVV=pd.DataFrame(), main_run=False
    ):
        self.exp_folder_search = exp_folder_search
        self.pre_load_RunOVV = pre_load_RunOVV
        self.get_PAR_files_run()
        if main_run:
            self.main_run()
        pass

    #        self.get_PAR_files_run()

    def get_PAR_files_run(self):
        par_folders, par_files, class_db_file_path = FindExpFolder(
            "VERSASTAT"
        ).EC_find_folders()
        if self.exp_folder_search:
            logger.info("Classify files only from exp folder")
            par_files_run = [i for i in par_files if self.exp_folder_search in str(i)]
        if not self.pre_load_RunOVV.empty:
            logger.info("Classify files only from pre loaded OVV")
            present_par_files = [
                str(i) for i in self.pre_load_RunOVV.PAR_file.to_list()
            ]
            par_files_run = [i for i in par_files if not str(i) in present_par_files]
        else:
            logger.info("Classify all .PAR files from folders")
            par_files_run = par_files
        logger.info(
            f"All PAR files: {len(par_files)}\nMissing pfiles: {len(par_files_run)}"
        )
        #        return par_files_run,exp_folder_search, pre_load_RunOVV
        #        self.pre_load_RunOVV = pre_load_RunOVV
        #        self.exp_folder_search = exp_folder_search
        self.par_files_run = par_files_run

    def main_run(self, _run_multi=True, out_path=""):
        __spec__ = None

        td = datetime.today().date()
        if not out_path:
            out_path = FindExpFolder("VERSASTAT").DestDir.joinpath(
                f"EXP_PAR_class_{td.year}-{td.month}-{td.day}.xlsx"
            )

        ### RUNNING of classifier
        #        _run_multi = False
        if _run_multi == True:
            with multiprocessing.Pool(os.cpu_count() - 2) as pool:
                try:
                    #                    results = pool.map(EC_classifier_multi_core.EC_PAR_file_check, self.par_files_run)
                    results = pool.map(PAR_file_parser, self.par_files_run)
                except Exception as e:
                    #                    print('FileHelper module not found:',e)
                    logger.error("Classifier multiprocessing error: {0}".format(e))
                    results = pool.map(PAR_file_parser, self.par_files_run)

        else:
            results = []
            for par_file in self.par_files_run:
                _res = PAR_file_parser(par_file)
                results.append(_res)

        out = pd.DataFrame([i.parse_result for i in results])
        logger.info(f"Classifier processing finished len:{len(out)}")
        #        for pf in par_files_run:
        #            EC_classifier_multi_core.EC_PAR_file_check(pf)

        if (
            "DataFrame" in str(type(self.pre_load_RunOVV))
            and not out.empty
            and not self.pre_load_RunOVV.empty
        ):
            try:
                out_runovv = pd.concat([self.pre_load_RunOVV, out], axis=0, sort=False)
                logger.info(
                    "Classifier multiprocessing succes merged with pre load runovv: {0}".format(
                        out_path
                    )
                )
            except:
                out_runovv = out
        elif self.pre_load_RunOVV.empty and not out.empty:
            out_runovv = out
            logger.warning(
                "Classifier pre load empty, out is good! : {0}".format(out_path)
            )
        elif out.empty and not self.pre_load_RunOVV.empty:
            out_runovv = self.pre_load_RunOVV
            logger.error(
                "Classifier ERROR out empty and not Preload empty: {0}".format(out_path)
            )
        else:
            out_runovv = pd.DataFrame()
            logger.error(
                "Classifier ERROR out and pre load runovv both empty: {0}".format(
                    out_path
                )
            )

        if not out_runovv.empty:
            grhash = out_runovv.groupby("PAR_hash")
            hash_dups = [(nm, gr) for nm, gr in grhash if len(gr) > 1]
            logger.info(
                "Classifier multiprocessing index contains {0} duplicates".format(
                    len(hash_dups)
                )
            )
        try:
            out_runovv.to_excel(out_path)
            logger.info(
                "Classifier multiprocessing succes out path: {0}".format(out_path)
            )
            final_out_path = out_path

        except Exception as e:
            print("error: Output PAR OVV, %s" % e)
            logger.warning(
                "Classifier multiprocessing succes out path error {0}: {1}".format(
                    e, out_path
                )
            )
            final_out_path = out_path.parent.joinpath(
                f"EXP_PAR_class_{td.year}-{td.month}-{td.day}.xlsx"
            )
            out_runovv.to_excel(final_out_path)
        print("RunEC classifier end")
        self.index_path, self.index = final_out_path, out_runovv

    @staticmethod
    def _deprc_EC_PAR_file_check(PARf):
        try:
            DestDir = FindExpFolder("VERSASTAT").DestDir
        except:
            DestDir = FindExpFolder("VERSASTAT").DestDir
        PARf = Path(PARf)

        if PARf.is_file():
            exp_dir = PARf.parent.parts[-1]
            exp_dir_path = PARf.parent

        elif PARf.is_dir():
            exp_dir = PARf.name
            exp_dir_path = PARf.parent

        else:
            logger.error(
                "EC PAR file check. PAR file: {0} is neither file nor dir?!".format(
                    PARf
                )
            )
            exp_dir, exp_dir_path = PARf.parent.parts[-1], PARf.parent
        #        print(PARf)
        try:
            method = GetSampleID(PARf)
            sID_file, sID_dict = method.SampleID, method.SampleID_dict
        except Exception as e:
            logger.error(
                "EC PAR file check. PAR file: {0} FindSampleID error {1}".format(
                    PARf, e
                )
            )
            method = GetSampleID(PARf.parent)
            sID_file, sID_dict = method.SampleID, method.SampleID_dict
        basepf = PARf.stem
        basepf_split = basepf.split("_")

        exp_date = getting_date_from_PARf(PARf, exp_dir)
        #        try:
        #            exp_date = pd.to_datetime(exp_dir.split('_')[0],format='%d.%m.%Y',errors='ignore')
        ##            exp_date = pd.to_datetime(exp_dir.split('_')[0],errors='ignore')
        #        except Exception as e1:
        #              try:
        #                  pd.to_datetime('_'.join(a.split('_')[0:3]),format='%Y_%m_%d',errors='ignore')
        #              except Exception as e2:
        #                  logger.warning('EXP date classifier error: 1st:{0}, 2nd: {1} for {2}'.format(e1, e2, exp_dir))
        #                  exp_date = pd.to_datetime(PARf.stat().st_ctime,unit='s')
        #            try:
        #                exp_date = pd.to_datetime(exp_dir.name.split('_')[0],format='%d.%m.%Y',errors='ignore')
        #            except:
        try:
            dl, el = DestDir.parts, exp_dir_path.parts
            diffeldl = len(el) - len(dl)
            if diffeldl > 0:
                Dest_dir = Path(*(dl + el[-diffeldl::]))
            else:
                Dest_dir = DestDir / "EC_plots" / el[-1]
        except Exception as e:
            print(e, exp_dir)
            Dest_dir = DestDir / "EC_plots" / exp_dir.parts[-1]
        out = []
        ol = []
        gas, tpnm, electrode = "NaN", "NaN", "NaN"
        #            for pf in par_files[::]+list(itertools.chain(*[a for a in sub_parf])):
        #            print(pf)
        #        try:
        #            SpID = FileHelper.FindSampleID(PARf).SampleID
        #        except:
        #            SpID = basepf
        #    #        with open(PARf) as N2f:
        #            N2_read = N2f.read()
        N2_read = PARf.read_text()
        N2_soup_par = BeautifulSoup(N2_read, "lxml")
        hash_f = FileOperations.take_hash(PARf)
        action0 = FindExpFolder.get_DF_soup_part(N2_soup_par.action0, hash_f)
        experiment = FindExpFolder.get_DF_soup_part(N2_soup_par.experiment, hash_f)
        try:
            comment_act0_DF = action0["Comment"].values[0]
        except:
            comment_act0_DF = "NoneDetermined"
        #        instr_DF = FileHelper.FindExpFolder.get_DF_soup_part(N2_soup_par.instrument,hash_f)
        #        date_df = FileHelper.FindExpFolder.get_DF_soup_part(N2_soup_par.experiment,hash_f)
        #        try:
        #            actions = [i.name for i in N2_soup_par.find_all(name=re.compile('action.'))][1::]
        #            action1_name = FileHelper.FindExpFolder.get_DF_soup_part(N2_soup_par.action1,N2_soup_par.action1.name)['Name'][0]
        #        except:
        #            actions, action1_name = 'NaN', 'NaN'
        try:
            experiment["DTC"] = pd.to_datetime(
                experiment["DateAcquired"] + " " + experiment["TimeAcquired"]
            )
            PAR_date = experiment["DTC"][0]

        except Exception as e:
            logger.warning("EC_PAR_file_check {0} error: {1}".format(PARf, e))
            #            print(e)
            PAR_date = pd.to_datetime(PARf.stat().st_ctime, unit="s")
        #        filesize = os.path.getsize(PARf)
        #        SampleID = FindSampleID(PARf).SampleID
        try:
            if exp_date.date() == PAR_date.date():
                exp_DATE_match = "yes"
            else:
                exp_DATE_match = "no"
        except Exception as e:
            logger.warning(
                "EC_PAR_file_check {0} match date error: {1}".format(PARf, e)
            )
            #            print(e)
            exp_DATE_match = "NaN"
        gas, tpnm = GetSampleID.determine_Gas_Exp_from_filename(
            basepf, comment_act0_DF, False
        )
        postAST = GetSampleID.determine_postAST_from_filename(basepf, verbose=False)
        #        pf.split('\\')[-1]
        pH = GetSampleID.determine_pH_from_filename(PARf, False)
        if pH["Electrolyte"][0] == "FALSE":
            pH = GetSampleID.determine_pH_from_filename(exp_dir, False)

        loading_name, loading_cm2 = GetSampleID.ink_loading_from_filename(
            basepf_split, PAR_date
        )
        #        if 'ring' in basepf or 'Ring' in basepf:
        #            electrode = 'Pt_ring'
        if "Pt-ring" in basepf or "ring" in basepf or "Ring" in basepf:
            electrode = "Pt_ring"
        elif "Disk" in basepf or "disk" in basepf:
            electrode = "disk"
        else:
            electrode = "Unknown"
        #        exps = ['RHE|HER_OCP|RHE_OCP|OCP_RHE','HPRR','N2','O2','_V3',' _V3F','_disk','_ring']
        #        exp_names = ['RHE','HPRR','N2','O2',' V3',' V3F','disk','ring']
        cr_date = pd.to_datetime(PARf.stat().st_ctime, unit="s")
        mod_date = pd.to_datetime(PARf.stat().st_mtime, unit="s")
        PAR_cr_Delta = (pd.to_datetime(mod_date) - PAR_date).seconds
        #        [exp_dir,Dest_Dir,exp_date,PAR_date,exp_DATE_match,exp_DATE_match,PARf,hash_f,gas,tpnm,PARf.stat().st_size,basepf,SampleID,
        #                   cr_date,mod_date,PAR_cr_Delta,electrode,pH['pH'][0],pH['Electrolyte'][0],
        #                   comment_act0_DF]
        dr = {
            "EXP_dir": exp_dir_path,
            "Dest_dir": Dest_dir,
            "EXP_date": exp_date,
            "PAR_date": PAR_date,
            "EXP_PAR_date_match": exp_DATE_match,
            "EXP_PAR_folder_match": sID_dict["sID_Match_file_dir"],
            "PAR_file": PARf,
            "PAR_hash": hash_f,
            "Gas": gas,
            "PAR_exp": tpnm,
            "postAST": postAST,
            "filesize": PARf.stat().st_size,
            "basename": basepf,
            "SampleID": sID_file,
            "SampleID_folder": sID_dict["sID_Folder"],
            "Creation_date": cr_date,
            "LastMod_date": mod_date,
            "Delta_PAR_LastMod": PAR_cr_Delta,
            "Electrode": electrode,
            "pH": pH["pH"][0],
            "Electrolyte": pH["Electrolyte"][0],
            "Comment": comment_act0_DF,
            "Loading_name": loading_name,
            "Loading_cm2": loading_cm2,
        }
        return dr


#        dtlst.append((PARf,exp_date))

#    out.groupby(by='EXP_dir').groups
if __name__ == "__main__":
    #    from FileHelper import *
    __spec__ = None
    print("==== NO EC CLASSIFICATION ===")
