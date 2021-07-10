## -*- coding: utf-8 -*-
# """
# Created on Sat Jul 22 15:05:57 2017
# @author: David Wallace
# @email: wallace@ese.tu-darmstadt.de
# @
# This module reads PAR files in a certain folder and tries to analyze the experimental data
# from N2,ORR,OER,HPRR measurements.
# """
#
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.lines as mlines
##from matplotlib import cm
# from matplotlib import gridspec
##from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
# import os
# import glob
# import sys
#
# from collections import namedtuple
#
# import numpy as np
# import pandas as pd
# import scipy.signal
# from scipy.stats import linregress
##import xarray as xr
##import statsmodels.api as sm
##from sklearn.linear_model import LinearRegression
##from math import isclose
##import shlex
##from PyQt import QtGui.QFileDialog
##import tkinter
##from tkinter.filedialog import *
# from bs4 import BeautifulSoup
##import io
##import errno
# import re
##import platform
##import socket
##import concurrent.futures
##from multiprocessing import Pool
##from functools import partial
# from itertools import repeat
# import multiprocessing
##from collections import Counter
##import scipy.optimize
# from pathlib import Path
# import logging
##import difflib
##import time
##import h5py
# from datetime import datetime
##import pytables
# matplotlib.rcParams.update({'font.size': 16})
# print('TODO: New runPAR fix imports and definitions...!!!')
# try:
#    import statsmodels.formula.api as smf
##    import statsmodels as sm
##    import seaborn as sns
# except Exception as e:
#    print('No modules: %s'%e)
# print('Name:', __name__)
# if __name__ == "__main__":
#    print(__file__)
##    import experiments.EIS.eis_run_ovv
#    from experiments.Loader import set_OCP_RHE
#    FH_path = Path(__file__).parent.parent.parent.joinpath('FileHelper')
#    sys.path.append(str(FH_path))
#    sys.path.append(str(Path(__file__).parent.parent.joinpath('indexer')))
#    sys.path.append(str(Path(__file__).parent.parent.parent))
#    sys.path.append("..")
##    print(sys.path)
#    import FileHelper
#    import prepare_input
#    from experiments.EIS import eis_run_ovv
##    import experiments
#
#
##    from experiments.EIS import *
##    import runEC
#
# else:
##    import FileHelper
##    import prepare_input
#    print(f'Do not run this:{__file__}')
#
##from functools import
##from FolderOrganizing.FileHelper import FindExpFolder
##from FolderOrganizing.RunEC_classifier import EC_classifier_multi_core
##import FolderOrganizing.FileHelper
##import FolderOrganizing.RunEC_classifier
##from FolderOrganizing import FileHelper
##from FolderOrganizing import RunEC_classifier
##from FolderOrganizing import PAR_EIS_fit_V2
##import PAR_postEC_v6
## from lmfit import minimize, Parameters, Parameter,Model, report_errors
#
#
##import scipy.fftpack
##from pyfinance.ols import OLS, RollingOLS, PandasRollingOLS
##%%
# EvRHE = 'E_AppV_RHE'
# def start_logging():
#    # Gets or creates a logger
#    logger = logging.getLogger(__name__)
#
#    # set log level
#    logger.setLevel(logging.INFO)
#
#    # define file handler and set formatter
#    file_handler = logging.FileHandler(FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_DW_logfile.log'))
#    formatter    = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s')
#    file_handler.setFormatter(formatter)
#
#    # add file handler to logger
#    logger.addHandler(file_handler)
#    logger.warning('=========== Started logging {0}...  ============='.format(__name__))
#    return logger
#
#    #%%
#
##ExpParOVV['EXP_dir'] = [str(find_CS_parts(Path.cwd())[0].joinpath(find_CS_parts(Path(i))[1])) for i in ExpParOVV.EXP_dir.values]
##def find_CS_parts(xyPath):
###    xyPath = ExpParOVV.EXP_dir.values[2]
###    [str(find_CS_parts(Path.cwd())[0].joinpath(find_CS_parts(Path(i)))) for i in ExpParOVV.EXP_dir.values]
##    CSroot = Path(*Path(xyPath).parts[0:[(i,a) for i,a in enumerate(Path.cwd().parts) if a == 'CloudStation' or a == 'Cloudstation'][0][0]+1])
###    Path(*Path(ExpParOVV.EXP_dir.values[2]).parts[0:[(i,a) for i,a in enumerate(Path.cwd().parts) if a == 'CloudStation'][0][0]+1])
###    CSafter = Path(*Path(ExpParOVV.EXP_dir.values[2]).parts[[(i,a) for i,a in enumerate(Path.cwd().parts) if a == 'CloudStation'][0][0]::])
##    CSafter = Path(*Path(xyPath).parts[[(i,a) for i,a in enumerate(Path(xyPath).parts) if a == 'Experimental data'][0][0]::])
##    return CSroot, CSafter
##def find_CS_List(xyPathLst):
##    xyPathLst = ExpParOVV.EXP_dir.values
##    for x in xyPathLst:
##        Pre = Path(*Path.cwd().parts[0:[(i,a) for i,a in enumerate(Path.cwd().parts) if a == 'CloudStation' or a == 'Cloudstation'][0][0]])
##        mid = [(i,a) for i,a in enumerate(Path(x).parts) if a == 'CloudStation' or a == 'Cloudstation']
##        Path(*Path(a).parts[0:[(i,a) for i,a in enumerate(Path.cwd().parts) if a == 'Experimental data'][0][0]])
##        CSroot = Path(*Path(a).parts[0:[(i,a) for i,a in enumerate(Path.cwd().parts) if a == 'CloudStation'][0][0]+1])
##        CSafter = Path(*Path(a).parts[[(i,a) for i,a in enumerate(Path(a).parts) if a == 'Experimental data'][0][0]::])
# #%%
# class ECRunOVV:
#    ''' Run over a DateFrame that is an Overview(OVV) of EC experiment in a folder.'''
#    def __init__(self,**kwargs):
#        pass
#        if 'load' in kwargs:
#            OnlyRecentMissingOVV = prepare_input.MainEC.EC_Analysis_Input(True)
#            FileHelper.FindExpFolder('VERSASTAT').IndexDir.mkdir(parents=True,exist_ok=True)
#            self.index = OnlyRecentMissingOVV
##        self.ovv = ovv
#    EvRHE = 'E_AppV_RHE'
##    PathDB = FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_Files_class.hdf5')
#
#    @staticmethod
#    def MakeOVVperExpDir(arg, **kwargs):
##%%
##        HPRR_pars, Cdl_fit_pars, Cdl_fit_data  = pd.DataFrame([]), pd.DataFrame([]), pd.DataFrame([])
##        outP1,outP2 = pd.DataFrame(),pd.DataFrame()
##        dest_dir,arg.EC_index
#        EvRHE = 'E_AppV_RHE'
##        PathDB = FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_Files_class.hdf5')
#        ovv, RHE_ovv,status, All_ovv,Samples_ovv = pd.DataFrame([]), pd.DataFrame([]),pd.DataFrame([]),pd.DataFrame([]),pd.DataFrame([])
##        dest_dir = Path(gr['Dest_dir'].unique()[0])
##        save_dir = Path(FileHelper.FindExpFolder('VERSASTAT').DestDir,dest_dir.parts[-2],dest_dir.parts[-1])
## ===== Compare with saved STATUS to skip or not  =====
#        '''Check status in an HDF5 file'''
##            if exp_dir.N2 == 1 and exp_dir.HPRR == 1 or exp_dir.ORR == 1:
##                print('Skipped %s'%basedir)
##                continue
##            else:
##                print('Perform analysis on: %s' %basedir)
## ===== Set default values =====
####### === Start with overview and classification of files in the Experimental Folder ==== #######
#        '''Filter gr of OVV for failed measurements'''
#        gr = arg.grp
#        exp_dir = arg.EXP_dir
#        ovv = gr.loc[(gr.PAR_file == [i for i in gr.PAR_file.values if not ('template|templ|fail|Kopie|Copy') in  str(i)])]
##        ovv = gr.loc[(~gr['PAR_file'].str.contains('template|templ|fail')),:]
##        SampleIDfolder = FileHelper.FindSampleID.match_SampleID(exp_dir)
#        if not ovv.Date_PAR_EXP.empty:
#            ovv = ovv.loc[((ovv.Date_PAR_EXP > pd.Timedelta('-2 days 0 hours')) | (ovv.PAR_exp == 'RHE' )),:]
#            if ovv.empty:
#                logger.warning('MakeOVVperExpDir files filter by DATE ignored')
#                ovv = gr.loc[(gr.PAR_file == [i for i in gr.PAR_file.values if not ('template|templ|fail|Kopie|Copy') in  str(i)])]
#        ovv_nfiles,ovv_nfls = ovv.PAR_hash.nunique(),ovv.PAR_date.nunique()
#        singles = ovv.loc[ovv.duplicated(subset=['PAR_hash','PAR_date'],keep=False) == False]
#        dups = ovv.loc[ovv.duplicated(subset=['PAR_hash','PAR_date'],keep=False) == True]
#        dups_folder_match = dups.loc[ dups.EXP_PAR_folder_match.str.contains('yes')]
#        already_found = pd.concat([singles,dups_folder_match])
#        left_overs = ovv[~ovv.PAR_hash.isin(already_found.PAR_hash)].drop_duplicates(subset='PAR_hash')
#        if not left_overs.empty:
#            ovv_filtered = pd.concat([already_found,left_overs])
#        else:
#            ovv_filtered = already_found
#        number_filtered = len(ovv)-len(ovv_filtered)
#        if ovv_filtered.PAR_hash.nunique() == ovv_nfiles:
#            logger.warning('MakeOVVperExpDir files ({0}) were filtered out as duplicates {1}'.format((),exp_dir))
#            ovv = ovv_filtered
##        dups_folder_match = ovv.loc[(ovv.duplicated(subset=['PAR_hash'],keep=False) == True & ovv.EXP_PAR_folder_match.str.contains('yes')) & (ovv.duplicated(subset=['PAR_hash'],keep=False) == False)]
##        ovv.loc[(ovv.duplicated(subset=['PAR_hash'],keep=False) == False)]
##        set(gr.PAR_file.values) - set(ovv.PAR_file.values)
#        OVV_GR_len = np.abs(len(gr)-len(ovv))
#        if OVV_GR_len > 3:
#            logger.warning('MakeOVVperExpDir Many files ({0}) were filtered out by date or name from {1}'.format(OVV_GR_len,exp_dir))
#        SkippedLst = []
#        if ovv.empty:
#            logger.warning(f'MakeOVVperExpDir Skipped, OVV EMPTY FOR: {exp_dir}')
#            SkippedLst.append(exp_dir)
#
#        for n,r in ovv.iterrows():
#            if 'ORR' in Path(r.PAR_file).stem.split('_') and r.PAR_exp != 'ORR':
##                    print(n,r)
#                print(n, Path(r.PAR_file).name,r.PAR_exp,r.Gas)
#                logger.info('{0}'.format([n, Path(r.PAR_file).name,r.PAR_exp,r.Gas]))
#                ovv.loc[n,['PAR_exp']] = 'ORR'
#                ovv.loc[n,['Gas']] = 'O2'
#                OnlyRecentMissingOVV.loc[n,['PAR_exp']] = 'ORR'
#                OnlyRecentMissingOVV.loc[n,['Gas']] = 'O2'
#
####### === Matching ORR scan with N2_activations in OVV ==== #######
#        ovv['ORR_act_N2_bg'] = 'None' # Make new empty columns for filename of ORR background scan
#        ovv_N2_act_local = ovv.query('PAR_exp == "N2_act"')
#
##        ovv_N2_act_local = OnlyRecentMissingOVV.loc[(OnlyRecentMissingOVV.EXP_date == ovv.EXP_date.unique()[0])].query('PAR_exp == "N2_act"')
#        for n,r in ovv_N2_act_local.iterrows():
#            N2_dest_dir = Path(r.Dest_dir).joinpath('N2_scans/{0}'.format('_'.join([r.Electrolyte,r.SampleID,r.Loading_name])))
#            N2_dest_dir.mkdir(parents=True,exist_ok=True)
#            N2_act_BG = Path(N2_dest_dir.parent).joinpath(r.basename+'_BG.xlsx')
##            print(r.basename+'_BG.xlsx')
#            ovv.loc[n,'ORR_act_N2_bg'] = N2_act_BG
##              ovv.loc[ovv['PAR_exp'] == "ORR" & ovv.SampleID == r.SampleID & ovv.postAST == r.postAST,'ORR-act_N2-bg'] = ('ORR;'+r.PAR_file)
#            if not 'ring' in N2_act_BG.stem and 'Pt' not in r.Comment:
#
#                if not ovv.loc[((ovv['PAR_exp'] != "N2_act") & (ovv.SampleID == r.SampleID) & (ovv.postAST == r.postAST) &
#                         (ovv.Electrode == r.Electrode) & (ovv.Loading_cm2 == r.Loading_cm2) & (ovv.Electrolyte == r.Electrolyte))].empty:
#                    ovv.loc[((ovv['PAR_exp'] != "N2_act") & (ovv.SampleID == r.SampleID) & (ovv.postAST == r.postAST) &
#                         (ovv.Electrode == r.Electrode) & (ovv.Loading_cm2 == r.Loading_cm2) & (ovv.Electrolyte == r.Electrolyte)),'ORR_act_N2_bg'] = str(N2_act_BG)
##                    print('reset', r.basename,'to ',(N2_act_BG.name))
#                else:
#                    pass
##                    print('ovv fail ',r.basename)
#            else:
#                pass
##                print('fail 2',r.basename)
#        faillst = []
#        for n,r in ovv.query('(ORR_act_N2_bg == "None") & (PAR_exp == "ORR")').iterrows():
##            n,r = list(ovv.query('(ORR_act_N2_bg == "None") & (PAR_exp == "ORR") & (Electrode != "Pt_ring")').iterrows())[-1]
##            n,r = faillst[0][0],faillst[0][-1]
##            print(Path(r.PAR_file))
##            Path(r.PAR_file).parent.name.split('_')[0]
##            "N2_act",r.SampleID,r.postAST,r.Electrolyte,r.EXP_date, r.Loading_name
##            N2_dest_dir = Path(gr_N2_ovv.Dest_dir.iloc[0]).joinpath('N2_scans/{0}'.format('_'.join([Electrolyte,SampleID,Loading_name])))
##            N2_dest_dir.mkdir(parents=True,exist_ok=True)
##            N2_act_BG = Path(N2_dest_dir.parent).joinpath(N2_fn+'_BG.xlsx')
#            N2exMatch_ovv_local = ovv.loc[((ovv['PAR_exp'] == "N2_act") & (ovv.SampleID == r.SampleID) & (ovv.postAST == r.postAST)
#                                                & (ovv.Electrolyte == r.Electrolyte) & (ovv.EXP_date == r.EXP_date)
#                                                & (ovv.Loading_name == r.Loading_name))]
#            if len(N2exMatch_ovv_local) == 1: # look for local N2 background inside of the experimental folder
#                 N2_bg_from_locl_N2match = N2exMatch_ovv_local.ORR_act_N2_bg.values[0]
#                 ovv.loc[n,'ORR_act_N2_bg'] = N2_bg_from_locl_N2match
#            else: # look for N2 background outside of the experimental folder
#                N2exMatch = OnlyRecentMissingOVV.loc[((OnlyRecentMissingOVV['PAR_exp'] == "N2_act") & (OnlyRecentMissingOVV.SampleID == r.SampleID) & (OnlyRecentMissingOVV.postAST == r.postAST)
#                                                    & (OnlyRecentMissingOVV.Electrolyte == r.Electrolyte) & (OnlyRecentMissingOVV.EXP_date == r.EXP_date)
#                                                    & (OnlyRecentMissingOVV.Loading_name == r.Loading_name))]
#
#                if len(N2exMatch) <= 1:
#                    N2exMatch = OnlyRecentMissingOVV.loc[((OnlyRecentMissingOVV['PAR_exp'] == "N2_act") & (OnlyRecentMissingOVV.SampleID == r.SampleID) & (OnlyRecentMissingOVV.postAST == r.postAST)
#                                                    & (OnlyRecentMissingOVV.Electrolyte == r.Electrolyte.split('+')[0]) & (OnlyRecentMissingOVV.EXP_date == r.EXP_date)
#                                                    & (OnlyRecentMissingOVV.Loading_name == r.Loading_name))]
#    #            if r.Electrode == 'Pt_ring' and len(N2exMatch) > 1:
#    #                N2exMatch = N2exMatch.iloc[0]
#                if len(N2exMatch) == 1:
#                    N2matchPath = Path(*Path(N2exMatch.PAR_file.values[0]).parts[-2::])
#    #                print('N2_bg succes! (%s) added in folder for: %s' %(N2matchPath,r.basename))
#                    logger.info('N2_bg succes! ({0}) added in folder {2} for: {1}'.format(N2matchPath,r.basename, r.Dest_dir.name))
#                    N2exMatch.Dest_dir = r.Dest_dir
#                    N2_ext_BG = Path(r.Dest_dir).joinpath(N2exMatch.basename.values[0]+'_BG.xlsx')
#                    N2exMatch = N2exMatch.assign(**{'ORR_act_N2_bg' : N2_ext_BG})
#                    ovv.loc[n,'ORR_act_N2_bg'] = N2_ext_BG
#                    ovv = ovv.append(N2exMatch)
#                else:
#                    N2_ext_BG = 'fail'
#                    faillst.append([n,r.basename,N2_ext_BG,r])
#                    if r.Electrode != 'Pt_ring':
#                        if r.PAR_exp == 'ORR':
#                            logger.warning('MakeOVVperExpDir N2 BG match fail, N2_act_BG ORR, {0}\n for folder {1}' .format(r.basename,dest_dir))
#    #            print(r.basename,'\n',Path(N2_ext_BG).name)
#        ovv = ovv.drop_duplicates(subset=['PAR_hash','PAR_file'])
#        if ovv.empty:
##            print('OVV empty skipped')
#            logger.warning('OVV empty skipped: {0}'.format(exp_dir))
#        RHE_index = get_RHE_OCP(ovv)
#        ovv = RHE_index[0]
#        #%%
#        return ovv
####### === EC Analysis Run over test1 in multiprocess ==== #######
#    @staticmethod
#    def EC_Analysis_Run(EC_run_slice,EC_index,multi_run=False):
#        SkippedLst, results,out = [], [], []
#        EC_run_groupby = ['EXP_date','EXP_dir']
#        grp_date_dir = EC_run_slice.groupby(by=EC_run_groupby)
#        grp_arg_lst = grp_date_dir.groups.items()
##        zip_grp_arg_lst = zip(grp_arg_lst,repeat(EC_index))
#        zip_grp_arg_lst = [(i,{'EC_index' : EC_index,'EC_run_groupby'  : EC_run_groupby}) for i in grp_arg_lst]
#        run_group = namedtuple('Run_groups' ,'EXP_date EXP_dir grp EC_index')
#        run_groups = [run_group(n[0],n[1],gr,EC_index) for n,gr in grp_date_dir]
#
#        grp_arg_lst
##        for exp_dir,gr in test3.groupby(by='EXP_dir'):
##        test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV['EXP_dir'].str.contains('06.03.2018_DW17_HPRR_0.1MHClO4_RRDE22775')] # TESTINGG!!
##        grp_arg_lst_ddd =[(i,OnlyRecentMissingOVV.to_dict()) for i in grp_arg_lst]
#        if multi_run == True:
##            grp_arg_lst = test1.groupby(by=['EXP_dir']).groups.items() # TODO check grouping EXP_date
##            grp_arg_lst_ovv = zip(grp_arg_lst,repeat(EC_index))
#            logger.info('PAR DW Classifier multiprocessing start: {0}, length {1}'.format(multi_run,len(grp_arg_lst)))
#            with multiprocessing.Pool(os.cpu_count()-1) as pool:
#                try:
#                    results = pool.starmap( ECRunOVV.EC_Analysis_Run, run_groups)
#                except Exception as e:
#                    print('PAR DW Classifier multiprocessing error: {0}'.format(e))
#                    logger.error('PAR DW Classifier multiprocessing error: {0}'.format(e))
#        #                    results = pool.map(EC_classifier_multi_core.EC_PAR_file_check, FolderOrganizing.FileHelper.FindExpFolder('VERSASTAT').PARfiles)
#        elif multi_run == False:
#            logger.info('PAR DW Classifier single core start: multi = False, length {0}'.format(len(grp_arg_lst)))
##            for nm,gr in test1.groupby(by=['EXP_dir']):
#            for arg in run_groups:
##                 arg = [nm,gr.index]
#                 res_arg = ECRunOVV.EC_Analysis_Run_groups(arg)
#                 results.append(res_arg)
#        out.append(results)
##        results_df = pd.concat(results)
#        return out
##        for dt,gr in test1.groupby(by=['EXP_date','EXP_dir']):
##            exp_date,exp_dir = dt
##        for DateSample,gr in test1.groupby(by=['EXP_date','SampleID']):
##            exp_dir = gr.EXP_dir.unique()[-1]
####### === EC Analysis Run over test1 grouped by Date ==== #######
#    @staticmethod
#    def EC_analysis_run_wrapper(arg):
#        args, kwargs = arg
#        return ECRunOVV.EC_Analysis_Run(*args, **kwargs)
#
#    @staticmethod
#    def EC_Analysis_Run_groups(arg, print_console=True, **kwargs):
##        EC_index = kwargs.get('EC_index',pd.DataFrame())
##        EC_run_groupby = kwargs.get('EC_run_groupby',[])
#        skip,results = False,[]
##        exp_date,
##         = arg[0]
##        gr_idx = arg[1]
##        OnlyRecentMissingOVV = Index
##        print('date',exp_date)
##        print('exp dir',exp_dir)
##        gr = OnlyRecentMissingOVV.loc[gr_idx]
#        dest_dir = arg.grp.Dest_dir.unique()[0]
#        exp_dir= arg.EXP_dir
#
#        if 'old' in Path(arg.EXP_dir).name:
#            logger.info('Skipped: {0}'.format(arg.EXP_dir))
#            skip = True
#            #                print('Skip old dir, %s' %Path(exp_dir).name)
#        if len(arg.grp) < 1:
#            logger.warning('Skipped for too few files {0} in: {1}'.format(len(arg.grp),arg.EXP_dir))
#            skip = True
##                print('Skipped, Too few files (%s) in %s. Skipped' %(len(gr),exp_dir))
##        gr,exp_dir = test1.groupby(by=['EXP_date']).get_group('2019-07-19 00:00:00'), test1.EXP_dir.unique()[1]
#        if skip == False:
#            try:
#                ovv = ECRunOVV.MakeOVVperExpDir(arg)
#            except Exception as e:
#                ovv = arg.grp
#                logger.error('EC_Analysis_Run_Date ovv preparation fail: {0}\n because {1}'.format(exp_dir,e))
#            try:
#                results = ECRunOVV().OVV_loop(ovv)
#            except Exception as e:
#                logger.warning('EC_Analysis_Run_Date OVV loop Fail: {0}\n because {1}'.format(exp_dir,e))
#    #                with multiprocessing.Pool(os.cpu_count()-1) as pool:
#    #                    try:
#    #                        results = pool.map( ECRunOVV.OVV_loop, ovv)
#    #                    except Exception as e:
#    #                        print('PAR DW Classifier multiprocessing error: {0}'.format(e))
#    #                        logger.error('PAR DW Classifier multiprocessing error: {0}'.format(e))
#    ##                    results = pool.map(EC_classifier_multi_core.EC_PAR_file_check, FolderOrganizing.FileHelper.FindExpFolder('VERSASTAT').PARfiles)
#    #                out.append(results)
#    #                print('Success: ', exp_dir)
#                results = []
#                if results:
#                    MainEC.MergeIndexOvv(results,ovv)
#                    logger.info('EC_Analysis_Run_Date Success Index: {0}'.format())
#                else:
#                    logger.warning('EC_Analysis_Run_Date Fail for empty results/Index: {0}, {1}'.format(results, exp_dir))
#            except Exception as e:
#                logger.warning('EC_Analysis_Run_Date Fail: {0}\n because {1}'.format(exp_dir,e))
#        print(results)
#        return results
#
####### === EC Analysis Run over test1 grouped by Date ovv and loop over PAR_exp types ==== #######
#    @staticmethod
#    def OVV_loop(ovv,plotEIS = True):
#        '''Start of Main Analysis Loop Where Each Type of Experiment from OVV is processed
#            and index files for output are collected ==== #######  '''
##        ovv = test1
#        if 'PAR_exp' in ovv.columns:
#            index_info = []
#            ovv_Dest_dir =  Path(ovv.Dest_dir.unique()[0])
#            ExpTypes_gr = ovv.groupby(by='PAR_exp')
#            for exp,gr in ExpTypes_gr:
#                logger.info('OVV loop starting: {0}'.format(exp))
#                if 'N2_act' in exp:
##                    print('N2_act_SKIP')
#                    index = ECRunOVV.N2_act(ovv)
#                    index['Script_run_date'] = pd.datetime.now()
#                    index_info.append(index)
#                elif 'EIS' in exp:
#                    EIS_skip_set,EIS_use_prepars_set, FreqLim_set = False,True, 10000
#                    print('Skip EIS: %s. Use prePars: %s\n Frequency limit: %.0f'%(EIS_skip_set,EIS_use_prepars_set,FreqLim_set))
##                    exp, gr = 'EIS', ExpTypes_gr.get_group('EIS')
#                    if EIS_skip_set == False:
#                        indexes, indexes_path_target = eis_run_ovv.eis_run_group_ovv(exp,gr,ovv,EIS_use_prepars = EIS_use_prepars_set, FreqLim = FreqLim_set,EIS_skip = EIS_skip_set,EIS_plot_combined=True )
##                        ovv_Dest_dir.joinpath('{})
#                        index_info.append(indexes)
#                elif 'HPRR' in exp:
##                    print('HPRR_skip')
##                    print('Skip EIS: %s. Use prePars: %s\n Frequency limit: %.0f'%(EIS_skip,EIS_use_prepars,FreqLim_set))
##                    if EIS_skip == False:
#                    index = ECRunOVV.HPRR(exp,ExpTypes_gr.get_group('HPRR'),ovv)
#                    index['Script_run_date'] = pd.datetime.now()
#                    index_info.append(index)
#                elif 'OER' in exp:
#                    index = ECRunOVV.OER(exp,ExpTypes_gr.get_group('OER'),ovv)
#                    index['Script_run_date'] = pd.datetime.now()
#                    index_info.append(index)
#                elif 'HER' in exp:
#                    index = ECRunOVV.HER(exp,ExpTypes_gr.get_group('HER'),ovv)
#                    index['Script_run_date'] = pd.datetime.now()
#                    index_info.append(index)
#
#                elif 'ORR' in exp or 'O2_nan' in exp:
#                    index = ECRunOVV.ORR(ovv)
#                    index['Script_run_date'] = pd.datetime.now()
#                    index_info.append(index)
#
#                elif 'RHE' in exp:
#                    pass
#                else:
##                    print('No run, unknown experiment type:', exp)
#                    logger.info('No run, unknown experiment type: {0}'.format(exp))
#                Index_exp = pd.concat([i for i in index_info],ignore_index=True,sort=False)
#                index_local_exp_path = ovv_Dest_dir.joinpath('index_{0}.xlsx'.format(exp))
#                index_local_exp_target = FileHelper.FileOperations.CompareHashDFexport(Index_exp,index_local_exp_path)
#                index_folder_exp_path = FileHelper.FindExpFolder('VERSASTAT').IndexDir.joinpath(ovv_Dest_dir.name + '_index_{0}_{1}.xlsx'.format(exp,len(ovv)))
#                index_folder_exp_target = FileHelper.FileOperations.CompareHashDFexport(Index_exp,index_folder_exp_path)
#
##            pd.concat([pd.DataFrame(i,index=[0]) for i in index_out],ignore_index=True,sort=False)
##            Index = pd.concat([i for i in index_info],ignore_index=True,sort=False)
##            index_local_path = ovv_Dest_dir.joinpath('index.xlsx')
#
##            if not Index.empty:
##            Index.to_excel(index_local_path)
##            Index.to_excel(index_folder_path)
#            logger.info('Main OVV loop,  Index  saved local to: {0}\n and in general Folder to '.format(index_local_exp_target,index_folder_exp_target))
#        else:
#            print('Your overview does not contain a PAR_exp column')
#            logger.error('Your overview does not contain a PAR_exp column: {0}'.format(Path(ovv.PAR_file.iloc[0]).parent))
#            index_folder_exp_target = []
#        return index_folder_exp_target
#
##%%
#### ============== ###
# def MultiMain():
#    __spec__ = None
#    with multiprocessing.Pool(os.cpu_count()-1) as pool:
#        try:
#            results = pool.map(prepare_input.MainEC().EC_Analysis_Input())
#        except Exception as e:
#            print('MultiprocessError module not found:',e)
##                results = pool.map(EC_classifier_multi_core.EC_PAR_file_check, FolderOrganizing.FileHelper.FindExpFolder('VERSASTAT').PARfiles)
#
##%%
# def index_run_selection(run:str,EC_index):
#    if 'y' in run:
#        if 'ya' in run or run == 'yes all' or run == 'y std':
#            test1 = EC_index
#    #        test3 = OnlyRecentMissingOVV.query('(EXP_date >= 20190301)  & ((SampleID >= "DW01") | (SampleID == "JOS12") | (SampleID == "JOS13") | (SampleID == "JOS14") | (SampleID == "JOS15")) & (PAR_exp == "EIS")')
##            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20190103)  & (SampleID == "JOS15") & (PAR_exp == "EIS")')
#        elif run == 'ys':
#            test1 = EC_index[EC_index.EXP_date.isin(EC_index.query('PAR_exp == "HPRR" ').EXP_date)]
#        elif 'yrecent' in run:
#             test1 = EC_index.query('(PAR_date >= 20200220)')
#             if 'slicedate' in run:
#                 test1 = test1.loc[(test1.Date_PAR_EXP > pd.Timedelta('-2 days 0 hours')),:]
#        elif 'test' in run:
##            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20191020)')
##            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190325)')
#            test1 = EC_index.loc[EC_index.PAR_file.str.contains('H2O2 Daten_BA')]
##            test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.PAR_file.str.contains('PTA1 KOH')]
##            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190912) & (pH < 7) & (SampleID == "JOS2")')
#        elif 'highload' in run:
#            test1 = EC_index.loc[EC_index.PAR_file.str.contains('high-load_267')].query('PAR_exp != "EIS" & EXP_date == 20190912')
##            'O2_ORR_JOS12_3rpm_258_#2_Disc_Parstat'
##            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20190101)')
#        elif 'serie' in run:
##            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20191020)')
##            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190325)')
#            serm = [i for i in run.split() if i in FileHelper.PostChar.SampleSelection.Series.index.values]
#            test1 = EC_index.query('(EXP_date >= 20190101) & (pH < 7) & (SampleID == "DW28")')
##            test1 = OnlyRecentMissingOVV.query('(EXP_date <= 20190228) & (EXP_date >= 20190101)')
##            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20180101) & (EXP_date <= 20180131)')
##            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20180123) & SampleID_folder == "DW28"')
##            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190125)')
##            PAR_date < pd.to_datetime('20190901') and PAR_date > pd.to_datetime('20190827')
#        elif 'missing'in run:
##            EIS_missing = FileHelper.FindExpFolder('VERSASTAT').EISmissing
##            EIS_missing = pd.read_excel(FileHelper.FindExpFolder('VERSASTAT').PostDir.joinpath('OVV_EIS_missing.xlsx'))
##            test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.EXP_date.isin(EIS_missing.EXP_date.unique())]
##            refits = 'DW16','DW19','DW17','DW28'
##            test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.SampleID.isin(refits) & OnlyRecentMissingOVV.pH == 1]
#             test1 = EC_index.loc[EC_index.basename.str.contains('SDFe2AL',case=False)].query('pH < 7' )
##            test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.SampleID.str.contains('SD')]
#        elif 'eisrpm' in run:
#            eisrpm_miss = pd.read_pickle(FileHelper.FindExpFolder('VERSASTAT').PostDir.joinpath('EIS_RPM_series.pkl.compress'))
#            test1 = EC_index.loc[EC_index.PAR_file.isin(eisrpm_miss.PAR_file.to_list())]
#        elif 'orrmiss' in run:
#            test1 = pd.read_pickle(PostDestDir.joinpath('ORR_missing.pkl.compress'))
#        elif 'eismiss' in run:
#            eis_refit = pd.read_pickle(PostDestDir.joinpath('EIS_ORR_refit_pars.pkl.compress')).PAR_file.unique()
#            test1 = EC_index.loc[EC_index.PAR_file.isin(eis_refit)].tail(1)
#        elif 'porph_refit' in run:
#            eis_metaf = pd.read_excel(list(FindExpFolder('PorphSiO2').compare.parent.rglob('EIS_Porph_SiO2\meta_data*EIS*origin.xlsx'))[0],index_col=[0])
#            test1 = EC_index.loc[EC_index.PAR_file.isin(eis_metaf.PAR_file.unique())]
#
#        elif 'now' in run:
#            today = pd.datetime.now().strftime('%Y%m%d')
#            print('run Today only')
#            test1 = EC_index.query('(EXP_date == @today)')
#            if test1.empty:
#                last_day = EC_index.sort_values(by='EXP_date').EXP_date.unique()[-1]
#                test1 = EC_index.query('(EXP_date == @last_day)')
##            & (Loading_name == "high") & (SampleID == "DW21")' )
##            test1 = OnlyRecentMissingOVV[OnlyRecentMissingOVV.EXP_date.isin(OnlyRecentMissingOVV.EXP_date.unique()[-3::])]
##            OnlyRecentMissingOVV.query('(EXP_date > 20190715)').tail(30)
##            test1 = OnlyRecentMissingOVV.query('(EXP_date > 20180101) & (EXP_date < 20181212)')
#    #        & (Gas == "O2")')
#    #        test1,OnlyRecentMissingOVV = MainEC().EC_Analysis_Input()
#        ExpTypes = test1.PAR_exp.unique()
#        RunTypes = [i for i in run.split() if i in ExpTypes]
#        Samples = [i for i in run.split() if i in test1.SampleID.unique()]
#        RunSeries = [i for i in run.split() if i in FileHelper.PostChar.SampleSelection.Series.index.values]
#        if RunTypes:
#            test1 = test1.loc[test1.PAR_exp.str.contains('|'.join(RunTypes))]
#        if Samples:
#            test1 = test1.loc[test1.SampleID.str.contains('|'.join(Samples))]
#        if RunSeries:
#            test1 = test1.loc[test1.SampleID.isin(FileHelper.PostChar.SampleSelection.Series.loc[RunSeries,'sIDs'].values[0])]
#        if 'short' in run:
#            test1 = test1.iloc[0:2]
#        MultiRun = [i for i in run.split() if 'multi' in i]
#        multi_set = True if MultiRun else False
#
#        if 'OER testing' in run:
#            OER_testing = pd.concat([(gr) for n,gr in test1.groupby('PAR_date') if len(gr) > 2]).query('(EXP_date == 20190308)')
#            test1 = OER_testing
#
#    elif 'sr' in run:
#        test1,EC_index = prepare_input.MainEC().EC_Analysis_Input()
#        serie = FileHelper.PostChar.SampleSelection('CB4','SeriesID',[])
#        serie_samples_ovv = pd.merge(serie.Prep_EA_BET.SampleID, EC_index, how='left', on='SampleID')
#        EC_run_Serie_OVV = pd.merge(serie_samples_ovv.EXP_dir.drop_duplicates(), EC_index, how='left', on='EXP_dir')
#
##        ECRunOVV.EC_Analysis_Run(EC_run_Serie_OVV,EC_index)
#    elif re.match('index(?!(force|.force))',run):
#        EC_index = prepare_input.MainEC.EC_Analysis_Input(TakeRecentList = False, exp_folder = None, force_recalc = 'index')
#    elif re.match('index(?=(force|.force))',run):
#        if 'folder' in run:
#            exp_folder = '25.07.2019_0.5MH2SO4_LM'
#            EC_index = prepare_input.MainEC.EC_Analysis_Input(False,exp_folder)
#        else:
#            EC_index = prepare_input.MainEC.EC_Analysis_Input(False,force_recalc='force')
#    return test1.drop_duplicates(subset=['PAR_hash','PAR_file'])
##print('Collection efficiency %s : %.3f' %collection_eff())
##data_filter =
##%%
#
# if __name__ == "__main__":
##    run = input('Want to start the fitting run?')
#    logger = start_logging()
#    run = 'y porph_refit'
#    try:
#        if not EC_index.empty:
#            pass
#        else:
#            EC_index = ECRunOVV(load=1).index.sort_values('PAR_date',ascending=False)
#    except Exception as e:
#        print('reload EC index')
#        EC_index = ECRunOVV(load=0).index.sort_values('PAR_date',ascending=False)
#
#    EC_run_slice = index_run_selection(run,EC_index)
#    ECRunOVV.EC_Analysis_Run(EC_run_slice,EC_index)

#    run = 'ytest ORR N2_act'
# run = 'ytest short EIS'
#        try:
#            MultiMain(
#        except:
#        pass
# from FolderOrganizing.FindExpFolder import *
# ExpDirs = FindExpFolder('VERSASTAT').EC_find_folders()[0]
# f=FindExpFolder('VERSASTAT').EC_find_folders()[0][-1][1][1]
# FindSampleID(FindExpFolder('VERSASTAT').EC_find_folders()[0][-1][1][1]).try_find_sampleID(f)
# [i[1] for i in FindExpFolder('VERSASTAT').EC_find_folders()[0]]
# RunEC(ExpDirs[12::])
# Jkin_gr1500 = Jkin1500_ovv.groupby(by=['SampleID','DATE'])
# for gr in Jkin_gr1500:
#            for name,O2 in O2gr:
#                O2.loc[:,'jcorr'] = O2['j A/cm2'].values - N2_scan['j A/cm2'].values
#            O2gr['j A/cm2']
#            O2_activity = pd.concat([O2_activity,O2_act])
#    O2_activity.to_csv(dest_plot_folder+'\\'+'O2_OVERVIEW.csv')
#  O2_out_ovv.to_csv(dest_plot_folder+'\\'+'ORR_1500_OVERVIEW_v6.csv')
# O2C = pd.read_excel(dest_plot_folder+'\\'+'ORR_JKIN_calc_OVV_compl_V2.xlsx')
# O2C = pd.read_excel(dest_plot_folder+'\\'+'ORR_JKIN_calc_OVV_compl_V3.xlsx')
#%%
