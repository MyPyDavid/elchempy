# @staticmethod

import sys
from pathlib import Path
from collections import namedtuple
from datetime import datetime
import numpy as np
import re
import os
import multiprocessing
from functools import partial
from itertools import repeat
import pandas as pd

from file_py_helper.find_folders import FindExpFolder

# from file_py_helper.file_functions import FileOperations

from .EC_DataLoader.CreateCV import create_CVs as create_CVs
from .N2.analyze_scans import N2_analyze_scan, N2_scans
from .run_Baserunner import base_Analyzer

if __name__ == "__main__":
    from run_experiment import BaseRunner

    pass

import logging

logger = logging.getLogger(__name__)


N2_meta = namedtuple("N2_meta", "file_ovv PAR_file")


def add_missing_create_CVs(
    N2act_ovv, index_fls={}, filter_out=["_1cls_", "_AST_20.000"]
):

    _df = index_fls["N2_CVs"]["data"]
    set_pf_ovv, set_df = set([Path(i) for i in N2act_ovv.PAR_file.unique()]), set(
        _df.PAR_file.unique()
    )
    if set_pf_ovv <= set_df:
        return index_fls
    else:
        _leftover = set_pf_ovv.difference(set_df)
        _leftfltrd = [
            i for i in _leftover if not re.search("|".join(filter_out), i.name)
        ]
        if _leftfltrd:
            _N2left = N2act_ovv.loc[
                N2act_ovv.PAR_file.isin([str(i) for i in _leftfltrd])
            ]
            N2_CVs_all, N2_actions_all = create_CVs(_N2left)
            _CVs_add = pd.concat([_df, N2_CVs_all], sort=False, ignore_index=True)
            _actions_add = pd.concat(
                [index_fls["N2_actions"]["data"], N2_actions_all],
                sort=False,
                ignore_index=True,
            )

            for k, val in index_fls.items():
                #                _df_add = pd.concat([_df,N2_CVs_all], sort=False,ignore_index=True)
                if "N2_actions" == k:
                    index_fls[k].update(
                        {
                            "data": _actions_add,
                            "PAR_files": set(_actions_add.PAR_file.unique()),
                        }
                    )
                    _actions_add.to_pickle(val["path"])
                elif "N2_Cvs" == k:
                    index_fls[k].update(
                        {"data": _CVs_add, "PAR_files": set(_CVs_add.PAR_file.unique())}
                    )
                    _CVs_add.to_pickle(val["path"])
                else:
                    print("read_create_CVs error")
            #                if val['path'].exists():
            return index_fls

        else:
            return index_fls


def eis_multi_wrap(fn, args, kwargs):
    return fn(args[0], args[1], args[-1])


def PF_fit_starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(PF_fit_apply_args_and_kwargs, args_for_starmap)


def PF_fit_apply_args_and_kwargs(fn, args, kwargs):
    return fn(args, **kwargs)


# logger = start_logging(__name__)
def CV_index(N2act_ovv, index_files={}):
    #    N2_CVs.to_pickle(FindExpFolder('VERSASTAT').PostDir.joinpath('N2_CVs_all.pkl'))
    #    N2_actions.to_pickle(FindExpFolder('VERSASTAT').PostDir.joinpath('N2_CVs_actions.pkl'))
    #    N2_CVs_all, N2_actions_all = create_CVs(N2act_ovv)

    for k, val in index_files.items():
        if val["path"].exists():
            _df = pd.read_pickle(val["path"])
            #            _df['PAR_file'] = [Path(i) for i in _df['PAR_file'].to_numpy()]
            #            EIS_ovv['PAR_file'] = [Path(i) for i in EIS_ovv['PAR_file'].to_numpy()]
            index_files[k].update(
                {"data": _df, "PAR_files": set(_df.PAR_file.unique())}
            )
    return index_files


#    set([Path(i) for i in N2act_ovv.PAR_file.unique()]) -   set(_df.PAR_file.unique())


def load_all_N2_CVs(N2act_ovv):
    CV_data_local = CV_index(N2act_ovv)
    CV_data_local = add_missing_create_CVs(
        N2act_ovv, index_fls=CV_data_local, filter_out=["_1cls_", "_AST_20.000"]
    )
    N2_CVs_local = CV_data_local["N2_CVs"]["data"]
    N2_CVs = N2_CVs_local.loc[N2_CVs_local.ActionId == 38]
    N2_actions = CV_data_local["N2_actions"]["data"]
    return N2_CVs.groupby("PAR_file"), N2_actions("PAR_file")


#    dest_files = []
#    for n,gr in N2act_ovv.groupby('PAR_file'):
#        N2_dest_dir = Path(gr.Dest_dir.unique()).joinpath('N2_scans')
#        N2_dest_dir.mkdir(parents=True, exist_ok=True)
#        dest_files.append({'PAR_file' : Path(n),'N2_dest_dir' : N2_dest_dir,
#         'N2_dest_Pars' : EIS_dest_dir.joinpath( Path(r.PAR_file).stem + '_pars.xlsx'),
#         'EIS_dest_spectra' :EIS_dest_dir.joinpath( Path(r.PAR_file).stem + '_Combined.xlsx')
#         })
#    EIS_ovv_destfiles = pd.DataFrame(dest_files).set_index('index')


class Analyze(base_Analyzer):  # TODO develop this class further
    """
    Main class for N2 experiment analysis
    """

    def __init__(
        self,
        index_slice: pd.DataFrame = pd.DataFrame(),
        multi_par_fit=False,
        **N2_kwargs,
    ):

        self.index_slice = index_slice
        self.multi_par_fit

    def run_parallel(self):
        pass


def N2_act(ovv_exp_grp, **N2_kwargs):
    ###### === Analyze the N2 experiments ==== #######
    """Filter out relevant N2 scans from All_ovv data and calculate Capacity values"""
    #            N2bgs = pd.DataFrame(),  N2bgs = pd.concat([N2_background,N2bgs])
    #        if 'N2_act' in Experiments: N2_bgs_lst = []
    index_out = []
    #                    All_N2 = N2_CVs.loc[(N2_CVs.EXP == "N2_act") & (~N2_CVs.File.astype(str).str.contains('Pt-ring')),:]
    #                    Cdl_fit_pars, Cdl_fit_data = pd.DataFrame(), pd.DataFrame()
    N2act_ovv = ovv_exp_grp.get_group("N2_act").drop_duplicates(
        subset=["PAR_file", "PAR_hash"], keep="first"
    )
    fit_run_args = [N2_meta(gr, pf) for pf, gr in N2act_ovv.groupby("PAR_file")]
    #    N2_skip = ['08-Aug-2017\\24.08.2017_DW05_0.1M-H2SO4\\N2_300-20cls.par' ,'PTA3\RingO2_vs_Ag_AgCl_1.par']
    #    if N2_file in N2_skip:
    #            logger.warning('Create CV N2 skipped because in skip list N2_skip')
    if len(N2act_ovv) > 1e3:
        load_all_N2_CVs(N2act_ovv)
    #    ovv.query('(PAR_exp == "N2_act")')
    #        N2act_ovv = ExpTypes_gr.get_group('N2_act').query('(PAR_exp == "N2_act")')
    multi_par_fit = N2_kwargs.get("multi_par_fit", True)
    if multi_par_fit:
        fit_export_EV_all = []
        try:
            #            EISgrEvRHE_data_grp = EISgrEvRHE_data_raw.groupby(['PAR_file',EvRHE, 'RPM_DAC'])
            fit_kwargs_iter = repeat(N2_kwargs)
            #            args_for_starmap = zip(repeat(eis_fit_PAR_file), t_run_args, t_kwargs_iter)
            #            EIS_run_grps = [(nm,gr,EIS_kwargs) for nm,gr in EIS_ovv.groupby(by='PAR_file')]
            logger.info(
                f"N2 run_group_ovv START multiprocessing {multi_par_fit} for len{len(fit_run_args)}"
            )
            pool_size = os.cpu_count() - 2
            if "linux" in sys.platform:
                #                os.system('taskset -cp 0-%d %s' % (pool_size, os.getpid()))
                os.system("taskset -p 0xff %d" % os.getpid())
            #                os.sched_setaffinity(0,{i for i in range(pool_size)})
            #            for chunk in fit_run_args[:pool_size if pool_size > 2 else 1]:
            with multiprocessing.Pool(pool_size) as pool:
                fit_export_EV_all_chunck = PF_fit_starmap_with_kwargs(
                    pool, N2_scans, fit_run_args, fit_kwargs_iter
                )
                fit_export_EV_all.append(fit_export_EV_all_chunck)
        #                results = pool.map(partial(eis_fit_PAR_file, 'This is arg1', **just_a_dict), list_of_args2)
        #                print('eis_run_group_ovv multiprocessing error: {0}'.format(e))
        #                logger.error('eis_run_group_ovv  multiprocessing error: {0}'.format(e))
        except Exception as e2:
            #        print('EIS eis_run_group_ovv multiprocessing error: {0}'.format(e2))
            logger.error(
                f"EIS eis_fit_PAR_file  multiprocessing erroe: {e2}, len out({len(fit_export_EV_all)})"
            )
    #            multi_par_fit = False
    if multi_par_fit == False:
        #                dest_dir = N2grF.Dest_dir.values[0]
        #        for N2_file, N2_ovv_file in N2act_ovv.groupby('PAR_file'):
        for fit_run_arg in fit_run_args:
            #            fit_run_arg = fit_run_args[-1]
            try:
                N2_background, Cdl_dt, Cdl_PARS = N2_scans(fit_run_arg, **N2_kwargs)
                if not N2_background.empty:
                    #                        N2_bgs_lst.append(N2_background)
                    logger.info(
                        f"N2 fit_run_arg background succes: {fit_run_arg.PAR_file}"
                    )
                else:
                    logger.warning(
                        f"N2 fit_run_arg background empty for: {fit_run_arg.PAR_file}"
                    )
            #                        Cdl_fit_pars = pd.concat([Cdl_fit_pars,Cdl_PARS])
            #                        Cdl_fit_data = pd.concat([Cdl_fit_data,Cdl_dt])
            #                    print('N2 background succes!')
            except Exception as e:
                #                print('N2 FAIL',e)
                logger.error(
                    f"ERROR N2 run_fit_arg FAIL: {fit_run_arg.PAR_file}.\nbecause  {e}"
                )
                N2_background, Cdl_dt, Cdl_PARS = (
                    pd.DataFrame([]),
                    pd.DataFrame([]),
                    pd.DataFrame([]),
                )


#            N2_BGs = pd.concat([i for i in N2_bgs_lst])
#        indexes = pd.concat([pd.DataFrame(i) for i in index_out],sort=False).reset_index(drop=True)
#    return indexes
#            print('=== no N2 in Experiments === EMPTY')
