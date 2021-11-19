# @staticmethod

import sys
from pathlib import Path
from collections import namedtuple
from dataclasses import dataclass, field
from typing import Dict, List

import logging

import datetime as dt
import numpy as np
import re
import os


import multiprocessing
from functools import partial
from itertools import repeat
import pandas as pd

import matplotlib.pyplot as plt

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations

from .EC_DataLoader.CreateCV import create_CVs
from .EC_DataLoader.set_OCP_RHE import get_RHE_OCP

from .ORR.ORR_analyze_scans import ORR_calculations, ORR_KL_loop

import logging

_logger = logging.getLogger(__name__)

Meta = namedtuple("Meta", "PAR_file data ovv N2_BG_f N2_BG_data")

# class ORR_scan_data():


@dataclass(order=True, frozen=False)
class ORR_scan_data:
    """ORR scan dataclass.\n
    Holds scan raw data in pd.DataFrame and metadata"""

    # _required_cols = ['Frequency(Hz)','Z Real','Z Imag']
    # _spectrum_grp_cols = ['PAR_file','Segment #',EvRHE, 'RPM_DAC']
    N2_scans_folder_version = f"N2_scans_v{FileOperations.version}"
    N2_BG_ext = ".pkl"
    # pf, grp,self.ovv_all
    PAR_file_disk: Path = field(default=Path(Path.cwd().joinpath("empty.txt")))
    ovv_disk: type(pd.DataFrame) = field(default=pd.DataFrame(), repr=False)
    ovv_all: type(pd.DataFrame) = field(default=pd.DataFrame(), repr=False)
    EC_index: type(pd.DataFrame) = field(default=pd.DataFrame(), repr=False)
    run_kwargs: Dict = field(default_factory=dict, repr=False)
    # E_dc_RHE: float = 0.
    # RPM_DAC : float = 0.

    def __post_init__(self):
        if type(self.PAR_file_disk) == str:
            self.PAR_file_disk = Path(self.PAR_file_disk)
        self._check_N2()
        self.get_N2_options()
        self.load_N2_data()

    # data : type(pd.DataFrame) = field(default=pd.DataFrame(),repr=False)
    # ovv : type(pd.DataFrame) = field(default=pd.DataFrame(),repr=False)

    def _check_N2(self):
        self.ORR_get_N2_BG(
            self.ovv_all, self.PAR_file_disk, self.ovv_disk, **self.run_kwargs
        )

    ### N2 testing ###
    # grps = list(gr_ovv_disk.groupby(by='PAR_file').groups)
    # ORR_file, ORR_ovv_file = grps[0], gr_ovv_disk.groupby(by='PAR_file').get_group(grps[0])

    def get_N2_options(self):
        _opts = []
        for n, r in self.ovv_disk.iterrows():
            try:
                N2exMatch_options = self.EC_index.loc[
                    (
                        (self.EC_index["PAR_exp"] == "N2_act")
                        & (self.EC_index.SampleID == r.SampleID)
                        & (self.EC_index.postAST == r.postAST)
                        & (self.EC_index.Electrolyte == r.Electrolyte)
                        # & (self.EC_index.PAR_date_day == r.PAR_date_day)
                        & (self.EC_index.Loading_name == r.Loading_name)
                    )
                ]
            except:
                N2exMatch_options = pd.DataFrame
            if not N2exMatch_options.empty:
                N2exMatch_options = N2exMatch_options.assign(
                    **{
                        "N2_PAR_date_diff": np.datetime64(r.PAR_date)
                        - N2exMatch_options.PAR_date,
                        "N2_PAR_date_diff_seconds": [
                            i.total_seconds()
                            for i in (
                                np.datetime64(r.PAR_date) - N2exMatch_options.PAR_date
                            )
                        ],
                    }
                )
                N2exMatch_options = N2exMatch_options.sort_values(
                    by="N2_PAR_date_diff_seconds", ascending=False
                )

                N2_ext_BG_opts = [
                    Path(r.Dest_dir).joinpath(
                        self.N2_scans_folder_version
                        + "/"
                        + r.basename
                        + f"_BG{self.N2_BG_ext}"
                    )
                    for n, r in N2exMatch_options.iterrows()
                ]
                N2exMatch_options = N2exMatch_options.assign(
                    **{"ORR_act_N2_bg": N2_ext_BG_opts}
                )

                # N2_scan_dd = Path(r.Dest_dir).joinpath(self.N2_scans_folder_version)
                # N2_ext_BG = N2_scan_dd.joinpath()
                # N2_ext_BG = Path(run_ORR.N2_BG_exclusion_list(N2_ext_BG))
                _opts.append(N2exMatch_options)
        if _opts:
            _N2_option_ovv = pd.concat(_opts)
            _N2_option_ovv = get_RHE_OCP(_N2_option_ovv)[0]
            # gr_ovv = RHE_index[0]
            self.N2_option_ovv = _N2_option_ovv

    def load_N2_data(self):
        if hasattr(self, "N2_option_ovv"):
            _N2coll = []
            for n, r in self.N2_option_ovv.iterrows():
                _N2read = pd.DataFrame()
                if r.ORR_act_N2_bg.is_file():

                    if "pkl" in self.N2_BG_ext:
                        _N2read = pd.read_pickle(r.ORR_act_N2_bg)
                    elif "xl" in self.N2_BG_ext:
                        _N2read = pd.read_excel(r.ORR_act_N2_bg)
                else:
                    _N2read = self.read_out_N2_data(r)

                if not _N2read.empty:
                    _N2read = _N2read.assign(
                        **{
                            "ORR_PAR_file": self.PAR_file_disk,
                            "N2_PAR_date_diff": r.N2_PAR_date_diff,
                            "N2_PAR_date_diff_seconds": r.N2_PAR_date_diff_seconds,
                            "N2_EC_index": n,
                        }
                    )
                    _N2coll.append(_N2read)

            if _N2coll:
                self.N2_BG_data = pd.concat(_N2coll, ignore_index=True)

    def read_out_N2_data(self, row):
        # TODO move
        _raw_data, _raw_actions = create_CVs(row.to_frame().T)

        _10mVs_data = _raw_data.loc[_raw_data["Scan Rate (V/s)"] == 0.01]
        if 1800 < len(_10mVs_data) < 2200:
            return _10mVs_data
        else:
            return pd.DataFrame()

    def ORR_get_N2_BG(self, ovv_all, ORR_file, ORR_ovv_file, **ORR_kwargs):
        N2_bg_file_used = ""
        N2_scans_folder_version = f"N2_scans_v{FileOperations.version}"
        N2_scan_dd = Path(ORR_ovv_file.iloc[0].Dest_dir).joinpath(
            N2_scans_folder_version
        )
        N2_bg_local = pd.DataFrame()
        N2_BG_ext = ORR_kwargs.get("N2_BG_ext", ".pkl")

        try:
            #        gr_ovv.loc[:,'ORR_act_N2_bg'] =
            #        gr_ovv.astype(['ORR_act_N2_bg'].astype(str)
            if len(ORR_ovv_file) == 1:
                ovv_N2_bg_file = str(ORR_ovv_file.ORR_act_N2_bg.values[0])
            #                                ovv_N2_bg_file  = Path([i.ORR_act_N2_bg for n,i in ovv.iterrows() if str(i.PAR_file) == str(ORR_file)][0])
            if not Path(ovv_N2_bg_file).exists():
                glob_N2file_opts = sorted(
                    [
                        (i.stem, i, i.stat().st_ctime)
                        for i in Path(ORR_ovv_file.iloc[0].Dest_dir).parent.rglob(
                            f"*{N2_scans_folder_version}/*_BG{N2_BG_ext}"
                        )
                        if i.stem == Path(ovv_N2_bg_file).stem
                    ],
                    key=lambda x: x[2],
                )
                if glob_N2file_opts:
                    ovv_N2_bg_file = glob_N2file_opts[-1][1]
                else:
                    _logger.warning(
                        "ORR ERROR ovv_N2_bg_fail file not found with rglob"
                    )
                if not ovv_N2_bg_file:
                    try:
                        ovv_N2_bg_file = (
                            ovv_all.query(
                                'SampleID == @ORR_ovv_file.SampleID.unique()[0] & PAR_exp == "N2_act"'
                            )
                            .sort_values("Date_PAR_EXP", ascending=False)
                            .head(1)
                            .ORR_act_N2_bg.values[0]
                        )
                    except Exception as e2:
                        _logger.warning(
                            "ORR ERROR ovv_N2_bg_fail file not found: {0}".format(e2)
                        )
            #        N2_BGs = pd.read_excel(ovv_N2_bg_file, index_col=[0])
            ovv_N2_bg_file = N2_BG_exclusion_list(ovv_N2_bg_file)
            N2_BGs = pd.read_pickle(ovv_N2_bg_file)
            N2_bg_file_used = ovv_N2_bg_file
            _logger.info("ORR ovv_N2_bg_succes: {0}".format(ovv_N2_bg_file))
        #                            Path(ovv.loc[ovv['PAR_file'] == str(ORR_file)].ORR_act_N2_bg.values[0])
        except Exception as e:
            #                            print('ERROR ovv_N2_bg_fail: {0}'.format(e))
            N2_BGs = pd.DataFrame()
            _logger.warning("ORR ERROR ovv_N2_bg_fail EMPTY: {0}".format(e))

        if not ORR_ovv_file.empty:

            if not N2_BGs.empty:
                try:
                    N2_bg_local = N2_BGs.loc[
                        N2_BGs.basename == Path(ovv_N2_bg_file).stem.split("_BG")[0]
                    ]
                except Exception as e:
                    _logger.warning(f"ORR ERROR ORR_get_N2_BG: {e}")
                    N2_bg_local = pd.DataFrame()

            if N2_bg_local.empty:
                ORR_ovv_all_slice = ovv_all.loc[ovv_all["PAR_file"] == Path(ORR_file)]
                if not ORR_ovv_all_slice.empty:
                    try:
                        N2_bg_local_file_try = ORR_ovv_all_slice.ORR_act_N2_bg.values[0]
                        N2_bg_local = pd.read_excel(N2_bg_local_file_try)
                        N2_bg_file_used = N2_bg_local_file_try
                    except:
                        N2_ext_match = []

                        N2_ext_match = list(
                            Path(ORR_ovv_all_slice.Dest_dir.values[0]).parent.rglob(
                                f"*{N2_scans_folder_version}/*{Path(ovv_N2_bg_file).stem}*BG"
                            )
                        )
                        if N2_ext_match:
                            N2_bg_local = pd.read_excel(N2_ext_match)
                            N2_bg_file_used = N2_ext_match
                else:
                    pass

        return N2_bg_file_used, N2_bg_local


def corrected_filenames():
    """/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT/2019-05-May
    /06.05.2019_0.1MH2SO4_cell3/O2_ORR_JOS3_pAST-sHA_285_#3_Disc_Parstat.par
    changed: N2_post_AST_JOS3
    """


class ORR_collection:
    def __init__(self, _obj):
        self.ORR_run_loop = _obj
        self.run_kwargs = self.ORR_run_loop.run_kwargs
        self.EC_index = self.ORR_run_loop.runOVV.EC_index
        self.ORR_prepare_ovv_dest()
        self.prep_collection()

    def ORR_prepare_ovv_dest(self):
        # ovv_exp_grp = self.ORR_run_loop.ovv_exp_grp,
        gr_ovv = self.ORR_run_loop.ovv_exp_grp.get_group("ORR")
        RHE_index = get_RHE_OCP(gr_ovv)
        gr_ovv = RHE_index[0]

        gr_ovv = gr_ovv.astype({"ORR_act_N2_bg": str})

        gr_ovv = gr_ovv.assign(
            **{
                "ORR_dest_dir": [
                    Path(i).joinpath(f"ORR_v{FileOperations.version}")
                    for i in gr_ovv.Dest_dir.to_numpy()
                ]
            }
        )
        #    for _dest_dir in gr_ovv.ORR_dest_dir.unique()
        _mkdirs = [
            i.mkdir(parents=True, exist_ok=True) for i in gr_ovv.ORR_dest_dir.unique()
        ]
        gr_ovv_disk = gr_ovv.loc[
            (gr_ovv._act0_Comment.str.contains("Pt|chrono") == False)
            & (gr_ovv.SampleID.str.contains("PT-RING") == False)
        ]

        _calc_dd = []
        for n, r in gr_ovv_disk.iterrows():
            _dd = r.ORR_dest_dir
            _ORR_calc_dest_dir = _dd.joinpath(
                f'{r["Electrolyte"]}_{r["SampleID"]}_{r["Loading_name"]}_{r["postAST"]}'
            )
            _ORR_calc_dest_dir.mkdir(parents=True, exist_ok=True)
            _fstem = FileOperations.Check_fstem(Path(r.PAR_file))
            _calc_dd.append((_ORR_calc_dest_dir, _fstem))

        gr_ovv_disk = gr_ovv_disk.assign(
            **{
                "ORR_ecexp_destdir": [i[0] for i in _calc_dd],
                "PAR_fstem": [i[1] for i in _calc_dd],
            }
        )

        self.ovv_all = gr_ovv
        self.gr_ovv_disk = gr_ovv_disk

    def prep_collection(self):
        # ovv_all, gr_ovv_disk = ORR_prepare_ovv_dest(self.ORR_run_loop.ovv_exp_grp, **self.ORR_run_loop.run_kwargs)
        # [Meta(Path(pf),grp,ovv_all.loc[ovv_all.PAR_date_day.isin(grp.PAR_date_day.unique())],*ORR_get_N2_BG(ovv_all, pf, grp, **self.run_kwargs))
        # for pf,grp in gr_ovv_disk.groupby(by='PAR_file')]
        _dict = {}
        for pf, grp in self.gr_ovv_disk.groupby(by="PAR_file"):
            pf, grp
            _orr = ORR_scan_data(pf, grp, self.ovv_all, self.EC_index, self.run_kwargs)
            _dict.update({pf: _orr})

        orr_run_args_raw = _dict
        # Path(pf),grp,ovv_all.loc[ovv_all.PAR_date_day.isin(grp.PAR_date_day.unique())],
        self.orr_run_args_raw = orr_run_args_raw

    def __iter__(self):
        return iter(self.orr_run_args_raw.values())

    def __len__(self):
        return len(self.orr_run_args_raw.values())

    def slice_N2_from_index(self, grp):
        _grpdict = grp.iloc[0].to_dict()
        _grpdict.get("SampleID")

        grp.SampleID.unique()
        self.EC_index

    #    EC_index = ORR_kwargs.get('EC_index', pd.DataFrame())
    #    if not EC_index.empty:
    #        _add_ring_ovv = EC_index.loc[EC_index.PAR_date.isin(gr_ovv_disk.PAR_date.unique())]
    #        gr_ovv
    # return gr_ovv, gr_ovv_disk


#    for ORR_file, ORR_gr_ovv in gr_ovv_disk.groupby(by='PAR_file'):
#        DestDir = Path(ORR_gr_ovv.Dest_dir.unique()[0])
#        ORR_dest_dir_gr = DestDir.joinpath(f'ORR_v{FileOperations.version}')
#        ORR_dest_dir_gr.mkdir(parents=True, exist_ok=True)


def eis_multi_wrap(fn, args, kwargs):
    return fn(args[0], args[1], args[-1])


def PF_fit_starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(PF_fit_apply_args_and_kwargs, args_for_starmap)


def PF_fit_apply_args_and_kwargs(fn, args, kwargs):
    return fn(args, **kwargs)


def PF_fit_starmap_with(pool, fn, args_iter):
    args_for_starmap = zip(repeat(fn), args_iter)
    return pool.starmap(PF_fit_apply_args_and_kwargs, args_for_starmap)


def PF_fit_apply_args_and_kwargs(fn, args):
    return fn(args)


class ORR_run_loop:
    def __init__(self, _obj, testing_mode=False):
        if "ECRunOVV" in _obj.__class__.__name__:
            self.runOVV = _obj
            self.run_kwargs = self.runOVV.run_kwargs
            self.ovv_exp_grp = self.runOVV.ovv_exp_grp
            self.run_prep()
            _logger.warning(
                f"{self.__class__.__name__} starting with {len(self.ovv_exp_grp.groups)/2}"
            )
        else:
            raise TypeError

        if (
            not "test" in self.runOVV.input_run
            and self.runOVV.input_run.endswith("y")
            and not testing_mode
        ):
            self.start_run()
        else:
            _logger.error(
                f'ORR run loop finished because {self.runOVV.input_run} has "test" or not endswith "y"'
            )

    #        self._gen_loop = _gen

    def run_prep(self):

        self.ORR_collection = ORR_collection(self)
        # self.orr_run_args = iter(self.ORR_collection)
        # ovv_all, gr_ovv_disk = ORR_prepare_ovv_dest(self.ovv_exp_grp, **self.run_kwargs)
        # orr_run_args_raw = [Meta(Path(pf),grp,ovv_all.loc[ovv_all.PAR_date_day.isin(grp.PAR_date_day.unique())],*ORR_get_N2_BG(ovv_all, pf, grp, **self.run_kwargs))
        #                     for pf,grp in gr_ovv_disk.groupby(by='PAR_file')]
        # _n2faillst = [i[0] for i in orr_run_args_raw if i[-1].empty]
        # orr_run_args = [i for i in orr_run_args_raw if not i[-1].empty]
        # self._dt_start = dt.datetime.now()

    def start_run(self, **ORR_kwargs):
        #%%
        ###### === Analyze the ORR experiments ==== #######
        #        exp,gr_ovv,ovv = 'ORR',ExpTypes_gr.get_group('ORR'),ovv
        #        if 'O2' in Gases and not 'N2_act' in Experiments:
        #            print('N2 BACKGROUND SCAN IS MISSING FOR THIS ORR EXPERIMENT.\nFailed: %s'%exp_dir)
        #            sORR = 1
        #        elif 'O2' in Gases and 'N2_act' in Experiments: # O2 and N2 changed
        #        O2_activity, O2_out_ovv, O2_Jcalc_ovv, overall_rows   = pd.DataFrame([]),pd.DataFrame([]),pd.DataFrame([]), pd.DataFrame([])
        index_ORR, index_out, faillst = [], [], []

        multi_par_fit = True
        _t = ORR_kwargs.get("multi_par_fit", True)
        # if gr_ovv_disk.empty:
        # return _logger.warning('ORR attemp failed for {0} because ovv empty'.format(gr_ovv.Dest_dir.unique()))
        #    _logger.warning('ORR disk OVV empty')
        _dt_start = dt.datetime.now()
        if multi_par_fit:
            fit_export_EV_all = []
            # orr_kwargs_iter = repeat(self.ORR_kwargs)
            try:
                _logger.info(
                    f"ORR orr_run_group_ovv START multiprocessing {multi_par_fit} for len{len(self.ORR_collection.gr_ovv_disk)}"
                )
                pool_size = os.cpu_count() - 2
                if "linux" in sys.platform:
                    #                os.system('taskset -cp 0-%d %s' % (pool_size, os.getpid()))
                    os.system("taskset -p 0xff %d" % os.getpid())
                #                os.sched_setaffinity(0,{i for i in range(pool_size)})
                #            for chunk in orr_run_args[:pool_size if pool_size > 2 else 1]:
                with multiprocessing.Pool(pool_size) as pool:
                    try:
                        PF_fit_starmap_with(
                            pool, self.ORR_calc_arg, iter(self.ORR_collection)
                        )
                    except Exception as e1:
                        _logger.error(f"ORR run multiprocessing pool error: {e1})")
            #                    fit_export_EV_all.append(fit_export_EV_all_chunck)
            except Exception as e2:
                _logger.error(f"ORR run multiprocessing erroe: {e2})")

        else:
            fit_export_EV_all = []
            #       fit_run_arg = orr_run_args[-1]
            #       [i[0] for i in orr_run_args]
            for fit_run_arg in iter(self.ORR_collection):
                try:
                    self.ORR_calc_arg(fit_run_arg)
                    # calc = ORR_calculations(fit_run_arg)
                    # ORR_KL_loop(calc, run_loop = True)
                    # Jkin_calculations(fit_run_arg, **self.ORR_kwargs)
                #                fit_export_EV_all.append([fit_export_EV])
                except Exception as e:
                    _logger.warning(f"ORR attemp failed for {fit_run_arg} because {e}")
        _dt_end = dt.datetime.now()
        _logger.warning(f"ORR run ovv finished in {_dt_end-_dt_start} ")

    @staticmethod
    def ORR_calc_arg(fit_run_arg):
        try:
            calc = ORR_calculations(fit_run_arg)
            ORR_KL_loop(calc, run_loop=True)
        except Exception as e:
            _logger.warning(f"ORR attemp failed for {fit_run_arg} because {e}")


#        for ORR_file, ORR_ovv_file in gr_ovv_disk.groupby(by='PAR_file'):
#            try:
##                N2_bg_file_used, N2_bg_local = ORR_get_N2_BG(gr_ovv,ORR_file, ORR_ovv_file)
##                if N2_bg_local.empty:
##                    _logger.warning('ORR ERROR N2_bg_local empty: {0}'.format(ORR_file))
##                    faillst.append((ORR_file,ORR_gr_ovv,N2_bg_file_used,N2_bg_local))
#                else:
#                    ###### === Run calculcations for ORR ==== #######
#                    index_ORR = Jkin_calculations(ORR_ovv_file, gr_ovv, N2_bg_file_used, N2_bg_local,**ORR_kwargs)
#
#                    index_out.append(index_ORR)
#                _logger.info('ORR attemp succes for {0}'.format(ORR_file))
#            except Exception as e:
#                    #                    print('ERROR ORR attemp failed, %s' %e)
#                _logger.warning('ORR attemp failed for {0} because {1}'.format(ORR_file, e))
#            DestDir = Path(ORR_gr_ovv.Dest_dir.unique()[0])
#            ORR_dest_dir_gr = DestDir.joinpath(f'ORR_v{FileOperations.version}')
#            O2_CVs.query('(Gas == "O2") & (Type_action == "Cyclic Voltammetry (Multiple Cycles)")' ).plot(x=EvRHE,y='j A/cm2',xlim=[0,1.2] )
#                O2_act = O2_act.assign(Sweep_Type=np.where(O2_act[EvRHE].diff() > 0, 'anodic',
#                                                           np.where(O2_act[EvRHE].diff() < 0, 'cathodic', 'NA')))
#            names_overall_row = ['SampleID','File','DATE','E_onset','E_half','J_diff_lim','Jkin_075','Tafel_slope','KL_slope']
#                ovv_ORR = ovv.query('(PAR_exp == "ORR")')
#                if ovv_ORR.empty:
#                    print('No ORR in ovv')
# ======= CALC ====#
###### === Prepare N2 background for ORR ==== #######
#                    if O2_act.empty:
#                        raise ValueError('O2_act empty')
#                    else:
#                    for nORR, grORR in O2_act.groupby(by='PAR_file'):
#                            if not 'O2_ORR_JOS5_257_#3_Disc_Parstat.par' in nORR.parts:
#                                continue
#                            pd.read_excel(ORR_dest_dir_gr.joinpath('ORR_pars_%s.xlsx'%grORR.SampleID.unique()[0]))
#                            nORR.name, [Path(i).name for i in N2bgs.File.unique()]
#                        coll_rows.to_excel(ORR_dest_dir.joinpath('ORR_pars_%s.xlsx'%coll_rows.SampleID.unique()[0]))
#                            if ORR_dest_dir.joinpath('ORR_pars.xlsx').is_file():
#                                pd.concat([coll_rows,pd.read_excel(ORR_dest_dir.joinpath('ORR_pars.xlsx'))])
#                            else:
#                                coll_rows.to_excel(ORR_dest_dir.joinpath('ORR_pars.xlsx'))
#                    print('ORR succes,%s'%exp_dir)
#            O2_out_ovv = pd.concat[O2_out_ovv,O2_out1]
#    indexes = pd.concat([pd.DataFrame(i) for i in index_out], sort=False).reset_index(drop=True)
#    return indexes
