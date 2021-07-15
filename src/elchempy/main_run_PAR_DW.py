# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 15:05:57 2017
@author: David Wallace
This module reads PAR files in a certain folder and tries to analyze the experimental data
from N2, ORR, OER, HER, HPRR and EIS measurements.
"""
# import os
import sys

# sys.setrecursionlimit(3000)

from collections import namedtuple

import numpy as np
import pandas as pd
import re

# import multiprocessing
from pathlib import Path

import datetime as dt

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper import PostChar

if __name__ == "__main__":

    from indexer import prepare_input

    from experiments import run_EIS
    from experiments import run_N2
    from experiments import run_ORR
    from experiments import run_OER
    from experiments import run_HER
    from experiments import run_HPRR

    from experiments.EC_DataLoader.set_OCP_RHE import get_RHE_OCP

    # from PostEC.collect_load import get_EIS_pars
    # from PostEC.EIS_export import EIS_all_check_redchi, export_to_Origin

else:
    from .indexer import prepare_input

import logging

_logger = logging.getLogger(__name__)
# _logger.propagate = False
#%%
# Globally used Variable in every module
globals()["EvRHE"] = "E_AppV_RHE"
#%%


def check_interactive_for_multi():
    if not sys.__stdin__.isatty():
        _multi = False
    else:
        _multi = True
    return _multi


#%% ==== ECRunOVV start =====
class ECRunOVV:
    global EvRHE

    POST_DIR = FindExpFolder("VERSASTAT").PostDir

    """ Run over a DateFrame that is an Overview(OVV) of EC experiment in a folder."""

    def __init__(self, **kwargs):
        self.run_kwargs = kwargs
        self.input_run = self.run_kwargs.get("input_run", "")
        self._start_time = dt.datetime.now()
        self.EC_run_slice = pd.DataFrame()
        if not self.run_kwargs.get("EC_index", pd.DataFrame()).empty:
            self.EC_index = self.run_kwargs.get("EC_index")
        else:
            #        if 'load_index' in self.run_kwargs:
            self.EC_index = prepare_input.MainEC.EC_Analysis_Input(
                self.run_kwargs.get("load", True)
            )

        if self.run_kwargs.get("input_run", ""):
            self.EC_run_slice, self.run_kwargs = index_run_selection(
                self.EC_index, **self.run_kwargs
            )

        if not self.EC_run_slice.empty:
            self.update_run_kwargs()
            self.prepare_ovv_EC_run()
            self.update_type_run_kwargs()
            if (
                not "000" in self.input_run
                and self.input_run.startswith("y")
                and not self.input_run.rstrip().endswith(" n")
            ):
                self.EC_Analysis_run()
            self._end_time = dt.datetime.now()
            self._finish()
        else:
            self._finish()

    #         if not '000' in self.input_run and\
    #                         self.input_run.startswith('y') and not self.input_run.rstrip().endswith(' n'):
    #             self.EC_Analysis_Run(self)
    # #        self.ovv = ovv

    #    PathDB = FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_Files_class.hdf5')
    def __repr__(self):
        _t = f'ECRunOVV, input: "{self.input_run}" from idx({len(self.EC_index):.0f}). '
        if not self.EC_run_slice.empty:
            _t += "\n"
            _t += f"Run slice: ({len(self.EC_run_slice)}),with {self.EC_run_slice.PAR_exp.value_counts().to_dict()}"
        else:
            _t += "Run slice empty"
        return _t

    def _timeit(self):
        if hasattr(self, "_start_time") and hasattr(self, "_end_time"):
            self._duration = self._end_time - self._start_time
        else:
            self._duration = 0

    def _finish(self):
        self._timeit()
        if self.EC_run_slice.empty:
            _logger.warning(
                f"{__name__} finished, input run slice was empty.\nlen(slice:{len(self.EC_run_slice)} and index:{len(self.EC_index)})\n"
            )
            _logger.info(
                f'EIS suggestions(bad): {len(self.run_kwargs.get("EIS_recheck_bad_fits", pd.DataFrame()))}.\nInput command was: {self.run_kwargs.get("input_run")}'
            )
        #            return print('Finished')
        else:
            _logger.warning(
                f'ECrunovv finished in {self._duration}, started at {self._start_time}.\nInput command was: {self.run_kwargs.get("input_run")}'
            )
            pass

    #    @staticmethod
    def MakeOVVperExpDir(self, arg):
        #%% MakeOVVperExpDir
        global EvRHE
        #        PathDB = FileHelper.FindExpFolder('VERSASTAT').DestDir.joinpath('PAR_Files_class.hdf5')
        ovv = pd.DataFrame([])
        #        dest_dir = Path(gr['Dest_dir'].unique()[0])
        #        save_dir = Path(FileHelper.FindExpFolder('VERSASTAT').DestDir,dest_dir.parts[-2],dest_dir.parts[-1])
        # ===== Compare with saved STATUS to skip or not  =====
        """Check status in an HDF5 file"""
        #            if exp_dir.N2 == 1 and exp_dir.HPRR == 1 or exp_dir.ORR == 1:
        #                print('Skipped %s'%basedir)
        #                continue
        #            else:
        #                print('Perform analysis on: %s' %basedir)
        # ===== Set default values =====
        ###### === Start with overview and classification of files in the Experimental Folder ==== #######
        """Filter gr of OVV for failed measurements"""
        gr = arg.grp
        exp_dir = arg.EXP_dir
        #        OnlyRecentMissingOVV = self.EC_index
        ovv = gr.loc[
            (
                gr.PAR_file
                == [
                    i
                    for i in gr.PAR_file.values
                    if not ("template|templ|fail|Kopie|Copy") in str(i)
                ]
            )
        ]
        N2_BG_ext = self.run_kwargs.get("N2_BG_ext", ".pkl")
        #        ovv = gr.loc[(~gr['PAR_file'].str.contains('template|templ|fail')),:]
        #        SampleIDfolder = FileHelper.FindSampleID.match_SampleID(exp_dir)
        if not ovv.Date_PAR_EXP.empty:
            ovv = ovv.loc[
                (
                    (ovv.Date_PAR_EXP > pd.Timedelta("-2 days 0 hours"))
                    | (ovv.PAR_exp == "RHE")
                ),
                :,
            ]
            if ovv.empty:
                _logger.warning("MakeOVVperExpDir files filter by DATE ignored")
                ovv = gr.loc[
                    (
                        gr.PAR_file
                        == [
                            i
                            for i in gr.PAR_file.values
                            if not ("template|templ|fail|Kopie|Copy") in str(i)
                        ]
                    )
                ]
        ovv_nfiles, ovv_nfls = ovv.PAR_hash.nunique(), ovv.PAR_date.nunique()
        singles = ovv.loc[
            ovv.duplicated(subset=["PAR_hash", "PAR_date"], keep=False) == False
        ]
        dups = ovv.loc[
            ovv.duplicated(subset=["PAR_hash", "PAR_date"], keep=False) == True
        ]
        dups_folder_match = dups.loc[dups.sID_Match_file_dir.str.contains("yes")]
        already_found = pd.concat([singles, dups_folder_match])
        left_overs = ovv[~ovv.PAR_hash.isin(already_found.PAR_hash)].drop_duplicates(
            subset="PAR_hash"
        )
        if not left_overs.empty:
            ovv_filtered = pd.concat([already_found, left_overs])
        else:
            ovv_filtered = already_found
        number_filtered = len(ovv) - len(ovv_filtered)
        if ovv_filtered.PAR_hash.nunique() == ovv_nfiles:
            #            _logger.warning('MakeOVVperExpDir files ({0}) were filtered out as duplicates {1}'.format((),exp_dir))
            ovv = ovv_filtered
        #        dups_folder_match = ovv.loc[(ovv.duplicated(subset=['PAR_hash'],keep=False) == True & ovv.sID_Match_file_dir.str.contains('yes')) & (ovv.duplicated(subset=['PAR_hash'],keep=False) == False)]
        #        ovv.loc[(ovv.duplicated(subset=['PAR_hash'],keep=False) == False)]
        #        set(gr.PAR_file.values) - set(ovv.PAR_file.values)
        OVV_GR_len = np.abs(len(gr) - len(ovv))
        if OVV_GR_len > 3:
            _logger.warning(
                "MakeOVVperExpDir Many files ({0}) were filtered out by date or name from {1}".format(
                    OVV_GR_len, exp_dir
                )
            )
        SkippedLst = []
        if ovv.empty:
            _logger.warning(f"MakeOVVperExpDir Skipped, OVV EMPTY FOR: {exp_dir}")
            SkippedLst.append(exp_dir)

        for n, r in ovv.iterrows():
            if "ORR" in Path(r.PAR_file).stem.split("_") and r.PAR_exp != "ORR":
                #                    print(n,r)
                print(n, Path(r.PAR_file).name, r.PAR_exp, r.Gas)
                _logger.info("{0}".format([n, Path(r.PAR_file).name, r.PAR_exp, r.Gas]))
                ovv.loc[n, ["PAR_exp"]] = "ORR"
                ovv.loc[n, ["Gas"]] = "O2"
                self.EC_index.loc[n, ["PAR_exp"]] = "ORR"
                self.EC_index.loc[n, ["Gas"]] = "O2"
        ###### === Matching ORR scan with N2_activations in OVV ==== #######
        ovv[
            "ORR_act_N2_bg"
        ] = "None"  # Make new empty columns for filename of ORR background scan
        ovv_N2_act_local = ovv.query('PAR_exp == "N2_act"')
        #        ovv_N2_act_local = OnlyRecentMissingOVV.loc[(OnlyRecentMissingOVV.EXP_date == ovv.EXP_date.unique()[0])].query('PAR_exp == "N2_act"')
        for n, r in ovv_N2_act_local.iterrows():
            N2_dir_name = "_".join(
                [r.Electrolyte, r.SampleID, r.Loading_name, r.postAST]
            )
            N2_dest_dir = Path(r.Dest_dir).joinpath(
                f"N2_scans_v{FileOperations.version}/{N2_dir_name}"
            )
            N2_dest_dir.mkdir(parents=True, exist_ok=True)
            N2_act_BG = Path(N2_dir_name).joinpath(r.basename + f"_BG{N2_BG_ext}")
            #            print(r.basename+'_BG.xlsx')
            ovv.loc[n, "ORR_act_N2_bg"] = N2_act_BG
            #              ovv.loc[ovv['PAR_exp'] == "ORR" & ovv.SampleID == r.SampleID & ovv.postAST == r.postAST,'ORR-act_N2-bg'] = ('ORR;'+r.PAR_file)
            if not "ring" in N2_act_BG.stem and "Pt" not in str(r._act0_Comment):
                if not ovv.loc[
                    (
                        (ovv["PAR_exp"] != "N2_act")
                        & (ovv.SampleID == r.SampleID)
                        & (ovv.postAST == r.postAST)
                        & (ovv.Electrode == r.Electrode)
                        & (ovv.Loading_cm2 == r.Loading_cm2)
                        & (ovv.Electrolyte == r.Electrolyte)
                    )
                ].empty:
                    ovv.loc[
                        (
                            (ovv["PAR_exp"] != "N2_act")
                            & (ovv.SampleID == r.SampleID)
                            & (ovv.postAST == r.postAST)
                            & (ovv.Electrode == r.Electrode)
                            & (ovv.Loading_cm2 == r.Loading_cm2)
                            & (ovv.Electrolyte == r.Electrolyte)
                        ),
                        "ORR_act_N2_bg",
                    ] = str(N2_act_BG)
                else:
                    pass
            else:
                pass
        faillst = []
        #%% MakeOVVperExpDir 2
        for n, r in ovv.query(
            '(ORR_act_N2_bg == "None") & (PAR_exp == "ORR")'
        ).iterrows():
            N2_scans_folder_version = f"N2_scans_v{FileOperations.version}"
            N2_scan_dd = Path(r.Dest_dir).joinpath(N2_scans_folder_version)
            #            n,r = list(ovv.query('(ORR_act_N2_bg == "None") & (PAR_exp == "ORR") & (Electrode != "Pt_ring")').iterrows())[-1]
            #            n,r = faillst[0][0],faillst[0][-1]
            #            print(Path(r.PAR_file))
            #            Path(r.PAR_file).parent.name.split('_')[0]
            #            "N2_act",r.SampleID,r.postAST,r.Electrolyte,r.EXP_date, r.Loading_name
            #            N2_dest_dir = Path(gr_N2_ovv.Dest_dir.iloc[0]).joinpath('N2_scans/{0}'.format('_'.join([Electrolyte,SampleID,Loading_name])))
            #            N2_dest_dir.mkdir(parents=True,exist_ok=True)
            #            N2_act_BG = Path(N2_dest_dir.parent).joinpath(N2_fn+'_BG.xlsx')
            N2exMatch_ovv_local = ovv.loc[
                (
                    (ovv["PAR_exp"] == "N2_act")
                    & (ovv.SampleID == r.SampleID)
                    & (ovv.postAST == r.postAST)
                    & (ovv.Electrolyte == r.Electrolyte)
                    & (ovv.PAR_date_day == r.PAR_date_day)
                    & (ovv.Loading_name == r.Loading_name)
                )
            ]
            if (
                len(N2exMatch_ovv_local) == 1
            ):  # look for local N2 background inside of the experimental folder
                N2_bg_from_locl_N2match = N2exMatch_ovv_local.ORR_act_N2_bg.values[0]
                ovv.loc[n, "ORR_act_N2_bg"] = N2_bg_from_locl_N2match
            else:  # look for N2 background outside of the experimental folder
                N2exMatch = self.EC_index.loc[
                    (
                        (self.EC_index["PAR_exp"] == "N2_act")
                        & (self.EC_index.SampleID == r.SampleID)
                        & (self.EC_index.postAST == r.postAST)
                        & (self.EC_index.Electrolyte == r.Electrolyte)
                        & (self.EC_index.PAR_date_day == r.PAR_date_day)
                        & (self.EC_index.Loading_name == r.Loading_name)
                    )
                ]

                if N2exMatch.empty:
                    N2exMatch = self.EC_index.loc[
                        (
                            (self.EC_index["PAR_exp"] == "N2_act")
                            & (self.EC_index.SampleID == r.SampleID)
                            & (self.EC_index.postAST == r.postAST)
                            & (self.EC_index.Electrolyte == r.Electrolyte.split("+")[0])
                            & (self.EC_index.PAR_date_day == r.PAR_date_day)
                            & (self.EC_index.Loading_name == r.Loading_name)
                        )
                    ]
                #            if r.Electrode == 'Pt_ring' and len(N2exMatch) > 1:
                #                N2exMatch = N2exMatch.iloc[0]
                if len(N2exMatch) == 1:
                    N2matchPath = Path(*Path(N2exMatch.PAR_file.values[0]).parts[-2::])
                    #                print('N2_bg succes! (%s) added in folder for: %s' %(N2matchPath,r.basename))
                    _logger.info(
                        "N2_bg succes! ({0}) added in folder {2} for: {1}".format(
                            N2matchPath, r.basename, Path(r.Dest_dir).name
                        )
                    )
                    N2exMatch.Dest_dir = N2_scan_dd
                    N2_ext_BG = N2_scan_dd.joinpath(
                        N2exMatch.basename.values[0] + f"_BG{N2_BG_ext}"
                    )
                    N2_ext_BG = Path(run_ORR.N2_BG_exclusion_list(N2_ext_BG))
                    N2exMatch = N2exMatch.assign(**{"ORR_act_N2_bg": N2_ext_BG})
                    ovv.loc[n, "ORR_act_N2_bg"] = N2_ext_BG
                    ovv = ovv.append(N2exMatch)
                elif len(N2exMatch) > 1:
                    if not len(
                        self.EC_index.loc[
                            self.EC_index.PAR_hash.isin(N2exMatch.PAR_hash.unique())
                        ]
                    ) == len(N2exMatch):
                        _logger.info(
                            "N2_bg error duplicate files in EC_index! ({0}) added in folder {2} for: {1}".format(
                                r.PAR_file, r.basename, Path(r.Dest_dir).name
                            )
                        )
                    else:
                        _recent_N2match = N2exMatch.sort_values(
                            "PAR_date", ascending=False
                        ).head(1)
                        _logger.info(
                            "N2_bg multiple options take most recent! ({0}) added in folder {2} for: {1}".format(
                                _recent_N2match.PAR_file.values[0],
                                r.basename,
                                Path(r.Dest_dir).name,
                            )
                        )
                        _recent_N2match.Dest_dir = N2_scan_dd
                        N2_ext_BG = N2_scan_dd.joinpath(
                            _recent_N2match.basename.values[0] + f"_BG{N2_BG_ext}"
                        )
                        N2_ext_BG = Path(run_ORR.N2_BG_exclusion_list(N2_ext_BG))
                        _recent_N2match = _recent_N2match.assign(
                            **{"ORR_act_N2_bg": N2_ext_BG}
                        )
                        ovv.loc[n, "ORR_act_N2_bg"] = N2_ext_BG
                        ovv = ovv.append(_recent_N2match)
                else:
                    N2_ext_BG = "fail"
                    faillst.append([n, r.basename, N2_ext_BG, r])
                    if r.Electrode != "Pt_ring":
                        if r.PAR_exp == "ORR":
                            _logger.warning(
                                "MakeOVVperExpDir N2 BG match EMPTY fail, N2_act_BG ORR, {0}\n for folder {1}".format(
                                    r.basename, exp_dir
                                )
                            )
        #            print(r.basename,'\n',Path(N2_ext_BG).name)
        ovv = ovv.drop_duplicates(subset=["PAR_hash", "PAR_file"])
        if ovv.empty:
            #            print('OVV empty skipped')
            _logger.warning("OVV empty skipped: {0}".format(exp_dir))
        RHE_index = get_RHE_OCP(ovv)
        ovv = RHE_index[0]
        #%% MakeOVVperExpDir end
        return ovv

    def update_run_kwargs(self, ttpfs):

        self.run_kwargs.update(
            {"N2_BG_ext": ".pkl", "multi_par_fit": check_interactive_for_multi()}
        )
        if "testing" in self.run_kwargs["input_run"]:
            tt4 = "N2_20cls_300_100_10_JOS3_high-load_269"

            #            tt4 = 'O2_ORR_DW16_data_PDX_Ch1'
            #            'N2_20cls_300_100_10_DW28_898'
            self.EC_run_slice = self.EC_run_slice.loc[
                self.EC_run_slice.basename.str.contains(tt4)
            ]
            self.EC_run_slice = self.EC_run_slice.loc[
                self.EC_run_slice.PAR_file.isin(ttpfs)
            ]
        #            EC_run_slice = EC_run_slice.iloc[-10:-7]
        if "ORR" in self.EC_run_slice.PAR_exp.unique():
            self.EC_run_slice = self.EC_index.loc[
                self.EC_index.PAR_date_min.isin(self.EC_run_slice.PAR_date_min.unique())
            ]

    def prepare_ovv_EC_run(self):

        results, out = [], []
        EC_run_groupby = ["PAR_exp", "EXP_date", "EXP_dir"]
        run_group_template = namedtuple("Run_groups", "PAR_exp EXP_date EXP_dir grp")
        EC_run_slice_grp_date_dir = self.EC_run_slice.groupby(by=EC_run_groupby)
        #        EC_run_slice_grp_date_dir = EC_run_slice.loc[EC_run_slice.PAR_file.isin(_n2faillst)].groupby(by=EC_run_groupby)
        run_group_lst = [
            run_group_template(n[0], n[1], n[2], gr)
            for n, gr in EC_run_slice_grp_date_dir
        ]
        ovv_all_date_dirs = pd.concat(
            [ECRunOVV.MakeOVVperExpDir(self, run_group) for run_group in run_group_lst]
        )

        ovv = ovv_all_date_dirs.loc[
            ovv_all_date_dirs.SampleID.isin(self.EC_run_slice.SampleID.unique())
        ]
        self.run_ovv = ovv
        ovv_exp_grp = ovv.groupby("PAR_exp")
        ovv_lst_exps = ovv.PAR_exp.unique()
        #        PAR_exp_grp_lst = [ovv.groupby(by='PAR_exp' ]
        _logger.info(
            f"=== PAR DW starting: groups {len(EC_run_slice_grp_date_dir)},\nLength slice {len(self.EC_run_slice)}, index {len(self.EC_index)}  ===="
        )
        self.ovv_exp_grp = ovv_exp_grp

    ###### === EC Analysis Run over test1 in multiprocess ==== #######
    def update_type_run_kwargs(self):
        #        if 'N2_act' in _exp_type:
        #        elif 'ORR' in _exp_type or 'O2_nan' in _exp_type:
        if "EIS" in self.ovv_exp_grp.groups.keys():
            # ==== UPDATING EIS KWARGS =====
            #                    EIS_kwargs = {'EIS_skip_set' : False, 'EIS_use_prepars_set' : True, 'FreqLim_set' : 30E3}
            EIS_kwargs = dict(
                EIS_skip_set=False,
                EIS_use_prepars_set=True,
                FreqLim=15e3,
                EIS_plot_combined=True,
                EIS_single_output="Text,Plot",
                perform_prefit=True,
                TrimData=False,
                FitOnlyTrimmedData=False,
                linKK_trimming_factor=3,
                export_raw_data=True,
                DP_DRT_fit=False,
                GP_DRT_fit=False,
            )
            EIS_kwargs.update(run_kwargs)
            print(
                f'Skip EIS: {EIS_kwargs["EIS_skip_set"]}. Use prePars: {EIS_kwargs["EIS_use_prepars_set"]}, Frequency limit: {EIS_kwargs["FreqLim"]:.0f}'
            )
            EIS_kwargs.update(run_EIS.get_eis_suggestions(), compare=self.EC_run_slice)
            self.run_kwargs.update(**EIS_kwargs)

    #        self.run_kwargs.update(**dict(ovv_exp_grp = self.ovv_exp_grp))
    #            self.run_kwargs.update({'EIS' : {**EIS_kwargs}})
    #        elif 'HPRR' in _exp_type:
    #        elif 'OER' in _exp_type:
    #        elif 'HER' in _exp_type:
    #        elif 'RHE' in _exp_type:
    def start_N2(self, run_kwargs, _exp_type="N2_act"):

        try:
            index_fls = {
                "N2_CVs": {"path": self.POST_DIR.joinpath("N2_CVs_all.pkl")},
                "N2_actions": {
                    "path": self.POST_DIR("VERSASTAT").PostDir.joinpath(
                        "N2_CVs_actions.pkl"
                    )
                },
            }
            run_N2.N2_act(self.ovv_exp_grp, index_fls=index_fls, **self.run_kwargs)
            pass
        except Exception as e:
            _logger.error(
                f'Run N2_act {run_kwargs["input_run"]} failed for {_exp_type}, because {e}'
            )

    def EC_Analysis_run(self):
        #%% EC_Analysis_run
        for _exp_type, _exp_grp in self.ovv_exp_grp:
            _logger.warning(
                f"Starting {_exp_type} len({len(_exp_grp)}) {self.__class__.__name__} from {self.input_run}. =========="
            )
            if "N2_act" in _exp_type:
                self.start_N2()

            elif "ORR" in _exp_type or "O2_nan" in _exp_type:
                try:
                    # run_ORR.ORR(self.ovv_exp_grp, **self.run_kwargs)
                    run_ORR.ORR_run_loop(self)
                except Exception as e:
                    _logger.error(
                        f'Run ORR {run_kwargs["input_run"]} failed for {_exp_type}, because {e}'
                    )
            elif "EIS" in _exp_type and not "HER" in _exp_type:
                #                    EIS_kwargs = {'EIS_skip_set' : False, 'EIS_use_prepars_set' : True, 'FreqLim_set' : 30E3}
                #                EIS_kwargs = dict(EIS_skip_set = False,EIS_use_prepars_set = True, FreqLim = 25E3,
                #                                  EIS_plot_combined=True, EIS_single_output = 'Text,Plot',
                #                                  perform_prefit = True, TrimData = False, FitOnlyTrimmedData = False, linKK_trimming_factor = 3,
                #                                  export_raw_data = True, DP_DRT_fit = False, GP_DRT_fit = False)
                #                EIS_kwargs.update(run_kwargs)
                #                print(f'Skip EIS: {EIS_kwargs["EIS_skip_set"]}. Use prePars: {EIS_kwargs["EIS_use_prepars_set"]}, Frequency limit: {EIS_kwargs["FreqLim"]:.0f}')
                #    #                    exp, gr = 'EIS', ExpTypes_gr.get_group('EIS')
                if self.run_kwargs.get("EIS_skip_set", False) == False:
                    try:
                        _logger.warning(
                            f"=====Starting EIS {self.__class__.__name__} from {self.input_run}. =========="
                        )
                        #                        eis_run_ovv.eis_run_group_ovv(self.ovv_exp_grp,self.run_kwargs)
                        #                        eis_run_ovv.EIS_Preparator(self.ovv_exp_grp, **self.run_kwargs)
                        self.run_kwargs.update({})
                        run.EIS_run_loop(self)
                        # eis_run_ovv.EIS_Preparator(self)
                    except Exception as e:
                        _logger.error(
                            f'EIS run failed command "{run_kwargs["input_run"]}" failed for {_exp_type}, because {e}'
                        )
                    # TODO Fix run EIS_HER files, only 2 with 0 rpm...
            #                        ovv_Dest_dir.joinpath('{})
            #                        index_info.append(index)
            elif "HPRR" in _exp_type:
                pass
            elif "OER" in _exp_type:
                try:
                    run_OER.OER(self.ovv_exp_grp, **self.run_kwargs)
                except Exception as e:
                    _logger.info(
                        f'Run {run_kwargs["input_run"]} failed for {_exp_type}, because {e}'
                    )

            elif "HER" in _exp_type and not "EIS" in _exp_type:
                try:
                    run_HER.HER(self.ovv_exp_grp, **self.run_kwargs)
                except Exception as e:
                    _logger.info(
                        f'Run {run_kwargs["input_run"]} failed for {_exp_type}, because {e}'
                    )
            elif "RHE" in _exp_type:
                pass
            else:
                #                    print('No run, unknown experiment type:', exp)
                _logger.info("No run, unknown experiment type: {0}".format(_exp_type))


###### === EC Analysis Run over test1 grouped by Date ==== #######
#%% Testing functions


def gas_filter(test1, run):
    gas_search = re.search("(O2|N2)", run)
    if gas_search:
        gas_select = gas_search.group()
        test1 = test1.loc[test1.Gas == gas_select]
    return test1


def PF_select(key):
    selector = {
        "eispASTsHA": [
            "O2_EIS-range_1500rpm_JOS1_pAST-sHA_285",
            "N2_EIS-range_1500rpm_pAST-sHA_JOS1_288",
            "N2_EIS-range_1500rpm_JOS3_288",
        ]
    }

    return selector.get(key, [])


def check_idx(EC_index, test_obj):
    # EIS_EC_index = EC_run_slice.reset_index()
    #    loc[EC_run_slice.PAR_exp.str.contains('EIS')].reset_index()
    test_obj = "N2_EIS-range_0rpm_JOS2_272"
    test_obj = "O2_EIS-range_1500rpm_JOS9_r1_283"
    test_obj = "O2_EIS-range_1500rpm_DW51_280"
    test_obj = "O2_EIS-range_1500rpm_JOS9_899"
    test_obj = "N2_EIS-range_0rpm_JOS7_264"
    tj = "N2_EIS-range_1500rpm_JOS4_268"
    tj = "O2_EIS-range_1500rpm_JOS3_270"
    tj = "O2_EIS-range_1500rpm_JOS3_285"
    tj = "O2_EIS-range_1500rpm_JOS3_pAST-sHA_285"
    tj = "O2_EIS-range_1500rpm_JOS1_pAST-sHA_285"
    tj = "O2_EIS-range_1500rpm_JOS1_pAST-sHA_285"
    tj = "N2_EIS-range_1500rpm_pAST-sHA_JOS1"
    tj = "O2_ORR_JOS4_257"

    test_idx = EC_index.loc[EC_index.basename.str.contains(tj)].index[0]
    #    test_idx = EIS_EC_index.loc[EIS_EC_index.basename.str.contains(test_obj)].index[0]
    print(
        f"index: {test_idx}, other {100*[n for n,i in enumerate(EC_index.index) if i == test_idx][0]/len(EC_index):.0f}"
    )
    return


#%% === Index Run Selection ===


def index_run_selection(EC_index, **run_kwargs):
    run = run_kwargs.get("input_run", "n empty")
    print(f"Run command: {run}")

    if run.startswith("y"):
        if "ya" in run or run == "yes all" or run == "y std":
            test1 = EC_index
        #        test3 = OnlyRecentMissingOVV.query('(EXP_date >= 20190301)  & ((SampleID >= "DW01") | (SampleID == "JOS12") | (SampleID == "JOS13") | (SampleID == "JOS14") | (SampleID == "JOS15")) & (PAR_exp == "EIS")')
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20190103)  & (SampleID == "JOS15") & (PAR_exp == "EIS")')
        elif run == "ys":
            test1 = EC_index[
                EC_index.EXP_date.isin(EC_index.query('PAR_exp == "HPRR" ').EXP_date)
            ]
        elif "yrecent" in run:

            test1 = EC_index.query("(PAR_date >= 20190506)").loc[
                EC_index.SampleID.str.contains("JOS4|DW")
            ]
            #             .head(20)
            if "slicedate" in run:
                test1 = test1.loc[
                    (test1.Date_PAR_EXP > pd.Timedelta("-2 days 0 hours")), :
                ]
        elif "testN2EIS" in run:
            #            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20191020)')
            #            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190325)')
            #            test1 = EC_index.loc[EC_index.PAR_file.str.contains('H2O2 Daten_BA')]
            _ll = [
                "/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT/2019-01-Jan/25.01.2019_0.1MH2SO4_cell2/N2_EIS-range_0rpm_HER_JOS4_257.par",
                "/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT/2019-01-Jan/25.01.2019_0.1MH2SO4_cell3/N2_EIS-range_0rpm_JOS4_postAST-LC_257.par",
            ]
            test1 = EC_index.loc[EC_index.PAR_file.isin(_ll)]
        #                                 str.contains('/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT/2019-03-Mrt/25.03.2019_0.1MH2SO4_cell2/O2_ORR_JOS3_270_2_#2_Disc_Parstat.par')]
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190912) & (pH < 7) & (SampleID == "JOS2")')
        elif "highload" in run:
            test1 = EC_index.loc[EC_index.PAR_file.str.contains("high-load_267")].query(
                'PAR_exp != "EIS" & EXP_date == 20190912'
            )
        #            'O2_ORR_JOS12_3rpm_258_#2_Disc_Parstat'
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20190101)')
        elif "serie" in run:
            #            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20191020)')
            #            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190325)')
            serm = [
                i
                for i in run.split()
                if i in PostChar.SampleSelection.Series.index.values
            ]
            test1 = EC_index.query(
                '(EXP_date >= 20190101) & (pH < 7) & (SampleID == "DW28") & (PAR_file.str.contains("3rpms"))'
            )
        #            EC_run_slice.loc[EC_run_slice.PAR_file.str.contains('3rpms')]
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date <= 20190228) & (EXP_date >= 20190101)')
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date >= 20180101) & (EXP_date <= 20180131)')
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20180123) & SampleID_folder == "DW28"')
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date == 20190125)')
        #            PAR_date < pd.to_datetime('20190901') and PAR_date > pd.to_datetime('20190827')
        elif "missing" in run:
            #            EIS_missing = FileHelper.FindExpFolder('VERSASTAT').EISmissing
            #            EIS_missing = pd.read_excel(FileHelper.FindExpFolder('VERSASTAT').PostDir.joinpath('OVV_EIS_missing.xlsx'))
            #            test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.EXP_date.isin(EIS_missing.EXP_date.unique())]
            #            refits = 'DW16','DW19','DW17','DW28'
            #            test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.SampleID.isin(refits) & OnlyRecentMissingOVV.pH == 1]
            test1 = EC_index.loc[
                EC_index.basename.str.contains("SDFe2AL", case=False)
            ].query("pH < 7")
        #            test1 = OnlyRecentMissingOVV.loc[OnlyRecentMissingOVV.SampleID.str.contains('SD')]

        elif "orrmiss" in run:
            PostDestDir = FindExpFolder("VERSASTAT").PostDir
            test1 = pd.read_pickle(PostDestDir.joinpath("ORR_missing.pkl.compress"))
        elif "eismiss" in run:
            PostDestDir = FindExpFolder("VERSASTAT").PostDir
            eis_refit = pd.read_pickle(
                PostDestDir.joinpath("EIS_ORR_refit_pars.pkl.compress")
            ).PAR_file.unique()
            test1 = EC_index.loc[EC_index.PAR_file.isin(eis_refit)].tail(1)
        elif "all eis" in run:
            EIS_EC_index = EC_index.loc[
                EC_index.PAR_exp.str.contains("EIS")
            ].reset_index()
            test1 = EIS_EC_index
            if "JOS" in run:
                test1 = EIS_EC_index.loc[
                    EIS_EC_index.SampleID.isin(["JOS1", "JOS2", "JOS3", "JOS4", "JOS5"])
                    & EIS_EC_index.basename.str.contains("1500")
                ]
                if "JOSn" in run:
                    test1 = EIS_EC_index.iloc[~EIS_EC_index.index.isin(test1.index)]
                if "JOS4" in run:
                    # pass
                    test1 = EIS_EC_index.loc[
                        (EIS_EC_index.SampleID.str.contains("JOS4"))
                        & EIS_EC_index.basename.str.contains("1500")
                        & (EIS_EC_index.Gas == "O2")
                    ].tail(1)

            if "continue" in run:
                test1 = test1.loc[
                    test1.PAR_date
                    >= test1.loc[
                        test1.PAR_file.str.contains("O2_EIS-range_0rpm_HER_JOS5_257"),
                        "PAR_date",
                    ].iloc[0]
                ]
                test1 = test1.loc[test1.PAR_date <= test1.loc[1734, "PAR_date"]]
            elif "rpm" in run:
                # print("RPM")
                eisrpm_miss = pd.read_pickle(
                    FindExpFolder("VERSASTAT").PostDir.joinpath(
                        "EIS_RPM_series.pkl.compress"
                    )
                )
                # print('eisrpm',len(eisrpm_miss))
                test1 = EC_index.loc[
                    EC_index.PAR_file.isin(
                        [Path(i) for i in eisrpm_miss.PAR_file.unique()]
                    )
                ]
                # print('test',len(test1))
            elif "recent" in run:
                eis_miss_recent = pd.read_pickle(
                    FindExpFolder("VERSASTAT").PostDir.joinpath("EIS_pars_nonrecent")
                )
                test1 = EC_index.loc[
                    EC_index.PAR_file.isin(
                        [Path(i) for i in eis_miss_recent.PAR_file.unique()]
                    )
                ]

            elif "rest" in run:
                eis_metaf = pd.read_excel(
                    list(
                        FindExpFolder("PorphSiO2").compare.parent.rglob(
                            "2020-24-03_EIS_Porph_SiO2/meta_data*EIS*origin.xlsx"
                        )
                    )[0],
                    index_col=[0],
                )
                porhp_refit = EIS_EC_index.loc[
                    EIS_EC_index.basename.isin(eis_metaf.basename.unique())
                ].index
                tt_inf = "2019-02-Feb/12.02.2019_0.1MKOH_cell2/O2_EIS-range_1500rpm_JOS6_postAST_899"
                idx_dbg = EIS_EC_index.loc[
                    EIS_EC_index.PAR_file.str.contains(tt_inf)
                ].index[0]
                not_idx = list(porhp_refit) + [idx_dbg + i for i in range(-7, 7)]
                test1 = EIS_EC_index.loc[
                    (~EIS_EC_index.index.isin(not_idx)) & (EIS_EC_index.index > 215)
                ]

        elif "porph" and "refit" in run:
            eis_metaf = pd.read_excel(
                list(
                    FindExpFolder("PorphSiO2").compare.parent.rglob(
                        "2020-24-03_EIS_Porph_SiO2/meta_data*EIS*origin.xlsx"
                    )
                )[0],
                index_col=[0],
            )
            test1 = EC_index.loc[EC_index.basename.isin(eis_metaf.basename.unique())]
            if "not" in run:
                _not_porph = EC_index.loc[
                    ~EC_index.SampleID.isin(eis_metaf.SampleID.unique())
                    & EC_index.PAR_exp.str.contains("EIS")
                ]
                test1 = _not_porph.loc[
                    _not_porph.SampleID.str.contains("DW|JOS12|JOS14|JOS15")
                ]

            elif "samples" in run:
                test1 = EC_index.loc[
                    EC_index.SampleID.isin(eis_metaf.SampleID.unique())
                    & ~EC_index.PAR_exp.str.contains("EIS")
                ]
                if "postAST" in run:
                    test1 = EC_index.loc[
                        EC_index.SampleID.isin(eis_metaf.SampleID.unique())
                        & EC_index.PAR_exp.str.contains("EIS")
                        & EC_index.basename.str.contains("AST")
                    ]
        #                & (~EC_index.PAR_file.isin(eis_metaf.PAR_file.unique()))]
        elif "JOS eismetaf" in run:
            eis_metaf = pd.read_excel(
                list(
                    FindExpFolder("PorphSiO2").compare.parent.rglob(
                        "2020-24-03_EIS_Porph_SiO2/meta_data*EIS*origin.xlsx"
                    )
                )[0],
                index_col=[0],
            )
            JOS2 = pd.DataFrame()
            if "ttN2" in run:

                JOS2 = eis_metaf.loc[
                    eis_metaf.basename.str.contains("N2_EIS-range_1500rpm_JOS2_288")
                ]
                JOS2 = eis_metaf.loc[
                    eis_metaf.basename.str.contains("N2_EIS-range_1500rpm_JOS4_288")
                ]
                JOS2 = eis_metaf.loc[
                    eis_metaf.basename.str.contains("N2_EIS-range_1500rpm_JOS1_288")
                ]

            elif "ttO2" in run:
                JOS4 = eis_metaf.loc[
                    eis_metaf.basename.str.contains("O2_EIS-range_1500rpm_JOS5_285")
                ]
                JOS2 = eis_metaf.loc[
                    eis_metaf.basename.str.contains("O2_EIS-range_1500rpm_JOS2_288")
                ]
                JOS2 = eis_metaf.loc[
                    eis_metaf.basename.str.contains("O2_EIS-range_1500rpm_JOS1_285")
                ]

            else:
                JOS4 = eis_metaf.loc[
                    eis_metaf.basename.str.contains("EIS-range_1500rpm_JOS2")
                ]
            #                JOS4 = eis_metaf.loc[eis_metaf.basename.str.contains('JOS5')]
            #                eis_metaf.loc[eis_metaf.basename.str.contains('O2_EIS-range_1500rpm_JOS5_285|N2_EIS-range_1500rpm_JOS5_285')]
            if not JOS2.empty:
                eis_metaf = JOS2
            test1 = EC_index.loc[EC_index.basename.isin(eis_metaf.basename.unique())]
        elif "pAST-sHA_JOS7" in run:
            tj = "EIS-range_1500rpm_pAST-sHA_JOS7_288"
            test1 = EC_index.loc[EC_index.basename.str.contains(tj)]

        elif "DW17DW29" in run:
            test1 = EC_index.loc[
                (EC_index.SampleID.str.contains("DW17|DW29|DW28"))
                & (EC_index.PAR_exp.str.contains("EIS"))
            ]

        elif "ORRselect JOS4" in run:
            tt4 = "O2_ORR_JOS4_257_#2_Disc_Parstat"
            tt4 = "O2_ORR_JOS4_postAST-LC_257_#3_Disc_Parstat"
            JOS4 = EC_index.loc[EC_index.basename.str.contains(tt4)]
            test1 = EC_index.loc[EC_index.basename.isin(JOS4.basename.unique())]
        elif "ORRporph":
            _smpls = ["JOS1", "JOS2", "JOS3", "JOS4", "JOS5"]
            test1 = EC_index.loc[
                (EC_index.SampleID.isin(_smpls))
                & (EC_index.PAR_exp.str.contains("ORR"))
            ]

        elif "eispASTsHA" in run:
            test1 = EC_index.loc[EC_index.basename.isin(PF_select("eispASTsHA"))]
        #            .query('SampleID == "JOS1"') # FIXME
        elif "sr" in run:
            test1, EC_index = prepare_input.MainEC().EC_Analysis_Input()
            serie = PostChar.SampleSelection("CB4", "SeriesID", [])
            serie_samples_ovv = pd.merge(
                serie.Prep_EA_BET.SampleID, EC_index, how="left", on="SampleID"
            )
            EC_run_Serie_OVV = pd.merge(
                serie_samples_ovv.EXP_dir.drop_duplicates(),
                EC_index,
                how="left",
                on="EXP_dir",
            )
        elif "debug" in run:
            tt_inf = "2019-02-Feb/12.02.2019_0.1MKOH_cell2/O2_EIS-range_1500rpm_JOS6_postAST_899"
            EIS_EC_index = EC_index.loc[
                EC_index.PAR_exp.str.contains("EIS")
            ].reset_index()
            idx_dbg = EIS_EC_index.loc[
                EIS_EC_index.PAR_file.str.contains(tt_inf)
            ].index[0]
            test1 = EIS_EC_index.iloc[idx_dbg - 7 : idx_dbg + 7]

        elif "Nicole" in run:
            test1 = EC_index.loc[EC_index.PAR_file.str.contains("Nicole")]

        elif "now" in run:
            today = pd.datetime.now().strftime("%Y%m%d")
            print("run Today only")
            test1 = EC_index.query("(EXP_date == @today)")
            if test1.empty:
                last_day = EC_index.sort_values(by="EXP_date").EXP_date.unique()[-1]
                test1 = EC_index.query("(EXP_date == @last_day)")
        else:
            test1 = EC_index.head(1)
        #            & (Loading_name == "high") & (SampleID == "DW21")' )
        #            test1 = OnlyRecentMissingOVV[OnlyRecentMissingOVV.EXP_date.isin(OnlyRecentMissingOVV.EXP_date.unique()[-3::])]
        #            OnlyRecentMissingOVV.query('(EXP_date > 20190715)').tail(30)
        #            test1 = OnlyRecentMissingOVV.query('(EXP_date > 20180101) & (EXP_date < 20181212)')
        #        & (Gas == "O2")')
        #        test1,OnlyRecentMissingOVV = MainEC().EC_Analysis_Input()
        test1 = gas_filter(test1, run)
        ExpTypes = test1.PAR_exp.unique()

        RunTypes = [i for i in run.split() if i in ExpTypes]
        if RunTypes:
            test1 = test1.loc[test1.PAR_exp.str.contains("|".join(RunTypes))]

        if "ORR" in test1.PAR_exp.unique():
            #            test1 = EC_index.loc[EC_index.PAR_date.isin(test1.PAR_date.values)] # adds the Ring ORR PAR files
            test1 = EC_index.loc[
                (EC_index.PAR_date.dt.date.isin(test1.PAR_date.dt.date.values))
                & (EC_index.PAR_exp.isin(["N2_act", "ORR"]))
                & (EC_index.SampleID.isin(test1.SampleID.unique()))
            ]  # adds the Ring ORR PAR files

        Samples = [i for i in run.split() if i in test1.SampleID.unique()]
        RunSeries = [
            i for i in run.split() if i in PostChar.SampleSelection.Series.index.values
        ]
        if Samples:
            test1 = test1.loc[test1.SampleID.str.contains("|".join(Samples))]
        if RunSeries:
            test1 = test1.loc[
                test1.SampleID.isin(
                    PostChar.SampleSelection.Series.loc[RunSeries, "sIDs"].values[0]
                )
            ]
        if "short" in run:
            test1 = test1.iloc[0:2]
        MultiRun = [i for i in run.split() if "multi" in i]
        multi_set = bool(MultiRun)

        if "OER testing" in run:
            OER_testing = pd.concat(
                [(gr) for n, gr in test1.groupby("PAR_date") if len(gr) > 2]
            ).query("(EXP_date == 20190308)")
            test1 = OER_testing

        if "only" and "bad" in run:
            bad_PFs = set(
                [
                    n
                    for n, gr in run_kwargs.get("EIS_recheck_bad_fits")
                    if n[0] in test1.PAR_file.unique()
                ]
            )
            run_kwargs.update({"BadOnly": bad_PFs})
            test1 = test1.loc[test1.PAR_file.isin([i[0] for i in bad_PFs])]
            test1.loc[~test1.PAR_file.isin(bad_PFs)]

    elif run.startswith("n"):
        #        ECRunOVV.EC_Analysis_Run(EC_run_Serie_OVV,EC_index)
        if re.match("index(?!(force|.force))", run.split("n ")[-1]):
            EC_index = prepare_input.MainEC.EC_Analysis_Input(
                TakeRecentList=False, exp_folder=None, force_recalc="force"
            )
            test1 = pd.DataFrame()
        elif re.match("index(?=(force|.force))", run.split("n ")[-1]):
            test1 = pd.DataFrame()
            if "folder" in run:
                exp_folder = "25.07.2019_0.5MH2SO4_LM"
                exp_folder = "2020-05-01_Nicole-EIS"
                EC_index = prepare_input.MainEC.EC_Analysis_Input(False, exp_folder)
            else:
                EC_index = prepare_input.MainEC.EC_Analysis_Input(
                    TakeRecentList=False, exp_folder=None, force_recalc="force"
                )
        elif "only post" in run:
            # Post_EIS_pars(reload=True)
            test1 = pd.DataFrame()
        else:
            test1 = pd.DataFrame()
    else:
        test1 = pd.DataFrame()
    return test1, run_kwargs


#%% == MAIN ==
if __name__ == "__main__":
    # run = input('Want to start the fitting run?')
    #    _logger = start_logging(__name__)
    #    run = 'y porph_refit'
    #    run = 'y all eis continue'
    run = [
        "y all eis recent y",
        "y JOS eismetaf ttO2 y",
        "y eispASTsHA y",
        "y JOS eismetaf y",
        "y pAST-sHA_JOS7 y",
        "y all eis JOS y",
        "y DW17DW29 y",
    ][3]
    run = ["y ORRselect JOS4 ", "y ORRporph"][0]
    run = "n "

    # run ='n '
    run = "n index force "
    run = "y N2"

    # run = []
    _start = "n"
    #    run = input('What to run for ECpy?')

    run_kwargs = dict(
        skip=False,
        input_run="y all eis JOS n",
        load_index=1,
        multi_par_fit=True,
        run_fit_mean="lmfit",
    )
    if run:
        run = run[:-1] + " " + _start
        run_kwargs.update({"input_run": run})

    #    try:
    #        if not EC_index.empty:
    #            pass
    #        else:
    #            EC_index = ECRunOVV(load=1).index.sort_values('PAR_date',ascending=False).reset_index()
    #    except Exception as e:
    #        print('reload EC index')
    ##        EC_index = prepare_input.MainEC.EC_Analysis_Input(False,force_recalc='force')
    #        EC_index = ECRunOVV(load=1).index.sort_values('PAR_date',ascending=False).reset_index()

    #    EC_run_slice,run_kwargs = index_run_selection(EC_index,**run_kwargs)
    test = ECRunOVV(**run_kwargs)

    # if not '000' in run_kwargs['input_run'] and run_kwargs['input_run'].startswith('y') and not run_kwargs['input_run'].rstrip().endswith(' n'):
    #     ECRunOVV.EC_Analysis_Run(EC_run_slice,EC_index,**run_kwargs)
    #     if not EC_run_slice.empty:
    #         if 'EIS' in EC_run_slice.PAR_exp.unique() and 'post' in run:
    #             Post_EIS_pars()
