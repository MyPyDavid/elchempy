# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:30:31 2020

@author: User
"""

import sys
import datetime as dt
from collections import Counter
import pprint

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

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper.PostChar import (
    SampleSelection,
    Characterization_TypeSetting,
    SampleCodesChar,
)


if __name__ == "__main__":
    print(f"Package: {__package__}, File: {__file__}")
    from ECpy.main_run_PAR_DW import ECRunOVV
    from ECpy.indexer.prepare_input import CleanUpCrew
    from ECpy.experiments.EIS.models import Model_Collection
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

    # logger = start_logging(__name__)

else:
    #    print('\n\n***** run_PAR_DW *****')
    print(f"File: {__file__}, Name:{__name__}, Package:{__package__}")
    #    FH_path = Path(__file__).parent.parent.parent
    #    sys.path.append(str(FH_path))
    #    import FileHelper
    from ECpy.main_run_PAR_DW import ECRunOVV
    from ECpy.indexer.prepare_input import CleanUpCrew
    from ECpy.runEC.EC_logging_config import start_logging
    from ECpy.PostEC import post_helper, merger
    from ECpy.experiments.EIS.models import Model_Collection

    # logger = start_logging(__name__)
_logger = logging.getLogger(__name__)
_logger.setLevel(20)

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


class EnterExitLog:
    def __init__(self, funcName):
        self.funcName = funcName

    def __enter__(self):
        _logger.info(f"Started: {self.funcName}")
        self.init_time = dt.datetime.now()
        return self

    def __exit__(self, type, value, tb):
        self.end_time = dt.datetime.now()
        self.duration = self.end_time - self.init_time
        _logger.info(f"Finished: {self.funcName} in {self.duration} seconds")


def func_timer_decorator(func):
    def func_wrapper(*args, **kwargs):
        with EnterExitLog(func.__name__):
            return func(*args, **kwargs)

    return func_wrapper


def get_daily_pickle(exp_type=""):
    today = dt.datetime.now().date()
    _result = {"today": today}
    if exp_type:
        daily_pickle_path = FindExpFolder("VERSASTAT").PostDir.joinpath(
            f"{today:%Y-%m-%d}_{exp_type}_{system()}.pkl.compress"
        )
        daily_pkl_options = list(
            FindExpFolder("VERSASTAT").PostDir.rglob(
                f"*_{exp_type}_{system()}.pkl.compress"
            )
        )
        daily_pkl_options = sorted(daily_pkl_options, key=lambda x: x.stat().st_ctime)

        _result.update(
            {
                "daily_path": daily_pickle_path,
                "_exists": daily_pickle_path.exists(),
                "daily_options": daily_pkl_options,
            }
        )

        daily_pickle_path_RAW = FindExpFolder("VERSASTAT").PostDir.joinpath(
            f"{today:%Y-%m-%d}_{exp_type}_{system()}_RAW.pkl.compress"
        )
        daily_pkl_options_RAW = list(
            FindExpFolder("VERSASTAT").PostDir.rglob(
                f"*_{exp_type}_{system()}_RAW.pkl.compress"
            )
        )
        daily_pkl_options_RAW = sorted(
            daily_pkl_options_RAW, key=lambda x: x.stat().st_ctime
        )
        _result.update(
            {
                "daily_path_RAW": daily_pickle_path_RAW,
                "_raw_exists": daily_pickle_path_RAW.exists(),
                "daily_options_RAW": daily_pkl_options_RAW,
            }
        )
        if "EIS" in exp_type:
            _result.update(
                {
                    "daily_path_BRUTE": FindExpFolder("VERSASTAT").PostDir.joinpath(
                        f"{today:%Y-%m-%d}_{exp_type}_BRUTE_{system()}.pkl.compress"
                    ),
                    "daily_path_RAW_WB": FindExpFolder("VERSASTAT").PostDir.joinpath(
                        f"{today:%Y-%m-%d}_{exp_type}_RAW_WB_{system()}.pkl.compress"
                    ),
                }
            )

    return _result


def _collect_test():
    tt = CollectLoadPars(load_type="fast")


class CollectLoadPars:
    def __init__(self, load_type="fast"):
        self.load_type = load_type
        self.load_pars()
        self.collect_dict()

    def load_pars(self):
        _BaseLoad = BaseLoadPars()
        _kws = {"EC_index": _BaseLoad.EC_index, "SampleCodes": _BaseLoad.SampleCodes}
        if "fast" in self.load_type:
            _kws.update(**{"reload": False, "reload_raw": False})
            self.EIS_load = EIS_LoadPars(**_kws)
            self.ORR_load = ORR_LoadPars(**_kws)
            self.N2_load = N2_LoadPars(**_kws)

    def collect_dict(self):
        _load_attrs = [i for i in self.__dict__.keys() if i.endswith("_load")]
        _collect = {}
        for _load_pars in _load_attrs:
            _pars_name = f'{_load_pars.split("_")[0]}_pars'
            if hasattr(getattr(self, _load_pars), _pars_name):
                _pars = getattr(getattr(self, _load_pars), _pars_name)
                _collect.update({_pars_name: _pars})
        self.pars_collection = _collect


class BaseLoadPars:

    _required_funcs = [
        "make_raw_pars_from_scratch",
        "edit_raw_columns",
        "search_pars_files",
        "read_in_pars_files",
        "extra_stuff_delegator",
    ]

    def __init__(
        self,
        EC_index=pd.DataFrame(),
        SampleCodes=pd.DataFrame(),
        exp_type="",
        reload=False,
        reload_raw=False,
    ):
        self.exp_type = exp_type
        self._auto_set_exp_type()
        self.EC_index = EC_index
        self.SampleCodes = SampleCodes

        self._check_class_req_functions()
        self.check_EC_index()
        self.set_OVV_exp_type()

        self._reload = reload
        self._reload_raw = reload_raw

        self.get_daily_pickle()

        if self.exp_type:
            self.load_delegator()

    def _auto_set_exp_type(self):
        _cls_name = self.__class__.__name__
        if "_" in _cls_name:
            _cls_exp_type = _cls_name.split("_")[0]
            _exp_type = f"{_cls_exp_type}_pars"
            self.exp_type = _exp_type

    def check_EC_index(self):
        if self.EC_index.empty:
            EC_index = ECRunOVV(load=1).EC_index
            EC_index = FileOperations.ChangeRoot_DF(EC_index, [])
            EC_index.PAR_file = EC_index.PAR_file.astype(str)
            EC_index["Loading_cm2"] = EC_index["Loading_cm2"].round(3)
            self.EC_index = EC_index

        if self.SampleCodes.empty:
            SampleCodes = FindExpFolder().LoadSampleCode()
            self.SampleCodes = SampleCodes
        # SampleCodesChar().load

    def set_OVV_exp_type(self):
        if not self.EC_index.empty and self.exp_type:
            PAR_exp_uniq = self.EC_index.PAR_exp.unique()

            PAR_match = [
                parexp
                for parexp in PAR_exp_uniq
                if self.exp_type.split("_")[0] in parexp
            ]
            self.exp_type_match = PAR_match
            # if PAR_match:
            EC_index_exp = self.EC_index.loc[self.EC_index.PAR_exp.isin(PAR_match)]
            self.EC_index_exp = EC_index_exp
            if EC_index_exp.empty:
                _logger.error(f'set_OVV_exp_type "{self.__class__.__name__}" empty')

            self.EC_index_exp_destdirs = EC_index_exp.Dest_dir.unique()

    def get_daily_pickle(self):
        exp_type = self.exp_type
        today = dt.datetime.now().date()
        _result = {"today": today}
        if exp_type:
            daily_pickle_path = FindExpFolder("VERSASTAT").PostDir.joinpath(
                f"{today:%Y-%m-%d}_{exp_type}_{system()}.pkl.compress"
            )
            daily_pkl_options = list(
                FindExpFolder("VERSASTAT").PostDir.rglob(
                    f"*_{exp_type}_{system()}.pkl.compress"
                )
            )
            daily_pkl_options = sorted(
                daily_pkl_options, key=lambda x: x.stat().st_ctime
            )

            _result.update(
                {
                    "daily_path": daily_pickle_path,
                    "_exists": daily_pickle_path.exists(),
                    "daily_options": daily_pkl_options,
                }
            )

            if not daily_pkl_options and not self._reload_raw:
                self._reload_raw = True

            daily_pickle_path_RAW = FindExpFolder("VERSASTAT").PostDir.joinpath(
                f"{today:%Y-%m-%d}_{exp_type}_{system()}_RAW.pkl.compress"
            )
            _pickle_path_RAW_read_in = FindExpFolder("VERSASTAT").PostDir.joinpath(
                f"{exp_type}_{system()}_RAW_read_in.pkl.compress"
            )
            daily_pkl_options_RAW = list(
                FindExpFolder("VERSASTAT").PostDir.rglob(
                    f"*_{exp_type}_{system()}_RAW.pkl.compress"
                )
            )
            daily_pkl_options_RAW = sorted(
                daily_pkl_options_RAW, key=lambda x: x.stat().st_ctime
            )
            _result.update(
                {
                    "daily_path_RAW": daily_pickle_path_RAW,
                    "_raw_exists": daily_pickle_path_RAW.exists(),
                    "daily_options_RAW": daily_pkl_options_RAW,
                    "pkl_path_RAW_read_in": _pickle_path_RAW_read_in,
                }
            )
            if "EIS" in exp_type:

                daily_pkl_options_RAW_WB = list(
                    FindExpFolder("VERSASTAT").PostDir.rglob(
                        f"*_{exp_type}_{system()}_RAW_WB.pkl.compress"
                    )
                )
                daily_pkl_options_RAW_WB = sorted(
                    daily_pkl_options_RAW_WB, key=lambda x: x.stat().st_ctime
                )

                _result.update(
                    {
                        "daily_path_BRUTE": FindExpFolder("VERSASTAT").PostDir.joinpath(
                            f"{today:%Y-%m-%d}_{exp_type}_{system()}_BRUTE.pkl.compress"
                        ),
                        "daily_path_RAW_WB": FindExpFolder(
                            "VERSASTAT"
                        ).PostDir.joinpath(
                            f"{today:%Y-%m-%d}_{exp_type}_{system()}_RAW_WB.pkl.compress"
                        ),
                        "daily_options_RAW_WB": daily_pkl_options_RAW_WB,
                    }
                )

        self.daily_pickle_path = _result

    def load_delegator(self):
        setattr(self, self.exp_type, pd.DataFrame())
        if self._reload:
            if self._reload_raw:
                self.make_raw_pars_from_scratch()
            else:
                self.read_in_daily_raw()

            if hasattr(self, "edit_raw_columns"):
                try:
                    self.edit_raw_columns()
                except Exception as e:
                    _logger.warning(
                        f'edit_raw_columns in load_delegator "{self.__class__.__name__}" {self.exp_type} failed because {e}'
                    )
            self.save_daily_pars()
        else:
            self.read_in_daily_pars()
        try:
            self.extra_stuff_delegator()
        except Exception as e:
            _logger.warning(
                f'extra_stuff_delegator "{self.__class__.__name__}" {self.exp_type} failed because {e}'
            )

    def _check_class_req_functions(self):

        for _f in self._required_funcs:
            if not hasattr(self, _f) and "BaseLoadPars" not in self.__class__.__name__:
                _logger.warning(
                    f'Class "{self.__class__.__name__}" is missing required func: "{_f}"'
                )

    def save_daily_pars(self):
        pars = getattr(self, self.exp_type)
        pars.to_pickle(self.daily_pickle_path["daily_path"])
        _logger.info(
            f'{self.exp_type} len({len(pars)}) OVV to daily pickle: {self.daily_pickle_path.get("daily_path")}'
        )

    def read_in_daily_pars(self):

        if self.daily_pickle_path.get("daily_options"):
            _pars_fp = self.daily_pickle_path.get("daily_options")[-1]

            _logger.info(
                f"start read_in_daily_pars {self.exp_type} pars OVV from  daily {_pars_fp} "
            )
            _pars = pd.read_pickle(_pars_fp)

            try:
                _pars = FileOperations.ChangeRoot_DF(_pars, [], coltype="string")
                setattr(self, self.exp_type, _pars)
                _logger.info(f"Loaded {self.exp_type} pars OVV from  daily {_pars_fp} ")
            except Exception as e:
                _pars = pd.DataFrame()
                _logger.error(
                    f" ERROR in Loaded {self.exp_type} pars OVV from  daily {_pars_fp} {e} "
                )
        else:
            _pars = pd.DataFrame()
            _pars_fp = "options empty list"

        if _pars.empty:
            _logger.error(
                f" ERROR in Loaded {self.exp_type} pars OVV from  daily {_pars_fp}: empty "
            )

    def reload_raw_df_delegator(self):
        _raw_read_fp = self.daily_pickle_path.get("pkl_path_RAW_read_in")
        if _raw_read_fp.exists() and not (self._reload or self._reload_raw):
            _pars_RAW_read_in = pd.read_pickle(_raw_read_fp)
            setattr(self, f"{self.exp_type}_RAW", _pars_RAW_read_in)
        else:
            self.generate_raw_df()
            self.reload_raw_df()
            _pars_RAW_read_in = getattr(self, f"{self.exp_type}_RAW")
            _pars_RAW_read_in.to_pickle(_raw_read_fp)

    def read_in_daily_raw(self):
        _raw_fp = self.daily_pickle_path.get("daily_options_RAW")[-1]
        _pars_RAW = pd.read_pickle(_raw_fp)
        _pars_RAW.sort_values("source_delta_mtime", inplace=True)
        if not "level_0" in _pars_RAW.columns:
            _pars_RAW = _pars_RAW.reset_index()
        setattr(self, f"{self.exp_type}_RAW", _pars_RAW)
        _logger.info(f"Loaded raw df {self.exp_type} from  daily {_raw_fp} ")

    def save_daily_raw(self):
        _pars_RAW = getattr(self, f"{self.exp_type}_RAW")
        _pars_RAW.to_pickle(self.daily_pickle_path.get("daily_path_RAW"))
        _logger.info(
            f'{self.exp_type} OVV to daily pickle: {self.daily_pickle_path.get("daily_path_RAW")}'
        )

    def set_gen_raw_fls(self):
        _par_files = [
            list(self.search_pars_files(d)) for d in self.EC_index_exp_destdirs
        ]
        self._par_files = _par_files
        if not _par_files:
            _logger.warning(f"{self.exp_type} set_gen_raw_fls: list empty ")
        self._par_fls_gen = (a for i in self._par_files for a in i)

    @func_timer_decorator
    def generate_raw_df(self):
        if not hasattr(self, "_par_fls_gen"):
            self.set_gen_raw_fls()

        _pars_lst = list(self.read_in_pars_files(self._par_fls_gen))
        try:
            _pars_RAW = pd.concat(_pars_lst, sort=False)
        except Exception as e:
            _pars_RAW = pd.DataFrame()
            _logger.warning(f"{self.exp_type} generate_raw_df: {e}")
        setattr(self, f"{self.exp_type}_RAW", _pars_RAW)

    @staticmethod
    def get_source_meta(filepath):
        i = filepath
        _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
        _delta_mtime = dt.datetime.now() - _source_mtime
        _meta_res = {
            "sourceFilename": i,
            "source_mtime": _source_mtime,
            "source_delta_mtime": _delta_mtime,
            "sourcebasename": i.stem,
        }
        return _meta_res

    def extra_stuff_delegator(self):
        _extra_funcs = [i for i in self.__dict__.keys() if i.startswith("_extra")]

        for _func in _extra_funcs:
            try:
                func = getattr(self, _func)
                func()
                # self._extra_plotting()
            except Exception as e:
                _logger.info(
                    f"{self.__class__.__name__} Extra stuff failed because {e}"
                )


def _testing():
    tt = EIS_LoadPars(reload=False, reload_raw=False)
    tt._reload_raw
    self = tt
    self.load_delegator()
    self.make_raw_pars_from_scratch()


class EIS_LoadPars(BaseLoadPars):

    col_names = ["File_SpecFit", "File_SpecRaw", "PAR_file"]

    def __init__(
        self,
        EC_index=pd.DataFrame(),
        SampleCodes=pd.DataFrame(),
        exp_type="EIS_pars",
        BRUTE_out=False,
        **kws,
    ):
        self.BRUTE_out = BRUTE_out
        super().__init__(
            EC_index=EC_index, SampleCodes=SampleCodes, exp_type=exp_type, **kws
        )

    def read_in_pars_files(self, _genlist):
        #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
        while True:
            try:
                i = next(_genlist)
                if i.name.endswith("xlsx"):
                    _pp = pd.read_excel(i, index_col=[0])
                elif i.name.endswith("pkl"):
                    _pp = pd.read_pickle(i)
                _pp = FileOperations.ChangeRoot_DF(_pp, [], coltype="string")
                _meta = self.get_source_meta(i)
                _pp = _pp.assign(**_meta)

                yield _pp
            except StopIteration:
                return "all done"
                print("gen empty")

    def search_pars_files(self, _dest_dir):
        return Path(_dest_dir.joinpath("EIS")).rglob(
            f"*_pars_v{FileOperations.EIS_version}.xlsx"
        )

    @func_timer_decorator
    def make_raw_pars_from_scratch(self):
        _logger.info(
            f'Reloading raw extra steps "{self.__class__.__name__}" {self.exp_type}'
        )
        self.reload_raw_df_delegator()
        self._load_WB_delegator()
        self._merge_WB_pars_raw()
        self._raw_finish_edit_columns()
        self.save_daily_raw()

    def reload_raw_df_delegator(self):
        _raw_read_fp = self.daily_pickle_path.get("pkl_path_RAW_read_in")
        if _raw_read_fp.exists() and not (self._reload or self._reload_raw):
            EIS_pars_RAW_read_in = pd.read_pickle(_raw_read_fp)
            setattr(self, f"{self.exp_type}_RAW", EIS_pars_RAW_read_in)
        else:
            self.generate_raw_df()
            self.reload_raw_df()
            EIS_pars_RAW_read_in = getattr(self, f"{self.exp_type}_RAW")
            EIS_pars_RAW_read_in.to_pickle(_raw_read_fp)

    def reload_raw_df(self):
        _pars_RAW = getattr(self, f"{self.exp_type}_RAW")
        _pars_RAW.sort_values("source_delta_mtime", inplace=True)
        _pars_RAW = _pars_RAW.reset_index()
        setattr(self, f"{self.exp_type}_RAW", _pars_RAW)
        self._raw_extra_steps()
        _logger.info(f'Reloading "{self.__class__.__name__}" {self.exp_type}')
        # self.EIS_pars_RAW = EIS_pars_RAW

    def _raw_extra_steps(self):
        _logger.info(
            f'Reloading raw extra steps "{self.__class__.__name__}" {self.exp_type}'
        )
        EIS_pars_all = getattr(self, f"{self.exp_type}_RAW")
        float_cols = set(
            [
                a
                for i in EIS_pars_all.lmfit_var_names.unique()
                if type(i) == str and not "(" in i
                for a in i.split(", ")
            ]
        )
        float_cols.update(
            set(
                [a for i in float_cols for a in EIS_pars_all.columns if a.startswith(i)]
            )
        )
        EIS_pars_all[list(float_cols)] = EIS_pars_all[list(float_cols)].fillna(0)
        # EIS_pars_all[list(float_cols)] = EIS_pars_all[list(float_cols)].astype(float)
        obj_flt_cols = [
            i
            for i in EIS_pars_all.columns
            if str(EIS_pars_all[i].dtype) == "object" and i in float_cols
        ]
        EIS_pars_all[obj_flt_cols] = EIS_pars_all[obj_flt_cols].replace("", 0)
        EIS_pars_all[list(float_cols)] = EIS_pars_all[list(float_cols)].astype(float)
        wrong_fls = [
            EIS_pars_all.loc[EIS_pars_all[i].astype(str).str.contains("Parameter")]
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
        setattr(self, f"{self.exp_type}_RAW", EIS_pars_all)

    def _load_WB_delegator(self):
        daily_options_WB = self.daily_pickle_path.get("daily_options_RAW_WB")
        if daily_options_WB:
            _WB_RAW_daily_path = daily_options_WB[-1]
            if _WB_RAW_daily_path.exists() and not (self._reload or self._reload_raw):
                _EIS_WB_pars_all = pd.read_pickle(_WB_RAW_daily_path)
                setattr(self, f"{self.exp_type}_WB", _EIS_WB_pars_all)
            else:
                self.reload_raw_WB_df()
        else:
            self.reload_raw_WB_df()

    def reload_raw_WB_df(self):
        _logger.info(f'Reloading "{self.__class__.__name__}" {self.exp_type} WB')
        _EIS_WB_files = [
            list(Path(d.joinpath("EIS/lin_Warburg")).rglob(f"lin_Warburg*.pkl"))
            for d in self.EC_index_exp_destdirs
        ]
        self._EIS_WB_files = _EIS_WB_files
        self._EIS_WB_fls = (a for i in _EIS_WB_files for a in i)

        _WB_lst = list(self.read_in_pars_files(self._EIS_WB_fls))
        _EIS_WB_pars_all = pd.concat(_WB_lst, sort=False, ignore_index=True)
        setattr(self, f"{self.exp_type}_WB", _EIS_WB_pars_all)
        _EIS_WB_pars_all.to_pickle(self.daily_pickle_path.get("daily_path_RAW_WB"))

    def _merge_WB_pars_raw(self):
        _EIS_WB_pars_all = getattr(self, f"{self.exp_type}_WB")
        EIS_pars_all = getattr(self, f"{self.exp_type}_RAW")
        _diffcols = set(EIS_pars_all.columns).difference(_EIS_WB_pars_all.columns)
        _mcols = [
            i
            for i in set(EIS_pars_all.columns).intersection(_EIS_WB_pars_all.columns)
            if i
            not in [
                "sourceFilename",
                "source_mtime",
                "source_delta_mtime",
                "sourcebasename",
            ]
        ]
        _dtype_mismatch = [
            (i, EIS_pars_all[i].dtype, _EIS_WB_pars_all[i].dtype)
            for i in _mcols
            if EIS_pars_all[i].dtype != _EIS_WB_pars_all[i].dtype
        ]
        if _dtype_mismatch:
            _excl = []
            for i in _dtype_mismatch:
                try:
                    _EIS_WB_pars_all[i[0]] = _EIS_WB_pars_all[i[0]].astype(i[1])
                except Exception as e:
                    _excl.append(i[0])
                    print(i, "\n", e)
            _mcols = [i for i in _mcols if i not in _excl]
            # EIS_pars_all[i[0]] = EIS_pars_all[i[0]].astype(i[2])
        _merge = pd.merge(
            EIS_pars_all, _EIS_WB_pars_all, on=_mcols, how="left", suffixes=("", "_WB")
        )
        if not _merge.empty:
            return _merge
        else:
            print("WB merge was empty")
            return EIS_pars_all
        setattr(self, f"{self.exp_type}_RAW", EIS_pars_all)

    def _raw_finish_edit_columns(self):
        # EIS_pars_all = self._merge_WB_pars_raw(EIS_pars_all)
        EIS_pars_all = getattr(self, f"{self.exp_type}_RAW")
        EIS_pars_all = EIS_pars_all.assign(
            **{
                "EIS_fake": [
                    "fakeZmean" in Path(i).name
                    for i in EIS_pars_all.PAR_file.to_numpy()
                ]
            }
        )
        _not_in_index = EIS_pars_all.loc[
            (
                ~(EIS_pars_all.PAR_file.isin(self.EC_index.PAR_file.values))
                & (~EIS_pars_all.EIS_fake == True)
            )
        ]
        CleanUpCrew(list_of_files=_not_in_index.sourceFilename.unique(), delete=True)
        EIS_pars_all = EIS_pars_all.iloc[
            ~(EIS_pars_all.index.isin(_not_in_index.index))
        ]
        EIS_pars_all = Load_from_Indexes.test_update_from_index(
            EIS_pars_all, self.EC_index
        )
        setattr(self, f"{self.exp_type}_RAW", EIS_pars_all)

    def edit_raw_columns(self):
        EIS_pars_all = getattr(self, f"{self.exp_type}_RAW")
        # EIS_pars_RAW = self._raw_extra_steps(EIS_pars_RAW)
        E_dc_RHE_cols = [
            (np.round(i, 3), np.round(i, 3) * 1e3) for i in EIS_pars_all[EvRHE].values
        ]
        EIS_pars_all = EIS_pars_all.assign(
            **{
                "E_dc_RHE": [i[0] for i in E_dc_RHE_cols],
                "E_dc_RHE_mV": [i[1] for i in E_dc_RHE_cols],
            }
        )
        EIS_pars_recent = EIS_pars_all.loc[
            (EIS_pars_all.source_mtime > pd.Timestamp(dt.date(2020, 11, 25)))
            & (EIS_pars_all.PAR_file.str.contains("None") == False)
        ]
        EIS_pars_undup = EIS_pars_recent.dropna(subset=self.col_names).drop_duplicates(
            keep="first"
        )
        # === POST EDITING OF LOADED PARS ===
        EIS_pars_undup = EIS_pars_undup.assign(
            **{"Loading_cm2": EIS_pars_undup["Loading_cm2"].round(3)}
        )
        EIS_pars_undup = post_helper.make_uniform_EvRHE(EIS_pars_undup)
        EIS_pars_undup = CollectPostOVV.MatchECconditions(EIS_pars_undup)
        #            EIS_pars_undup = Load_from_Indexes.add_missing_ECindex_cols(EC_index, EIS_pars_undup)
        _oc_OVV = list(EIS_pars_undup.columns.intersection(self.EC_index_exp.columns))
        if not set(self.EC_index_exp.groupby(_oc_OVV).groups.keys()).intersection(
            EIS_pars_undup.groupby(_oc_OVV).groups.keys()
        ):
            _drpcols = [
                a
                for a in EIS_pars_undup.columns
                if (
                    a in [i for i in _oc_OVV if i not in "PAR_file"]
                    or "_".join(a.split("_")[0:-1])
                    in [i for i in _oc_OVV if i not in "PAR_file"]
                )
            ]
            #                EIS_pars_undup.drop(columns =_drpcols)
            EIS_pars_undup = Load_from_Indexes.add_missing_ECindex_cols(
                self.EC_index, EIS_pars_undup.drop(columns=_drpcols)
            )
        #            EIS_pars_undup = pd.merge(EIS_pars_undup,EIS_OVV,on=_oc_OVV, how='left')
        _oc_SC = list(EIS_pars_undup.columns.intersection(self.SampleCodes.columns))
        EIS_pars_undup = pd.merge(
            EIS_pars_undup, self.SampleCodes, how="left", on=_oc_SC
        )
        EIS_pars_BRUTE = EIS_pars_undup.loc[
            (EIS_pars_undup.BRUTE_FIT == 1) | (EIS_pars_undup.FINAL_FIT == 0)
        ]
        if self.BRUTE_out:
            EIS_pars_BRUTE.to_pickle(eis_daily["daily_path_BRUTE"])

        EIS_pars = EIS_pars_undup.loc[(EIS_pars_undup.FINAL_FIT == 1)]
        EIS_pars = EIS_extra_methods.add_best_model_per_spectrum(EIS_pars)

        setattr(self, self.exp_type, EIS_pars)

    # def extra_stuff_delegator(self):
    #     try:
    #         self._extra_best_models()
    #         self._extra_plotting()
    #     except Exception as e:
    #         _logger.info(f'{self.__class__.__name__} Extra stuff failed because {e}')

    def _extra_best_models(self):
        _err_type = "lmfit_MSE"
        _filter = "(EIS_pars.lmfit_MSE < 65E4) & (EIS_pars.Rct < 2E3) & (EIS_pars.Rct > 2E-2) \
                                   & (EIS_pars.Rs > 0.01)  & (EIS_pars.Rs < 200) & (EIS_pars.Cdlp < 0.075)\
                                   & (EIS_pars.lmfit_redchi < 1E3)  & (EIS_pars.Aw < 10E3) & (EIS_pars.Aw > 10E-2)\
                                   & (EIS_pars.Qad < 1) & (EIS_pars.tau < 1E3)"
        _filter += '& (EIS_pars.SampleID.str.contains("JOS1|JOS2|JOS3|JOS4|JOS5"))'
        _filter += "& (EIS_pars.EIS_fake == False)"
        _grps = ["Model_EEC", "Gas", "lmfit_var_names"][0:2]

        EIS_pars = self.EIS_pars

        best_models = (
            EIS_pars.loc[eval(_filter)]
            .dropna(axis=0, subset=[_err_type])
            .groupby(_grps)[_err_type]
            .agg(["count", "mean", "std"])
            .sort_values(["Gas", "mean"], ascending=True)
        )
        print(best_models)

        keep_models = (
            best_models.loc[(best_models["count"] > 5) & (best_models["std"] > 0)]
            .index.get_level_values(0)
            .unique()
        )
        EIS_pars = EIS_pars.loc[EIS_pars.Model_EEC.isin(keep_models)]

        if hasattr(EIS_pars, "best_mod_name"):
            # EIS_best_mods = EIS_pars.loc[EIS_pars.Model_EEC_name.isin([i for i in EIS_pars.best_mod_name.unique() if not pd.isna(i)])]
            EIS_best_mods = EIS_pars.loc[
                EIS_pars.index.isin(
                    [i for i in EIS_pars.best_mod_n.unique() if not pd.isna(i)]
                )
            ]
            self.EIS_pars_best_mods = EIS_best_mods
            _agg = (
                EIS_best_mods.dropna(subset=[_err_type])
                .groupby(_grps + ["E_RHE"])[_err_type]
                .agg(["count", "mean", "std"])
            )

            _agg_best = _agg.loc[_agg["count"] > 3].sort_values(
                ["Gas", "E_RHE", "mean"], ascending=True
            )

    def _extra_plotting(self):
        if hasattr(self, "EIS_pars_best_mods"):
            self.EIS_pars_best_mods.query("pH < 15").plot(
                y="Qad",
                x="E_RHE",
                c="pH",
                colormap="rainbow_r",
                kind="scatter",
                ylim=(0, 0.05),
            )
            self.EIS_pars_best_mods.query("pH < 15").plot(
                y="Rs",
                x="E_RHE",
                c="pH",
                colormap="rainbow_r",
                kind="scatter",
                ylim=(0, 80),
            )
            self.EIS_pars_best_mods.query("pH < 15").plot(
                y="Rs",
                x="R_ion",
                c="E_RHE",
                colormap="rainbow_r",
                kind="scatter",
                ylim=(0, 80),
                xlim=(0.1, 2e3),
                logx=True,
            )


def _testing():
    t2 = ORR_LoadPars(reload=True, reload_raw=True)
    tf2 = ORR_LoadPars(reload=False, reload_raw=False)
    t2._reload_raw
    self = tf2
    self.load_delegator()
    self.make_raw_pars_from_scratch()


class ORR_LoadPars(BaseLoadPars):

    read_types = ["ORR_pars", "KL_pars"]

    def __init__(
        self,
        EC_index=pd.DataFrame(),
        SampleCodes=pd.DataFrame(),
        exp_type="ORR_pars",
        BRUTE_out=False,
        **kws,
    ):
        self.BRUTE_out = BRUTE_out
        super().__init__(
            EC_index=EC_index, SampleCodes=SampleCodes, exp_type=exp_type, **kws
        )

    def read_in_pars_files(self, _genlist):
        #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
        while True:
            try:
                i = next(_genlist)
                # _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
                # _delta_mtime = dt.datetime.now() - _source_mtime
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

                _meta = self.get_source_meta(i)
                _meta.update({"source_type": _type})

                if _type in self.read_types:
                    _pp = pd.read_excel(i, index_col=[0])
                    _pp = FileOperations.ChangeRoot_DF(_pp, [], coltype="string")
                    _pp = _pp.assign(**_meta)

                else:
                    _pp = pd.DataFrame(_meta, index=[0])
                # _meta.update({'DF' : _pp})
                yield _pp
            except StopIteration:
                return "all done"
                print("gen empty")

    @func_timer_decorator
    def make_raw_pars_from_scratch(self):
        _logger.info(
            f'Reloading raw extra steps "{self.__class__.__name__}" {self.exp_type}'
        )
        self.reload_raw_df_delegator()
        if hasattr(self, "_raw_finish_edit_columns"):
            self._raw_finish_edit_columns()
        self.save_daily_raw()

    def search_pars_files(self, dest_dir):
        return Path(dest_dir.joinpath(f"ORR_v{FileOperations.version}")).rglob("*xlsx")

    def reload_raw_df(self):
        _pars_RAW = getattr(self, f"{self.exp_type}_RAW")
        _pars_RAW.sort_values("source_delta_mtime", inplace=True)
        _pars_RAW = _pars_RAW.reset_index()
        setattr(self, f"{self.exp_type}_RAW", _pars_RAW)
        # self._raw_extra_steps()
        _logger.info(f'Reloading "{self.__class__.__name__}" {self.exp_type}')
        # self.EIS_pars_RAW = EIS_pars_RAW

    def edit_raw_columns(self):
        ### Fixing the pars after loading...
        # TODO : taking out duplicates based on time_since_run....
        ORR_pars_char = getattr(self, f"{self.exp_type}_RAW")
        #         Load_na = ORR_pars_char.loc[(ORR_pars_char.Loading_cm2.isna()) & (ORR_pars_char.PAR_file.isna() == False)]
        #         if not Load_na.empty:
        #             Load_na_missingvalues =[(n,*GetSampleID.ink_loading_from_filename(i.PAR_file)) for n,i in Load_na.iterrows()]
        #             Load_na_vals = pd.DataFrame(Load_na_missingvalues).rename(columns={1 : 'Loading_name',2 : 'Loading_cm2'}).set_index([0])
        #             ORR_pars_char.Loading_cm2.fillna(value=Load_na_vals.Loading_cm2,inplace=True)
        # #            ORR_char_merge_cols = [i for i in ORR_pars.columns if i in SampleCodes.columns]
        ORR_pars_char = ORR_pars_char.drop(
            columns=[i for i in ORR_pars_char.columns if "Unnamed" in i]
        )
        if not ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna()].empty:
            _loading_cols = ["Loading_cm2", "Loading_name", "Loading_date"]
            ORR_pars_char = ORR_pars_char.drop(columns=_loading_cols)
            ORR_pars_char = pd.merge(
                ORR_pars_char,
                self.EC_index[["PAR_file"] + _loading_cols],
                on="PAR_file",
                how="left",
            )
            ORR_pars_char.Loading_cm2 = ORR_pars_char.Loading_cm2.fillna(
                value=0.379
            )  # fillna for Loading_cm2
        ORR_pars_char.Loading_cm2 = ORR_pars_char.Loading_cm2.round(3)

        if ORR_pars_char.postAST.dropna().empty:
            ORR_pars_char = ORR_pars_char.drop(columns="postAST")
            #                 _int = list(set(ORR_pars_char.columns).intersection(set(EC_index.columns)))
            ORR_pars_char = pd.merge(
                ORR_pars_char,
                self.EC_index[["PAR_file", "postAST"]],
                on="PAR_file",
                suffixes=("", ""),
            )

        ORR_pars_char = make_uniform_RPM_DAC(ORR_pars_char)
        setattr(self, f"{self.exp_type}", ORR_pars_char)

    # def extra_stuff_delegator(self):
    #     try:
    #         self._extra_plotting()
    #     except Exception as e:
    #         _logger.info(f'{self.__class__.__name__} Extra stuff failed because {e}')

    def _extra_plotting(self):
        ORR_pars_char = getattr(self, f"{self.exp_type}")
        for swp, swgrp in ORR_pars_char.query("(pH < 14) & (RPM_DAC > 900)").groupby(
            "Sweep_Type"
        ):
            fig, (ax1, ax2) = plt.subplots(figsize=(10, 4), ncols=2)
            #                plt.figure()
            swgrp.plot(
                y="ORR_Jkin_min_750",
                x="ORR_E_onset",
                c="pH",
                title=f"{swp}",
                kind="scatter",
                logy=True,
                colormap="rainbow_r",
                ylim=[0.1, 50],
                xlim=(0.5, 1),
                ax=ax1,
            )
            ax1.set_xlabel("E onset / mV_RHE")
            swgrp.plot(
                y="ORR_Frac_H2O2_600",
                x="ORR_E_onset",
                c="pH",
                title=f"{swp}",
                kind="scatter",
                logy=True,
                colormap="rainbow_r",
                ylim=[0.1, 100],
                xlim=(0.5, 1),
                ax=ax2,
            )
            # ax2.set_xlabel('E onset / mV_RHE')
            plt.suptitle("ORR with E_onset")
            plt.show()

            fig, (ax1, ax2) = plt.subplots(figsize=(10, 4), ncols=2)
            swgrp.plot(
                y="ORR_E_onset",
                x="N2_BG_lin_slope",
                c="pH",
                title=f"{swp}",
                kind="scatter",
                logy=True,
                logx=True,
                colormap="rainbow_r",
                xlim=[0.01, 4],
                ylim=(0.5, 1),
                ax=ax1,
            )
            swgrp.plot(
                y="ORR_Jkin_min_750",
                x="N2_BG_lin_slope",
                c="pH",
                title=f"{swp}",
                kind="scatter",
                logy=True,
                logx=True,
                colormap="rainbow_r",
                xlim=[0.01, 4],
                ylim=(0.001, 50),
                ax=ax2,
            )
            # ax2.set_xlabel('E onset / mV_RHE')
            plt.suptitle("ORR with N2_BG lin slope")
            plt.show()

        plt.close()


def _N2_testing():
    n2 = N2_LoadPars(reload=True, reload_raw=True)
    n2r = N2_LoadPars(reload=True, reload_raw=False)


class N2_LoadPars(BaseLoadPars):
    def __init__(
        self,
        EC_index=pd.DataFrame(),
        SampleCodes=pd.DataFrame(),
        exp_type="",
        BRUTE_out=False,
        **kws,
    ):
        self.BRUTE_out = BRUTE_out
        super().__init__(
            EC_index=EC_index, SampleCodes=SampleCodes, exp_type=exp_type, **kws
        )

    @func_timer_decorator
    def make_raw_pars_from_scratch(self):
        _logger.info(
            f'Reloading raw extra steps "{self.__class__.__name__}" {self.exp_type}'
        )
        self.reload_raw_df_delegator()
        if hasattr(self, "_raw_finish_edit_columns"):
            self._raw_finish_edit_columns()
        self.save_daily_raw()

    def _old(self):
        IndexOVV_N2_pars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "N2Cdl_pars_IndexOVV_v{0}.pkl.compress".format(FileOperations.version)
        )
        n2_daily = get_daily_pickle(exp_type="N2_all")

        if n2_daily.get("_exists", False) and reload != True:
            #            Cdl_pars_char = pd.read_excel(IndexOVV_N2_pars_fn,index_col=[0])
            Cdl_pars_char = pd.read_pickle(n2_daily.get("daily_path"))
            Cdl_pars_char = FileOperations.ChangeRoot_DF(
                Cdl_pars_char, [], coltype="string"
            )
        else:
            # @@ Check POST_AST status from OVV and PRM
            _logger.info(
                f'START reloading N2_pars OVV from daily {n2_daily["today"]:%Y-%m-%d}'
            )

    #            EC_index = ECRunOVV(load=1).index
    #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
    #            EC_index = FileOperations.ChangeRoot_DF(OnlyRecentMissingOVV,[])
    #            OnlyRecentMissingOVV.PAR_file = OnlyRecentMissingOVV.PAR_file.astype(str)
    #            OnlyRecentMissingOVV['Loading_cm2'] = OnlyRecentMissingOVV['Loading_cm2'].round(3)
    #            SampleCodes = SampleCodesChar().load
    # EC_index, SampleCodes = Load_from_Indexes.get_EC_index()
    # def read_df(_par_fls, ):
    #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
    def search_pars_files(self, destdir):
        return Path(destdir.joinpath(f"N2_scans_v{FileOperations.version}")).rglob(
            "*.xlsx"
        )

    def read_in_pars_files(self, _genlist, read_types=["Cdl_data", "Cdl_pars"]):
        while True:
            try:
                i = next(_genlist)

                _i_stem = i.stem

                _meta = self.get_source_meta(i)
                if _i_stem.endswith("_BG"):
                    _N2_type = "BG"
                else:
                    if _i_stem.startswith("CV_"):
                        _N2_type = "CV"
                        if _i_stem.endswith(f"_first_v{FileOperations.version}"):
                            _N2_type = "CV_first"
                    #                                if not 'Scan Rate' in _pp.columns:
                    #                                    'N2_CV_raw =  N2_CV_raw.assign(**{'ScanRate' : [i.split(f'_v{FileOperations.version}')[0].split('_')[-1] for i in N2_CV_raw.basename.to_numpy()]})
                    elif _i_stem.startswith("Cdl_data_"):
                        _N2_type = "Cdl_data"
                    elif _i_stem.startswith("Cdl_pars"):
                        _N2_type = "Cdl_pars"
                    else:
                        _N2_type = "N2_unknown"
                _meta.update({"N2_type": _N2_type})

                if _N2_type in read_types:
                    _pp = pd.read_excel(i, index_col=[0])
                    _pp = FileOperations.ChangeRoot_DF(_pp, [], coltype="string")
                    _pp = _pp.assign(**_meta)
                else:
                    _pp = pd.DataFrame(_meta, index=[0])
                # _meta.update({'DF' :  _pp})
                yield _pp
            except StopIteration:
                return "all done"
                print("gen empty")

    def reload_raw_df(self):
        _pars_RAW = getattr(self, f"{self.exp_type}_RAW")
        if not _pars_RAW.empty:
            _pars_RAW.sort_values("source_delta_mtime", inplace=True)
            _pars_RAW = _pars_RAW.reset_index()
        setattr(self, f"{self.exp_type}_RAW", _pars_RAW)
        _logger.info(
            f'Reloading "{self.__class__.__name__}" {self.exp_type} len({len(_pars_RAW)}'
        )

    def _old_stuff():
        if n2_daily.get("_raw_exists", False) and use_daily is True:
            N2_pars_all = pd.read_pickle(n2_daily.get("daily_path_RAW"))
        elif n2_daily.get("daily_options_RAW", False) and use_daily is True:
            if n2_daily.get("daily_options_RAW")[-1]:
                N2_pars_all = pd.read_pickle(n2_daily.get("daily_options_RAW")[-1])
        else:  # Construct new N2 pars ovv from reading in files
            N2_OVV = EC_index.loc[EC_index.PAR_exp == "N2_act"]
            _par_files = [
                list(Path(d.joinpath("N2_scans_v30")).rglob("*.xlsx"))
                for d in N2_OVV.Dest_dir.unique()
            ]
            _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
            _par_reads = read_df(_par_fls, read_types=["Cdl_data", "Cdl_pars"])
            N2_pars_all = pd.concat([i["DF"] for i in _par_reads], sort=False)

            for n, gr in N2_pars_all.groupby("PAR_file"):
                print(
                    n,
                    f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}',
                    ",".join(gr.N2_type.unique()),
                )

            N2_pars_all, _missing_index = Load_from_Indexes.check_missing_ECindex(
                EC_index, N2_pars_all, clean_up=True
            )
            N2_pars_all.to_pickle(n2_daily["daily_path_RAW"])

    def _extra_pivot_CV(self):
        N2_type_grps = N2_pars_all.groupby("N2_type")

        if "CV" in N2_type_grps.groups.keys():
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
    def edit_raw_columns(self):

        N2Cdl_pars_raw = getattr(self, f"{self.exp_type}_RAW")

        N2_type_grps = N2Cdl_pars_raw.groupby("N2_type")
        N2Cdl_pars_raw = N2_type_grps.get_group("Cdl_pars").dropna(axis=1, how="all")
        N2Cdl_pars_raw.drop_duplicates(
            subset=N2Cdl_pars_raw.columns[0:19], keep="first", inplace=True
        )
        N2Cdl_pars_raw = FileOperations.ChangeRoot_DF(
            N2Cdl_pars_raw, [], coltype="string"
        )
        Cdl_pars = post_helper.make_uniform_EvRHE(N2Cdl_pars_raw)
        Cdl_pars.drop_duplicates(subset=Cdl_pars.columns[0:19], inplace=True)
        # Cdl_pars_merge_cols = [i for i in Cdl_pars.columns if i in SampleCodes.columns and not 'Unnamed' in i]
        # Cdl_pars_char = pd.merge(Cdl_pars,SampleCodes,on=Cdl_pars_merge_cols,how='left')
        # Cdl_pars_char.drop_duplicates(subset=Cdl_pars_char.columns[0:19],inplace=True)
        _int = list(set(Cdl_pars.columns).intersection(set(self.EC_index.columns)))
        if Cdl_pars.postAST.dropna().empty and len(self.EC_index.columns) != len(_int):
            Cdl_pars = Cdl_pars.drop(columns="postAST")
            #                 _int = list(set(Cdl_pars_char.columns).intersection(set(EC_index.columns)))
            Cdl_pars = pd.merge(
                Cdl_pars,
                self.EC_index[["PAR_file", "postAST"]],
                on="PAR_file",
                suffixes=("", ""),
            )

        Cdl_pars = Load_from_Indexes.add_missing_ECindex_cols(self.EC_index, Cdl_pars)
        setattr(self, f"{self.exp_type}", Cdl_pars)

    def _extra_xls_out(self):
        if xls_out:
            new_N2_pars_char_target = FileOperations.CompareHashDFexport(
                Cdl_pars_char, IndexOVV_N2_pars_fn
            )
            _logger.info(
                "PostEC Cdl N2 CVs re-indexed and saved: {0}".format(
                    new_N2_pars_char_target
                )
            )
        Cdl_pars_char.to_pickle(IndexOVV_N2_pars_fn)

    def _extra_plotting(self):
        try:
            Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)').plot(
                y="Cdl",
                x="E_RHE",
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
                acidCdl = Cdl_pars_char.query(
                    '(Sweep_Type_N2 == "cathodic") & (pH < 7)'
                )
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
        except Exception as e:
            _logger.warning(f"PostEC Cdl N2 CVs extra plotting fail:\n{e}")


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
            _logger.info("PostEC loaded IndexOVV from recent: {0}".format(IndexOVV_fn))
        else:
            _logger.info(
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
                _logger.info(
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
class Load_from_Indexes:
    """This class loads the parameters of Electrochemical Data files and merge it with the Overview"""

    SampleCodes = FindExpFolder().LoadSampleCode()
    #    EC_label_cols = ['SampleID','pH','Electrolyte','Loading_cm2','postAST','PAR_date_day']
    EC_label_cols = [
        "PAR_file",
        "SampleID",
        "postAST",
        "Loading_cm2",
        "Electrolyte",
        "pH",
        "Gas",
        "RPM_DAC",
        "E_RHE",
    ]
    PostDestDir = FindExpFolder("VERSASTAT").PostDir

    def __init__(self, **kwargs):
        if "reload" in kwargs:
            #            self.postOVVout  = CollectPostOVV.LoadPostOVV(kwargs['reload'])
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
        EC_index = ECRunOVV(load=1).EC_index
        #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
        EC_index = FileOperations.ChangeRoot_DF(EC_index, [])
        EC_index.PAR_file = EC_index.PAR_file.astype(str)
        EC_index["Loading_cm2"] = EC_index["Loading_cm2"].round(3)
        SampleCodes = FindExpFolder().LoadSampleCode()
        # SampleCodesChar().load
        return EC_index, SampleCodes

    @staticmethod
    def check_missing_ECindex(OnlyRecentMissingOVV, DF_pars, clean_up=False):
        not_in_index = DF_pars.loc[
            ~DF_pars.PAR_file.isin(OnlyRecentMissingOVV.PAR_file.values)
        ]
        CleanUpCrew(list_of_files=not_in_index.sourceFilename.unique(), delete=clean_up)
        return (
            DF_pars.loc[DF_pars.PAR_file.isin(OnlyRecentMissingOVV.PAR_file.values)],
            not_in_index,
        )

    @staticmethod
    def add_missing_ECindex_cols(EC_index, DF):
        if list(EC_index.columns.difference(DF.columns)):
            DF = pd.merge(
                DF,
                EC_index[["PAR_file"] + list(EC_index.columns.difference(DF.columns))],
                on="PAR_file",
                how="left",
            )
        return DF

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

    def test_update_from_index(pars, EC_index):

        _olap = pars.columns.intersection(EC_index.columns)
        _olap_minus = [i for i in _olap if not "PAR_file" == i]
        _mtime = [i for i in pars.columns if i.endswith("delta_mtime")]
        if _mtime:
            _idx = pars[_mtime[0]].idxmin()
        else:
            _idx = 0

        _ECidx = (
            EC_index.loc[EC_index.PAR_file == pars.iloc[_idx].PAR_file][_olap]
            .iloc[0]
            .to_dict()
        )
        _prsx = pars.iloc[_idx][_olap].to_dict()
        _check = {
            key: {"pars": val, "EC_index": _ECidx.get(key, "xx")}
            for key, val in _prsx.items()
            if _ECidx.get(key, "xx") != val
        }
        _pars_bad = False
        if _check:
            _pars_bad = any(
                "error" in str(i) for i in [i["pars"] for i in _check.values()]
            )
        if _pars_bad:
            _logger.info(f"Overwriting columns in Pars from EC_index")
            _new_pars = pd.merge(
                pars[[i for i in pars.columns if i not in _olap_minus]],
                EC_index[_olap],
                on="PAR_file",
                how="left",
            )
        else:
            _new_pars = pars
        return _new_pars

    @staticmethod
    def EIS_pars_OVV(
        reload=False,
        extra_plotting=False,
        xls_out=False,
        BRUTE_out=False,
        use_daily=True,
        use_latest=False,
        **kwargs,
    ):

        #        IndexOVV_EISpars_fn_xls = PostDestDir.joinpath('EIS_pars_IndexOVV_v{0}.xlsx'.format(FileOperations.version))
        #        IndexOVV_EISpars_fn = PostDestDir.joinpath('EIS_pars_IndexOVV_v{0}.pkl.compress'.format(FileOperations.version))
        #        PostDestDir = Load_from_Indexes.PostDestDir
        #        FindExpFolder('VERSASTAT').PostDir
        eis_daily = get_daily_pickle(exp_type="EIS_pars")

        #        today = dt.datetime.now().date()
        #        eis_daily_pickle_path = PostDestDir.joinpath(f'{today.year}-{today.month}-{today.day}_EIS_pars_{system()}.pkl.compress')
        #        eis_daily_pickle_path_RAW = PostDestDir.joinpath(f'{today.year}-{today.month}-{today.day}_EIS_pars_{system()}_RAW.pkl.compress')

        if eis_daily.get("_exists", False) and not reload and use_daily:
            EIS_pars = pd.read_pickle(eis_daily.get("daily_path"))
            EIS_pars = FileOperations.ChangeRoot_DF(EIS_pars, [], coltype="string")
            _logger.info(
                f'Loaded EIS_pars OVV from  daily {eis_daily["today"]} pickle: {eis_daily.get("daily_path","")}'
            )
        elif (
            eis_daily.get("daily_options", [])
            and not reload
            and (use_latest or use_daily)
        ):
            EIS_pars = pd.read_pickle(eis_daily.get("daily_options")[-1])
            EIS_pars = FileOperations.ChangeRoot_DF(EIS_pars, [], coltype="string")
            _logger.info(
                f'Loaded EIS_pars OVV from  daily {eis_daily.get("daily_options")[-1]} '
            )
        else:
            # @@ Read EIS pars files and extend with columns from Samples
            # try other way::   idx_files_EIS = [list(Path(i).rglob('**/EIS/*pars_v20.xlsx')) for i in OnlyRecentMissingOVV.Dest_dir.unique() if list(Path(i).rglob('**/EIS/*pars_v20.xlsx'))]
            _logger.info(
                f'START reloading EIS_pars OVV from  daily {eis_daily["today"]}'
            )
            #            OnlyRecentMissingOVV = ECRunOVV(load=1).index
            ##            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
            #            OnlyRecentMissingOVV = FileOperations.ChangeRoot_DF(OnlyRecentMissingOVV,[])
            #            OnlyRecentMissingOVV.PAR_file = OnlyRecentMissingOVV.PAR_file.astype(str)
            #            OnlyRecentMissingOVV['Loading_cm2'] = OnlyRecentMissingOVV['Loading_cm2'].round(3)
            #            SampleCodes = SampleCodesChar().load
            EC_index, SampleCodes = Load_from_Indexes.get_EC_index()

            def read_df(_par_fls):
                #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
                while True:
                    try:
                        i = next(_par_fls)
                        if i.name.endswith("xlsx"):
                            _pp = pd.read_excel(i, index_col=[0])
                        elif i.name.endswith("pkl"):
                            _pp = pd.read_pickle(i)
                        _pp = FileOperations.ChangeRoot_DF(_pp, [], coltype="string")
                        _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
                        _delta_mtime = dt.datetime.now() - _source_mtime
                        _pp = _pp.assign(
                            **{
                                "sourceFilename": i,
                                "source_mtime": _source_mtime,
                                "source_delta_mtime": _delta_mtime,
                                "sourcebasename": i.stem,
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
            EIS_OVV = EC_index.loc[EC_index.PAR_exp == "EIS"]

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
            _EIS_WB_files = [
                list(Path(d.joinpath("EIS/lin_Warburg")).rglob(f"lin_Warburg*.pkl"))
                for d in EIS_OVV.Dest_dir.unique()
            ]
            _EIS_WB_fls = (a for i in _EIS_WB_files for a in i)
            _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
            #            tt = (i for i in _par_fls if bn in i.name)

            #            __ttp = list(read_df(tt, col_names))
            if eis_daily.get("_raw_exists", False) and use_daily == True:
                EIS_pars_all = pd.read_pickle(eis_daily.get("daily_path_RAW"))
            elif (
                not eis_daily.get("_raw_exists", False)
                and use_daily == True
                and eis_daily.get("daily_options_RAW")
            ):
                EIS_pars_all = pd.read_pickle(eis_daily.get("daily_options_RAW")[-1])
            else:
                _pars_lst = list(read_df(_par_fls))
                EIS_pars_RAW = pd.concat(_pars_lst, sort=False)
                EIS_pars_RAW.sort_values("source_delta_mtime", inplace=True)
                EIS_pars_RAW = EIS_pars_RAW.reset_index()
                EIS_pars_all = EIS_pars_RAW
                float_cols = set(
                    [
                        a
                        for i in EIS_pars_all.lmfit_var_names.unique()
                        if type(i) == str and not "(" in i
                        for a in i.split(", ")
                    ]
                )
                float_cols.update(
                    set(
                        [
                            a
                            for i in float_cols
                            for a in EIS_pars_all.columns
                            if a.startswith(i)
                        ]
                    )
                )
                EIS_pars_all[list(float_cols)] = EIS_pars_all[list(float_cols)].fillna(
                    0
                )
                # EIS_pars_all[list(float_cols)] = EIS_pars_all[list(float_cols)].astype(float)
                obj_flt_cols = [
                    i
                    for i in EIS_pars_all.columns
                    if str(EIS_pars_all[i].dtype) == "object" and i in float_cols
                ]
                EIS_pars_all[obj_flt_cols] = EIS_pars_all[obj_flt_cols].replace("", 0)
                EIS_pars_all[list(float_cols)] = EIS_pars_all[list(float_cols)].astype(
                    float
                )
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

                def _add_WB_pars(EIS_pars_all):
                    _WB_RAW_daily_path = eis_daily.get("daily_path_RAW_WB")
                    if _WB_RAW_daily_path.exists():
                        _EIS_WB_pars_all = pd.read_pickle(_WB_RAW_daily_path)
                    else:
                        _WB_lst = list(read_df(_EIS_WB_fls))
                        _EIS_WB_pars_all = pd.concat(
                            _WB_lst, sort=False, ignore_index=True
                        )
                        _EIS_WB_pars_all.to_pickle(_WB_RAW_daily_path)

                    _diffcols = set(EIS_pars_all.columns).difference(
                        _EIS_WB_pars_all.columns
                    )
                    _mcols = [
                        i
                        for i in set(EIS_pars_all.columns).intersection(
                            _EIS_WB_pars_all.columns
                        )
                        if i
                        not in [
                            "sourceFilename",
                            "source_mtime",
                            "source_delta_mtime",
                            "sourcebasename",
                        ]
                    ]
                    _dtype_mismatch = [
                        (i, EIS_pars_all[i].dtype, _EIS_WB_pars_all[i].dtype)
                        for i in _mcols
                        if EIS_pars_all[i].dtype != _EIS_WB_pars_all[i].dtype
                    ]
                    if _dtype_mismatch:
                        _excl = []
                        for i in _dtype_mismatch:
                            try:
                                _EIS_WB_pars_all[i[0]] = _EIS_WB_pars_all[i[0]].astype(
                                    i[1]
                                )
                            except Exception as e:
                                _excl.append(i[0])
                                print(i, "\n", e)
                        _mcols = [i for i in _mcols if i not in _excl]
                        # EIS_pars_all[i[0]] = EIS_pars_all[i[0]].astype(i[2])
                    _merge = pd.merge(
                        EIS_pars_all,
                        _EIS_WB_pars_all,
                        on=_mcols,
                        how="left",
                        suffixes=("", "_WB"),
                    )
                    if not _merge.empty:
                        return _merge
                    else:
                        print("WB merge was empty")
                        return EIS_pars_all

                EIS_pars_all = _add_WB_pars(EIS_pars_all)

                EIS_pars_all = EIS_pars_all.assign(
                    **{
                        "EIS_fake": [
                            "fakeZmean" in Path(i).name
                            for i in EIS_pars_all.PAR_file.to_numpy()
                        ]
                    }
                )
                _not_in_index = EIS_pars_all.loc[
                    (
                        ~(EIS_pars_all.PAR_file.isin(EC_index.PAR_file.values))
                        & (~EIS_pars_all.EIS_fake == True)
                    )
                ]
                CleanUpCrew(
                    list_of_files=_not_in_index.sourceFilename.unique(), delete=True
                )

                EIS_pars_all = EIS_pars_all.iloc[
                    ~(EIS_pars_all.index.isin(_not_in_index.index))
                ]
                EIS_pars_all = Load_from_Indexes.test_update_from_index(
                    EIS_pars_all, EC_index
                )
                EIS_pars_all.to_pickle(eis_daily.get("daily_path_RAW"))
            #                EIS_pars_all = pd.read_pickle(eis_daily.get('daily_path_RAW'))
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
                (EIS_pars_all.source_mtime > pd.Timestamp(dt.date(2020, 11, 25)))
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

            #            EIS_pars_undup = Load_from_Indexes.add_missing_ECindex_cols(EC_index, EIS_pars_undup)
            _oc_OVV = list(EIS_pars_undup.columns.intersection(EIS_OVV.columns))
            if not set(EIS_OVV.groupby(_oc_OVV).groups.keys()).intersection(
                EIS_pars_undup.groupby(_oc_OVV).groups.keys()
            ):
                _drpcols = [
                    a
                    for a in EIS_pars_undup.columns
                    if (
                        a in [i for i in _oc_OVV if i not in "PAR_file"]
                        or "_".join(a.split("_")[0:-1])
                        in [i for i in _oc_OVV if i not in "PAR_file"]
                    )
                ]
                #                EIS_pars_undup.drop(columns =_drpcols)
                EIS_pars_undup = Load_from_Indexes.add_missing_ECindex_cols(
                    EC_index, EIS_pars_undup.drop(columns=_drpcols)
                )
            #            EIS_pars_undup = pd.merge(EIS_pars_undup,EIS_OVV,on=_oc_OVV, how='left')

            _oc_SC = list(EIS_pars_undup.columns.intersection(SampleCodes.columns))
            EIS_pars_undup = pd.merge(
                EIS_pars_undup, SampleCodes, how="left", on=_oc_SC
            )

            EIS_pars_BRUTE = EIS_pars_undup.loc[
                (EIS_pars_undup.BRUTE_FIT == 1) | (EIS_pars_undup.FINAL_FIT == 0)
            ]
            if BRUTE_out:
                EIS_pars_BRUTE.to_pickle(eis_daily["daily_path_BRUTE"])

            EIS_pars = EIS_pars_undup.loc[(EIS_pars_undup.FINAL_FIT == 1)]
            EIS_pars = EIS_extra_methods.add_best_model_per_spectrum(EIS_pars)
            EIS_pars.to_pickle(eis_daily["daily_path"])

            _logger.info(f'EIS_pars OVV to daily pickle: {eis_daily.get("daily_path")}')

        _err_type = "lmfit_MSE"
        _filter = "(EIS_pars.lmfit_MSE < 65E4) & (EIS_pars.Rct < 2E3) & (EIS_pars.Rct > 2E-2) \
                                   & (EIS_pars.Rs > 0.01)  & (EIS_pars.Rs < 200) & (EIS_pars.Cdlp < 0.075)\
                                   & (EIS_pars.lmfit_redchi < 1E3)  & (EIS_pars.Aw < 10E3) & (EIS_pars.Aw > 10E-2)\
                                   & (EIS_pars.Qad < 1) & (EIS_pars.tau < 1E3)"
        _filter += '& (EIS_pars.SampleID.str.contains("JOS1|JOS2|JOS3|JOS4|JOS5"))'
        _filter += "& (EIS_pars.EIS_fake == False)"
        _grps = ["Model_EEC", "Gas", "lmfit_var_names"][0:2]
        best_models = (
            EIS_pars.loc[eval(_filter)]
            .dropna(axis=0, subset=[_err_type])
            .groupby(_grps)[_err_type]
            .agg(["count", "mean", "std"])
            .sort_values("mean", ascending=True)
        )
        print(best_models)
        keep_models = (
            best_models.loc[(best_models["count"] > 5) & (best_models["std"] > 0)]
            .index.get_level_values(0)
            .unique()
        )
        EIS_pars = EIS_pars.loc[EIS_pars.Model_EEC.isin(keep_models)]
        best_models = (
            EIS_pars.loc[eval(_filter)]
            .dropna(axis=0, subset=[_err_type])
            .groupby(_grps)[_err_type]
            .agg(["count", "mean", "std"])
            .sort_values(["Gas", "mean"], ascending=True)
        )
        print(best_models)

        if hasattr(EIS_pars, "best_mod_name"):
            # EIS_best_mods = EIS_pars.loc[EIS_pars.Model_EEC_name.isin([i for i in EIS_pars.best_mod_name.unique() if not pd.isna(i)])]
            EIS_best_mods = EIS_pars.loc[
                EIS_pars.index.isin(
                    [i for i in EIS_pars.best_mod_n.unique() if not pd.isna(i)]
                )
            ]

            _agg = (
                EIS_best_mods.dropna(subset=[_err_type])
                .groupby(_grps + ["E_RHE"])[_err_type]
                .agg(["count", "mean", "std"])
            )

            _agg_best = _agg.loc[_agg["count"] > 3].sort_values(
                ["Gas", "E_RHE", "mean"], ascending=True
            )

        #             fast_checking_EEC_models =[(2, 'EEC_2CPEpRW',50),
        #                                        (3, 'EEC_2CPEpW',120),(4,'EEC_2CPE_W',100),
        #                                        (5, 'EEC_2CPE',100), (6,'EEC_Randles_RWpCPE_CPE',60)]
        # #                ['Model(Singh2015_RQRQR)', 'Model(Singh2015_RQRWR)', 'Model(Singh2015_R3RQ)', 'Model(Bandarenka_2011_RQRQR)' ]
        if extra_plotting == "blocked":
            for n, r in best_models.head(1).iterrows():
                modname = r.name[0]
                varnames = [
                    a
                    for i in EIS_pars.loc[
                        EIS_pars["Model_EEC"] == modname
                    ].lmfit_var_names.unique()
                    for a in i.split(", ")
                ]
                #            [1]]+[fast_checking_EEC_models[4]]:
                # modname = f'Model({_modname})'
                EIS_pars_fltr = EIS_pars.loc[
                    (EIS_pars["Model_EEC"] == modname) & eval(_filter)
                ]
                for var in varnames:
                    EIS_pars_fltr.query("pH < 7 & Rct < 2E3").plot(
                        y=var,
                        x="E_RHE",
                        c="BET_cat_agg",
                        colormap="rainbow_r",
                        kind="scatter",
                        title=modname,
                        logy=0,
                    )

                # .query('pH < 15').plot(y='Rs',x='E_RHE',c='pH',colormap='rainbow_r',kind='scatter',ylim=(0,100),title=modname)
                EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH < 15").plot(
                    y="Qad",
                    x="E_RHE",
                    c="pH",
                    colormap="rainbow_r",
                    kind="scatter",
                    ylim=(0, 0.05),
                    title=modname,
                )
                EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH < 7").plot(
                    y="R_ion",
                    x="E_RHE",
                    c="BET_cat_agg",
                    colormap="rainbow_r",
                    kind="scatter",
                    title=modname,
                )
                EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH < 7").plot(
                    y="tau",
                    x="E_RHE",
                    c="BET_cat_agg",
                    colormap="rainbow_r",
                    kind="scatter",
                    ylim=(0, 100),
                    title=modname,
                )
                EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH < 7").plot(
                    y="Rct",
                    x="E_RHE",
                    c="BET_cat_agg",
                    colormap="rainbow_r",
                    kind="scatter",
                    ylim=(0.1, 1e4),
                    logy=True,
                    title=modname,
                )
                if (
                    not EIS_pars.loc[EIS_pars["Model_EEC"] == modname]
                    .query("pH > 7")
                    .empty
                ):
                    EIS_pars.loc[EIS_pars["Model_EEC"] == modname].query("pH > 7").plot(
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
        if "update_index" in kwargs.keys():
            pass

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
            _logger.info(
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
    def HER_pars_OVV(reload=False, use_daily=True, extra_plotting=False, xls_out=False):
        #        exp_type = 'H
        #        PostDestDir = Load_from_Indexes.PostDestDir
        her_daily = get_daily_pickle(exp_type="HER_pars")
        #        IndexOVV_HER_pars_fn = FindExpFolder('VERSASTAT').PostDir.joinpath('Pars_IndexOVV_HER_v{0}.pkl.compress'.format(FileOperations.version))

        if her_daily.get("_exists", False) and reload != True:
            #            Cdl_pars_char = pd.read_excel(IndexOVV_N2_pars_fn,index_col=[0])
            HER_pars_char = pd.read_pickle(her_daily.get("daily_path"))
            HER_pars_char = FileOperations.ChangeRoot_DF(
                HER_pars_char, [], coltype="string"
            )

        else:
            # @@ Check POST_AST status from OVV and PRM
            EC_index, SampleCodes = Load_from_Indexes.get_EC_index()

            def read_df(_par_fls, read_types=["HER_pars"]):
                #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
                while True:
                    try:
                        i = next(_par_fls)
                        _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
                        _delta_mtime = dt.datetime.now() - _source_mtime
                        _i_stem = i.stem
                        _pparts = i.parent.parts

                        if f"HER_v{FileOperations.version}" in _pparts[-2]:
                            if _i_stem.startswith("HER") or "HER" in _i_stem.split("_"):
                                #                            any([_i_stem.startswith(_p) for _p in ['N2_HER|N2_EIS']]):
                                _type = "HER_pars"
                            else:
                                _type = "HER_unknown"
                        else:
                            _type = "_unknown"

                        _meta = {
                            "sourceFilename": i,
                            "source_mtime": _source_mtime,
                            "source_delta_mtime": _delta_mtime,
                            "source_basename": _i_stem,
                            "source_type": _type,
                        }
                        if _type in read_types:
                            _pp = pd.read_excel(i, index_col=[0])
                            _pp = FileOperations.ChangeRoot_DF(
                                _pp, [], coltype="string"
                            )
                            _pp = _pp.assign(**_meta)
                        else:
                            _pp = pd.DataFrame(_meta, index=[0])

                        if not "Analysis_date" in _pp.columns:
                            _pp = _pp.assign(
                                **{
                                    "Analysis_date": dt.datetime.fromtimestamp(
                                        i.stat().st_ctime
                                    )
                                }
                            )

                        _meta.update({"DF": _pp})
                        yield _meta
                    except StopIteration:
                        return "all done"
                        print("gen empty")

            if her_daily.get("_raw_exists", False) and use_daily:
                HER_pars_all = pd.read_pickle(her_daily.get("daily_path_RAW"))
            elif her_daily.get("daily_options_RAW", False) and use_daily:
                HER_pars_all = pd.read_pickle(her_daily.get("daily_options_RAW")[-1])
            else:  # Construct new N2 pars ovv from reading in files
                HER_OVV = EC_index.loc[EC_index.PAR_exp.str.contains("HER")]
                _par_files = [
                    list(
                        Path(d.joinpath(f"HER_v{FileOperations.version}")).rglob(
                            "*xlsx"
                        )
                    )
                    for d in HER_OVV.Dest_dir.unique()
                ]
                _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
                _par_reads = read_df(_par_fls, read_types=["HER_pars"])
                _reads_out = [i for i in _par_reads]
                HER_pars_all = pd.concat(
                    [i["DF"] for i in _reads_out], sort=False, ignore_index=True
                )
                not_in_index = HER_pars_all.loc[
                    ~HER_pars_all.PAR_file.isin(EC_index.PAR_file.values)
                ]
                if not_in_index.empty:
                    print("HER pars, not-in-index is empty... success!")
                else:
                    print("HER pars, not-in-index is NOT empty... delete wrong pars??")
                #                    CleanUpCrew(list_of_files = not_in_index.SourceFilename.unique(), delete = True)
                HER_pars_all = HER_pars_all.loc[
                    HER_pars_all.PAR_file.isin(EC_index.PAR_file.values)
                ]

            HER_pars_recent = HER_pars_all.loc[
                HER_pars_all.Analysis_date > dt.datetime.fromisoformat("2020-07-15")
            ]
            for n, gr in HER_pars_recent.groupby("_type"):
                print(
                    n,
                    f" len {len(gr)}",
                    f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}',
                )
            HER_pars_recent.to_pickle(her_daily["daily_path_RAW"])
            #            ORR_merge_cols = [i for i in ORR_pars.columns if i in ORR_pars_index.columns and not 'Segment' in i]
            #            p2,ovv2 = ORR_pars.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols), ORR_pars_index.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols)
            #            ORR_pars_ovv  = p2.join(ovv2,rsuffix='_ovv').reset_index()
            #            ORR_pars_ovv.query('(pH < 7)').plot(y='E_onset',x='Loading_cm2',kind='scatter',logy=False)
            #        ORR_pars_ovv = pd.merge(ORR_pars,ORR_pars_index,on=ORR_merge_cols,suffixes=('','_ovv'),how='left')
            #        ORR_pars = pd.merge(ORR_pars,postOVVout,on=['PAR_file','SampleID','Electrolyte','pH','postAST'],how='left',suffixes=('','_ovv'))
            #            print('Leftover SampleIDs: {0}'.format(set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())))
            HER_pars_char = pd.merge(
                HER_pars_recent, SampleCodes, on="SampleID", how="left"
            )
            HER_pars_char = pd.merge(
                HER_pars_char, EC_index, on="PAR_file", suffixes=("", "_index")
            )
            ### Fixing the pars after loading...
            # TODO : taking out duplicates based on time_since_run....
            Load_na = HER_pars_char.loc[HER_pars_char.Loading_cm2.isna()]
            if not Load_na.empty:
                Load_na_missingvalues = [
                    (n, *GetSampleID.ink_loading_from_filename(i.PAR_file))
                    for n, i in Load_na.iterrows()
                ]
                Load_na_vals = (
                    pd.DataFrame(Load_na_missingvalues)
                    .rename(columns={1: "Loading_name", 2: "Loading_cm2"})
                    .set_index([0])
                )
                HER_pars_char.Loading_cm2.fillna(
                    value=Load_na_vals.Loading_cm2, inplace=True
                )
            #            ORR_char_merge_cols = [i for i in ORR_pars.columns if i in SampleCodes.columns]
            #            ORR_pars_char = pd.merge(ORR_pars,SampleCodes,on=ORR_char_merge_cols,how='left')
            HER_pars_char = HER_pars_char.drop(
                columns=[i for i in HER_pars_char.columns if "Unnamed" in i]
            )
            if HER_pars_char.loc[HER_pars_char.Loading_cm2.isna()].empty == False:
                HER_pars_char.Loading_cm2 = HER_pars_char.Loading_cm2.fillna(
                    value=0.379
                )  # fillna for Loading_cm2
            HER_pars_char.Loading_cm2 = HER_pars_char.Loading_cm2.round(3)
            HER_pars_char.HER_at_E_slice = HER_pars_char.HER_at_E_slice.round(3)

            if HER_pars_char.postAST.dropna().empty:
                HER_pars_char = HER_pars_char.drop(columns="postAST")
                #                 _int = list(set(ORR_pars_char.columns).intersection(set(EC_index.columns)))
                HER_pars_char = pd.merge(
                    HER_pars_char, EC_index, on="PAR_file", suffixes=("", "_index")
                )

            HER_pars_char = make_uniform_RPM_DAC(HER_pars_char)
            #                 ORR_pars_char = pd.merge(ORR_pars_char, EC_index[['PAR_file', 'postAST']], on = 'PAR_file')
            _sgdct = []
            for pf, pfgrp in HER_pars_char.groupby("PAR_file"):
                _segs = pfgrp["Segment #"].unique()
                for _n, _seg in enumerate(_segs):
                    _sgdct.append({"PAR_file": pf, "Segment #": _seg, "HER_Segnum": _n})
            _HER_segnums = pd.DataFrame(_sgdct)
            HER_pars_char = pd.merge(
                HER_pars_char, _HER_segnums, on=["PAR_file", "Segment #"]
            )

            #            ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna() == True]
            #            if xls_out:
            #                IndexOVV_HER_pars_fn = FileOperations.CompareHashDFexport(HER_pars_char,IndexOVV_HER_pars_fn)
            HER_pars_char.to_pickle(her_daily["daily_path"])

        if extra_plotting:
            jmA2_slice = HER_pars_char.loc[(HER_pars_char["Segment #"] > 1)].query(
                '(HER_type == "j_slice_onset") & (HER_at_J_slice == -2)'
            )
            jmA2_slice.plot(
                x="Metal_wt", y="HER_Tafel_slope", kind="scatter", ylim=(0, 1e3)
            )
            jmA2_slice.plot(
                x="N_content",
                y="HER_Tafel_slope",
                s=50,
                c="g",
                kind="scatter",
                ylim=(0, 1e3),
            )
            #        HER_atE = HER_pars_char.loc[(HER_pars_char['Segment #'] > 1) & np.isclose(HER_pars_char[EvRHE+'_upper'],-0.3,atol=0.02)].query('(E_type == "E_slice")')
            if extra_plotting:
                E_350mV_slice = HER_pars_char.loc[
                    (HER_pars_char["Segment #"] > 1)
                ].query(
                    '(HER_type == "E_slice") & (HER_at_E_slice < -0.29) & (HER_at_E_slice > -0.33)'
                )
                fig, ax = plt.subplots()
                for n, Hgr in E_350mV_slice.groupby(["postAST", "RPM"]):
                    c_set = "g" if "no" in n else "r"
                    _ms_set = "o" if n[-1] < 100 else "*"
                    Hgr.plot(
                        x="N_content",
                        y="HER_J_upper",
                        s=50,
                        c=c_set,
                        kind="scatter",
                        label=n,
                        title="HER at -0.3 Vrhe, j vs N_content",
                        ax=ax,
                        **{"marker": _ms_set},
                    )
                E_350mV_slice.plot(
                    x="N_content",
                    y="HER_J_upper",
                    kind="bar",
                    title="HER, j vs N_content at",
                )
                E_350mV_slice.plot(
                    x="BET_cat_agg",
                    y="HER_J_upper",
                    s=50,
                    c="g",
                    kind="scatter",
                    title="HER, j vs N_content at",
                )
        return HER_pars_char

    def old_HER():
        if IndexOVV_HER_pars_fn.exists() and reload is not True:
            HER_pars_char = pd.read_pickle(IndexOVV_HER_pars_fn)

            if HER_pars_char.SourceFilename.iloc[0].exists() == False:
                HER_pars_char = FileOperations.ChangeRoot_DF(
                    HER_pars_char, [], coltype="string"
                )
        #            ORR_pars_char = ORR_pars_char.drop_duplicates(subset=ORR_pars_char.columns[0:19])
        elif reload == "pickle":
            IndexOVV_HER_pars_fn_pkl = list(
                PostDestDir.rglob(
                    f"{today.year}-{today.month}-*_HER_pars_{system()}.pkl.compress"
                )
            )[-1]
            HER_pars_char = pd.read_pickle(IndexOVV_HER_pars_fn_pkl)

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
            _logger.info(
                "PostEC HPRR re-indexed and saved: {0}".format(
                    new_IndexOVV_HERpars_target
                )
            )

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
            _logger.info(
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
    def ORR_pars_OVV(reload=False, extra_plotting=False, xls_out=False, use_daily=True):
        #        exp_type = 'H
        PostDestDir = Load_from_Indexes.PostDestDir
        orr_daily = get_daily_pickle(exp_type="ORR_pars")
        IndexOVV_ORRpars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "Pars_IndexOVV_ORR_v{0}.pkl.compress".format(FileOperations.version)
        )

        if IndexOVV_ORRpars_fn.exists() and reload is not True:
            #            ORR_pars_char = pd.read_excel(IndexOVV_ORRpars_fn,index_col=[0])
            ORR_pars_char = pd.read_pickle(IndexOVV_ORRpars_fn)

            if ORR_pars_char.sourceFilename.iloc[0].exists() == False:
                ORR_pars_char = FileOperations.ChangeRoot_DF(
                    ORR_pars_char, [], coltype="string"
                )

        #            ORR_pars_char = ORR_pars_char.drop_duplicates(subset=ORR_pars_char.columns[0:19])
        elif reload == "pickle":
            IndexOVV_ORRpars_fn_pkl = list(
                PostDestDir.rglob(
                    f"{today.year}-{today.month}-*_ORR_pars_{system()}.pkl.compress"
                )
            )[-1]
            ORR_pars_char = pd.read_pickle(IndexOVV_ORRpars_fn_pkl)
        else:
            # @@ Check POST_AST status from OVV and PRM
            #            ORR_pars_index = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_Pars')
            #            ORR_pars_index_RRDE = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_RRDE')
            #            ORR_pars_index_RRDE_Chrono = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_RRDE_Chrono').drop_duplicates(subset=['PAR_file','DestFile','Time_since_run']) # cathodic
            #            ORR_Pars_files = [i for i in ORR_pars_index['SourceFilename'].unique() if re.search('(?i)(_pars|_v20)', Path(i).stem) and Path(i).exists()]
            #            ORR_pars_raw = pd.concat([pd.read_excel(i,index_col=[0]) for i in ORR_Pars_files],sort=False)
            #            ORR_pars_raw.PAR_file.fillna(value=ORR_pars_raw.File,inplace=True)
            #            ORR_pars = ORR_pars_raw.drop(columns=['File'],axis=1)
            #            .rename(columns={'File' : 'PAR_file'})
            #            ORR_pars =  FileOperations.ChangeRoot_DF(ORR_pars,[i for i in ORR_pars.columns if re.search('([F-f]ile)',i)],coltype='string')
            #            ORR_pars.PAR_file = ORR_pars.PAR_file.astype(str)
            #            ORR_pars_index.PAR1_file = ORR_pars_index.PAR_file.astype(str)

            EC_index, SampleCodes = Load_from_Indexes.get_EC_index()

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
                            "sourceFilename": i,
                            "source_mtime": _source_mtime,
                            "source_delta_mtime": _delta_mtime,
                            "source_basename": _i_stem,
                            "source_type": _type,
                        }
                        if _type in read_types:
                            _pp = pd.read_excel(i, index_col=[0])
                            _pp = FileOperations.ChangeRoot_DF(
                                _pp, [], coltype="string"
                            )
                            _pp = _pp.assign(**_meta)

                        else:
                            _pp = pd.DataFrame(_meta, index=[0])
                        _meta.update({"DF": _pp})
                        yield _meta
                    except StopIteration:
                        return "all done"
                        print("gen empty")

            if orr_daily.get("_raw_exists", False) and use_daily:
                ORR_pars_all = pd.read_pickle(orr_daily.get("daily_path_RAW"))
            elif orr_daily.get("daily_options_RAW", False) and use_daily:
                ORR_pars_all = pd.read_pickle(orr_daily.get("daily_options_RAW")[-1])
            else:  # Construct new N2 pars ovv from reading in files
                ORR_OVV = EC_index.loc[EC_index.PAR_exp == "ORR"]
                _par_files = [
                    list(
                        Path(d.joinpath(f"ORR_v{FileOperations.version}")).rglob(
                            "*xlsx"
                        )
                    )
                    for d in ORR_OVV.Dest_dir.unique()
                ]
                _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
                _par_reads = read_df(_par_fls, read_types=["ORR_pars", "KL_pars"])
                _reads_out = [i for i in _par_reads]

                ORR_pars_all = pd.concat(
                    [i["DF"] for i in _reads_out], sort=False, ignore_index=True
                )
                ORR_pars_all.to_pickle(orr_daily["daily_path_RAW"])
                _ORR_type_grp = ORR_pars_all.groupby("source_type")

                pprint.pprint(
                    dict(
                        sorted(Counter([i["source_type"] for i in _reads_out]).items())
                    )
                )
                pprint.pprint(
                    {n: gr.sourceFilename.nunique() for n, gr in _ORR_type_grp}
                )
                not_in_index = ORR_pars_all.loc[
                    ~ORR_pars_all.PAR_file.isin(EC_index.PAR_file.unique())
                ]
                if len(not_in_index) > 10:
                    pprint.pprint(
                        {
                            n: gr.sourceFilename.nunique()
                            for n, gr in not_in_index.groupby("source_type")
                        }
                    )
                    _logger.warning(
                        f"PostEC ORR not-in-index is too large: {len(not_in_index)}"
                    )
                CleanUpCrew(
                    list_of_files=not_in_index.sourceFilename.unique(), delete=False
                )
                ORR_pars_all_index = ORR_pars_all.loc[
                    ORR_pars_all.PAR_file.isin(EC_index.PAR_file.values)
                ]
                pprint.pprint(
                    {n: gr.sourceFilename.nunique() for n, gr in _ORR_type_grp}
                )

                ORR_pars_recent = ORR_pars_all_index.loc[
                    ORR_pars_all_index.source_mtime
                    > dt.datetime.fromisoformat("2020-09-09")
                ]
                for n, gr in ORR_pars_recent.groupby("source_type"):
                    print(
                        n,
                        f" len {len(gr)}",
                        f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}',
                    )
                ORR_pars_char = pd.merge(
                    ORR_pars_recent, SampleCodes, on="SampleID", how="left"
                )

            #                _1 = ORR_pars_all.loc[~ORR_pars_all.PAR_file_x.isna()].PAR_file.unique()
            #                _2 = ORR_pars_all.loc[~ORR_pars_all.PAR_file_y.isna()].PAR_file.unique()
            #                ttpfs = set(_1).intersection(_2)

            #                _kl_reads = read_df(_par_fls, read_types = ['KL_pars'])
            #                ORR_KL_all = pd.concat(_kl_reads,sort=False,ignore_index=True)

            #            ORR_merge_cols = [i for i in ORR_pars.columns if i in ORR_pars_index.columns and not 'Segment' in i]
            #            p2,ovv2 = ORR_pars.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols), ORR_pars_index.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols)
            #            ORR_pars_ovv  = p2.join(ovv2,rsuffix='_ovv').reset_index()
            #            ORR_pars_ovv.query('(pH < 7)').plot(y='E_onset',x='Loading_cm2',kind='scatter',logy=False)
            #        ORR_pars_ovv = pd.merge(ORR_pars,ORR_pars_index,on=ORR_merge_cols,suffixes=('','_ovv'),how='left')
            #        ORR_pars = pd.merge(ORR_pars,postOVVout,on=['PAR_file','SampleID','Electrolyte','pH','postAST'],how='left',suffixes=('','_ovv'))
            #            print('Leftover SampleIDs: {0}'.format(set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())))

            ### Fixing the pars after loading...
            # TODO : taking out duplicates based on time_since_run....
            Load_na = ORR_pars_char.loc[
                (ORR_pars_char.Loading_cm2.isna())
                & (ORR_pars_char.PAR_file.isna() == False)
            ]
            if not Load_na.empty:
                Load_na_missingvalues = [
                    (n, *GetSampleID.ink_loading_from_filename(i.PAR_file))
                    for n, i in Load_na.iterrows()
                ]
                Load_na_vals = (
                    pd.DataFrame(Load_na_missingvalues)
                    .rename(columns={1: "Loading_name", 2: "Loading_cm2"})
                    .set_index([0])
                )
                ORR_pars_char.Loading_cm2.fillna(
                    value=Load_na_vals.Loading_cm2, inplace=True
                )
            #            ORR_char_merge_cols = [i for i in ORR_pars.columns if i in SampleCodes.columns]
            #            ORR_pars_char = pd.merge(ORR_pars,SampleCodes,on=ORR_char_merge_cols,how='left')
            ORR_pars_char = ORR_pars_char.drop(
                columns=[i for i in ORR_pars_char.columns if "Unnamed" in i]
            )
            if ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna()].empty == False:
                ORR_pars_char.Loading_cm2 = ORR_pars_char.Loading_cm2.fillna(
                    value=0.379
                )  # fillna for Loading_cm2
            ORR_pars_char.Loading_cm2 = ORR_pars_char.Loading_cm2.round(3)

            if ORR_pars_char.postAST.dropna().empty:
                ORR_pars_char = ORR_pars_char.drop(columns="postAST")
                #                 _int = list(set(ORR_pars_char.columns).intersection(set(EC_index.columns)))
                ORR_pars_char = pd.merge(
                    ORR_pars_char, EC_index, on="PAR_file", suffixes=("", "_index")
                )

            ORR_pars_char = make_uniform_RPM_DAC(ORR_pars_char)
            #                 ORR_pars_char = pd.merge(ORR_pars_char, EC_index[['PAR_file', 'postAST']], on = 'PAR_file')

            #            ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna() == True]
            if xls_out:
                IndexOVV_ORRpars_fn = FileOperations.CompareHashDFexport(
                    ORR_pars_char, IndexOVV_ORRpars_fn
                )
            ORR_pars_char.to_pickle(IndexOVV_ORRpars_fn)
            _logger.info(
                "PostEC ORR re-indexed and saved: {0}".format(IndexOVV_ORRpars_fn)
            )
        if extra_plotting:
            for swp, swgrp in ORR_pars_char.query(
                "(pH < 14) & (RPM_DAC > 900)"
            ).groupby("Sweep_Type"):
                fig, (ax1, ax2) = plt.subplots(figsize=(10, 4), ncols=2)
                #                plt.figure()
                swgrp.plot(
                    y="ORR_Jkin_min_750",
                    x="ORR_E_onset",
                    c="pH",
                    title=f"{swp}",
                    kind="scatter",
                    logy=True,
                    colormap="rainbow_r",
                    ylim=[0.1, 50],
                    xlim=(0.5, 1),
                    ax=ax1,
                )
                ax1.set_xlabel("E onset / mV_RHE")
                swgrp.plot(
                    y="ORR_Frac_H2O2_600",
                    x="ORR_E_onset",
                    c="pH",
                    title=f"{swp}",
                    kind="scatter",
                    logy=True,
                    colormap="rainbow_r",
                    ylim=[0.1, 100],
                    xlim=(0.5, 1),
                    ax=ax2,
                )
                # ax2.set_xlabel('E onset / mV_RHE')
                plt.suptitle("ORR with E_onset")
                plt.show()

                fig, (ax1, ax2) = plt.subplots(figsize=(10, 4), ncols=2)
                swgrp.plot(
                    y="ORR_E_onset",
                    x="N2_BG_lin_slope",
                    c="pH",
                    title=f"{swp}",
                    kind="scatter",
                    logy=True,
                    logx=True,
                    colormap="rainbow_r",
                    xlim=[0.01, 4],
                    ylim=(0.5, 1),
                    ax=ax1,
                )
                swgrp.plot(
                    y="ORR_Jkin_min_750",
                    x="N2_BG_lin_slope",
                    c="pH",
                    title=f"{swp}",
                    kind="scatter",
                    logy=True,
                    logx=True,
                    colormap="rainbow_r",
                    xlim=[0.01, 4],
                    ylim=(0.001, 50),
                    ax=ax2,
                )
                # ax2.set_xlabel('E onset / mV_RHE')
                plt.suptitle("ORR with N2_BG lin slope")
                plt.show()

            plt.close()

        #            ORR_pars_char.query('(pH < 14) & (RPM > 900)').plot(y='Jkin_075',x='E_onset',c='pH',kind='scatter',logy=True,colormap='rainbow_r',xlim=(0.5,1))
        return ORR_pars_char

    #    @staticmethod
    #    def ORR_KL_pars_OVV(reload=False, extra_plotting=False, xls_out = False):
    ##        exp_type = 'H
    #        PostDestDir = Load_from_Indexes.PostDestDir
    #        IndexOVV_ORR_KLpars_fn = FindExpFolder('VERSASTAT').PostDir.joinpath('Pars_IndexOVV_ORR-KL_v{0}.pkl.compress'.format(FileOperations.version))
    #
    #        if IndexOVV_ORRpars_fn.exists() and reload is not True:
    #            ORR_KL_pars = pd.read_excel(IndexOVV_ORR_KLpars_fn,index_col=[0])
    #            ORR_pars_char =  FileOperations.ChangeRoot_DF(ORR_pars_char,[],coltype='string')
    #            ORR_pars_char = ORR_pars_char.drop_duplicates(subset=ORR_pars_char.columns[0:19])
    #        elif reload == 'pickle':
    #            IndexOVV_ORRpars_fn_pkl =  list(PostDestDir.rglob(f'{today.year}-{today.month}-*_ORR_pars_{system()}.pkl.compress'))[-1]
    #            ORR_pars_char = pd.read_pickle(IndexOVV_ORRpars_fn_pkl)
    #
    #        else:
    #            #@@ Check POST_AST status from OVV and PRM
    ##            ORR_index_KL_pars = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_KL_pars')
    ##            ORR_index_KL_data = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_KL_data')
    ##            ORR_pars_index_RRDE_Chrono = postOVVout.groupby('Type_output').get_group('ORR_Jkin_calc_RRDE_Chrono').drop_duplicates(subset=['PAR_file','DestFile','Time_since_run']) # cathodic
    ##            ORR_Pars_files = [i for i in ORR_pars_index['SourceFilename'].unique() if re.search('(?i)(_pars|_v20)', Path(i).stem) and Path(i).exists()]
    ##            ORR_pars_raw = pd.concat([pd.read_excel(i,index_col=[0]) for i in ORR_Pars_files],sort=False)
    #
    #            if orr_daily_pickle_path_RAW.exists():
    #                N2_pars_all = pd.read_pickle(orr_daily_pickle_path_RAW)
    #            elif orr_daily_pickle_path_RAW:
    #                if orr_daily_pickle_path_RAW[-1].exists():
    #                    N2_pars_all = pd.read_pickle(orr_daily_pickle_path_RAW[-1])
    #            else: # Construct new N2 pars ovv from reading in files
    #                N2_OVV = EC_index.loc[OnlyRecentMissingOVV.PAR_exp == 'N2_act']
    #                _par_files = [list(Path(d.joinpath('N2_scans_v30')).rglob(f'*.xlsx')) for d in N2_OVV.Dest_dir.unique()]
    #                _par_fls = (a for i in _par_files for a in i) #if 'EIS' in a.name)
    #                _par_reads = read_df(_par_fls)
    #                N2_pars_all = pd.concat(_par_reads,sort=False)
    #
    #                for n,gr in N2_pars_all.groupby('PAR_file'):
    #                    print(n,f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}', ','.join(gr.N2_type.unique()))
    #                N2_pars_all.to_pickle(n2_daily_pickle_path_RAW)
    #
    #            ORR_pars_raw.PAR_file.fillna(value=ORR_pars_raw.File,inplace=True)
    #            ORR_pars = ORR_pars_raw.drop(columns=['File'],axis=1)
    ##            .rename(columns={'File' : 'PAR_file'})
    #            ORR_pars =  FileOperations.ChangeRoot_DF(ORR_pars,[i for i in ORR_pars.columns if re.search('([F-f]ile)',i)],coltype='string')
    #            ORR_pars.PAR_file = ORR_pars.PAR_file.astype(str)
    #            ORR_pars_index.PAR_file = ORR_pars_index.PAR_file.astype(str)
    #
    #            rrde_fls,emptylst = [],[]
    #            for fl in ORR_pars_index.PAR_file.values:
    #                rrde_df_slice = ORR_pars_index_RRDE_Chrono.loc[(ORR_pars_index_RRDE_Chrono.PAR_file == fl)]
    #                if not rrde_df_slice.empty:
    #                    rrde_df =  rrde_df_slice.loc[(rrde_df_slice.Time_since_run.idxmin())]
    #                    if 'Series' in str(type(rrde_df)):
    #                        spf = rrde_df.DestFile
    #                    else:
    #                        if len(rrde_df) == 1: # pd.read_excel(spectra,index_col=[0])
    #                            spf = rrde_df.DestFile.unique()[0]
    #                        elif len(rrde_df) > 1: # pd.read_excel(spectra,index_col=[0])
    #                            spf = rrde_df.DestFile.unique()[0]
    #                            two = rrde_Df
    #        #                   print('ORR prep Took 1st spectra file: {0}'.format(rrde_df))
    #                        else:
    #                            print('ORR prep Missing spectra file: {0}'.format(rrde_df))
    #                            miss = rrde_df
    #                            spf = None
    #                    rrde_fls.append(spf)
    #                else:
    #                    emptylst.append(fl)
    #                    rrde_fls.append(None)
    #            if len(ORR_pars_index.PAR_file.values) != len(rrde_fls):
    #                print('ORR mismatch length for adding data')
    #            else:
    #                print('ORR pars length matches RRDE Datafiles... adding column')
    #                ORR_pars_index = ORR_pars_index.assign(**{'RRDE_DataFile' : rrde_fls})
    #
    #            ORR_merge_cols = [i for i in ORR_pars.columns if i in ORR_pars_index.columns and not 'Segment' in i]
    #            p2,ovv2 = ORR_pars.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols), ORR_pars_index.dropna(subset=ORR_merge_cols).set_index(ORR_merge_cols)
    #            ORR_pars_ovv  = p2.join(ovv2,rsuffix='_ovv').reset_index()
    ##            ORR_pars_ovv.query('(pH < 7)').plot(y='E_onset',x='Loading_cm2',kind='scatter',logy=False)
    #    #        ORR_pars_ovv = pd.merge(ORR_pars,ORR_pars_index,on=ORR_merge_cols,suffixes=('','_ovv'),how='left')
    #    #        ORR_pars = pd.merge(ORR_pars,postOVVout,on=['PAR_file','SampleID','Electrolyte','pH','postAST'],how='left',suffixes=('','_ovv'))
    #            print('Leftover SampleIDs: {0}'.format(set(ORR_pars.SampleID.unique()) - set(SampleCodes.SampleID.unique())))
    #            ORR_pars = pd.merge(ORR_pars_ovv,SampleCodes,on='SampleID',how='left')
    #            # TODO : taking out duplicates based on time_since_run....
    #            Load_na = ORR_pars.loc[ORR_pars.Loading_cm2.isna()]
    #            Load_na_missingvalues =[(n,*GetSampleID.ink_loading_from_filename(i.PAR_file)) for n,i in Load_na.iterrows()]
    #            Load_na_vals = pd.DataFrame(Load_na_missingvalues).rename(columns={1 : 'Loading_name',2 : 'Loading_cm2'}).set_index([0])
    #            ORR_pars.Loading_cm2.fillna(value=Load_na_vals.Loading_cm2,inplace=True)
    #
    #            ORR_char_merge_cols = [i for i in ORR_pars.columns if i in SampleCodes.columns]
    #            ORR_pars_char = pd.merge(ORR_pars,SampleCodes,on=ORR_char_merge_cols,how='left')
    #            ORR_pars_char = ORR_pars_char.drop(columns=[i for i in ORR_pars_char.columns if 'Unnamed' in i])
    #            if ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna()].empty == False:
    #                ORR_pars_char.Loading_cm2 = ORR_pars_char.Loading_cm2.fillna(value=0.379) # fillna for Loading_cm2
    ##            ORR_pars_char.loc[ORR_pars_char.Loading_cm2.isna() == True]
    #            if xls_out:
    #                IndexOVV_ORRpars_fn = FileOperations.CompareHashDFexport(ORR_pars_char,IndexOVV_ORRpars_fn)
    #            ORR_pars_char.to_pickle(IndexOVV_ORRpars_fn)
    #            _logger.info('PostEC ORR re-indexed and saved: {0}'.format(IndexOVV_ORRpars_fn))
    ##        ORR_pars_char.query('(pH < 7) & (RPM > 900)').plot(y='Jkin_075',x='AD/AG',kind='scatter',logy=False)
    #        return ORR_pars_char

    @staticmethod
    def N2_pars_OVV(reload=False, use_daily=True, extra_plotting=False, xls_out=False):
        #        exp_type = 'H
        #        PostDestDir = Load_from_Indexes.PostDestDir
        #        IndexOVV_N2_pars_fn_xls = FindExpFolder('VERSASTAT').PostDir.joinpath('Pars_IndexOVV_CdlN2_v{0}.xlsx'.format(FileOperations.version))
        IndexOVV_N2_pars_fn = FindExpFolder("VERSASTAT").PostDir.joinpath(
            "N2Cdl_pars_IndexOVV_v{0}.pkl.compress".format(FileOperations.version)
        )
        n2_daily = get_daily_pickle(exp_type="N2_all")

        if n2_daily.get("_exists", False) and reload != True:
            #            Cdl_pars_char = pd.read_excel(IndexOVV_N2_pars_fn,index_col=[0])
            Cdl_pars_char = pd.read_pickle(n2_daily.get("daily_path"))
            Cdl_pars_char = FileOperations.ChangeRoot_DF(
                Cdl_pars_char, [], coltype="string"
            )
        else:
            # @@ Check POST_AST status from OVV and PRM
            _logger.info(
                f'START reloading N2_pars OVV from daily {n2_daily["today"]:%Y-%m-%d}'
            )
            #            EC_index = ECRunOVV(load=1).index
            #            ['EXP_dir','Dest_dir','PAR_file','PAR_file_Ring', 'ORR_act_N2_bg','DestFile']
            #            EC_index = FileOperations.ChangeRoot_DF(OnlyRecentMissingOVV,[])
            #            OnlyRecentMissingOVV.PAR_file = OnlyRecentMissingOVV.PAR_file.astype(str)
            #            OnlyRecentMissingOVV['Loading_cm2'] = OnlyRecentMissingOVV['Loading_cm2'].round(3)
            #            SampleCodes = SampleCodesChar().load
            EC_index, SampleCodes = Load_from_Indexes.get_EC_index()

            def read_df(_par_fls, read_types=["Cdl_data", "Cdl_pars"]):
                #                _ps = Path(d).rglob(f'*_pars_v{FileOperations.version}.xlsx' )
                while True:
                    try:
                        i = next(_par_fls)
                        _source_mtime = dt.datetime.fromtimestamp(i.stat().st_mtime)
                        _delta_mtime = dt.datetime.now() - _source_mtime
                        _i_stem = i.stem
                        _meta = {
                            "sourceFilename": i,
                            "source_mtime": _source_mtime,
                            "source_delta_mtime": _delta_mtime,
                            "source_basename": _i_stem,
                        }

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
                        _meta.update({"N2_type": _N2_type})

                        if _N2_type in read_types:
                            _pp = pd.read_excel(i, index_col=[0])
                            _pp = FileOperations.ChangeRoot_DF(
                                _pp, [], coltype="string"
                            )
                            _pp = _pp.assign(**_meta)
                        else:
                            _pp = pd.DataFrame()
                        _meta.update({"DF": _pp})

                        yield _meta
                    except StopIteration:
                        return "all done"
                        print("gen empty")

            if n2_daily.get("_raw_exists", False) and use_daily is True:
                N2_pars_all = pd.read_pickle(n2_daily.get("daily_path_RAW"))
            elif n2_daily.get("daily_options_RAW", False) and use_daily is True:
                if n2_daily.get("daily_options_RAW")[-1]:
                    N2_pars_all = pd.read_pickle(n2_daily.get("daily_options_RAW")[-1])
            else:  # Construct new N2 pars ovv from reading in files
                N2_OVV = EC_index.loc[EC_index.PAR_exp == "N2_act"]
                _par_files = [
                    list(Path(d.joinpath("N2_scans_v30")).rglob("*.xlsx"))
                    for d in N2_OVV.Dest_dir.unique()
                ]
                _par_fls = (a for i in _par_files for a in i)  # if 'EIS' in a.name)
                _par_reads = read_df(_par_fls, read_types=["Cdl_data", "Cdl_pars"])
                N2_pars_all = pd.concat([i["DF"] for i in _par_reads], sort=False)

                for n, gr in N2_pars_all.groupby("PAR_file"):
                    print(
                        n,
                        f'\nSamples: {", ".join([str(i) for i in gr.SampleID.unique()])}',
                        ",".join(gr.N2_type.unique()),
                    )

                N2_pars_all, _missing_index = Load_from_Indexes.check_missing_ECindex(
                    EC_index, N2_pars_all, clean_up=True
                )
                N2_pars_all.to_pickle(n2_daily["daily_path_RAW"])

            #
            N2_type_grps = N2_pars_all.groupby("N2_type")

            if "CV" in N2_type_grps.groups.keys():
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

            _int = list(set(Cdl_pars_char.columns).intersection(set(EC_index.columns)))
            if Cdl_pars_char.postAST.dropna().empty and len(EC_index.columns) != len(
                _int
            ):
                Cdl_pars_char = Cdl_pars_char.drop(columns="postAST")
                #                 _int = list(set(Cdl_pars_char.columns).intersection(set(EC_index.columns)))
                Cdl_pars_char = pd.merge(
                    Cdl_pars_char, EC_index, on="PAR_file", suffixes=("", "_index")
                )

            Cdl_pars_char = Load_from_Indexes.add_missing_ECindex_cols(
                EC_index, Cdl_pars_char
            )

            if xls_out:
                new_N2_pars_char_target = FileOperations.CompareHashDFexport(
                    Cdl_pars_char, IndexOVV_N2_pars_fn
                )
                _logger.info(
                    "PostEC Cdl N2 CVs re-indexed and saved: {0}".format(
                        new_N2_pars_char_target
                    )
                )
            Cdl_pars_char.to_pickle(IndexOVV_N2_pars_fn)

        try:
            Cdl_pars_char.query('(Sweep_Type_N2 == "cathodic") & (pH < 7)').plot(
                y="Cdl",
                x="E_RHE",
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
                acidCdl = Cdl_pars_char.query(
                    '(Sweep_Type_N2 == "cathodic") & (pH < 7)'
                )
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
        except Exception as e:
            _logger.warning(f"PostEC Cdl N2 CVs extra plotting fail:\n{e}")
            # ==== #
        return Cdl_pars_char

    #        .plot(y='Cdl',x='E_RHE',c='BET_cat_agg',colormap='viridis',kind='scatter',ylim=(0,0.03),title='Cdl in alkaline')

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
    #                _logger.info('PostEC EIS re-indexed and saved: {0}'.format(new_IndexOVV_EISpars_target))
    #            except:
    #                pass
    #            EIS_pars_char.query('pH < 17').groupby('Model_EEC').plot(y='RedChisqr',x='E_RHE',colormap='viridis',kind='scatter',yscale='log')


def make_uniform_RPM_DAC(DF):

    _template = [0, 200, 400, 900, 1500, 2000, 2500]
    _out = []
    if "RPM_DAC" in DF.columns and not "RPM_DAC_uni" in DF.columns:
        for _rpm in DF.RPM_DAC.to_numpy():
            _rpmtest = [_t for _t in _template if np.isclose(_t, _rpm, rtol=0.15)]
            if _rpmtest:
                _out.append(_rpmtest[0])
            else:
                _out.append(_rpm)
        DF = DF.assign(**{"RPM_DAC_uni": _out})
    return DF


class EIS_extra_methods:

    _best_mod_cols = [
        "best_mod_n",
        "best_mod_name",
        "best_mod_badvars",
        "best_mod_rank",
        "best_mod_sortn",
    ]

    @staticmethod
    def add_best_model_per_spectrum(EIS_pars):
        _grpkeys = ["PAR_file", "Segment #"]
        EIS_pars_in = EIS_pars
        EIS_Models = Model_Collection()
        _EIS_pars_mod_cols = [
            i for i in EIS_pars.columns if i in EIS_extra_methods._best_mod_cols
        ]
        EIS_pars = EIS_pars.loc[
            EIS_pars.Model_EEC_name.isin([i.name for i in EIS_Models.lmfit_models])
        ]
        EIS_pars = (
            EIS_pars.drop_duplicates(subset=[*_grpkeys, "Model_EEC"])
            .drop(columns=_EIS_pars_mod_cols)
            .dropna(subset=["basename"])
        )
        EIS_pars_grp = EIS_pars.groupby(_grpkeys)
        _tt = EIS_pars.query('SampleID == "JOS4" & Gas == "O2"')
        _tt = EIS_pars.loc[
            EIS_pars.basename.str.contains("O2_EIS-range_1500rpm_JOS4_285")
            & (EIS_pars.E_RHE == 0.7)
        ]
        _tt = EIS_pars.loc[
            (EIS_pars.basename.str.contains("O2_EIS-range_1500rpm_JOS5_285"))
            & (EIS_pars.E_RHE == 0.3)
        ]
        _ttgrp = _tt.groupby(_grpkeys)
        _mod_results = []
        _failed = []
        for pfseg, PF_pars in EIS_pars_grp:  # EIS_pars_grp:
            # for pfseg,PF_pars in _ttgrp: #EIS_pars_grp: # FIXME
            # (pf,seg),PF_pars
            # pfseg = _failed[-2]
            # PF_pars = EIS_pars_grp.get_group(pfseg)
            try:
                _best_result = EIS_extra_methods.find_best_model_per_spectrum(
                    PF_pars, EIS_Models
                )
            except Exception as e:
                _failed.append(pfseg)
                print("Add best mod error:", e, "\n", pfseg)
            _mod_results.append([*pfseg, *_best_result])
        Best_Models = pd.DataFrame(
            _mod_results, columns=[*_grpkeys, *EIS_extra_methods._best_mod_cols]
        )
        _Pars_BestMods_merged = pd.merge(
            EIS_pars, Best_Models, on=[*_grpkeys], how="left"
        )
        _EIS_pars = EIS_extra_methods.add_best_mod_index_col(_Pars_BestMods_merged)
        return _EIS_pars

    @staticmethod
    def add_best_mod_index_col(EIS_pars):
        _grpkeys = ["PAR_file", "Segment #"]
        _result = []
        for pfseg, PF_pars in EIS_pars.groupby(_grpkeys):
            pfseg, PF_pars
            _bestmod = PF_pars.loc[
                PF_pars.Model_EEC_name.isin(PF_pars.best_mod_name.unique())
            ]
            _result.append([*pfseg, int(*_bestmod.index)])
        Best_Models_index = pd.DataFrame(_result, columns=[*_grpkeys, "best_mod_index"])
        _EIS_pars = pd.merge(EIS_pars, Best_Models_index, on=[*_grpkeys], how="left")
        return _EIS_pars

    @staticmethod
    def find_best_model_per_spectrum(
        PF_pars, EIS_Models, var_lim_max=1e5, var_lim_min=1e-8, var_err_lim=5e3
    ):

        PF_pars = PF_pars.loc[PF_pars.lmfit_message.str.contains("satisfied") == True]

        _lmfitcols = [i for i in PF_pars.columns if i.startswith("lmfit")]
        _best_result = []
        if PF_pars.Model_EEC.nunique() >= 2:
            # and PF_pars.Model_EEC.nunique() == len(PF_pars):
            _res = []
            aic_mean1 = PF_pars.lmfit_aic.mean()
            aic_mean = (
                PF_pars.loc[PF_pars.lmfit_aic < aic_mean1].lmfit_aic.mean()
                - 0.01 * aic_mean1
            )
            aic_std = PF_pars.lmfit_aic.std()

            for n, r in PF_pars.iterrows():
                _rowres = []
                _vars = r.lmfit_var_names.split(", ")
                _varserr = [i for i in [i + "_stderr" for i in _vars] if i in r.index]
                _vsum, _vmax = r[_vars].sum(), r[_vars].max()
                _bad_vars = set(
                    [i for i in _vars if r[i] > var_lim_max and r[i] < var_lim_min]
                )
                _verrsum, _verrmax = 0, 0
                _aic_mean_smaller = r.lmfit_aic <= aic_mean

                if _varserr:
                    _verrsum, _verrmax, _verrmean = (
                        r[_varserr].sum(),
                        r[_varserr].max(),
                        r[_varserr].mean(),
                    )
                    _bad_varserr = {
                        i: {
                            "val": r[i],
                            "rel_val": r[i] / r[(i.split("_stderr")[0])],
                            "name": (i.split("_stderr")[0]),
                            "name_val": r[(i.split("_stderr")[0])],
                        }
                        for i in _varserr
                    }
                    _bad_varserr_val = [
                        val["name"]
                        for k, val in _bad_varserr.items()
                        if val["val"] > var_err_lim
                    ]
                    _bad_varserr_perc = [
                        val["name"]
                        for k, val in _bad_varserr.items()
                        if val["rel_val"] > var_err_lim
                    ]
                    _bad_vars_cont = EIS_extra_methods.get_context(r, _vars)
                    _bad_vars_lims = EIS_extra_methods.get_badvars_from_lim(
                        r, _vars, EIS_Models
                    )
                    _bad_vars_lims_bool = bool(_bad_vars_lims)
                    _bad_vars_err = set(
                        _bad_varserr_val
                        + _bad_varserr_perc
                        + _bad_vars_cont
                        + _bad_vars_lims
                    )
                    _bad_vars = _bad_vars.union(_bad_vars_err)
                    _testlow = r.test_low

                    _rowres = [
                        n,
                        r.Model_EEC_name,
                        len(_vars),
                        _vsum,
                        _vmax,
                        ", ".join(_bad_vars),
                        len(_bad_vars),
                        _verrsum,
                        _verrmax,
                        _verrmean,
                        r.lmfit_aic,
                        r.lmfit_redchi,
                        r.lmfit_chiqsr,
                        _aic_mean_smaller,
                        ", ".join(_bad_vars_lims),
                        _bad_vars_lims_bool,
                        _testlow,
                    ]
                _res.append(_rowres)
            var_res_raw = pd.DataFrame(
                _res,
                columns=[
                    "pf_index",
                    "Model_EEC_name",
                    "len_vars",
                    "varsum",
                    "varmax",
                    "bad_vars",
                    "bad_vars_len",
                    "err_varsum",
                    "err_varmax",
                    "err_varmean",
                    "lmfit_aic",
                    "lmfit_redchi",
                    "lmfit_chiqsr",
                    "_aic_mean_smaller",
                    "_bad_vars_lims",
                    "_bad_vars_lims_bool",
                    "_testlow",
                ],
            )
            # var_res_err = var_res.loc[var_res.err_varmax > 0]

            var_res_fltr = var_res_raw.loc[
                ~var_res_raw._bad_vars_lims.str.contains("Rs|Qad")
            ]
            if var_res_fltr.empty:
                var_res_fltr = var_res_raw.loc[
                    ~var_res_raw._bad_vars_lims.str.contains("Rs")
                ]
            if var_res_fltr.empty:
                var_res_fltr = var_res_raw

            var_res = var_res_fltr.loc[
                (var_res_fltr._aic_mean_smaller == True)
                & (var_res_fltr._bad_vars_lims_bool == False)
            ]
            if var_res.empty:
                var_res = var_res_fltr.loc[(var_res_fltr._bad_vars_lims_bool == False)]

            if var_res.empty:
                var_res = var_res_fltr.loc[(var_res_fltr._aic_mean_smaller == True)]

            if var_res.empty:
                var_res = var_res_fltr

            # if
            _rankby = ["_testlow", "bad_vars_len", "lmfit_aic", "len_vars"]
            var_res = var_res.sort_values(by=_rankby)
            best_mod_row = var_res.head(1)
            _sorted_rank = ", ".join([str(i) for i in var_res.pf_index.values])
            _best_result = [
                best_mod_row.pf_index.iloc[0],
                best_mod_row.Model_EEC_name.iloc[0],
                best_mod_row.bad_vars.iloc[0],
                _sorted_rank,
                ", ".join(_rankby),
            ]
            # var_res.bad_vars_len.unique()
        return _best_result

    @staticmethod
    def get_context(r, _vars):
        _bad_vars_cont = []
        if "Gas" in r.index:

            if pd.isna(r.Gas):
                _testgas = Path(r.PAR_file).stem
            else:
                _testgas = r.Gas

            if "Rorr" in _vars:
                if "N2" in _testgas:
                    _bad_vars_cont.append("Rorr")
                elif "O2" in _testgas:
                    if r.Rorr > 5e3:
                        _bad_vars_cont.append("Rorr")

        return _bad_vars_cont

    def get_badvars_from_lim(r, _vars, EIS_Models):
        _modr = [i for i in EIS_Models.lmfit_models if i.name == r.Model_EEC_name]
        _badvarlims = []
        if _modr:
            _modinit = _modr[0]
            for v in _vars:
                vval = r[v]
                if vval < _modinit.parameters_guesses[v].min * 1.5:
                    if not _modinit.name == "RL-TLM(Rct-Qad-W" and v == "Cdlp":
                        pass
                    else:
                        _badvarlims.append(v)

                if vval > _modinit.parameters_guesses[v].max * 0.66:
                    _badvarlims.append(v)

                if v == "Aw" and vval > 3000:
                    _badvarlims.append(v)

                if v == "R_ion" and vval > 2000:
                    _badvarlims.append(v)

        return _badvarlims
        # r.Model_EEC_name


def get_EIS_pars(kwargs):
    EIS_pars = Load_from_Indexes.EIS_pars_OVV(**kwargs)  # EIS_Pars2 6745, 17994
    return EIS_pars


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


#%% === IF MAIN block====
if __name__ == "__main__":
    reloadOVV = False
    # EIS_pars = Load_from_Indexes.EIS_pars_OVV(reload=reloadOVV, extra_plotting=True)
    if reloadOVV:

        PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        today = dt.today()
        postOVVout = Load_from_Indexes.PreparePostOVV(fastload=False)  # len(22965)
        SampleCodes = FindExpFolder().LoadSampleCode()

    #        postOVVout.PAR_file = postOVVout.PAR_file.astype(str)
    # === Loading preparation overview of Samples and merging with the data from Characterization techniques === #
    reloading, Reload_set, _use_daily = True, True, False
    if 0:
        #        logger = start_logger()

        EIS_pars = Load_from_Indexes.EIS_pars_OVV(
            reload=Reload_set,
            extra_plotting=True,
            xls_out=False,
            BRUTE_out=False,
            use_daily=_use_daily,
        )
        #        (reload= Reload_set, source = 'ExpDirs') # EIS_Pars2 6745, 22813
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
            reload=True, use_daily=_use_daily
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
