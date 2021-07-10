# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:42:36 2020

@author: User
"""

import sys
from pathlib import Path

# from collections import namedtuple
# from datetime import datetime
import numpy as np

# import itertools
from typing import Dict, List
from cmath import phase
from dataclasses import dataclass, field
import os
import multiprocessing
import logging

# from functools import partial
from itertools import repeat
import pandas as pd


if __name__ == "__main__":
    from fitting import (
        fit_EEC,
        prepare_rand_models,
        Fit_Spectrum,
        Fit_Spectra_Collection,
    )
    import models

    FH_path = Path(__file__).parent.parent.parent.parent.joinpath("FileHelper")
    sys.path.append(str(FH_path))
    EC_path = Path(__file__).parent.parent.parent
    sys.path.append(str(FH_path))
    from ECpy.experiments.Loader.CreateCV import create_CVs as create_CVs
    from ECpy.main_run_PAR_DW import ECRunOVV
    from FileHelper.FileFunctions import FileOperations
    from FileHelper.FindFolders import FindExpFolder

    sys.path.append(Path(__file__).parent.parent.parent)
    #    from EC_logging_config import start_logging
    from ECpy.runEC.EC_logging_config import start_logging

    from ECpy.experiments.EIS.plotting import (
        EIS_plotting_per_EV,
        EIS_Trimming_plot,
        EIS_plotting_EvRHE,
    )

#    logger = start_logging(__name__)
#    pass
#    FH_path = Path(__file__).parent.parent.parent
#    sys.path.append(str(FH_path))
#    from FileHelper import FileOperations

else:
    from .fitting import prepare_rand_models, Fit_Spectrum, Fit_Spectra_Collection
    from . import models
    from .plotting import EIS_plotting_per_EV, EIS_Trimming_plot, EIS_plotting_EvRHE
    from ..Loader.CreateCV import create_CVs as create_CVs

    FH_path = Path(__file__).parent.parent.parent.parent.joinpath("FileHelper")
    sys.path.append(str(FH_path))
    from FileHelper.FileFunctions import FileOperations
    from FileHelper.FindFolders import FindExpFolder

    sys.path.append(Path(__file__).parent.parent.joinpath("runEC"))
    from ECpy.runEC.EC_logging_config import start_logging

    #    from ECpy.experiments.EIS.plotting import EIS_plotting_per_EV, EIS_Trimming_plot, EIS_plotting_EvRHE
    from ECpy.main_run_PAR_DW import ECRunOVV
#    logger = start_logging(__name__)


_logger = logging.getLogger(__name__)
_logger.propagate = False

globals()["EvRHE"] = "E_AppV_RHE"

# print('TODO: eis run fix MULTIPROCESSING!')
# TODO: do it
# Meta = namedtuple('Meta', 'PAR_file Segment E_dc_RHE E_dc_RHE_mV RPM_DAC data ovv')
@dataclass(order=True, frozen=False)
class EIS_Spectrum:
    """EIS Spectrum dataclass.\n
    Holds spectrum raw data in pd.DataFrame and metadata"""

    spectrum_grp_cols = ["PAR_file", "Segment #", EvRHE, "RPM_DAC"]

    PAR_file: Path = field(default=Path)
    Segment: int = 0
    E_dc_RHE: float = 0.0
    RPM_DAC: float = 0.0

    data: type(pd.DataFrame) = field(default=pd.DataFrame(), repr=False)
    ovv: type(pd.DataFrame) = field(default=pd.DataFrame(), repr=False)
    EIS_kwargs: Dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        #        self.E_dc_RHE = np.round(self.E_dc_RHE, 3)
        self.E_dc_RHE_mV = np.round(self.E_dc_RHE * 1e3, 3)
        self.EvRHE = self.E_dc_RHE
        self.add_complex_impedance_columns()
        self.check_freqlim()

    def __repr__(self):
        try:
            _file = f'File: "{self.PAR_file.relative_to(FindExpFolder("VERSASTAT").DataDir)}"'
        except Exception as e:
            _file = f'File:"{self.PAR_file.parent.name}/{self.PAR_file.name}"'
        _data = f"at {self.E_dc_RHE_mV} mV with {self.RPM_DAC} rpm, data({len(self.data)}) and ovv({len(self.ovv)})"
        _keys = f"attrs: {self.E_dc_RHE_mV} mV with {self.RPM_DAC} rpm, data({len(self.data)}) and ovv({len(self.ovv)})"
        return _file + "\n" + _data

    def add_complex_impedance_columns(
        self, required_cols=["Frequency(Hz)", "Z Real", "Z Imag"]
    ):
        #        ['Frequency(Hz)'
        _check_cols = []
        if hasattr(self.data, "columns"):
            _check_cols = [i in self.data.columns for i in required_cols]

        EIS_data_raw = pd.DataFrame()
        if all(_check_cols):
            _add_cols = {}
            freq = self.data["Frequency(Hz)"].values
            _add_cols.update(
                {
                    "Frequency(Hz)": freq,
                    "Angular": freq * 2 * np.pi,
                    "Ang_Warburg": 1 / (np.sqrt(freq * 2 * np.pi)),
                }
            )
            Zre, Zim = self.data["Z Real"].values, self.data["Z Imag"].values
            Zdata = Zre + 1j * Zim
            Ydata = Zdata ** -1
            Yre, Yim = Ydata.real, Ydata.imag
            DataWeights_modulus_Z = Zre ** 2 + Zim ** 2
            DataWeights_modulus_Y = (Zre ** 2 + Zim ** 2) ** -1

            _add_cols.update(
                {
                    "DATA_Z": Zdata,
                    "DATA_Zre": Zre,
                    "DATA_Zim": Zim,
                    "DATA_-Zim": -1 * Zim,
                    "DATA_Zphase": [phase(i) for i in Zdata],
                    "DATA_Zmod": abs(Zdata),
                    "DATA_Zangle": np.angle(Zdata, deg=True),
                    "DATA_-Zangle": -1 * np.angle(Zdata, deg=True),
                    "DATA_Y": Ydata,
                    "DATA_Yre": Yre,
                    "DATA_Yim": Yim,
                    "Valid": True,
                    "lmfit_weights_mod_Z": DataWeights_modulus_Z,
                    "lmfit_weights_mod_Y": DataWeights_modulus_Y,
                    "lmfit_weights_unit": Zre / Zre,
                }
            )

            _meta_info_cols = {
                key: val
                for key, val in vars(self).items()
                if any([i in type(val).__name__ for i in ["int", "float", "Path"]])
            }
            _add_cols.update(_meta_info_cols)
            EIS_data_raw = pd.DataFrame(_add_cols, index=self.data.index)
        elif self.data.empty:
            pass
        else:
            print(
                "Error in EIS_spectrum, missing columns:",
                ", ".join([i for i in self.data.columns if i not in required_cols]),
            )
        self.EIS_data = EIS_data_raw

    def check_freqlim(self):
        if not self.EIS_data.empty:
            FreqLim = self.EIS_kwargs.get("FreqLim", 30e3)
            self.EIS_data.loc[
                (self.EIS_data["Frequency(Hz)"] >= FreqLim), "Valid"
            ] = False  # TODO
            #            & (self.EIS_data['DATA_-Zim'] < -15
            self.EIS_data_freqlim = self.EIS_data.query("Valid == True")

    def EIS_exp_fit_cols(self):
        _check = ["PAR_file", "Segment", EvRHE, "RPM_DAC", "Model_EEC"]
        return list(self.__dataclass_fields__.keys())[:-3] + ["Model_EEC"]


@dataclass(order=True)
class EIS_spectra_collection:
    global EvRHE

    PAR_file: Path = field(default=Path)
    spectra: List[EIS_Spectrum] = field(default=list)
    data: type(pd.DataFrame) = field(default=pd.DataFrame(), repr=False)
    ovv: type(pd.DataFrame) = field(default=pd.DataFrame(), repr=False)
    EIS_kwargs: Dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.check_if_list_contains_EIS_spectra()
        self.make_mean_spectrum()

    def check_if_list_contains_EIS_spectra(self):
        _check = [isinstance(i, EIS_Spectrum) for i in self.spectra]
        if not all(_check) and self.spectra:
            raise ValueError

    def make_mean_spectrum(self):

        _new_PF_mean = self.PAR_file.with_name(
            self.PAR_file.stem + "_fakeZmean" + self.PAR_file.suffix
        )
        _PF_mean_ovv = self.ovv.loc[self.ovv.PAR_file == self.PAR_file]
        #    _PF_mean_ovv['Measured_OCP'] =  [i[0] for i in _PF_mean_ovv['_act0_Measured Open Circuit'].str.split()]
        #    _PF_mean_ovv['PAR_file'] = _new_PF_mean
        _PF_mean_ovv = _PF_mean_ovv.assign(**{"PAR_file": _new_PF_mean})
        _PF_data = pd.concat(i.data for i in self.spectra)
        _numcols = [
            i
            for i in _PF_data.columns
            if _PF_data[i].dtype != "O" and not "Frequency" in i
        ]

        _PF_data_mean = _PF_data.groupby("Frequency(Hz)")[_numcols].mean().reset_index()
        #    _PF_data_mean.plot(x='Z Real',y='Z Imag', kind='scatter',c='Segment #') # test plot
        _PF_data_mean = _PF_data_mean.assign(**_PF_mean_ovv.iloc[0].to_dict())
        _PF_data_mean_meta = (
            _PF_data[[i for i in _PF_data.columns if _PF_data[i].nunique() == 1]]
            .iloc[0]
            .to_dict()
        )
        _PF_mean_ovv_meta = _PF_mean_ovv.iloc[0].to_dict()
        _PF_meta_mean = {**_PF_data_mean_meta, **_PF_mean_ovv_meta}
        #        _join_cols = list(_PF_data_mean.columns.intersection(_PF_mean_ovv.columns))
        #    _PF_data_mean = _PF_data_mean.join(_PF_mean_ovv.set_index('PAR_file'),on=_join_cols,how='left')
        #    _merge = pd.concat([_PF_data_mean, _PF_mean_ovv],axis=1)
        #                       .set_index('PAR_file'),on=_join_cols,how='outer')
        #        _PF_data_mean[[ 'Segment #', EvRHE, 'RPM_DAC']]
        _PF_data_mean_grp = _PF_data_mean.groupby(EIS_Spectrum.spectrum_grp_cols)
        self.mean_spectrum = [
            EIS_Spectrum(
                Path(PF),
                int(Seg),
                np.round(float(E_V), 3),
                int(RPM_DAC),
                gr,
                _PF_mean_ovv,
            )
            for (PF, Seg, E_V, RPM_DAC), gr in _PF_data_mean_grp
        ][0]
        self.meta_data = _PF_meta_mean


def _apply_df(args):
    df, func = args
    return df.groupby(level=0).apply(func)


def mp_apply(df, func):
    workers = 4
    pool = mp.Pool(processes=workers)
    split_dfs = np.array_split(df, workers, axis=1)
    result = pool.map(_apply_df, [(d, func) for d in split_dfs])
    pool.close()
    # result = sorted(result, key=lambda x: x[0])
    return pd.concat(result, axis=1)


# if __name__ == '__mainee__':
#    df = pd.DataFrame([[1, 2, 3, 1], [1, 2, 3, 1], [4, 5, 6, 2], [7, 8, 9, 2]], columns=['A', 'B', 'C', 'cluster_id'])
#    df = df.set_index('cluster_id')
#    out = mp_apply(df, my_func)
#    print(out)


def test_func(arg1, arg2, **kwargs):
    print(arg1)
    print(arg2)
    print(kwargs)
    return arg2


# if __name__ == '__main__':
#    list_of_args2 = [1, 2, 3]
#    just_a_dict = {'key1': 'Some value'}
#    with multiprocessing.Pool(processes=3) as pool:
#        results = pool.map(partial(test_func, 'This is arg1', **just_a_dict), list_of_args2)
#    print(results)


def get_eis_suggestions(how="grouped"):
    #    PostDestDir =
    EvRHE = "E_AppV_RHE"
    EIS_version = FileOperations.EIS_version
    PostDestDir = FindExpFolder("VERSASTAT").PostDir
    PDD_eischeck = PostDestDir.joinpath("EIS/redchi_check")
    _out_check = {}
    if "group" in how:
        _out_check = {"EIS_recheck_grpcols": EIS_Spectrum().EIS_exp_fit_cols()}
    for recheck in [
        "EIS_recheck_good_fits",
        "EIS_recheck_bad_fits",
        "EIS_recheck_bad_fits_suggestions",
    ]:
        suggestions_path = PDD_eischeck.joinpath(f"{recheck}.pkl")
        if suggestions_path.is_file():
            try:
                EIS_suggestions = pd.read_pickle(suggestions_path)
                EIS_suggestions[EvRHE] = EIS_suggestions[EvRHE].round(3)

                if not EIS_suggestions.empty:
                    _df = EIS_suggestions
                else:
                    _df = pd.DataFrame()
            except Exception as e:
                _df = pd.DataFrame()
        else:
            _df = pd.DataFrame()
        if "group" in how and not _df.empty:
            _out_check.update({recheck: _df.groupby(EIS_Spectrum().EIS_exp_fit_cols())})
        else:
            _out_check.update({recheck: _df})
    return _out_check


def starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)


def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)


##
# args_iter = zip(repeat(project_name), api_extensions)
# run_groups = [run_group(n,n[1],gr,EC_index) for n,gr in grp_date_dir]
# branches = starmap_with_kwargs(pool, eis_fit_PAR_file, args_iter, kwargs_iter)


def get_recheck_parfiles_only(ovv, **EIS_kwargs):
    pfs = ovv.PAR_file.unique()
    for recheck in [
        "EIS_recheck_good_fits",
        "EIS_recheck_bad_fits",
        "EIS_recheck_bad_fits_suggestions",
    ]:
        kwgrps = EIS_kwargs.get(recheck)

        if kwgrps.groups:
            _pfgrps = [i for i in kwgrps.groups if i[0] in pfs]

            #        _df = EIS_kwargs.get(recheck,pd.DataFrame())
            if _pfgrps:
                _dfgrp = pd.concat([kwgrps.get_group(i) for i in _pfgrps]).groupby(
                    EIS_Spectrum().EIS_exp_fit_cols()
                )
                if "good_fits" in recheck:
                    _good_grp = pd.concat(
                        [kwgrps.get_group(i) for i in _pfgrps]
                    ).groupby(
                        [i for i in EIS_Spectrum().EIS_exp_fit_cols() if "Model" in i]
                    )
                    _good_pars = {}
                    for n, gr in _good_grp:
                        try:
                            _good_pars.update(
                                {
                                    n: [
                                        list(zip(r.dtype.names, r))
                                        for r in gr.loc[
                                            :,
                                            models.EEC_models_index(select=n)[0][
                                                1
                                            ].param_names,
                                        ].to_records(index=False)
                                    ]
                                }
                            )
                        except:
                            pass
                    EIS_kwargs.update({"good_fit_pars": _good_pars})

                EIS_kwargs.update({recheck: _dfgrp})
    return EIS_kwargs


# TODO continue means:
def means_EIS(fit_run_args):
    EIS_data_pf = pd.concat([i.data for i in fit_run_args])
    EIS_data_pf = EIS_data_pf.assign(
        **{
            "Y Real": ((EIS_data_pf["Z Real"] + 1j * EIS_data_pf["Z Imag"]) ** -1)
            .to_numpy()
            .real,
            "Y Imag": ((EIS_data_pf["Z Real"] + 1j * EIS_data_pf["Z Imag"]) ** -1)
            .to_numpy()
            .imag,
            "Z-Imag": -1 * EIS_data_pf["Z Imag"],
        }
    )
    EIS_data_pf.groupby("Frequency(Hz)").mean().plot(
        x="Z Real", y="Z-Imag", kind="scatter", c=EvRHE, cmap="viridis"
    )
    EIS_data_pf.plot(x="Y Real", y="Y Imag", kind="scatter", c=EvRHE, cmap="viridis")
    for n, gr in EIS_data_pf.groupby("Frequency(Hz)"):
        #    for n,gr in  EIS_data_pf.groupby(EvRHE):
        #        plt.clf()

        fig, ax = plt.subplots()
        gr.plot(
            x="Y Real",
            y="Y Imag",
            kind="scatter",
            c=EvRHE,
            cmap="viridis",
            ax=ax,
            title=f"Freq: {n}",
        )
        plt.show()
        plt.close()


#


def read_in_EIS_data(EIS_ovv):
    EvRHE = "E_AppV_RHE"
    #    EIS_ovv = ovv_exp_grp.get_group('EIS')
    #    Parsout_lst, nm2 = [], Path(nm2)
    #    index_per_file_EV_lst = []
    #    OutDataCombi, AllData_E_file = [], pd.DataFrame(EIS_ovv)
    EISgrEvRHE_data_raw, EISgrEvRHE_actions = create_CVs(EIS_ovv, multi_run=True)
    EISgrEvRHE_data = EISgrEvRHE_data_raw.loc[
        (EISgrEvRHE_data_raw.ActionId.astype(int) == 21)
        | (EISgrEvRHE_data_raw["Z Real"] != 0)
    ]
    EISgrEvRHE_data = EISgrEvRHE_data.round({EvRHE: 3})
    EISgrEvRHE_data["PAR_file"] = [
        Path(i) for i in EISgrEvRHE_data["PAR_file"].to_numpy()
    ]
    EIS_ovv["PAR_file"] = [Path(i) for i in EIS_ovv["PAR_file"].to_numpy()]

    read_data_diff = len(EISgrEvRHE_data_raw) - len(EISgrEvRHE_data)
    if read_data_diff > 0:
        logger.info(
            f"EIS filtered data {read_data_diff} because ActionId != 21 but in {EISgrEvRHE_data_raw.ActionId.unique()}"
        )
    #    for ddir in EIS_ovv.Dest_dir.unique():
    #        EIS_dest_dir = Path(ddir).joinpath('EIS')
    #        EIS_dest_dir.mkdir(parents=True, exist_ok=True)
    dest_files = []
    for n, r in EIS_ovv.iterrows():
        EIS_dest_dir = Path(r.Dest_dir).joinpath("EIS")
        EIS_dest_dir.mkdir(parents=True, exist_ok=True)
        dest_files.append(
            {
                "index": n,
                "PAR_file": Path(r.PAR_file),
                "EIS_dest_dir": EIS_dest_dir,
                "EIS_dest_Pars": EIS_dest_dir.joinpath(
                    Path(r.PAR_file).stem + "_pars.xlsx"
                ),
                "EIS_dest_spectra": EIS_dest_dir.joinpath(
                    Path(r.PAR_file).stem + "_Combined.xlsx"
                ),
            }
        )
    EIS_ovv_destfiles = pd.DataFrame(dest_files).set_index("index")

    EIS_ovv_dest_out = pd.merge(
        EIS_ovv, EIS_ovv_destfiles, left_index=True, on="PAR_file", how="inner"
    )
    #    pd.concat([EIS_ovv, EIS_ovv_destfiles],axis=1)
    #    EISgrEvRHE_data.PAR_file = EISgrEvRHE_data.PAR_file.astype(str)
    EISgrEvRHE_data = EISgrEvRHE_data.join(
        EIS_ovv_destfiles.set_index("PAR_file"), on="PAR_file", how="left"
    )
    return EISgrEvRHE_data, EISgrEvRHE_data_raw, EISgrEvRHE_actions, EIS_ovv_dest_out


#    EIS_fit_kwargs.update({'EIS_dest_dir' : EIS_dest_dir})
#    PF_fit_pars_target = EIS_dest_dir.joinpath(nm2.stem + '_pars.xlsx')
#    Parsout_path_target = FileOperations.CompareHashDFexport(pd.DataFrame(), Parsout_path)
#    PF_fit_spectra_target = EIS_dest_dir.joinpath(nm2.stem + '_Combined.xlsx')
#    index_out = []
# def eis_multi_wrap(fn, args,kwargs):
#    return fn(args[0],args[1],args[-1])


def PF_fit_starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(PF_fit_apply_args_and_kwargs, args_for_starmap)


def PF_fit_apply_args_and_kwargs(fn, args, kwargs):
    return fn(args, **kwargs)


def filter_error_debug():
    "12.02.2019_0.1MKOH_cell2/O2_EIS-range_1500rpm_JOS6_postAST_899 at 0.999 (1500 rpm)"
    #    FileOperations.CompareHashDFexport(EISgrEvRHE_actions.T,EIS_ovv_destfiles.EIS_dest_dir.iloc[0].joinpath('PAR_actions_T.xlsx'))
    bad_fit = ("NSS_0103-EIS-7", 3)
    tt = "N2_EIS-range_1500rpm_JOS4_288"
    fit_run_arg = [i for i in fit_run_args if tt in i[0].name and 558 == i.E_dc_RHE_mV]
    fit_run_args_tests = [i for i in fit_run_args if tt in i[0].name]


class EIS_Preparator:
    global EvRHE

    def __init__(self, *args, **kwargs):
        #        print(type(args[0]))
        _t = [
            n
            for n, i in enumerate([a for a in args if hasattr(type(a), "__name__")])
            if "ECRunOVV" in type(i).__name__
        ]
        #        print(_t)
        if _t:
            instance_attrs = vars(args[_t[0]])
            self.run_kwargs = instance_attrs.get("run_kwargs")
            self.ovv_exp_grp = instance_attrs.get("ovv_exp_grp")
        else:
            self.run_kwargs = kwargs
            self.ovv_exp_grp = kwargs.get(
                "ovv_exp_grp",
                pd.DataFrame({"PAR_exp": "EIS"}, index=[0]).groupby("PAR_exp"),
            )
        #        super().__init__(**kwargs)
        #        self.ovv_exp_grp = ovv_exp_grp
        self.EIS_ovv = pd.concat([gr for n, gr in self.ovv_exp_grp if "EIS" in n])
        self.EIS_kwargs = self.run_kwargs
        #        if self.EIS_kwargs.get()
        self.update_EIS_kwargs()
        self.prepare_EIS_data()

    #        self.prepare_fit_run_args()
    #        self.iter = self.prepare_fit_run_args_iter()
    #    def __getattr__(self, attr):
    #        return getattr(self._ECRunOVV, attr)

    def update_EIS_kwargs(self):
        print("starting prep_random")
        _prepare_rand_params = prepare_rand_models(
            rand_params_dest_dir=FindExpFolder("VERSASTAT").PostDir.joinpath("EIS")
        )
        self.EIS_kwargs.update(**dict(random_params=_prepare_rand_params))
        print("finished prep_random")
        print("start get_eis_suggestions")
        self.EIS_kwargs.update(get_eis_suggestions())
        print("start get_recheck_parfiles_only")
        if not True:  # FIXME
            self.EIS_kwargs = get_recheck_parfiles_only(self.EIS_ovv, **self.EIS_kwargs)

    def prepare_EIS_data(self):
        (
            EISgrEvRHE_data,
            EISgrEvRHE_data_raw,
            EISgrEvRHE_actions,
            EIS_ovv_destfiles,
        ) = read_in_EIS_data(self.EIS_ovv)
        self.EISgrEvRHE_data, self.EISgrEvRHE_data_raw = (
            EISgrEvRHE_data,
            EISgrEvRHE_data_raw,
        )
        self.EISgrEvRHE_actions, self.EIS_ovv_destfiles = (
            EISgrEvRHE_actions,
            EIS_ovv_destfiles,
        )
        self.ovv_exp_grp_PF = self.EIS_ovv_destfiles.groupby("PAR_file")
        self.EISgrEvRHE_data_grp = EISgrEvRHE_data.groupby(["PAR_file"])

    #    def prepare_fit_run_args(self):
    #        ovv_exp_grp_PF = self.EIS_ovv_destfiles.groupby('PAR_file')
    #        EISgrEvRHE_data_grp = self.EISgrEvRHE_data.groupby(['PAR_file','Segment #',EvRHE, 'RPM_DAC'])

    #        fit_run_args = [EIS_spectrum_arg(Path(PF),int(Seg), np.round(float(E_V),3),\
    #                        int(RPM_DAC), gr, ovv_exp_grp_PF.get_group(PF), self.EIS_kwargs)
    #                        for (PF, Seg, E_V, RPM_DAC),gr in EISgrEvRHE_data_grp]
    #        self.fit_run_args = fit_run_args
    def create_fit_rung_args(self):
        pass

    def collection_fit_run_args(self):
        pass

    def EIS_collector_gen(self):

        _EIS_fit_run = (
            EIS_spectra_collection(
                PF,
                [
                    EIS_Spectrum(
                        Path(PF),
                        int(Seg),
                        np.round(float(E_V), 3),
                        int(RPM_DAC),
                        gr,
                        self.ovv_exp_grp_PF.get_group(PF),
                        self.EIS_kwargs,
                    )
                    for (PF, Seg, E_V, RPM_DAC), gr in PFgrp.groupby(
                        ["PAR_file", "Segment #", EvRHE, "RPM_DAC"]
                    )
                ],
                PFgrp,
                self.ovv_exp_grp_PF.get_group(PF),
                self.EIS_kwargs,
            )
            for PF, PFgrp in self.EISgrEvRHE_data_grp
        )
        return _EIS_fit_run

    def PF_group_gen(self):
        #        ovv_exp_grp_PF = self.EIS_ovv_destfiles.groupby('PAR_file')
        #        EISgrEvRHE_data_grp = self.EISgrEvRHE_data.groupby(['PAR_file','Segment #',EvRHE, 'RPM_DAC'])
        #        fit_run_args = [Meta(Path(PF),int(Seg), np.round(float(E_V),3),\
        #                        np.round(float(E_V)*1E3,3), int(RPM_DAC), gr, ovv_exp_grp_PF.get_group(PF))
        #                        for (PF, Seg, E_V, RPM_DAC),gr in EISgrEvRHE_data_grp]
        for PF, PFgrp in self.EISgrEvRHE_data_grp:
            yield EIS_spectra_collection(
                PF,
                [
                    EIS_Spectrum(
                        Path(PF),
                        int(Seg),
                        np.round(float(E_V), 3),
                        int(RPM_DAC),
                        gr,
                        self.ovv_exp_grp_PF.get_group(PF),
                        self.EIS_kwargs,
                    )
                    for (PF, Seg, E_V, RPM_DAC), gr in PFgrp.groupby(
                        ["PAR_file", "Segment #", EvRHE, "RPM_DAC"]
                    )
                ],
                PFgrp,
                self.ovv_exp_grp_PF.get_group(PF),
                self.EIS_kwargs,
            )


#    def __next__(self):
#        return next(self.__iter__())
#        fit_run_args = [EIS_spectrum_arg(Path(PF),int(Seg), np.round(float(E_V),3),\
#                        int(RPM_DAC), gr, ovv_exp_grp_PF.get_group(PF), self.EIS_kwargs)
#                        for (PF, Seg, E_V, RPM_DAC),gr in EISgrEvRHE_data_grp]
#        self.fit_run_args = fit_run_args


class EIS_run_loop:
    def __init__(self, _obj):
        if "ECRunOVV" in type(_obj).__name__:
            self.EISrunOVV = _obj
            self.EIS_Prep = EIS_Preparator(self.EISrunOVV)

    #        self._gen_loop = _gen

    def run_loop(self):
        for EIS_spec_collection in self.EIS_Prep.EIS_collector_gen():
            EIS_spec_collection
            _mean_spec = Fit_Spectrum(EIS_spec_collection.mean_spectrum)
            _mean_spec.linKK()

            EIS_spec_colection
            print(EIS_spec_colection)


#
# class EIS_input_run(ECRunOVV):
#
#    def __init__(self, *args, **kwargs):
#
#        _t = [n for n,i in enumerate(args) if isinstance(i, ECRunOVV)]
#        if _t:
#            self._ECRun = args[_t[0]]
#            instance_attrs = vars(self._ECRun)
#            super().__init__(**instance_attrs)
#

#        pass
#        if any(isin)
#        super.()


def eis_run_group_ovv(ovv_exp_grp, **EIS_kwargs):
    # EIS_use_prepars = False
    #        EIS_use_prepars_set, FreqLim_set,EIS_skip, EIS_plot_combined = True , 50000, False, True
    ###### === Analyze the EIS experiments ==== #######
    #                exportEIS_option = 'Text,Plot'
    #    EIS_skip = EIS_kwargs.get('EIS_skip',False)

    #        EIS_plot_combined = False
    # %%
    #    EIS_ovv = ovv.query('PAR_exp == "EIS"')
    #    if not ovv.empty or EIS_skip == False:
    #                EIS_scans = All_ovv.query('EXP == "EIS"')
    #                EIS_scans = ovv.query('PAR_exp == "EIS"')
    #                for nm2 ,EISgrEvRHE in EIS_scans.groupby('File'):
    #        index_out = []
    #        Meta = namedtuple('Meta', 'PAR_file E_dc_RHE E_dc_RHE_mV RPM_DAC data ovv')
    #        EIS_fit_kwargs.update({'meta_template' : Meta})
    #    ovv_exp_grp.groups.keys()
    #    EIS_ovv = pd.concat([gr for n,gr in ovv_exp_grp if 'EIS' in n])

    if "BadOnly" in EIS_fit_kwargs.keys():
        _bak_fit_run_args = fit_run_args
        filter_run_args = EIS_fit_kwargs.get("BadOnly", [])
        fltr_fit_args = [
            (Path(i[0]), int(i[1]), *i[2:4], int(i[4])) for i in filter_run_args
        ]
        fit_run_args = list(filter(lambda x: x[0:5] in fltr_fit_args, fit_run_args))
        print(f"bak:{len(_bak_fit_run_args)}, len(fit):{len(fit_run_args)}")
    #        [i for i in fit_run_args if (str(i[0]),i[1:-1]) in [i[:-1] for i in filter_run_args]]

    fit_kwargs_iter = repeat(EIS_fit_kwargs)
    fit_export_EV_all = []
    multi_par_fit = True
    _test_fitargs, fit_testing = (), False
    if "error_file" in EIS_fit_kwargs.get("input_run", "n"):
        #        '12.02.2019_0.1MKOH_cell2/O2_EIS-range_1500rpm_JOS6_postAST_899 at 0.999 (1500 rpm)'
        #        _test_fitargs = [i for i in fit_run_args if np.round(i[1],1) in [1] ] #fit_run_args[0]fit_run_args
        #        _test_fitargs =  #fit_run_args[0]fit_run_args
        _test_fitargs = fit_run_args[0:4]
        #        [i for i in fit_run_args if 'NSS-0103-EIS-7' in i[0].name and i[1] == 3]
        "O2_EIS-range_1500rpm_JOS2_285_705mV_1500rpm_4_spectrumfit_v20"
        "O2_EIS-range_1500rpm_JOS5_285_305mV_1500rpm_10_linkK_v20"
        fit_run_arg = fit_run_args[18]  # _test_fitargs[0]
        fit_run_arg = [i for i in fit_run_args if "O2" in i[0].name and 12 == i[1]][0]
    if "1500rpm" in EIS_fit_kwargs.get("input_run", "n"):
        _fit_run_args_bak = fit_run_args
        fit_run_args = [i for i in fit_run_args if i.RPM_DAC > 1000]
    #        multi_par_fit, fit_testing = False, True
    #%%

    for PF, PF_run_arg in PF_grp_fit_run_args:
        a, b = PF, list(PF_run_arg)

    if multi_par_fit:
        fit_export_EV_all = []
        try:
            #            EISgrEvRHE_data_grp = EISgrEvRHE_data_raw.groupby(['PAR_file',EvRHE, 'RPM_DAC'])
            #            fit_run_args = [Meta(Path(PF), np.round(float(E_V),3), np.round(float(E_V)*1E3,3), int(RPM_DAC), gr, ovv_exp_grp_PF.get_group(str(PF))) for (PF,E_V,RPM_DAC),gr in EISgrEvRHE_data_grp]
            #            fit_kwargs_iter = repeat(EIS_fit_kwargs)
            #            args_for_starmap = zip(repeat(eis_fit_PAR_file), t_run_args, t_kwargs_iter)
            #            EIS_run_grps = [(nm,gr,EIS_kwargs) for nm,gr in EIS_ovv.groupby(by='PAR_file')]
            logger.info(
                f"EIS eis_run_group_ovv START multiprocessing {multi_par_fit} for len{len(fit_run_args)}"
            )
            pool_size = os.cpu_count() - 2
            if "linux" in sys.platform:
                #                os.system('taskset -cp 0-%d %s' % (pool_size, os.getpid()))
                os.system("taskset -p 0xff %d" % os.getpid())
            #                os.sched_setaffinity(0,{i for i in range(pool_size)})

            #            for chunk in fit_run_args[::pool_size if pool_size > 2 else 1]:
            #                print(len(chunk))
            with multiprocessing.Pool(pool_size) as pool:
                fit_export_EV_all_chunck = PF_fit_starmap_with_kwargs(
                    pool, fit_EEC, fit_run_args, fit_kwargs_iter
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
        #    for (PF,E_V,RPM_DAC),EISgr_data_EV in EISgrEvRHE_data_raw.groupby(['PAR_file',EvRHE, 'RPM_DAC']):
        #            pass#test
        #            fit_run_arg  = Meta(Path(PF), np.round(float(E_V),3), np.round(float(E_V)*1E3,3), int(RPM_DAC), EISgr_data_EV, ovv_exp_grp_PF.get_group(str(PF)))
        #            fit_run_arg = Meta(Path(nm2), np.round(float(E_V),3), np.round(float(E_V)*1E3,3), int(RPM_DAC), EISgr_data_EV, gr_EIS_ovv)
        if fit_testing:
            if not _test_fitargs:
                print("E_dc : \n", [np.round(i[2], 2) for i in fit_run_args])
                _test_fitargs = [
                    i
                    for i in fit_run_args
                    if np.round(i[2], 2) in [0.7] and Path(i[0]).name in ["JOS1", "O2"]
                ]  # fit_run_args[0]
                _test_fitargs = [
                    i
                    for i in fit_run_args
                    if "JOS" in i[0].name
                    and "O2" in i[0].name
                    and np.round(i[2], 1) in [0.1]
                ]
                _test_fitargs = [
                    i
                    for i in fit_run_args
                    if "JOS3" in i[0].name
                    and "O2" in i[0].name
                    and np.round(i[2], 3) in [0.805]
                ]
                _test_fitargs = [
                    i
                    for i in fit_run_args
                    if "JOS4" in i[0].name
                    and "N2" in i[0].name
                    and np.round(i[2], 3) in [0.758]
                ]
                _test_fitargs = [
                    i
                    for i in fit_run_args
                    if "JOS3" in i[0].name
                    and "N2" in i[0].name
                    and np.round(i[2], 3) in [0.758]
                ]
                _test_fitargs = [
                    i
                    for i in fit_run_args
                    if "JOS4" in i[0].name and "O2" in i[0].name and i.RPM_DAC == 1500
                ]

            fit_run_arg = _test_fitargs[25]
        #             # [0.21,0.76,0.91]
        #        a = fit_EEC(_test_fitargs[0], **EIS_fit_kwargs)
        logger.warning(f"EIS starting single core loop fitting {multi_par_fit}")
        fit_export_EV_all = []
        for fit_run_arg in fit_run_args_tests:
            fit_export_EV = fit_EEC(fit_run_arg, **EIS_fit_kwargs)
            fit_export_EV_all.append([fit_export_EV])
    #    else:
    #         logger.error(f'EIS no fitting {multi_par_fit}')
    try:
        fit_export_EV_all_flat = [a for i in fit_export_EV_all for a in i]
        PF_fit_spectra_all_grp = pd.concat(
            [i[0] for i in fit_export_EV_all_flat], sort=False
        ).groupby(["PAR_file"])
        PF_fit_pars_all_grp = pd.concat(
            [i[1] for i in fit_export_EV_all_flat], sort=False
        ).groupby(["PAR_file"])

        for PF, PF_pars in PF_fit_pars_all_grp:
            dest_PF = ovv_exp_grp_PF.get_group(PF)
            destPars = dest_PF["EIS_dest_Pars"].iloc[0]
            Parsout_path_target = FileOperations.CompareHashDFexport(PF_pars, destPars)
            if "models_pars" in EIS_fit_kwargs.get("input_run", ""):
                for mod, mgrp in PF_pars.groupby("Model_EEC"):
                    #                 mod_dest_pf = destPars.with_name(f'{destPars.stem}_{mod}.xlsx').name
                    mod_dir = destPars.parent.joinpath(mod)
                    mod_dir.mkdir(parents=True, exist_ok=True)
                    mod_dest_pf = mod_dir.joinpath(
                        destPars.with_name(f"{destPars.stem}_{mod}.xlsx").name
                    )
                    mod_target = FileOperations.CompareHashDFexport(mgrp, mod_dest_pf)
                    var_names = mgrp.lmfit_var_names.iloc[0].split(", ")
            #                 for var in var_names:
            #                     fig,ax = plt.subplots()
            #                     mgrp.plot(x='Segment' , y= var,ax=ax,kind='scatter')
            #                     plt.savefig(mod_dest_pf.with_name(f'{var}_{mod_dest_pf.stem}.png'),dpi=100,bbox_inches='tight')
            ##                     plt.show()
            #                     plt.close()
            E_data_combined_path_target = FileOperations.CompareHashDFexport(
                PF_fit_spectra_all_grp.get_group(PF),
                dest_PF["EIS_dest_spectra"].iloc[0],
            )
            EIS_plotting_EvRHE(
                PF_fit_spectra_all_grp.get_group(PF),
                ovv_exp_grp_PF.get_group(PF),
                PF_pars,
                dest_PF["EIS_dest_spectra"].iloc[0],
            )
    except Exception as e:
        logger.error(
            f"EIS finish plotting error: {e}, PAR files {len(ovv_exp_grp_PF)} , spectra: {len(PF_fit_spectra_all_grp)}"
        )

    return logger.info(
        f"EIS FINISHED: PAR files {len(ovv_exp_grp_PF)} , spectra: {len(EISgrEvRHE_data_grp)}"
    )


#            with multiprocessing.Pool(os.cpu_count()-2) as pool:
#                try:
##                    index_out = pool.starmap(eis_multi_wrap(eis_fit_PAR_file, EIS_run_grps))
#                    starmap_with_kwargs(pool, eis_fit_PAR_file, t_run_args, t_kwargs_iter)
#    #                results = pool.map(partial(eis_fit_PAR_file, 'This is arg1', **just_a_dict), list_of_args2)
#                except Exception as e:
#                    print('eis_run_group_ovv multiprocessing error: {0}'.format(e))
#                    logger.error('eis_run_group_ovv  multiprocessing error: {0}'.format(e))
#
#        except Exception as e2:
#            print('EIS eis_run_group_ovv multiprocessing error: {0}'.format(e2))
#            logger.error('EIS eis_run_group_ovv  multiprocessing erroe: {0}'.format(e2))
#    else:
#        for nm2, gr_EIS_ovv in ovv.groupby(by='PAR_file'):
##                nm2
##                index_out_pf =
#            eis_fit_PAR_file(nm2,gr_EIS_ovv, **EIS_fit_kwargs)
##                index_out.append(index_out_pf)
#            index_out = [a for i in index_out for a in i]
#         EIS_ovv.groupby(by='PAR_file')
#         with multiprocessing.Pool(os.cpu_count()-1) as pool:
#            try:
#                results = pool.starmap(partial(eis_fit_PAR_file, run_groups)
#                results = pool.map(partial(eis_fit_PAR_file, 'This is arg1', **just_a_dict), list_of_args2)
#            except Exception as e:
#                print('PAR DW Classifier multiprocessing error: {0}'.format(e))
#                logger.error('PAR DW Classifier multiprocessing error: {0}'.format(e))
# %%
#    try:
#        indexes = pd.concat(index_out, sort=False, ignore_index=True)
##        [a for i in index_out for a in i]
##        breakpoint()
##        print('Indexes:',indexes)
#        if indexes.empty:
#            logger.error(f'EIS eis_run_group_ovv indexes empty')
#        indexes['Script_run_date'] = datetime.now()
#        indexes_path = Path(gr_ovv.Dest_dir.unique()[0]).joinpath(
#            '{0}_index_{1}.xlsx'.format(exp, indexes.DestFile.nunique()))
#        indexes_path_target = FileOperations.CompareHashDFexport(indexes, indexes_path)
#    except Exception as e:
#        logger.error(f'EIS eis_run_group_ovv indexes out: {e}')
#    return indexes, indexes_path_target


#
# def eis_fit_PAR_file(nm2,gr_EIS_ovv, **EIS_fit_kwargs):
##    print(EIS_fit_kwargs)
#
#    EvRHE = 'E_AppV_RHE'
#
#    Parsout_lst, nm2 = [], Path(nm2)
#    index_per_file_EV_lst = []
#    OutDataCombi, AllData_E_file = [], pd.DataFrame()
#    EISgrEvRHE_data_raw, EISgrEvRHE_actions = create_CVs(gr_EIS_ovv)
#    EISgrEvRHE_data = EISgrEvRHE_data_raw.query('ActionId == 21')
#    EISgrEvRHE_data[EvRHE] = EISgrEvRHE_data[EvRHE].round(3)
#
#    if len(EISgrEvRHE_data_raw) - len(EISgrEvRHE_data) > 0:
#        print('EIS filtered data because ActionId != 21 but in {0}'.format(EISgrEvRHE_data_raw.ActionId.unique()))
#    EIS_dest_dir = Path(gr_EIS_ovv.Dest_dir.iloc[0]).joinpath('EIS')
#    EIS_dest_dir.mkdir(parents=True, exist_ok=True)
#    EIS_fit_kwargs.update({'EIS_dest_dir' : EIS_dest_dir})
#    #                    dest_dir.joinpath('EIS')
#    PF_fit_pars_target = EIS_dest_dir.joinpath(nm2.stem + '_pars.xlsx')
##    Parsout_path_target = FileOperations.CompareHashDFexport(pd.DataFrame(), Parsout_path)
#    PF_fit_spectra_target = EIS_dest_dir.joinpath(nm2.stem + '_Combined.xlsx')
#    index_out = []
##    fit_kwargs = dict(**dict(perform_prefit=True, exportEIS=EIS_single_output)+**EIS_kwargs)
##    Meta = namedtuple('Meta', 'PAR_file E_dc_RHE E_dc_RHE_mV RPM_DAC data ovv')
##    templ = EIS_fit_kwargs.get('meta_template')
#    fit_export_EV_all = [] # becomes list of namedtuples
#    multi_par_fit = False
#    if multi_par_fit:
#        try:
#            EISgrEvRHE_data_grp = EISgrEvRHE_data.groupby([EvRHE, 'RPM_DAC'])
#            fit_run_args = [Meta(Path(nm2), np.round(float(E_V),3), np.round(float(E_V)*1E3,3), int(RPM_DAC), gr, gr_EIS_ovv) for (E_V,RPM_DAC),gr in EISgrEvRHE_data_grp]
#            fit_kwargs_iter = repeat(EIS_fit_kwargs)
#
#    #            args_for_starmap = zip(repeat(eis_fit_PAR_file), t_run_args, t_kwargs_iter)
#    #            EIS_run_grps = [(nm,gr,EIS_kwargs) for nm,gr in EIS_ovv.groupby(by='PAR_file')]
#            logger.info(f'EIS eis_fit_PAR_file START multiprocessing for {nm2}')
#            with multiprocessing.Pool(len(EISgrEvRHE_data_grp)) as pool:
#                fit_export_EV_all = PF_fit_starmap_with_kwargs(pool, fit_EEC, fit_run_args, fit_kwargs_iter)
#    #                results = pool.map(partial(eis_fit_PAR_file, 'This is arg1', **just_a_dict), list_of_args2)
#    #                print('eis_run_group_ovv multiprocessing error: {0}'.format(e))
#    #                logger.error('eis_run_group_ovv  multiprocessing error: {0}'.format(e))
#
#        except Exception as e2:
#    #        print('EIS eis_run_group_ovv multiprocessing error: {0}'.format(e2))
#            logger.error('EIS eis_fit_PAR_file  multiprocessing erroe: {0}'.format(e2))
#            multi_par_fit = False
#
#    if multi_par_fit == False:
#        for (E_V,RPM_DAC),EISgr_data_EV in EISgrEvRHE_data.groupby([EvRHE, 'RPM_DAC']):
##            if EISgr_data_EV.empty or gr_EIS_ovv.empty:
##                print('ERROR: EIS DATA or OVV EMPTY')
##            continue
##    PAR_file = E_dc_RHE E_dc_RHE_mV = E_dc_RHE * 1000 RPM_DAC =
#            fit_run_arg = Meta(Path(nm2), np.round(float(E_V),3), np.round(float(E_V)*1E3,3), int(RPM_DAC), EISgr_data_EV, gr_EIS_ovv)
##        meta_info = Meta(Path(nm2), np.round(float(E_V),3), np.round(float(E_V)*1E3,3), int(RPM_DAC))
#            fit_export_EV = fit_EEC(fit_run_arg, **EIS_fit_kwargs)
#            fit_export_EV_all.append(fit_export_EV)
#
#    PF_fit_pars = pd.concat([i.fit_pars for i in fit_export_EV_all], sort=False)
#    PF_fit_spectra = pd.concat([i.fit_spectra for i in fit_export_EV_all], sort=False)
#
#    Parsout_path_target = FileOperations.CompareHashDFexport(PF_fit_pars, PF_fit_pars_target)
#    E_data_combined_path_target = FileOperations.CompareHashDFexport(PF_fit_spectra, PF_fit_spectra_target)
#    EIS_plotting_EvRHE(PF_fit_spectra, gr_EIS_ovv, PF_fit_pars, PF_fit_spectra_target)

#                    Parsout_path_target = Parsout_path_target.with_suffix('.TESTING.xlsx' )
#    AllData_E_file = pd.concat(OutDataCombi, sort=False)
#                        AllData_E_fn = EIS_dest_dir.joinpath(nm2.stem+'_Combined.xlsx')


#            Parsout_lst.append(pars_models), OutDataCombi.append(EIS_fit_data), index_per_file_EV_lst.append(indexes_per_EV_out)

#            if 'Plot' in EIS_fit_kwargs.get('EIS_single_output','Text, Plot'):
#                EIS_outPath_target_png = pars_models.FitDataFile.unique()[0].with_suffix('.png').with_suffix('.png')
#                index_dataplot_output = {'PAR_file': nm2,'Type_output' : 'EIS_fit_data_png', 'Type_exp' : 'EIS', 'DestFile' : EIS_outPath_target_png, 'E_V' : meta_info.E_dc_RHE,'RPM_DAC' : meta_info.RPM_DAC}
#                EIS_plotting_per_EV(EIS_fit_data,pars_models, EIS_outPath_target_png,plot_show = False )
#                indexes_per_EV_out.append(index_dataplot_output)
#            breakpoint()
#                            (EISgr_data_EV,gr_EIS_ovv,InitP,EIS_dest_dir,perform_prefit = True,TrimData = False,FitOnlyTrimmedData = False, FreqLim = 10000,exportEIS = 'Text,Plot'):
#                                print('PAR EIS Succes:','\n', InitP1, InitP2,'\n %.2f'%nm1,nm2)
#        except Exception as e:
#                                print(gr_EIS_ovv,'\n', InitP1,'\n', InitP2,'\n %.2f'%nm1,nm2)
#            print('PAR EIS fit failed: %s \n %s' % (e, nm2))
#        1    EIS_fit_data, pars_models, indexes_per_EV_out = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
#                                if plotEIS:
#                                    try:
#                                        PAR_EIS_fit.EIS_plotting(EISgr_data_EV,gr_EIS_ovv,EIS_dest_dir,outData,out1,out2)
#                                    except Exception as e:
#                                        pass
#        Parsout_lst.append(pars_models), OutDataCombi.append(EIS_fit_data), index_per_file_EV_lst.append(
#            indexes_per_EV_out)


#    try:
#        index_out.append(pd.concat([pd.DataFrame(i) for i in index_per_file_EV_lst], sort=False, ignore_index=True))
#
#        Parsout = pd.concat(Parsout_lst, sort=False)
#        #                    Parsout_path_target = Parsout_path_target.with_suffix('.TESTING.xlsx' )
#        Parsout_path_target = FileOperations.CompareHashDFexport(Parsout, Parsout_path)
##        Parsout.to_excel(Parsout_path_target)
#        index_info_pars = pd.DataFrame(
#            {'PAR_file': nm2, 'DestFile': Parsout_path_target, 'Type_output': 'EIS_Pars', 'Type_exp': 'EIS',
#             'nModels': Parsout.Model_EEC.nunique(), 'Models_used': Parsout.Model_EEC.unique()})
#        index_out.append(index_info_pars)
#
#        #                    Pars2out_path_target = FolderOps.FileOperations.CompareHashDFexport(Pars2 ,Pars2out_path)
#        #                    index_info_p2 = {'PAR_file'  : nm2,'DestFile' : Pars2out_path_target,'Type_output' : 'EIS_Pars2', 'Type_exp' : 'EIS'}
#        #                    index_out.append(index_info_p2)
#        #                                Pars2.to_excel(EIS_dest_dir.joinpath(nm2.stem+'_pars2.xlsx'))
#        AllData_E_file = pd.concat(OutDataCombi, sort=False)
#        #                        AllData_E_fn = EIS_dest_dir.joinpath(nm2.stem+'_Combined.xlsx')
#        E_data_combined_path_target = FileOperations.CompareHashDFexport(AllData_E_file, E_data_combined_path)
#        index_info_all = pd.DataFrame(
#            {'PAR_file': nm2, 'DestFile': E_data_combined_path_target, 'Type_output': 'EIS_AllData_combined',
#             'Type_exp': 'EIS', 'nModels': AllData_E_file.Model_EEC.nunique(),
#             'Models_used': AllData_E_file.Model_EEC.unique()})
#        index_out.append(index_info_all)
#        #                                AllData_E_file.to_excel(AllData_E_fn)
#
#        if EIS_fit_kwargs.get('EIS_plot_combined',True) == True:
#            try:
#                EIS_plotting_EvRHE(AllData_E_file, gr_EIS_ovv, Parsout, E_data_combined_path_target)
#            except Exception as e:
#                logger.warning('ERROR: EIS Plot combined output for %s: %s' % (nm2, e))
#        #                    print('EIS exported success %s' %Path(nm2.parent.name).joinpath(nm2.name))
#        logger.info('EIS exported success %s' % Path(nm2.parent.name).joinpath(nm2.name))
#    except Exception as e:
#        #                    print('ERROR: EIS Pars output for %s: %s' %(nm2,e))
#        logger.warning('ERROR: EIS Pars output for %s: %s' % (nm2, e))
#    return index_out

if __name__ == "__main__":
    pass
#    try:
#        t_run_args = [(nm,gr) for nm,gr in EIS_ovv.groupby(by='PAR_file')]
#        t_kwargs_iter = repeat(EIS_kwargs)
#    #            args_for_starmap = zip(repeat(eis_fit_PAR_file), t_run_args, t_kwargs_iter)
#        EIS_run_grps = [(nm,gr,EIS_kwargs) for nm,gr in EIS_ovv.groupby(by='PAR_file')]
#        with multiprocessing.Pool(os.cpu_count()-2) as pool:
#            try:
#    #                    index_out = pool.starmap(eis_multi_wrap(eis_fit_PAR_file, EIS_run_grps))
#                index_out = starmap_with_kwargs(pool, eis_fit_PAR_file, t_run_args, t_kwargs_iter)
#    #                results = pool.map(partial(eis_fit_PAR_file, 'This is arg1', **just_a_dict), list_of_args2)
#            except Exception as e:
#                print('PAR DW Classifier multiprocessing error: {0}'.format(e))
#                logger.error('PAR DW Classifier multiprocessing error: {0}'.format(e))
#
#    except Exception as e2:
#        print('EIS multiprocessing error: {0}'.format(e2))
#        logger.error('EIS multiprocessing erroe: {0}'.format(e2))

#    Pars1, Pars2 = Pars1.append(pars1),Pars2.append(pars2)
#                                outP1,outP2 = outP1.append(pars1),outP2.append(pars2)
#                   outData.to_hdf(PathDB,Path(dest_dir.parts[-2],dest_dir.parts[-1]).joinpath(('%.2f_Vrhe'%nm1)).as_posix(),table=True,mode='a')
