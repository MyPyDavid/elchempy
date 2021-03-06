# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:42:36 2020

@author: User
"""

# import sys
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

# from functools import partial
from itertools import repeat
import pandas as pd


from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations


from EIS.fitting import Fit_Spectrum, Fit_Spectra_Collection
from EIS.models import Model_Collection

# prepare_rand_models

from EC_DataLoader.CreateCV import create_CVs as create_CVs

# from EIS.plotting import (
#     EIS_plotting_per_EV,
#     EIS_Trimming_plot,
#     EIS_plotting_EvRHE,
# )


import logging

_logger = logging.getLogger(__name__)

globals()["EvRHE"] = "E_AppV_RHE"

# print('TODO: eis run fix MULTIPROCESSING!')
# TODO: do it
# Meta = namedtuple('Meta', 'PAR_file Segment E_dc_RHE E_dc_RHE_mV RPM_DAC data ovv')
@dataclass(order=True, frozen=False)
class EIS_Spectrum:
    """EIS Spectrum dataclass.\n
    Holds spectrum raw data in pd.DataFrame and metadata"""

    _required_cols = ["Frequency(Hz)", "Z Real", "Z Imag"]
    _spectrum_grp_cols = ["PAR_file", "Segment #", EvRHE, "RPM_DAC"]

    PAR_file: Path = field(default=Path(Path.cwd().joinpath("empty.txt")))
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
        _file = f'File: "{self.PAR_file}"'
        _data = f"at {self.E_dc_RHE_mV} mV with {self.RPM_DAC} rpm, data({len(self.data)}) and ovv({len(self.ovv)})"
        _keys = f"attrs: {self.E_dc_RHE_mV} mV with {self.RPM_DAC} rpm, data({len(self.data)}) and ovv({len(self.ovv)})"
        return _file + "\n" + _data + "\n" + _keys

    def add_complex_impedance_columns(self):
        #        ['Frequency(Hz)'
        _check_cols = []
        if hasattr(self.data, "columns"):
            _check_cols = [i in self.data.columns for i in self._required_cols]

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
            # DataWeights_modulus_Z = np.sqrt((Zre**2+Zim**2))
            # 'lmfit_weights_mod_Z' : DataWeights_modulus_Z
            DataWeights_modulus_Y = np.sqrt((Zre ** 2 + Zim ** 2)) ** -1

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
                    "lmfit_weights_mod_Y": DataWeights_modulus_Y,
                    "lmfit_weights_unit": Zre / Zre,
                    "lmfit_weights_prop": Ydata,
                }
            )

            _meta_info_cols = {
                key: val
                for key, val in vars(self).items()
                if any([i in type(val).__name__ for i in ["int", "float", "Path"]])
            }
            _add_cols.update(_meta_info_cols)
            EIS_data_raw = pd.DataFrame(_add_cols, index=self.data.index)
            EIS_data_raw = EIS_data_raw.sort_values(by="Frequency(Hz)", ascending=True)
        elif self.data.empty:
            if not "empty" in self.PAR_file.name:
                _logger.error(f"Error in EIS_spectrum, empty data for {self}")
        else:
            _logger.error(
                "Error in EIS_spectrum, missing columns:",
                ", ".join(
                    [i for i in self.data.columns if i not in self._required_cols]
                ),
            )
            raise ValueError
        self.EIS_data = EIS_data_raw

    def check_freqlim(self):
        if not self.EIS_data.empty:
            FreqLim = self.EIS_kwargs.get("FreqLim", 30e3)
            self.EIS_data.loc[
                (
                    (self.EIS_data["Frequency(Hz)"] >= FreqLim)
                    | (self.EIS_data["DATA_Zre"] < 0)
                    | (np.abs(self.EIS_data["DATA_Zre"]) > 1e4)
                    | (np.abs(self.EIS_data["DATA_Zim"]) > 1e4)
                ),
                "Valid",
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

        _mean_E_selection = [
            i
            for i in self.spectra
            if i.E_dc_RHE > 0.4 and i.E_dc_RHE < 0.85 and "O2" in self.PAR_file.name
        ]
        if _mean_E_selection:
            _PF_data = pd.concat(i.data for i in _mean_E_selection)
        else:
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
        _mean_EIS_kwargs = {"linKK_trimming_factor": 3, "res_max": 0.2}
        self.EIS_kwargs.update(_mean_EIS_kwargs)

        _PF_data_mean_grp = _PF_data_mean.groupby(EIS_Spectrum._spectrum_grp_cols)
        self.mean_spectra = [
            EIS_Spectrum(
                Path(PF),
                int(Seg),
                np.round(float(E_V), 3),
                int(RPM_DAC),
                gr,
                _PF_mean_ovv,
                self.EIS_kwargs,
            )
            for (PF, Seg, E_V, RPM_DAC), gr in _PF_data_mean_grp
        ]
        self.mean_spectrum = self.mean_spectra[0]
        self.meta_data = _PF_meta_mean

    def __repr__(self):
        _sp = f"EIS_spectra_collection ({len(self.spectra)}), {self.spectra[0]}"
        return _sp


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
                    models = {}  # FIXME
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
        _logger.info(
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
        # print('starting prep_random')
        # print('finished prep_random')
        # print("start get_eis_suggestions")
        self.EIS_kwargs.update(get_eis_suggestions())
        print("start get_recheck_parfiles_only")
        if False:  # FIXME
            # _prepare_rand_params = prepare_rand_models(
            #     rand_params_dest_dir=FindExpFolder("VERSASTAT").PostDir.joinpath("EIS")
            # )
            # self.EIS_kwargs.update(**dict(random_params=_prepare_rand_params))
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

    # def EIS_collector_gen(self):
    def __iter__(self):
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

    def __len__(self):
        return len(self.EISgrEvRHE_data_grp)


#    def __next__(self):
#        return next(self.__iter__())
#        fit_run_args = [EIS_spectrum_arg(Path(PF),int(Seg), np.round(float(E_V),3),\
#                        int(RPM_DAC), gr, ovv_exp_grp_PF.get_group(PF), self.EIS_kwargs)
#                        for (PF, Seg, E_V, RPM_DAC),gr in EISgrEvRHE_data_grp]
#        self.fit_run_args = fit_run_args


class EIS_run_loop:
    def __init__(self, _obj):
        if "ECRunOVV" in _obj.__class__.__name__:
            self.EISrunOVV = _obj
            self.EIS_Prep = EIS_Preparator(self.EISrunOVV)
            _logger.warning(
                f"{self.__class__.__name__} starting with {len(self.EIS_Prep)}"
            )
        else:
            raise TypeError

        if (
            not "test" in self.EISrunOVV.input_run
            and self.EISrunOVV.input_run.endswith("y")
        ):
            self.run_loop()
        else:
            _logger.error(
                f'EIS run loop finished because {self.EISrunOVV.input_run} has "test" or not endswith "y"'
            )

    #        self._gen_loop = _gen
    def run_loop(self):

        if self.EISrunOVV.run_kwargs.get("multi_par_fit", False):
            _logger.error(
                f"Starting multi {self.__class__.__name__}  {len(self.EIS_Prep)}"
            )
            try:
                pool_size = os.cpu_count() - 2
                with multiprocessing.Pool(pool_size) as pool:
                    # _fit_spec_coll_result = pool.map(Fit_Spectra_Collection, iter(self.EIS_Prep))
                    for n, r in enumerate(
                        pool.imap_unordered(Fit_Spectra_Collection, self.EIS_Prep)
                    ):
                        _logger.error(
                            f"Progression {n+1}/{len(self.EIS_Prep)} {r} {self.__class__.__name__} "
                        )
                        # _logger.error(f'multi finished {self.__class__.__name__} {len(self.EIS_Prep)}')
                    # for n,EIS_spec_collection in enumerate(self.EIS_Prep):

                    # pool.apply_async((Fit_Spectra_Collection(EIS_spec_collection)))
                    # _fit_spec_coll_result = pool.apply(Fit_Spectra_Collection, self.EIS_Prep)
                _logger.error(
                    f"multi finished {self.__class__.__name__} {len(self.EIS_Prep)}"
                )
            except Exception as e:
                _logger.error(
                    f"Error multi {self.__class__.__name__} {len(self.EIS_Prep)}, {e}"
                )
        else:
            _fit_spec_coll_result = []
            for EIS_spec_collection in self.EIS_Prep:
                EIS_spec_collection
                try:
                    _logger.error(
                        f"Starting {self.__class__.__name__} {EIS_spec_collection}"
                    )
                    _fit_spec_coll = Fit_Spectra_Collection(
                        EIS_spec_collection, run_fit_mean="lmfit"
                    )
                    _fit_spec_coll_result.append(_fit_spec_coll)

                except Exception as e:
                    _logger.error(
                        f"Error  {self.__class__.__name__} {EIS_spec_collection} {len(self.EIS_Prep)}, {e}"
                    )
