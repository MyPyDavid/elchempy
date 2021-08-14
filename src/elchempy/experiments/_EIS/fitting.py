# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:00:16 2020

@author: DWXMG
"""

# import importlib.util
# import sys
# from collections import namedtuple, OrderedDict
from pathlib import Path
import itertools

# from itertools import combinations
# import operator
import datetime as dt
import random
from math import pi

# import copy
# import pickle
import pandas as pd
import numpy as np
from cmath import phase
import matplotlib.pyplot as plt
from matplotlib.offsetbox import (
    AnchoredOffsetbox,
    DrawingArea,
    HPacker,
    TextArea,
)  # TODO
from scipy import stats
from scipy.stats import linregress
from scipy.optimize import curve_fit

# from lmfit import Parameters, conf_interval, Minimizer, minimize

# import numdifftools
# import corner
from file_py_helper.file_functions import FileOperations

from .validation import get_KKvalid, prep_GP_DRT_raw_data
from .models import Model_Collection
from .plotting import (
    plot_linKK,
    EIS_Trimming_plot,
    EIS_plotting_per_EV,
    plot_lin_Warburg,
)
from .DRT_DP_fitting import DP_DRT_analysis
from .GP_DRT_fitting import run_GP_DRT_fit

#    from .plotting import EIS_plotting_per_EV, EIS_Trimming_plot, EIS_plotting_EvRHE
#    _logger = start_logging(__name__)
import logging

_logger = logging.getLogger(__name__)

# fit_export_templ = namedtuple('fit_export_templ', 'fit_spectra fit_pars meta_index')
# Meta = namedtuple('Meta', 'PAR_file Segment E_dc_RHE E_dc_RHE_mV RPM_DAC data ovv')
globals()["EvRHE"] = "E_AppV_RHE"


def func_lin(a):
    def func(x, b):
        return a * x + b

    return func


def fitting_recheck_params(fit_run_arg, modname, params_model, **EIS_fit_kwargs):
    #    PF,E,RPM = str(fit_run_arg.PAR_file), fit_run_arg.E_dc_RHE, fit_run_arg.RPM_DAC
    _key = (
        str(fit_run_arg[0]),
        int(fit_run_arg[1]),
        *[float(i) for i in fit_run_arg[2:4]],
        int(fit_run_arg[4]),
        modname,
    )
    #    (PF,E,RPM,modname)
    _get_params = pd.DataFrame()
    #     ('/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT/2019-05-May/06.05.2019_0.1MH2SO4_cell2/O2_EIS-range_1500rpm_JOS2_288.par',
    #  4,  0.708,  708.0,  1500,  'Model(Randles_RQRQ)')
    # ['PAR_file',EvRHE,'RPM_DAC','Model_EEC']
    bad_grp, good_grp = EIS_fit_kwargs.get("EIS_recheck_bad_fits"), EIS_fit_kwargs.get(
        "EIS_recheck_good_fits"
    )
    sugg_grp = EIS_fit_kwargs.get("EIS_recheck_bad_fits_suggestions")
    recheck_msg = ""
    if all([len(i.groups) > 0 for i in [bad_grp, good_grp, sugg_grp]]):
        if [i for i in good_grp.groups if _key == i]:
            recheck_msg += "Prefit recheck in good keys"
            _get_params = good_grp.get_group(_key)
        elif [i for i in bad_grp.groups if _key == i]:
            recheck_msg += f"Prefit recheck in bad keys {_key} and"
            _sugg_match = [i for i in sugg_grp.groups if _key == i]
            if _sugg_match:
                recheck_msg += " taking suggestions."
                _get_params = sugg_grp.get_group(_key)
            else:
                recheck_msg += " not in suggestions."
        #                         _logger.warning(f'Prefit recheck bad key {_key} not in suggestions')
        else:
            recheck_msg += f"Prefit recheck keys not in good or bad {_key}"
    else:
        recheck_msg += f"Prefit recheck empty frames"
    #    _logger.warning(recheck_msg)
    return _get_params, recheck_msg


#%%
def make_prefit_frame(
    EIS_data_KKvalid,
    lmfitting,
    prefix="pp",
    plot=False,
    check_err=True,
    norm=np.array([]),
    get_frame=False,
):
    #            norm = 1/(Z_KKv.real**2+Z_KKv.imag**2)
    #            abs(Z_KKv)
    #            norm = 1
    #            norm = np.sqrt(EIS_data_KKvalid.DATA_weightsmod_Z.values)
    #            norm = 1/abs(Z_KKv)
    #            norm = 1/np.sqrt(EIS_data_KKvalid.DATA_weightsmod_Z.values)
    #    lmfitting = best_trial
    # EIS_data_KKvalid, lmfitting = _spectrum.EIS_data_KKvalid,best_trial_weights # FIXME
    # make_prefit_frame(EIS_data_KKvalid, out, plot = 'Y')
    if np.array([]).any() == False:
        norm = np.array([1] * len(lmfitting.best_fit.real))
    if "DataFrame" in type(EIS_data_KKvalid).__name__:
        Z_KKv = EIS_data_KKvalid.DATA_Z.values
    elif "array" in type(EIS_data_KKvalid).__name__:
        Z_KKv = EIS_data_KKvalid
        EIS_data_KKvalid = pd.DataFrame(Z_KKv)
    EIS_data_KKvalid = EIS_data_KKvalid.loc[
        EIS_data_KKvalid.DATA_Zre.isin(lmfitting.data.real)
        & EIS_data_KKvalid.DATA_Zim.isin(lmfitting.data.imag)
    ]

    if norm.size == 0:
        norm = (
            lmfitting.data.real / lmfitting.data.real
        )  # (Z_KKv.real**2+Z_KKv.imag**2)**-1

    if not "float" in str(type(lmfitting.residual)) and lmfitting.success:
        resIm, resRe = lmfitting.residual[1::2], lmfitting.residual[0::2]
    else:
        resIm, resRe = 1e9, 1e9

    pp_errRE, pp_errIM = (lmfitting.best_fit.real - lmfitting.data.real) ** 2, (
        lmfitting.best_fit.imag - lmfitting.data.imag
    ) ** 2
    pp_errRE_mean, pp_errIM_mean = pp_errRE.mean(), pp_errIM.mean()

    #    MSE_Re,MSE_Im= (lmfitting.best_fit.real-Z_KKv.real)**2, (lmfitting.best_fit.imag-Z_KKv.imag)**2
    MSE = np.sqrt(sum(pp_errRE) + sum(pp_errIM))
    #    pp_errRE_std,pp_errIM_std  = np.abs(pp_errRE).std(), np.abs(pp_errIM).std()
    pp_Z = lmfitting.best_fit
    pp_Y = pp_Z ** -1

    prefit_data = EIS_data_KKvalid.assign(
        **{
            f"{prefix}_err_Re": pp_errRE,
            f"{prefix}_err_Im": pp_errIM,
            f"{prefix}_Z_Re": pp_Z.real,
            f"{prefix}_Z_Im": -1 * pp_Z.imag,
            f"{prefix}_Y_Re": pp_Y.real,
            f"{prefix}_Y_Im": pp_Y.imag,
            f"{prefix}_res_Re": resRe,
            f"{prefix}_res_Im": resIm,
            f"{prefix}_norm": norm,
        }
    )
    if get_frame:
        return prefit_data

    ext_ang = np.logspace(-3, 4, endpoint=True) * 2.0 * pi
    ext_model = lmfitting.eval(lmfitting.params, ang=ext_ang)
    extended_data = pd.DataFrame(
        {
            f"{prefix}_ext_Z_Re": ext_model.real,
            f"{prefix}_ext_Z_Im": -1 * ext_model.imag,
            f"{prefix}_ext_Y_Re": (ext_model ** -1).real,
            f"{prefix}_ext_Y_Im": (ext_model ** -1).imag,
            f"{prefix}_ext_freq": ext_ang / (2.0 * pi),
        }
    )

    if plot:
        fig, ax = plt.subplots(3, 1, figsize=(4, 8))
        #                if 'Y' in str(plot_pp):
        prefit_data.plot(x=f"DATA_Yre", y=f"DATA_Yim", kind="scatter", ax=ax[0])
        prefit_data.plot(x=f"{prefix}_Y_Re", y=f"{prefix}_Y_Im", c="r", ax=ax[0])
        #                else:
        prefit_data.plot(
            x=f"DATA_Zre",
            y=f"DATA_-Zim",
            kind="scatter",
            ax=ax[1],
            logy=True,
            logx=True,
        )
        prefit_data.plot(
            x=f"{prefix}_Z_Re",
            y=f"{prefix}_Z_Im",
            c="r",
            ax=ax[1],
            logy=True,
            logx=True,
        )

        if not extended_data.empty:
            extended_data.plot(
                x=f"{prefix}_ext_Z_Re",
                y=f"{prefix}_ext_Z_Im",
                c="g",
                ax=ax[1],
                logy=True,
                logx=True,
            )
            extended_data.plot(
                x=f"{prefix}_ext_Y_Re", y=f"{prefix}_ext_Y_Im", c="g", ax=ax[0]
            )
        #            extended_data.plot(x=f'{prefix}_ext_Z_Re', y=f'{prefix}_ext_Z_Im',c='g',xlim=(0,500),ylim=(0,500))
        #
        prefit_data.plot(x="Frequency(Hz)", y=f"{prefix}_err_Re", c="g", ax=ax[2])
        prefit_data.plot(x="Frequency(Hz)", y=f"{prefix}_err_Im", c="k", ax=ax[2])

        box1 = TextArea(
            lmfitting.fit_report(min_correl=0.45), textprops=dict(color="k")
        )

        box = HPacker(children=[box1], align="center", pad=0, sep=5)

        anchored_box = AnchoredOffsetbox(
            loc="lower left",
            child=box,
            pad=0.0,
            frameon=True,
            bbox_to_anchor=(1.1, 0.02),
            bbox_transform=ax[2].transAxes,
            borderpad=0.0,
        )
        ax[0].add_artist(anchored_box)
        #        axbox = plt.axes([1.1, 0.05, 0.8, 0.075])
        #        text_box = TextBox(axbox, 'Evaluate', initial=initial_text)
        #        text_box.on_submit(submit)
        #        ax[0].text(print(lmfitting.fit_report()))
        #        ax22 = ax[2].twinx()
        #        prefit_data.plot(x='Frequency(Hz)', y=f'{prefix}_res_Re',c='g', ax= ax[2])
        #        prefit_data.plot(x='Frequency(Hz)', y=f'{prefix}_res_Im',c='k', ax= ax[2])
        #                f'{prefix}_norm'
        #                prefit_data.plot(x='Frequency(Hz)', y='DATA_weightsmod_Z' ,c='orange', ax= ax22)
        ax[2].set_xscale("log")

        ax[0].set_xlim(0, prefit_data[f"{prefix}_Y_Re"].max() * 1.5)
        ax[0].set_ylim(0, prefit_data[f"{prefix}_Y_Im"].max() * 1.5)

        ax[1].set_xlim(0, prefit_data[f"{prefix}_Z_Im"].max() * 4)
        ax[1].set_ylim(0, prefit_data[f"{prefix}_Z_Im"].max() * 2)

        #                ax[2].set_yscale('log')
        plt.show()
        plt.close()
    if check_err:
        #        test_out =''
        #        test_out = [False,False]
        #        n_std = 1E-3
        #        while all(test_out) == False:
        #            for name,freqlim in [('low freq',20),('high freq',500)]:
        lf_data = prefit_data.loc[prefit_data["Frequency(Hz)"] < 20]
        hf_data = prefit_data.loc[prefit_data["Frequency(Hz)"] > 500]
        #                    , prefit_data.loc[prefit_data['Frequency(Hz)'] > freqlim]
        lf_errRe_mean, lf_errIm_mean = sum(lf_data[f"{prefix}_err_Re"] ** 2), sum(
            lf_data[f"{prefix}_err_Im"] ** 2
        )
        hf_errRe_mean, hf_errIm_mean = sum(hf_data[f"{prefix}_err_Re"] ** 2), sum(
            hf_data[f"{prefix}_err_Im"] ** 2
        )

        lf_ratio_Re, hf_ratio_Re = (
            lf_errRe_mean / pp_errRE_mean,
            hf_errRe_mean / pp_errRE_mean,
        )
        lf_ratio_Im, hf_ratio_Im = (
            lf_errIm_mean / pp_errIM_mean,
            hf_errIm_mean / pp_errIM_mean,
        )
        #            if all([lf_errRe_mean > n_std*pp_errRE_std + pp_errRE_mean, lf_errIm_mean > n_std*pp_errIM_std + pp_errIM_mean]):
        #                if test_out[0] == False:
        #                    test_out[0] = n_std
        #            if all([hf_errRe_mean > n_std*pp_errRE_std + pp_errRE_mean, hf_errIm_mean > n_std*pp_errIM_std + pp_errIM_mean]):
        #                if test_out[1] == False:
        #                    test_out[1] = n_std
        ##                =test_out + f'bad fit {name} (lim {freqlim} Hz)\n'
        #            n_std += 0.01
        test_out = [lf_errRe_mean + lf_errIm_mean, hf_errRe_mean + hf_errIm_mean]
        good_fit_test = True if any(i < 0.5 for i in test_out) else False
        return good_fit_test, test_out, MSE


def residual(params, Z_KKv, ang, model_set, weights=None):
    model = model_set.eval(params, ang=ang)
    MSE_re = (model.real - Z_KKv.real) ** 2
    MSE_im = (model.imag - Z_KKv.imag) ** 2
    MSE = MSE_re + MSE_im
    resid = model - Z_KKv
    return resid.view(np.float)


class Fit_Spectra_per_file:
    """This class will take the fit_run_arg and
    run the steps for fitting the EIS spectrum"""

    def __init__(self, _eis_run):
        pass

    #        self. =
    def __getattr__(self, attr):
        return getattr(self.eis_run, attr)

    def fit_mean_PAR_file(self):
        for PF in itertools.groupby(self.fit_run_args, lambda x: x.PAR_file):
            yield self._get_mean_EIS_from_args(*PF)

    def _get_mean_EIS_from_args(self, PF, PF_run_args_gen):
        global EvRHE

        _new_PF_mean = PF.with_name(PF.stem + "_Zmean" + PF.suffix)

        _prep_mean_args = list(PF_run_args_gen)
        _PF_mean_ovv = _prep_mean_args[0].ovv
        #    _PF_mean_ovv['Measured_OCP'] =  [i[0] for i in _PF_mean_ovv['_act0_Measured Open Circuit'].str.split()]
        #    _PF_mean_ovv['PAR_file'] = _new_PF_mean
        _PF_mean_ovv = _PF_mean_ovv.assign(**{"PAR_file": _new_PF_mean})
        _PF_data = pd.concat(i.data for i in _prep_mean_args)
        _numcols = [
            i
            for i in _PF_data.columns
            if _PF_data[i].dtype != "O" and not "Frequency" in i
        ]

        _PF_data_mean = _PF_data.groupby("Frequency(Hz)")[_numcols].mean().reset_index()
        #    _PF_data_mean.plot(x='Z Real',y='Z Imag', kind='scatter') # test plot

        _PF_data_mean = _PF_data_mean.assign(**_PF_mean_ovv.iloc[0].to_dict())
        # _PF_data_mean = _PF_data_mean.sort_values('Frequency(Hz)',ascending=False)

        # _join_cols = list(_PF_data_mean.columns.intersection(_PF_mean_ovv.columns))
        #    _PF_data_mean = _PF_data_mean.join(_PF_mean_ovv.set_index('PAR_file'),on=_join_cols,how='left')
        #    _merge = pd.concat([_PF_data_mean, _PF_mean_ovv],axis=1)
        #                       .set_index('PAR_file'),on=_join_cols,how='outer')
        _PF_data_mean[["Segment #", EvRHE, "RPM_DAC"]]
        _PF_data_mean_grp = _PF_data_mean.groupby(
            ["PAR_file", "Segment #", EvRHE, "RPM_DAC"]
        )
        fit_run_arg_mean = [
            Fit_Spectrum(
                Path(PF),
                int(Seg),
                np.round(float(E_V), 3),
                int(RPM_DAC),
                gr,
                _PF_mean_ovv,
            )
            for (PF, Seg, E_V, RPM_DAC), gr in _PF_data_mean_grp
        ][0]

        self.fit_run_arg_mean = fit_run_arg_mean


class Fit_Spectra_Collection:
    """This class will take the EIS_spectra_collection and
    run the steps for fitting the EIS spectrum"""

    global EvRHE

    def __init__(self, _EIS_spectra_pf, **kwargs):
        assert type(_EIS_spectra_pf).__name__ == "EIS_spectra_collection"
        #        print(isinstance(_EIS_spectrum_arg,EIS_Spectrum))

        self._spectra = _EIS_spectra_pf
        self._kwargs = kwargs
        _logger.warning(
            f"Starting {dt.datetime.now():%Y-%m-%d %H:%M:%S}; {len(self._spectra.spectra)} {self._spectra} {self.__class__.__name__}"
        )
        self._dest_pars = self._spectra.ovv.EIS_dest_Pars.iloc[0].with_suffix(".pkl")
        self.Model_Collection = Model_Collection()
        # _startswith='Q_'
        self.lmfit_models = self.Model_Collection

        try:
            self.fit_mean_spec = Fit_Spectrum(
                self._spectra.mean_spectrum,
                **{**self._kwargs, **dict(linKK_trimming_factor=4, res_max=0.16)},
            )
            self.load_pars_models_set_pretrials()
            self.fit_mean_spec.Model_Collection = self.Model_Collection

            if "lmfit" in self._spectra.EIS_kwargs.get("run_fit_mean"):
                self.lmfit_mean = LMfit_method(self.fit_mean_spec, run_prefit=True)

            self.results = {}

            self.load_pars_models_set_pretrials()
            for _spectrum in self._spectra.spectra:
                self._results = {}
                try:
                    _fit_spectrum = Fit_Spectrum(_spectrum, **self._spectra.EIS_kwargs)
                    self._results.update({"fit": _fit_spectrum})
                    if hasattr(
                        self, "lmfit_mean"
                    ) and not "mean" in self._spectra.EIS_kwargs.get("run_fit_mean"):
                        _fit_spectrum.Model_Collection = self.Model_Collection
                        _fit_spectrum.PRETRIALS_weights = (
                            self.lmfit_mean.PRETRIALS_weights
                        )
                        _lmfit_spectrum = LMfit_method(_fit_spectrum)
                        self._results.update({"lmfit": _lmfit_spectrum})
                    self.results.update({str(_spectrum): self._results})

                except Exception as e:
                    _logger.error(
                        f"Errror in trying fit: {_spectrum} {self.__class__.__name__}, {e} "
                    )
                    _fit_spectrum = f"fail: {e}"

            self.save_pars()
        except Exception as e:
            _logger.error(
                f"Errror in lmfit: {len(self._spectra.spectra)} {_EIS_spectra_pf} {self.__class__.__name__}, {e} "
            )

    #        for _spectrum in self._spectra.spectra:
    #            self.lmfit_results.update({str(_spectrum) : LMfit_method(self.results.get(str(_spectrum)),
    #                                           _extra_init_params = self.lmfit_mean.PRETRIALS_weights)})
    def load_pars_models_set_pretrials(self):

        if self._dest_pars.is_file():
            try:

                if not hasattr(self, "loaded_pars"):
                    self.loaded_pars = pd.DataFrame()
                    try:
                        self.loaded_pars = pd.read_pickle(self._dest_pars)
                    except AttributeError:
                        self._dest_pars.unlink()
                        _logger.warning(
                            f"AttributeError Error load pars models {self._spectra}, removed file: {self._dest_pars}"
                        )

                if hasattr(self, "lmfit_mean") and not self.loaded_pars.empty:
                    _load_pars_cols = self.loaded_pars[
                        self.lmfit_mean.PRETRIALS_weights.columns.intersection(
                            self.loaded_pars.columns
                        )
                    ]
                    self.load_mean = pd.concat(
                        [_load_pars_cols, self.lmfit_mean.PRETRIALS_weights], axis=0
                    ).reset_index()
                else:
                    self.load_mean = self.loaded_pars

                if not self.load_mean.empty:
                    self.load_mean = self.load_mean.query("Rs > 1")
                    self._best_lpars = self.load_mean.loc[
                        self.load_mean.groupby("Model_EEC")
                        .lmfit_MSE.transform("idxmin")
                        .unique()
                    ]
                    self.fit_mean_spec.PRETRIALS_weights = self._best_lpars
                # for _spectrum in self._spectra.spectra:
                #     _spectrum.PRETRIALS_weights = self._best_lpars
            except Exception as e:
                _logger.warning(
                    f"{e.__class__} Error load pars models {self._spectra}, {e}"
                )
        else:
            pass
            # for _spectrum in self._spectra.spectra:
            #     _spectrum.PRETRIALS_weights = self.lmfit_mean.PRETRIALS_weights

    def save_pars(self):
        if self.results:
            if all(["lmfit" in i.keys() for i in self.results.values()]):
                _pars = pd.concat(
                    [
                        val["lmfit"].pars_models.assign(**{"spectrum_name": key})
                        for key, val in self.results.items()
                    ],
                    sort=False,
                    ignore_index=True,
                )
                self._pars = _pars
            if not _pars.empty:
                self._pars.to_pickle(self._dest_pars)

    def check_results(self):
        """performs a re-run with best fiting as input for initial params"""
        _best_pars = self._pars.iloc[
            self._pars.groupby("Model_EEC").lmfit_aic.transform("idxmin").unique()
        ]

        for _spec, _specres in self.results.items():
            _specres["fit"].PRETRIALS_weights = _best_pars
            _lmfit_spectrum_recheck = LMfit_method(_specres["fit"])
            self.results[_spec].update(
                {
                    "recheck_fit": _specres["fit"],
                    "recheck_lmfit": _lmfit_spectrum_recheck,
                }
            )
        self.save_pars()

    def run_lmfit_over_spectra(self):
        pass

    #        self.fit_mean.linKK(linKK_trimming_factor = 4, res_max = 0.2)

    def fit_mean(self):
        pass


class Fit_Spectrum:
    """This class will take the fit_run_arg and
    run the steps for fitting the EIS spectrum"""

    global EvRHE

    def __init__(self, _EIS_spectrum_arg, **kwargs):
        assert type(_EIS_spectrum_arg).__name__ in ["EIS_Spectrum", "Fit_Spectrum"]
        self._spectrum = _EIS_spectrum_arg
        self._spectrum.EIS_kwargs.update(**kwargs)

        self.prepare_meta_data()
        self.linKK(**kwargs)
        self.Rs_guess()
        self.Rion_guess_from_lowfreq()
        self.Warburg_type_fitting()
        self.plot_linkK()
        self._add_weights_col()
        if hasattr(_EIS_spectrum_arg, "PRETRIALS_weights"):
            self.PRETRIALS_weights = _EIS_spectrum_arg.PRETRIALS_weights

    def __repr__(self):
        return f"Fit_Spectrum: {self._spectrum.__repr__()}"

    def prepare_meta_data(self):
        #        self.dataself.ovv
        PF_EIS_dest_dir = self._spectrum.data.EIS_dest_dir.unique()[0]
        #    it_run_arg.ovv.Dest_dir.iloc[0]
        # FIXME preparing meta dicts
        EISgr_data_meta = (
            self._spectrum.data[
                [
                    i
                    for i in self._spectrum.data.columns
                    if self._spectrum.data[i].nunique() == 1
                ]
            ]
            .iloc[0]
            .to_dict()
        )
        gr_EIS_ovv_meta = self._spectrum.ovv.iloc[0].to_dict()
        EISgr_meta_combined = {**gr_EIS_ovv_meta, **EISgr_data_meta}
        overlap_keys = [
            i for i in EISgr_data_meta.keys() if i in gr_EIS_ovv_meta.keys()
        ]
        matching_keys = [
            i for i in overlap_keys if EISgr_data_meta[i] == gr_EIS_ovv_meta[i]
        ]
        unmatching_keys = [i for i in overlap_keys if i not in matching_keys]
        self.spectrum_meta_info = EISgr_meta_combined
        EIS_outPath = PF_EIS_dest_dir.joinpath(
            Path(
                str(self._spectrum.PAR_file.stem.replace(".", "-"))
                + f"_{self._spectrum.E_dc_RHE_mV:.0f}mV_{self._spectrum.RPM_DAC:.0f}rpm_{self._spectrum.Segment}"
            )
        )
        self.EIS_outPath = EIS_outPath
        #    EIS_outPath_OLD = EIS_kwargs.get('EIS_dest_dir').joinpath(Path(str(Path(EISgr_data_meta['PAR_file']).stem)+'_{:.0f}mV'.format(fit_run_arg.E_dc_RHE_mV)))
        EIS_outPath.parent.mkdir(parents=True, exist_ok=True)

        _msg = f"{str(self._spectrum.PAR_file.name):s} seg({self._spectrum.Segment}) at {self._spectrum.E_dc_RHE:G} ({self._spectrum.RPM_DAC} rpm)"
        self._errmsg = _msg

    @staticmethod
    def _path_and_subdir(_EIS_path, _dir_type: str):
        _path = _EIS_path.parent.joinpath(_dir_type).joinpath(
            _dir_type + "_" + _EIS_path.name
        )
        _path.parent.mkdir(parents=True, exist_ok=True)
        return _path

    def linKK(self, **linKK_kwargs):
        #        restype_set = 'Z',  res_scale = 1, res_ylim = 0.25, res_max = 0.1)
        EIS_outpath_linKK_target = FileOperations.CompareHashDFexport(
            pd.DataFrame(), Path(str(self.EIS_outPath) + "_linkK")
        ).with_suffix(".png")
        self.EIS_outpath_linKK_target = EIS_outpath_linKK_target
        self._spectrum.EIS_kwargs.update(**linKK_kwargs)
        #        linKK_kwargs = dict(restype_set = restype_set,  res_scale = res_scale, res_ylim = res_ylim, res_max = res_max)
        try:
            if (
                self._spectrum.EIS_data_freqlim.loc[
                    self._spectrum.EIS_data_freqlim["Frequency(Hz)"].idxmin(),
                    "DATA_Zre",
                ]
                > 10
            ):
                self._spectrum.EIS_kwargs.update(dict(add_cap=True, fit_type="complex"))
            # FIXME EIS_data_freqlim
            EIS_data_KKvalid, linkKK_invalid_prefit, linKK_pars = get_KKvalid(
                self._spectrum.EIS_data, **self._spectrum.EIS_kwargs
            )
            self.EIS_data_KKvalid, self.linkKK_invalid_prefit, self.linKK_pars = (
                EIS_data_KKvalid,
                linkKK_invalid_prefit,
                linKK_pars,
            )
        except Exception as e:
            self.EIS_data_KKvalid = self.EIS_data_freqlim
            _logger.error(f"EIS fit_EEC error in validation, {self._errmsg}.\n {e}")
        self.freq_KKv = self.EIS_data_KKvalid["Frequency(Hz)"].to_numpy()
        self.ang_KKv = self.EIS_data_KKvalid["Angular"].to_numpy()
        self.Z_KKv = self.EIS_data_KKvalid["DATA_Z"].to_numpy()
        #        self.EIS_data = self.EIS_data.loc[~self.EIS_data['Frequency(Hz)'].isin(EIS_data_KKvalid['Frequency(Hz)']),'Valid'] = False
        #        linKK_pars.update({'linKK_invalid_freqs' : ', '.join(EIS_data_KKvalid.loc[~EIS_data_KKvalid['Frequency(Hz)'].isin(EIS_data_KKvalid['Frequency(Hz)']),'Frequency(Hz)'].astype(str).to_list())})

        if self._spectrum.EIS_kwargs.get("export_raw_testing", False):
            testing_path = Path.cwd().joinpath(
                "testing_data", self.EIS_outPath.name + "_spectrum_raw"
            )
            self.EIS_data_KKvalid.to_excel(testing_path.with_suffix(".xlsx"))

    def plot_linkK(self):
        try:
            plot_linKK(
                self.EIS_data_KKvalid,
                self._spectrum.EIS_data,
                save_target=self.EIS_outpath_linKK_target,
                plot_prefit=True,
                linkKK_invalid_prefit=self.linkKK_invalid_prefit,
                **{**self.linKK_pars, **self.spectrum_meta_info},
            )
        except Exception as e:
            _logger.warning(
                f"EIS in plotting linKK, {str(self._spectrum.PAR_file):s} at {self._spectrum.E_dc_RHE:G} ({self._spectrum.RPM_DAC} rpm).\n {e}"
            )

    def _add_weights_col(self):
        lin_angle = stats.linregress(
            np.log(self.freq_KKv),
            np.angle(self.EIS_data_KKvalid.linKK_Z.values, deg=True),
        )
        lin_subt_Zangle = np.angle(self.EIS_data_KKvalid.linKK_Z.values, deg=True) - (
            lin_angle.slope * np.log(self.freq_KKv) + lin_angle.intercept
        )
        lmfit_weights_linZang = (lin_subt_Zangle) ** 2 / np.mean(
            (lin_subt_Zangle) ** 2
        ) + self.EIS_data_KKvalid.lmfit_weights_unit
        self.EIS_data_KKvalid = self.EIS_data_KKvalid.assign(
            **{"lmfit_weights_linZangle": lmfit_weights_linZang}
        )

    def Rs_guess(self, Rs_limits=(0, 200)):
        Rs_guess_data = self.EIS_data_KKvalid.loc[
            np.abs(self.EIS_data_KKvalid.DATA_Zmod).idxmin(), "DATA_Zre"
        ]
        Rs_guess_linKK = self.EIS_data_KKvalid.linKK_Z.to_numpy().real.min()
        Rs_guess = np.mean(
            [
                i
                for i in [Rs_guess_data, Rs_guess_linKK]
                if Rs_limits[0] < i < Rs_limits[1]
            ]
        )
        if not Rs_guess:
            Rs_guess = Rs_guess_data
        self.Rs_guess = Rs_guess

    def Rion_guess_from_lowfreq(self):
        _lintangent = {}
        _lint_spec = {}
        try:

            _x, _y = self.EIS_data_KKvalid.DATA_Zre, self.EIS_data_KKvalid["DATA_-Zim"]
            _xmax = _x.max()
            xrange = np.arange(0, _xmax, _xmax / len(_y))
            if len(xrange) > len(_x):
                xrange = xrange[0 : len(_x)]

            for curve in ["hf", "lf", "uf"][::-1]:
                if curve == "uf":
                    _freq = self.EIS_data_KKvalid.loc[
                        (self.EIS_data_KKvalid["Frequency(Hz)"] < 2)
                    ]
                    _range = (0.85, 1.2, 2, 3)

                if curve == "lf":
                    _freq = self.EIS_data_KKvalid.loc[
                        (self.EIS_data_KKvalid["Frequency(Hz)"] < 1e3)
                        & (self.EIS_data_KKvalid["Frequency(Hz)"] > 3)
                    ]
                    _range = (0.85, 1.2, 2, 4)
                if curve == "hf":
                    _freq = self.EIS_data_KKvalid.loc[
                        self.EIS_data_KKvalid["Frequency(Hz)"] > 1e3
                    ]
                    _range = (0.45, 1.57, 2, 2)

                _xc, _yc = _freq.DATA_Zre, _freq["DATA_-Zim"]
                x0 = 0
                _dydx = _yc.diff() / _xc.diff()
                _dydx = _dydx.rolling(_range[3]).mean()
                _tang = _dydx[(_range[0] < _dydx) & (_dydx < _range[1])]
                _tangclean = []
                if len(_tang) > 2:
                    _tangclean = [
                        n for n, i in enumerate(np.diff(_tang.index)) if np.abs(i) < 3
                    ]
                if _tangclean and curve == "hf":
                    _tang = _tang.iloc[_tangclean]
                _idx = []
                if len(_tang) >= _range[2] and np.all(np.diff(_tang.index) < 3):
                    _idx = list(_tang.index)
                elif len(_tang) <= _range[2] and len(_tang) > 0:
                    for n in _tang.index:
                        if _x.min() < _xc[n] and _xc[n] < _x.min() * 1.3:
                            _idx.append(n)

                if _idx:
                    x0 = _xc[_tang.index].mean() - _yc[_tang.index].mean()
                    _dataslice = self.EIS_data_KKvalid.loc[_tang.index]
                    _dtfreq = _dataslice["Frequency(Hz)"]
                    _lintangx = np.array([x + x0 for x in xrange])
                    _lintangy = _lintangx - _lintangx.min()
                    _lint_spec.update(
                        {
                            f"_lintangent_{curve}_re": _lintangx,
                            f"_lintangent_{curve}_im": _lintangy,
                        }
                    )
                    _fmin, _fmax, _fmean = _dtfreq.min(), _dtfreq.max(), _dtfreq.mean()
                    _lintangent.update(
                        {
                            f"_lintangent_{curve}": {
                                "x0": x0,
                                "freq_min": _fmin,
                                "len_n": len(_tang),
                                "freq_max": _fmax,
                                "freq_mean": _fmean,
                            }
                        }
                    )
                    self.EIS_data_KKvalid = self.EIS_data_KKvalid.assign(**_lint_spec)

            self._lintangent_pars = _lintangent
            _flatpars = {
                f"{k1}_{k0}": v0
                for k1, val1 in _lintangent.items()
                for k0, v0 in val1.items()
            }
            self.spectrum_meta_info.update(**_flatpars)

            def testplot():
                fig, ax = plt.subplots()
                ax.scatter(_xc, _yc, label="data")
                ax.plot(_x, _x - _x.min(), label="x")
                for freq in set([i[0:-3] for i in _lint_spec.keys()]):
                    ax.plot(
                        _lint_spec[f"{freq}_re"],
                        _lint_spec[f"{freq}_im"],
                        ls="--",
                        label=freq,
                    )
                ax.legend()

            # testplot()

        except Exception as e:
            _logger.warning(
                f"Rion_guess_from_lowfreq plotting failed for: {self._spectrum}, because {e}"
            )

        _Rion_intrcept = _lintangent.get("_lintangent_lf", {"x0": 0}).get("x0", 0)
        Z_Rion_diff = _Rion_intrcept - self.Rs_guess
        if Z_Rion_diff > 0:
            self.R_ion_guess = Z_Rion_diff * 3
            self.spectrum_meta_info.update(**{"R_ion_guess": self.R_ion_guess})

    def Warburg_type_fitting(self, _lin_window_size=15):
        _lin = {}
        lin_slopes = [(0.25, "lightgreen"), (0.5, "grey"), (1, "orange")]
        #        for yax in ['Zre','-Zim']:
        #        _lin.update({yax : linregress(spec.query('ang_Warburg > 0.3').ang_Warburg,spec.query('ang_Warburg > 0.3')['DATA_'+yax]),
        #                     yax+'_low' : linregress(spec.query('ang_Warburg < 0.045').ang_Warburg,spec.query('ang_Warburg < 0.045')['DATA_'+yax])})
        #        spec['W_lin_'+yax] = _lin[yax].slope * spec.ang_Warburg + _lin[yax].intercept
        #        spec['W_lin_'+yax+'_low'] = _lin[yax+'_low'].slope * spec.ang_Warburg + _lin[yax+'_low'].intercept
        try:
            _anglim = 30
            _lin_WB_freq_loq = self.EIS_data_KKvalid.query("Angular < 30")
            while (
                len(_lin_WB_freq_loq) < 3
                and _anglim < 0.9 * self.EIS_data_KKvalid.Angular.max()
            ):
                _anglim += 10
                _lin_WB_freq_loq = self.EIS_data_KKvalid.query("Angular < @_anglim")

            _angWhighlim = 0.3
            _lin_WB_angW_high = self.EIS_data_KKvalid.query("Ang_Warburg > 0.3")
            while (
                len(_lin_WB_angW_high) < 3
                and _angWhighlim > 1.5 * self.EIS_data_KKvalid.Ang_Warburg.min()
            ):
                _angWhighlim -= 0.05
                _lin_WB_angW_high = self.EIS_data_KKvalid.query(
                    "Ang_Warburg > @_angWhighlim"
                )

            _lin_WB_angW_low = self.EIS_data_KKvalid.query("Ang_Warburg < 0.045")
            _yax_cols = {}
            for yax in ["Zre", "-Zim"]:
                _DATA_Z_col = f"DATA_{yax}"
                for _range in ["high", "low"]:
                    _rvar = eval(f"_lin_WB_angW_{_range}")
                    if len(_rvar) > 2:
                        _lrgess = linregress(_rvar.Ang_Warburg, _rvar[_DATA_Z_col])
                        _lr_res = dict(zip(_lrgess._fields, _lrgess))
                        _lin.update({f"{yax}_angW_lin_{_range}": _lr_res})
                        _yax_cols.update(
                            {
                                f"{yax}_lin_WB_angW_{_range}": _lrgess.slope
                                * self.EIS_data_KKvalid.Ang_Warburg
                                + _lrgess.intercept
                            }
                        )

                if yax == "Zre":
                    if len(_lin_WB_freq_loq) > 3:
                        _lrgress = linregress(
                            _lin_WB_freq_loq["DATA_Zre"], _lin_WB_freq_loq["DATA_-Zim"]
                        )
                        _lr_ress = dict(zip(_lrgress._fields, _lrgress))
                        _lin.update({f"{yax}_lin_slopes": _lr_ress})
                        _yax_cols.update(
                            {
                                f"{yax}_lin_slopes": _lin[f"{yax}_lin_slopes"]["slope"]
                                * self.EIS_data_KKvalid[_DATA_Z_col]
                                + _lin[f"{yax}_lin_slopes"]["intercept"]
                            }
                        )
            for _range in ["high", "low"]:  # add rel. difference of slopes between axes
                _slRe = _lin.get(f"Zre_angW_lin_{_range}", {}).get("slope", 0)
                _slIm = _lin.get(f"-Zim_angW_lin_{_range}", {}).get("slope", 0)
                if _slIm != _slRe != 0:
                    _sl_d = _slRe - _slIm
                    _sl_rel = 100 * _sl_d / _slRe
                    _lin.update(
                        {
                            f"ZreZim_angW_lin_{_range}": {
                                "slopes_diff_abs": _sl_d,
                                "slope_diff_rel": _sl_rel,
                            }
                        }
                    )
                    # f'angW_lin_ZreZim_diff_slope_{_range}' : _sl_d})
                #                self.EIS_data_KKvalid['W_lin_'+yax] = _lin[yax].slope * self.EIS_data_KKvalid.ang_Warburg + _lin[yax].intercept
            #            self.EIS_data_KKvalid['W_lin_'+yax+'_low'] = _lin[yax+'_low'].slope * self.EIS_data_KKvalid.ang_Warburg + _lin[yax+'_low'].intercept
            for _slope, _ in lin_slopes:
                perr_set = 1000
                #            for _win_size in [7,10,15,25]:
                #            for win in spec.rolling(_lin_window_size):
                best = ()
                for i in range((len(self.EIS_data_KKvalid) - _lin_window_size)):
                    win = self.EIS_data_KKvalid.iloc[i : i + _lin_window_size]
                    popt, pcov = curve_fit(
                        func_lin(_slope), win.DATA_Zre, win["DATA_-Zim"]
                    )
                    perr = float(np.sqrt(np.diag(pcov)))
                    #                print(win.index,popt,pcov,perr)
                    if perr < perr_set:
                        perr_set = perr
                        best = (_slope, win, popt, perr)
                if best:
                    _yax_cols.update(
                        {
                            f"WB_lin_-Zim_a{_slope}": func_lin(best[0])(
                                self.EIS_data_KKvalid.DATA_Zre, best[2][0]
                            )
                        }
                    )
                    _lin.update(
                        {
                            f"lin_slope_{_slope}": {
                                "popt": best[2][0],
                                "win_size": len(best[1]),
                                "perr": best[-1],
                            }
                        }
                    )

            self.EIS_data_KKvalid = self.EIS_data_KKvalid.assign(**_yax_cols)
        except Exception as e:
            _logger.error(f"Warburg fitting failed for: {self._spectrum}, because {e}")

        try:

            _WB_lin_pars = pd.DataFrame(
                dict(
                    [
                        (f"WB_{key}_{k2}", v2)
                        for key, val in _lin.items()
                        for k2, v2 in val.items()
                    ]
                ),
                index=self._spectrum.ovv.index,
            )
            _WB_lin_pars = _WB_lin_pars.assign(**self.spectrum_meta_info)
            # pd.concat([self.ovv,_WB_lin_pars],axis=1)
            self.WB_lin_pars = _WB_lin_pars
            plot_lin_Warburg(
                self.EIS_data_KKvalid,
                _lin,
                lin_slopes,
                pars=_WB_lin_pars,
                dest_path=self._path_and_subdir(self.EIS_outPath, "lin_Warburg"),
            )

        except Exception as e:
            _logger.warning(
                f"Warburg plotting failed for: {self._spectrum}, because {e}"
            )

    def DRT_fitting(self):
        if self._spectrum.EIS_fit_kwargs.get("DP_DRT_fit", False):
            try:
                DP_DRT_destpath = self._path_and_subdir(self.EIS_outPath, "DP_DRT")
                #                self.EIS_outPath.parent.joinpath('DP_DRT').joinpath('DP_DRT_'+self.EIS_outPath.name)
                #                DP_DRT_destpath.parent.mkdir(parents=True,exist_ok=True)
                #                DP_DRT_destdir = self.EIS_outPath.parent.joinpath('DP_DRT')
                DP_DRT_out = DP_DRT_analysis(
                    self.EIS_data_KKvalid.DATA_Z.to_numpy()
                )  # TODO FIX EXPORTING
                self.spectrum_meta_info.update(
                    {"DP_DRT_fit": DP_DRT_destpath, "DP_DRT_run_success": True}
                )
            except Exception as e:
                self.spectrum_meta_info.update(
                    {
                        "DP_DRT_fit": DP_DRT_destpath,
                        "DP_DRT_run_success": False,
                        "DP_DRT_run_msg": e,
                    }
                )
                _logger.error(
                    f"EIS fit_EEC error in GP_DRT_fit, {self._errmsg}).\n {e}"
                )
        #        EIS_outPath.name+'_DP_DRT_'
        if self._spectrum.EIS_kwargs.get("GP_DRT_fit", False):
            try:
                GP_DRT_destpath = self._path_and_subdir(self.EIS_outPath, "GP_DRT")
                #                self.EIS_outPath.parent.joinpath('GP_DRT').joinpath('GP_DRT_'+self.EIS_outPath.name)
                #                GP_DRT_destpath.parent.mkdir(parents=True,exist_ok=True)
                DRT_data_KKvalid = self._spectrum.data.sort_values(
                    by="Frequency(Hz)", ascending=True
                )
                #            prep_GP_DRT_raw_data(EIS_data_raw)
                _Z_exp, _Z_star, res_fit_params, res_fit_arrays = run_GP_DRT_fit(
                    DRT_data_KKvalid, **{"GP_DRT_savefig_path": GP_DRT_destpath}
                )
                self.spectrum_meta_info.update(
                    {"GP_DRT_fit": GP_DRT_destpath, "GP_DRT_run_success": True}
                )
            except Exception as e:
                self.spectrum_meta_info.update(
                    {
                        "GP_DRT_fit": DP_DRT_destpath,
                        "GP_DRT_run_success": False,
                        "GP_DRT_run_msg": e,
                    }
                )
                _logger.error(f"EIS fit_EEC error in GP_DRT_fit, {self._errmsg}.\n {e}")

        # ====!!!==== #


# class Fit_Model_spectrum():
#     # TODO install model fitting
#     def __init__(self):
#         pass
# self.validate_params()
# if self._missing_pars:
#     self.set_standard_init_from_guesses()
# self.validate_params()
# if not self._missing_pars:


class LMfit_method:
    """Take a Fit_Spectrum class and performs the lmfit on models from Model_Collection.
    Performs prefit when PRETRIALS_weights is not given.
    returns
    """

    _standard_Rs_guess = 20

    def __init__(self, _obj, **kwargs):
        self._obj = _obj
        self.kwargs = kwargs

        if not hasattr(self._obj, "Model_Collection"):
            self.Model_Collection = Model_Collection()
        else:
            self.Model_Collection = self._obj.Model_Collection
        self.standard_init_params = self.Model_Collection.standard_init_params
        if hasattr(self._obj, "PRETRIALS_weights"):
            self.PRETRIALS_weights = self._obj.PRETRIALS_weights

        if type(self._obj).__name__ == "Fit_Spectrum":
            self.spectrum = self._obj
            self.Rs_guess_set()
            self._init_params()
            self.set_weight_opts()
            try:
                if not hasattr(self, "PRETRIALS_weights") or self.kwargs.get(
                    "run_prefit", False
                ):
                    self.run_prefit_models()
                #                self.PRETRIALS_weights = self.lmfit_mean.PRETRIALS_weights
                self.run_final_fit_models()
                self.final_export()
            except Exception as e:
                _logger.error(
                    f"LMFit failed for {self._obj}, because:" + "\n " + str(e)
                )

    def Rs_guess_set(self):
        if hasattr(self._obj, "Rs_guess"):
            self.Rs_guess = self._obj.Rs_guess
        else:
            self.Rs_guess = self._standard_Rs_guess

    def set_weight_opts(self):
        weight_opts = {
            key: self.spectrum.EIS_data_KKvalid[key]
            for key in self.spectrum.EIS_data_KKvalid.columns
            if key.startswith("lmfit_weights")
        }
        self.weight_opts = weight_opts

    def _init_params(self):

        lower_bound, max_ub = 1e-9, 1e10
        max_C_ub = 0.5
        n_lb = 0.0

        self.standard_init_params.add(
            "Rs",
            value=self.Rs_guess,
            min=0.65 * self.Rs_guess,
            max=2 * self.Rs_guess,
            brute_step=1,
            vary=True,
        )
        if hasattr(self._obj, "R_ion_guess"):
            self.standard_init_params["R_ion"].set(value=self._obj.R_ion_guess)
        else:
            self.standard_init_params["R_ion"].set(value=self.Rs_guess * 1.7)
        self.lower_bound = lower_bound
        self.upper_bound = max_ub

    #        standard_init_params = params_extra_setting(standard_init_params, EISgr_data_EV, EIS_fit_kwargs) # skip extra settings
    #

    def run_prefit_models(self):
        _prefit_methods = ["least_squares", "nelder", "ampgo", "differential_evolution"]
        default_prefit_methods = [
            _prefit_methods[0],
            _prefit_methods[1],
            _prefit_methods[-1],
        ][0:1]
        _prefit_args = (
            self.spectrum,
            self.standard_init_params,
            self.weight_opts,
            self.Model_Collection,
        )
        if hasattr(self, "PRETRIALS_weights"):
            _extr_init = self.PRETRIALS_weights
        else:
            _extr_init = pd.DataFrame()

        (
            PRETRIALS_weights,
            _bad_models_out,
            PREFITS_weights,
            PREFITS_spectra,
        ) = prefit_test(
            *_prefit_args,
            _extra_init_params=_extr_init,
            prefit_methods=default_prefit_methods,
            pickle_spectra=False,
            **self.kwargs,
        )

        # PRETRIALS_weights,_bad_models, all_PRETRIALS,all_PRE_spectra
        self.PRETRIALS_weights = PRETRIALS_weights
        self._bad_models = _bad_models_out
        self.all_PRETRIALS = PREFITS_weights
        self.all_PRE_spectra = PREFITS_spectra

    def run_final_fit_models(self):
        _finalfit_args = (
            self.spectrum.ang_KKv,
            self.spectrum.Z_KKv,
            self.spectrum.EIS_data_KKvalid,
            self.weight_opts,
            self.PRETRIALS_weights,
            self.Model_Collection,
            self.standard_init_params,
        )

        _weights = [
            "lmfit_weights_unit",
            "lmfit_weights_mod_Z",
            "lmfit_weights_linZangle",
            "lmfit_weights_mod_Y",
            "lmfit_weights_prop",
        ]
        EIS_fit_data, pars_models = final_fit(
            *_finalfit_args,
            default_method="least_squares",
            default_weights=_weights[-2],
            use_default_wts=True,
            use_default_method=True,
        )
        self.EIS_fit_data = EIS_fit_data
        self.pars_models = pars_models

    def final_export(self):
        vers = FileOperations.EIS_version
        spectra_fit_outpath = self.spectrum.EIS_outPath.with_name(
            self.spectrum.EIS_outPath.stem + f"_spectrumfit_v{vers}"
        ).with_suffix(".xlsx")
        spectra_raw_outpath = self.spectrum.EIS_outPath.with_name(
            self.spectrum.EIS_outPath.stem + f"_spectrumraw_v{vers}"
        ).with_suffix(".xlsx")
        pars_outpath = self.spectrum.EIS_outPath.with_name(
            self.spectrum.EIS_outPath.stem + f"_pars_v{vers}"
        ).with_suffix(".xlsx")

        self.EIS_fit_data = self.EIS_fit_data.assign(
            **{"File_Pars": pars_outpath, "File_SpecRaw": spectra_raw_outpath}
        )
        self.EIS_fit_data = self.EIS_fit_data.assign(**self.spectrum.spectrum_meta_info)
        spectra_fit_outpath_target = FileOperations.CompareHashDFexport(
            self.EIS_fit_data, spectra_fit_outpath
        )

        self.spectrum._spectrum.data = self.spectrum._spectrum.data.assign(
            **{"File_Pars": pars_outpath, "File_SpecFit": spectra_fit_outpath}
        )
        spectra_raw_outpath_target = FileOperations.CompareHashDFexport(
            self.spectrum._spectrum.data, spectra_raw_outpath
        )

        self.pars_models = self.pars_models.assign(
            **{
                "File_SpecFit": spectra_fit_outpath_target,
                "File_SpecRaw": spectra_raw_outpath_target,
            }
        )
        self.pars_models = self.pars_models.assign(**self.spectrum.spectrum_meta_info)
        pars_outpath_target = FileOperations.CompareHashDFexport(
            self.pars_models, pars_outpath
        )

        if "Plot" in self.spectrum._spectrum.EIS_kwargs.get(
            "EIS_single_output", "Text, Plot"
        ):
            spectra_fit_outpath_png = spectra_fit_outpath.with_suffix(
                ".png"
            ).with_suffix(".png")
            try:
                EIS_plotting_per_EV(
                    self.EIS_fit_data,
                    self.pars_models.query("FINAL_FIT == 1"),
                    spectra_fit_outpath_png,
                    plot_show=False,
                    std_print_model="M_RC_CPE_W",
                )
            except Exception as e:
                _logger.error(f"LM fit export plotting error: {e}")


#    def __repr__(self):
#        _spec = f'LMfit: {str(self.spectrum)}'
#        = prefit_test(fit_run_arg, EIS_fit_kwargs, weight_opts, EIS_data_KKvalid,
#                                  Z_KKv, ang_KKv, standard_init_params, MTHDS,
#                                  prefit_methods= ['leastsq'], prefit_models = mod_lst_filter)
#    for n,r in PRETRIALS_weights.iterrows(): # TODO
#        make_prefit_frame(EIS_data_KKvalid, r['lmfit_best_trial'],plot = True) # TODO
#        all_PRETRIALS_grp = all_PRETRIALS.groupby('Model_EEC')
# ang_KKv,Z_KKv,EIS_data_KKvalid, weight_opts, \
# PRETRIALS_weights, lmfit_models, standard_init_param= self.spectrum.ang_KKv, self.spectrum.Z_KKv, \
#                                                        self.spectrum.EIS_data_KKvalid, self.weight_opts,\
#                                                        self.PRETRIALS_weights, self.lmfit_models, self.standard_init_params
#                                                        default_method='least_squares'


def prefit_test(
    *_prefit_args,
    _extra_init_params=pd.DataFrame(),
    prefit_methods=["least_squares"],
    pickle_spectra=False,
    **kwargs,
):
    # prefit_models=Model_Collection(),
    #%%
    _spectrum, standard_init_params, weight_opts, prefit_models = _prefit_args
    # self.spectrum, self.standard_init_params, self.weight_opts
    # prefit_methods= ['least_squares']
    # prefit_models=self.lmfit_models
    #    pars_lst, fit_data_valid_lst,outP = [], [], {}
    #    EIS_data_valid =  _spectrum.EIS_data_KKval.id
    # weight_opts = {key : _spectrum.EIS_data_KKvalid[key] for key in _spectrum.EIS_data_KKvalid.columns if key.startswith('lmfit_weights')}
    unit_weight = weight_opts["lmfit_weights_unit"]
    PREFITS_spectra = pd.DataFrame()

    wt_check, wt_spectra_lst = {}, []
    # prefit_models = [i for i in prefit_models if 'Q_TLM_W_pQadRct' in str(i.classname)] # FIXME
    for mod_instance in prefit_models:  #
        # mod_indx, model_set,rand_n,_color = prefit_models[-3] # FIXME
        mod_indx = mod_instance.classname
        model_set = mod_instance.model
        modname = model_set.name
        mod_color = mod_instance.color

        if mod_instance.parameters_guesses:
            params_model = mod_instance.parameters_guesses
            if "Rs" in params_model.keys():
                params_model["Rs"].set(value=standard_init_params["Rs"].value)
            if (
                "R_ion" in params_model.keys()
                and "R_ion" in standard_init_params.keys()
            ):
                params_model["R_ion"].set(value=standard_init_params["R_ion"].value)
                # max=_spectrum.Z_KKv.real.max()*3)
        else:
            params_model = model_set.make_params()

            for pn in params_model:
                params_model[pn].set(
                    value=standard_init_params[pn].value,
                    min=standard_init_params[pn].min,
                    max=standard_init_params[pn].max,
                )
        par_options = {}
        par_options.update({"standard": params_model})

        try:
            _get_params = pd.DataFrame()
            recheck_msg = "not used"
            #            _get_params, recheck_msg = fitting_recheck_params(fit_run_arg,modname,params_model, **EIS_fit_kwargs) # TODO
            _extr_params_mod = pd.DataFrame()
            if not _extra_init_params.empty and hasattr(
                _extra_init_params, "Model_EEC"
            ):
                #            _logger.warning(recheck_msg)
                _extr_params_mod = _extra_init_params.loc[
                    _extra_init_params["Model_EEC"] == modname
                ]

            if not _extr_params_mod.empty:
                pre_params = params_model.copy()
                for pn in pre_params:
                    _pn_ext = _extr_params_mod[pn].iloc[0]
                    # print('pn_ext:',_pn_ext)
                    if (
                        pre_params[pn].vary == True
                        and pd.isna(_pn_ext) == False
                        and pn != "Rs"
                    ):
                        # print('pn_ext:',_pn_ext)
                        pre_params[pn].set(
                            value=_pn_ext * random.randint(75, 110) / 100
                        )
                par_options.update({"pre_params": pre_params})
            #                _logger.info(f'Prefit used recheked params for {modname}')
            else:
                _logger.info(
                    f'Prefit _extra_init_params empty {_spectrum}\n message:"{recheck_msg}'
                )
        except Exception as e:
            _logger.error(
                f"fitting_recheck_params error: {e}" + "\n" + f"for model: {modname}"
            )

        #        std_prefit_DE = model_set.fit(Z_KKv,params_model, ang=ang_KKv, weights= wval, method= MTHDS[5])
        # TODO get_from_good_fits
        #            if 'Rs' in _good_params.keys():
        try:
            for parnm, par_obj in par_options.items():
                # parnm,par_obj = 'standard',par_options['standard']
                prefit_params = model_set.make_params()
                for pn in prefit_params:
                    prefit_params.add(
                        pn,
                        value=par_obj[pn].value,
                        min=par_obj[pn].min,
                        max=par_obj[pn].max,
                    )
                for wname, wval in weight_opts.items():
                    wname, wval
                    for method in prefit_methods:
                        method
                        _method_kwargs = {}
                        # if 'evolution' in method:
                        #     _method_kwargs = {'strategy' : 'rand1bin'}
                        # , 'workers' : -1}
                        #            make_prefit_frame(EIS_data_KKvalid, pre_prefit,plot = True, norm= wval) #TODO
                        try:

                            best_trial_weights = model_set.fit(
                                _spectrum.Z_KKv,
                                prefit_params,
                                ang=_spectrum.ang_KKv,
                                weights=wval,
                                method=method,
                            )
                            # fit_kws = _method_kwargs
                            # best_trial_weights.plot(ax_fit_kws={'xscale' : 'log'}) # for fast plotting
                            #                    best_trial_weights = model_set.fit(Z_KKv,pre_prefit.params,ang=ang_KKv, weights = unit_weight , method= pre_prefit.method)
                            if best_trial_weights.success:
                                (
                                    pp_good_low_high,
                                    best_test_msg,
                                    MSE,
                                ) = make_prefit_frame(
                                    _spectrum.EIS_data_KKvalid,
                                    best_trial_weights,
                                    prefix="pp",
                                    plot=0,
                                    norm=unit_weight,
                                )
                                #                     = make_prefit_frame(EIS_data_KKvalid,best_trial_weights, prefix = 'pp',plot = 0, norm= unit_weight)
                                if pickle_spectra:
                                    post_br_frame = make_prefit_frame(
                                        _spectrum.EIS_data_KKvalid,
                                        best_trial_weights,
                                        get_frame=1,
                                    )
                                    post_br_frame = post_br_frame.assign(
                                        **{
                                            "Model_EEC": modname,
                                            "lmfit_weights_nm": wname,
                                            "lmfit_Rs_setting": parnm,
                                            "lmfit_method": method,
                                            "mod_index": ", ".join(
                                                (modname, wname, parnm, method)
                                            ),
                                        }
                                    )
                                    wt_spectra_lst.append(post_br_frame)
                                #                    merge_cols = list(post_br_frame.columns.intersection(EIS_data_KKvalid.columns))
                                #                    EIS_data_KKvalid_BR_specs = pd.merge(EIS_data_KKvalid_BR_specs,post_br_frame,on =merge_cols,how='inner',left_index=True,right_index=True)

                                wt_check.update(
                                    {
                                        (modname, wname, parnm, method): {
                                            "Model_EEC": modname,
                                            "lmfit_weights_nm": wname,
                                            "lmfit_Rs_setting": parnm,
                                            "lmfit_method": method,
                                            "lmfit_MSE": MSE,
                                            "best_test_msg": best_test_msg,
                                            "sum_test_MSE": sum(MSE + best_test_msg),
                                            "lmfit_method": best_trial_weights.method,
                                            "lmfit_aic": best_trial_weights.aic,
                                            "lmfit_succes": best_trial_weights.success,
                                            "lmfit_best_trial": best_trial_weights,
                                            "plot_color": mod_color,
                                            "BRUTE_FIT": 1,
                                            "FINAL_FIT": 0,
                                            "PREFIT": 0,
                                            **best_trial_weights.params.valuesdict(),
                                        }
                                    }
                                )
                        except RuntimeError:
                            _logger.error(
                                f'RunTimeError for {_spectrum}\n message:"{mod_indx}'
                            )
        except Exception as e:
            _logger.error(
                f'PREfit Model fit error {modname},{parnm},{wname},{method} for {_spectrum}\n message for "{mod_indx}" :{e}'
            )
    #            wt_check.append({modname : {'weights_name' : (wname,MSE,best_test_msg,sum(MSE+best_test_msg),best_trial_weights))
    #        model_set.eval(pre_prefit.params,ang=np.logspace(-3,4,endpoint=True)*2.*pi, weights = statiscial_weights, method= pre_prefit.method)
    _wt_ch_idx = ["_".join(i for i in a) for a in wt_check.keys()]
    PREFITS_weights = pd.DataFrame(
        data=wt_check.values(), index=_wt_ch_idx
    ).sort_index()
    # .sort_values(by='lmfit_aic')
    # PRE_min = PREFITS_weights.iloc[PREFITS_weights.reset_index().groupby('level_0')['lmfit_aic'].transform('idxmin').unique()]
    if wt_spectra_lst:
        PREFITS_spectra = pd.concat(wt_spectra_lst, sort=False, ignore_index=True)
    #    [make_prefit_frame(_spectrum.EIS_data_KKvalid, r['lmfit_best_trial'],plot = True, norm= wval) for n,r in PREFITS_weights.iterrows()] #TODO
    min_idxs = (
        PREFITS_weights.groupby("Model_EEC")["lmfit_MSE"].transform("idxmin").unique()
    )
    # post_check = {}
    PREFITS_weights.loc[min_idxs, "PREFIT"] = 1
    PRETRIALS_weights = PREFITS_weights.query("PREFIT == 1")
    # loc[min_idxs]
    #%%
    #    [make_prefit_frame(EIS_data_KKvalid, r['lmfit_best_trial'],plot = True, norm= wval) for n,r in PRETRIALS_weights.iterrows()] #TODO
    #    PP_out = pd.DataFrame(data=post_check.values(),index=post_check.keys()).sort_values(by='MSE')
    #    PRETRIALS_weights = PP_out.iloc[PP_out.reset_index().groupby('level_0')['MSE'].transform('idxmin').unique()]
    _bad_models = PRETRIALS_weights.loc[
        PRETRIALS_weights.lmfit_MSE
        > PRETRIALS_weights.lmfit_MSE.mean() + PRETRIALS_weights.lmfit_MSE.std()
    ]
    _bad_models_out = [i for i in _bad_models.index]
    #    [gr['MSE'].idxmin() for n,gr in PREFITS_weights.reset_index().groupby('level_0')]
    #    PREFITS_weights['MSE'].idxmin()
    return PRETRIALS_weights, _bad_models_out, PREFITS_weights, PREFITS_spectra


def final_fit(
    *_finalfit_args,
    default_method="least_squares",
    default_weights="lmfit_weights_unit",
    use_default_wts=True,
    use_default_method=True,
    **kwargs,
):

    (
        ang_KKv,
        Z_KKv,
        EIS_data_KKvalid,
        weight_opts,
        PRETRIALS_weights,
        Model_Collection,
        standard_init_params,
    ) = _finalfit_args

    pars_lst, fit_data_valid_lst, outP = [], [], {}
    if use_default_method:
        _best_trial_method = default_method
    else:
        _best_trial_method = kwargs.get("best_trial_method", default_method)
    # lstsq_method = kwargs.get('best_trial_method' ,default_method)
    # weights_used_out = 'lmfit_weights_unit'
    pretrials_lmfit_models = [
        mi
        for mi in Model_Collection
        if mi.model.name in PRETRIALS_weights.groupby("Model_EEC").groups.keys()
    ]
    # for mod_indx, model_set, rand_n,_color in pretrials_lmfit_models: #
    for mod_instance in Model_Collection:  #
        # mod_indx, model_set,rand_n,_color = prefit_models[-3] # FIXME
        mod_indx = mod_instance.classname
        mod_inst_name = mod_instance.name
        model_set = mod_instance.model
        modname = model_set.name
        mod_color = mod_instance.color
        modname = model_set.name
        params_model = model_set.make_params()
        mod_par_guess = mod_instance.parameters_guesses
        outP = {}
        #            print('use previous')
        #            params = InitP1.query('Model_EEC')
        par_options = {}

        for pn in params_model:

            if pn != "Rs" and pn in mod_par_guess.keys():
                params_model[pn].set(
                    value=mod_par_guess[pn].value,
                    min=mod_par_guess[pn].min,
                    max=mod_par_guess[pn].max,
                )
            else:
                params_model[pn].set(
                    value=standard_init_params[pn].value,
                    min=standard_init_params[pn].min,
                    max=standard_init_params[pn].max,
                )

        # if mod_instance.parameters_guesses:
        #     params_model = mod_instance.parameters_guesses
        # else:

        #     for pn in params_model:

        par_options.update(
            {"standard": {"params": params_model, "weights": default_weights}}
        )
        if not PRETRIALS_weights.empty:
            # PRETRIALS_weights.loc[PRETRIALS_weights(modname),:]*PRETRIALS_weights.loc[modname]['lmfit_aic'].idxmin()
            pretrial_best_row = PRETRIALS_weights.query("Model_EEC == @modname")

            if not pretrial_best_row.empty:
                pretrial_best_row = pretrial_best_row.sort_values(
                    "lmfit_MSE", ascending=True
                )
                pretrial_lmfit = pretrial_best_row["lmfit_best_trial"].iloc[0]
                pretrial_weights = pretrial_best_row["lmfit_weights_nm"].iloc[0]
                if not pretrial_weights in weight_opts.keys():
                    pretrial_weights = default_weights

                pretrial_params = params_model.copy()

                for pn in pretrial_params:
                    if "ModelResult" in str(type(pretrial_lmfit)):
                        _pn_newval = (
                            pretrial_lmfit.params[pn].value
                            * random.randint(75, 110)
                            / 100
                        )
                    else:
                        _pn_newval = pretrial_best_row[pn].iloc[0]

                    if _pn_newval > 1e5:
                        _pn_newval = 1e3
                    if (pn.startswith("Q") or pn.startswith("C")) and _pn_newval > 0.1:
                        _pn_newval = standard_init_params[pn].value

                    if pn == "tau" and _pn_newval > 150:
                        _pn_newval = 5
                    if pn == "Aw" and _pn_newval > 500:
                        _pn_newval = 150

                    if (
                        pn == "Rs"
                        and not 0.3 < standard_init_params[pn].value / _pn_newval < 8
                    ):
                        _pn_newval = standard_init_params[pn].value

                    pretrial_params[pn].set(
                        value=_pn_newval,
                        min=params_model[pn].min,
                        max=params_model[pn].max,
                        vary=params_model[pn].vary,
                    )

                par_options.update(
                    {
                        "pretrial": {
                            "lmfit": pretrial_lmfit,
                            "params": pretrial_params,
                            "weights": pretrial_weights,
                        }
                    }
                )
        #            for pn in params_model:
        #                params_model[pn].set(value=pretrial_best.params[pn].value)

        pre_options_select = par_options.get("pretrial", par_options["standard"])

        best_trial = ""
        if "pretrial" in par_options.keys():
            if "lmfit.model" in str(type(pre_options_select.get("lmfit", False))):
                best_trial = pre_options_select["lmfit"]
        if not "lmfit.model" in str(type(best_trial)):
            best_trial = model_set.fit(
                Z_KKv,
                pre_options_select["params"],
                weights=weight_opts[pre_options_select["weights"]],
                ang=ang_KKv,
                method=_best_trial_method,
            )
        #        make_prefit_frame(EIS_data_KKvalid, best_trial,plot = True) #TODO
        init = model_set.eval(best_trial.params, ang=ang_KKv)
        #        _logger.info(f'''Prefit finished with method {best_trial.method} for {fit_run_arg.PAR_file.stem}
        #                    at {fit_run_arg.E_dc_RHE:0.2G} V on model {model_set.name},RedChiSqr = {best_trial.redchi:.3G} {best_trial.aic:.4G}''')
        # ADD RANDOM TO PRE-OUT params
        errRE_init, errIM_init = (init.real - Z_KKv.real) / abs(Z_KKv), (
            init.imag - Z_KKv.imag
        ) / abs(Z_KKv)

        # pre_options_select['weights']
        # weight_opts[pre_options_select['weights']]
        #        === Actual fitting step ===
        if use_default_wts:
            weights_name_select = default_weights
        else:
            weights_name_select = pre_options_select["weights"]

        weigths_select = weight_opts[weights_name_select]
        out = model_set.fit(
            Z_KKv,
            best_trial.params,
            ang=ang_KKv,
            weights=weigths_select,
            method=default_method,
            nan_policy="omit",
        )
        #        =====   ADDING WEIGHTS NAME TO FIT REPORT       =====
        fitrep_splt = out.fit_report(min_correl=0.6).split("\n[[Fit Statistics]]")
        fr_wghts = (
            fitrep_splt[0]
            + "\n"
            + f"    weights: {weights_name_select}"
            + "\n[[Fit Statistics]]"
        )
        out.fit_report_wts = fr_wghts + fitrep_splt[1]
        # ================ ===============#
        #        make_prefit_frame(EIS_data_KKvalid, out, prefix = 'pp',plot = 'Z') # TODO
        bf_good_low_high, (out_low_err, out_high_err), MSE_out = make_prefit_frame(
            EIS_data_KKvalid, out
        )
        #        _logger.info(f'''Fit out with method {out.method} for {fit_run_arg.PAR_file.stem} at {fit_run_arg.E_dc_RHE:0.2G} on model \
        #                    {model_set.name},ChiSqr = {out.chisqr:.3G}, RedChiSqr = {out.redchi:.3G}, {out.aic:.4G}''')
        ### == RETRY SECTION: when first fit is not good enough.. some initial values are re-set and fit is run again ... ####
        retry_attempt = "no"
        fit = out.eval(ang=ang_KKv)
        # === EIS REFIT USING RESIDUAL ===#
        errRE, errIM = (fit.real - Z_KKv.real) / abs(Z_KKv), (
            fit.imag - Z_KKv.imag
        ) / abs(Z_KKv)
        EIS_data_KKvalid_fit = pd.DataFrame()
        _mod_info = {
            "errRe": errRE,
            "errIm": errIM,
            "Model_EEC": str(modname),
            "Model_index": str(mod_indx),
            "Model_EEC_name": mod_inst_name,
            "lmfit_weights": out.weights,
            "lmfit_weights_nm": weights_name_select,
        }

        EIS_data_KKvalid_fit = EIS_data_KKvalid.assign(
            **{
                "INIT_Yim": (1 / init).imag,
                "INIT_Yre": (1 / init).real,
                "INIT_Zphase": [phase(i) for i in init],
                "INIT_errRe": errRE_init,
                "INIT_errIm": errIM_init,
                "FIT_Zre": fit.real,
                "FIT_Zim": fit.imag,
                "FIT_-Zim": -1 * fit.imag,
                "FIT_Yre": (1 / fit).real,
                "FIT_Yim": (1 / fit).imag,
                "FIT_Zmod": abs(fit),
                "FIT_Zphase": [phase(i) for i in fit],
                "FIT_Zangle": np.angle(fit, deg=True),
                "FIT_-Zangle": -1 * np.angle(fit, deg=True),
            }
        )

        EIS_data_KKvalid_fit = EIS_data_KKvalid_fit.assign(**_mod_info)

        # === Output organizing ===
        outP = out.best_values
        out_params_stderss = [
            (i + "_stderr", out.params.__getitem__(i).stderr) for i in out.params
        ]
        #        out_params_correl = [(i+'_correl', out.params.__getitem__(i).correl)  for i in out.params]
        outP.update(
            dict(
                zip(
                    [i[0] for i in out_params_stderss],
                    [i[1] for i in out_params_stderss],
                )
            )
        )
        #        outP.update(dict(zip([i[0] for i in out_params_correl],[i[1] for i in out_params_correl])))
        #        outP.update(EISgr_meta_add) # from beginning of function
        #        outP.update(linKK_pars) # from beginning of function after linKK validation
        outP.update(_mod_info)
        #        fit_run_arg_add = {EvRHE : fit_run_arg.E_dc_RHE, 'PAR_file' : fit_run_arg.PAR_file, 'PAR_date': EISgr_meta_combined['PAR_date'],
        #                     'RPM_DAC' : fit_run_arg.RPM_DAC,'Segment' : fit_run_arg.Segment, 'Model_EEC' : modname,'Model_index' : mod_indx}
        #        outP.update(fit_run_arg_add)
        #        outP.update(extraP)
        if "Qad" not in outP.keys():
            outP.update({"Qad": 0, "nAd": 0})
        if "Cdlp" not in outP.keys():
            outP.update({"Cdlp": 0, "nDL": 0})
        if "Rct" in outP.keys():
            outP.update({"Rct_kin": outP["Rct"] ** -1})
        xtra = {
            "Qad+Cdlp": outP["Qad"] + outP["Cdlp"],
            "lmfit_chiqsr": out.chisqr,
            "lmfit_redchi": out.redchi,
            "lmfit_aic": out.aic,
            "lmfit_bic": out.bic,
            "lmfit_method": out.method,
            "lmfit_message": out.message,
            "lmfit_out": out,
            "retry_attempt": retry_attempt,
            "lmfit_var_names": ", ".join(out.var_names),
            "lmfit_MSE": MSE_out,
            "lmfit_weights_nm": weights_name_select,
            "test_errcheck_msg": bf_good_low_high,
            "test_low": out_low_err,
            "test_high": out_high_err,
            "test_sum": out_high_err + out_low_err,
            "plot_color": mod_color,
            "fit_report": out.fit_report_wts,
            "BRUTE_FIT": 0,
            "FINAL_FIT": 1,
        }
        outP.update(xtra)
        outP = {
            key: val
            if not any(i in type(val).__name__ for i in ["list", "array"])
            else ", ".join([str(i) for i in val])
            for key, val in outP.items()
        }
        pars_mod_out = pd.DataFrame(outP, index=[EIS_data_KKvalid_fit.index[0]])
        #        if modname in all_PRETRIALS_grp.groups.keys():
        #            if not all_PRETRIALS_grp.get_group(modname).empty:
        #                BR_best_fit_weight_opt_out = all_PRETRIALS_grp.get_group(modname)
        #            pars_diff_cols = pars_mod_out.columns.difference(BR_best_fit_weight_opt_out.columns)
        #            pars_mod_out = pd.concat([pars_mod_out,BR_best_fit_weight_opt_out],sort=False,ignore_index=True)
        pars_lst.append(pars_mod_out)
        fit_data_valid_lst.append(EIS_data_KKvalid_fit)
    #        lmfit_out_lst.append({'Model_EEC' : modname,'Model_index' : mod_indx,'lmfit_out' : out})
    pars_models = pd.concat(pars_lst, sort=False, ignore_index=True)
    EIS_fit_data = pd.concat(fit_data_valid_lst, sort=False, ignore_index=True)
    return EIS_fit_data, pars_models
