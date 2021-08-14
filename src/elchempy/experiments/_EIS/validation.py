# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:43:56 2020

@author: DW
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy import stats
import logging
from pathlib import Path

from .impedancepy.impedance import validation as imppy_validation
from .impedancepy.impedance.models.circuits.fitting import rmse
from .impedancepy.impedance.models.circuits.elements import (
    circuit_elements,
    K,
)  # noqa

# freq,Zdata,Z_linKK,res_real,res_imag = EIS_data_valid['Frequency(Hz)'].values,EIS_data_valid.DATA_Z.values,EIS_data_valid.linKK_Z.values,EIS_data_valid.linKK_resRe.values,EIS_data_valid.linKK_resIm.values

# if __name__ == "__main__":
# from eis_run_ovv import EIS_Spectrum
# from repos.imppy.impedance import validation as imppy_validation
# from repos.imppy.impedance.models.circuits.fitting import rmse
import logging

_logger = logging.getLogger(__name__)


#%%
def rmse(a, b):
    """
    Code copied from : https://github.com/ECSHackWeek/impedance.py [2019-12-02]
    A function which calculates the root mean squared error
    between two vectors.
    Notes
    ---------
    .. math::
        RMSE = \\sqrt{\\frac{1}{n}(a-b)^2}
    """

    n = len(a)
    return np.linalg.norm(a - b) / np.sqrt(n)


def typeChecker(p, f, name, length):
    assert isinstance(p, list), "in {}, input must be of type list".format(name)
    for i in p:
        assert isinstance(
            i, (float, int, np.int32, np.float64)
        ), "in {}, value {} in {} is not a number".format(name, i, p)
    for i in f:
        assert isinstance(
            i, (float, int, np.int32, np.float64)
        ), "in {}, value {} in {} is not a number".format(name, i, f)
    assert len(p) == length, "in {}, input list must be length {}".format(name, length)
    return


def element_metadata(num_params, units):
    """decorator to store metadata for a circuit element
    Parameters
    ----------
    num_params : int
        number of parameters for an element
    units : list of str
        list of units for the element parameters
    """

    def decorator(func):
        def wrapper(p, f):
            typeChecker(p, f, func.__name__, num_params)
            return func(p, f)

        wrapper.num_params = num_params
        wrapper.units = units
        wrapper.__name__ = func.__name__
        wrapper.__doc__ = func.__doc__

        return wrapper

    return decorator


# === DATA VALIDATION Kramers-Kronig ====
class KramersKronigValidation:
    def __init__(self, spectrum, **kwargs):
        self.spectrum = spectrum
        #        if isinstance(spectrum, EIS_Spectrum):
        if type(spectrum).__name__ == "EIS_Spectrum":
            self._freqs = self.spectrum.EIS_data_freqlim["Frequency(Hz)"].to_numpy()
            self._Z = self.spectrum.EIS_data_freqlim["DATA_Z"].to_numpy()
            self._df_in = self.spectrum.EIS_data_freqlim
        elif type(spectrum).__name__ == "DataFrame":
            if all(hasattr(spectrum, i) for i in ["Frequency(Hz)", "DATA_Z"]):
                self._freqs = self.spectrum["Frequency(Hz)"].to_numpy()
                self._Z = self.spectrum["DATA_Z"].to_numpy()
                self._df_in = self.spectrum

        self.linKK(**kwargs)
        self.assigncols(**kwargs)

    #            print('is EIS_Spectrum')
    #        self.f = f
    #        self.Z = Z

    def _assigncols(df_in, Z_linKK, res_real, res_imag, *args, **kwargs):

        linKK = pd.DataFrame(
            {
                "linKK_Z": Z_linKK,
                "linKK_Zreal": Z_linKK.real,
                "linKK_Zimag": Z_linKK.imag,
                "linKK_Y": Z_linKK ** -1,
                "linKK_Yreal": (Z_linKK ** -1).real,
                "linKK_Yimag": (Z_linKK ** -1).imag,
                "linKK_resRe": res_real,
                "linKK_resIm": res_imag,
            },
            index=df_in.index,
        )
        if "suffix" in kwargs.keys():
            suffix = kwargs["suffix"]
            linKK.rename(
                columns=dict(zip(linKK.columns, [i + suffix for i in linKK.columns])),
                inplace=True,
            )

        return pd.concat([df_in, linKK], axis=1)

    def assigncols(self, **kwargs):
        Z_fit = self.linKK_result_dict["Z_fit"]
        linKK = pd.DataFrame(
            {
                "linKK_Z": Z_fit,
                "linKK_Zreal": Z_fit.real,
                "linKK_Zimag": Z_fit.imag,
                "linKK_-Zimag": -1 * Z_fit.imag,
                "linKK_Y": Z_fit ** -1,
                "linKK_Yreal": (Z_fit ** -1).real,
                "linKK_Yimag": (Z_fit ** -1).imag,
                "linKK_resRe": self.linKK_result_dict["resids_real"],
                "linKK_resIm": self.linKK_result_dict["resids_imag"],
            },
            index=self._df_in.index,
        )
        if "suffix" in kwargs.keys():
            suffix = kwargs["suffix"]
            linKK.rename(
                columns=dict(zip(linKK.columns, [i + suffix for i in linKK.columns])),
                inplace=True,
            )
        _linKK_df = pd.concat([self._df_in, linKK], axis=1).sort_values(
            by="Frequency(Hz)"
        )
        self.linKK_df = _linKK_df
        # pd.concat([self._df_in,linKK],axis=1)

    def linKK(self, *args, **kwargs):
        #        _freqs, _Z = self.spectrum.EIS_data_freqlim['Frequency(Hz)'].to_numpy(), self.spectrum.EIS_data_freqlim['DATA_Z'].to_numpy()
        _linKKkws = {}
        _kwskeys = ["c", "max_M", "fit_type", "add_cap"]
        if any(i in kwargs.keys() for i in _kwskeys):
            _linKKkws = {key: val for key, val in kwargs.items() if key in _kwskeys}

        self.linKK_result = imppy_validation.linKK(self._freqs, self._Z, **_linKKkws)
        self.linKK_result_dict = dict(
            zip(
                [
                    "M",
                    "mu",
                    "Z_fit",
                    "resids_real",
                    "resids_imag",
                    "fit_type",
                    "add_cap",
                ],
                self.linKK_result,
            )
        )

    def calc_meta_info(
        self, linKK_trimming_std_factor, res_max, prefit_suffix, _prefit
    ):
        linKK_limit_Re = np.min(
            [
                res_max,
                np.abs(_prefit["linKK_resRe" + prefit_suffix]).mean()
                + linKK_trimming_std_factor
                * np.abs(_prefit["linKK_resRe" + prefit_suffix]).std(),
            ]
        )
        linKK_limit_Im = np.min(
            [
                res_max,
                np.abs(_prefit["linKK_resIm" + prefit_suffix]).mean()
                + linKK_trimming_std_factor
                * np.abs(_prefit["linKK_resIm" + prefit_suffix]).std(),
            ]
        )
        linKK_limit_min_prefit = np.min([linKK_limit_Re, linKK_limit_Im])
        #    linKK_limit_mean = np.mean([linKK_limit_Re,linKK_limit_Im])
        #    linKK_limit_min = np.min([linKK_limit_Re,linKK_limit_Im])
        _prefit["linKK_res_mean" + prefit_suffix] = _prefit[
            ["linKK_resRe" + prefit_suffix, "linKK_resIm" + prefit_suffix]
        ].mean(axis=1)
        #    linKK_res_mean_lim = np.min([res_max, np.abs(_prefit['linKK_res_mean'+prefit_suffix]).mean()+linKK_trimming_std_factor*np.abs(_prefit['linKK_res_mean'+prefit_suffix]).std()])
        linKK_res_mean_lim = np.min(
            [
                res_max,
                np.min([linKK_limit_Re, linKK_limit_Im])
                + np.abs(_prefit["linKK_res_mean" + prefit_suffix]).std(),
            ]
        )
        #    np.abs(_prefit['linKK_res_mean'].mean()+linKK_trimming_std_factor*_prefit['linKK_res_mean'].std()
        #    _prefit['linKK_res_mean'].mean()+_prefit['linKK_res_mean'].std()
        resRe_Z, resIm_Z = (
            "linKK_resRe" + prefit_suffix + "_Zscore",
            "linKK_resIm" + prefit_suffix + "_Zscore",
        )
        _prefit[resRe_Z] = np.abs(stats.zscore(_prefit["linKK_resRe" + prefit_suffix]))
        _prefit[resIm_Z] = np.abs(stats.zscore(_prefit["linKK_resIm" + prefit_suffix]))
        resRe_Z_lim = _prefit[resRe_Z].mean() + 1.33 * _prefit[resRe_Z].std()
        resIm_Z_lim = _prefit[resIm_Z].mean() + 1.33 * _prefit[resIm_Z].std()


def get_KKvalid(
    EIS_data_freqlim,
    restype_set="Z",
    res_scale=1,
    res_ylim=0.25,
    res_max=0.1,
    **EIS_fit_kwargs,
):
    ### Performing lin KK prefit on freqlim data
    linKK_trimming_std_factor = EIS_fit_kwargs.get("linKK_trimming_factor", 1.60)
    try:
        prefit_suffix = "_prefit"
        #        _freqs, _Z = EIS_data_freqlim['Frequency(Hz)'].to_numpy(), EIS_data_freqlim['DATA_Z'].to_numpy()
        #        p_M, p_mu, p_Z_linKK, p_res_real, p_res_imag = KramersKronigValidation.linKK(_freqs, _Z,type_res=restype_set)
        #        _linKK_prefit = KramersKronigValidation.assigncols(EIS_data_freqlim,p_Z_linKK,p_res_real,p_res_imag,suffix=prefit_suffix)
        _prefit_result = KramersKronigValidation(
            EIS_data_freqlim, suffix=prefit_suffix, **EIS_fit_kwargs
        )
        _prefit = _prefit_result.linKK_df

        ### Slicing of lin KK prefit

        linKK_limit_Re = np.min(
            [
                res_max,
                np.abs(_prefit["linKK_resRe" + prefit_suffix]).mean()
                + linKK_trimming_std_factor
                * np.abs(_prefit["linKK_resRe" + prefit_suffix]).std(),
            ]
        )
        linKK_limit_Im = np.min(
            [
                res_max,
                np.abs(_prefit["linKK_resIm" + prefit_suffix]).mean()
                + linKK_trimming_std_factor
                * np.abs(_prefit["linKK_resIm" + prefit_suffix]).std(),
            ]
        )
        linKK_limit_min_prefit = np.min([linKK_limit_Re, linKK_limit_Im])
        #    linKK_limit_mean = np.mean([linKK_limit_Re,linKK_limit_Im])
        #    linKK_limit_min = np.min([linKK_limit_Re,linKK_limit_Im])
        _prefit["linKK_res_mean" + prefit_suffix] = _prefit[
            ["linKK_resRe" + prefit_suffix, "linKK_resIm" + prefit_suffix]
        ].mean(axis=1)
        #    linKK_res_mean_lim = np.min([res_max, np.abs(_prefit['linKK_res_mean'+prefit_suffix]).mean()+linKK_trimming_std_factor*np.abs(_prefit['linKK_res_mean'+prefit_suffix]).std()])
        linKK_res_mean_lim = np.min(
            [
                res_max,
                np.min([linKK_limit_Re, linKK_limit_Im])
                + np.abs(_prefit["linKK_res_mean" + prefit_suffix]).std(),
            ]
        )
        #    np.abs(_prefit['linKK_res_mean'].mean()+linKK_trimming_std_factor*_prefit['linKK_res_mean'].std()
        #    _prefit['linKK_res_mean'].mean()+_prefit['linKK_res_mean'].std()
        resRe_Z, resIm_Z = (
            "linKK_resRe" + prefit_suffix + "_Zscore",
            "linKK_resIm" + prefit_suffix + "_Zscore",
        )
        _prefit[resRe_Z] = np.abs(stats.zscore(_prefit["linKK_resRe" + prefit_suffix]))
        _prefit[resIm_Z] = np.abs(stats.zscore(_prefit["linKK_resIm" + prefit_suffix]))
        resRe_Z_lim = _prefit[resRe_Z].mean() + 1.33 * _prefit[resRe_Z].std()
        resIm_Z_lim = _prefit[resIm_Z].mean() + 1.33 * _prefit[resIm_Z].std()

        KKvalid_res = _prefit.loc[
            (
                (
                    (
                        (
                            np.abs(_prefit["linKK_resRe" + prefit_suffix])
                            < linKK_limit_Re
                        )
                        & (
                            np.abs(_prefit["linKK_resIm" + prefit_suffix])
                            < linKK_limit_Im
                        )
                        & (
                            np.abs(_prefit["linKK_res_mean" + prefit_suffix])
                            < linKK_res_mean_lim
                        )
                    )
                    & (_prefit["Frequency(Hz)"] > 2)
                )
                | (_prefit["Frequency(Hz)"] <= 2)
            )
        ]

        KKvalid_Zscore = _prefit.loc[
            (
                (
                    (
                        (_prefit[resRe_Z] < resRe_Z_lim)
                        & (_prefit[resIm_Z] < resIm_Z_lim)
                        & (_prefit["Frequency(Hz)"] > 2)
                    )
                    | (
                        (_prefit[resRe_Z] < 3.5 * resRe_Z_lim)
                        & (_prefit[resIm_Z] < 3.5 * resIm_Z_lim)
                        & (_prefit["Frequency(Hz)"] < 2)
                    )
                )
                & (_prefit.DATA_Yre < (1 / _prefit.linKK_Zreal_prefit.min()) * 1.1)
                | (
                    np.abs(_prefit["linKK_res_mean" + prefit_suffix])
                    < linKK_res_mean_lim * 0.2
                )
            )
        ]
        #    KKvalid_Zscore_2 = KKvalid_res.loc[((KKvalid_res[resRe_Z] < resRe_Z_lim)
        #                                & (KKvalid_res[resIm_Z] < resIm_Z_lim)
        #                                & (KKvalid_res.DATA_Yre < (1/KKvalid_res.linKK_Zreal_prefit.min())*1.01))]
        KKvalid_Zscore_3 = _prefit.loc[
            (
                (
                    np.abs(_prefit["linKK_resRe" + prefit_suffix])
                    < linKK_limit_min_prefit
                )
                & (
                    np.abs(_prefit["linKK_resIm" + prefit_suffix])
                    < linKK_limit_min_prefit
                )
            )
        ]
        #     & (_prefit[resRe_Z] < resRe_Z_lim)
        KKvalid = _prefit.loc[
            ((KKvalid_Zscore.index) & (KKvalid_res.index)) | (KKvalid_Zscore_3.index)
        ]

        # _first_KK_valid = KramersKronigValidation(KKvalid, suffix='first_valid',add_cap=True)
        # KKvalid = _first_KK_valid.linKK_df
        # linKK_valid_result = KramersKronigValidation(KKvalid)
        def test_plots():
            fig, ax = plt.subplots()
            KKvalid.plot(x="linKK_Yreal_prefit", y="linKK_Yimag_prefit", ax=ax)
            KKvalid.plot(x="DATA_Yre", y="DATA_Yim", ax=ax, kind="scatter")

            fig, ax = plt.subplots()
            KKvalid.plot(x="linKK_Zreal_prefit", y="linKK_-Zimag_prefit", ax=ax)
            KKvalid.plot(x="DATA_Zre", y="DATA_-Zim", ax=ax, kind="scatter")

            KKvalid.plot(x="Frequency(Hz)", y=[resIm_Z, resRe_Z], logx=True)
            KKvalid.plot(
                x="Frequency(Hz)",
                y=["linKK_resRe" + prefit_suffix, "linKK_resIm" + prefit_suffix],
                logx=True,
            )

            KKvalid.plot(x="Frequency(Hz)", y=[resIm_Z, resRe_Z])

        # One more slice for extremely bad Z scores
        KKvalid = KKvalid.loc[
            (
                (
                    (
                        KKvalid[resIm_Z]
                        < (KKvalid[resIm_Z].mean() + KKvalid[resIm_Z].std() * 5)
                    )
                    | (
                        np.abs(KKvalid["linKK_resIm" + prefit_suffix])
                        < linKK_limit_min_prefit * 1.5
                    )
                )
                & (
                    (
                        KKvalid[resRe_Z]
                        < (KKvalid[resRe_Z].mean() + KKvalid[resRe_Z].std() * 5)
                    )
                    | (
                        np.abs(KKvalid["linKK_resRe" + prefit_suffix])
                        < linKK_limit_min_prefit * 2.5
                    )
                )
            )
        ]
        _low_freq_Zmod_mean = KKvalid[KKvalid["Frequency(Hz)"] < 5].DATA_Zmod.mean()
        if _low_freq_Zmod_mean:
            KKvalid = KKvalid.loc[KKvalid.DATA_Zmod < _low_freq_Zmod_mean * 1.66]

        linKK_invalid_prefit = _prefit.loc[
            (~_prefit["Frequency(Hz)"].isin(KKvalid["Frequency(Hz)"]))
        ]

        #    _prefit.loc[(KKvalid_Zscore.index) & (~_prefit['Frequency(Hz)'].isin(KKvalid['Frequency(Hz)']))]
        # Setting Valid == False for lin KK on raw data
        #    plot_linKK(_prefit['Frequency(Hz)'].values,_prefit.DATA_Z.values,
        #                                     _prefit['linKK_Z_prefit'],_prefit['linKK_resRe'+prefit_suffix],
        #                                     _prefit['linKK_resIm'+prefit_suffix],
        #                                     meta=EISgr_data_meta,linKK=[p_M,p_mu],type_res=restype_set)
        # plot 1st KK validation for testing purposes
        #    linKK = pd.DataFrame({'linKK_Z' : Z_linKK, 'linKK_Zreal' : Z_linKK.real,'linKK_Zimag' : Z_linKK.imag,
        #                          'linKK_Y' : Z_linKK**-1, 'linKK_Yreal' : (Z_linKK**-1).real,'linKK_Yimag' : (Z_linKK**-1).imag,
        #                          'linKK_resRe' : res_real, 'linKK_resIm' : res_imag }, index=EIS_data_freqlim.index)
        #    EIS_data_linKK = pd.concat([EIS_data_freqlim,linKK],axis=1)
        #    linKK_limit = 0.0475
        #        M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(KKvalid['Frequency(Hz)'].values,KKvalid['DATA_Z'].values,type_res=restype_set,max_M=50)
        linKK_valid_result = KramersKronigValidation(KKvalid, **EIS_fit_kwargs)
        _linKK_valid_fit = linKK_valid_result.linKK_df
        (
            M,
            mu,
            Z_linKK,
            res_real,
            res_imag,
            fit_type,
            add_cap,
        ) = linKK_valid_result.linKK_result

    except Exception as e:
        _logger.error(
            f'Lin KK fit error {EIS_data_freqlim.iloc[0].to_dict().get("PAR_file","PF look up error"), {e}}'
        )
        M, mu = 0, 0
        fit_type, add_cap = "complex", False
        KKvalid, Z_linKK = EIS_data_freqlim, EIS_data_freqlim["DATA_Z"].to_numpy()
        _linKK_valid_fit = KKvalid
        res_real, res_imag = 2 * [np.array([1] * EIS_data_freqlim["DATA_Z"].shape[0])]

    EIS_data_valid_linKK = _linKK_valid_fit
    #    KramersKronigValidation.assigncols(KKvalid,Z_linKK,res_real,res_imag)
    linKK_pars = {
        "linKK_M": M,
        "linKK_mu": mu,
        "linKK_trimming_factor": linKK_trimming_std_factor,
        "linKK_limit_Im": linKK_limit_Im,
        "linKK_limit_Re": linKK_limit_Re,
        "type_res": restype_set,
        "res_scale": res_scale,
        "res_ylim": res_ylim,
        "linKK_invalid_freqs": linKK_invalid_prefit["Frequency(Hz)"].to_numpy(),
        "fit_type": fit_type,
        "add_cap": add_cap,
    }
    linKK_pars = {
        key if key.startswith("linKK_") else f"linKK_{key}": val
        for key, val in linKK_pars.items()
    }

    return EIS_data_valid_linKK, linKK_invalid_prefit, linKK_pars


#    EIS_data_valid_KK = KramersKronigValidation.assigncols(EIS_data,Z_linKK,res_real,res_imag)
#    plot_linKK(EIS_data_valid, EIS_data, res_scale=1,meta=EISgr_data_meta,linKK=[M,mu],type_res=restype_set,save_target=EIS_outpath_linKK_target,plot_prefit=True,KK_valid_limit=(linKK_limit_Re, linKK_limit_Im))
#    return
def prep_GP_DRT_raw_data(EIS_data_raw):
    _idxs = EIS_data_raw["Frequency(Hz)"].to_numpy()
    _invald_idxs = EIS_data_raw.loc[EIS_data_raw.Valid == False][
        "Frequency(Hz)"
    ].to_numpy()
    _rem = []
    for _idx in _invald_idxs:
        if _idx == _idxs[0]:
            _rem.append(_idx)
            _idxs = _idxs[1:]

    for _idx in _invald_idxs[::-1]:
        if _idx == _idxs[-1]:
            _rem.append(_idx)
            _idxs = _idxs[:-1]
    GP_DRT_valid = EIS_data_raw.loc[
        ~EIS_data_raw["Frequency(Hz)"].isin(_rem)
    ].sort_values(by="Frequency(Hz)", ascending=True)
    return GP_DRT_valid


def linKK(f, Z, c=0.85, max_M=50, fit_type="real", add_cap=False):
    """A method for implementing the Lin-KK test for validating linearity [1]

    Parameters
    ----------
    f: np.ndarray
        measured frequencies
    Z: np.ndarray of complex numbers
        measured impedances
    c: np.float
        cutoff for mu
    max_M: int
        the maximum number of RC elements
    fit_type: str
        selects which components of data are fit ('real', 'imag', or
        'complex')
    add_cap: bool
        option to add a serial capacitance that helps validate data with no
        low-frequency intercept

    Returns
    -------
    M: int
        number of RC elements used
    mu: np.float
        under- or over-fitting measure
    Z_fit: np.ndarray of complex numbers
        impedance of fit at input frequencies
    resids_real: np.ndarray
        real component of the residuals of the fit at input frequencies
    resids_imag: np.ndarray
        imaginary component of the residuals of the fit at input frequencies


    Notes
    -----

    The lin-KK method from Schönleber et al. [1] is a quick test for checking
    the
    validity of EIS data. The validity of an impedance spectrum is analyzed by
    its reproducibility by a Kramers-Kronig (KK) compliant equivalent circuit.
    In particular, the model used in the lin-KK test is an ohmic resistor,
    :math:`R_{Ohm}`, and :math:`M` RC elements.

    .. math::

        \\hat Z = R_{Ohm} + \\sum_{k=1}^{M} \\frac{R_k}{1 + j \\omega \\tau_k}

    The :math:`M` time constants, :math:`\\tau_k`, are distributed
    logarithmically,

    .. math::
        \\tau_1 = \\frac{1}{\\omega_{max}} ; \\tau_M = \\frac{1}{\\omega_{min}}
        ; \\tau_k = 10^{\\log{(\\tau_{min}) + \\frac{k-1}{M-1}\\log{{(
            \\frac{\\tau_{max}}{\\tau_{min}}}})}}

    and are not fit during the test (only :math:`R_{Ohm}` and :math:`R_{k}`
    are free parameters).

    In order to prevent under- or over-fitting, Schönleber et al. propose using
    the ratio of positive resistor mass to negative resistor mass as a metric
    for finding the optimal number of RC elements.

    .. math::

        \\mu = 1 - \\frac{\\sum_{R_k \\ge 0} |R_k|}{\\sum_{R_k < 0} |R_k|}

    The argument :code:`c` defines the cutoff value for :math:`\\mu`. The
    algorithm starts at :code:`M = 3` and iterates up to :code:`max_M` until a
    :math:`\\mu < c` is reached. The default of 0.85 is simply a heuristic
    value based off of the experience of Schönleber et al., but a lower value
    may give better results.

    If the argument :code:`c` is :code:`None`, then the automatic determination
    of RC elements is turned off and the solution is calculated for
    :code:`max_M` RC elements. This manual mode should be used with caution as
    under- and over-fitting should be avoided.

    [1] Schönleber, M. et al. A Method for Improving the Robustness of
    linear Kramers-Kronig Validity Tests. Electrochimica Acta 131, 20–27 (2014)
    `doi: 10.1016/j.electacta.2014.01.034
    <https://doi.org/10.1016/j.electacta.2014.01.034>`_.

    """

    if c is not None:
        M = 0
        mu = 1
        while mu > c and M <= max_M:
            M += 1
            ts = get_tc_distribution(f, M)
            elements, mu = fit_linKK(f, ts, M, Z, fit_type, add_cap)

            if M % 10 == 0:
                print(M, mu, rmse(eval_linKK(elements, ts, f), Z))
    else:
        M = max_M
        ts = get_tc_distribution(f, M)
        elements, mu = fit_linKK(f, ts, M, Z, fit_type, add_cap)

    Z_fit = eval_linKK(elements, ts, f)
    resids_real = residuals_linKK(elements, ts, Z, f, residuals="real")
    resids_imag = residuals_linKK(elements, ts, Z, f, residuals="imag")
    return M, mu, Z_fit, resids_real, resids_imag


def get_tc_distribution(f, M):
    """Returns the distribution of time constants for the linKK method"""

    t_max = 1 / (2 * np.pi * np.min(f))
    t_min = 1 / (2 * np.pi * np.max(f))
    ts = np.zeros(shape=(M,))
    ts[0] = t_min
    ts[-1] = t_max
    if M > 1:
        for k in range(2, M):
            ts[k - 1] = 10 ** (
                np.log10(t_min) + ((k - 1) / (M - 1)) * np.log10(t_max / t_min)
            )
    return ts


def fit_linKK(f, ts, M, Z, fit_type="real", add_cap=False):
    """Fits the linKK model using linear regression

    Parameters
    ----------
    f: np.ndarray
        measured frequencies
    ts: np.ndarray
        logarithmically spaced time constants of RC elements
    M: int
        the number of RC elements
    Z: np.ndarray of complex numbers
        measured impedances
    fit_type: str
        selects which components of data are fit ('real', 'imag', or
        'complex')
    add_cap: bool
        option to add a serial capacitance that helps validate data with no
        low-frequency intercept

    Returns
    -------
    elements: np.ndarray
        values of fit :math:`R_k` in RC elements and series :math:`R_0`,
        L, and optionally C.
    mu: np.float
        under- or over-fitting measure

    Notes
    -----
    Since we have a system of equations, :math:`Ax ~= b`, that's linear wrt
    :math:`R_k`, we can fit the model by calculating the pseudo-inverse of A.
    :math:`Ax` is our model fit, :math:`\\hat{Z}`, and :math:`b` is the
    normalized real or imaginary component of the impedance data,
    :math:`Re(Z)/|Z|` or :math:`Im(Z)/|Z|`, respectively.

    :math:`\\hat{Z} = R_0 + \\sum^M_{k=1}(R_k / |Z|(1 + j * w * \\tau_k))`.
    :math:`x` is an (M+1) :math:`\\times` 1 matrix where the first row
    contains :math:`R_0` and subsequent rows contain :math:`R_k` values.
    A is an N :math:`\\times` (M+1) matrix, where N is the number of data
    points, and M is the number of RC elements.

    Examples
    --------

    Fitting the real part of data, the first column of A contains
    values of :math:`\\frac{1}{|Z|}`, the second column contains
    :math:`Re(1 / |Z| (1 + j * w * \\tau_1))`, the third contains
    :math:`Re(1 / |Z| (1 + j * w * \\tau_2))` and so on. The :math:`R_k` values
    within the x matrix are found using :code:`numpy.linalg.pinv` when
    fit_type = 'real' or 'imag'. When fit_type = 'complex' the coefficients are
    found "manually" using :math:`r = ||A'x - b'||^2 + ||A''x - b'||^2`
    according to Eq 14 of Schonleber [1].

    [1] Schönleber, M. et al. A Method for Improving the Robustness of
    linear Kramers-Kronig Validity Tests. Electrochimica Acta 131, 20–27 (2014)
    `doi: 10.1016/j.electacta.2014.01.034
    <https://doi.org/10.1016/j.electacta.2014.01.034>`_.
    """

    w = 2 * np.pi * f

    # Fitting model has M RC elements plus 1 series resistance and 1 series
    # inductance
    a_re = np.zeros((f.size, M + 2))
    a_im = np.zeros((f.size, M + 2))

    if add_cap:
        a_re = np.zeros((f.size, M + 3))
        a_im = np.zeros((f.size, M + 3))

        # Column for series capacitance. Real part = 0.
        a_im[:, -2] = -1 / (w * np.abs(Z))

    # Column for series resistance, R_o in model. Imaginary part = 0.
    a_re[:, 0] = 1 / np.abs(Z)

    # Column for series inductance to capture inevitable contributions from
    # the measurement system. Real part = 0.
    a_im[:, -1] = w / np.abs(Z)

    # Columns for series RC elements
    for i, tau in enumerate(ts):
        a_re[:, i + 1] = K([1, tau], f).real / np.abs(Z)
        a_im[:, i + 1] = K([1, tau], f).imag / np.abs(Z)

    if fit_type == "real":
        elements = np.linalg.pinv(a_re).dot(Z.real / np.abs(Z))

        # After fitting real part, need to use imaginary component of fit to
        # find values of series inductance and capacitance
        a_im = np.zeros((f.size, 2))
        a_im[:, -1] = w / np.abs(Z)
        if add_cap:
            a_im[:, -2] = -1 / (w * np.abs(Z))
            elements[-2] = 1e-18  # nullifies series C without dividing by 0

        Z_fit_re = eval_linKK(elements, ts, f)
        coefs = np.linalg.pinv(a_im).dot((Z.imag - Z_fit_re.imag) / np.abs(Z))

        if add_cap:
            elements[-2:] = coefs
        else:
            elements[-1] = coefs[-1]
    elif fit_type == "imag":
        elements = np.linalg.pinv(a_im).dot(Z.imag / np.abs(Z))

        # Calculates real part of impedance from fitting to imaginary parts
        # without ohmic resistance, i.e. only the real parts of series RC
        # elements.
        z_re = eval_linKK(elements, ts, f)

        # Weighting used in Boukamp et al - "A Linear Kronig-Kramers
        # Transform for Immittance Data Validation" 1995, J. Electrochem
        # Soc. 142 (6)
        ws = 1 / (Z.real ** 2 + Z.imag ** 2)

        # Finds ohmic resistance for imaginary part fit according to Eq 7
        # of Boukamp et al.
        elements[0] = np.sum(ws * (Z.real - z_re.real)) / np.sum(ws)
    elif fit_type == "complex":
        # x = (A*•A)^-1
        # y = A*•b
        # Pseudoinsverse, A^+ = (A*•A)^-1•A* and R_k values = A^+•b
        x = np.linalg.inv(a_re.T.dot(a_re) + a_im.T.dot(a_im))
        y = a_re.T.dot(Z.real / np.abs(Z)) + a_im.T.dot(Z.imag / np.abs(Z))
        elements = x.dot(y)
    else:
        raise ValueError(
            "Invalid choice of fit_type, please choose from "
            "'real', 'imag', or 'complex'"
        )

    if add_cap:
        mu = calc_mu(elements[1:-2])
    else:
        mu = calc_mu(elements[1:-1])

    return elements, mu


def eval_linKK(elements, ts, f):
    """Builds a circuit of RC elements to be used in LinKK"""
    circuit_string = "s([R({},{}),".format([elements[0]], f.tolist())

    for (Rk, tk) in zip(elements[1:], ts):
        circuit_string += f"K({[Rk, tk]},{f.tolist()}),"

    circuit_string += "L({},{}),".format([elements[-1]], f.tolist())
    if elements.size == (ts.size + 3):
        circuit_string += "C({},{}),".format([1 / elements[-2]], f.tolist())

    circuit_string = circuit_string.strip(",")
    circuit_string += "])"

    return eval(circuit_string, circuit_elements)


def residuals_linKK(elements, ts, Z, f, residuals="real"):
    """Calculates the residual between the data and a LinKK fit"""

    err = (Z - eval_linKK(elements, ts, f)) / np.abs(Z)

    if residuals == "real":
        return err.real
    elif residuals == "imag":
        return err.imag
    elif residuals == "both":
        z1d = np.zeros(Z.size * 2, dtype=np.float64)
        z1d[0 : z1d.size : 2] = err.real
        z1d[1 : z1d.size : 2] = err.imag

        return z1d


def calc_mu(Rs):
    """Calculates mu for use in LinKK"""

    neg_sum = sum(abs(x) for x in Rs if x < 0)
    pos_sum = sum(abs(x) for x in Rs if x >= 0)

    return 1 - neg_sum / pos_sum


def plot_prefit(_linKK_prefit):

    fig, ax = plt.subplot()


def old_linKK(f, Z, c=0.85, max_M=50, type_res="Z"):
    """Code copied from https://github.com/ECSHackWeek/impedance.py/blob/master/impedance/validation.py [2019-12-02]
    A method for implementing the Lin-KK test for validating linearity [1]
    Parameters
    ----------
    f: np.ndarray
        measured frequencies
    Z: np.ndarray of complex numbers
        measured impedances
    c: np.float
        cutoff for mu
    max_M: int
        the maximum number of RC elements
    Returns
    -------
    mu: np.float
        under- or over-fitting measure
    residuals: np.ndarray of complex numbers
        the residuals of the fit at input frequencies
    Z_fit: np.ndarray of complex numbers
        impedance of fit at input frequencies
    Notes
    -----
    The lin-KK method from Schönleber et al. [1] is a quick test for checking
    the
    validity of EIS data. The validity of an impedance spectrum is analyzed by
    its reproducibility by a Kramers-Kronig (KK) compliant equivalent circuit.
    In particular, the model used in the lin-KK test is an ohmic resistor,
    :math:`R_{Ohm}`, and :math:`M` RC elements.
    .. math::
        \\hat Z = R_{Ohm} + \\sum_{k=1}^{M} \\frac{R_k}{1 + j \\omega \\tau_k}
    The :math:`M` time constants, :math:`\\tau_k`, are distributed
    logarithmically,
    .. math::
        \\tau_1 = \\frac{1}{\\omega_{max}} ; \\tau_M = \\frac{1}{\\omega_{min}}
        ; \\tau_k = 10^{\\log{(\\tau_{min}) + \\frac{k-1}{M-1}\\log{{(
            \\frac{\\tau_{max}}{\\tau_{min}}}})}}
    and are not fit during the test (only :math:`R_{Ohm}` and :math:`R_{k}`
    are free parameters).
    In order to prevent under- or over-fitting, Schönleber et al. propose using
    the ratio of positive resistor mass to negative resistor mass as a metric
    for finding the optimal number of RC elements.
    .. math::
        \\mu = 1 - \\frac{\\sum_{R_k \\ge 0} |R_k|}{\\sum_{R_k < 0} |R_k|}
    The argument :code:`c` defines the cutoff value for :math:`\\mu`. The
    algorithm starts at :code:`M = 3` and iterates up to :code:`max_M` until a
    :math:`\\mu < c` is reached. The default of 0.85 is simply a heuristic
    value based off of the experience of Schönleber et al.
    If the argument :code:`c` is :code:`None`, then the automatic determination
    of RC elements is turned off and the solution is calculated for
    :code:`max_M` RC elements. This manual mode should be used with caution as
    under- and over-fitting should be avoided.
    [1] Schönleber, M. et al. A Method for Improving the Robustness of
    linear Kramers-Kronig Validity Tests. Electrochimica Acta 131, 20–27 (2014)
    `doi: 10.1016/j.electacta.2014.01.034
    <https://doi.org/10.1016/j.electacta.2014.01.034>`_.
    """

    def calc_mu(Rs):
        """Calculates mu for use in LinKK"""

        neg_sum = sum(abs(x) for x in Rs if x < 0)
        pos_sum = sum(abs(x) for x in Rs if x >= 0)
        return 1 - neg_sum / pos_sum

    def fitLinKK(f, ts, M, Z):
        """Fits the linKK model using scipy.optimize.least_squares"""
        initial_guess = np.append(
            min(np.real(Z)),
            np.ones(shape=(M,)) * ((max(np.real(Z)) - min(np.real(Z))) / M),
        )

        result = least_squares(
            residuals_linKK,
            initial_guess,
            method="lm",
            args=(ts, Z, f, "both", type_res),
            ftol=1e-13,
            gtol=1e-10,
        )
        p_values = result["x"]
        mu = calc_mu(p_values[1:])
        return p_values, mu

    def circuit_element_s(series):
        """sums elements in series
        Notes
        ---------
        .. math::
            Z = Z_1 + Z_2 + ... + Z_n
        """
        z = len(series[0]) * [0 + 0 * 1j]
        for elem in series:
            z += elem
        return z

    @element_metadata(num_params=1, units=["Ohm"])
    def circuit_element_R(p, f):
        """defines a resistor
        Notes
        ---------
        .. math::
            Z = R
        """
        return np.array(len(f) * [p[0]])

    @element_metadata(num_params=2, units=["Ohm", "sec"])
    def circuit_element_K(p, f):
        """An RC element for use in lin-KK model
        Notes
        -----
        .. math::
            Z = \\frac{R}{1 + j \\omega \\tau_k}
        """
        omega = np.array(f)
        return p[0] / (1 + 1j * omega * p[1])

    def eval_linKK(Rs, ts, f):
        """Builds a circuit of RC elements to be used in LinKK"""
        s, R, K = circuit_element_s, circuit_element_R, circuit_element_K
        circuit_string = "s([R({},{}),".format([Rs[0]], f.tolist())

        for i, (Rk, tk) in enumerate(zip(Rs[1:], ts)):
            circuit_string += "K({},{}),".format([Rk, tk], f.tolist())

        circuit_string = circuit_string.strip(",")
        circuit_string += "])"
        return eval(circuit_string)

    def residuals_linKK(Rs, ts, Z, f, residuals="both", type_residual="Z"):
        """Calculates the residual between the data and a LinKK fit"""

        err = Z - eval_linKK(Rs, ts, f)
        Y = Z ** -1
        errY = Y - 1 / eval_linKK(Rs, ts, f)

        if "Z" in type_residual:
            err_take = 0
        elif "Y" in type_residual:
            err_take = 1

        if residuals == "real":
            err_real = [err.real / (np.abs(Z)), errY.real / (np.abs(Y))]
            return err_real[err_take]
        elif residuals == "imag":
            err_imag = err.imag / (np.abs(Z)), errY.imag / (np.abs(Y))
            return err_imag[err_take]
        elif residuals == "both":
            z1d = np.zeros(Z.size * 2, dtype=np.float64)
            z1d[0 : z1d.size : 2] = err.real / (np.abs(Z))
            z1d[1 : z1d.size : 2] = err.imag / (np.abs(Z))

            y1d = np.zeros(Y.size * 2, dtype=np.float64)
            y1d[0 : y1d.size : 2] = errY.real / (np.abs(Y))
            y1d[1 : y1d.size : 2] = errY.imag / (np.abs(Y))
            err_both = [z1d, y1d]
            return err_both[err_take]

    def get_tc_distribution(f, M):
        """Returns the distribution of time constants for the linKK method"""

        t_max = 1 / np.min(f)
        t_min = 1 / np.max(f)

        ts = np.zeros(shape=(M,))
        ts[0] = t_min
        ts[-1] = t_max
        if M > 1:
            for k in range(2, M):
                ts[k - 1] = 10 ** (
                    np.log10(t_min) + ((k - 1) / (M - 1)) * np.log10(t_max / t_min)
                )

        ts *= 2 * np.pi
        return ts

    if c is not None:
        M = 0
        mu = 1
        while mu > c and M <= max_M:
            M += 1
            ts = get_tc_distribution(f, M)
            p_values, mu = fitLinKK(f, ts, M, Z)

            if M % 10 == 0:
                pass  # print(M, mu, rmse(eval_linKK(p_values, ts, f), Z))
    else:
        M = max_M
        ts = get_tc_distribution(f, M)
        p_values, mu = fitLinKK(f, ts, M, Z)

    return (
        M,
        mu,
        eval_linKK(p_values, ts, f),
        residuals_linKK(p_values, ts, Z, f, residuals="real"),
        residuals_linKK(p_values, ts, Z, f, residuals="imag"),
    )
