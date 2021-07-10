# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:00:16 2020

@author: DWXMG
"""

import sys
from collections import namedtuple, OrderedDict
from pathlib import Path
from itertools import combinations
import random

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

from lmfit import Parameters

print(f"Name: {__name__} for file {__file__}")

if __name__ == "__main__":
    from validation import get_KKvalid
    from models import EEC_models_index, params_extra_setting
    from plotting import plot_linKK, EIS_Trimming_plot

    sys.path.append(str(Path(__file__).parent.parent.parent.joinpath("runEC")))
    #    sys.path.append(Path(__file__).parent.parent.parent.parent)
    from EC_logging_config import start_logging

    #    from EC_logging_config import start_logging
    #    sys.path.append(Path(__file__).parent.parent.joinpath('runEC'))
    #    from runEC.EC_logging_config import start_logging
    logger = start_logging(__name__)


elif __name__ in "fitting":
    from validation import get_KKvalid
    from models import EEC_models_index, params_extra_setting
    from plotting import plot_linKK, EIS_Trimming_plot, EIS_plotting_per_EV

    FH_path = Path(__file__).parent.parent.parent.parent.joinpath("FileHelper")
    sys.path.append(str(FH_path))
    from FileHelper.FileFunctions import FileOperations

    sys.path.append(Path(__file__).parent.parent.parent)
    from ECpy.runEC.EC_logging_config import start_logging

    #    from EC_logging_config import start_logging
    #    sys.path.append(Path(__file__).parent.parent.joinpath('runEC'))
    #    from runEC.EC_logging_config import start_logging
    logger = start_logging(__name__)


else:
    print(__name__)
    from .validation import get_KKvalid
    from .models import EEC_models_index, params_extra_setting
    from .plotting import plot_linKK, EIS_Trimming_plot, EIS_plotting_per_EV

    #    from .plotting import EIS_plotting_per_EV, EIS_Trimming_plot, EIS_plotting_EvRHE

    FH_path = Path(__file__).parent.parent.parent.parent.joinpath("FileHelper")
    sys.path.append(str(FH_path))
    from FileHelper.FileFunctions import FileOperations

    sys.path.append(Path(__file__).parent.parent.joinpath("runEC"))
    from runEC.EC_logging_config import start_logging

    logger = start_logging(__name__)


def fitting_recheck_params(fit_run_arg, modname, params_model, **EIS_fit_kwargs):
    PF, E, RPM = str(fit_run_arg.PAR_file), fit_run_arg.E_dc_RHE, fit_run_arg.RPM_DAC
    _key = (PF, E, RPM, modname)
    _get_params = pd.DataFrame()
    # ['PAR_file',EvRHE,'RPM_DAC','Model_EEC']
    bad_grp, good_grp = EIS_fit_kwargs.get("EIS_recheck_bad_fits"), EIS_fit_kwargs.get(
        "EIS_recheck_good_fits"
    )
    sugg_grp = EIS_fit_kwargs.get("EIS_recheck_bad_fits_suggestions")
    recheck_msg = ""
    if [i for i in good_grp.groups if _key == i]:
        recheck_msg += f"Prefit recheck in good keys"
        _get_params = good_grp.get_group(_key)
    elif [i for i in bad_grp.groups if _key == i]:
        recheck_msg += f"Prefit recheck in bad keys {_key} and"
        _sugg_match = [i for i in sugg_grp.groups if _key == i]
        if _sugg_match:
            recheck_msg += f" taking suggestions."
            _get_params = sugg_grp.get_group(_key)
        else:
            recheck_msg += f" not in suggestions."
    #                         logger.warning(f'Prefit recheck bad key {_key} not in suggestions')
    else:
        recheck_msg += f"Prefit recheck keys not in good or bad {_key}"
    #    logger.warning(recheck_msg)
    return _get_params, recheck_msg


#    if all([i in params_model for i in ['Rct','Rorr']]):
#        if params_model['Rct'] < params_model['Rorr']:
#            _Rorr = params_model['Rct']
#            params_model['Rct'].set(value=params_model['Rorr'])
#            params_model['Rorr'].set(value=_Rorr )
# TODO finish brute preparation


def random_combination(iterable, r):
    i = 0
    pool = tuple(iterable)
    n = len(pool)
    rng = range(n)
    while i < r:
        i += 1
        yield [pool[j] for j in random.sample(rng, r)]


def random_combo_select(_brutepardict, number_guesses):
    _guesses = []
    while len(_guesses) < number_guesses:
        _combo = [
            (key, random.sample(val["options"], 1)[0])
            for key, val in _brutepardict.items()
        ]
        if _combo in _guesses:
            continue
        if "Rct" and "Rorr" in _brutepardict.keys():
            if (
                _combo[_brutepardict.get("Rct")["index"]][1]
                < _combo[_brutepardict.get("Rorr")["index"]][1]
            ):
                _guesses.append(_combo)
        #                    brutepar_lst = [i for i in brutepar_lst if i[_brutepardict.get('Rct')][1] < i[_brutepardict.get('Rorr')][1]]
        elif "Qad" and "Cdlp" in _brutepardict.keys():
            if (
                _combo[_brutepardict.get("Qad")["index"]][1]
                > _combo[_brutepardict.get("Cdlp")["index"]][1]
            ):
                if "Q3" in _brutepardict.keys():
                    if (
                        _combo[_brutepardict.get("Q3")["index"]][1]
                        < _combo[_brutepardict.get("Cdlp")["index"]][1]
                    ):
                        _guesses.append(_combo)
                    else:
                        pass
                else:
                    _guesses.append(_combo)
            else:
                pass
    return _guesses


def prepare_brute_pars(params_model, number_guess):
    #    mod_combs = {}
    #    for mod_indx, model_set in EEC_models_index(): #
    #    modname = model_set.name
    #    params_model = model_set.make_params()
    _t, _tl = [], 0
    _brutepardict = OrderedDict()
    if "Rct" in params_model:
        #                    _brutepars.update({'Rct' : [5,5E2,1E3,10E3,100E3]})
        _opts = [1, 5, 10, 30, 50, 1e3, 10e3, 100e3]
        _t += [("Rct", i) for i in _opts]
        _brutepardict.update({"Rct": {"index": _tl, "options": _opts}})
        _tl += 1
    if "Rorr" in params_model:
        _opts = [10, 100, 500, 1e3, 3e3, 5e3, 20e3, 200e3]
        _t += [("Rorr", i) for i in _opts]
        _brutepardict.update({"Rorr": {"index": _tl, "options": _opts}})
        _tl += 1

    if "nAd" in params_model:
        _opts = [0.4, 0.6, 0.9][::-1]
        _t += [("nAd", i) for i in _opts]
        _brutepardict.update({"nAd": {"index": _tl, "options": _opts}})
        _tl += 1
    if "nDL" in params_model:
        _opts = [0.4, 0.6, 0.9][::-1]
        _t += [("nDL", i) for i in _opts]
        _brutepardict.update({"nDL": {"index": _tl, "options": _opts}})
        _tl += 1

    if "Qad" in params_model:
        _opts = [1e-3, 5e-3, 25e-3]
        _t += [("Qad", i) for i in _opts]
        _brutepardict.update({"Qad": {"index": _tl, "options": _opts}})
        _tl += 1
    if "Cdlp" in params_model:
        _opts = [5e-5, 1e-3, 2e-2]
        _t += [("Cdlp", i) for i in _opts]
        _brutepardict.update({"Cdlp": {"index": _tl, "options": _opts}})
        _tl += 1

    if "R3" in params_model:
        _opts = [5, 10, 50, 150, 500, 5e3, 5e5]
        _t += [("R3", i) for i in _opts]
        _brutepardict.update({"R3": {"index": _tl, "options": _opts}})
        _tl += 1

    if "Q3" in params_model:
        _opts = [5e-5, 1e-3, 2e-2]
        _t += [("Q3", i) for i in _opts]
        _brutepardict.update({"Q3": {"index": _tl, "options": _opts}})
        _tl += 1
    if "n3" in params_model:
        _opts = [0.4, 0.6, 0.9]
        _t += [("n3", i) for i in _opts]
        _brutepardict.update({"n3": {"index": _tl, "options": _opts}})
        _tl += 1
    #    _tset = set([i[0] for i in _t])
    #       random_combination( combinations(_t,_tl),5)
    selected_combis = random_combo_select(_brutepardict, number_guess)
    return selected_combis


def make_prefit_frame(
    EIS_data_KKvalid,
    lmfitting,
    prefix="pp",
    plot=False,
    check_err=True,
    norm=np.array([]),
):
    #            norm = 1/(Z_KKv.real**2+Z_KKv.imag**2)
    #            abs(Z_KKv)
    #            norm = 1
    #            norm = np.sqrt(EIS_data_KKvalid.DATA_weightsmod_Z.values)
    #            norm = 1/abs(Z_KKv)
    #            norm = 1/np.sqrt(EIS_data_KKvalid.DATA_weightsmod_Z.values)
    # make_prefit_frame(EIS_data_KKvalid, out, plot = 'Y')
    if np.array([]).any() == False:
        np.array([1] * len(lmfitting.best_fit.real))
    Z_KKv = EIS_data_KKvalid.DATA_Z.values
    if norm.size == 0:
        norm = (Z_KKv.real ** 2 + Z_KKv.imag ** 2) ** -1

    resIm, resRe = lmfitting.residual[1::2], lmfitting.residual[0::2]
    pp_errRE, pp_errIM = (lmfitting.best_fit.real - Z_KKv.real) * norm, (
        lmfitting.best_fit.imag - Z_KKv.imag
    ) * norm
    pp_errRE_mean, pp_errIM_mean = np.abs(pp_errRE).mean(), np.abs(pp_errIM).mean()
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
    if plot:
        fig, ax = plt.subplots(3, 1, figsize=(4, 8))
        #                if 'Y' in str(plot_pp):
        prefit_data.plot(x=f"DATA_Yre", y=f"DATA_Yim", kind="scatter", ax=ax[0])
        prefit_data.plot(x=f"{prefix}_Y_Re", y=f"{prefix}_Y_Im", c="r", ax=ax[0])
        #                else:
        prefit_data.plot(x=f"DATA_Zre", y=f"DATA_-Zim", kind="scatter", ax=ax[1])
        prefit_data.plot(x=f"{prefix}_Z_Re", y=f"{prefix}_Z_Im", c="r", ax=ax[1])

        prefit_data.plot(x="Frequency(Hz)", y=f"{prefix}_err_Re", c="g", ax=ax[2])
        prefit_data.plot(x="Frequency(Hz)", y=f"{prefix}_err_Im", c="k", ax=ax[2])

        box1 = TextArea(lmfitting.fit_report(), textprops=dict(color="k"))

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
        hf_data = prefit_data.loc[prefit_data["Frequency(Hz)"] < 500]
        #                    , prefit_data.loc[prefit_data['Frequency(Hz)'] > freqlim]
        lf_errRe_mean, lf_errIm_mean = (
            lf_data[f"{prefix}_err_Re"].abs().mean(),
            lf_data[f"{prefix}_err_Im"].abs().mean(),
        )

        hf_errRe_mean, hf_errIm_mean = (
            hf_data[f"{prefix}_err_Re"].abs().mean(),
            hf_data[f"{prefix}_err_Im"].abs().mean(),
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
        #
        #            if all([hf_errRe_mean > n_std*pp_errRE_std + pp_errRE_mean, hf_errIm_mean > n_std*pp_errIM_std + pp_errIM_mean]):
        #                if test_out[1] == False:
        #                    test_out[1] = n_std
        ##                =test_out + f'bad fit {name} (lim {freqlim} Hz)\n'
        #            n_std += 0.01
        test_out = [lf_ratio_Re + lf_ratio_Im, hf_ratio_Im + hf_ratio_Re]
        good_fit_test = True if any(i < 0.5 for i in test_out) else False
        return good_fit_test, test_out


#        mod_combs.update({modname : {'tset' : _tset, 'combs' : combinations(_t,_tl)}})
#        [(i,[].append(a[1])) for i in _tset for a in _t if i in a[0]]


#        for k,val in _brutepardict.items():
#            random.sample(val['options'],1)
#
#        brutepar_lst = [i for i in combinations(_t,_tl) if set([a[0] for a in i]) == _tset]
#        selected_combis = []
#        for i in :
#            if set([a[0] for a in i]) == _tset:
#                if 'Rct' and 'Rorr' in _tset:
#                    if i[_brutepardict.get('Rct')][1] < i[_brutepardict.get('Rorr')][1]:
#                        selected_combis.append(i)
##                    brutepar_lst = [i for i in brutepar_lst if i[_brutepardict.get('Rct')][1] < i[_brutepardict.get('Rorr')][1]]
#                elif 'Qad' and 'Cdlp' in _tset:
#                    if i[_brutepardict.get('Qad')][1] > i[_brutepardict.get('Cdlp')][1]:
#                        selected_combis.append(i)
##                    brutepar_lst = [i for i in brutepar_lst if i[_brutepardict.get('Qad')][1] > i[_brutepardict.get('Cdlp')][1]]
#                        if 'Q3' in _tset:
#                            if i[_brutepardict.get('Q3')][1] < i[_brutepardict.get('Cdlp')][1]:
#                                selected_combis.append(i)
#                else:
#                    selected_combis.append(i)
#                        brutepar_lst = [i for i in brutepar_lst if i[_brutepardict.get('Q3')][1] < i[_brutepardict.get('Cdlp')][1]]


def fit_EEC(fit_run_arg, **EIS_fit_kwargs):
    #%%
    EISgr_data_EV = fit_run_arg.data
    gr_EIS_ovv = fit_run_arg.ovv
    PF_EIS_dest_dir = EISgr_data_EV.EIS_dest_dir.unique()[0]
    #    it_run_arg.ovv.Dest_dir.iloc[0]
    indexes_per_EV = []
    EvRHE = "E_AppV_RHE"
    # FIXME preparing meta dicts
    EISgr_data_meta = (
        EISgr_data_EV[
            [i for i in EISgr_data_EV.columns if EISgr_data_EV[i].nunique() == 1]
        ]
        .iloc[0]
        .to_dict()
    )
    gr_EIS_ovv_meta = gr_EIS_ovv.iloc[0].to_dict()
    EISgr_meta_combined = {**gr_EIS_ovv_meta, **EISgr_data_meta}
    overlap_keys = [i for i in EISgr_data_meta.keys() if i in gr_EIS_ovv_meta.keys()]
    matching_keys = [
        i for i in overlap_keys if EISgr_data_meta[i] == gr_EIS_ovv_meta[i]
    ]
    unmatching_keys = [i for i in overlap_keys if i not in matching_keys]
    # Important variables defined for shortened
    #    EIS_kwargs.get(')

    if len(unmatching_keys) == 0:
        pass
    elif len(unmatching_keys) == 1 and unmatching_keys[0] == "PAR_file":
        pass
        if (len(EISgr_data_meta) + len(gr_EIS_ovv_meta)) - (
            len(EISgr_meta_combined) + len(matching_keys) + len(unmatching_keys)
        ) != 0:
            logger.error(
                "EIS missing keys in EISgr_data_meta and ovv_meta dicts, {:s} at {:G}".format(
                    fit_run_arg.PAR_file.stem, fit_run_arg.E_dc_RHE
                )
            )
    else:
        logger.error(
            "EIS non-matching keys in EISgr_data_meta and ovv_meta dicts, {:s} at {:G}".format(
                fit_run_arg.PAR_file.stem, fit_run_arg.E_dc_RHE
            )
        )
    # PAR_file = Path(gr_EIS_ovv['PAR_file'].unique()[0]) # FALSE GROUP, contains different PAR files !!!!

    combined_meta_slice = [
        "postAST",
        "Electrolyte",
        "pH",
        "Gas",
        "RPM",
        "Loading_cm2",
        "Loading_name",
    ] + [
        "SampleID",
        "Measured_OCP",
        "Instrument",
        "Instrument_SN",
        "Electrode",
        "Scanrate",
        "RPM_DAC",
        "Type_action",
    ]

    EISgr_meta_add = {}
    if all([(i in EISgr_meta_combined) for i in combined_meta_slice]):
        # EISgr_meta_combined =  EISgr_data_meta[['pH','RPM','Electrolyte','SampleID','Gas','Measured_OCP','Instrument']].to_dict()
        EISgr_meta_add.update({k: EISgr_meta_combined[k] for k in combined_meta_slice})
    else:
        [i for i in combined_meta_slice if i not in EISgr_meta_combined]
        EISgr_meta_add.update(
            {
                k: EISgr_meta_combined[k]
                for k in combined_meta_slice
                if k in EISgr_meta_combined
            }
        )
        logger.error(
            "EIS outP missing columns in EISgr_data_meta, {:s} at {:G}  mV".format(
                fit_run_arg.PAR_file.stem, fit_run_arg.E_dc_RHE_mV
            )
        )

    EIS_outPath = PF_EIS_dest_dir.joinpath(
        Path(
            str(fit_run_arg.PAR_file.stem)
            + "_{:.0f}mV_{:.0f}rpm".format(fit_run_arg.E_dc_RHE_mV, fit_run_arg.RPM_DAC)
        )
    )
    #    EIS_outPath_OLD = EIS_kwargs.get('EIS_dest_dir').joinpath(Path(str(Path(EISgr_data_meta['PAR_file']).stem)+'_{:.0f}mV'.format(fit_run_arg.E_dc_RHE_mV)))
    EIS_outPath.parent.mkdir(parents=True, exist_ok=True)
    EIS_outpath_linKK_target = FileOperations.CompareHashDFexport(
        pd.DataFrame(), Path(str(EIS_outPath) + "_linkK")
    ).with_suffix(".png")
    #    EISgr = EISgr.loc[EISgr["Frequency(Hz)"] < FreqLim,:]
    freq = EISgr_data_EV["Frequency(Hz)"].values
    ang = freq * 2 * np.pi
    Zre, Zim = EISgr_data_EV["Z Real"].values, EISgr_data_EV["Z Imag"].values
    Zdata = Zre + 1j * Zim
    Ydata = Zdata ** -1
    Yre, Yim = Ydata.real, Ydata.imag

    DataWeights_modulus = (Zre ** 2 + Zim ** 2) ** -1
    #    DataWeights_modulus_Y = (Yre**2+Yim**2)
    DataWeights_unitary = Zre / Zre

    EIS_data_raw = pd.DataFrame(
        {
            "Frequency(Hz)": freq,
            "Angular": ang,
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
            EvRHE: fit_run_arg.E_dc_RHE,
            "DATA_weightsmod_Z": DataWeights_modulus,
            "PAR_file": fit_run_arg.PAR_file,
            "RPM_DAC": fit_run_arg.RPM_DAC,
        },
        index=EISgr_data_EV.index,
    )

    # ====Important Frequency limit slice on EIS data!!!==== #
    FreqLim = EIS_fit_kwargs.get("FreqLim", 30e3)
    EIS_data_raw.loc[
        (EIS_data_raw["Frequency(Hz)"] >= FreqLim) | (EIS_data_raw["DATA_Zim"] > 0),
        "Valid",
    ] = False
    EIS_data_freqlim = EIS_data_raw.query("Valid == True")
    # ====!!!==== #
    """ include here KK-test!! """
    EIS_data_KKvalid, linkKK_invalid_prefit, linKK_pars = get_KKvalid(
        EIS_data_freqlim, restype_set="Z", res_scale=1, **EIS_fit_kwargs
    )
    EIS_data_raw.loc[
        ~EIS_data_raw["Frequency(Hz)"].isin(EIS_data_KKvalid["Frequency(Hz)"]), "Valid"
    ] = False
    linKK_pars.update(
        {
            "Unvalid_freqs": ", ".join(
                EIS_data_raw.loc[
                    ~EIS_data_raw["Frequency(Hz)"].isin(
                        EIS_data_KKvalid["Frequency(Hz)"]
                    ),
                    "Frequency(Hz)",
                ]
                .astype(str)
                .to_list()
            )
        }
    )
    plot_linKK(
        EIS_data_KKvalid,
        EIS_data_raw,
        meta=EISgr_data_meta,
        save_target=EIS_outpath_linKK_target,
        plot_prefit=True,
        linkKK_invalid_prefit=linkKK_invalid_prefit,
        **linKK_pars,
    )
    #    restype_set = 'Z'

    #    M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(EIS_data_freqlim['Frequency(Hz)'].values,EIS_data_freqlim['DATA_Z'].values,type_res=restype_set)
    ##    plot_linKK(EIS_data['Frequency(Hz)'].values,EIS_data.DATA_Z.values,Z_linKK,res_real,res_imag,meta=EISgr_data_meta,linKK=[M,mu],type_res=restype_set)
    #    # plot 1st KK validation for testing purposes
    #    prefit_suffix = '_prefit'
    #    EIS_data_linKK = KramersKronigValidation.assigncols(EIS_data_freqlim,Z_linKK,res_real,res_imag,suffix=prefit_suffix)
    ##    linKK = pd.DataFrame({'linKK_Z' : Z_linKK, 'linKK_Zreal' : Z_linKK.real,'linKK_Zimag' : Z_linKK.imag,
    ##                          'linKK_Y' : Z_linKK**-1, 'linKK_Yreal' : (Z_linKK**-1).real,'linKK_Yimag' : (Z_linKK**-1).imag,
    ##                          'linKK_resRe' : res_real, 'linKK_resIm' : res_imag }, index=EIS_data_freqlim.index)
    ##    EIS_data_linKK = pd.concat([EIS_data_freqlim,linKK],axis=1)
    ##    linKK_limit = 0.0475
    #    linKK_trimming_std_factor = EIS_fit_kwargs.get('linKK_trimming_factor',1.550)
    #    linKK_limit_Re = np.abs(EIS_data_linKK['linKK_resRe'+prefit_suffix]).mean()+linKK_trimming_std_factor*np.abs(EIS_data_linKK['linKK_resRe'+prefit_suffix]).std()
    #    linKK_limit_Im = np.abs(EIS_data_linKK['linKK_resIm'+prefit_suffix]).mean()+linKK_trimming_std_factor*np.abs(EIS_data_linKK['linKK_resIm'+prefit_suffix]).std()
    #    KKvalid = EIS_data_linKK.loc[((((np.abs(EIS_data_linKK['linKK_resRe'+prefit_suffix]) < linKK_limit_Re) & (np.abs(EIS_data_linKK['linKK_resIm'+prefit_suffix]) < linKK_limit_Im)) & (EIS_data_linKK['Frequency(Hz)'] > 2)) | (EIS_data_linKK['Frequency(Hz)'] <= 2))]
    #    EIS_data.loc[~EIS_data['Frequency(Hz)'].isin(KKvalid['Frequency(Hz)']),'Valid'] = False
    #    M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(KKvalid['Frequency(Hz)'].values,KKvalid['DATA_Z'].values,type_res=restype_set,max_M=50)
    #
    #    EIS_data_valid = KramersKronigValidation.assigncols(KKvalid,Z_linKK,res_real,res_imag)
    ##    EIS_data_valid_KK = KramersKronigValidation.assigncols(EIS_data,Z_linKK,res_real,res_imag)
    #    EIS_data.loc[~EIS_data['Frequency(Hz)'].isin(KKvalid['Frequency(Hz)']),'Valid'] = False
    #    outData.loc[outData['Frequency(Hz)'].isin([i for i in outData['Frequency(Hz)'].values if i not in KKvalid['Frequency(Hz)'].values])]
    #    EIS_data.loc[~EIS_data.index.isin(EIS_data_valid.index),'Valid'] = False

    #        index_raw = {'PAR_file'  : nm2,'DestFile' : Parsout_path_target,'Type_output' : 'EIS_Pars', 'Type_exp' : 'EIS','nModels' : Pars.Model_EEC.nunique(), 'Models_used' : Pars.Model_EEC.unique()}

    #    M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(KKvalid['Frequency(Hz)'].values,KKvalid['DATA_Z'].values,type_res=restype_set,max_M=2)
    #    outData = pd.DataFrame(data=EIS_data,index=EISgr.index)
    #    outData = pd.DataFrame(data=EIS_data,index=EISgr.index)
    #    'DATA_Zmod' : [cmath.phase(i) for i in Zdata],'DATA_Zmod' : [cmath.polar(i)[0] for i in Zdata][cmath.polar(i)[1]*180/np.pi for i in Zdata],
    #%%
    # if True:
    #    EIS_data_KKvalid = EIS_data.query('Valid == True')

    Z_KKv, ang_KKv = EIS_data_KKvalid.DATA_Z.values, EIS_data_KKvalid.Angular.values
    freq_KKv = EIS_data_KKvalid["Frequency(Hz)"].values
    statiscial_weights = EIS_data_KKvalid.DATA_weightsmod_Z.values

    pre_prefit_linkKK_data = "DATA_Z"
    if "linKK_Z" in EIS_data_KKvalid.columns:
        if (
            np.sum(
                EIS_data_KKvalid.linKK_resRe ** 2 + EIS_data_KKvalid.linKK_resIm ** 2
            )
            < 5e-3
        ):
            pre_prefit_linkKK_data = "linKK_Z"

    # create a set of Parameters

    Rs_guess = EIS_data_KKvalid.loc[
        np.abs(EIS_data_KKvalid.DATA_Zmod).idxmin(), "DATA_Zre"
    ]
    if Rs_guess > 5 and Rs_guess < 200:
        pass
    else:
        Rs_guess = EIS_data_KKvalid.DATA_Zre.min()

    standard_init_params = Parameters()

    standard_init_params.add(
        "Rs", value=Rs_guess, min=0.5 * Rs_guess, max=1.5 * Rs_guess
    )
    standard_init_params.add("Cdlp", value=0.003, min=1e-08, max=5e-2)
    standard_init_params.add("nDL", value=0.90, min=0.3, max=1)
    standard_init_params.add(
        "Rct", value=10, min=1e-3, max=1e5, vary=True, brute_step=1000
    )
    #    params.add('Aw', value= 3000, min=5.0, max=1E6)
    #    params.add('Rad', value= 5E5, min=0.01, max=1E7)
    standard_init_params.add("Qad", value=1e-02, min=1e-04, max=0.5)
    standard_init_params.add("nAd", value=0.6, min=0.2, max=1)
    standard_init_params.add("Rorr", value=3e3, min=0.01, max=1e11, brute_step=1000)
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    standard_init_params.add_many(
        ("R3", 100, True, 0.01, 1e8, None, None),
        ("Q3", 1e-04, True, 1e-06, 0.1, None, None),
        ("n3", 0.6, True, 0.2, 1, None, None),
    )

    standard_init_params = params_extra_setting(standard_init_params, EISgr_data_EV)
    # Old style:  mod,mod2 = Model(EEC_Singh2015_RQRQR,name='Singh2015_RQRQR'),Model(EEC_ORRpOx,name='Bandarenka_2011_RQRQR')
    #   skip(5,Model(EEC_Bandarenka_Ads,name='Bandarenka2011_RQRW'))
    mod_lst_filter = EEC_models_index(exclude=["(Singh2015_RQR)"])

    # TODO check whats wrong with Bandarenka Fit model

    MTHDS = [
        "leastsq",
        "nelder",
        "lbfgsb",
        "anneal",
        "powell",
        "cobyla",
        "slsqp",
        "differential-evolution",
        "ampgo",
        "basinhopping",
        "dual_annealing",
        "brute",
    ]
    pars_lst, fit_data_valid_lst, outP = [], [], {}
    for mod_indx, model_set in mod_lst_filter:  #
        modname = model_set.name
        params_model = model_set.make_params()
        lstsq_method = MTHDS[0]
        #            print('use previous')
        #            params = InitP1.query('Model_EEC')
        for pn in params_model:
            params_model.add(standard_init_params[pn])

        try:
            _get_params = pd.DataFrame()
            _get_params, recheck_msg = fitting_recheck_params(
                fit_run_arg, modname, params_model, **EIS_fit_kwargs
            )
            #            logger.warning(recheck_msg)
            if not _get_params.empty:
                for pn in params_model:
                    if (
                        params_model[pn].vary == True
                        and pd.isna(_get_params[pn].iloc[0]) == False
                    ):
                        params_model[pn].set(value=_get_params[pn].iloc[0])
                logger.info(f"""Prefit used recheked params""")
                run_brute_prefit = False
            else:
                logger.warning(
                    f'Prefit _get_params empty {(str(fit_run_arg.PAR_file), fit_run_arg.E_dc_RHE, fit_run_arg.RPM_DAC)}\n message:"{recheck_msg}'
                )
                run_brute_prefit = True
        except Exception as e:
            run_brute_prefit = True
            logger.error(f"fitting_recheck_params error: {e}")

        if "nAd" in params_model:
            if modname in ["Model(Singh2015_RQRWR)", "Model(Bandarenka2011_RQRW)"]:
                params_model["nAd"].set(value=0.5, vary=False)
            else:
                params_model["nAd"].set(value=0.9, vary=True)
        # params_model.pretty_print()
        #    Weights = np.array([i/i for i in range(1,len(ang)+1)][::-1])
        #    result = minimize(model_ORR, params, args=(ang,Zdata),method=methods[0])
        #    result2 = minimize(model_ORRpOx, params, args=(ang,Zdata),method=methods[0])
        #    z == z.real + z.imag*1j
        #    if Initp1 is not {} or Initp2 and not {}:
        #        params = InitP1

        #        errDF = pd.DataFrame()
        #        pre_init1,pre_init2 =  mod.eval(pre_prefit1.params, ang=ang),mod2.eval(pre_prefit2.params, ang=ang)
        #        pre_prefit1 , pre_prefit2

        #        if EIS_fit_kwargs.get('perform_prefit',False) == True:

        #            print('Prefit on model "{:s}" with method {:s}' .format(modname, prefit_method))
        #            for method in (MTHDS[0],MTHDS[-5]):
        # TODO            # PrePrefit on Z_linKK !! check
        ### Making preprefit on linKK Z values NOT on Zdata

        pre_prefit_leastsq = model_set.fit(
            EIS_data_KKvalid[pre_prefit_linkKK_data],
            params_model,
            ang=ang_KKv,
            weights=statiscial_weights,
            method=MTHDS[0],
        )
        pre_prefit_DE = model_set.fit(
            EIS_data_KKvalid[pre_prefit_linkKK_data],
            params_model,
            ang=ang_KKv,
            weights=statiscial_weights,
            method=MTHDS[-5],
        )

        logger.info(
            f"Prefit test method: {MTHDS[0]} ({pre_prefit_leastsq.aic:.0f}), {MTHDS[-5]} ({pre_prefit_DE.aic:.0f})"
        )
        #        pre_prefit_leastsq.redchi,pre_prefit_DE.redchi
        if pre_prefit_leastsq.aic < pre_prefit_DE.aic:
            pre_prefit = pre_prefit_leastsq
            #            make_prefit_frame(EIS_data_KKvalid, pre_prefit_leastsq,plot = True, norm= statiscial_weights)
            succes_msg = "Fit succeeded"
            logger.info(f"Best prefit method: {MTHDS[0]}")
        else:
            logger.info(f"Best prefit method: {MTHDS[-5]}")
            #            make_prefit_frame(EIS_data_KKvalid, pre_prefit_DE,plot = True, norm= statiscial_weights)
            pre_prefit = pre_prefit_DE
            succes_msg = "successfully"
        best_trial = model_set.fit(
            Z_KKv,
            pre_prefit.params,
            ang=ang_KKv,
            weights=statiscial_weights,
            method=pre_prefit.method,
        )
        make_prefit_frame(EIS_data_KKvalid, best_trial, norm=statiscial_weights)
        #        make_prefit_frame(EIS_data_KKvalid, pre_prefit_DE, prefix = 'pp',plot = True, norm= statiscial_weights)
        pp_good_low_high, best_test_msg = make_prefit_frame(
            EIS_data_KKvalid, best_trial, norm=statiscial_weights, check_err=True
        )
        run_brute_prefit = (
            False
            if all(
                [
                    succes_msg in best_trial.message,
                    best_trial.redchi < 5e-06,
                    best_trial.aic < -1000,
                    sum(best_test_msg) < 1.6,
                ]
            )
            else True
        )
        if (
            all([i in params_model for i in ["Rct", "Rorr"]])
            and run_brute_prefit == False
        ):
            run_brute_prefit = (
                False
                if all([best_trial.params["Rct"] < best_trial.params["Rorr"]])
                else True
            )
            run_brute_prefit = True if all([best_trial.params["Rct"] > 5000]) else False

        #        run_brute_prefit, test_msg = make_prefit_frame(EIS_data_KKvalid, best_trial, prefix = 'pp', check_err = True)
        brute_prefit_method = MTHDS[-5]
        #        run_brute_prefit = True
        if run_brute_prefit == True:
            brutepar_lst = prepare_brute_pars(
                params_model, EIS_fit_kwargs.get("EISfit_brute_random_n", 90)
            )
            logger.info(
                f"""Prefit brute run starting with {brute_prefit_method}({len(brutepar_lst)} tests) for {fit_run_arg.PAR_file.stem}
             at {fit_run_arg.E_dc_RHE:0.2G} on model {model_set.name},RedChiSqr = {best_trial.redchi:.3G},{best_trial.aic:.4G}"""
            )
            br_tests = []
            for bpar in brutepar_lst:
                for i in bpar:
                    best_trial.params[i[0]].set(
                        value=i[1], vary=True, min=0.1 * i[1], max=10 * i[1]
                    )
                trial_prefit_brute = model_set.fit(
                    Z_KKv,
                    best_trial.params,
                    ang=ang_KKv,
                    weights=statiscial_weights,
                    method=brute_prefit_method,
                )
                #                    print(f'Test Rct{bpar}, redchi ({trial_prefit_brute.redchi:.4G}), aic {trial_prefit_brute.aic:.4G}')
                bf_good_low_high, brute_test_msg = make_prefit_frame(
                    EIS_data_KKvalid, trial_prefit_brute, norm=statiscial_weights
                )
                _br_tests = {}
                _br_tests = {
                    **best_trial.params.valuesdict(),
                    **{
                        "aic": trial_prefit_brute.aic,
                        "redchi": trial_prefit_brute.redchi,
                        "test_low": brute_test_msg[0],
                        "test_high": brute_test_msg[1],
                    },
                }
                br_tests.append(_br_tests)
                if (
                    trial_prefit_brute.aic < best_trial.aic
                    and sum(brute_test_msg) < sum(best_test_msg)
                    and brute_test_msg[0] < best_test_msg[0]
                    and brute_test_msg[1] < best_test_msg[1]
                ):
                    best_trial = trial_prefit_brute
                    best_test_msg = brute_test_msg
            #                    make_prefit_frame(EIS_data_KKvalid, best_trial, norm= statiscial_weights,plot=1)
            #            BR_test_ovv = pd.DataFrame(br_tests)
            #            BR_test_ovv.plot(x='aic' ,y='Rct',kind='scatter',logy=True)
            #            BR_test_ovv.plot(x='aic' ,y='Rs',kind='scatter',logy=True)
            #            BR_test_ovv.plot(x='aic' ,y='Rorr',kind='scatter',logy=True)
            #            BR_test_ovv.plot(x='aic' ,y='Cdlp',kind='scatter')
            #            BR_test_ovv.plot(x='aic' ,y='Qad',kind='scatter')

            for i in params_model.keys():
                best_trial.params[i].set(
                    min=standard_init_params[i].min, max=standard_init_params[i].max
                )
            best_trial = model_set.fit(
                Z_KKv,
                best_trial.params,
                ang=ang_KKv,
                weights=statiscial_weights,
                method="leastsq",
            )
            #            make_prefit_frame(EIS_data_KKvalid, best_trial, prefix = 'pp',plot = 'Y')
            logger.info(
                f"""Prefit brute finished with {best_trial.method} for {fit_run_arg.PAR_file.stem}
             at {fit_run_arg.E_dc_RHE:0.2G} V on model {model_set.name},RedChiSqr = {best_trial.redchi:.3G}, {pre_prefit.aic:.4G}"""
            )

        prefit = model_set.fit(
            Z_KKv,
            best_trial.params,
            ang=ang_KKv,
            weights=statiscial_weights,
            method="leastsq",
        )
        #            make_prefit_frame(EIS_data_KKvalid, prefit, prefix = 'pp',plot = 'Y')
        init = model_set.eval(prefit.params, ang=ang_KKv)
        logger.info(
            f"""Prefit finished with method {prefit.method} for {fit_run_arg.PAR_file.stem}
                    at {fit_run_arg.E_dc_RHE:0.2G} V on model {model_set.name},RedChiSqr = {prefit.redchi:.3G} {pre_prefit.aic:.4G}"""
        )
        # ADD RANDOM TO PRE-OUT params
        InitParam = prefit.params
        for par in InitParam.keys():
            if InitParam[par].vary:
                _parval = InitParam[par].value
                _rpar = random.gauss(_parval, 0.2 * _parval)
                while _rpar > InitParam[par].max and _rpar < InitParam[par].min:
                    _rpar = random.gauss(_parval, 0.2 * _parval)
                InitParam[par].set(_rpar)
        #        InitParam = prefit.params

        errRE_init, errIM_init = (init.real - Z_KKv.real) / abs(Z_KKv), (
            init.imag - Z_KKv.imag
        ) / abs(Z_KKv)
        #        Weights = 1/abs(Z_KKv)
        if modname in ["Model(Singh2015_RQRWR)", "Model(Bandarenka2011_RQRW)"]:
            InitParam["nAd"].set(value=0.5, vary=False)

        for i in InitParam.keys():
            InitParam[i].set(
                min=standard_init_params[i].min, max=standard_init_params[i].max
            )

        out = model_set.fit(
            Z_KKv,
            InitParam,
            ang=ang_KKv,
            weights=statiscial_weights,
            method=lstsq_method,
        )
        #        make_prefit_frame(EIS_data_KKvalid, out, prefix = 'pp',plot = 'Y')
        bf_good_low_high, (out_low_err, out_high_err) = make_prefit_frame(
            EIS_data_KKvalid, out
        )
        logger.info(
            f"""Fit out with method {out.method} for {fit_run_arg.PAR_file.stem} at {fit_run_arg.E_dc_RHE:0.2G} on model \
                    {model_set.name},ChiSqr = {out.chisqr:.3G}, RedChiSqr = {out.redchi:.3G}, {out.aic:.4G}"""
        )
        ### == RETRY SECTION: when first fit is not good enough.. some initial values are re-set and fit is run again ... ####
        retry_attempt = "no"

        #                print(pre_prefit.fit_report())
        #                print(best_trial.fit_report())
        #                print(trial_prefit_brute.fit_report())
        #                    for Rct_guess in [5,5E2,1E3,10E3,100E3]:
        #                        best_trial.params['Rct'].set(value= Rct_guess,vary=True, min= 0.01*Rct_guess, max=20*Rct_guess)
        #                        if 'Rorr' in params_model:
        #                            for Rorr_guess in [500,1E3,3E3,1E4,3E4,1E5,5E8][::-1]:
        #                                best_trial.params['Rorr'].set(value= Rorr_guess, min= 0.01*Rorr_guess, max=20*Rorr_guess)
        #
        #                                if 'R3' in params_model:
        #                                    for R3_guess in [500,1E3,1E4,1E6][::-1]:
        #                                        if Rct_guess < Rorr_guess:
        #                                            best_trial.params['R3'].set(value= R3_guess, min= 0.01*R3_guess, max=20*R3_guess)
        #                                            trial_prefit_brute = model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv,weights= statiscial_weights, method= prefit_method)
        #                                            if trial_prefit_brute.redchi < best_trial.redchi:
        #                                                best_trial = trial_prefit_brute
        #                                    best_trial.params['R3'].set(min= 1, max=1E10)
        #                                else:
        #                                    if Rct_guess < Rorr_guess:
        #                                        trial_prefit_brute = model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv,weights= statiscial_weights, method= prefit_method)
        ##                                        print(f'Test Rct{Rct_guess}, Rorr{Rorr_guess}: redchi ({trial_prefit_brute.redchi})')
        #                                        if trial_prefit_brute.redchi < best_trial.redchi:
        #                                            best_trial = trial_prefit_brute
        #                            best_trial.params['Rorr'].set(min=1, max=1E10)
        #    #                                print('best trial updated: R_ct[{:.3G}] -> {:.3G}  & R_orr[{:.3G}] -> {:.3G}, improvement of {:.2%} ({:.2E})'.format(Rct_guess, trial_prefit_brute.params['Rct'].value,
        #    #                                      Rorr_guess, trial_prefit_brute.params['Rorr'].value,  1E-2*(best_trial.chisqr-trial_prefit_brute.chisqr) /best_trial.chisqr, trial_prefit_brute.redchi))
        #    #                                trial_prefit_brute.params.pretty_print()
        #                        else:
        #                            trial_prefit_brute = model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv,weights= statiscial_weights, method= prefit_method)
        #                            if trial_prefit_brute.redchi < best_trial.redchi:
        #    #                            print('best trial updated: Rct[{:.3G}] -> {:.3G} improvement of {:.2%} ({:.2E})'.format(Rct_guess,trial_prefit_brute.params['Rct'].value, 1E-2*(best_trial.chisqr-trial_prefit_brute.chisqr) /best_trial.chisqr, trial_prefit_brute.redchi))
        #                                best_trial = trial_prefit_brute
        #                    best_trial.params['Rct'].set(min=1E-2, max=1E6)
        #                            trial_prefit_brute.params.pretty_print()
        #            print(pre_prefit.fit_report())
        #            else:

        #        WeightsRes = ((DataRes2['index']*DataRes2[0]).values[1::2]**2 + (DataRes2['index']*DataRes2[0]).values[0::2]**2)
        #        WeightsRes =   abs(errRE_init+1j*errIM_init)
        #        Weights = WeightsRes
        #        DataRes2[(DataRes2 < Rmaxlim) & (DataRes2 > Rminlim)]
        #        DateResFiltered = DataRes2.loc[(DataRes2[0] < Rmaxlim) & (DataRes2[0] > Rminlim)]
        #            InitParam = prefit.params
        #            if EIS_fit_kwargs.get('TrimData',False) == True:
        #                errRE_init,errIM_init = (init.real-Z_KKv.real)/abs(Z_KKv), (init.imag-Z_KKv.imag)/abs(Z_KKv)
        #                errRE_initSTD,errIM_initSTD =  errRE_init.std(),errIM_init.std()
        #                errRE_initMEAN,errIM_initMEAN =  errRE_init.mean(),errIM_init.mean()
        #                trimSTDmax = 1.75
        #                errRmaxlim,errRminlim = (errRE_initMEAN+trimSTDmax*np.abs(errRE_initSTD)),(errRE_initMEAN-trimSTDmax*np.abs(errRE_initSTD))
        #                errRmaxlimIm,errRminlimIm = (errIM_initMEAN+trimSTDmax*np.abs(errIM_initSTD)),(errIM_initMEAN-trimSTDmax*np.abs(errIM_initSTD))
        #        #        pd.DataFrame([freq,prefit2.residual])
        #        #        abs(errRE+1j*errIM)
        #                errDF = pd.DataFrame([errRE_init,errIM_init,freq]).T
        #        #        errDF.loc[~((errDF[0] < errRmaxlim) & (errDF[0] > errRminlim)),0] = errRE_initMEAN
        #        #        errDF.loc[~((errDF[1] < errRmaxlimIm) & (errDF[1] > errRminlimIm)),1] = errIM_initMEAN
        #        #        badFreq = errDF.loc[~((errDF[1] < errRmaxlimIm) & (errDF[1] > errRminlimIm)),2]
        #                badFreq = errDF.loc[~(((errDF[0] < errRmaxlim) & (errDF[0] > errRminlim)) & ((errDF[1] < errRmaxlimIm) & (errDF[1] > errRminlimIm))),2]
        #        #        | ((errDF[0] < errRmaxlim) & (errDF[0] > errRminlim))),2]
        #                TrimmedOutData = EIS_data_KKvalid.loc[~EIS_data_KKvalid['Frequency(Hz)'].isin(badFreq.values)]
        #                trimAng = TrimmedOutData['Frequency(Hz)'].values*2*np.pi
        #                TrimmedWeightsResOutlier = abs(errDF.loc[~errDF[2].isin(badFreq.values),0]+errDF.loc[~errDF[2].isin(badFreq.values),1]*1j)
        #                WeightsResOutlier = abs(errDF[0]+errDF[1]*1j)
        #                prefit_res_mean,prefit_res_std = prefit.residual.mean(),prefit.residual.std()
        #                Rmaxlim,Rminlim = (prefit_res_mean+trimSTDmax*np.abs(prefit_res_std)),(prefit_res_mean-trimSTDmax*np.abs(prefit_res_std))
        #                prefit.residual[(prefit.residual < Rmaxlim) & (prefit.residual > Rminlim)]
        #
        #                DataRes = pd.DataFrame(prefit.residual,Zdata.view(np.float)).reset_index()
        #                DataRes.loc[~((DataRes[0] < Rmaxlim) & (DataRes[0] > Rminlim)),0] = 1E-12
        #
        #                trim_ang = TrimmedOutData['Frequency(Hz)'].values*2*np.pi
        #    #            Zre, Zim = EISgr_data_EV['Z Real'].values, EISgr_data_EV['Z Imag'].values
        #                Z_TrimmedData = TrimmedOutData['DATA_Z'].values
        #
        #    #            Ydata = 1/Z_TrimmedData
        #    #            WeightsRes = TrimmedWeightsResOutlier
        #                EISgr_trimmed_out = EISgr_data_EV.loc[~EISgr_data_EV['Frequency(Hz)'].isin(TrimmedOutData['Frequency(Hz)'].values)]
        #                EISgr_data_EV = EISgr_data_EV.loc[EISgr_data_EV['Frequency(Hz)'].isin(TrimmedOutData['Frequency(Hz)'].values)]
        #
        #                prefit = model_set.fit(Z_TrimmedData,pre_prefit.params,ang=trim_ang,weights= 1/abs(Z_TrimmedData), method= prefit_method )
        ##                prefit2 = mod2.fit(Z_TrimmedData,pre_prefit2.params,ang=trim_ang,weights= 1/abs(Z_TrimmedData), method='differential-evolution')
        #                init =  model_set.eval(prefit.params, ang=trim_ang)
        #                InitParam = prefit.params
        #                EIS_Trimming_plot(EISgr_data_EV,gr_EIS_ovv,EIS_dest_dir,outData,TrimmedOutData,fit_run_arg.E_dc_RHE)
        #    #            WeightsRes =   abs(errRE_init+1j*errIM_init)
        #    #            1/(Zre**2+1/Zim**2)
        #            if EIS_fit_kwargs.get('FitOnlyTrimmedData',False) == True:
        #                ang,Zdata = trim_ang,Z_TrimmedData
        #                outData = TrimmedOutData
        #        else:
        #            logger.warning('No prefit for {:s} at {:0.2f} with model {:s}'.format(PAR_file.stem,fit_run_arg.E_dc_RHE,modname))
        #            init = model_set.eval(params_model, ang=ang_KKv)
        #            InitParam   = params_model

        #        out1 = mod.fit(Z_TrimmedData,InitParam1,ang=trim_ang,weights= 1/Z_TrimmedData, method=MTHDS[0])
        #        out2 = mod2.fit(Z_TrimmedData,InitParam2,ang=trim_ang,weights=1/Z_TrimmedData, method=MTHDS[0])
        ##       out2_initW = mod2.fit(Zdata,InitParam2,ang=ang,weights=DataWeights, method=MTHDS[-2])
        #        fit1,fit2 = out1.eval(ang=trim_ang),out2.eval(ang=trim_ang)
        ##=== EIS REFIT USING RESIDUAL ===#
        #        errRE,errIM = (fit2.real-Z_TrimmedData.real)/abs(Z_TrimmedData), (fit2.imag-Z_TrimmedData.imag)/abs(Z_TrimmedData)
        #        errRE1,errIM1 = (fit1.real-Z_TrimmedData.real)/abs(Z_TrimmedData), (fit1.imag-Z_TrimmedData.imag)/abs(Z_TrimmedData)
        #        outData =

        #        if out.chisqr > 1E-4:
        #            large_chisqr = out.chisqr
        #            retry_method = MTHDS[0]
        #            logger.warning('Restarting Fit >>> Chisqr larger than 0.07 (={:0.3G}) for model {:s}.\n  fit with method {:s} for {:s} at {:0.2G}'.format(large_chisqr, modname, retry_method, PAR_file.stem, E_dc_RHE))
        #            retry_params = out.params
        #            if 'O2' in gr_EIS_ovv['Gas'].unique()[0]:
        #                if 'Rct' in retry_params.valuesdict().keys():
        #                    retry_params['Rct'].value = 20
        #                if 'nAd' in retry_params.valuesdict().keys():
        #                    retry_params['nAd'].value = 0.77
        #                if 'Rorr' in retry_params.valuesdict().keys():
        #                    retry_params['Rorr'].value = 1000
        #            if 'N2' in gr_EIS_ovv['Gas'].unique()[0] and EISgr_data_EV.pH.unique()[0] < 7:
        #                if 'Rct' in retry_params.valuesdict().keys():
        #                    retry_params['Rct'].value = 1
        #                if 'nAd' in retry_params.valuesdict().keys():
        #                    retry_params['nAd'].value = 1
        #                if 'Rorr' in retry_params.valuesdict().keys():
        #                    retry_params['Rorr'].value = 5E5
        #            if 'nAd' in retry_params:
        #                if modname in ['Model(Singh2015_RQRWR)','Model(Bandarenka2011_RQRW)']:
        #                    retry_params['nAd'].set(value= 0.5, vary=False)
        #                else:
        #                    retry_params['nAd'].set(value= 0.8, vary=True)
        #            retry_out = model_set.fit(Z_KKv,retry_params,ang=ang_KKv,weights= statiscial_weights , method = retry_method)
        #            retry_chisqr_diff = np.abs(large_chisqr - retry_out.chisqr)
        #            retry_attempt = 'tried'
        #            if retry_out.chisqr < large_chisqr:
        #                out = retry_out
        #                retry_attempt = 'tried and result'
        #                logger.warning('Results of Re-fit Chisqr {:0.3f}.\n used method {:s} for {:s} at {:G}'.format(retry_chisqr_diff,MTHDS[0],PAR_file.stem,E_dc_RHE))
        #            else:
        #                pass
        fit = out.eval(ang=ang_KKv)
        #    out2_initW = mod2.fit(Zdata,InitParam2,ang=ang,weights=DataWeights, method=MTHDS[-2])
        # === EIS REFIT USING RESIDUAL ===#
        errRE, errIM = (fit.real - Z_KKv.real) / abs(Z_KKv), (
            fit.imag - Z_KKv.imag
        ) / abs(Z_KKv)
        #    print('Prefit Weights: DataWeights)
        # 'Init_Z1' : init1,'Init_Z2' : init2, 'FIT1_Z' : fit1,
        EIS_data_KKvalid_fit = pd.DataFrame()
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
                "errRe": errRE,
                "errIm": errIM,
                "Model_EEC": modname,
                "Model_index": mod_indx,
                "PAR_file": fit_run_arg.PAR_file,
            }
        )
        # add the metadata to the EIS_data_KKvalid_fit DF !!
        EISgr_meta_add_to_fit = pd.DataFrame(
            [EISgr_meta_add] * len(EIS_data_KKvalid_fit),
            index=EIS_data_KKvalid_fit.index,
        )
        #        EISgr_meta_add_to_fit.drop()
        EIS_data_KKvalid_fit = pd.concat(
            [EIS_data_KKvalid_fit, EISgr_meta_add_to_fit], axis=1
        )
        # EIS_data_KKvalid_fit = EIS_data_KKvalid_fit.assign(**[{i[0] : [i[1]]*len(EIS_data_KKvalid_fit)} for i in EISgr_meta_add.items()])
        # EIS_data_KKvalid_fit = EIS_data_KKvalid_fit.set_index(EIS_set_index_columns())
        #    metadata to OutData
        #'FIT2_Zmod' : [cmath.polar(i)[0] for i in fit2], FIT2_Zphase'[cmath.polar(i)[1]*180/np.pi for i in fit2]
        #    print('%s\n'%out.message,out.fit_report(min_correl=0.50))
        #    print('%s\n'%out2.message,out2.fit_report(min_correl=0.90)) EvRHE : [EISgr[EvRHE].unique()[0]]*len(outData)
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

        outP.update(EISgr_meta_add)  # from beginning of function
        outP.update(linKK_pars)  # from beginning of function after linKK validation
        outP.update(
            {
                EvRHE: fit_run_arg.E_dc_RHE,
                "PAR_file": fit_run_arg.PAR_file,
                "PAR_date": EISgr_meta_combined["PAR_date"],
                "RPM_DAC": fit_run_arg.RPM_DAC,
            }
        )
        #        outP.update(extraP)
        if "Qad" not in outP.keys():
            outP.update({"Qad": 0, "nAd": None})
        xtra = {
            "Rct_kin": outP["Rct"] ** -1,
            "Qad+Cdlp": outP["Qad"] + outP["Cdlp"],
            "Chisqr": out.chisqr,
            "RedChisqr": out.redchi,
            "aic": out.aic,
            "bic": out.bic,
            "Model_EEC": modname,
            "Model_index": mod_indx,
            "lmfit_method": out.method,
            "lmfit_message": out.message,
            "lmfit_out": out,
            "retry_attempt": retry_attempt,
            "lmfit_var_names": ", ".join(out.var_names),
            "lmfit_errcheck_msg": bf_good_low_high,
            "lmfit_errcheck_low": out_low_err,
            "lmfit_errcheck_high": out_high_err,
        }
        outP.update(xtra)
        pars_lst.append(pd.DataFrame(outP, index=[EIS_data_KKvalid_fit.index[0]]))
        fit_data_valid_lst.append(EIS_data_KKvalid_fit)
    #        lmfit_out_lst.append({'Model_EEC' : modname,'Model_index' : mod_indx,'lmfit_out' : out})
    pars_models = pd.concat(pars_lst, sort=False, ignore_index=True)
    EIS_fit_data = pd.concat(fit_data_valid_lst, sort=False, ignore_index=True)
    #%%
    #    EIS_fit_data.groupby('Model_index').plot(x=['DATA_Yre','FIT_Yre'],y=['DATA_Yim','FIT_Yim'],kind='scatter')
    # +++ PLOTTING of EIS spectrum +++
    vers = FileOperations.version
    spectra_fit_outpath = EIS_outPath.with_name(
        EIS_outPath.stem + f"_spectrumfit_v{vers}"
    ).with_suffix(".xlsx")
    spectra_raw_outpath = EIS_outPath.with_name(
        EIS_outPath.stem + f"_spectrumraw_v{vers}"
    ).with_suffix(".xlsx")
    pars_outpath = EIS_outPath.with_name(
        EIS_outPath.stem + f"_pars_v{vers}"
    ).with_suffix(".xlsx")

    EIS_fit_data = EIS_fit_data.assign(
        **{"File_Pars": pars_outpath, "File_SpecRaw": spectra_raw_outpath}
    )
    EIS_data_raw = EIS_data_raw.assign(
        **{"File_Pars": pars_outpath, "File_SpecFit": spectra_fit_outpath}
    )

    spectra_fit_outpath_target = FileOperations.CompareHashDFexport(
        EIS_data_raw, spectra_fit_outpath
    )
    spectra_raw_outpath_target = FileOperations.CompareHashDFexport(
        EIS_fit_data, spectra_raw_outpath
    )

    pars_models = pars_models.assign(
        **{
            "File_SpecFit": spectra_fit_outpath_target,
            "File_SpecRaw": spectra_raw_outpath_target,
        }
    )
    pars_outpath_target = FileOperations.CompareHashDFexport(pars_models, pars_outpath)

    if "Plot" in EIS_fit_kwargs.get("EIS_single_output", "Text, Plot"):
        spectra_fit_outpath_png = spectra_fit_outpath.with_suffix(".png").with_suffix(
            ".png"
        )
        EIS_plotting_per_EV(
            EIS_fit_data, pars_models, spectra_fit_outpath_png, plot_show=False
        )
    #        index_dataplot_output = {'PAR_file': nm2,'Type_output' : 'EIS_fit_data_png', 'Type_exp' : 'EIS', 'DestFile' : EIS_outPath_target_png, 'E_V' : fit_run_arg.E_dc_RHE,'RPM_DAC' : fit_run_arg.RPM_DAC}
    #        indexes_per_EV_out.append(index_dataplot_output)
    #    index_data_output = {'PAR_file': fit_run_arg.PAR_file,'Type_output' : 'EIS_fit_data', 'Type_exp' : 'EIS', 'DestFile' : EIS_fit_data_outPath_target, 'E_V' : fit_run_arg.E_dc_RHE,'RPM_DAC' : fit_run_arg.RPM_DAC}
    #    indexes_per_EV.append(index_data_output)
    #    if EIS_fit_kwargs.get('export_raw_data',False):

    #        index_raw_output = {'PAR_file': fit_run_arg.PAR_file,'Type_output' : 'EIS_raw_data', 'Type_exp' : 'EIS', 'DestFile' : EIS_outpath_target_raw,'E_V' : fit_run_arg.E_dc_RHE}
    #        indexes_per_EV.append(index_raw_output)
    meta_export_templ = namedtuple(
        "meta_export",
        fit_run_arg._fields + ("File_SpecFit", "File_SpecRaw", "File_Pars"),
    )
    meta_export = meta_export_templ(
        *fit_run_arg,
        spectra_fit_outpath_target,
        spectra_raw_outpath_target,
        pars_outpath_target,
    )
    fit_export_templ = namedtuple("fit_export", "fit_spectra fit_pars meta_index")
    fit_export = fit_export_templ(EIS_fit_data, pars_models, meta_export)
    #    indexes_per_EV_out = pd.DataFrame(indexes_per_EV)
    return fit_export
