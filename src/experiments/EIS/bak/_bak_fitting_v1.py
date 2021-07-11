# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:00:16 2020

@author: DWXMG
"""

import sys
from collections import namedtuple, OrderedDict
from pathlib import Path
import itertools
from itertools import combinations
import operator

import random
from math import pi
import copy
import pickle
import logging
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

from lmfit import Parameters, conf_interval, Minimizer, minimize

# import numdifftools
import corner


print(f"Name: {__name__} for file {__file__}")

if __name__ == "__main__":
    sys.path.append(str(Path(__file__).parent.parent.parent.joinpath("runEC")))
    sys.path.append(str(Path(__file__).parent.parent.parent))
    from validation import get_KKvalid
    from models import EEC_models_index, params_extra_setting

    #    from .plotting import plot_linKK, EIS_Trimming_plot
    #    import ECpy
    from ECpy.experiments.EIS.plotting import (
        plot_linKK,
        EIS_Trimming_plot,
        EIS_plotting_per_EV,
        plot_lin_Warburg,
    )
    from ECpy.experiments.EIS.eis_run_ovv import EIS_Spectrum
    from DRT_DP_fitting import DP_DRT_analysis
    from GP_DRT_fitting import run_GP_DRT_fit

    #    sys.path.append(Path(__file__).parent.parent.parent.parent)
    from EC_logging_config import start_logging

#    from EC_logging_config import start_logging
#    sys.path.append(Path(__file__).parent.parent.joinpath('runEC'))
#    from runEC.EC_logging_config import start_logging
#    _logger = start_logging(__name__)
elif __name__ in "fitting":
    from validation import get_KKvalid, prep_GP_DRT_raw_data
    from models import EEC_models_index, params_extra_setting
    from ECpy.experiments.EIS.plotting import (
        plot_linKK,
        EIS_Trimming_plot,
        EIS_plotting_per_EV,
        plot_lin_Warburg,
    )
    from ECpy.experiments.EIS.eis_run_ovv import EIS_Spectrum
    from DRT_DP_fitting import DP_DRT_analysis
    from GP_DRT_fitting import run_GP_DRT_fit

    FH_path = Path(__file__).parent.parent.parent.parent.joinpath("FileHelper")
    sys.path.append(str(FH_path))
    from FileHelper.FileFunctions import FileOperations

    sys.path.append(Path(__file__).parent.parent.parent)
    from ECpy.runEC.EC_logging_config import start_logging

#    from EC_logging_config import start_logging
#    sys.path.append(Path(__file__).parent.parent.joinpath('runEC'))
#    from runEC.EC_logging_config import start_logging
#    _logger = start_logging(__name__)


else:
    print(__name__)
    from .validation import get_KKvalid, prep_GP_DRT_raw_data
    from .models import EEC_models_index, params_extra_setting
    from .plotting import plot_linKK, EIS_Trimming_plot, EIS_plotting_per_EV
    from .DRT_DP_fitting import DP_DRT_analysis
    from .GP_DRT_fitting import run_GP_DRT_fit

    #    from .plotting import EIS_plotting_per_EV, EIS_Trimming_plot, EIS_plotting_EvRHE

    FH_path = Path(__file__).parent.parent.parent.parent.joinpath("FileHelper")
    sys.path.append(str(FH_path))
    from FileHelper.FileFunctions import FileOperations

    sys.path.append(Path(__file__).parent.parent.joinpath("runEC"))
    from ECpy.runEC.EC_logging_config import start_logging
#    _logger = start_logging(__name__)


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


#    if all([i in params_model for i in ['Rct','Rorr']]):
#        if params_model['Rct'] < params_model['Rorr']:
#            _Rorr = params_model['Rct']
#            params_model['Rct'].set(value=params_model['Rorr'])
#            params_model['Rorr'].set(value=_Rorr )
# TODO finish brute preparation
def prepare_rand_models(total_combis=5e3, rand_params_dest_dir=""):
    mod_lst_filter = EEC_models_index(exclude=["(Singh2015_RQR)"])

    _hash = hash(str(mod_lst_filter))
    dump_pkl = False
    _prep_params = {}
    if rand_params_dest_dir:
        _pkl_file = Path(rand_params_dest_dir).joinpath(
            f"randon_model_pars_{_hash}.pkl"
        )
        if _pkl_file.is_file():
            with open(_pkl_file, "rb") as handle:
                _prep_params = pickle.load(handle)
        else:
            dump_pkl = True

    if not _prep_params:
        _pp = []
        for i, mod, _ in mod_lst_filter:
            modname = mod.name
            #        print(modname)
            params_model = mod.make_params()
            if "nAd" in params_model:
                if modname == "Model(Bandarenka2011_RQRW)":
                    params_model["nAd"].set(value=0.5, vary=False)

                elif modname == "Model(Singh2015_RQRWR)":
                    params_model["nDL"].set(value=0.5, vary=False)
                    params_model["nAd"].set(value=0.9, vary=True)
            #        else:
            #            params_model['nDL'].set(value= 0.5, vary=True)
            #            params_model['nAd'].set(value= 0.9, vary=True)
            _pp.append((modname, prepare_brute_pars(params_model, total_combis)))
        _prep_params = dict(_pp)
        if dump_pkl:
            with open(_pkl_file, "wb") as handle:
                pickle.dump(_prep_params, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return _prep_params


#     brutepar_lst_big = prepare_brute_pars(params_model,20E3)


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
    _default = {"index": 0}
    number_guesses = np.min(
        [
            0.1 * np.prod([len(k["options"]) for k in _brutepardict.values()]),
            number_guesses,
        ]
    )
    _counter = 0
    while len(_guesses) < number_guesses and _counter < 3 * number_guesses:
        _combo = [
            (key, random.sample(val["options"], 1)[0])
            for key, val in _brutepardict.items()
        ]

        if _combo in _guesses:
            continue
        if "Rct" and "Rorr" in _brutepardict.keys():
            _Rct, _Rorr = (
                _combo[_brutepardict.get("Rct", _default)["index"]][1],
                _combo[_brutepardict.get("Rorr", _default)["index"]][1],
            )
            if "Qad" and "Cdlp" in _brutepardict.keys():
                _Cdlp, _Qad = (
                    _combo[_brutepardict.get("Cdlp", _default)["index"]][1],
                    _combo[_brutepardict.get("Qad", _default)["index"]][1],
                )
                ##                 if 'Q3' and 'R3' in _brutepardict.keys():
                ##                     _R3, _Q3 =  _combo[_brutepardict.get('R3',_default)['index']][1], _combo[_brutepardict.get('3',_default)['index']][1]
                #                 if _Rct < _Rorr and _Cdlp < _Qad and _Rorr < _R3 and _Qad < _Q3:
                #                     _guesses.append(_combo)
                #                 else:
                if _Rct < _Rorr and _Cdlp < _Qad:
                    _guesses.append(_combo)
            else:
                if _Rct < _Rorr:
                    _guesses.append(_combo)
        else:
            _guesses.append(_combo)
        _counter += 1
    #        if not _counter % 1000:
    #            print(_counter)
    #                    brutepar_lst = [i for i in brutepar_lst if i[_brutepardict.get('Rct')][1] < i[_brutepardict.get('Rorr')][1]]
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
        _opts = [1, 2, 5, 8, 10, 30, 50, 2e2, 1e3, 3e3]
        _t += [("Rct", i) for i in _opts]
        _brutepardict.update({"Rct": {"index": _tl, "options": _opts}})
        _tl += 1
    if "Rorr" in params_model:
        _opts = [i * 3 + 1 for i in [1, 2, 5, 8, 10, 30, 50, 2e2, 5e2, 1e3, 3e3, 5e4]]
        #        _opts = [10, 100, 500,750,1E3,1.5E3,2E3,3.5E3,5E3]
        _t += [("Rorr", i) for i in _opts]
        _brutepardict.update({"Rorr": {"index": _tl, "options": _opts}})
        _tl += 1

    if "nAd" in params_model:
        _opts = [0.4, 0.5, 0.51, 0.6, 0.9][::-1]
        _t += [("nAd", i) for i in _opts]
        _brutepardict.update({"nAd": {"index": _tl, "options": _opts}})
        _tl += 1
    if "nDL" in params_model:
        _opts = [0.5, 0.6, 0.9, 1][::-1]
        _t += [("nDL", i) for i in _opts]
        _brutepardict.update({"nDL": {"index": _tl, "options": _opts}})
        _tl += 1

    if "Qad" in params_model:
        _opts = [1e-5, 3.5e-05, 1e-4, 1e-3, 2e-2]
        #        _opts = [1E-3,5E-3,25E-3, 70E-3]
        _t += [("Qad", i) for i in _opts]
        _brutepardict.update({"Qad": {"index": _tl, "options": _opts}})
        _tl += 1

    if "Cad" in params_model:
        _opts = [1e-5, 3.5e-05, 1e-4, 1e-3, 2e-2]
        #        _opts = [1E-3,5E-3,25E-3, 70E-3]
        _t += [("Cad", i) for i in _opts]
        _brutepardict.update({"Cad": {"index": _tl, "options": _opts}})
        _tl += 1

    if "Cdlp" in params_model:
        _opts = [i * 0.7 for i in [1e-6, 5e-06, 1e-5, 7e-5, 7e-04]]
        _t += [("Cdlp", i) for i in _opts]
        _brutepardict.update({"Cdlp": {"index": _tl, "options": _opts}})
        _tl += 1

    if "R3" in params_model:
        _opts = np.logspace(0, 4.5, 7)
        #        [i*2+1 for i in [1, 2, 5, 8, 10, 30, 50,2E2,5E2, 1E3,3E3,5E4]]
        #        [5,10,17,50,150, 500, 1600, 2.5E3, 5E3,5E5]
        _t += [("R3", i) for i in _opts]
        _brutepardict.update({"R3": {"index": _tl, "options": _opts}})
        _tl += 1

    if "Q3" in params_model:
        _opts = [1e-5, 1e-4, 25e-3]
        _t += [("Q3", i) for i in _opts]
        _brutepardict.update({"Q3": {"index": _tl, "options": _opts}})
        _tl += 1
    if "n3" in params_model:
        _opts = [0.4, 0.6, 0.75, 1]
        _t += [("n3", i) for i in _opts]
        _brutepardict.update({"n3": {"index": _tl, "options": _opts}})
        _tl += 1

    if "Aw" in params_model:
        _opts = list(np.logspace(0.1, 3.5, 7))
        _t += [("Aw", i) for i in _opts]
        _brutepardict.update({"Aw": {"index": _tl, "options": _opts}})
        _tl += 1
    if "sigmaW" in params_model:
        _opts = list(np.logspace(-2, 3, 6))
        _t += [("sigmaW", i) for i in _opts]
        _brutepardict.update({"sigmaW": {"index": _tl, "options": _opts}})
        _tl += 1

    #    _tset = set([i[0] for i in _t])

    for key in _brutepardict.keys():
        if not params_model[key].vary:
            _brutepardict[key]["options"] = [params_model[key].value]
        if params_model[key].min:
            _brutepardict[key]["options"] = [
                i
                if i > params_model[key].min
                else np.min(_brutepardict[key]["options"])
                for i in _brutepardict[key]["options"]
            ]
        if params_model[key].max:
            _brutepardict[key]["options"] = [
                i
                if i < params_model[key].max
                else np.max(_brutepardict[key]["options"])
                for i in _brutepardict[key]["options"]
            ]
    #       random_combination( combinations(_t,_tl),5)
    selected_combis = random_combo_select(_brutepardict, number_guess)
    return selected_combis


def run_brute_fit():

    (
        best_weights,
        best_MSE,
        best_higlow,
        best_sum,
        _bmethod,
        _b_pretest,
        best_trial_weight,
    ) = sorted(wt_check, key=lambda MSE: MSE[1])[0]
    best_MSE_norm = best_MSE * weight_opts["lmfit_weights_mod_Y"].mean()
    statiscial_weights = weight_opts[best_weights]
    best_trial = model_set.fit(
        Z_KKv,
        best_trial_weight.params,
        ang=ang_KKv,
        weights=statiscial_weights,
        method="leastsq",
    )
    _logger.info(
        f"""Prefit test method: {MTHDS[0]} ({best_trial.aic:.0f}), {MTHDS[-5]} ({best_trial.aic:.0f}), took : {best_trial.method}.
    best_weight ({best_weights}), MSE({best_MSE:.3F}) """
    )
    #        make_prefit_frame(EIS_data_KKvalid, best_trial, prefix = 'pp',plot = True, norm= statiscial_weights) #TODO
    pp_good_low_high, best_test_msg, MSE = make_prefit_frame(
        EIS_data_KKvalid, best_trial, norm=statiscial_weights, check_err=True
    )
    run_brute_prefit = (
        False if all([best_trial.success, best_MSE_norm < 50e-03]) else True
    )
    if all([i in params_model for i in ["Rct", "Rorr"]]) and run_brute_prefit == False:
        run_brute_prefit = (
            False
            if all([best_trial.params["Rct"].value < best_trial.params["Rorr"].value])
            else True
        )
        run_brute_prefit = (
            True if all([best_trial.params["Rct"].value > 5000]) else False
        )
    if (
        run_brute_prefit == True
        and "leastsq" in best_trial.method
        and best_trial.redchi < 2e-01
    ):
        run_brute_prefit = False

    #        run_brute_prefit, test_msg = make_prefit_frame(EIS_data_KKvalid, best_trial, prefix = 'pp', check_err = True)
    brute_prefit_method = MTHDS[-5]
    run_brute_prefit = False  # TODO OVERWRITE BRUTE FIT
    BR_best_fit_weight_opt_out, EIS_data_KKvalid_BR_specs = (
        pd.DataFrame(),
        pd.DataFrame(),
    )

    #%%
    _brutepars_lookup = []
    _brutepars_lookup = EIS_fit_kwargs.get("random_params", {}).get(modname, [])
    good_pars_grp = EIS_fit_kwargs.get("good_fit_pars", {})
    if good_pars_grp:
        good_pars_mod = good_pars_grp.get(modname, [])
        _brutepars_lookup += good_pars_mod
    if wt_check:
        _brutepars_lookup += [
            list(
                zip(
                    i[-1].params.valuesdict().keys(), i[-1].params.valuesdict().values()
                )
            )
            for i in wt_check
        ]

    if _brutepars_lookup:
        brutepar_lst_big = _brutepars_lookup
    else:
        brutepar_lst_big = prepare_brute_pars(params_model, 10e3)
    #            brutepar_lst = prepare_brute_pars(params_model, np.min([rand_n,EIS_fit_kwargs.get('EISfit_brute_random_n',150)]))
    brutepar_lst_big.append(
        [
            (k, val * (1 + np.random.normal(0, 0.1)))
            for k, val in best_trial.params.valuesdict().items()
        ]
    )

    _logger.info(
        f"""Prefit brute run starting with {brute_prefit_method}({len(brutepar_lst_big)} tests) for {fit_run_arg.PAR_file.stem}
     at {fit_run_arg.E_dc_RHE:0.2G} on model {model_set.name},RedChiSqr = {best_trial.redchi:.3G},{best_trial.aic:.4G}"""
    )
    br_eval_test = []
    #            brutepar_lst_big = prepare_brute_pars(params_model,20E3)
    for bpar_eval in brutepar_lst_big:
        #                brute_prefit_method = random.sample((MTHDS[-5],MTHDS[0]),1)[0]
        list(filter(lambda x: x[0] in standard_init_params.keys(), bpar_eval))
        for pname, pval in bpar_eval:
            if standard_init_params[pname].vary:
                best_trial.params[pname].set(
                    value=pval,
                    vary=standard_init_params[pname].vary,
                    min=np.max([standard_init_params[pname].min, 0.1 * pval]),
                    max=np.min([standard_init_params[pname].max, 5 * pval]),
                )
            else:
                best_trial.params[pname].set(
                    value=standard_init_params[pname].value,
                    vary=standard_init_params[pname].vary,
                )
        best_trial.params["Rs"].set(
            value=Rs_guess, vary=True, min=0.8 * Rs_guess, max=1.2 * Rs_guess
        )
        #                trial_prefit_brute = model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv,weights= statiscial_weights, method= brute_prefit_method)
        trial_prefit_brute_fit = model_set.eval(best_trial.params, ang=ang_KKv)
        MSE_re = (trial_prefit_brute_fit.real - Z_KKv.real) ** 2
        MSE_im = (trial_prefit_brute_fit.imag - Z_KKv.imag) ** 2
        MSE = sum(MSE_im + MSE_re)
        MSE_high = sum(MSE_im[-25:] + MSE_re[-25:])
        br_eval_test.append((MSE, MSE_high, bpar_eval))
    br_eval_test_all = sorted(br_eval_test, key=lambda MSE: MSE[0])
    br_eval_test = (
        br_eval_test_all[0:5]
        + sorted(br_eval_test_all, key=lambda MSEhf: MSEhf[1])[0:3]
    )
    br_tests = []
    for br_eval_MSE, br_MSE_hf, bpar in br_eval_test + random.sample(
        br_eval_test_all, 2
    ):
        #                brute_prefit_method = random.sample((MTHDS[-5],MTHDS[0]),1)[0]
        #                br_wt_opt = random.sample(weight_opts.keys(),1)[0]
        for br_wt_opt in weight_opts.keys():
            try:
                for pkey, pval in bpar:
                    if standard_init_params[pkey].vary:
                        best_trial.params[pkey].set(
                            value=pval,
                            vary=standard_init_params[pkey].vary,
                            min=np.max([standard_init_params[pkey].min, 0.1 * pval]),
                            max=np.min([standard_init_params[pkey].max, 5 * pval]),
                        )
                    else:
                        best_trial.params[pkey].set(
                            value=standard_init_params[pkey].value,
                            vary=standard_init_params[pkey].vary,
                        )
                best_trial.params["Rs"].set(
                    value=Rs_guess, vary=True, min=0.8 * Rs_guess, max=1.2 * Rs_guess
                )
                trial_prefit_brute = model_set.fit(
                    Z_KKv,
                    best_trial.params,
                    ang=ang_KKv,
                    weights=weight_opts[br_wt_opt],
                    method="leastsq",
                )

                #                trial_prefit_brute_fit = model_set.eval(best_trial.params,ang=ang_KKv,weights= statiscial_weights, method= brute_prefit_method)
                #                MSE_re = (trial_prefit_brute_fit.real-Z_KKv.real)**2
                #                MSE_im = (trial_prefit_brute_fit.imag-Z_KKv.imag)**2
                #                MSE = sum(MSE_im + MSE_re)
                #                    print(f'Test Rct{bpar}, redchi ({trial_prefit_brute.redchi:.4G}), aic {trial_prefit_brute.aic:.4G}')
                bf_good_low_high, brute_test_msg, br_fit_MSE = make_prefit_frame(
                    EIS_data_KKvalid, trial_prefit_brute, norm=statiscial_weights
                )
                _br_tests = {}
                _br_tests = {
                    **trial_prefit_brute.params.valuesdict(),
                    **{
                        "lmfit_aic": trial_prefit_brute.aic,
                        "lmfit_redchi": trial_prefit_brute.redchi,
                        "test_low": brute_test_msg[0],
                        "test_high": brute_test_msg[1],
                        "test_sum": sum(brute_test_msg),
                        "br_fit_MSE": br_fit_MSE,
                        "br_fit_weight": br_wt_opt,
                        "br_eval_MSE": br_eval_MSE,
                        "lmfit_var_names": ", ".join(trial_prefit_brute.var_names),
                    },
                }
                tr_bf_params_stderss = [
                    (i + "_stderr", trial_prefit_brute.params.__getitem__(i).stderr)
                    for i in trial_prefit_brute.params
                ]
                _br_tests.update(
                    dict(
                        zip(
                            [i[0] for i in tr_bf_params_stderss],
                            [i[1] for i in tr_bf_params_stderss],
                        )
                    )
                )

                br_tests.append(_br_tests)
                if (
                    trial_prefit_brute.aic < best_trial.aic
                    and sum(brute_test_msg) < sum(best_test_msg)
                    and brute_test_msg[0] < best_test_msg[0]
                    and brute_test_msg[1] < best_test_msg[1]
                ):
                    best_trial = trial_prefit_brute
                    best_test_msg = brute_test_msg
            except Exception as ebr:
                _logger.error(
                    f"""Prefit brute run for {fit_run_arg.PAR_file.stem} at {fit_run_arg.E_dc_RHE:0.2G} on model {model_set.name}.
                                error: {ebr}"""
                )

    #                    make_prefit_frame(EIS_data_KKvalid, trial_prefit_brute, norm= statiscial_weights,plot=1)
    BR_test_ovv = pd.DataFrame(br_tests).sort_values(by="br_fit_MSE")
    #            BR_test_ovv['test_sum'] = BR_test_ovv['test_low']+BR_test_ovv['test_high']
    #            BR_test_ovv.sort_values(by='test_sum')
    _br_gr = []
    for n, gr in BR_test_ovv.groupby("br_fit_weight"):
        BR_best_gr = gr.loc[
            (gr.lmfit_aic < gr.lmfit_aic.mean() - 0.5 * gr.lmfit_aic.std())
            & (gr.lmfit_redchi < 0.5 * gr.lmfit_redchi.mean())
            & (gr.test_sum < gr.test_sum.mean())
        ].sort_values(by="br_fit_MSE")
        if BR_best_gr.empty:
            BR_best_gr = gr.loc[
                (gr.lmfit_aic < gr.lmfit_aic.mean())
                & (gr.lmfit_redchi < 0.5 * gr.lmfit_redchi.mean())
                & (gr.test_sum < gr.test_sum.mean())
            ].sort_values(by="br_fit_MSE")

        _br_gr.append(BR_best_gr)
    BR_best = pd.concat(_br_gr).sort_values(by="br_fit_MSE")
    BR_best_fit = pd.concat(
        [
            BR_test_ovv.sort_values("test_high").head(10),
            BR_test_ovv.sort_values("test_low").head(10),
            BR_best,
        ]
    )
    BR_best_fit.drop_duplicates(
        subset=["test_high", "br_fit_weight", "br_fit_MSE"], inplace=True
    )
    #             BR_best = BR_test_ovv.loc[(BR_test_ovv.aic < BR_test_ovv.aic.mean()-0.5*BR_test_ovv.aic.std()) & (BR_test_ovv.redchi < 0.5*BR_test_ovv.redchi.mean())
    #                                        & (BR_test_ovv.test_sum < BR_test_ovv.test_sum.mean())].sort_values(by='redchi')
    #            pd.DataFrame(_BR_best_out.values(),index=_BR_best_out.keys()) =
    BR_best_fit_weight_opt = BR_best_fit.loc[
        [gr["br_fit_MSE"].idxmin() for n, gr in BR_best_fit.groupby("br_fit_weight")]
    ]
    BR_best_fit_weight_opt["BRUTE_FIT"] = 1
    BR_best_fit_weight_opt["FINAL_FIT"] = 0
    #            BR_best_5 = BR_best.head(3)
    _BR_best_out, _BR_specs = {}, []
    EIS_data_KKvalid_BR_specs = EIS_data_KKvalid.copy()
    for n, r in BR_best_fit_weight_opt.iterrows():
        for p in params_model:
            if standard_init_params[p].vary:
                _min = np.min([standard_init_params[p].value, 0.95 * BR_best[p].min()])
                _max = np.max([standard_init_params[p].value, 1.05 * BR_best[p].max()])
                if _min == _max:
                    _min = 0.9 * _min
                best_trial.params[p].set(
                    value=r[p], vary=standard_init_params[p].vary, min=_min, max=_max
                )
            else:
                best_trial.params[p].set(
                    value=standard_init_params[p].value,
                    vary=standard_init_params[p].vary,
                )
        #                    best_trial.params[p].set(value= r[p],vary=standard_init_params[p].vary,
        #                                     min=standard_init_params[p].min, max = standard_init_params[p].max)
        post_brute = model_set.fit(
            Z_KKv,
            best_trial.params,
            ang=ang_KKv,
            weights=weight_opts[r["br_fit_weight"]],
            method="leastsq",
        )
        spec_prefix = "br_" + r["br_fit_weight"]
        post_br_frame = make_prefit_frame(
            EIS_data_KKvalid, post_brute, get_frame=1, prefix=spec_prefix
        )
        merge_cols = list(post_br_frame.columns.intersection(EIS_data_KKvalid.columns))
        EIS_data_KKvalid_BR_specs = pd.merge(
            EIS_data_KKvalid_BR_specs,
            post_br_frame,
            on=merge_cols,
            how="inner",
            left_index=True,
            right_index=True,
        )
        #                list(post_br_frame.columns.difference(EIS_data_KKvalid.columns))
        #                _BR_specs.append(post_br_frame)
        #                BR_best_spec_out post_br_frame

        #                EIS_data_KKvalid.DATA_weightsmod_Z.values
        pb_ow_high, pb_test_msg, MSE = make_prefit_frame(
            EIS_data_KKvalid, post_brute, norm=weight_opts[r["br_fit_weight"]]
        )
        _BR_best_out.update(
            {
                n: {
                    "brute_fit_method": post_brute.method,
                    "brute_fit_obj": post_brute,
                    "brute_spec_prefix": spec_prefix,
                }
            }
        )
    #                make_prefit_frame(EIS_data_KKvalid, post_brute, norm= statiscial_weights,plot=1)
    BR_best_fit_weight_opt_out = pd.concat(
        [
            BR_best_fit_weight_opt,
            pd.DataFrame(_BR_best_out.values(), index=_BR_best_out.keys()),
        ],
        axis=1,
    )
    #            BR_best_spec_out = pd.concat(_BR_specs,axis=1,join='inner')
    #            BR_best_spec_out.columns = sorted(list(set(list(BR_best_spec_out))))
    BR_best_wt_idx = BR_best_fit_weight_opt_out.br_fit_MSE.idxmin()
    BR_best_trial = BR_best_fit_weight_opt_out.loc[BR_best_wt_idx, "brute_fit_obj"]
    BR_best_wtname = BR_best_fit_weight_opt_out.loc[BR_best_wt_idx, "br_fit_weight"]
    BR_best_wtMSE = BR_best_fit_weight_opt_out.loc[BR_best_wt_idx, "br_fit_MSE"]

    if best_weights != BR_best_wtname:
        _logger.info(
            f"""Prefit brute finished with changing ({best_weights}) to {BR_best_wtname} for {fit_run_arg.PAR_file.stem}
     at {fit_run_arg.E_dc_RHE:0.2G} V on model {model_set.name},RedChiSqr = {BR_best_trial.redchi:.3G}, {BR_best_trial.aic:.4G}"""
        )
        best_weights, best_MSE = (BR_best_wtname, BR_best_wtMSE)
        statiscial_weights = weight_opts[best_weights]
    #            sorted(_BR_best_out, key=lambda fit: fit[0])[0][1]
    #            make_prefit_frame(EIS_data_KKvalid, best_trial, norm= statiscial_weights,plot=1)
    #            BR_test_ovv.plot(x='lmfit_aic' ,y='Rct',kind='scatter',logy=True)
    #            BR_test_ovv.plot(x='aic' ,y='Rs',kind='scatter',logy=True)
    #            BR_test_ovv.plot(x='aic' ,y='Rorr',kind='scatter',logy=True)
    #            BR_test_ovv.plot(x='aic' ,y='Cdlp',kind='scatter')
    #            BR_test_ovv.plot(x='aic' ,y='nDL',kind='scatter')
    #            BR_test_ovv.plot(x='aic' ,y='Qad',kind='scatter')
    for i in params_model.keys():
        if best_trial.params[i].vary:
            best_trial.params[i].set(
                value=BR_best_trial.params[i].value * (1 + np.random.normal(0, 0.05)),
                min=standard_init_params[i].min,
                max=standard_init_params[i].max,
                vary=standard_init_params[i].vary,
            )
        else:
            best_trial.params[i].set(value=standard_init_params[i].value)
    best_trial = model_set.fit(
        Z_KKv,
        best_trial.params,
        ang=ang_KKv,
        weights=statiscial_weights,
        method="leastsq",
    )
    #            EIS_data_KKvalid.DATA_weightsmod_Z.values**-1
    #            make_prefit_frame(EIS_data_KKvalid, best_trial, prefix = 'pp',plot = 'Y') # TODO
    _logger.info(
        f"""Prefit brute finished with {best_trial.method} for {fit_run_arg.PAR_file.stem}
     at {fit_run_arg.E_dc_RHE:0.2G} V on model {model_set.name},RedChiSqr = {best_trial.redchi:.3G}, {pre_prefit.aic:.4G}"""
    )


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
    # make_prefit_frame(EIS_data_KKvalid, out, plot = 'Y')
    if np.array([]).any() == False:
        norm = np.array([1] * len(lmfitting.best_fit.real))
    Z_KKv = EIS_data_KKvalid.DATA_Z.values
    if norm.size == 0:
        norm = Z_KKv.real / Z_KKv.real  # (Z_KKv.real**2+Z_KKv.imag**2)**-1

    resIm, resRe = lmfitting.residual[1::2], lmfitting.residual[0::2]
    pp_errRE, pp_errIM = (lmfitting.best_fit.real - Z_KKv.real) ** 2, (
        lmfitting.best_fit.imag - Z_KKv.imag
    ) ** 2
    pp_errRE_mean, pp_errIM_mean = pp_errRE.mean(), pp_errIM.mean()

    #    MSE_Re,MSE_Im= (lmfitting.best_fit.real-Z_KKv.real)**2, (lmfitting.best_fit.imag-Z_KKv.imag)**2
    MSE = sum(pp_errRE) + sum(pp_errIM)
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


def get_best_pars_from_good_fits(
    EIS_fit_kwargs,
    Z_KKv,
    ang_KKv,
    params_model,
    standard_init_params,
    MTHDS,
    modname,
    model_set,
):
    br_eval_test = []
    good_pars_grp = EIS_fit_kwargs.get("good_fit_pars", {}).get(modname, [])
    _brutepars_lookup = EIS_fit_kwargs.get("random_params", {}).get(modname, [])
    _stdpar = standard_init_params
    brute_prefit_method = MTHDS[-5]
    #    params_model
    # random.sample((MTHDS[-5],MTHDS[0]),1)[0]
    #    _stdpar[el[0]].min,_stdpar[el[0]].max,
    if good_pars_grp:
        #         good_pars_mod = good_pars_grp.get(modname,[])
        #         good_pars_minmax_raw = {list(set(el[0]))[0] : {'min' : np.max([ np.min(el[1])*0.7]), 'max' : np.min([np.max(el[1])*1.25])}
        #         for el in[list(zip(*i))
        #         for i in  [map(operator.itemgetter(i), good_pars_grp)
        #         for i in range(0,len(good_pars_grp[0]))]]
        #         if np.mean(el[1]) > 0}
        good_pars_minmax_raw = {
            k: {
                "min": min(
                    [
                        dict(i).get(k)
                        for i in good_pars_grp
                        if set(dict(i).keys()) == set(params_model.valuesdict().keys())
                    ]
                ),
                "max": max(
                    [
                        dict(i).get(k)
                        for i in good_pars_grp
                        if set(dict(i).keys()) == set(params_model.valuesdict().keys())
                    ]
                ),
            }
            for k in params_model.valuesdict().keys()
            if not k == "Rs"
        }
        good_pars_minmax = {
            key: {
                "min": np.max([val["min"] * 0.1, _stdpar[key].min, 1e-12]),
                "max": np.min([val["max"] * 5, _stdpar[key].max, 1e12]),
            }
            for key, val in good_pars_minmax_raw.items()
        }
    else:
        good_pars_minmax = {
            k: {"min": params_model[k].min, "max": params_model[k].max}
            for k in params_model.valuesdict().keys()
            if not k == "Rs"
        }
    #            brutepar_lst_big = prepare_brute_pars(params_model,20E3)
    if _brutepars_lookup:
        good_pars_mod = good_pars_grp + _brutepars_lookup
    for bpar_eval in good_pars_mod:
        #        brute_prefit_method = random.sample((MTHDS[-5],MTHDS[0]),1)[0]
        _new_params = Parameters()
        bpar_fltr = dict((filter(lambda x: x[0] in params_model.keys(), bpar_eval)))
        for parnm in params_model.keys():
            if (
                params_model[parnm].vary
                and params_model[parnm].value > 0
                and parnm != "Rs"
                and parnm in bpar_fltr.keys()
            ):
                #                print(f'adapting {parnm,parv} {bpar_eval[0:2]}')
                parv = bpar_fltr[parnm]
                _min = good_pars_minmax.get(parnm, {"min": params_model[parnm].min})[
                    "min"
                ]
                _max = good_pars_minmax.get(parnm, {"max": params_model[parnm].max})[
                    "max"
                ]
                _new_params.add(
                    parnm, value=parv, vary=params_model[parnm].vary, min=_min, max=_max
                )
            elif parnm == "Rs":
                parv = (
                    _stdpar[parnm].value
                    if not "Rs" in bpar_fltr.keys()
                    else bpar_fltr[parnm]
                )
                _new_params.add(
                    parnm,
                    value=parv,
                    vary=[False, _stdpar[parnm].vary][1],
                    min=_stdpar[parnm].value * 0.00,
                    max=_stdpar[parnm].value * 1.5,
                )
            else:
                _new_params.add(
                    parnm,
                    value=params_model[parnm].value,
                    vary=params_model[parnm].vary,
                    min=_stdpar[parnm].min,
                    max=_stdpar[parnm].max,
                )
        #        _new_params.pretty_print()
        #        print(f'\n')
        #        _new_params['Rs'].set(value=Rs_guess,vary=True, min=0.8*Rs_guess, max = 1.2*Rs_guess)
        #                trial_prefit_brute = model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv,weights= statiscial_weights, method= brute_prefit_method)
        trial_prefit_brute_eval = model_set.eval(_new_params, ang=ang_KKv)
        MSE_re = (trial_prefit_brute_eval.real - Z_KKv.real) ** 2
        MSE_im = (trial_prefit_brute_eval.imag - Z_KKv.imag) ** 2
        MSE = sum(MSE_im + MSE_re)
        MSE_high = sum(MSE_im[-25:] + MSE_re[-25:])
        MSE_low = sum(MSE_im[:10] + MSE_re[:10])
        MSE_highlow = MSE_high + MSE_low
        br_eval_test.append((MSE, MSE_highlow, bpar_eval, _new_params))
    #    br_eval_test_all = sorted(br_eval_test, key=lambda MSE: MSE[0])
    if br_eval_test:
        br_eval_test_all = sorted(br_eval_test, key=lambda MSE: MSE[1])
        br_eval_test_all_MSE = sorted(br_eval_test, key=lambda MSE: MSE[0])[0:3]
        br_eval_len = np.min(
            [
                int(len(br_eval_test_all) * 0.3 + 1),
                10 if len(br_eval_test_all) > 10 else int(len(br_eval_test_all) / 2),
            ]
        )

        #        br_eval_test_all[0:br_eval_len]+br_eval_test_all_MSE
        br_eval_test_fit = []
        _Rs_test = {"Rs_vary": True, "Rs_fixed": False}
        for MSE, MSE_high, bpar_eval, _test_params in (
            br_eval_test_all[0:br_eval_len] + br_eval_test_all_MSE
        ):
            #            DE_kwgs = {'strategy' : 'rand1bin','popsize' : 20,'recombination' : 0.7, 'maxiter' : 100,'disp' : False}
            bpars = dict(bpar_eval)
            for _Rs_nm, _Rs_set in _Rs_test.items():
                if _Rs_set:
                    parv = (
                        _stdpar["Rs"].value if not "Rs" in bpars.keys() else bpars["Rs"]
                    )
                    _test_params["Rs"].set(vary=_Rs_set, value=parv)
                else:
                    _test_params["Rs"].set(vary=_Rs_set, value=_stdpar["Rs"].value)
                trial_prefit_brute_fit = model_set.fit(
                    Z_KKv, _test_params, ang=ang_KKv, method="least_squares"
                )
                #            make_prefit_frame(EIS_data_KKvalid, trial_prefit_brute_fit,plot = True, norm= wval) #TODO
                _aic, _redchi, _chi = (
                    trial_prefit_brute_fit.aic,
                    trial_prefit_brute_fit.redchi,
                    trial_prefit_brute_fit.chisqr,
                )
                br_eval_test_fit.append(
                    {
                        "_chi": _chi,
                        "_aic": _aic,
                        "_redchi": _redchi,
                        "lmfit": trial_prefit_brute_fit,
                        "lmfit_params": trial_prefit_brute_fit.params,
                        "_init_params": trial_prefit_brute_fit.init_params,
                        "delta_Rs": np.abs(
                            _stdpar["Rs"].value
                            - trial_prefit_brute_fit.params["Rs"].value
                        ),
                        "Rs_setting": _Rs_nm,
                    }
                )
        br_eval_fit_all = sorted(br_eval_test_fit, key=lambda MSE: MSE["_redchi"])
        br_eval_fit_DF = pd.DataFrame(br_eval_fit_all)
        #        br_eval_fit_DF.loc[br_eval_fit_DF[0] < br_eval_fit_DF[0].mean()].iloc[0]
        #        br_eval_fit_DF.loc[br_eval_fit_DF[0] < br_eval_fit_DF[0].mean()].iloc[0].to_numpy()
        #        for n,r in  br_eval_fit_DF.iterrows(): #TODO
        #            make_prefit_frame(EIS_data_KKvalid, r['lmfit'],plot = True) #TODO
        return br_eval_fit_DF
    else:
        return pd.DataFrame()


def old_prefit():
    #            else:
    #                params_model['nDL'].set(value= 0.5, vary=False)
    #                params_model['nAd'].set(value= 1, vary= False)
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
    wt_check = []
    for wname, wval in weight_opts.items():
        for methodnm in [MTHDS[0], MTHDS[-5]]:
            for pre_par_setnm, pre_par_set in par_options.items():
                pre_prefit = model_set.fit(
                    EIS_data_KKvalid["DATA_Z"].to_numpy(),
                    pre_par_set,
                    ang=ang_KKv,
                    weights=wval,
                    method=methodnm,
                )
                pp_good_low_high, best_test_msg, MSE = make_prefit_frame(
                    EIS_data_KKvalid, pre_prefit, prefix="pp", plot=0, norm=wval
                )
                wt_check.append(
                    [
                        wname,
                        MSE,
                        best_test_msg,
                        sum(MSE + best_test_msg),
                        methodnm,
                        pre_par_setnm,
                        pre_prefit,
                    ]
                )
    #                pre_prefit_pretrial = model_set.fit(EIS_data_KKvalid[pre_prefit_linkKK_data],pretrial_best,ang=ang_KKv, weights = wval, method= methodnm)
    #                pretrial_best
    #        _logger.info(f'Prefit test method: {MTHDS[0]} ({pre_prefit_leastsq.aic:.0f}), {MTHDS[-5]} ({pre_prefit_DE.aic:.0f})')
    #        pre_prefit_leastsq.redchi,pre_prefit_DE.redchi
    [
        make_prefit_frame(EIS_data_KKvalid, wt[-1], prefix="pp", plot=True)
        for wt in wt_check
    ]


#    #            wt_check.update({(modname,wname) : {'MSE' :MSE, 'best_test_msg' : best_test_msg,
#    #                                        'sum_test_MSE' : sum(MSE+best_test_msg),'lmfit_best_trial' : best_trial_weights}})
#                wt_check.append([wname,methodnm MSE,best_test_msg, sum(MSE+best_test_msg), best_trial_weights])
#            wt_check.append({modname : {'weights_name' : (wname,MSE,best_test_msg,sum(MSE+best_test_msg),best_trial_weights))
#        model_set.eval(pre_prefit.params,ang=np.logspace(-3,4,endpoint=True)*2.*pi, weights = statiscial_weights, method= pre_prefit.method)
#        PREFITS_weights = pd.DataFrame(data=wt_check.values(),index=wt_check.keys()).sort_values(by='MSE')
#%%
#        best_trial2 = model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv, weights = wval, method= 'leastsq')
#        run_brute_fit() #TODO
#        if run_brute_prefit == True and modname not in _bad_models:
#%%
#        prefit = model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv,weights= statiscial_weights, method='leastsq')
#            make_prefit_frame(EIS_data_KKvalid, prefit, prefix = 'pp',plot = 'Y')


def test_DE_emcee():

    DE_kwgs = {
        "strategy": "rand1bin",
        "popsize": 200,
        "recombination": 0.9,
        "maxiter": 1000,
        "disp": True,
    }
    trial_prefit_brute_fit = model_set.fit(
        Z_KKv,
        par_options.get("pretrial", par_options["standard"]),
        ang=ang_KKv,
        method="differential_evolution",
        fit_kws=DE_kwgs,
        weights=weight_opts["lmfit_weights_mod_Y"],
    )

    make_prefit_frame(EIS_data_KKvalid, trial_prefit_brute_fit, plot=True)  # TODO

    emcee_kws = dict(steps=5000, burn=300, thin=10, is_weighted=False, progress=False)
    emcee_params = result.params.copy()
    emcee_params.add("__lnsigma", value=np.log(0.1), min=np.log(0.001), max=np.log(2.0))
    result_emcee = model_set.fit(
        Z_KKv,
        emcee_params,
        ang=ang_KKv,
        method="emcee",
        nan_policy="omit",
        fit_kws=emcee_kws,
    )
    from pprint import pprint

    pprint(result_emcee.fit_report())
    make_prefit_frame(EIS_data_KKvalid, result_emcee, prefix="pp", plot=True)  # TODO
    best_trial = result_emcee
    emcee_corner = corner.corner(
        result_emcee.flatchain,
        labels=result_emcee.var_names,
        truths=list(result_emcee.params.valuesdict().values()),
    )


def test_emcee():
    mini = Minimizer(residual, _test_params, fcn_args=(Z_KKv, ang_KKv, model_set))
    trial_prefit_brute_fit = mini.minimize(
        method="differential-evolution",
        **{"strategy": "rand1bin", "popsize": 700, "recombination": 0.9},
    )
    #            method='leastsq')
    #            trial_prefit_brute_fit.plot()
    try:
        ci = conf_interval(mini, trial_prefit_brute_fit)
        res = minimize(
            residual,
            method="emcee",
            nan_policy="omit",
            burn=300,
            steps=1000,
            thin=50,
            params=mini.params,
            args=(Z_KKv, ang_KKv, model_set),
            is_weighted=False,
            progress=False,
        )
    except:
        ci = None

    result = model_set.fit(Z_KKv, _test_params, ang=ang_KKv, method="Nelder")
    make_prefit_frame(EIS_data_KKvalid, result, plot=True)  # TODO
    emcee_kws = dict(steps=1000, burn=300, thin=20, is_weighted=False, progress=False)
    emcee_params = result.params.copy()
    #    emcee_params.add('__lnsigma', value=np.log(0.1), min=np.log(0.001), max=np.log(2.0))
    result_emcee = model_set.fit(
        Z_KKv,
        emcee_params,
        ang=ang_KKv,
        method="emcee",
        nan_policy="omit",
        fit_kws=emcee_kws,
    )
    pprint(result_emcee.fit_report())

    make_prefit_frame(EIS_data_KKvalid, result_emcee, plot=True)  # TODO
    plt.plot(result_emcee.acceptance_fraction)
    plt.xlabel("walker")
    plt.ylabel("acceptance fraction")
    plt.show()

    if hasattr(result_emcee, "acor"):
        print("Autocorrelation time for the parameters:")
        print("----------------------------------------")
        for i, p in enumerate(result.params):
            print(p, result.acor[i])
    emcee_corner = corner.corner(
        result_emcee.flatchain,
        labels=result_emcee.var_names,
        truths=list(result_emcee.params.valuesdict().values()),
    )


# params = _test_params
def residual(params, Z_KKv, ang, model, weights=None):
    model = model_set.eval(params, ang=ang)
    MSE_re = (model.real - Z_KKv.real) ** 2
    MSE_im = (model.imag - Z_KKv.imag) ** 2
    MSE = MSE_re + MSE_im
    resid = model - Z_KKv
    return resid.view(np.float)


def prefit_test(
    fit_run_arg,
    EIS_fit_kwargs,
    weight_opts,
    EIS_data_KKvalid,
    Z_KKv,
    ang_KKv,
    standard_init_params,
    MTHDS,
    prefit_methods=["leastsq"],
    prefit_models=EEC_models_index(),
    pickle_spectra=False,
):
    pars_lst, fit_data_valid_lst, outP = [], [], {}
    unit_weight = weight_opts["lmfit_weights_unit"]
    PREFITS_spectra = pd.DataFrame()

    wt_check, wt_spectra_lst = {}, []
    for mod_indx, model_set, rand_n in prefit_models:  #
        modname = model_set.name
        params_model = model_set.make_params()
        lstsq_method = MTHDS[0]
        par_options = {}
        for pn in params_model:
            params_model[pn].set(
                value=standard_init_params[pn].value,
                min=standard_init_params[pn].min,
                max=standard_init_params[pn].max,
            )
        try:
            _get_params = pd.DataFrame()
            _get_params, recheck_msg = fitting_recheck_params(
                fit_run_arg, modname, params_model, **EIS_fit_kwargs
            )
            #            _logger.warning(recheck_msg)
            if not _get_params.empty:
                for pn in params_model:
                    if (
                        params_model[pn].vary == True
                        and pd.isna(_get_params[pn].iloc[0]) == False
                        and pn != "Rs"
                    ):
                        params_model[pn].set(value=_get_params[pn].iloc[0])
            #                _logger.info(f'Prefit used recheked params for {modname}')

            else:
                _logger.info(
                    f'Prefit _get_params empty {(Path(fit_run_arg.PAR_file).name, fit_run_arg.E_dc_RHE, fit_run_arg.RPM_DAC)}\n message:"{recheck_msg}'
                )
        except Exception as e:
            _logger.error(f"fitting_recheck_params error: {e}")
        par_options.update({"standard": params_model})
        #        std_prefit_DE = model_set.fit(Z_KKv,params_model, ang=ang_KKv, weights= wval, method= MTHDS[5])
        _good_fits_par_fit = get_best_pars_from_good_fits(
            EIS_fit_kwargs,
            Z_KKv,
            ang_KKv,
            params_model,
            standard_init_params,
            MTHDS,
            modname,
            model_set,
        )
        #         make_prefit_frame(EIS_data_KKvalid, _good_fits_par_fit[3],plot = True)
        #        br_eval_fit_DF.loc[br_eval_fit_DF[0] < br_eval_fit_DF[0].mean()].iloc[0].to_numpy()
        if not _good_fits_par_fit.empty:
            #            _good_params = _good_fits_par_fit[-2]
            for _Rsnm, _Rsgrp in _good_fits_par_fit.groupby("Rs_setting"):
                _Rsnm, _Rsgrp
                _Rsgrp = _Rsgrp.sort_values(by=["delta_Rs", "_chi"])
                par_options.update({_Rsnm: _Rsgrp.iloc[0]["lmfit_params"]})

        #            if 'Rs' in _good_params.keys():
        #                _good_params['Rs'].set(value = standard_init_params['Rs'].value,
        #                                         min=standard_init_params['Rs'].min,
        #                                         max=standard_init_params['Rs'].max)
        #                par_options.update({'good_params' : _good_params})
        #        if not PRETRIALS_weights.empty:
        #            pretrial_best = PRETRIALS_weights.loc[(mo  params_model = model_set.make_params()dname, PRETRIALS_weights.loc[modname]['MSE'].idxmin()),'lmfit_best_trial']
        #            par_options.update({'pretrial' : pretrial_best.params})
        ##        for pre_par_setnm, pre_par_set in  par_options.items():
        for parnm, par_obj in par_options.items():
            prefit_params = model_set.make_params()
            for pn in prefit_params:
                prefit_params.add(
                    pn,
                    value=par_obj[pn].value,
                    min=standard_init_params[pn].min,
                    max=standard_init_params[pn].max,
                )
            for wname, wval in weight_opts.items():
                for method in prefit_methods:
                    #                    pre_prefit_leastsq = model_set.fit(Z_KKv,par_obj, ang=ang_KKv, weights= wval, method= MTHDS[0])
                    #                pre_prefit_DE = model_set.fit(Z_KKv,par_obj,ang=ang_KKv, weights = wval, method= MTHDS[1])
                    #        _logger.info(f'Prefit test method: {MTHDS[0]} ({pre_prefit_leastsq.aic:.0f}), {MTHDS[-5]} ({pre_prefit_DE.aic:.0f})')
                    #        pre_prefit_leastsq.redchi,pre_prefit_DE.redchi
                    #                if pre_prefit_leastsq.aic < pre_prefit_DE.aic:
                    #                    pre_prefit = pre_prefit_leastsq
                    #    #                make_prefit_frame(EIS_data_KKvalid, pre_prefit_leastsq,plot = True, norm= statiscial_weights)
                    #    #                succes_msg = 'Fit succeeded'
                    #        #            _logger.info(f'Best prefit method: {MTHDS[0]}')
                    #                else:
                    #        #            _logger.info(f'Best prefit method: {MTHDS[-5]}')
                    #    #                make_prefit_frame(EIS_data_KKvalid, pre_prefit_DE,plot = True, norm= wval)
                    #                    pre_prefit = pre_prefit_DE
                    #    #                succes_msg = 'successfully'
                    #            make_prefit_frame(EIS_data_KKvalid, pre_prefit,plot = True, norm= wval) #TODO
                    best_trial_weights = model_set.fit(
                        Z_KKv, prefit_params, ang=ang_KKv, weights=wval, method=method
                    )
                    #                    best_trial_weights = model_set.fit(Z_KKv,pre_prefit.params,ang=ang_KKv, weights = unit_weight , method= pre_prefit.method)
                    pp_good_low_high, best_test_msg, MSE = make_prefit_frame(
                        EIS_data_KKvalid,
                        best_trial_weights,
                        prefix="pp",
                        plot=0,
                        norm=unit_weight,
                    )

                    #                     = make_prefit_frame(EIS_data_KKvalid,best_trial_weights, prefix = 'pp',plot = 0, norm= unit_weight)
                    if pickle_spectra:
                        post_br_frame = make_prefit_frame(
                            EIS_data_KKvalid, best_trial_weights, get_frame=1
                        )
                        post_br_frame = post_br_frame.assign(
                            **{
                                "Model_EEC": modname,
                                "lmfit_Weights": wname,
                                "lmfit_Rs_setting": parnm,
                                "lmfit_method": method,
                                "mod_index": ", ".join((modname, wname, parnm, method)),
                            }
                        )
                        wt_spectra_lst.append(post_br_frame)
                    #                    merge_cols = list(post_br_frame.columns.intersection(EIS_data_KKvalid.columns))
                    #                    EIS_data_KKvalid_BR_specs = pd.merge(EIS_data_KKvalid_BR_specs,post_br_frame,on =merge_cols,how='inner',left_index=True,right_index=True)

                    wt_check.update(
                        {
                            (modname, wname, parnm, method): {
                                "Model_EEC": modname,
                                "lmfit_Weights": wname,
                                "lmfit_Rs_setting": parnm,
                                "lmfit_method": method,
                                "MSE": MSE,
                                "best_test_msg": best_test_msg,
                                "sum_test_MSE": sum(MSE + best_test_msg),
                                "lmfit_method": best_trial_weights.method,
                                "lmfit_aic": best_trial_weights.aic,
                                "lmfit_succes": best_trial_weights.success,
                                "lmfit_best_trial": best_trial_weights,
                                "BRUTE_FIT": 1,
                                "FINAL_FIT": 0,
                                "PREFIT": 0,
                                **best_trial_weights.params.valuesdict(),
                            }
                        }
                    )
    #            wt_check.append({modname : {'weights_name' : (wname,MSE,best_test_msg,sum(MSE+best_test_msg),best_trial_weights))
    #        model_set.eval(pre_prefit.params,ang=np.logspace(-3,4,endpoint=True)*2.*pi, weights = statiscial_weights, method= pre_prefit.method)
    PREFITS_weights = pd.DataFrame(
        data=wt_check.values(), index=wt_check.keys()
    ).sort_values(by="MSE")
    PRE_min = PREFITS_weights.iloc[
        PREFITS_weights.reset_index()
        .groupby("level_0")["MSE"]
        .transform("idxmin")
        .unique()
    ]
    if wt_spectra_lst:
        PREFITS_spectra = pd.concat(wt_spectra_lst, sort=False, ignore_index=True)
    #    [make_prefit_frame(EIS_data_KKvalid, r['lmfit_best_trial'],plot = True, norm= wval) for n,r in PREFITS_weights.iterrows()] #TODO
    min_idxs = (
        PREFITS_weights.groupby("Model_EEC")["sum_test_MSE"]
        .transform("idxmin")
        .unique()
    )
    post_check = {}
    PREFITS_weights.loc[min_idxs, "PREFIT"] = 1
    PRETRIALS_weights = PREFITS_weights.loc[min_idxs]
    #    [make_prefit_frame(EIS_data_KKvalid, r['lmfit_best_trial'],plot = True, norm= wval) for n,r in PRETRIALS_weights.iterrows()] #TODO
    #    PP_out = pd.DataFrame(data=post_check.values(),index=post_check.keys()).sort_values(by='MSE')
    #    PRETRIALS_weights = PP_out.iloc[PP_out.reset_index().groupby('level_0')['MSE'].transform('idxmin').unique()]
    _bad_models = PRETRIALS_weights.loc[
        PRETRIALS_weights.MSE
        > PRETRIALS_weights.MSE.mean() + PRETRIALS_weights.MSE.std()
    ]
    _bad_models_out = [i[0] for i in _bad_models.index]
    #    [gr['MSE'].idxmin() for n,gr in PREFITS_weights.reset_index().groupby('level_0')]
    #    PREFITS_weights['MSE'].idxmin()
    return PRETRIALS_weights, _bad_models_out, PREFITS_weights, PREFITS_spectra


def prefit_loop_over_models():
    for idx in min_idxs[::]:
        prefit_row = PREFITS_weights.iloc[idx]
        prefit_params = prefit_row.lmfit_best_trial.params
        prefit_R_pars = [
            (k, val)
            for k, val in prefit_params.valuesdict().items()
            if k in ["Rct", "Rorr", "R3"]
        ]
        for mod_indx, model_set, rand_n in EEC_models_index():  #
            modname = model_set.name
            params_model_test = model_set.make_params()
            lstsq_method = MTHDS[0]

            test_R_pars = [
                (k, val)
                for k, val in params_model_test.valuesdict().items()
                if k in ["Rct", "Rorr", "R3"]
            ]
            test_R_pars = sorted(test_R_pars, key=lambda x: x[1], reverse=True)
            for pn in params_model_test:
                if not np.isnan(prefit_row[pn]) and pn != "Rs":
                    params_model_test.add(prefit_params[pn])
                else:
                    params_model_test.add(standard_init_params[pn])
            #            if len(prefit_R_pars) < len(test_R_pars):
            #                for R in list(filter(lambda x: x in [i[0] for i in prefit_R_pars], ['Rct','Rorr','R3'])):
            #                    params_model_test[R].value = test_R_pars.pop()[1]
            #                filter(['Rct','Rorr','R3'],prefit_R_pars):
            for pn in params_model_test:
                params_model_test[pn].set(
                    value=random.gauss(
                        params_model_test[pn].value, 0.07 * params_model_test[pn].value
                    )
                )

            post_pre_lmfit = model_set.fit(
                Z_KKv,
                params_model_test,
                ang=ang_KKv,
                weights=unit_weight,
                method="leastsq",
            )
            pp_good_low_high, best_test_msg, MSE = make_prefit_frame(
                EIS_data_KKvalid, post_pre_lmfit, prefix="pp", plot=0, norm=unit_weight
            )
            post_check.update(
                {
                    (modname, idx): {
                        "MSE": MSE,
                        "best_test_msg": best_test_msg,
                        "lmfit_aic": post_pre_lmfit.aic,
                        "sum_test_MSE": sum(MSE + best_test_msg),
                        "lmfit_best_trial": post_pre_lmfit,
                        **post_pre_lmfit.params.valuesdict(),
                    }
                }
            )
    PP_out = pd.DataFrame(
        data=post_check.values(), index=post_check.keys()
    ).sort_values(by="MSE")
    PRETRIALS_weights = PP_out.iloc[
        PP_out.reset_index().groupby("level_0")["MSE"].transform("idxmin").unique()
    ]


def fit_mean_PAR_file(fit_run_args):

    PF_grp_fit_run_args = itertools.groupby(fit_run_args, lambda x: x.PAR_file)

    PF, PF_run_args_gen = next(PF_grp_fit_run_args)


def _get_mean_EIS_from_args(PF, PF_run_args_gen):
    global EvRHE

    _new_PF_mean = PF.with_name(PF.stem + "_Zmean" + PF.suffix)

    _prep_mean_args = list(PF_run_args_gen)
    _PF_mean_ovv = _prep_mean_args[0].ovv
    #    _PF_mean_ovv['Measured_OCP'] =  [i[0] for i in _PF_mean_ovv['_act0_Measured Open Circuit'].str.split()]
    #    _PF_mean_ovv['PAR_file'] = _new_PF_mean
    _PF_mean_ovv = _PF_mean_ovv.assign(**{"PAR_file": _new_PF_mean})
    _PF_data = pd.concat(i.data for i in _prep_mean_args)
    _numcols = [
        i for i in _PF_data.columns if _PF_data[i].dtype != "O" and not "Frequency" in i
    ]

    _PF_data_mean = _PF_data.groupby("Frequency(Hz)")[_numcols].mean().reset_index()
    #    _PF_data_mean.plot(x='Z Real',y='Z Imag', kind='scatter') # test plot

    _PF_data_mean = _PF_data_mean.assign(**_PF_mean_ovv.iloc[0].to_dict())
    _join_cols = list(_PF_data_mean.columns.intersection(_PF_mean_ovv.columns))
    #    _PF_data_mean = _PF_data_mean.join(_PF_mean_ovv.set_index('PAR_file'),on=_join_cols,how='left')
    #    _merge = pd.concat([_PF_data_mean, _PF_mean_ovv],axis=1)
    #                       .set_index('PAR_file'),on=_join_cols,how='outer')
    _PF_data_mean[["Segment #", EvRHE, "RPM_DAC"]]
    _PF_data_mean_grp = _PF_data_mean.groupby(
        ["PAR_file", "Segment #", EvRHE, "RPM_DAC"]
    )
    fit_run_arg_mean = [
        Meta(
            Path(PF),
            int(Seg),
            np.round(float(E_V), 3),
            np.round(float(E_V) * 1e3, 3),
            int(RPM_DAC),
            gr,
            _PF_mean_ovv,
        )
        for (PF, Seg, E_V, RPM_DAC), gr in _PF_data_mean_grp
    ][0]
    #    _PF_mean_ovv[0]._fields
    try:
        a = fit_EEC(fit_run_arg_mean, **EIS_fit_kwargs)
        EIS_fit_data, pars_models = a
    except Exception as e:
        print(e)


#    fit_run_arg=fit_run_arg_mean
class EIS_loop:
    def __init__(self, fit_run_args):
        self.fit_run_args = fit_run_args


#    def fit_mean_PAR_file(self):
#        for PF in itertools.groupby(self.fit_run_args, lambda x: x.PAR_file):


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
        _join_cols = list(_PF_data_mean.columns.intersection(_PF_mean_ovv.columns))
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

    def __init__(self, _EIS_spectra_pf):
        assert type(_EIS_spectra_pf).__name__ == "EIS_spectra_collection"
        #        print(isinstance(_EIS_spectrum_arg,EIS_Spectrum))
        self._spectra = _EIS_spectra_pf
        self.fit_mean = Fit_Spectrum(
            self._spectra.mean_spectrum, linKK_trimming_factor=4, res_max=0.2
        )

    #        self.fit_mean.linKK(linKK_trimming_factor = 4, res_max = 0.2)
    #        self.fit_mean()

    def fit_mean(self):
        pass


#        self._spectra.EIS_kwargs.update({'linKK_trimming_factor' : 3})
#        self.fit_mean = Fit_Spectrum(self._EIS_spectra_pf.mean_spectrum)

#    def linKK_mean(self):
#        self.fit_mean.linKK(linKK_trimming_factor = 4, res_max = 0.2)


class Fit_Spectrum:
    """This class will take the fit_run_arg and
    run the steps for fitting the EIS spectrum"""

    global EvRHE

    def __init__(self, _EIS_spectrum_arg, **kwargs):
        assert type(_EIS_spectrum_arg).__name__ == "EIS_Spectrum"
        #        print(isinstance(_EIS_spectrum_arg,EIS_Spectrum))
        self._spectrum = _EIS_spectrum_arg
        self._spectrum.EIS_kwargs.update(**kwargs)

        self.prepare_meta_data()
        self.linKK(**kwargs)
        self._add_weights_col()
        self.Rs_guess()
        self.Warburg_type_fitting()

    #        self.fit_run_arg = fit_run_arg
    def __getattr__(self, attr):
        return getattr(self._spectrum, attr)

    def __repr__(self):
        return f"Fit_Spectrum: {self._spectrum.__repr__()}"

    def prepare_meta_data(self):
        #        self.dataself.ovv
        PF_EIS_dest_dir = self.data.EIS_dest_dir.unique()[0]
        #    it_run_arg.ovv.Dest_dir.iloc[0]
        # FIXME preparing meta dicts
        EISgr_data_meta = (
            self.data[[i for i in self.data.columns if self.data[i].nunique() == 1]]
            .iloc[0]
            .to_dict()
        )
        gr_EIS_ovv_meta = self.ovv.iloc[0].to_dict()
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
                str(self.PAR_file.stem)
                + f"_{self.E_dc_RHE_mV:.0f}mV_{self.RPM_DAC:.0f}rpm_{self.Segment}"
            )
        )
        self.EIS_outPath = EIS_outPath
        #    EIS_outPath_OLD = EIS_kwargs.get('EIS_dest_dir').joinpath(Path(str(Path(EISgr_data_meta['PAR_file']).stem)+'_{:.0f}mV'.format(fit_run_arg.E_dc_RHE_mV)))
        EIS_outPath.parent.mkdir(parents=True, exist_ok=True)

        _msg = f"{str(self.PAR_file.name):s} seg({self.Segment}) at {self.E_dc_RHE:G} ({self.RPM_DAC} rpm)"
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
            EIS_data_KKvalid, linkKK_invalid_prefit, linKK_pars = get_KKvalid(
                self.EIS_data_freqlim, **self._spectrum.EIS_kwargs
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
        #        self.EIS_data = self.EIS_data.loc[~self.EIS_data['Frequency(Hz)'].isin(EIS_data_KKvalid['Frequency(Hz)']),'Valid'] = False
        #        linKK_pars.update({'linKK_invalid_freqs' : ', '.join(EIS_data_KKvalid.loc[~EIS_data_KKvalid['Frequency(Hz)'].isin(EIS_data_KKvalid['Frequency(Hz)']),'Frequency(Hz)'].astype(str).to_list())})
        try:
            plot_linKK(
                EIS_data_KKvalid,
                self.EIS_data_freqlim,
                save_target=EIS_outpath_linKK_target,
                plot_prefit=True,
                linkKK_invalid_prefit=linkKK_invalid_prefit,
                **{**linKK_pars, **self.spectrum_meta_info},
            )
        except Exception as e:
            _logger.warning(
                f"EIS in plotting linKK, {str(self.PAR_file):s} at {self.E_dc_RHE:G} ({self.RPM_DAC} rpm).\n {e}"
            )

        if self._spectrum.EIS_kwargs.get("export_raw_testing", False):
            testing_path = Path.cwd().joinpath(
                "testing_data", self.EIS_outPath.name + "_spectrum_raw"
            )
            self.EIS_data_KKvalid.to_excel(testing_path.with_suffix(".xlsx"))

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

    def Rs_guess(self, Rs_limits=(2, 200)):
        Rs_guess_data = self.EIS_data_KKvalid.loc[
            np.abs(self.EIS_data_KKvalid.DATA_Zmod).idxmin(), "DATA_Zre"
        ]
        Rs_guess_linKK = self.EIS_data_KKvalid.linKK_Z.to_numpy().real.min()
        Rs_guess = np.mean(
            [
                i
                for i in [Rs_guess_data, Rs_guess_linKK]
                if i > Rs_limits[0] and i < Rs_limits[1]
            ]
        )
        if not Rs_guess:
            Rs_guess_data = Rs_guess
        self.Rs_guess = Rs_guess

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
                _logger.error(f"EIS fit_EEC error in GP_DRT_fit, {self._errmsg}.\n {e}")

    def Warburg_type_fitting(self, _lin_window_size=11):
        _lin = {}
        lin_slopes = [(0.25, "lightgreen"), (0.5, "grey"), (1, "orange")]
        #    zip([0.25, 0.5, 1], ['lightgreen','grey','orange'])
        #        for yax in ['Zre','-Zim']:
        #        _lin.update({yax : linregress(spec.query('ang_Warburg > 0.3').ang_Warburg,spec.query('ang_Warburg > 0.3')['DATA_'+yax]),
        #                     yax+'_low' : linregress(spec.query('ang_Warburg < 0.045').ang_Warburg,spec.query('ang_Warburg < 0.045')['DATA_'+yax])})
        #        spec['W_lin_'+yax] = _lin[yax].slope * spec.ang_Warburg + _lin[yax].intercept
        #        spec['W_lin_'+yax+'_low'] = _lin[yax+'_low'].slope * spec.ang_Warburg + _lin[yax+'_low'].intercept
        _lin_WB_freq_loq = self.EIS_data_KKvalid.query("Angular < 30")
        _lin_WB_angW_high = self.EIS_data_KKvalid.query("Ang_Warburg > 0.3")
        _lin_WB_angW_low = self.EIS_data_KKvalid.query("Ang_Warburg < 0.045")

        for yax in ["Zre", "-Zim"]:
            _yax_cols = {}
            _DATA_Z_col = f"DATA_{yax}"
            _lin.update(
                {
                    f"{yax}_angW_lin_high": linregress(
                        _lin_WB_angW_high.Ang_Warburg, _lin_WB_angW_high[_DATA_Z_col]
                    ),
                    f"{yax}_angW_lin_low": linregress(
                        _lin_WB_angW_low.Ang_Warburg, _lin_WB_angW_low[_DATA_Z_col]
                    ),
                }
            )
            if yax == "Zre":
                _lin.update(
                    {
                        f"{yax}_lin_slopes": linregress(
                            _lin_WB_freq_loq["DATA_Zre"], _lin_WB_freq_loq["DATA_-Zim"]
                        )
                    }
                )
                _yax_cols.update(
                    {
                        f"{yax}_lin_slopes": _lin[f"{yax}_lin_slopes"].slope
                        * self.EIS_data_KKvalid[_DATA_Z_col]
                        + _lin[f"{yax}_lin_slopes"].intercept
                    }
                )

            for _WB_slice_ax in [
                f"{yax}{i}" for i in ["_angW_lin_high", "_angW_lin_low"]
            ]:
                _yax_cols.update(
                    {
                        _WB_slice_ax: _lin[_WB_slice_ax].slope
                        * self.EIS_data_KKvalid.Ang_Warburg
                        + _lin[_WB_slice_ax].intercept
                    }
                )
            self.EIS_data_KKvalid = self.EIS_data_KKvalid.assign(**_yax_cols)

            #                self.EIS_data_KKvalid['W_lin_'+yax] = _lin[yax].slope * self.EIS_data_KKvalid.ang_Warburg + _lin[yax].intercept
            #            self.EIS_data_KKvalid['W_lin_'+yax+'_low'] = _lin[yax+'_low'].slope * self.EIS_data_KKvalid.ang_Warburg + _lin[yax+'_low'].intercept
            for _slope, _ in lin_slopes:
                perr_set = 1000
                #            for _win_size in [7,10,15,25]:
                #            for win in spec.rolling(_lin_window_size):
                for i in range((len(self.EIS_data_KKvalid) - _lin_window_size)):
                    win = self.EIS_data_KKvalid.iloc[i : i + _lin_window_size]

                    popt, pcov = curve_fit(
                        func_lin(_slope), win.DATA_Zre, win["DATA_" + yax]
                    )
                    perr = np.sqrt(np.diag(pcov))
                    #                print(win.index,popt,pcov,perr)
                    if perr < perr_set:
                        perr_set = perr
                        best = (_slope, win, popt, perr)
                #            popt, pcov = curve_fit(func_lin(_slope), spec.query('Angular > 30').DATA_Zre, spec.query('Angular > 30')['DATA_'+yax])
                self.EIS_data_KKvalid[f"Z_lin_WB_a{_slope}"] = func_lin(best[0])(
                    self.EIS_data_KKvalid.DATA_Zre, best[2][0]
                )
                _lin.update(
                    {
                        _slope: {
                            "popt": best[2][0],
                            "win_size": len(best[1]),
                            "perr": best[-1],
                        }
                    }
                )

        plot_lin_Warburg(
            self.EIS_data_KKvalid,
            _lin,
            lin_slopes,
            dest_path=self._path_and_subdir(self.EIS_outPath, "lin_Warburg"),
        )
        self.lin_fitting_Warburg = _lin
        # ====!!!==== #


def fit_EEC(fit_run_arg, **EIS_fit_kwargs):
    #%%
    EISgr_data_EV = fit_run_arg.data.sort_values("Frequency(Hz)", ascending=True)
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
            _logger.error(
                "EIS missing keys in EISgr_data_meta and ovv_meta dicts, {:s} at {:G}".format(
                    fit_run_arg.PAR_file.stem, fit_run_arg.E_dc_RHE
                )
            )
    else:
        _logger.error(
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
        _logger.error(
            "EIS outP missing columns in EISgr_data_meta, {:s} at {:G}  mV".format(
                fit_run_arg.PAR_file.stem, fit_run_arg.E_dc_RHE_mV
            )
        )

    EIS_outPath = PF_EIS_dest_dir.joinpath(
        Path(
            str(fit_run_arg.PAR_file.stem)
            + f"_{fit_run_arg.E_dc_RHE_mV:.0f}mV_{fit_run_arg.RPM_DAC:.0f}rpm_{fit_run_arg.Segment}"
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

    DataWeights_modulus_Z = Zre ** 2 + Zim ** 2
    DataWeights_modulus_Y = (Zre ** 2 + Zim ** 2) ** -1
    #    DataWeights_modulus_Y = (Yre**2+Yim**2)
    DataWeights_unitary = Zre / Zre

    EIS_data_raw = pd.DataFrame(
        {
            "Frequency(Hz)": freq,
            "Angular": ang,
            "Ang_Warburg": 1 / (np.sqrt(ang)),
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
            "lmfit_weights_mod_Z": DataWeights_modulus_Z,
            "lmfit_weights_mod_Y": DataWeights_modulus_Y,
            "lmfit_weights_unit": DataWeights_unitary,
            "PAR_file": fit_run_arg.PAR_file,
            "Segment": fit_run_arg.Segment,
        },
        index=EISgr_data_EV.index,
    )

    # ====Important Frequency limit slice on EIS data!!!==== #
    FreqLim = EIS_fit_kwargs.get("FreqLim", 30e3)
    EIS_data_raw.loc[
        (EIS_data_raw["Frequency(Hz)"] >= FreqLim) & (EIS_data_raw["DATA_-Zim"] < -5),
        "Valid",
    ] = False  # TODO
    EIS_data_freqlim = EIS_data_raw.query("Valid == True")
    # ====!!!==== #
    """ include here KK-test!! """
    #%%
    try:
        EIS_data_KKvalid, linkKK_invalid_prefit, linKK_pars = get_KKvalid(
            EIS_data_freqlim, restype_set="Z", res_scale=1, **EIS_fit_kwargs
        )
    except Exception as e:
        _logger.error(
            f"EIS fit_EEC error in validation, {str(fit_run_arg.PAR_file):s} at {fit_run_arg.E_dc_RHE:G} ({fit_run_arg.RPM_DAC} rpm).\n {e}"
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
    if EIS_fit_kwargs.get("export_raw_testing", False):
        testing_path = Path.cwd().joinpath(
            "testing_data", EIS_outPath.name + "_spectrum_raw"
        )
        EIS_data_KKvalid.to_excel(testing_path.with_suffix(".xlsx"))
    #
    if EIS_fit_kwargs.get("DP_DRT_fit", False):
        DP_DRT_destdir = EIS_outPath.parent.joinpath("DP_DRT")
        DP_DRT_destdir.joinpath("DP_DRT").mkdir(parents=True, exist_ok=True)
        DP_DRT_out = DP_DRT_analysis(EIS_data_KKvalid.DATA_Z.to_numpy())
    #        EIS_outPath.name+'_DP_DRT_'

    if EIS_fit_kwargs.get("GP_DRT_fit", False):
        try:
            GP_DRT_destpath = EIS_outPath.parent.joinpath("GP_DRT").joinpath(
                EIS_outPath.name
            )
            GP_DRT_destpath.parent.mkdir(parents=True, exist_ok=True)
            DRT_data_KKvalid = EIS_data_raw.sort_values(
                by="Frequency(Hz)", ascending=True
            )
            #            prep_GP_DRT_raw_data(EIS_data_raw)
            _Z_exp, _Z_star, res_fit_params, res_fit_arrays = run_GP_DRT_fit(
                DRT_data_KKvalid, **{"GP_DRT_savefig_path": GP_DRT_destpath}
            )
            EISgr_meta_add.update(
                {"GP_DRT_fit": GP_DRT_destpath, "GP_DRT_run_success": True}
            )
        except Exception as e:
            _logger.error(
                f"EIS fit_EEC error in GP_DRT_fit, {str(fit_run_arg.PAR_file):s} at {fit_run_arg.E_dc_RHE:G} ({fit_run_arg.RPM_DAC} rpm).\n {e}"
            )

    #    M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(EIS_data_freqlim['Frequency(Hz)'].values,EIS_data_freqlim['DATA_Z'].values,type_res=restype_set)
    #%%
    # if True:
    #    EIS_data_KKvalid = EIS_data.query('Valid == True')
    Z_KKv, ang_KKv = (
        EIS_data_KKvalid.DATA_Z.to_numpy(),
        EIS_data_KKvalid.Angular.to_numpy(),
    )
    freq_KKv = EIS_data_KKvalid["Frequency(Hz)"].to_numpy()
    ### Warburg check
    #    lin_Warburg= stats.linregress(np.log10(freq_KKv),np.log10(np.abs(Z_KKv)))
    #    lin_Warburg= stats.linregress(Z_KKv.real[0:20],-1*Z_KKv.imag[0:20])
    ##                                  np.angle(EIS_data_KKvalid.linKK_Z.values,deg=True))
    #    plt.plot(Z_KKv.real,-1*Z_KKv.imag)
    #    plt.plot(Z_KKv.real,lin_Warburg.slope*Z_KKv.real+lin_Warburg.intercept)
    #    plt.plot(Z_KKv.real,1*Z_KKv.real+lin_Warburg.intercept)
    ### Adding Zangle lin weights
    lin_angle = stats.linregress(
        np.log(EIS_data_KKvalid["Frequency(Hz)"].to_numpy()),
        np.angle(EIS_data_KKvalid.linKK_Z.values, deg=True),
    )
    lin_subt_Zangle = np.angle(EIS_data_KKvalid.linKK_Z.values, deg=True) - (
        lin_angle.slope * np.log(freq_KKv) + lin_angle.intercept
    )
    lmfit_weights_linZang = (lin_subt_Zangle) ** 2 / np.mean(
        (lin_subt_Zangle) ** 2
    ) + EIS_data_KKvalid.lmfit_weights_unit
    EIS_data_KKvalid = EIS_data_KKvalid.assign(
        **{"lmfit_weights_linZangle": lmfit_weights_linZang}
    )
    #    plt.plot(np.log10(ang_KKv),np.log10(np.abs(Z_KKv)))
    #    plt.plot(np.log(freq_KKv),np.angle(EIS_data_KKvalid.linKK_Z.values,deg=True))
    #    plt.plot(np.log(freq_KKv),(lin_angle.slope*np.log(freq_KKv)+ lin_angle.intercept))
    #    plt.plot(np.log(freq_KKv),lmfit_weights_linZang+EIS_data_KKvalid.lmfit_weights_unit)
    weight_opts = {
        key: EIS_data_KKvalid[key].values
        for key in EIS_data_KKvalid.columns
        if key.startswith("lmfit_weights_")
    }
    #    del weight_opts['lmfit_weights_mod_Z']
    #    statiscial_weights = EIS_data_KKvalid.DATA_weightsmod_Z.values
    #    unitary_weights = EIS_data_KKvalid.DATA_weightsmod_Z.values / EIS_data_KKvalid.DATA_weightsmod_Z.values
    #    pre_prefit_linkKK_data = 'DATA_Z'
    #    if 'linKK_Z' in EIS_data_KKvalid.columns:
    #        if np.sum(EIS_data_KKvalid.linKK_resRe**2 + EIS_data_KKvalid.linKK_resIm**2) < 5E-3:
    #            pre_prefit_linkKK_data = 'linKK_Z'

    # create a set of Parameters
    Rs_guess_data = EIS_data_KKvalid.loc[
        np.abs(EIS_data_KKvalid.DATA_Zmod).idxmin(), "DATA_Zre"
    ]
    Rs_guess_linKK = EIS_data_KKvalid.linKK_Z.to_numpy().real.min()
    Rs_guess = np.mean([Rs_guess_data, Rs_guess_linKK])
    if Rs_guess > 5 and Rs_guess < 200:
        pass
    else:
        Rs_guess = EIS_data_KKvalid.linKK_Z.to_numpy().real.min()

    standard_init_params = Parameters()
    lower_bound, max_ub = 1e-8, 1e10
    max_C_ub = 10
    n_lb = 0.0
    standard_init_params.add(
        "Rs",
        value=Rs_guess,
        min=0 * Rs_guess,
        max=1.5 * Rs_guess,
        brute_step=1,
        vary=True,
    )
    standard_init_params.add(
        "Ls", value=2e-05, min=lower_bound, max=1, brute_step=1e-03
    )
    standard_init_params.add(
        "Cdlp", value=1e-05, min=lower_bound, max=max_C_ub, brute_step=0.5e-06
    )
    standard_init_params.add("nDL", value=0.50, min=n_lb, max=1, brute_step=0.05)
    standard_init_params.add(
        "Rct", value=25, min=lower_bound, max=max_ub, vary=True, brute_step=10
    )
    #    params.add('Aw', value= 3000, min=5.0, max=1E6)
    #    params.add('Rad', value= 5E5, min=0.01, max=1E7)
    standard_init_params.add(
        "Qad", value=1e-02, min=lower_bound, max=max_C_ub, brute_step=1e-05
    )

    standard_init_params.add(
        "Cad", value=1e-03, min=lower_bound, max=max_C_ub, brute_step=1e-05
    )

    standard_init_params.add("nAd", value=0.9, min=n_lb, max=1, brute_step=0.05)
    standard_init_params.add(
        "Rorr", value=150, min=lower_bound, max=max_ub, brute_step=1000
    )
    standard_init_params.add(
        "Aw", value=200, min=lower_bound, max=max_ub, brute_step=1000
    )
    standard_init_params.add(
        "Z0", value=500, min=lower_bound, max=max_ub, brute_step=1000
    )
    standard_init_params.add(
        "tau", value=0.5, min=lower_bound, max=max_ub, brute_step=1000
    )
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    standard_init_params.add_many(
        ("R3", 50, True, lower_bound, max_ub, None, 5),
        ("Q3", 1e-02, True, lower_bound, max_C_ub, None, 1e-05),
        ("n3", 0.5, True, n_lb, 1, None, 0.05),
    )
    #                                  ('sigmaW',200,True, 0.1, 3E3, None, 10))
    standard_init_params = params_extra_setting(
        standard_init_params, EISgr_data_EV, EIS_fit_kwargs
    )
    # Old style:  mod,mod2 = Model(EEC_Singh2015_RQRQR,name='Singh2015_RQRQR'),Model(EEC_ORRpOx,name='Bandarenka_2011_RQRQR')
    #   skip(5,Model(EEC_Bandarenka_Ads,name='Bandarenka2011_RQRW'))

    exclude_models_list = [
        "(Singh2015_RQR)",
        "EEC_RQ_RQ_RQ",
        "Model(EEC_RQ_RQ_RW)",
        "(Randles_RQRQ)",
    ]
    mod_lst_filter = EEC_models_index(exclude=exclude_models_list)
    if "Nicole" in EIS_fit_kwargs.get("input_run", False):
        mod_lst_filter = EEC_models_index(select=["(Randles_RQRQ)"])

        mod_lst_filter = EEC_models_index(
            select=["EEC_RQ_RQ_RQ", "Model(EEC_RQ_RQ_RW)"][-1]
        )
        mod_lst_filter = EEC_models_index(
            select=["Model(EEC_2CPE)"], exclude=["EEC_2CPEpW"]
        )
        mod_lst_filter = EEC_models_index(
            select=["EEC_2CPEpRWs"], exclude=["EEC_2CPEpW"]
        )

        if len(mod_lst_filter) == 1:
            model_set, modname = mod_lst_filter[0][1], mod_lst_filter[0][1].name
        EEC_models_index(select=["Model(EEC_Randles_RWpCPE"])
    #    'EEC_Randles_RWpCPE'
    # TODO check whats wrong with Bandarenka Fit model
    MTHDS = [
        "leastsq",
        "least_squares",
        "nelder",
        "lbfgsb",
        "anneal",
        "powell",
        "cobyla",
        "slsqp",
        "differential_evolution",
        "ampgo",
        "basinhopping",
        "dual_annealing",
        "brute",
    ]
    _best_trial_method = MTHDS[1]
    lstsq_method = MTHDS[1]  #'leastsq'
    _prefit_methods = MTHDS[1:2]
    _lmfit_kwargs = {"least_squares": {"loss": "soft_l1"}}
    wt_check, PRETRIALS_weights, all_PRETRIALS = {}, pd.DataFrame(), pd.DataFrame()
    #%%
    PRETRIALS_weights, _bad_models, all_PRETRIALS, all_PRE_spectra = prefit_test(
        fit_run_arg,
        EIS_fit_kwargs,
        weight_opts,
        EIS_data_KKvalid,
        Z_KKv,
        ang_KKv,
        standard_init_params,
        MTHDS,
        prefit_methods=["leastsq"],
        prefit_models=mod_lst_filter,
    )
    #    for n,r in PRETRIALS_weights.iterrows(): # TODO
    #        make_prefit_frame(EIS_data_KKvalid, r['lmfit_best_trial'],plot = True) # TODO
    all_PRETRIALS_grp = all_PRETRIALS.groupby("Model_EEC")
    #    prefit_test(mod_lst_filter,standard_init_params)
    #%%
    pars_lst, fit_data_valid_lst, outP = [], [], {}
    for mod_indx, model_set, rand_n in mod_lst_filter:  #
        modname = model_set.name
        params_model = model_set.make_params()

        #            print('use previous')
        #            params = InitP1.query('Model_EEC')
        par_options = {}
        for pn in params_model:
            params_model[pn].set(
                value=standard_init_params[pn].value,
                min=standard_init_params[pn].min,
                max=standard_init_params[pn].max,
            )
        try:
            _get_params = pd.DataFrame()
            _get_params, recheck_msg = fitting_recheck_params(
                fit_run_arg, modname, params_model, **EIS_fit_kwargs
            )
            #            _logger.warning(recheck_msg)
            if not _get_params.empty:
                for pn in params_model:
                    if (
                        params_model[pn].vary == True
                        and pd.isna(_get_params[pn].iloc[0]) == False
                        and pn != "Rs"
                    ):
                        params_model[pn].set(value=_get_params[pn].iloc[0])
                _logger.info("Prefit used recheked params")
                run_brute_prefit = False
            else:
                _logger.info(
                    f'Prefit _get_params empty {(str(fit_run_arg.PAR_file), fit_run_arg.E_dc_RHE, fit_run_arg.RPM_DAC)}\n message:"{recheck_msg}'
                )
                run_brute_prefit = True
        except Exception as e:
            run_brute_prefit = True
            _logger.error(f"fitting_recheck_params error: {e}")

        par_options.update(
            {"standard": {"params": params_model, "weights": "lmfit_weights_mod_Y"}}
        )
        if not PRETRIALS_weights.empty:
            pretrial_best_row = PRETRIALS_weights.loc[
                (modname, *PRETRIALS_weights.loc[modname]["MSE"].idxmin()), :
            ]

            pretrial_lmfit = pretrial_best_row["lmfit_best_trial"]
            pretrial_weights = pretrial_best_row["lmfit_Weights"]

            for pn in pretrial_lmfit.params:
                pretrial_lmfit.params[pn].set(
                    min=standard_init_params[pn].min, max=standard_init_params[pn].max
                )

            par_options.update(
                {
                    "pretrial": {
                        "lmfit": pretrial_lmfit,
                        "params": pretrial_lmfit.params,
                        "weights": pretrial_weights,
                    }
                }
            )
        #            for pn in params_model:
        #                params_model[pn].set(value=pretrial_best.params[pn].value)

        pre_options_select = par_options.get("pretrial", par_options["standard"])
        if "pretrial" in par_options.keys():
            best_trial = pre_options_select["lmfit"]
        else:
            best_trial = model_set.fit(
                Z_KKv,
                pre_options_select["params"],
                weights=weight_opts[pre_options_select["weights"]],
                ang=ang_KKv,
                method=_best_trial_method,
            )
        #        make_prefit_frame(EIS_data_KKvalid, best_trial,plot = True) #TODO

        init = model_set.eval(best_trial.params, ang=ang_KKv)
        _logger.info(
            f"""Prefit finished with method {best_trial.method} for {fit_run_arg.PAR_file.stem}
                    at {fit_run_arg.E_dc_RHE:0.2G} V on model {model_set.name},RedChiSqr = {best_trial.redchi:.3G} {best_trial.aic:.4G}"""
        )
        # ADD RANDOM TO PRE-OUT params
        #        InitParam = best_trial.params.copy()
        #        if '__lnsigma' in InitParam.keys():
        #            del InitParam['__lnsigma']
        ##
        #        for par in InitParam.keys():
        #            if InitParam[par].vary and not 'Rs' in par:
        #                _parval = InitParam[par].value
        #                _rpar = random.gauss(_parval,0.02*_parval)
        ##                bounds_mean = (InitParam[par].max-InitParam[par].min)/2
        ##                _rguess = random.gauss(-0.1,0.05) if _parval > bounds_mean else random.gauss(0.1,0.05)
        ##                _rpar = _parval * (1 + _rguess)
        #                while _rpar > 0.9*InitParam[par].max and _rpar < 0.1+0.1*InitParam[par].min:
        #                    _rpar = random.gauss(_parval,0.02*_parval)
        #                InitParam[par].set(_rpar)
        #            elif 'Rs' in par:
        #                InitParam[par].set(value=Rs_guess)
        ##        InitParam = prefit.params
        errRE_init, errIM_init = (init.real - Z_KKv.real) / abs(Z_KKv), (
            init.imag - Z_KKv.imag
        ) / abs(Z_KKv)
        #        Weights = 1/abs(Z_KKv)
        #        if modname in ['Model(Singh2015_RQRWR)','Model(Bandarenka2011_RQRW)']:
        #            InitParam['nAd'].set(value= 0.5, vary=False)
        #            InitParam['Rorr'].set(value= 3E3, vary=False)
        #        for i in InitParam.keys():
        #            InitParam[i].set(value= InitParam[i].value , min= standard_init_params[i].min, max= standard_init_params[i].max)
        weights_used_out = pre_options_select["weights"]
        weight_opts[pre_options_select["weights"]]
        out = model_set.fit(
            Z_KKv,
            best_trial.params,
            ang=ang_KKv,
            weights=best_trial.weights,
            method=lstsq_method,
            **_lmfit_kwargs.get(lstsq_method),
        )
        #        make_prefit_frame(EIS_data_KKvalid, out, prefix = 'pp',plot = 'Y') # TODO
        #        model_set.fit(Z_KKv,best_trial.params,ang=ang_KKv,weights= statiscial_weights, method='leastsq')

        #        fitter = Minimizer(model_set, InitParam, fcn_args=(ang_KKv, Z_KKv))
        #        result_brute = fitter.minimize(method='brute', Ns=25, keep=25)
        #        out = model_set.fit(Z_KKv,InitParam,ang=ang_KKv,weights= statiscial_weights , method= 'brute')

        bf_good_low_high, (out_low_err, out_high_err), MSE_out = make_prefit_frame(
            EIS_data_KKvalid, out
        )
        _logger.info(
            f"""Fit out with method {out.method} for {fit_run_arg.PAR_file.stem} at {fit_run_arg.E_dc_RHE:0.2G} on model \
                    {model_set.name},ChiSqr = {out.chisqr:.3G}, RedChiSqr = {out.redchi:.3G}, {out.aic:.4G}"""
        )
        ### == RETRY SECTION: when first fit is not good enough.. some initial values are re-set and fit is run again ... ####
        retry_attempt = "no"

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
                "lmfit_weights": out.weights,
                "lmfit_weights_nm": weights_used_out,
                "PAR_file": fit_run_arg.PAR_file,
            }
        )
        # add the metadata to the EIS_data_KKvalid_fit DF !!
        EISgr_meta_add_to_fit = pd.DataFrame(
            [EISgr_meta_add] * len(EIS_data_KKvalid_fit),
            index=EIS_data_KKvalid_fit.index,
        )
        #        EISgr_meta_add_to_fit.drop()
        merge_cols = list(
            EISgr_meta_add_to_fit.columns.intersection(EIS_data_KKvalid.columns)
        )
        EIS_data_KKvalid_fit = pd.concat(
            [EIS_data_KKvalid_fit, EISgr_meta_add_to_fit], axis=1
        )
        #        if not EIS_data_KKvalid_BR_specs.empty:
        #            merge_cols = list(EIS_data_KKvalid_BR_specs.columns.intersection(EIS_data_KKvalid.columns))
        #            EIS_data_KKvalid_fit = pd.merge(EIS_data_KKvalid_fit,EIS_data_KKvalid_BR_specs,on =merge_cols,how='inner')
        #            EIS_data_KKvalid_fit =  pd.merge([EIS_data_KKvalid_fit,EIS_data_KKvalid_BR_specs],axis=1)
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
        fit_run_arg_add = {
            EvRHE: fit_run_arg.E_dc_RHE,
            "PAR_file": fit_run_arg.PAR_file,
            "PAR_date": EISgr_meta_combined["PAR_date"],
            "RPM_DAC": fit_run_arg.RPM_DAC,
            "Segment": fit_run_arg.Segment,
            "Model_EEC": modname,
            "Model_index": mod_indx,
        }
        outP.update(fit_run_arg_add)
        #        outP.update(extraP)
        if "Qad" not in outP.keys():
            outP.update({"Qad": 0, "nAd": None})
        xtra = {
            "Rct_kin": outP["Rct"] ** -1,
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
            "test_errcheck_msg": bf_good_low_high,
            "test_low": out_low_err,
            "test_high": out_high_err,
            "test_sum": out_high_err + out_low_err,
            "BRUTE_FIT": 0,
            "FINAL_FIT": 1,
        }
        outP.update(xtra)
        pars_mod_out = pd.DataFrame(outP, index=[EIS_data_KKvalid_fit.index[0]])
        if modname in all_PRETRIALS_grp.groups.keys():
            if not all_PRETRIALS_grp.get_group(modname).empty:
                BR_best_fit_weight_opt_out = all_PRETRIALS_grp.get_group(
                    modname
                ).assign(**{**EISgr_meta_add, **fit_run_arg_add})
            pars_diff_cols = pars_mod_out.columns.difference(
                BR_best_fit_weight_opt_out.columns
            )
            pars_mod_out = pd.concat(
                [pars_mod_out, BR_best_fit_weight_opt_out],
                sort=False,
                ignore_index=True,
            )
        pars_lst.append(pars_mod_out)
        fit_data_valid_lst.append(EIS_data_KKvalid_fit)
    #        lmfit_out_lst.append({'Model_EEC' : modname,'Model_index' : mod_indx,'lmfit_out' : out})
    pars_models = pd.concat(pars_lst, sort=False, ignore_index=True)
    EIS_fit_data = pd.concat(fit_data_valid_lst, sort=False, ignore_index=True)
    #%%
    #    EIS_fit_data.groupby('Model_index').plot(x=['DATA_Yre','FIT_Yre'],y=['DATA_Yim','FIT_Yim'],kind='scatter')
    # +++ PLOTTING of EIS spectrum +++
    vers = FileOperations.EIS_version
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
        EIS_fit_data, spectra_fit_outpath
    )
    spectra_raw_outpath_target = FileOperations.CompareHashDFexport(
        EIS_data_raw, spectra_raw_outpath
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
            EIS_fit_data,
            pars_models.loc[pars_models.BRUTE_FIT == 0],
            spectra_fit_outpath_png,
            plot_show=False,
            std_print_model="EEC_2CPEpRW",
        )
    return (EIS_fit_data, pars_models)


#        index_dataplot_output = {'PAR_file': nm2,'Type_output' : 'EIS_fit_data_png', 'Type_exp' : 'EIS', 'DestFile' : EIS_outPath_target_png, 'E_V' : fit_run_arg.E_dc_RHE,'RPM_DAC' : fit_run_arg.RPM_DAC}
#        indexes_per_EV_out.append(index_dataplot_output)
#    index_data_output = {'PAR_file': fit_run_arg.PAR_file,'Type_output' : 'EIS_fit_data', 'Type_exp' : 'EIS', 'DestFile' : EIS_fit_data_outPath_target, 'E_V' : fit_run_arg.E_dc_RHE,'RPM_DAC' : fit_run_arg.RPM_DAC}
#    indexes_per_EV.append(index_data_output)
#    if EIS_fit_kwargs.get('export_raw_data',False):


#        index_raw_output = {'PAR_file': fit_run_arg.PAR_file,'Type_output' : 'EIS_raw_data', 'Type_exp' : 'EIS', 'DestFile' : EIS_outpath_target_raw,'E_V' : fit_run_arg.E_dc_RHE}
#        indexes_per_EV.append(index_raw_output)

#    meta_export = meta_export_templ(*fit_run_arg, spectra_fit_outpath_target, spectra_raw_outpath_target,pars_outpath_target)
#    meta_export_templ = namedtuple('meta_export_templ', fit_run_arg._fields + ('File_SpecFit','File_SpecRaw','File_Pars'))
#    fit_export = fit_export_templ(EIS_fit_data,pars_models)
#    indexes_per_EV_out = pd.DataFrame(indexes_per_EV)

# OLD MODELS
#  mod_lst_filter = EEC_models_index(select=['Singh2015_RQRQR'])
#        mod_lst_filter = EEC_models_index(select=['Singh2015_RQRWR'])
#        mod_lst_filter = EEC_models_index(select=['Singh2015_R3RQ'])
#        mod_lst_filter = EEC_models_index(select=['Singh2015_RQRQR'])
#        mod_lst_filter = EEC_models_index(select=['EEC_Singh2015_RQRQRW'])
#        mod_lst_filter = EEC_models_index(select=['EEC_Singh2015_RQRQRserW'])
#        mod_lst_filter = EEC_models_index(select=['EEC_Singh2015_RQRWQR'])
#        mod_lst_filter = EEC_models_index(select=['EEC_Singh2015_RQRQRW',
#                                                  'EEC_Singh2015_RQRQRserW',
#                                                  'EEC_Singh2015_RQRWQR'])
