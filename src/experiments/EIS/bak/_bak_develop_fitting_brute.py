"""
old functions for finding initial parameter values in a 'brute' way for fitting
"""


# TODO Classify get_best_pars
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


#    if all([i in params_model for i in ['Rct','Rorr']]):
#        if params_model['Rct'] < params_model['Rorr']:
#            _Rorr = params_model['Rct']
#            params_model['Rct'].set(value=params_model['Rorr'])
#            params_model['Rorr'].set(value=_Rorr )
# TODO finish brute preparation
def prepare_rand_models(total_combis=5e3, rand_params_dest_dir=""):
    mod_lst_filter = Model_Collection().lmfit_models

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
        for i, mod, *_ in mod_lst_filter:
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
