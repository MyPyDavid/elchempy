#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:29:50 2020

@author: zmg
"""


def fit_EEC(fit_run_arg, **EIS_fit_kwargs):
    #%%

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
        "Qad", value=1e-04, min=lower_bound, max=max_C_ub, brute_step=1e-05
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
            # _get_params, recheck_msg = fitting_recheck_params(fit_run_arg,modname,params_model, **EIS_fit_kwargs) # TODO
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
