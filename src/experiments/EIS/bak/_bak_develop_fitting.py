"""
Created on Sun Jul 11 10:31:48 2021

@author: DW
"""

# print(f"Name: {__name__} for file {__file__}")

# if __name__ == "__main__":
#     from validation import get_KKvalid
#     from models import EEC_models_index, Model_Collection

#     #    from .plotting import plot_linKK, EIS_Trimming_plot
#     #    import ECpy
#     from ECpy.experiments.EIS.plotting import (
#         plot_linKK,
#         EIS_Trimming_plot,
#         EIS_plotting_per_EV,
#         plot_lin_Warburg,
#     )
#     from ECpy.experiments.EIS.eis_run_ovv import EIS_Spectrum
#     from DRT_DP_fitting import DP_DRT_analysis
#     from GP_DRT_fitting import run_GP_DRT_fit

#     #    sys.path.append(Path(__file__).parent.parent.parent.parent)
#     from EC_logging_config import start_logging

# #    spec = importlib.util.spec_from_file_location('validation',_fp)
# #    module = importlib.util.module_from_spec(spec)
# #    spec.loader.exec_module(module)

# #    from impedance import validation as imppy_validaton

#     #    from ECpy.experiments.EIS.eis_run_ovv import EIS_Spectrum
#     from DRT_DP_fitting import DP_DRT_analysis
#     from GP_DRT_fitting import run_GP_DRT_fit

#    from EC_logging_config import start_logging
#    sys.path.append(Path(__file__).parent.parent.joinpath('runEC'))
#    from runEC.EC_logging_config import start_logging
#    _logger = start_logging(__name__)


def get_from_good_fits():
    _good_fits_par_fit = get_best_pars_from_good_fits(
        Z_KKv, ang_KKv, params_model, standard_init_params, MTHDS, modname, model_set
    )  # TODO
    make_prefit_frame(EIS_data_KKvalid, _good_fits_par_fit[3], plot=True)
    br_eval_fit_DF.loc[br_eval_fit_DF[0] < br_eval_fit_DF[0].mean()].iloc[0].to_numpy()
    if not _good_fits_par_fit.empty:
        #            _good_params = _good_fits_par_fit[-2]
        for _Rsnm, _Rsgrp in _good_fits_par_fit.groupby("Rs_setting"):
            _Rsnm, _Rsgrp
            _Rsgrp = _Rsgrp.sort_values(by=["delta_Rs", "_chi"])
            par_options.update({_Rsnm: _Rsgrp.iloc[0]["lmfit_params"]})


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
# class EIS_loop():
#      def __init__(self,fit_run_args):
#         self.fit_run_args = fit_run_args
#    def fit_mean_PAR_file(self):
#        for PF in itertools.groupby(self.fit_run_args, lambda x: x.PAR_file):


def plot_overE():
    E_data_combined_path_target = FileOperations.CompareHashDFexport(
        PF_fit_spectra_all_grp.get_group(PF), dest_PF["EIS_dest_spectra"].iloc[0]
    )
    EIS_plotting_EvRHE(
        PF_fit_spectra_all_grp.get_group(PF),
        self._pars,
        PF_pars,
        dest_PF["EIS_dest_spectra"].iloc[0],
    )


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


#%% OLD STUFF
def _init_Models(self):
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


def old_stuff():
    pass
    # self.standard_init_params['Rion'].value = self.Rion_guess
    # self.Rs_guess = self._obj.Rs_guess
    #     standard_init_params.add('Ls', value= 2E-05, min= 1E-10, max=1,brute_step = 1E-03)
    #     standard_init_params.add('Cdlp', value= 1E-05,  min=lower_bound, max=max_C_ub, brute_step = 0.5E-06)
    #     standard_init_params.add('nDL', value= 0.60, min=n_lb, max=1, brute_step = 0.05)
    #     standard_init_params.add('Rct', value= 25, min=lower_bound, max=max_ub,vary=True,brute_step=10)
    # #    params.add('Aw', value= 3000, min=5.0, max=1E6)
    # #    params.add('Rad', value= 5E5, min=0.01, max=1E7)
    #     standard_init_params.add('Qad', value= 1E-03,  min=lower_bound, max=max_C_ub, brute_step = 1E-05)
    #     standard_init_params.add('Cad', value= 1E-04,  min=lower_bound, max=max_C_ub, brute_step = 1E-05)
    #     standard_init_params.add('nAd', value= 0.9, min= n_lb, max=1, brute_step = 0.05)
    #     standard_init_params.add('Rorr', value= 150, min=lower_bound, max=max_ub, brute_step=1000)
    #     standard_init_params.add('Aw', value= 100, min=lower_bound, max=max_ub, brute_step=1000)
    #     standard_init_params.add('tau', value= 0.5, min=lower_bound, max=2E2, brute_step=1000)
    #     standard_init_params.add('rho_el', value= 10, min=lower_bound, max=1E8, brute_step=1000)
    #     standard_init_params.add('n_pores', value= 1E10, min=1, max=1E20, brute_step=1000)
    #     standard_init_params.add('l_pore', value= 1E-3, min=lower_bound, max=10, brute_step=1000)
    #     standard_init_params.add('r_pore', value= 1E-6, min=1E-10, max=1E-1, brute_step=1000)
    #     standard_init_params.add('Rel_p', value= 50, min=lower_bound, max=1E8, brute_step=1000)
    #     standard_init_params.add('t_G', value= 20, min=lower_bound, max=1E8, brute_step=1000)
    #     standard_init_params.add('R_G', value= 2E3, min=lower_bound, max=1E8, brute_step=1000)
    #     standard_init_params.add('phi_G', value= 1, min=lower_bound, max=1E8, brute_step=1000)
    # standard_init_params.add('Aws', value= 200, min=lower_bound, max=max_ub, brute_step=1000)
    # standard_init_params.add('Awo', value= 200, min=lower_bound, max=max_ub, brute_step=1000)
    # standard_init_params.add('Z0', value= 500, min=lower_bound, max=max_ub, brute_step=1000)
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    # standard_init_params.add_many(('R3', 50,True, lower_bound, max_ub, None, 5),
    #                               ('Q3', 1E-02,True, lower_bound, max_C_ub, None, 1E-05),
    #                               ('n3', 0.5,True, n_lb, 1, None, 0.05))
    #                                  ('sigmaW',200,True, 0.1, 3E3, None, 10))
    # self.standard_init_params = standard_init_params
    # def validate_params(self):
    #     _missing_pars = []
    #     if not self.Model_Collection.params_set.issubset(set(self.standard_init_params.keys())):
    #         _missing_pars = [i for i in self.Model_Collection.params_set if i not in set(self.standard_init_params.keys())]
    #         # print(f'missing params, {_missing_pars}')
    #     self._missing_pars = _missing_pars
    # def set_standard_init_from_guesses(self):

    #     _guesses = [i.guesses for i in self.Model_Collection.model_selection]
    #     _check = [mp for mp in self._missing_pars if any(mp in _g.keys() for _g in _guesses)]
    #     if _check:
    #         for mp in self._missing_pars:
    #             mp_g = [i[mp] for i in _guesses if mp in i.keys()]
    #             if mp_g:
    #                 # print(mp_g)
    #                 self.standard_init_params.add(mp, value=mp_g[0],min = 1E-02*mp_g[0], max = 1E2*mp_g[0] )
    def lin_Rion():

        ax.plot(np.arange(0, x1), np.arange(0, x1))
        xrange = np.arange(0, x1)[0 : len(_y)]
        tang = np.array([x + x0 for x in xrange])

        _step = 0
        _index = []
        # fig,ax =plt.subplots1()
        while bool(_index) == False and x0 < _x.max() * 2:
            tang = np.array([x + x0 for x in xrange])
            _y_tang = tang - tang.min()
            ax.plot(tang, _y_tang, alpha=0.2)
            _diff = np.array([i[0] - i[1] for i in zip(_y.values, _y_tang)])

            _sumdiff = np.sqrt(np.abs(_diff)).sum()
            if x0 == 0:
                _newdiff = _sumdiff
            _close_diff = _diff[np.abs(_diff) < 1]
            if _sumdiff < _newdiff:
                _newdiff = _sumdiff
                _xbest = x0

            if _close_diff.size > 0:
                print(f"x0: {_x}", _close_diff)
            # np.all(np.diff(_close_diff.index) < 3)
            if not _close_diff.size > 0 and len(_close_diff) >= 2:
                _index = list(_close_diff.index)
                _tangentres.append([curve, _index])

            else:
                x0 += 0.05
                # if curve == 'lf':
                # _xadd = 0.5
                # elif curve == 'hf':
                # _xadd = +0.5

    def _func(x0, y, xrange):
        tang = np.array([x + x0 for x in xrange])
