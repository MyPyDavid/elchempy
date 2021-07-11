"""
Created on Sun Jul 11 11:35:02 2021

@author: DW
"""


# TODO continue means:
# finished alread in other class??
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


def filter_error_debug():
    "12.02.2019_0.1MKOH_cell2/O2_EIS-range_1500rpm_JOS6_postAST_899 at 0.999 (1500 rpm)"
    #    FileOperations.CompareHashDFexport(EISgrEvRHE_actions.T,EIS_ovv_destfiles.EIS_dest_dir.iloc[0].joinpath('PAR_actions_T.xlsx'))
    bad_fit = ("NSS_0103-EIS-7", 3)
    tt = "N2_EIS-range_1500rpm_JOS4_288"
    fit_run_arg = [i for i in fit_run_args if tt in i[0].name and 558 == i.E_dc_RHE_mV]
    fit_run_args_tests = [i for i in fit_run_args if tt in i[0].name]


def eis_run_group_ovv(ovv_exp_grp, **EIS_kwargs):

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
            # pool_size = os.cpu_count()-2
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


def testing_mode():
    er = EIS_Preparator(test)
    erter = iter(er)
    es = next(erter)
    fsp = Fit_Spectra_Collection(es, run_fit_mean="lmfit mean")
    lm = LMfit_method(fsp.fit_mean)
    fsp_all = Fit_Spectra_Collection(es, run_fit_mean="lmfit")
