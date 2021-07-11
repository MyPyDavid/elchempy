"""
Created on Sun Jul 11 11:21:52 2021

@author: DW
"""


# if __name__ == '__mainee__':
#    df = pd.DataFrame([[1, 2, 3, 1], [1, 2, 3, 1], [4, 5, 6, 2], [7, 8, 9, 2]], columns=['A', 'B', 'C', 'cluster_id'])
#    df = df.set_index('cluster_id')
#    out = mp_apply(df, my_func)
#    print(out)

# if __name__ == '__main__':
#    list_of_args2 = [1, 2, 3]
#    just_a_dict = {'key1': 'Some value'}
#    with multiprocessing.Pool(processes=3) as pool:
#        results = pool.map(partial(test_func, 'This is arg1', **just_a_dict), list_of_args2)
#    print(results)


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


def test_func(arg1, arg2, **kwargs):
    print(arg1)
    print(arg2)
    print(kwargs)
    return arg2


def run_OER(ovv_exp_grp, ORR_get_N2_BG, **OER_kwargs):
    #%%
    ###### === Analyze the ORR experiments ==== #######
    #        exp,gr_ovv,ovv = 'ORR',ExpTypes_gr.get_group('ORR'),ovv
    #        if 'O2' in Gases and not 'N2_act' in Experiments:
    #            print('N2 BACKGROUND SCAN IS MISSING FOR THIS ORR EXPERIMENT.\nFailed: %s'%exp_dir)
    #            sORR = 1
    #        elif 'O2' in Gases and 'N2_act' in Experiments: # O2 and N2 changed
    #        O2_activity, O2_out_ovv, O2_Jcalc_ovv, overall_rows   = pd.DataFrame([]),pd.DataFrame([]),pd.DataFrame([]), pd.DataFrame([])
    #    index_ORR, index_out, faillst = [], [], []
    ovv_all, gr_ovv_disk = ORR_prepare_ovv_dest(ovv_exp_grp, **OER_kwargs)

    run_args_raw = [
        Meta(
            Path(pf),
            grp,
            ovv_all.loc[ovv_all.PAR_date_day.isin(grp.PAR_date_day.unique())],
            *ORR_get_N2_BG(ovv_all, pf, grp),
        )
        for pf, grp in gr_ovv_disk.groupby(by="PAR_file")
    ]

    _n2faillst = [i[0] for i in run_args_raw if i[-1].empty]

    orr_run_args = [i for i in run_args_raw if not i[-1].empty]

    #%%
    if gr_ovv_disk.empty:
        return logger.warning(
            "OER attemp failed for {0} because ovv empty".format(
                ovv_exp_grp.Dest_dir.unique()
            )
        )
    #    logger.warning('ORR disk OVV empty')
    #%%
    multi_par_fit = True
    if multi_par_fit:
        fit_export_EV_all = []
        orr_kwargs_iter = repeat(OER_kwargs)
        try:
            logger.info(
                f"ORR orr_run_group_ovv START multiprocessing {multi_par_fit} for len{len(orr_run_args)}"
            )
            pool_size = os.cpu_count() - 2
            if "linux" in sys.platform:
                #                os.system('taskset -cp 0-%d %s' % (pool_size, os.getpid()))
                os.system("taskset -p 0xff %d" % os.getpid())
            #                os.sched_setaffinity(0,{i for i in range(pool_size)})
            #            for chunk in orr_run_args[:pool_size if pool_size > 2 else 1]:
            with multiprocessing.Pool(pool_size) as pool:
                PF_fit_starmap_with_kwargs(
                    pool, Jkin_calc_multi, orr_run_args, orr_kwargs_iter
                )
        #                    fit_export_EV_all.append(fit_export_EV_all_chunck)
        except Exception as e2:
            logger.error(
                f"ORR run multiprocessing erroe: {e2}, len out({len(fit_export_EV_all)})"
            )

    else:
        fit_export_EV_all = []
        for fit_run_arg in orr_run_args:

            try:
                #                fit_run_arg = orr_run_args[1]
                Jkin_calculations(fit_run_arg, **ORR_kwargs)
            #                fit_export_EV_all.append([fit_export_EV])
            except Exception as e:
                logger.warning(f"ORR attemp failed for {fit_run_arg[0]} because {e}")
