"""
Created on Sun Jul 11 10:56:48 2021

@author: DW
"""


def CreateCV_multi(ovv, OutPutRaw=False):

    try:
        logger.info(
            f"EIS eis_run_group_ovv START multiprocessing {multi_par_fit} for len{len(fit_run_args)}"
        )
        pool_size = os.cpu_count() - 2
        #            if 'linux' in sys.platform:
        #                os.system("taskset -p 0xff %d" % os.getpid())
        with multiprocessing.Pool(pool_size) as pool:
            fit_export_EV_all_chunck = PF_fit_starmap_with_kwargs(
                pool, fit_EEC, fit_run_args, fit_kwargs_iter
            )
            fit_export_EV_all.append(fit_export_EV_all_chunck)
    except Exception as e2:
        #        print('EIS eis_run_group_ovv multiprocessing error: {0}'.format(e2))
        logger.error(
            f"EIS eis_fit_PAR_file  multiprocessing erroe: {e2}, len out({len(fit_export_EV_all)})"
        )
