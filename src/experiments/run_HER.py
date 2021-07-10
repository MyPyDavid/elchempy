import sys
from pathlib import Path
from collections import namedtuple
import datetime as dt
import numpy as np
import re
import os
import multiprocessing
from functools import partial
from itertools import repeat
import pandas as pd


from .EC_DataLoader.CreateCV import create_CVs
from .HER.HER_analyze_scans import HER_scan

if __name__ == "__main__":
    pass
#    from EC_logging_config import start_logging
#    from .analyze_scans import N2_scans
import logging

logger = logging.getLogger(__name__)


Meta = namedtuple("Meta", "PAR_file data ovv")


def PF_fit_starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(PF_fit_apply_args_and_kwargs, args_for_starmap)


def PF_fit_apply_args_and_kwargs(fn, args, kwargs):
    return fn(args, **kwargs)


def HER_prepare_ovv_dest(ovv_exp_grp, dest_dir_name="HER_dest_dir", **HER_kwargs):
    #    gr_ovv = ovv_exp_grp.get_group('HER')
    gr_ovv = pd.concat([gr for n, gr in ovv_exp_grp if "HER" in n])
    #    gr_ovv = gr_ovv.astype({'ORR_act_N2_bg' : str})
    gr_ovv = gr_ovv.assign(
        **{
            dest_dir_name: [
                Path(i).joinpath(f"HER_v{FileOperations.version}")
                for i in gr_ovv.Dest_dir.to_numpy()
            ]
        }
    )
    #    for _dest_dir in gr_ovv.ORR_dest_dir.unique()
    _mkdirs = [
        i.mkdir(parents=True, exist_ok=True) for i in gr_ovv[dest_dir_name].unique()
    ]
    gr_ovv_disk = gr_ovv.loc[
        (gr_ovv.Comment.str.contains("Pt|chrono") == False)
        & (gr_ovv.SampleID.str.contains("PT-RING") == False)
    ]
    _calc_dd = []
    for n, r in gr_ovv_disk.iterrows():
        _dd = r[dest_dir_name]
        _calc_dest_dir = _dd.joinpath(
            f'{r["Electrolyte"]}_{r["SampleID"]}_{r["Loading_name"]}_{r["postAST"]}'
        )
        _calc_dest_dir.mkdir(parents=True, exist_ok=True)
        _fstem = FileOperations.Check_fstem(Path(r.PAR_file))
        _calc_dd.append((_calc_dest_dir, _fstem))

    gr_ovv_disk = gr_ovv_disk.assign(
        **{
            f"{dest_dir_name}_ecexp": [i[0] for i in _calc_dd],
            "PAR_fstem": [i[1] for i in _calc_dd],
        }
    )

    HER_kwargs.update({"dest_dir_col": dest_dir_name})
    return gr_ovv, gr_ovv_disk, HER_kwargs


def HER(ovv_exp_grp, **HER_kwargs):
    #%%
    ###### === Analyze the ORR experiments ==== #######
    #        exp,gr_ovv,ovv = 'ORR',ExpTypes_gr.get_group('ORR'),ovv
    #        if 'O2' in Gases and not 'N2_act' in Experiments:
    #            print('N2 BACKGROUND SCAN IS MISSING FOR THIS ORR EXPERIMENT.\nFailed: %s'%exp_dir)
    #            sORR = 1
    #        elif 'O2' in Gases and 'N2_act' in Experiments: # O2 and N2 changed
    #        O2_activity, O2_out_ovv, O2_Jcalc_ovv, overall_rows   = pd.DataFrame([]),pd.DataFrame([]),pd.DataFrame([]), pd.DataFrame([])
    #    index_ORR, index_out, faillst = [], [], []
    #    ovv_all, gr_ovv_disk = ORR_prepare_ovv_dest(ovv_exp_grp, **ORR_kwargs)
    ovv_all, gr_ovv_disk, HER_kwargs = HER_prepare_ovv_dest(ovv_exp_grp, **HER_kwargs)
    run_args_raw = [
        Meta(
            Path(pf),
            grp,
            ovv_all.loc[ovv_all.PAR_date_day.isin(grp.PAR_date_day.unique())],
        )
        for pf, grp in gr_ovv_disk.groupby(by="PAR_file")
    ]

    #    _n2faillst = [i[0] for i in orr_run_args_raw if i[-1].empty]

    #    orr_run_args = [i for i in orr_run_args_raw if not i[-1].empty]

    #%%
    if gr_ovv_disk.empty:
        return logger.warning(
            "ORR attemp failed for {0} because ovv empty".format(
                gr_ovv.Dest_dir.unique()
            )
        )
    #    logger.warning('ORR disk OVV empty')
    #%%
    multi_par_fit = True
    if multi_par_fit:
        fit_export_EV_all = []
        kwargs_iter = repeat(HER_kwargs)
        try:
            logger.info(
                f"ORR orr_run_group_ovv START multiprocessing {multi_par_fit} for len{len(run_args_raw)}"
            )
            pool_size = os.cpu_count() - 2
            if "linux" in sys.platform:
                #                os.system('taskset -cp 0-%d %s' % (pool_size, os.getpid()))
                os.system("taskset -p 0xff %d" % os.getpid())
            #                os.sched_setaffinity(0,{i for i in range(pool_size)})
            #            for chunk in orr_run_args[:pool_size if pool_size > 2 else 1]:
            with multiprocessing.Pool(pool_size) as pool:
                PF_fit_starmap_with_kwargs(pool, HER_scan, run_args_raw, kwargs_iter)
        #                    fit_export_EV_all.append(fit_export_EV_all_chunck)
        except Exception as e2:
            logger.error(
                f"HER_scan run multiprocessing error: {e2}, len out({len(fit_export_EV_all)})"
            )

    else:
        fit_export_EV_all = []
        for fit_run_arg in run_args_raw:
            try:
                #                fit_run_arg = [i for i in run_args_raw if i[0].stem.startswith('N2_')][1]
                HER_scan(fit_run_arg, **HER_kwargs)
            #                fit_export_EV_all.append([fit_export_EV])
            except Exception as e:
                logger.warning(f"HER attemp failed for {fit_run_arg[0]} because {e}")


def old_HER_run(ovv_exp_grp, **HER_kwargs):
    ####### === Analyze the HER experiments ==== #######
    #        DestDir = Path(gr_ovv.Dest_dir.unique()[0])
    #        print('OER skipped')
    #        return pd.DataFrame()
    gr_ovv = ovv.groupby(by="PAR_exp").get_group("HER")
    DestDir = Path(gr_ovv.Dest_dir.unique()[0])
    #        gr_ovv = ExpTypes_gr
    HER_index_lst, HER_pars, faillst = [], [], []
    for HER_file, HER_ovv_file in gr_ovv.groupby(by="PAR_file"):
        try:
            HER_CVs, HER_actions = create_CVs(HER_ovv_file)
            if not HER_CVs.empty:
                Samples_ovv_cv = HER_CVs[
                    (~HER_CVs["SampleID"].str.contains("Pt_ring|Pt-ring"))
                ]
                if len(HER_CVs) - len(Samples_ovv_cv) > 0:
                    logger.info(
                        "HER filtered samples from Pt-ring files {0}".format(
                            len(HER_CVs) - len(Samples_ovv)
                        )
                    )
                if not Samples_ovv_cv.empty:
                    HER_file_index = HER_scan(
                        Samples_ovv_cv,
                        HER_ovv_file,
                        Path(HER_ovv_file.Dest_dir.iloc[0]),
                    )
                    HER_index_lst.append(HER_file_index)
                else:
                    logger.error(
                        "HER === ERROR HER_scan Samples ovv empty === {0}\n because {1}".format(
                            HER_file
                        )
                    )
                    faillst.append([HER_file, "Samples_ovv_cv.empty"])
            else:
                logger.error(
                    "HER === ERROR HER_CVs empty === {0}\n because {1}".format(HER_file)
                )
                faillst.append([OER_file, "OER_CV.empty"])
        except Exception as e:
            logger.error(
                "HER === ERROR in HER_scan === {0}\n because {1}".format(HER_file, e)
            )
    HER_index = pd.concat([i for i in HER_index_lst], ignore_index=True)
    HER_pars_target = FolderOps.FileOperations.CompareHashDFexport(
        HER_index, DestDir.joinpath("HER_Pars_index.xlsx")
    )
    return HER_index


#    @staticmethod
#    def HER(exp,gr_ovv,ovv):

##        if 'HER' in Experiments:
##            print('HER')
#        HER_CVs,HER_info = pd.DataFrame(), pd.DataFrame()
#        HER_ovv = ovv.loc[((ovv['PAR_exp'] == "HER") & (ovv.basename.str.contains('HER'))),:]
#        if HER_ovv.empty:
#            HER_ovv = ovv.loc[((ovv['PAR_exp'] == "EIS") & (ovv.basename.str.contains('HER'))),:]
##                    HER_CVs,HER_info = create_CVs(HER_ovv,PathDB,True)
#        HER_dest_dir = dest_dir.joinpath('HER')
#        HER_dest_dir.mkdir(parents=True,exist_ok=True)
#        HER_out_fn = HER_dest_dir.joinpath('HER_Pars.xlsx')
#        HER_pars = HER_calc(HER_ovv,HER_out_fn,PathDB)
#                if not HER_CVs.empty:
#                    HER_pars = H$R_scan(HER_CVs,dest_dir)
#                else:
#                    print('OER failed: %s'%exp_dir)
#                OER_pars,OER_action = OER_scan(create_CVs(ovv.query('PAR_exp == "OER"'),PathDB,False),dest_dir)

#    else:
#        plt.close()
#        print('No ORR')
#        sORR = 0
