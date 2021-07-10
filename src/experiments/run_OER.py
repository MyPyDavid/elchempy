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


if __name__ == "__main__":
    pass

import logging

logger = logging.getLogger(__name__)


Meta = namedtuple("Meta", "PAR_file data ovv N2_BG_f N2_BG_data")


def ORR_prepare_ovv_dest(ovv_exp_grp, **OER_kwargs):
    gr_ovv = ovv_exp_grp.get_group("OER")

    #    gr_ovv = gr_ovv.astype({'ORR_act_N2_bg' : str})

    gr_ovv = gr_ovv.assign(
        **{
            "OER_dest_dir": [
                Path(i).joinpath(f"OER_v{FileOperations.version}")
                for i in gr_ovv.Dest_dir.to_numpy()
            ]
        }
    )
    #    for _dest_dir in gr_ovv.ORR_dest_dir.unique()
    _mkdirs = [
        i.mkdir(parents=True, exist_ok=True) for i in gr_ovv.ORR_dest_dir.unique()
    ]

    gr_ovv_disk = gr_ovv.loc[
        (gr_ovv.Comment.str.contains("Pt|chrono") == False)
        & (gr_ovv.SampleID.str.contains("PT-RING") == False)
    ]

    _calc_dd = []
    for n, r in gr_ovv_disk.iterrows():
        _dd = r.ORR_dest_dir
        _ORR_calc_dest_dir = _dd.joinpath(
            f'{r["Electrolyte"]}_{r["SampleID"]}_{r["Loading_name"]}_{r["postAST"]}'
        )
        _ORR_calc_dest_dir.mkdir(parents=True, exist_ok=True)
        _fstem = FileOperations.Check_fstem(Path(r.PAR_file))
        _calc_dd.append((_ORR_calc_dest_dir, _fstem))

    gr_ovv_disk = gr_ovv_disk.assign(
        **{
            "ORR_ecexp_destdir": [i[0] for i in _calc_dd],
            "PAR_fstem": [i[1] for i in _calc_dd],
        }
    )
    #    EC_index = ORR_kwargs.get('EC_index', pd.DataFrame())
    #    if not EC_index.empty:
    #        _add_ring_ovv = EC_index.loc[EC_index.PAR_date.isin(gr_ovv_disk.PAR_date.unique())]
    #        gr_ovv

    return gr_ovv, gr_ovv_disk


def OER(ovv_exp_grp, **OER_kwargs):
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

    run_args_raw = [
        Meta(
            Path(pf),
            grp,
            ovv_all.loc[ovv_all.PAR_date_day.isin(grp.PAR_date_day.unique())],
            *ORR_get_N2_BG(ovv_all, pf, grp),
        )
        for pf, grp in gr_ovv_disk.groupby(by="PAR_file")
    ]

    _n2faillst = [i[0] for i in orr_run_args_raw if i[-1].empty]

    orr_run_args = [i for i in orr_run_args_raw if not i[-1].empty]

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
        orr_kwargs_iter = repeat(ORR_kwargs)
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


def OER(exp, gr_ovv, ovv):
    ###### === Analyze the OER experiments ==== #######
    #        DestDir = Path(gr_ovv.Dest_dir.unique()[0])
    #        print('OER skipped')
    #        return pd.DataFrame()
    gr_ovv = ovv.groupby(by="PAR_exp").get_group("OER").query('Electrode != "Pt_ring"')
    DestDir = Path(gr_ovv.Dest_dir.unique()[0])
    #        gr_ovv = ExpTypes_gr
    OER_index_lst, OER_pars, faillst = [], [], []
    for OER_file, OER_ovv_file in failed_ovv.groupby(by="PAR_file"):
        try:
            #                ring = ovv.loc[ovv.PAR_date.isin(OER_ovv_file.PAR_date.values) & ovv.Electrode.str.contains('Pt_ring' )]
            ring = ovv.loc[
                (ovv.PAR_date - OER_ovv_file.PAR_date.values[0] < pd.Timedelta("10s"))
                & (ovv.Electrode.str.contains("Pt_ring"))
            ]
            if not ring.empty:
                rr = ring

            OER_CVs, OER_actions = create_CVs(OER_ovv_file)
            if not OER_CVs.empty:
                Samples_ovv_cv = OER_CVs[
                    (~OER_CVs["SampleID"].str.contains("Pt_ring|Pt-ring"))
                    & OER_CVs.Type.str.contains("Voltammetry")
                ]
                if len(OER_CVs) - len(Samples_ovv_cv) > 0:
                    logger.info(
                        "OER filtered samples from CV files {0}".format(
                            len(OER_CVs) - len(Samples_ovv_cv)
                        )
                    )
                if not Samples_ovv_cv.empty:
                    OER_file_index = OER_scan(
                        Samples_ovv_cv,
                        OER_ovv_file,
                        Path(OER_ovv_file.Dest_dir.iloc[0]),
                    )
                    OER_index_lst.append(OER_file_index)
                else:
                    logger.error(
                        "OER === ERROR OER_scan Samples ovv empty === {0}\n because {1}".format(
                            OER_file
                        )
                    )
                    faillst.append([OER_file, "Samples_ovv_cv.empty"])
            else:
                logger.error(
                    "OER === ERROR OER_CVs empty === {0}\n because {1}".format(OER_file)
                )
                faillst.append([OER_file, "OER_CV.empty"])
        except Exception as e:
            logger.error(
                "OER === ERROR in OER_scan === {0}\n because {1}".format(OER_file, e)
            )
            faillst.append([OER_file, "OER exception {0}".format(e)])
    failed_ovv = gr_ovv.loc[gr_ovv.PAR_file.isin([i[0] for i in faillst])]
    OER_index = pd.concat([i for i in OER_index_lst], ignore_index=True)
    OER_pars_target = FolderOps.FileOperations.CompareHashDFexport(
        OER_index, DestDir.joinpath("OER_Pars_index.xlsx")
    )
    return OER_index


#        OER_pars.to_excel(DestDir.joinpath('OER_Pars.xlsx'))
#        index = pd.DataFrame([['OER',OER_pars_target,i] for i in OER_index.PAR_file.unique()],columns=['Type_output','DestFile','PAR_file'])
#        if 'OERtest' in Experiments:
#            print('OER')
##                OER_CVs,OER_actions = create_CVs(ovv.query('PAR_exp == "OER"'),PathDB,False)
#                if not OER_CVs.empty:
#                    OER_pars = OER_scan(OER_CVs,dest_dir)
#                else:
#                    print('OER failed: %s'%exp_dir)
#        OER_CVs,HER_info = pd.DataFrame(), pd.DataFrame()
#        OER_ovv = ovv.loc[((ovv['PAR_exp'] == "OER") & (ovv.basename.str.contains('OER'))),:]
##                    HER_CVs,HER_info = create_CVs(HER_ovv,PathDB,True)
#        OER_dest_dir = dest_dir.joinpath('OER')
#        OER_dest_dir.mkdir(parents=True,exist_ok=True)
#        OER_out_fn = OER_dest_dir.joinpath('OER_Pars.xlsx')
#        OER_pars = OER_calc(OER_ovv,OER_out_fn,PathDB)

#                OER_pars,OER_action = OER_scan(create_CVs(ovv.query('PAR_exp == "OER"'),PathDB,False),dest_dir)
