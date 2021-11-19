# import sys
from pathlib import Path
from collections import namedtuple
from datetime import datetime
import numpy as np
import logging

from scipy.stats import linregress

# from scipy.interpolate import UnivariateSpline
from scipy import constants

# from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

# import os
# import multiprocessing
# from functools import partial
# from itertools import repeat
import pandas as pd

print("File", __file__, "\nName;", __name__)
# from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations


if __name__ == "__main__":

    from elchempy.experiments.Loader.CreateCV import create_CVs

    # from elchempy.experiments.EC_conditions.electrode import WE_SA_collection_eff
    from O2_chrono import O2_Chrono
    from plotting import ORR_plot_ring_disk
    from KL_calc import KL_plots
#    from EC_logging_config import start_logging
#    pass
else:
    from ..EC_DataLoader.CreateCV import create_CVs as create_CVs

    # from .EC_conditions.electrode import WE_SA_collection_eff
    from .O2_chrono import O2_Chrono
    from .plotting import ORR_plot_ring_disk
    from .KL_calc import KL_plots
#    from ..Loader.CreateCV import create_CVs as create_CVs


_logger = logging.getLogger(__name__)
#
# All_ovv.loc[(All_ovv['Elapsed Time(s)'].isin(Te['Elapsed Time(s)'])) &  (All_ovv['Type'] =='Chronoamperometry') & (All_ovv['DATE'].isin(O2_act.loc[O2_act['hash'] == h1,'DATE'].unique()))]
#    n= 4*ID/ (ID +IR/N)


class ORR_calculations:
    """
    Prepares the ORR data. Subtracts background N2 scan from ORR current:
    [ N2_jcorr, N2_Icorr ]
    """

    E_slices = {"mid": (0.25, 0.5), "high": (0.85, 0.95)}
    N2_corr_E_slice = {"mid": (0.35, 0.65)}
    mA = 1000

    N2_testcols = ["jmAcm-2", "N2_jcorr"]

    globals()["EvRHE"] = "E_AppV_RHE"

    def __init__(self, fit_run_arg):
        assert "ORR_scan_data" in str(type(fit_run_arg))

        self.ORR_data = fit_run_arg
        self.add_attributes()

        self.prepare_N2_BG_data()

        self.read_CV_disk_data()
        self.check_ORR_capacity()

        self.read_ring_data()

    def add_attributes(self):

        for att in self.ORR_data.__dict__.keys():
            setattr(self, att, getattr(self.ORR_data, att))

        self.DestDir = self.ovv_disk.ORR_ecexp_destdir.unique()[0]
        self._set_electrode_info()

    def _doublecheck_O2_act(self):
        assert self.O2_act.PAR_file.nunique() == 1
        assert self.O2_act.PAR_file.unique()[0] == self.PAR_file_disk
        # logger.warning(f'ORR difference in PAR_files from OVV and CreateCV:\n{PAR_file_ovv_file} vs {PAR_file}')

    def read_CV_disk_data(self):
        O2_CVs, O2_action = create_CVs(self.ovv_disk)
        self.O2_disk_raw_data, self.O2_disk_raw_actions = O2_CVs, O2_action

        if not O2_CVs.empty:
            _logger.info(f"Starting ORR for {self.DestDir}/{self.PAR_file_disk.stem}")
            O2_act = self.O2_disk_raw_data.query(
                '(Gas == "O2") & (Type_action == "Cyclic Voltammetry (Multiple Cycles)") & (Scanrate == 0.01) & (SampleID != "Pt_ring")'
            )
            if not O2_act.empty:
                self.O2_act = O2_act
                self.O2_act_seggrp = self.O2_act.groupby("Segment #")
                self._doublecheck_O2_act()
        else:
            _logger.warning(
                f"Not starting ORR -> ERROR empty for {self.DestDir}/{self.PAR_file_disk.stem}"
            )

    def add_rpm_list(self):
        # O2_segs = self.O2_act['Segment #'].unique()
        # O2_act_seggrp = self.O2_act.groupby('Segment #')
        _longsegs = [n for n, gr in self.O2_act_seggrp if len(gr) > 1000]
        _rpm_ref_lists = {
            5: [0, 200, 400, 900, 1500],
            3: [0, 900, 1500],
            2: [1500, 1500],
            4: [0, 200, 400, 900],
            6: [0, 200, 400, 900, 1500, 2500],
        }
        _rpm_list = _rpm_ref_lists.get(len(_longsegs), [])
        self.rpm_list = _rpm_list
        _logger.debug("Jkin calculation rpm list used: {0}".format(_rpm_list))

    def _set_electrode_info(self):
        # TODO move out and input as kwargs
        # electrode_name, collN, SA_disk, SA_ring
        for k, val in self.kwargs["WE_SA_collection_eff"]("PINE").items():
            setattr(self, k, val)

    def read_ring_data(self):
        O2_chrono_Ring, O2_chrono_actions = ORR_read_Ring_file(
            self.ovv_all, self.ovv_disk
        )
        self.O2_ring_raw_data, self.O2_ring_raw_actions = (
            O2_chrono_Ring,
            O2_chrono_actions,
        )

        if not self.O2_ring_raw_data.empty:
            self.PAR_file_ring = O2_chrono_Ring.PAR_file.unique()[0]
        if (
            "DAC Control" not in self.O2_disk_raw_actions.Type_action.unique()
            and "DAC Control" in O2_chrono_actions.Type_action.unique()
            and sum(self.O2_act.RPM_DAC.unique()) == 0
        ):
            O2_chrono_Ring[["Segment #", "RPM_DAC"]]
            _seg_rpm_cols = ["Segment #", "RPM_DAC"]
            _seg_rpms_chrono = set(
                [(i[0], i[1]) for i in O2_chrono_Ring[_seg_rpm_cols].values]
            )
            _seg_rpm_DF = pd.DataFrame(_seg_rpms_chrono, columns=_seg_rpm_cols)
            O2_act = pd.merge(
                self.O2_act.drop(columns="RPM_DAC"),
                _seg_rpm_DF,
                on=_seg_rpm_cols[0],
                how="left",
            )
            self.O2_act = O2_act
            _logger.warning(f"Jkin calculation, RPM_DAC used from Ring self.fstem")
            # TODO Match Segments with RPM_DACs of Ring and Merge with Disk

    def _check_capacity_sweeps(self, dataDF):

        all_res = []
        E_res = {}
        for Ename, Eval in self.E_slices.items():
            Elow, Ehigh = Eval
            _ORRcheck = dataDF.loc[(dataDF[EvRHE] < Ehigh) & (dataDF[EvRHE] > Elow)]
            _name = f"{Ename}_{Elow}_{Ehigh}"

            if not _ORRcheck.empty:
                _segres = {}
                # _segr = []
                for seg, sgrp in _ORRcheck.groupby("Segm"):
                    _res = {}
                    for swp, swpgrp in sgrp.groupby("Sweep_Type"):
                        jmean = swpgrp["jmAcm-2"].mean()
                        _res = {**_res, **{swp: jmean}}
                    if "anodic" and "cathodic" in _res.keys():
                        _diff = _res["anodic"] - _res["cathodic"]
                        _res = {**_res, **{"diff": _diff}}
                    # _segr.append({seg : _res})
                    _segres = {**_segres, **{seg: _res}}
                E_res = {**E_res, **{_name: _segres}}
                # all_res.append({_name : _segr})

        _E_results = [
            {"name": f"{k}", "seg": f"{k1}_{k3}", "val": v2["diff"]}
            for k, val in E_res.items()
            for k1, v2 in val.items()
            for k3, v3 in v2.items()
            if "diff" in k3
        ]
        _check_cap = pd.DataFrame(_E_results)
        return _check_cap

    def check_ORR_capacity(self):
        ORR_check_cap = self._check_capacity_sweeps(
            self.O2_act.loc[(self.O2_act.Segm > 1)]
        )
        self.ORR_check_cap = ORR_check_cap

    def compare_cap(self):

        if hasattr(self, "ORR_check_cap") and hasattr(self, "N2_check_cap"):
            _ORR_cap_agg = self.ORR_check_cap.groupby("name").agg("val").mean()
            _N2_diff_lst = [
                gr.assign(
                    **{
                        "val_diff": gr.val - _ORR_cap_agg.loc[n],
                        "O2_val": _ORR_cap_agg.loc[n],
                    }
                )
                for n, gr in self.N2_check_cap.groupby("name")
            ]
            N2_cap_diff = pd.concat(_N2_diff_lst, ignore_index=True).sort_values(
                ["name", "val_diff"]
            )

    def __iter__(self):

        for N2_cols in self.N2_testcols:
            for N2_pf, N2_BG in self.N2_BG_data_grp:
                N2_cols
                N2_pf, N2_BG
                # O2_act_slice = O2_act.loc[(O2_act['Segment #'] == seg)]
                for seg, O2_disk_seg in self.O2_act_seggrp:
                    seg, O2_disk_seg

                    yield self, O2_disk_seg, N2_BG, N2_cols


#      def _test():
# #         ORR_disk_pars = ORR_disk_pars.assign(**_N2BG_pars).drop_duplicates()
# # #            (PAR_file,rpm_n, seg, ORR_dest_dir_file, dest_file, O2_join)
# # #                                'Analysis_date': datetime.now()}

# #             ### +++++++++-   ++++++++++++++++++++++    -+++++++++++ ###
# #     ### +++++++++- PLOTTING -+++++++++++ ###
# #             # https://stackoverflow.com/questions/37737538/merge-matplotlib-subplots-with-shared-x-axis
# #             #                fig, (axRing,axJ) = plt.subplots(2, sharex=True)
# #             #                SegGr.loc[:,[x]+y].to_csv(ORR_dest_dir.joinpath(Path(f).stem+'_%s'%rpm+'.csv'))
# #             #                                gr.plot(x=x,y=y,xlim=(0,1.2),ylim=(-6,0.5),title=Path(f).stem)
# #             #                                plt.savefig(ORR_dest_dir.joinath(Path(f).stem)+'_O2.png',bbox_inches='tight')
# #             #                                plt.close()
# #         ORR_plot_ring_disk( O2_join, ORR_disk_pars , Iring_Chrono, _ORR_RRDE_pars,
# #                            ORR_dest_dir_file, dest_file)
# # ### +++++++++- FINALIZING -+++++++++++ ###
# #         #                O2_out_ovv,O2_CalcOut = pd.DataFrame([]), pd.DataFrame([])

#         O2_CalcOut = pd.merge(O2_join, Iring_Chrono, on=[EvRHE,'Sweep_Type'],how='left')
# #            O2_join.plot(x=EvRHE, y=['Jcorr']),  Iring_Chrono.plot(x=EvRHE, y=['J_ring'])
# #            O2_CalcOut.groupby('Sweep_Type').plot(x=EvRHE, y=['Jcorr','J_ring'])
#         yOut = ['E_AppV_RHE', 'jmAcm-2', 'jmAcm-2_N2', 'Jcorr', 'Jkin_min', 'Jkin_max', 'J_ring',  'Frac_H2O2', 'n_ORR','Elapsed Time(s)','Sweep_Type','RPM_DAC']
#         #                O2_CalcOut[yOut].to_excel(ORR_dest_dir.joinpath(dest_file+'_RRDE.xlsx'))

#         ORR_Excel_dest_base = ORR_dest_dir_file.joinpath(f'{dest_file}_RRDE.xlsx')
#         ORR_Excel_dest = FileOperations.CompareHashDFexport(O2_CalcOut[yOut], ORR_Excel_dest_base)

#         O2_ParsOut = pd.merge(ORR_disk_pars,_ORR_RRDE_pars).assign(**{'RRDE_merged_data_file' : ORR_Excel_dest}).drop_duplicates()
#         _ORR_parslst.append(O2_ParsOut )
#         ORR_KL_data_rpm = ORR_select_KL_data(O2_CalcOut)
#         _KL_data_all_rpms.append(ORR_KL_data_rpm)


# def add_mean_Jcorr_col(_DF):
#    cath = _DF.groupby('Sweep_Type').get_group('cathodic')
#    anod = _DF.groupby('Sweep_Type').get_group('anodic')
#    _n1cols = [i for i in _DF.columns if _DF[i].nunique() ==1]
#    swp_merge = pd.merge_asof(cath, anod[[i for i in anod.columns if i not in _n1cols]],on=EvRHE,suffixes=['_cath','_anod'])
#    swp_merge = swp_merge.assign(**{'Sweep_Type' : 'mean','Jcorr_raw': swp_merge[['Jcorr_raw_cath', 'Jcorr_raw_anod']].mean(axis=1),
#                                    'Jcorr_raw': swp_merge[['Jcorr_raw_cath', 'Jcorr_raw_anod']].mean(axis=1)})
##                _n1cols_an = [i for i in swp_merge if swp_merge[i].nunique() ==1 if i.endswith('_anod')]
##                _n1cols_cath = [i for i in swp_merge if swp_merge[i].nunique() ==1 if i.endswith('_cath')]
##                [i for i in _n1cols_an if swp_merge[i].unique()[0] == swp_merge[i.replace('_anod','_cath')].unique()[0] ]
##                [swp_merge[i].unique() for i in _n1cols_an]
##                swp_merge.plot(x='E_AppV_RHE',y=['Jcorr_raw_cath', 'Jcorr_raw_anod','Jcorr_raw'],xlim=(0,1.1),ylim=(-3E0,1E0))
#    O2_join_Jmean = pd.concat([_DF, swp_merge],ignore_index=True)
#    return O2_join_Jmean
#
