# import sys
from pathlib import Path

# from collections import namedtuple
# from datetime import datetime
import decimal

# import os
# import multiprocessing
# from functools import partial
# from itertools import repeat

import numpy as np
import pandas as pd


import scipy

# from scipy.stats import linregress
import matplotlib.pyplot as plt

# from file_py_helper.find_folders import FindExpFolder
# from file_py_helper.file_functions import FileOperations

# print('File', __file__,'\nName;',__name__)
if __name__ == "__main__":
    pass

import logging

logger = logging.getLogger(__name__)
# __all__ = ['O2_Chrono']


def round_ElapsedTime_decimal_DR(Disk, ChronoI):
    tMin, tMax = Disk["Elapsed Time(s)"].min(), Disk["Elapsed Time(s)"].max()
    _tMin_dec = decimal.Decimal(str(tMin))
    _tDisk_len = len(str(_tMin_dec).split(".")[1])

    _R_tMin, _R_tMax = (
        ChronoI["Elapsed Time(s)"].min(),
        ChronoI["Elapsed Time(s)"].max(),
    )
    _R_tMin_dec = decimal.Decimal(str(_R_tMin))
    _tRing_len = len(str(_R_tMin_dec).split(".")[1])
    if _tDisk_len < _tRing_len:
        if np.round(_R_tMin, _tDisk_len) == tMin:
            _ET_round = ChronoI["Elapsed Time(s)"].round(_tDisk_len)
            if _ET_round.nunique() == len(_ET_round):
                ChronoI = ChronoI.assign(
                    **{
                        "Elapsed_Time_raw_ring": ChronoI["Elapsed Time(s)"],
                        "Elapsed Time(s)": _ET_round,
                    }
                )
    elif _tDisk_len > _tRing_len:
        if np.round(tMin, _tRing_len) == tMin:
            _ET_round = Disk["Elapsed Time(s)"].round(_tRing_len)
            if _ET_round.nunique() == len(_ET_round):
                Disk = Disk.assign(
                    **{
                        "Elapsed_Time_raw_disk": Disk["Elapsed Time(s)"],
                        "Elapsed Time(s)": _ET_round,
                    }
                )
    else:
        pass
    return Disk.sort_values("Elapsed Time(s)"), ChronoI.sort_values("Elapsed Time(s)")


# sweep_col= 'Sweep_Type_disk'
# E_col = 'E_AppV_RHE_disk'
#  mean_cols = ['I(A)_disk','I(A)_ring']
def add_mean_Jcorr_col(
    O2DF,
    sweep_col="Sweep_Type",
    E_col="E_AppV_RHE",
    mean_cols=["I(A)_disk", "I(A)_ring", "Jcorr"],
):
    #    O2DF = ORR_Both
    O2DF_res = O2DF.reset_index()
    if "Sweep_Type" not in O2DF_res.columns and "Sweep_Type_disk" in O2DF_res.columns:
        O2DF_res = O2DF_res.assign(**{"Sweep_Type": O2DF_res.Sweep_Type_disk})

    cath = O2DF_res.groupby(sweep_col).get_group("cathodic").sort_values(E_col)
    anod = O2DF_res.groupby(sweep_col).get_group("anodic").sort_values(E_col)
    _n1cols = [i for i in O2DF_res.columns if O2DF_res[i].nunique() <= 2]
    _nXcols = [i for i in O2DF_res.columns if O2DF_res[i].nunique() > 2]
    mean_cols_used = [i for i in mean_cols if i in O2DF.columns and i not in _n1cols]
    swp_merge = pd.merge_asof(
        cath,
        anod[[i for i in anod.columns if i not in _n1cols]],
        on=E_col,
        suffixes=["_cath", "_anod"],
    )
    _mean_dct = {
        i: swp_merge[[f"{i}_cath", f"{i}_anod"]].mean(axis=1) for i in mean_cols_used
    }
    swp_merge = swp_merge.assign(
        **{**{"Sweep_Type": "mean", "Sweep_Type_disk": "mean"}, **_mean_dct}
    )
    #    swp_merge = swp_merge.assign(**{'Jcorr_O2diff' : (swp_merge['Jcorr_raw'].rolling(31).mean().diff() *1E3),
    #                'Jcorr_O2diffdiff' : (swp_merge['Jcorr_raw'].rolling(31).mean().diff().rolling(31).mean().diff().rolling(31).mean())*1E5})
    #                _n1cols_an = [i for i in swp_merge if swp_merge[i].nunique() ==1 if i.endswith('_anod')]
    #                _n1cols_cath = [i for i in swp_merge if swp_merge[i].nunique() ==1 if i.endswith('_cath')]
    #                [i for i in _n1cols_an if swp_merge[i].unique()[0] == swp_merge[i.replace('_anod','_cath')].unique()[0] ]
    #                [swp_merge[i].unique() for i in _n1cols_an]
    #                swp_merge.plot(x='E_AppV_RHE',y=['Jcorr_raw_cath', 'Jcorr_raw_anod','Jcorr_raw'],xlim=(0,1.1),ylim=(-3E0,1E0))
    swp_merge = swp_merge.drop(
        columns=[
            i for i in swp_merge.columns if i.endswith("anod") or i.endswith("cath")
        ]
    )
    Te_mean = pd.concat([O2DF_res, swp_merge], ignore_index=True, axis=0).dropna(
        subset=[sweep_col], axis=0
    )
    #    Te_mean =  O2DF_res.append(swp_merge)
    #    Te_mean.plot(x=E_col,y=mean_cols)
    #    Te_mean.loc[Te_mean['Sweep_Type_disk'] == 'mean'].plot(x=E_col,y=mean_cols)
    #    Te_mean.loc[Te_mean['Sweep_Type_disk'] == 'mean']
    #    tt2 = Te_mean.loc[Te_mean['Sweep_Type_disk'] == 'mean']
    return Te_mean


def _debugging(self):
    seg, rpm_n, O2_join, O2_chrono_Ring, ORR_dest_dir_file = (
        self.segment,
        self.RPM_DAC,
        self.O2_disk_N2BG_mean_Jkin,
        self.ORR_calc.O2_ring_raw_data,
        self.ORR_dest_dir_file,
    )
    # self.ORR_calc.O2_ring_raw_data
    plot_BipotChrono, export_excel, ring_filter = True, True, True


def O2_Chrono(
    seg,
    rpm_n,
    O2_join,
    O2_chrono_Ring,
    ORR_dest_dir_file,
    plot_BipotChrono=True,
    export_excel=True,
    ring_filter=True,
    electrode_properties={},
):
    """Calculates the ring-disk properties for the RRDE measurement"""
    EvRHE = "E_AppV_RHE"
    # TODO fix WE_SA_collection_eff('PINE')

    electrode_name, CollEff, SA_disk, SA_ring = electrode_properties.values()
    mA, Fltr_Win = 1000, 71
    #    index_info_ORR_O2Chrono = {}
    RingDiskDir = ORR_dest_dir_file.joinpath("RingDisk")
    RingDiskDir.mkdir(parents=True, exist_ok=True)

    O2_join_seg = O2_join.loc[(O2_join["Segment #"] == seg)]

    _meta_O2_join = (
        O2_join_seg[[i for i in O2_join_seg.columns if O2_join_seg[i].nunique() == 1]]
        .iloc[0]
        .to_dict()
    )
    _meta_O2_join.update({"ring_filter": ring_filter})
    #    ====
    RRDE_cols = [
        "Elapsed Time(s)",
        "E_AppV_RHE",
        "I(A)",
        "Segment #",
        "Sweep_Type",
        "Scanrate",
        "I(A)_raw",
        "PAR_file",
    ]
    if "Elapsed Time(s)" not in O2_join:
        O2_join.reset_index(inplace=True)
    #    Disk = O2_act.loc[(O2_act['hash'] == h1 ) & (O2_act['Segment #'] == seg),RRDE_cols]
    #    ====
    Disk = O2_join_seg
    #    DiskFile = O2_join_seg.loc[(O2_join['Segment #'] == seg), 'PAR_file'].unique()[0]

    #    tMin, tMax = Disk['Elapsed Time(s)'].min(), Disk['Elapsed Time(s)'].max()
    #    _tMin_dec = decimal.Decimal(str(tMin))
    #    _tDisk_len = len(str(_tMin_dec).split(".")[1])
    #
    #    ORRdisk, ORRring = OER_CV.set_index('Elapsed Time(s)'), OER_Chrono.set_index('Elapsed Time(s)')
    #   OER_Both = OER1.join(OER2,lsuffix='_disk',rsuffix='_ring',how='outer')
    #   OER_Both.loc[:,['Segment #_disk',EvRHE+'_disk']] = OER_Both.loc[:,['Segment #_disk',EvRHE+'_disk']].interpolate()
    #    ChronoI = All_ovv.loc[(All_ovv['Elapsed Time(s)'].isin(Disk['Elapsed Time(s)'])) &  (All_ovv['Type'] =='Chronoamperometry') & (All_ovv['DATE'].isin(O2_act.loc[O2_act['hash'] == h1,'DATE'].unique())),['I(A)','Elapsed Time(s)']]
    ChronoI = O2_chrono_Ring.loc[
        (O2_chrono_Ring["Type_action"] == "Chronoamperometry")
        & (O2_chrono_Ring["Segment #"] == seg)
        & (O2_chrono_Ring["PAR_date_min"].isin(O2_join.PAR_date_min.unique())),
        :,
    ]
    ChronoI = ChronoI.assign(**{"I(A)_raw": ChronoI["I(A)"].values})

    Disk, ChronoI = round_ElapsedTime_decimal_DR(Disk, ChronoI)
    _chrono_meta = {
        f"Ring_{k}": val
        for k, val in ChronoI[[i for i in ChronoI.columns if ChronoI[i].nunique() == 1]]
        .iloc[0]
        .to_dict()
        .items()
    }

    RingFile = ChronoI["PAR_file"].unique()[0]
    ChronoI = ChronoI[
        RRDE_cols + [i for i in ChronoI.columns if "_raw" in i and i not in RRDE_cols]
    ]
    #    ChronoI['I(A)'] = ChronoI['I(A)_raw'].rolling(30).mean()
    if ChronoI.empty:
        logger.warning("ORR O2 Chrono Ring chrono empty, {0}".format(ORR_dest_dir_file))
    ### Apply filter to Ring Current to reduce noise ###
    if ring_filter:
        try:
            ChronoI["I(A)"] = scipy.signal.savgol_filter(
                ChronoI["I(A)_raw"].values, Fltr_Win, 3
            )
        except:
            ChronoI["I(A)"] = ChronoI["I(A)_raw"].rolling(30).mean()
    #    ChronoI.plot(x='Elapsed Time(s)',y=['I(A)_raw','I(A)'])
    #    ChronoI = All_ovv.loc[(All_ovv['Elapsed Time(s)'].isin(Disk['Elapsed Time(s)'])) &  (All_ovv['Type'] =='Chronoamperometry') & (All_ovv['DATE'].isin(O2_act.loc[O2_act['hash'] == h1,'DATE'].unique())),['I(A)','Elapsed Time(s)']]

    ORR_Both = pd.merge_asof(
        Disk.dropna(subset=["Elapsed Time(s)"]),
        ChronoI,
        on="Elapsed Time(s)",
        suffixes=["_disk", "_ring"],
    )
    ORR_Both = ORR_Both.dropna(axis=1, how="all")
    #   'OLD way with join'
    #    Disk.set_index('Elapsed Time(s)', inplace=True), ChronoI.set_index('Elapsed Time(s)', inplace=True)
    #    ORR_Both = Disk.join(ChronoI, lsuffix='_disk', rsuffix='_ring', how='outer')
    _ipcols = [
        "Segment #_disk",
        EvRHE + "_disk",
        "Scanrate_disk",
        "I(A)_disk",
        "I(A)_ring",
    ]
    ORR_Both.loc[:, _ipcols] = ORR_Both.loc[:, _ipcols].interpolate()
    ORR_Both = ORR_Both.assign(**{EvRHE: ORR_Both[f"{EvRHE}_disk"]})
    mean_J_cols = [
        i
        for i in ORR_Both.columns
        if (i.upper().startswith("J") or "I" in i) and "float" in str(ORR_Both[i].dtype)
    ]
    mean_J_cols += ["index", f"{EvRHE}_ring", f"{EvRHE}_disk"]
    ORR_Both = add_mean_Jcorr_col(
        ORR_Both, sweep_col="Sweep_Type_disk", E_col=EvRHE, mean_cols=mean_J_cols
    )

    Te = ORR_Both[
        [f"{EvRHE}_disk", "I(A)_ring", "I(A)_disk", "Sweep_Type_disk", f"{EvRHE}_ring"]
    ].rename(columns={f"{EvRHE}_disk": EvRHE, "Sweep_Type_disk": "Sweep_Type"})
    if "mean" not in Te.Sweep_Type.unique():
        Te = add_mean_Jcorr_col(
            Te, sweep_col="Sweep_Type", E_col=EvRHE, mean_cols=mean_J_cols
        )

    #    Te['I(A)_ring'].rolling_mean(10)
    #    Te = pd.merge_asof(Disk.sort_values(by='Elapsed Time(s)'),ChronoI.sort_values(by='Elapsed Time(s)').iloc[::2],on='Elapsed Time(s)',suffixes=('_disk', '_ring'))

    ###### === Calculations of Selectivity H2O2 yield, n electron from I(Ring and I(Disk) ==== #######
    def RRDE_calc_H2O2(_DF, CollEff, SA_ring, mA):
        """Calculations of Selectivity H2O2 yield, n electron from I(Ring and I(Disk)"""
        _I_disk, _I_ring = np.abs(_DF["I(A)_disk"]).values, np.abs(
            _DF["I(A)_ring"].values
        )
        _FracH2O2 = 200 * (_I_ring / CollEff) / (_I_disk + _I_ring / CollEff)
        _J_ring = mA * _I_ring / (CollEff * SA_ring)
        n_ORR = 4 * _I_disk / (_I_disk + _I_ring / CollEff)
        _RRDE = {"Frac_H2O2": _FracH2O2, "J_ring": _J_ring, "n_ORR": n_ORR}
        _DF = _DF.assign(**_RRDE)
        return _DF

    Te = RRDE_calc_H2O2(Te, CollEff, SA_ring, mA)
    ORR_Both = RRDE_calc_H2O2(ORR_Both, CollEff, SA_ring, mA)
    #    '''Plotting for testing quick if data is in the frames ...'''
    #    fig,ax =plt.subplots()
    #    Te.groupby('Sweep_Type').plot(x=EvRHE,y='Frac_H2O2',ax=ax,ylim=(0,50))
    #    [gr.plot(x=EvRHE,y='Jcorr',ax=ax,ylim=(-6,1), label=n) for n,gr in ORR_Both.groupby('Sweep_Type')]
    #    fig,ax =plt.subplots()
    #    [gr.plot(x=EvRHE,y='Frac_H2O2',ax=ax,ylim=(0,20), label=n) for n,gr in ORR_Both.groupby('Sweep_Type')]
    #    fig,ax =plt.subplots()
    #    [gr.plot(x=EvRHE,y='Jkin_min',ax=ax,ylim=(0,10), xlim=(0.9,0.5), label=n) for n,gr in ORR_Both.groupby('Sweep_Type')]
    #    Te = Te.assign(**{'Frac_H2O2': 200 * np.abs(Te['I(A)_ring'].values / CollEff) / (
    #                np.abs(Te['I(A)_disk']).values + np.abs(Te['I(A)_ring'].values / CollEff)),
    #                      'J_ring': mA * Te['I(A)_ring'].values / (CollEff * SA_ring),
    #                      'n_ORR': (4 * np.abs(Te['I(A)_disk'].values)) / (
    #                                  np.abs(Te['I(A)_disk']).values + np.abs(Te['I(A)_ring']).values / CollEff)})
    #    ORR_Both = ORR_Both.assign(**{'Frac_H2O2': 200 * np.abs(ORR_Both['I(A)_ring'].values / CollEff) / (
    #                np.abs(ORR_Both['I(A)_disk']).values + np.abs(ORR_Both['I(A)_ring'].values / CollEff)),
    #                                  'J_ring': mA * ORR_Both['I(A)_ring'].values / (CollEff * SA_ring),
    #                                  'n_ORR': (4 * np.abs(ORR_Both['I(A)_disk'].values)) / (
    #                                              np.abs(ORR_Both['I(A)_disk']).values + np.abs(
    #                                          ORR_Both['I(A)_ring']).values / CollEff)})
    try:
        Segment_Chrono = int(
            np.unique(
                [i for i in ORR_Both["Segment #_disk"].unique() if "nan" not in str(i)]
            )[0]
        )
    except:
        Segment_Chrono = str(
            [i for i in ORR_Both["Segment #_disk"].unique() if "nan" not in str(i)][0]
        )

    #    if export_excel:
    #        ORR_Both = drop_cols('ORR', ORR_Both)
    #        _meta_O2_join = ORR_Both[[i for i in ORR_Both.columns if ORR_Both[i].nunique() == 1]].iloc[0].to_dict()
    #        _meta_O2_join, _chrono_meta
    DiskStem = Path(_meta_O2_join["PAR_file"]).stem
    _pars_lst = []
    for swp, gr in ORR_Both.groupby(by="Sweep_Type_disk"):
        RingDisk_out_path = RingDiskDir.joinpath(
            f"{DiskStem}_{Segment_Chrono}_{swp}.xlsx"
        )
        gr.to_excel(RingDisk_out_path)
        _Chrono_pars = {"Sweep_Type": swp, "RRDE_swp_data_file": RingDisk_out_path}
        _Chrono_pars.update(_meta_O2_join)
        _Chrono_pars.update(_chrono_meta)

        for E_mV in range(500, 850, 50):
            _EV = E_mV * 1e-3
            for _parnm in ["J_ring", "Frac_H2O2", "n_ORR"]:
                _parv = gr.loc[np.isclose(gr[EvRHE], _EV, atol=0.005), _parnm].mean()
                _Chrono_pars.update({f"ORR_{_parnm}_{E_mV}": _parv})
        _pars_lst.append(_Chrono_pars)

    _ORR_RRDE_pars = pd.DataFrame(_pars_lst)
    #            index_info_ORR_O2Chrono = {'PAR_file': DiskFile, 'DestFile': RingDisk_out_path,
    #                                       'Type_output': 'ORR_Jkin_calc_RRDE_Chrono', 'Type_exp': 'ORR',
    #                                       'PAR_file_Ring': RingFile, 'Sweep_Type': swp}
    #    ORR_Both.plot(x=EvRHE+'_disk',y=['I(A)_disk','I(A)_ring'])
    if plot_BipotChrono:
        fig, axORR = plt.subplots(1, figsize=(8, 8))
        #        axORR.set_title(nmBase+'(%s,%s)' %(nmBoth[0],nmBoth[1]))
        axRing = axORR.twinx()
        for swp, swpgr in ORR_Both.groupby(by="Sweep_Type_disk"):
            swpgr[[f"{EvRHE}", "I(A)_disk"]].dropna().iloc[3:-3].plot(
                x=EvRHE, y=["I(A)_disk"], ax=axORR, label=[f"{swp} Disk"]
            )
            try:
                swpgr[[f"{EvRHE}", "I(A)_ring"]].dropna().plot(
                    x=EvRHE,
                    y=["I(A)_ring"],
                    label=[f"{swp} Ring"],
                    ax=axRing,
                    **{"ls": "--"},
                )
            except Exception as e:
                logger.warning("ORR O2 Chrono ERROR in O2_Chrono with plotting", e)
        #        ORR_Both[[EvRHE+'_disk','I(A)_disk']].dropna().plot(x=EvRHE+'_disk',y=['I(A)_disk'],ax=axORR,label=['Disk'])
        #        axORR.legend(loc=0)
        #        axRing.legend(loc=1)
        axORR.legend(
            loc="upper right",
            bbox_to_anchor=(0.95, 1.1),
            ncol=3,
            fancybox=True,
            shadow=True,
        )
        axRing.legend(
            loc="upper left",
            bbox_to_anchor=(0.05, 1.20),
            ncol=3,
            fancybox=True,
            shadow=True,
        )
        axORR.set_ylabel("I (disk) / A")
        axRing.set_ylabel("I (ring) / A")

        #        plt.show()

        RingDisk_plot_path = RingDiskDir.joinpath(f"{DiskStem}_{Segment_Chrono}.png")
        #        index_info_ORR_O2Chrono.update({'DestPlot': RingDisk_plot_path})
        #        int(Segment_Chrono)
        plt.savefig(RingDisk_plot_path, dpi=100, bbox_inches="tight")
        plt.close()
    #    Te = pd.merge(Disk,ChronoI.iloc[::1],on='Elapsed Time(s)',how='outer',suffixes=('_disk', '_ring'))
    #    .dropna(axis=0, how='any')
    #    ORR_Both['I(A)_disk']
    return Te, _ORR_RRDE_pars.dropna(axis=1, how="all")


#    'Sweep_Type' : np.where(Te[EvRHE].diff() > 0, 'anodic', np.where(Te[EvRHE].diff() < 0, 'cathodic', 'NA'))
#    Te.query('Sweep_Type == "cathodic"').plot(x=EvRHE,y='J_ring',xlim=(0,1.2),kind='scatter')
##    ,ylim=(0,10E-03))
#    Te.query('Sweep_Type == "anodic"').plot(x=EvRHE,y='J_ring',xlim=(0,1.2),kind='scatter')
##    ,ylim=(0,1E-02))
#    Te.query('Sweep_Type == "cathodic"').plot(x=EvRHE,y='Frac_H2O2',xlim=(0,0.86),kind='scatter',ylim=(0,20))
# out_Jring = gr.loc[np.isclose(gr[EvRHE], _E_mV, atol=0.005), 'J_ring'].mean()
# out_FracH2O2 = gr.loc[np.isclose(gr[EvRHE], _E_mV, atol=0.005), 'Frac_H2O2'].mean()
# out_nORR = gr.loc[np.isclose(gr[EvRHE], _E_mV, atol=0.005), 'Frac_H2O2'].mean()
