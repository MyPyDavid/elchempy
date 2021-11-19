#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:07:59 2020

@author: zmg
"""

# import sys
from pathlib import Path

# from collections import namedtuple
# from datetime import datetime
import numpy as np

# from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines

# import os
# import multiprocessing
# from functools import partial
# from itertools import repeat
import pandas as pd

# from file_py_helper.find_folders import FindExpFolder
# from file_py_helper.file_functions import FileOperations

if __name__ == "__main__":
    pass


import logging

logger = logging.getLogger(__name__)


def ORR_plot_ring_disk(
    O2_join, ORR_disk_pars, Iring_Chrono, _ORR_RRDE_pars, ORR_dest_dir_file, dest_file
):
    #        (O2_join, Iring_Chrono, Diff_lim, E_onset,E_half,Jkin_075,Jkin_080, N2_bg_file_used,
    #                       Ring_ovv, RHE_OCP_0,rpm_n, ORR_dest_dir, dest_file)
    #%%
    EvRHE = "E_AppV_RHE"
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3], hspace=0.01)
    #                fig.add_subplot()
    #                axJ = O2_join.plot(x=EvRHE,y='Jcorr',xlim=(0,1.2),ylim=(-6,0.3),label=O2_join['SampleID'].unique()[0]+'_'+str(rpm_list[rpm])+' rpm')
    #                fig.set_title(O2_join['SampleID'].unique()[0]+'_'+str(rpm_list[rpm])+' rpm')
    axJ = plt.subplot(gs[1])

    #    axJ.plot(O2_join.iloc[:-1][EvRHE], O2_join.iloc[:-1]['jmAcm-2_N2'], 'k--', label='$N_{2}\/scan$')
    O2_join_swgrp = O2_join.groupby("Sweep_Type")

    _cath_O2join = O2_join_swgrp.get_group("cathodic")
    _cath_O2meta = (
        _cath_O2join[
            [i for i in _cath_O2join.columns if _cath_O2join[i].nunique() == 1]
        ]
        .iloc[0]
        .to_dict()
    )
    _anod_O2join = O2_join_swgrp.get_group("anodic")
    if "mean" in O2_join_swgrp.groups.keys():
        _mean_O2join = O2_join_swgrp.get_group("mean")
    else:
        _mean_O2join = pd.DataFrame()

    if not Iring_Chrono.empty and not _ORR_RRDE_pars.empty:
        Iring_Chrono_grp = Iring_Chrono.groupby("Sweep_Type")
        _cath_Iring = Iring_Chrono_grp.get_group("cathodic")
        _anod_Iring = Iring_Chrono_grp.get_group("anodic")

        if "mean" in Iring_Chrono_grp.groups.keys():
            _mean_Iring = Iring_Chrono_grp.get_group("mean")
        else:
            _mean_Iring = pd.DataFrame()

        axRing = plt.subplot(gs[0])
        #                , sharex = axJ)
        try:
            axRing.set_title(
                O2_join["SampleID"].unique()[0]
                + "_"
                + str(_cath_O2meta.get("RPM_DAC", -99))
                + " rpm"
            )
            #                    axRing = axJ.twinx()
            #                    axRing.set_ylabel('$j_{R}\//\/mA/cm^{2}$')
            axRing.set_ylabel("$H_{2}O_{2}\/yield\//\/(\%)$")
            try:
                ymax = (
                    Iring_Chrono.query(" (E_AppV_RHE < 0.60)")
                    .iloc[:-1]["Frac_H2O2"]
                    .max()
                    * 1.5
                )
                int(ymax)
                if ymax > 100:
                    ymax = 100
            except Exception as e:
                logger.error(
                    "Jkin calculation ERROR Iring_Chrono ymax: Iring_Chrono {0}".format(
                        e
                    )
                )
                ymax = 100

            axRing.set_ylim(0, ymax)
            axRing.set_xlim([0, 1])
            E_AppV_ring = Iring_Chrono[EvRHE + "_ring"].dropna().unique().mean()
            _OCP_ring = _ORR_RRDE_pars.Ring_Measured_OCP.unique().mean()

            #            Iring_Chrono.query('Sweep_Type == "cathodic"')[EvRHE + '_ring'].dropna().unique()[0]
            axRing.annotate(
                f"E_ring: {E_AppV_ring:.2f} V_rhe",
                xy=(0.1, 0.9),
                xycoords="axes fraction",
                verticalalignment="bottom",
                fontsize=8,
            )
            axRing.annotate(
                f"OCP_ring vs RE: {_OCP_ring:.2f} mV",
                xy=(0.1, 0.8),
                xycoords="axes fraction",
                verticalalignment="bottom",
                fontsize=8,
            )

            axRing.plot(
                _cath_Iring.iloc[:-1]["E_AppV_RHE"],
                _cath_Iring.iloc[:-1]["Frac_H2O2"],
                "b",
            )
            axRing.plot(
                _anod_Iring.iloc[:-1]["E_AppV_RHE"],
                _anod_Iring.iloc[:-1]["Frac_H2O2"],
                "r--",
            )
            if not _mean_Iring.empty:
                axRing.plot(
                    _mean_Iring.iloc[:-1]["E_AppV_RHE"],
                    _mean_Iring.iloc[:-1]["Frac_H2O2"],
                    "g-.",
                )
        except Exception as e:
            logger.error("Jkin calculation ERROR in  Chrono Ring plotting: %s" % e)

    axJ.plot(
        _cath_O2join.iloc[:-1][EvRHE],
        _cath_O2join.iloc[:-1]["Jcorr"],
        "b",
        label="cathodic",
    )
    axJ.plot(
        _cath_O2join.iloc[:-1][EvRHE],
        _cath_O2join.iloc[:-1]["Jcorr_raw"],
        "b",
        label="Jcorr raw cathodic",
        alpha=0.4,
    )
    axJ.plot(
        _cath_O2join.iloc[:-1][EvRHE],
        _cath_O2join.iloc[:-1]["jmAcm-2"],
        "k",
        label="J raw cathodic",
        alpha=0.2,
    )

    axJ.plot(
        _anod_O2join.iloc[:-1][EvRHE],
        _anod_O2join.iloc[:-1]["Jcorr"],
        "r--",
        label="anodic",
    )
    axJ.plot(
        _anod_O2join.iloc[:-1][EvRHE],
        _anod_O2join.iloc[:-1]["Jcorr_raw"],
        "r--",
        label="Jcorr raw anodic",
        alpha=0.4,
    )
    axJ.plot(
        _anod_O2join.iloc[:-1][EvRHE],
        _anod_O2join.iloc[:-1]["jmAcm-2"],
        "k--",
        label="J raw anodic",
        alpha=0.2,
    )

    if not _mean_O2join.empty:
        axJ.plot(
            _mean_O2join.iloc[:-1][EvRHE],
            _mean_O2join.iloc[:-1]["Jcorr"],
            "g-.",
            label="mean",
        )
        axJ.plot(
            _mean_O2join.iloc[:-1][EvRHE],
            _mean_O2join.iloc[:-1]["Jcorr_raw"],
            "g-.",
            label="Jcorr raw mean",
            alpha=0.4,
        )
        axJ.plot(
            _mean_O2join.iloc[:-1][EvRHE],
            _mean_O2join.iloc[:-1]["jmAcm-2"],
            "k-.",
            label="J raw mean",
            alpha=0.2,
        )

    try:
        #        _cath_O2join, _anod_O2join
        #        axJ.plot(O2_join.iloc[:-1][EvRHE], O2_join.iloc[:-1]['jmAcm-2_N2'], 'k--', label='$N_{2}\/scan$')
        _N2_use_col = ORR_disk_pars.N2_BG_usecols_for_BG_correction.unique()[0]
        axJ.plot(
            _cath_O2join.iloc[:-1][EvRHE],
            _cath_O2join.iloc[:-1][_N2_use_col],
            "k",
            label=f"{_N2_use_col}",
            alpha=0.6,
        )
        axJ.plot(
            _anod_O2join.iloc[:-1][EvRHE],
            _anod_O2join.iloc[:-1][_N2_use_col],
            "k--",
            label=f"{_N2_use_col}",
            alpha=0.6,
        )
    except Exception as e:
        logger.error(
            "Jkin calculation ERROR in  Chrono Ring plotting on N2 SCAN: %s" % e
        )

    #                axJ.set_title(O2_join['SampleID'].unique()[0]+'_'+str(rpm_list[rpm])+' rpm')
    #                label=O2_join['SampleID'].unique()[0]+'_'+str(rpm_list[rpm])+' rpm'
    axJ.set_xlim([0, 1])
    axJ.set_ylim([-6.2, 1])

    if not ORR_disk_pars.empty:
        _cath_ORR_pars = ORR_disk_pars.query('Sweep_Type == "mean"')
        cath_pars = (
            _cath_ORR_pars[
                [i for i in _cath_ORR_pars.columns if _cath_ORR_pars[i].nunique() == 1]
            ]
            .iloc[0]
            .to_dict()
        )
        try:
            _diff_lim = cath_pars.get("ORR_J_diff_lim", -6)
            if _diff_lim < -6:
                axJ.set_ylim([_diff_lim, 1])
            #                xlim=(0,1.2),ylim=(-6,0.3),
            l = mlines.Line2D(
                [O2_join[EvRHE].values],
                [len(O2_join[EvRHE]) * [0]],
                ls="dashdot",
                lw=0.5,
                color="grey",
            )
            axJ.add_line(l)
            _fs = 14
            _diff_lim_E = _cath_O2join.loc[
                np.isclose(_cath_O2join["Jcorr"], _diff_lim, atol=0.01), EvRHE
            ].mean()
            axJ.annotate(
                "$J_{diff} = %.2f\/mA/cm^{2}$" % (_diff_lim),
                (_diff_lim_E, _diff_lim),
                xytext=(0, 40),
                textcoords="offset points",
                fontsize=_fs,
                arrowprops=dict(arrowstyle="-|>"),
            )

            _E_half = cath_pars.get("ORR_E_half", 0.5)
            _E_half_J = _cath_O2join.loc[
                np.isclose(_cath_O2join[EvRHE], _E_half, atol=0.005), "Jcorr"
            ].mean()
            axJ.annotate(
                "$E_{1/2}$ = %.2f $V_{RHE}$" % (_E_half),
                (_E_half, _E_half_J),
                xytext=(-120, 20),
                textcoords="offset points",
                fontsize=_fs,
                arrowprops=dict(arrowstyle="-|>"),
            )

            _Jkin0750 = cath_pars.get("ORR_Jkin_min_750", 0.500)
            _Jkin080 = cath_pars.get("ORR_Jkin_min_800", 0.499)
            _Jkin0750_Jcorr = _cath_O2join.loc[
                np.isclose(_cath_O2join["Jkin_min"], _Jkin0750, atol=0.005), "Jcorr"
            ].mean()

            axJ.annotate(
                "$J_{kin(0.75 V_{RHE})}$ \n = %.2f $\/mA/cm^{2}$ \n $@0.8V_{RHE}$=%.2f"
                % (np.abs(_Jkin0750), np.abs(_Jkin080)),
                (0.75, _Jkin0750_Jcorr),
                xytext=(20, -100),
                textcoords="offset points",
                fontsize=_fs,
                arrowprops=dict(arrowstyle="-|>"),
            )

            _E_onset = cath_pars.get("ORR_E_onset", 0.5)
            _E_onset_J = _cath_O2join.loc[
                np.isclose(_cath_O2join[EvRHE], _E_onset, atol=0.005), "Jcorr"
            ].mean()
            axJ.annotate(
                "$E_{onset}$: %.2f $V_{RHE}$" % (_E_onset),
                (_E_onset, _E_onset_J),
                xytext=(30, -20),
                textcoords="offset points",
                fontsize=_fs,
                arrowprops=dict(arrowstyle="-|>"),
            )
            axJ.text(
                0.05,
                -5.2,
                f'OCP Disk vs RE: {cath_pars.get("Measured_OCP", 0):.0f} mV',
                fontsize=10,
            )
            axJ.text(
                0.05,
                -5.5,
                f'OCP RE vs RHE: {cath_pars.get("RHE_OCP_mV", 0):.0f} mV',
                fontsize=10,
            )
            #                     '$\mathrm{RE vs RHE:}$ %.0f mV' % (cath_pars.get('RHE_OCP_mV', 0)), fontsize=10)

            #    cath_pars['PAR_file_N2']

            # _N2_use_col= cath_pars.get('N2_BG_usecols_for_BG_correction')
            axJ.annotate(
                f'N2 bg: {Path(cath_pars["PAR_file_N2"]).name} with {_N2_use_col} at {cath_pars.get("Scan Rate (V/s)_N2", 0):.2f} V/s ',
                xy=(0.025, 0.022),
                xycoords="figure fraction",
                horizontalalignment="left",
                verticalalignment="bottom",
                fontsize=12,
            )
            axJ.annotate(
                f"O2 ring: {_ORR_RRDE_pars.Ring_PAR_file.unique()[0].name}",
                xy=(0.025, 0),
                xycoords="figure fraction",
                horizontalalignment="left",
                verticalalignment="bottom",
                fontsize=12,
            )
        except Exception as e:
            logger.error(f"Jkin calc PLOTTING ERROR in ORR Pars:{e}")

    axJ.set_ylabel("$j_{Disk}\//\/mA/cm^{2}$")
    axJ.set_xlabel("$E_{Disk} \//\/V_{RHE}$")
    axJ.legend(loc=4, fontsize=8)
    #                axJ.text(-0.1,-7.5,'$\mathrm{RE vs RHE:}$ %.0f mV' %(RHE_OCP_0),fontsize=10)
    #                ORR_dest_dir.mkdir(parents=True,exist_ok=True)
    #                O2_join.to_csv(ORR_dest_dir.joinpath(dest_file+'.csv'))
    #                plt.show()
    #                plt.tight_layout()
    #%%
    plt.savefig(ORR_dest_dir_file.joinpath(dest_file + "_RRDE.png"), dpi=100)
    #                plt.savefig(ORR_dest_dir.joinpath(dest_file+'_RRDE.png'),dpi=100,bbox_inches='tight')
    plt.close()
