"""
plotting functions for N2 experiments.
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import logging

logger = logging.getLogger(__name__)


def N2_plot_raw_scans_scanrate(raw_data, EvRHE="E_vs_RHE", savepath=None):
    fig, ax = plt.subplots(figsize=(8, 8))

    for sr, sgr in raw_data.groupby("scanrate"):
        _seg = int(sgr["Segment #"].unique()[0])
        ax.plot(sgr[EvRHE], sgr["j_mA_cm2"], label=f"{sr} mV, seg {_seg}")
    #                ax.legend(True)
    # fig.suptitle(f"{Path(Cdl_scans.PAR_file.unique()[0]).name}")
    ax.legend()
    ax.set_ylabel("$j \//\/mA/cm^{2}$")
    ax.set_xlabel("$E \//\/V_{RHE}$")
    if savepath:
        plt.savefig(savepath, dpi=100, bbox_inches="tight")
        plt.close()
        return


def N2_plot_Cdl_sweeptype_scatter(Cdl_pars, **kwargs):
    # SampleID, ScanRates, Cdl_fit, Cdl_cath_slice, Cdl_an_slice, N2_dest_dir, N2_fn

    EvRHE = "E_AppV_RHE"
    fig, ax = plt.subplots()

    Cdl_cath_slice = Cdl_pars.query('SweepType == "cathodic"')
    Cdl_an_slice = Cdl_pars.query('SweepType == "anodic"')

    plt.title(
        "%s made with linear fit of\n %s (R=%.3f)"
        % (
            kwargs.get("SampleID"),
            kwargs.get("scanrates"),
            Cdl_pars["lin_rvalue"].mean(),
        )
    )

    ylim = (0, 1.2 * Cdl_pars.Cdl.max())

    Cdl_cath_slice.plot(
        x=EvRHE,
        y="Cdl_corr",
        kind="scatter",
        ylim=ylim,
        color="orange",
        ax=ax,
        label="Cdl_Cath_corr",
    )
    Cdl_cath_slice.plot(
        x=EvRHE,
        y="Cdl",
        kind="scatter",
        ylim=ylim,
        color="r",
        ax=ax,
        label="Cdl_Cath",
    )
    Cdl_an_slice.plot(
        x=EvRHE,
        y="Cdl_corr",
        kind="scatter",
        ylim=ylim,
        color="c",
        ax=ax,
        label="Cdl_Anod_corr",
    )
    Cdl_an_slice.plot(
        x=EvRHE,
        y="Cdl",
        kind="scatter",
        ylim=ylim,
        ax=ax,
        label="Cdl_Anod",
    )

    if kwargs.get("savepath"):
        plt.savefig(
            N2_dest_dir.joinpath(f"Cdl_{N2_fn}.png"), dpi=100, bbox_inches="tight"
        )
        pd.concat([Cdl_an_slice, Cdl_cath_slice], sort=False, axis=0).to_csv(
            N2_dest_dir.joinpath(f"Cdl_FIT_{N2_fn}.csv")
        )
        logger.info(f"Cdl fit saved to: Cdl_FIT_{N2_fn}.csv")
    plt.close()
