"""
plotting functions for N2 experiments.
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import logging

logger = logging.getLogger(__name__)


def N2_plot_Cdl_scans_scanrate(Cdl_scans, _savepath):
    EvRHE = "E_AppV_RHE"
    fig, ax = plt.subplots(figsize=(8, 8))

    for sr, sgr in Cdl_scans.groupby("ScanRate_mVs"):
        _seg = int(sgr["Segment #"].unique()[0])
        ax.plot(sgr[EvRHE], sgr["jmAcm-2"], label=f"{sr} mV, seg {_seg}")
    #                ax.legend(True)
    fig.suptitle(f"{Path(Cdl_scans.PAR_file.unique()[0]).name}")
    ax.legend()
    ax.set_ylabel("$j \//\/mA/cm^{2}$")
    ax.set_xlabel("$E \//\/V_{RHE}$")
    plt.savefig(_savepath, dpi=100, bbox_inches="tight")
    plt.close()


def N2_plot_Cdl_sweeptype_scatter(
    SampleID, ScanRates, Cdl_fit, Cdl_cath_slice, Cdl_an_slice, N2_dest_dir, N2_fn
):
    EvRHE = "E_AppV_RHE"
    fig, ax = plt.subplots()
    plt.title(
        "%s made with linear fit of\n %s (R=%.3f)"
        % (SampleID, ScanRates, Cdl_fit["Cdl_R"].mean())
    )

    Cdl_cath_slice.plot(
        x=EvRHE,
        y="Cdl_corr",
        kind="scatter",
        ylim=(0, 1.2 * Cdl_fit.Cdl.max()),
        color="orange",
        ax=ax,
        label="Cdl_Cath_corr",
    )
    Cdl_cath_slice.plot(
        x=EvRHE,
        y="Cdl",
        kind="scatter",
        ylim=(0, 1.2 * Cdl_fit.Cdl.max()),
        color="r",
        ax=ax,
        label="Cdl_Cath",
    )
    Cdl_an_slice.plot(
        x=EvRHE,
        y="Cdl_corr",
        kind="scatter",
        ylim=(0, 1.2 * Cdl_fit.Cdl.max()),
        color="c",
        ax=ax,
        label="Cdl_Anod_corr",
    )
    Cdl_an_slice.plot(
        x=EvRHE,
        y="Cdl",
        kind="scatter",
        ylim=(0, 1.2 * Cdl_fit.Cdl.max()),
        ax=ax,
        label="Cdl_Anod",
    )

    plt.savefig(N2_dest_dir.joinpath(f"Cdl_{N2_fn}.png"), dpi=100, bbox_inches="tight")
    pd.concat([Cdl_an_slice, Cdl_cath_slice], sort=False, axis=0).to_csv(
        N2_dest_dir.joinpath(f"Cdl_FIT_{N2_fn}.csv")
    )
    logger.info(f"Cdl fit saved to: Cdl_FIT_{N2_fn}.csv")
    plt.close()
