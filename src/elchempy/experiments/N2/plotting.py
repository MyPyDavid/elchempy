"""
plotting functions for N2 experiments.
"""


# std lib

import logging

logger = logging.getLogger(__name__)

## local
from elchempy.plotters.plot_helpers import PlotterMixin

## for developing and testing
# from elchempy.experiments._dev_datafiles._dev_fetcher import get_files

## constants
from elchempy.constants import EvRHE

### 3rd party
import pandas as pd
import matplotlib.pyplot as plt


#%%
class N2_Plotter(PlotterMixin):

    # sweep_type_mapper = {'ls' : {'anodic' : '-.', 'cathodic' : '-', 'chrono' : ':'}}

    def plot_all_scans_scanrate(
        self,
        plot_data=None,
        xlabel=EvRHE,
        ylabel="j_mA_cm2",
        max_seg_only=True,
        savepath=None,
    ):

        if not isinstance(plot_data, pd.DataFrame):
            plot_data = self.data_selection

        fig, ax = plt.subplots(figsize=(8, 8))

        for sr, sgr in plot_data.groupby("scanrate"):
            max_seg = max(sgr["Segment #"].unique())
            if max_seg_only:
                sgr = sgr.loc[sgr["Segment #"] == max_seg]
            for swp, swpgrp in sgr.groupby("SweepType"):
                ax.plot(
                    swpgrp[xlabel],
                    swpgrp[ylabel],
                    label=f"{swp} {sr} mV, seg {max_seg}",
                    linestyle=self.sweep_type_mapper["ls"][swp],
                )
        #                ax.legend(True)
        fig.suptitle(
            f"{self.filepath.parent.name}\n{self.filepath.name}\nData({len(plot_data)})"
        )
        ax.legend()
        ax.set_ylabel("$j \//\/mA/cm^{2}$")
        ax.set_xlabel("$E \//\/V_{RHE}$")
        if savepath:
            plt.savefig(savepath, dpi=100, bbox_inches="tight")
            plt.close()
            return

    def plot_Cdl_sweeptype_scatter(Cdl_pars, **kwargs):
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
