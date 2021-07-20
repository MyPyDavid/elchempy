#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 11:04:07 2020

@author: zmg
"""

import itertools
from pathlib import Path
from datetime import datetime


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import statsmodels.formula.api as smf
import seaborn as sns

try:
    import SampleSelection
except ModuleNotFoundError:
    SampleSelection = lambda x: x


class PlotAnalysis:
    """Some properties of electrochemical experiments collected and reference values"""

    def __init__():
        pass

    @staticmethod
    def LinearTrend(df, x_col, y_col):
        data = pd.DataFrame([])
        #    colX,colY = OCP_DATA['Elapsed Time(s)'],OCP_DATA['E(V)']
        data["x"] = x_col
        data["y"] = y_col
        res = smf.ols("data.y ~ data.x", data=data).fit()
        #    X = sm.add_constant(x), mod = sm.OLS(y,X),res = mod.fit()
        intercept, slope, rsquared = res.params[0], res.params[1], res.rsquared

        if np.abs(slope) < 1e-5:
            check = True
        else:
            check = False
        return intercept, slope, rsquared, check

    def heatmap(x, y, **kwargs):
        if "color" in kwargs:
            color = kwargs["color"]
        else:
            color = [1] * len(x)

        if "palette" in kwargs:
            palette = kwargs["palette"]
            n_colors = len(palette)
        else:
            n_colors = 256  # Use 256 colors for the diverging color palette
            palette = sns.color_palette("Blues", n_colors)

        if "color_range" in kwargs:
            color_min, color_max = kwargs["color_range"]
        else:
            color_min, color_max = min(color), max(
                color
            )  # Range of values that will be mapped to the palette, i.e. min and max possible correlation

        def value_to_color(val):
            if color_min == color_max:
                return palette[-1]
            else:
                val_position = float((val - color_min)) / (
                    color_max - color_min
                )  # position of value in the input range, relative to the length of the input range
                val_position = min(
                    max(val_position, 0), 1
                )  # bound the position betwen 0 and 1
                ind = int(
                    val_position * (n_colors - 1)
                )  # target index in the color palette
                return palette[ind]

        if "size" in kwargs:
            size = kwargs["size"]
        else:
            size = [1] * len(x)

        if "size_range" in kwargs:
            size_min, size_max = kwargs["size_range"][0], kwargs["size_range"][1]
        else:
            size_min, size_max = min(size), max(size)

        size_scale = kwargs.get("size_scale", 500)

        def value_to_size(val):
            if size_min == size_max:
                return 1 * size_scale
            else:
                val_position = (val - size_min) * 0.99 / (
                    size_max - size_min
                ) + 0.01  # position of value in the input range, relative to the length of the input range
                val_position = min(
                    max(val_position, 0), 1
                )  # bound the position betwen 0 and 1
                return val_position * size_scale

        if "x_order" in kwargs:
            x_names = [t for t in kwargs["x_order"]]
        else:
            x_names = [t for t in sorted(set([v for v in x]))]
        x_to_num = {p[1]: p[0] for p in enumerate(x_names)}

        if "y_order" in kwargs:
            y_names = [t for t in kwargs["y_order"]]
        else:
            y_names = [t for t in sorted(set([v for v in y]))]
        y_to_num = {p[1]: p[0] for p in enumerate(y_names)}

        plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1)  # Setup a 1x10 grid
        ax = plt.subplot(
            plot_grid[:, :-1]
        )  # Use the left 14/15ths of the grid for the main plot

        marker = kwargs.get("marker", "s")

        kwargs_pass_on = {
            k: v
            for k, v in kwargs.items()
            if k
            not in [
                "color",
                "palette",
                "color_range",
                "size",
                "size_range",
                "size_scale",
                "marker",
                "x_order",
                "y_order",
            ]
        }

        ax.scatter(
            x=[x_to_num[v] for v in x],
            y=[y_to_num[v] for v in y],
            marker=marker,
            s=[value_to_size(v) for v in size],
            c=[value_to_color(v) for v in color],
            **kwargs_pass_on,
        )
        ax.set_xticks([v for k, v in x_to_num.items()])
        ax.set_xticklabels(
            [k for k in x_to_num], rotation=45, horizontalalignment="right"
        )
        ax.set_yticks([v for k, v in y_to_num.items()])
        ax.set_yticklabels([k for k in y_to_num])

        ax.grid(False, "major")
        ax.grid(True, "minor")
        ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
        ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

        ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
        ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
        ax.set_facecolor("#F1F1F1")

        # Add color legend on the right side of the plot
        if color_min < color_max:
            ax = plt.subplot(plot_grid[:, -1])  # Use the rightmost column of the plot

            col_x = [0] * len(palette)  # Fixed x coordinate for the bars
            bar_y = np.linspace(
                color_min, color_max, n_colors
            )  # y coordinates for each of the n_colors bars

            bar_height = bar_y[1] - bar_y[0]
            ax.barh(
                y=bar_y,
                width=[5] * len(palette),  # Make bars 5 units wide
                left=col_x,  # Make bars start at 0
                height=bar_height,
                color=palette,
                linewidth=0,
            )
            ax.set_xlim(
                1, 2
            )  # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
            ax.grid(False)  # Hide grid
            ax.set_facecolor("white")  # Make background white
            ax.set_xticks([])  # Remove horizontal ticks
            ax.set_yticks(
                np.linspace(min(bar_y), max(bar_y), 3)
            )  # Show vertical ticks for min, middle and max
            ax.yaxis.tick_right()  # Show vertical ticks on the right

    def corrplot(data, size_scale=500, marker="s"):
        #        data = rcorr
        corr = pd.melt(data.reset_index(), id_vars="index").dropna()
        #        corr = corr.loc[corr.x.isin(SampleSelection.EC_EIS_par_cols)]
        corr.columns = ["x", "y", "value"]
        PlotAnalysis.heatmap(
            corr["x"],
            corr["y"],
            color=corr["value"],
            color_range=[-1, 1],
            palette=sns.diverging_palette(20, 220, n=256),
            size=corr["value"].abs(),
            size_range=[0, 1],
            marker=marker,
            x_order=data.columns,
            y_order=data.columns[::-1],
            size_scale=size_scale,
        )


def plot_spectra_SampleIDs(xc, yc, spectras, destdir):
    #    stgr = EIS_O2_065_acid1_no
    #    stgr = grE
    #     fig,ax = plt.subplots()
    #     spdf_char.plot(x=['DATA_Yre'],y=['DATA_Yim'],kind='scatter',ax=ax)
    #     spdf_char.plot(x=['FIT2_Yre'],y=['FIT2_Yim'],kind='scatter',ax=ax)
    spdf_char = spectras
    x_col, y_col = "DATA_Yre", "DATA_Yim"
    x_fit2, y_fit2 = "FIT2_Yre", "FIT2_Yim"
    #    Y_adm = {'x_data' : }
    AdmImp = eisplot.AmdImp()
    #    x_col, y_col = 'N_content','Rct_kin'
    #    x_col, y_col = 'Qad','Rct'
    #    'BET_Area_RPT', 'BET_SA m2/g_RAW'
    #%%
    for ZY in AdmImp:
        ZYinfo = AdmImp[ZY]
        fig, ax = plt.subplots()
        EC_exp_nuniq = [
            i
            for i in SampleSelection.EC_exp_cols + ["E_RHE"]
            if spdf_char[i].nunique() == 1
        ]
        EC_exp_nuniq_non = [
            i
            for i in SampleSelection.EC_exp_cols + ["E_RHE"]
            if spdf_char[i].nunique() != 1
        ]
        EC_exp_cond1 = [
            "{0} : {1}".format(i, spdf_char[i].unique()[0]) for i in EC_exp_nuniq[0:3]
        ]
        EC_exp_cond_png = "_".join(
            ["{0}".format(spdf_char[i].unique()[0]) for i in EC_exp_nuniq]
        )
        EC_exp_cond2 = [
            "{0} : {1}".format(i, spdf_char[i].unique()[0]) for i in EC_exp_nuniq[3::]
        ]
        if "pH" in EC_exp_nuniq_non:
            extra_pH = "pH : {0}".format(
                " & ".join([str(i) for i in list(spdf_char["pH"].unique())])
            )
            EC_exp_cond2.append(extra_pH)

        ax.text(
            0.4,
            1.1,
            "{0} with {1} \n in {2}\n {3}".format(
                xc, yc, ", ".join(EC_exp_cond1), ", ".join(EC_exp_cond2)
            ),
            horizontalalignment="center",
            transform=ax.transAxes,
        )
        #    fig.suptitle('{0} with {1}\n in {2}'.format(x_col,y_col,', '.join(EC_exp_cond)))
        #                                    ym = stgr[y_col].max()*1.3 if stgr[y_col].max() < 2*stgr[y_col].mean() else stgr[y_col].mean()*1.5
        x_col, y_col = ZYinfo["x_data"], ZYinfo["y_data"]
        x_fit2, y_fit2 = ZYinfo["x_fit"], ZYinfo["y_fit"]

        ymean, ystd = spdf_char[y_col].mean(), spdf_char[y_col].std()
        y_fltr = spdf_char.loc[spdf_char[y_col] < ymean + 4 * ystd, y_col]
        yfltrmean, yfltrstd, yfltrmax = y_fltr.mean(), y_fltr.std(), y_fltr.max()
        xmean, xstd = spdf_char[x_col].mean(), spdf_char[x_col].std()
        x_fltr = spdf_char.loc[spdf_char[x_col] < xmean + 4 * xstd, x_col]
        xfltrmean, xfltrstd, xfltrmax, xfltrmin = (
            x_fltr.mean(),
            x_fltr.std(),
            x_fltr.max(),
            x_fltr.min(),
        )

        ymin, ymax = spdf_char[y_col].min(), spdf_char[y_col].max()
        ymax_set = (
            yfltrmax * 1.1
        )  # if not yfltrmax < yfltrmean+3*yfltrstd else yfltrmean*2
        ymin_set = ymin * 0.9
        if np.abs(xfltrmean) - 5 * np.abs(xfltrstd) < 0:
            if spdf_char.loc[spdf_char[x_col] < 0].empty:
                xmin_set = 0
            else:
                xmin_set = xfltrmin * 0.9
        else:
            xmin_set = xfltrmin * 0.9
        #    xmin_set = xfltrmin*0.9 if not stgr.loc[stgr[x_col] < 0].empty and not xfltrmean-6*xfltrstd < 0 else 0
        xmax_set = xfltrmax * 1.1
        #    MLmax = stgr.ML.max()*1.2
        EC_st_out = []
        OriginColors = Characterization_TypeSetting.OriginColorList()
        xls_path = "{0}_{1}_{2}.xlsx".format(x_col, y_col, ZY).replace("/", "")
        png_path = "{0}_{1}_{2}_{3}_spectra.png".format(
            xc, yc, EC_exp_cond_png, ZY
        ).replace("/", "")
        print(xls_path)
        for Parf, IDgr in spdf_char.groupby(by=["PAR_file"]):
            ID_cl1 = [i for i in IDgr.Colorcode.unique() if str(i) != "nan"]
            IDnm = IDgr.SampleID.unique()[0]
            if ID_cl1:
                ID_colorcode = ID_cl1[0]
            else:
                ID_colorcode = 1

            RGB = [
                float(i) / 255
                for i in (
                    OriginColors.loc[
                        OriginColors.OriginCode == ID_colorcode, "RGB"
                    ].values[0]
                ).split(",")
            ]
            x, y = IDgr[x_col].values, IDgr[y_col].values
            xfit2, yfit2 = IDgr[x_fit2].values, IDgr[y_fit2].values
            alpha_val_ML = IDgr.ML.unique()[0] / 21.4 + 0.2
            alpha_val_ML_set = alpha_val_ML if not alpha_val_ML > 1 else 1
            alpha_val_ML_set = alpha_val_ML_set if not alpha_val_ML_set < 0 else 0.1
            try:
                ax.scatter(
                    x,
                    y,
                    label="{0} ({1})".format(
                        IDgr.SampleLabel.unique()[0], IDgr.SampleID.unique()[0]
                    ),
                    color=RGB,
                    s=50,
                    alpha=alpha_val_ML_set,
                )
            except Exception as e:
                print("No plot for: {0}, because {1}".format(IDnm, e))
                if "RGB" in e:
                    print("RGB: {}, alpha = {:.3f} ".format(RGB, alpha_val_ML))
            try:
                alpha_val_ML_fit_set = (
                    alpha_val_ML_set + 0.2 if alpha_val_ML < 0.7 else alpha_val_ML_set
                )
                fit_label = "{0} ({1})".format(
                    IDgr.SampleLabel.unique()[0], IDgr.SampleID.unique()[0]
                )
                ax.plot(xfit2, yfit2, color=RGB, lw=2, alpha=alpha_val_ML_fit_set)
            except Exception as e:
                print("No plot for: {0}, because {1}".format(IDnm, e))
                if "RGB" in e:
                    print("RGB: {}, alpha = {:.3f} ".format(RGB, alpha_val_ML))
        #        EC_st_out.append(IDgr)
        # XLS file already there, otherwise: pd.concat([i for i in EC_st_out]).to_excel(destdir.joinpath(xls_path))
        ax.set_ylim(0, ymax_set)
        ax.set_xlim(0, xmax_set)
        ax.set_xlabel(ZYinfo["x_label"])
        ax.set_ylabel(ZYinfo["y_label"])
        ax.legend(
            bbox_to_anchor=(-0.35, 1.45, 1.7, 0.102),
            loc="lower left",
            ncol=2,
            mode="expand",
            borderaxespad=0.0,
        )
        #%%
        plt.savefig(destdir.joinpath(png_path), dpi=300, bbox_inches="tight")
        #    plt.show()
        plt.close()


def plot_pd_SampleIDs(stgr, x_col, y_col, corr_val, ddir, **kwargs):
    #    stgr = EIS_O2_065_acid1_no
    #    stgr,kwargs = grE, {}
    #    x_col,y_col = xc, yc
    #    x_col, y_col = 'N_content','Rct_kin'
    #    x_col, y_col = 'Qad','Rct'
    #    'BET_Area_RPT', 'BET_SA m2/g_RAW'
    before, b_samples = len(stgr), stgr.SampleID.unique()
    stgr = stgr.dropna(how="any", subset=[x_col, y_col])
    after, after_samples = len(stgr), stgr.SampleID.unique()
    missing_samples = [i for i in b_samples if i not in after_samples]
    ddir.mkdir(exist_ok=True, parents=True)
    if corr_val == None or corr_val == 0:
        try:
            corr_val = stgr[[x_col, y_col]].corr().iloc[1, 0]
        except Exception as e:
            print("Error trying to set corr val, {0}".format(e))

    fig, ax = plt.subplots()
    if any([i in SampleSelection.EC_exp_cols + ["E_RHE"] for i in stgr.columns]):
        EC_exp_nuniq = [
            i for i in SampleSelection.EC_exp_cols + ["E_RHE"] if stgr[i].nunique() == 1
        ]
        EC_exp_nuniq_non = [
            i for i in SampleSelection.EC_exp_cols + ["E_RHE"] if stgr[i].nunique() != 1
        ]
        EC_exp_cond1 = [
            "{0} : {1}".format(i, stgr[i].unique()[0]) for i in EC_exp_nuniq[0:3]
        ]
        EC_exp_cond_png = "_".join(
            ["{0}".format(stgr[i].unique()[0]) for i in EC_exp_nuniq]
        )
        EC_exp_cond2 = [
            "{0} : {1}".format(i, stgr[i].unique()[0]) for i in EC_exp_nuniq[3::]
        ]
        if "pH" in EC_exp_nuniq_non:
            extra_pH = "pH : {0}".format(
                " & ".join([str(i) for i in list(stgr["pH"].unique())])
            )
            EC_exp_cond2.append(extra_pH)
    else:
        EC_exp_cond1, EC_exp_cond2, EC_exp_cond_png = "", "", ""
    extra_text = "{0} with {1} (c={2})\n in {3}\n {4}".format(
        x_col,
        y_col,
        np.round(corr_val, 2),
        ", ".join(EC_exp_cond1),
        ", ".join(EC_exp_cond2),
    )
    #    if corr_val == 0:
    #        extra_text = '{0} with {1} \n in {2}\n {3}'.format(x_col,y_col,', '.join(EC_exp_cond1),', '.join(EC_exp_cond2))
    ax.text(0.4, 1.1, extra_text, horizontalalignment="center", transform=ax.transAxes)
    #    fig.suptitle('{0} with {1}\n in {2}'.format(x_col,y_col,', '.join(EC_exp_cond)))
    #                                    ym = stgr[y_col].max()*1.3 if stgr[y_col].max() < 2*stgr[y_col].mean() else stgr[y_col].mean()*1.5
    ymean, ystd = stgr[y_col].mean(), stgr[y_col].std()
    y_fltr = stgr.loc[stgr[y_col] < ymean + 4 * ystd, y_col]
    yfltrmean, yfltrstd, yfltrmax, yfltrmin = (
        y_fltr.mean(),
        y_fltr.std(),
        y_fltr.max(),
        y_fltr.min(),
    )
    xmean, xstd = stgr[x_col].mean(), stgr[x_col].std()
    x_fltr = stgr.loc[stgr[x_col] < xmean + 4 * xstd, x_col]
    xfltrmean, xfltrstd, xfltrmax, xfltrmin = (
        x_fltr.mean(),
        x_fltr.std(),
        x_fltr.max(),
        x_fltr.min(),
    )

    ymin, ymax = stgr[y_col].min(), stgr[y_col].max()
    ymax_set = (
        yfltrmax + 2 * yfltrstd
    )  # if not yfltrmax < yfltrmean+3*yfltrstd else yfltrmean*2
    ymin_set = yfltrmin - 2 * yfltrstd
    if ymin_set < 0:
        if stgr.loc[stgr[y_col] < 0].empty:
            ymin_set = 0
    if np.abs(xfltrmean) - 5 * np.abs(xfltrstd) < 0:
        if stgr.loc[stgr[x_col] < 0].empty:
            xmin_set = 0
        else:
            xmin_set = xfltrmin - 2 * xfltrstd
    else:
        xmin_set = xfltrmin - 2 * xfltrstd
    #    xmin_set = xfltrmin*0.9 if not stgr.loc[stgr[x_col] < 0].empty and not xfltrmean-6*xfltrstd < 0 else 0
    xmax_set = xfltrmax + 2 * xfltrstd
    #    MLmax = stgr.ML.max()*1.2
    EC_st_out = []
    OriginColors = Characterization_TypeSetting.OriginColorList()
    xls_path = "{0}_{1}.xlsx".format(x_col, y_col).replace("/", "")
    png_path = "{0}_{1}_{2}_{3}.png".format(
        x_col, y_col, EC_exp_cond_png, kwargs.get("corr_method", "")
    ).replace("/", "")
    print(xls_path)
    for IDnm, IDgr in stgr.groupby(by=["SampleID"]):
        ID_cl1 = [i for i in IDgr.Colorcode.unique() if str(i) != "nan"]
        if ID_cl1:
            ID_colorcode = ID_cl1[0]
        else:
            ID_colorcode = 1
        RGB = [
            float(i) / 255
            for i in (
                OriginColors.loc[OriginColors.OriginCode == ID_colorcode, "RGB"].values[
                    0
                ]
            ).split(",")
        ]
        if len(IDgr) == 1:
            x_mean = IDgr[x_col].values[0]
            y_mean, y_std = IDgr[y_col].values[0], 0.0001
        else:
            x_mean = IDgr[x_col].mean()
            y_mean, y_std = IDgr[y_col].mean(), IDgr[y_col].std()
        if "E_RHE" in x_col:
            x_mean, y_mean = IDgr[x_col].values, IDgr[y_col].values
            xmax_set, xmin_set = 1.1, 0
        if pd.isna(IDgr.ML.unique()[0]):
            alpha_val_ML_set = 0.8
        else:
            alpha_val_ML = IDgr.ML.unique()[0] / 21.4 + 0.2
            alpha_val_ML_set = alpha_val_ML if not alpha_val_ML > 1 else 1
            alpha_val_ML_set = alpha_val_ML_set if not alpha_val_ML_set < 0 else 0.1
        try:
            ax.scatter(
                x_mean,
                y_mean,
                label="{0} ({1})".format(
                    IDgr.SampleLabel.unique()[0], IDgr.SampleID.unique()[0]
                ),
                color=RGB,
                s=200,
                alpha=alpha_val_ML_set,
                marker=Characterization_TypeSetting.Series_Colors(
                    IDgr.SeriesID.unique()[0]
                )["marker"],
            )
        #            ax.errorbar(x_mean,y_mean,yerr=y_std, label='{0} ({1})'.format(IDgr.SampleLabel.unique()[0], IDgr.SampleID.unique()[0]),color=RGB,
        #                       alpha=alpha_val_ML_set)
        except Exception as e:
            print("No plot for: {0}, because {1}".format(IDnm, e))
            if "RGB" in str(e):
                print("RGB: {}, alpha = {:.3f} ".format(RGB, alpha_val_ML))
    # XLS file already there, otherwise: pd.concat([i for i in EC_st_out]).to_excel(ddir.joinpath(xls_path))
    ax.set_ylim(ymin_set, ymax_set)
    if np.log10(np.abs(ymax_set) - np.abs(ymin_set)) > 2.5:
        ax.set_yscale("log")
        ax.autoscale(True)
    ax.set_xlim(xmin_set, xmax_set)
    if np.log10(np.abs(xmax_set) - np.abs(xmin_set)) > 2.5:
        pass
    #        ax.set_xscale('log')
    #        ax.set_xlim(xmax_set-700,xmax_set)
    #                                    if log_reg.rsquared > reg.rsquared:
    #                                        ax.autoscale(True)
    #                                        ax.set_yscale('log')
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    #                                    ax.autoscale(True)
    if "logx" in kwargs.keys():
        ax.set_xscale("log")

    ax.legend(
        bbox_to_anchor=(-0.35, 1.45, 1.80, 0.122),
        loc="lower left",
        ncol=2,
        mode="expand",
        fontsize=14,
        borderaxespad=0.0,
    )
    plt.savefig(ddir.joinpath(png_path), dpi=300, bbox_inches="tight")
    #    plt.show()
    plt.close()


class Characterization_TypeSetting:
    def __init__(self):
        self.jlbl = r"$\mathrm{j\/ / \/ mAcm^{-2}}$"
        self.Erhe = r"$\mathrm{E\/ / \/ V_{RHE}$"

    @staticmethod
    def Series_Colors(serie, type=str):
        series = [
            "CB3",
            "CB4",
            "CB6",
            "CB7",
            "CB_pyr",
            "CB_bim",
            "Porph_SiO2",
            "Co_PANI",
            "C3N4-Plu",
            "nan",
            "all",
        ]
        colors = [
            "aqua",
            "cadetblue",
            "deepskyblue",
            "steelblue",
            "grey",
            "darkviolet",
            "sandybrown",
            "lime",
            "gold",
            "black",
            "black",
        ]
        mt = ["o", "^", "s", "P", "p", "X", "*", "d", "H", "D", "D"]
        series_formatting_ref = {
            x: {"color": y, "marker": z} for x, y, z in zip(series, colors, mt)
        }
        return series_formatting_ref.get(serie, {"color": "pink", "marker": "x"})

    @staticmethod
    def Sample_Metal_Color(label, type=str):
        metals = ("FeTMPP", "FeTPP", "Mn", "Co", "Cu", "H2")
        mc = dict(
            zip(metals, ("red", "orange", "magenta", "dodgerblue", "orangered", "lime"))
        )
        match = [i for i in metals if i in label and not "Free" in label]
        fill_style = "full"
        if len(match) == 1:
            c = mc[match[0]]
            c_alt = c
        elif len(match) > 1:
            c = mc[match[0]]
            c_alt = mc[match[1]]
            fill_style = "left"
        elif len(match) < 1:
            c = "grey"
            c_alt = c
        return c, c_alt, fill_style

    def Marker_Style(label):
        c, c_alt, fillst = Characterization_TypeSetting.Sample_Metal_Color(label)
        marker_style = dict(
            color=c,
            linestyle=":",
            marker="o",
            markersize=80,
            markerfacecoloralt=c_alt,
            fillstyle=fillst,
        )
        return marker_style

    def Scatter_Marker_Style(label, **kwargs):
        c, c_alt, fillst = Characterization_TypeSetting.Sample_Metal_Color(label)
        mt = "<"
        if "serie" in kwargs.keys():
            mt = Characterization_TypeSetting.Series_Colors(kwargs["serie"])[1][
                kwargs["serie"]
            ]
        marker_style = dict(color=c, linestyle=":", marker=mt, s=200)
        return marker_style


#    def OriginColorList():
#        try:
#            OriginColorList= FileHelper.FindExpFolder().LoadOriginColor()
#        except:
#            print('Print no Orginicolorlist')
#        return OriginColorList


class eisplot:
    parlst = ["Rct", "Rs", "Rorr", "Rct_kin", "Cdlp", "Qad", "nDL", "nAd", "Qad+Cdlp"]
    extrapars = ["R3", "Q3", "n3", "R3_kin"]
    fitparlst = ["Rct", "Rs", "Rorr", "Cdlp", "Qad", "nDL", "nAd"]

    def __init__(self, par, **kwargs):
        self.par = par
        if "Gas" in par:
            self.ms = eisplot.gas_plot(self)
        #       ELIF  par in eisplot.parlst + eisplot.extrapars:
        else:
            self.ylim, self.logy, self.logyscale = eisplot.pars_plot_lims(self)
            if "yvalues" in kwargs.keys():
                ymax = kwargs["yvalues"].max()
                print(ymax)
                if np.max(self.ylim) * 3 < ymax:
                    self.ylim = (np.min(self.ylim), ymax * 1.2)
                if ymax * 1.5 < np.max(self.ylim):
                    self.ylim = (np.min(self.ylim), ymax * 1.2)

    def pars_plot_lims(self):
        yPar = self.par
        if yPar == "Rs":
            ylim, logy = (0, 60), False
        elif yPar == "Rct":
            ylim, logy = (1e-2, 5e4), True
        elif yPar == "Rct_kin":
            ylim, logy = (0, 2), False
        elif yPar == "Rorr":
            ylim, logy = (1, 5e9), True
        elif yPar == "Qad":
            ylim, logy = (0, 300e-4), False
        elif yPar == "Q3":
            ylim, logy = (0, 1e-1), False
        elif yPar == "nAd":
            ylim, logy = (0.2, 1), False
        elif yPar == "n3":
            ylim, logy = (0.1, 1), False
        elif yPar == "Cdlp":
            ylim, logy = (0, 1e-4), False
        elif yPar == "Ls":
            ylim, logy = (1e-06, 1e-4), False
        elif yPar == "nDL":
            ylim, logy = (0.2, 1), False
        elif yPar == "Qad+Cdlp":
            ylim, logy = (0, 1e-1), False
        elif yPar == "Rorr_kin":
            ylim, logy = (0, 5e-3), False
        elif yPar == "Aw":
            ylim, logy = (1e-1, 5e4), True
        elif yPar == "tau":
            ylim, logy = (1e-02, 1e3), True
        elif yPar == "R3":
            ylim, logy = (1, 5e3), True
        elif yPar == "R_ion":
            ylim, logy = (1, 3e3), True
        elif yPar == "R3_kin":
            ylim, logy = (0, 0.1), False
        elif yPar == "J_diff_lim":
            ylim, logy = (0, -6.2), False
        elif yPar == "Jkin_075":
            ylim, logy = (0.1, 100), True
        elif yPar == "Jkin_080":
            ylim, logy = (0.1, 100), True
        elif yPar == "E_onset":
            ylim, logy = (0.4, 1), False
        elif yPar == "E_half":
            ylim, logy = (0.4, 1), False
        elif yPar == "TSa_l":
            ylim, logy = (10, 200), False
        elif yPar == "TSb_l":
            ylim, logy = (0, 2), False
        elif yPar == "TSa_h":
            ylim, logy = (10, 400), False
        elif yPar == "TSb_h":
            ylim, logy = (0, 2), False
        elif yPar == "Jring_050":
            ylim, logy = (0, 0.1), False
        elif yPar == "FracH2O2_050":
            ylim, logy = (0, 50), False
        else:
            ylim, logy, logyscale = (0, 1e3), True, "log"
        if logy == False:
            logyscale = "linear"
        else:
            logyscale = "log"
        return ylim, logy, logyscale

    def dict_plot_lims(tt):
        _dict = {
            i.split("== ")[-1][1:-2]: {
                "ylim": tt.split("\n")[n + 1].split("=")[-1].split()[0:2]
            }
            for (n, i) in enumerate(tt.split("\n"))
            if "==" in i
        }

    def gas_plot(self):
        gas = self.par
        if "N2" in gas:
            ms = "^"
            alpha = 0.6
        elif "O2" in gas:
            ms = "o"
            alpha = 0.9
        else:
            ms = "*"
            alpha = 0.4
        #        return ms, alpha
        return {"marker": ms, "alpha": alpha}

    def AmdImp():
        x_col, y_col = "DATA_Yre", "DATA_Yim"
        x_fit2, y_fit2 = "FIT2_Yre", "FIT2_Yim"
        Y_adm = {
            "x_data": x_col,
            "y_data": y_col,
            "x_fit": x_fit2,
            "y_fit": y_fit2,
            "x_label": "$\mathrm{Y_{Re}}$",
            "y_label": "$\mathrm{Y_{Im}}$",
        }
        Z_imp = {
            "x_data": "DATA_Zre",
            "y_data": "DATA_-Zim",
            "x_fit": "FIT2_Zre",
            "y_fit": "FIT2_Zim",
            "x_label": "$\mathrm{Z_{Re}}$",
            "y_label": "$\mathrm{Z_{Im}}$",
        }
        return {"Y": Y_adm, "Z": Z_imp}

    def ReFit():
        bad_fittings = ["N2_EIS-range_1500rpm_JOS3_high-load_267"]

    def StandardMeasurements():
        files = []

    def read_varnames(grp):
        lmfit_var_names_cols = [i for i in grp.columns if "lmfit_var_names" in i][0]
        _varsgrp = [
            a for i in grp[lmfit_var_names_cols].unique() for a in i.split(", ")
        ]
        _varsgrp += list(set(["Rct_kin" for i in _varsgrp if "Rct" in i])) + list(
            set(
                [
                    "Qad+Cdlp"
                    for i in _varsgrp
                    if all([i in _varsgrp for i in ["Qad", "Cdlp"]])
                ]
            )
        )
        _vars_in_grp = [
            i
            for i in grp.columns
            if any(v in i for v in _varsgrp) and not "stderr" in i
        ]
        return _varsgrp, _vars_in_grp

    def EIS_ParsPlotting_Rs_per_Sample(SampleData, PDDirEIScom):
        for yPar in ["Rct", "Rs", "Rorr", "Rct_kin", "Cdlp"]:
            fig, ax = plt.subplots(1, 1)
            SampleLabel, SampleID = (
                SampleData.SampleCode.unique()[0],
                SampleData.SampleID.unique()[0],
            )
            for Gas, gasGr in SampleData.groupby(by="Gas"):
                for Status, stGr in gasGr.groupby(by="postAST"):
                    if yPar == "Cdlp" or yPar == "Qad":
                        ax2 = ax.twinx()
                        sc1 = ax.scatter(
                            stGr["E_AppV_RHE"].values,
                            stGr["Cdlp"].values,
                            label=str("Cdl" + ": " + "%s, %s" % (Gas, Status)),
                            s=80,
                            marker="D",
                        )
                        sc2 = ax.scatter(
                            stGr["E_AppV_RHE"].values,
                            stGr["Qad"].values,
                            label=str("Qad" + ": " "%s, %s" % (Gas, Status)),
                            s=80,
                            marker="h",
                        )
                        sc_12 = ax.scatter(
                            stGr["E_AppV_RHE"].values,
                            stGr["Qad"].values + stGr["Cdlp"].values,
                            label=str("Cdl+Qad" + ": " + "%s, %s" % (Gas, Status)),
                            s=80,
                            marker="h",
                        )

                        #            sc1 = ax1.scatter(Pgr['E_AppV_RHE'].values,Pgr[yPar[0]].values,label=str(yPar[0]+': '+Path(fn).stem),s=80)
                        sc3 = ax2.scatter(
                            stGr["E_AppV_RHE"].values,
                            stGr["nDL"].values,
                            label=str("n_Cdl" + ": " + "%s, %s" % (Gas, Status)),
                            s=40,
                            alpha=0.6,
                            marker="o",
                        )
                        sc4 = ax2.scatter(
                            stGr["E_AppV_RHE"].values,
                            stGr["nAd"].values,
                            label=str("n_Qad" + ": " + "%s, %s" % (Gas, Status)),
                            s=40,
                            alpha=0.6,
                            marker="h",
                        )
                        ax.set_ylim([0, 0.03])
                        ax2.set_ylim([0, 1])
                        ax.legend(bbox_to_anchor=(0.5, 1.5), ncol=2, loc="upper center")
                    else:
                        ax.scatter(
                            stGr["E_AppV_RHE"].values,
                            stGr[yPar].values,
                            label="%s, %s" % (Gas, Status),
                            s=80,
                        )

            ax.legend(ncol=1, loc="upper left", fontsize=10)
            ax.set_ylabel(yPar)
            ax.set_xlabel("E / V v RHE")
            ax.set_title(SampleLabel)
            ax.grid(True)
            DestFile = PDDirEIScom.joinpath(
                "_".join([SampleID, yPar, SampleLabel])
            ).with_suffix(".png")
            plt.savefig(DestFile, dpi=300, bbox_inches="tight")
            plt.close()

    #
    def EIS_ParsPlotting(AllData_E_file, Cdl, DestFile, SaveFigsC=False):
        SampleLabel = AllData_E_file.Sample.unique()[0]
        fig1, ax1 = plt.subplots(1, 1)
        ax2 = ax1.twinx()
        minl, maxl = (
            AllData_E_file[["Cdlp", "Qad"]].min().min(),
            1.1 * AllData_E_file[["Cdlp", "Qad"]].max().max(),
        )
        #     for yPar in [('Cdlp','nDL'),('Qad','nAd')]:
        #        minl, maxl = 0.5*AllData_E_file[yPar[0]].min(),1.1*AllData_E_file[yPar[0]].max()
        scts = []
        for fn, Pgr in AllData_E_file.groupby("File"):
            sc1 = ax1.scatter(
                Pgr["E_AppV_RHE"].values,
                Pgr["Cdlp"].values,
                label=str("Cdl" + ": " + Path(fn).stem),
                s=80,
                c="royalblue",
                marker="D",
            )
            sc2 = ax1.scatter(
                Pgr["E_AppV_RHE"].values,
                Pgr["Qad"].values,
                label=str("Qad" + ": " + Path(fn).stem),
                s=80,
                c="tomato",
                marker="h",
            )
            sc_12 = ax1.scatter(
                Pgr["E_AppV_RHE"].values,
                Pgr["Qad"].values + Pgr["Cdlp"].values,
                label=str("Cdl+Qad" + ": " + Path(fn).stem),
                s=80,
                c="darkviolet",
                marker="h",
            )
            #            sc1 = ax1.scatter(Pgr['E_AppV_RHE'].values,Pgr[yPar[0]].values,label=str(yPar[0]+': '+Path(fn).stem),s=80)
            sc3 = ax2.scatter(
                Pgr["E_AppV_RHE"].values,
                Pgr["nDL"].values,
                label=str("n_Cdl" + ": " + Path(fn).stem),
                s=40,
                c="lightgrey",
                alpha=0.6,
                marker="o",
            )
            sc4 = ax2.scatter(
                Pgr["E_AppV_RHE"].values,
                Pgr["nAd"].values,
                label=str("n_Qad" + ": " + Path(fn).stem),
                s=40,
                c="black",
                alpha=0.6,
                marker="h",
            )
        #            sc2 = ax2.scatter(Pgr['E_AppV_RHE'].values,Pgr[yPar[1]].values,label=str(yPar[1]+': '+Path(fn).stem),s=40,c='lightgrey',alpha=0.5, marker='^')
        scts.append([sc1, sc2, sc_12, sc3, sc4])
        if not Cdl.empty:
            sc5 = ax1.scatter(
                Cdl["E_AppV_RHE"].values,
                Cdl["Cdl"].values,
                label=str([Path(i).stem + "\n" for i in Cdl.Filename.unique()]),
                s=20,
                alpha=0.8,
                marker="o",
                c="aqua",
            )
            Cdl_CV_max = Cdl["Cdl"].mean() + 2 * Cdl["Cdl"].std()
            if Cdl_CV_max > maxl:
                maxl = Cdl_CV_max
            scts.append([sc5])
        #        ax1.set_ylabel(yPar[0])
        #        ax2.set_ylabel(yPar[1])
        ax1.set_ylabel("Cdl / Qad")
        ax2.set_ylabel("n_CdL / n_Qad")
        ax1.set_xlabel("E / V v RHE")
        ax1.set_ylim([0, maxl])
        ax2.set_ylim([0, 1])
        ax1.grid(True)
        #            lns = lns1+lns2+lns3
        flat_scts = list(itertools.chain.from_iterable(scts))
        labs = [l.get_label() for l in flat_scts]
        ax1.legend(
            flat_scts,
            labs,
            bbox_to_anchor=(0.4, 1.1 + 0.1 * len(labs)),
            ncol=1,
            loc="upper center",
            fontsize=10,
        )
        #    ax1.set_title('%s (%s)'%(SampleLabel,DestFile.stem))
        fig1.suptitle(("%s (%s)" % (SampleLabel, DestFile.stem)))
        #            ax2.legend(bbox_to_anchor=(0.4,1.24),ncol=1, loc="upper center",fontsize=10)
        if SaveFigsC:
            #            plt.savefig(dest_dir.joinpath('EIS_%s.png'%(yPar[0])),dpi=300,bbox_inches='tight')
            plt.savefig(
                DestFile.parent.joinpath("Cdl" + "_" + DestFile.name),
                dpi=300,
                bbox_inches="tight",
            )
        plt.show()
        plt.close()

    def PlotParsEIS(AllData_E_file, EISovv, SampleCode, DestFile, xEIS="Y", yEIS="Rct"):
        #    maxLim = (AllData_E_file[['DATA_Yim','DATA_Yre']].max()).max()
        #%%
        EvRHE = "E_AppV_RHE"
        maxYim = np.abs(AllData_E_file["DATA_%sim" % xEIS]).max()
        maxYre = np.abs(AllData_E_file["DATA_%sre" % xEIS]).max()
        #    AllData_E_file['DATA_Yre'].max()
        #    maxZim = AllData_E_file['DATA_Yim'].max()
        #    maxZre = AllData_E_file['DATA_Yre'].max()
        Lenrows = len(AllData_E_file[EvRHE].unique())
        #    fig,axes = plt.subplots(nrows=Lenrows ,sharex=True,sharey=True)
        ht, wd = 10, 15
        fig, ax = plt.subplots(figsize=(ht, wd))

        if SampleCode.empty:
            Scode = EISovv["SampleID"].unique()[0]
        else:
            Scode = SampleCode.Sample.values[0]

        fig.suptitle(
            "%s %s, in \n %s saturated %s \n %s \n %s"
            % (
                Scode,
                EISovv.postAST.values[0],
                EISovv["Gas"].unique()[0],
                EISovv["Electrolyte"].unique()[0],
                EISovv.EXP_date.values[0],
                Path(EISovv["SourceFilename"].unique()[0]).stem,
            )
        )
        dataC, fitC, extraC, initC = "tab:blue", "tab:red", "gold", "gray"
        #    ax.set_xlim(0,maxYre)
        #    ax.set_ylim(0,maxYre+maxYim*Lenrows)
        ax.grid(True)
        #    ax.axis('equal')
        for En, Ev in enumerate(AllData_E_file[EvRHE].unique()):
            Edata = AllData_E_file.loc[(AllData_E_file[EvRHE] == Ev)].sort_values(
                by="Frequency(Hz)"
            )
            if xEIS == "Y":
                FIT2_Im = Edata["FIT2_%sim" % xEIS].values + len(Edata) * [maxYim * En]
                FIT1_Im = Edata["FIT1_%sim" % xEIS].values + len(Edata) * [maxYim * En]
                DATA_Im = Edata["DATA_%sim" % xEIS].values + len(Edata) * [maxYim * En]
                xText = 0.06
            elif xEIS == "Z":
                FIT2_Im = np.abs(Edata["FIT2_%sim" % xEIS].values) + len(Edata) * [
                    maxYim * En
                ]
                FIT1_Im = np.abs(Edata["FIT1_%sim" % xEIS].values) + len(Edata) * [
                    maxYim * En
                ]
                DATA_Im = np.abs(Edata["DATA_%sim" % xEIS].values) + len(Edata) * [
                    maxYim * En
                ]
                xText = maxYre * 0.9

            ax.plot(Edata["FIT2_%sre" % xEIS].values, FIT2_Im, c=fitC, lw=2.5)
            ax.plot(
                Edata["FIT1_%sre" % xEIS].values,
                FIT1_Im,
                c="lightgrey",
                lw=2.5,
                alpha=0.5,
            )
            ax.scatter(Edata["DATA_%sre" % xEIS].values, DATA_Im, c=dataC, s=150)
            ax.annotate(
                "$\mathrm{%.2f \/ V_{RHE} }$" % Ev,
                xy=(xText, 0.001 + maxYim * En),
                xycoords="data",
            )
        ax.set_ylim(0, maxYim * (En + 2))
        if xEIS == "Y":
            ax.set_xlim(0, 0.07)
            xunit = "mS"
        elif xEIS == "Z":
            ax.set_xlim(0, maxYre * 1.1)
            xunit = "\Omega"
        ax.set_ylabel("$\mathrm{%s_{Im}\/ offset \//\/ %s}$" % (xEIS, xunit))
        ax.set_xlabel("$\mathrm{%s_{Re} \//\/ %s}$" % (xEIS, xunit))
        #    *(ht/wd)
        #    fig_path = EIS_dest_dir.with_suffix('.png')
        #    plt.show()
        #%%
        plt.savefig(DestFile, bbox_inches="tight", dpi=200)
        plt.close()

    def PlotCombinedEIS(
        AllData_E_file_spectras, EISovv, SampleCode, DestFile, xEIS="Y"
    ):
        #    maxLim = (AllData_E_file[['DATA_Yim','DATA_Yre']].max()).max()
        # AllData_E_file_spectras,EISovv,SampleCode,DestFile = spectras_comb_bestmods, bgrp,SampleCodes, DestFile
        #%%
        EvRHE = "E_AppV_RHE"

        if EISovv.Model_EEC.nunique() == 1:
            plot_Model_EEC = EISovv.Model_EEC.unique()[0]
        else:
            plot_Model_EEC = f"{EISovv.Model_EEC.unique()[0]} and more"
            # '; '.join(str(i) for i in EISovv.Model_EEC.unique())

        AllData_E_file = AllData_E_file_spectras
        # .loc[AllData_E_file_spectras.Model_EEC.isin(EISovv.Model_EEC.unique())]
        if xEIS in ["Y", "Z"]:
            colIm, colRe = f"{xEIS}im", f"{xEIS}re"
            maxYre = np.abs(AllData_E_file[f"DATA_{colRe}"]).max()

        elif "-Zangle" in xEIS:
            colIm, colRe = f"{xEIS}", "Frequency(Hz)"
            maxYre = np.abs(AllData_E_file[f"{colRe}"]).max()

        #    AllData_E_file['DATA_Yre'].max()
        #    maxZim = AllData_E_file['DATA_Yim'].max()
        #    maxZre = AllData_E_file['DATA_Yre'].max()
        # Lenrows = len(AllData_E_file[EvRHE].unique())
        #    fig,axes = plt.subplots(nrows=Lenrows ,sharex=True,sharey=True)
        ht, wd = 15, 20
        fig, ax = plt.subplots(figsize=(ht, wd))

        sID = EISovv["SampleID"].unique()[0]
        if not EISovv.empty:
            Scode = EISovv.query("SampleID == @sID").SampleCode.values[0]
        else:
            Scode = ""
        sID = f"{sID}({Scode})"

        def ___convert_to_datetime(d):
            return datetime.strptime(
                np.datetime_as_string(d, unit="s"), "%Y-%m-%dT%H:%M:%S"
            )

        fig.suptitle(
            f"{sID} {EISovv.postAST.unique()[0]}, in"
            + " \n"
            + f"{EISovv.Gas.unique()[0]} saturated {EISovv.Electrolyte.unique()[0]}"
            + " \n"
            + f"{___convert_to_datetime(EISovv.PAR_date.values[0]):%c}"
            + "\n"
            + f"{Path(EISovv.PAR_file.unique()[0]).stem}"
            + "\n"
            + f"{plot_Model_EEC}"
        )
        dataC, fitC, extraC, initC = "tab:blue", "tab:red", "gold", "gray"
        #    ax.set_xlim(0,maxYre)
        #    ax.set_ylim(0,maxYre+maxYim*Lenrows)
        ax.grid(True)
        #    ax.axis('equal')
        maxYim = np.abs(AllData_E_file[f"DATA_{colIm}"]).max()
        _splitter = "Segment #"
        _uniqsegs = AllData_E_file["Segment #"].nunique()
        _uniq = [
            k
            for k, val in AllData_E_file.nunique().to_dict().items()
            if val == _uniqsegs and any([k in c for c in ["E_AppV_RHE", "RPM_DAC"]])
        ]
        # E_AppV_RHE'

        for En, Ev in enumerate(AllData_E_file[_splitter].unique()):
            Edata = AllData_E_file.loc[(AllData_E_file[_splitter] == Ev)].sort_values(
                by="Frequency(Hz)"
            )
            _ERHE = Edata[EvRHE].unique()[0]
            _RPM = Edata["RPM_DAC"].unique()[0]
            _MOD = Edata["Model_EEC"].unique()[0]
            # print(_ERHE,_MOD)
            if xEIS == "Y":
                #                FIT2_Im = Edata['FIT2_%sim'%xEIS].values+len(Edata)*[maxYim*En]
                FIT1_Im = Edata[f"FIT_{colIm}"].values + len(Edata) * [maxYim * En]
                DATA_Im = Edata[f"DATA_{colIm}"].values + len(Edata) * [maxYim * En]
                xText = 0.06
                _xtra = 0.02
            elif xEIS == "Z":
                #                FIT2_Im = np.abs(Edata['FIT2_%sim'%xEIS].values)+len(Edata)*[maxYim*En]
                FIT1_Im = np.abs(Edata[f"FIT_{colIm}"].values) + len(Edata) * [
                    maxYim * En
                ]
                DATA_Im = np.abs(Edata[f"DATA_{colIm}"].values) + len(Edata) * [
                    maxYim * En
                ]
                xText = maxYre * 0.9
                _xtra = maxYre * 0.3
            elif xEIS == "-Zangle":
                FIT1_Im = np.abs(Edata[f"FIT_{colIm}"].values) + len(Edata) * [
                    maxYim * En
                ]
                DATA_Im = np.abs(Edata[f"DATA_{colIm}"].values) + len(Edata) * [
                    maxYim * En
                ]
                xText = maxYre * 0.9
                _xtra = maxYre * 0.3

            #            ax.plot(Edata['FIT2_%sre' %xEIS].values,FIT2_Im,c=fitC,lw=2.5)
            if not "Frequency" in colRe:
                ax.plot(
                    Edata[f"FIT_{colRe}"].values, FIT1_Im, c=fitC, lw=3.5, alpha=0.7
                )
                ax.scatter(Edata[f"DATA_{colRe}"].values, DATA_Im, c=dataC, s=150)
            else:
                ax.plot(Edata[f"{colRe}"].values, FIT1_Im, c=fitC, lw=3.5, alpha=0.7)
                ax.scatter(Edata[f"{colRe}"].values, DATA_Im, c=dataC, s=150)
            ax.annotate(
                "$\mathrm{%.2f \/ V_{RHE}, %.0f }$  %s" % (_ERHE, _RPM, _MOD),
                xy=(xText, 0.01 + maxYim * En),
                xycoords="data",
            )
            # ax.annotate(f'{_MOD}',xy=(xText-(0.5*_xtra), 0.1+maxYim*En), xycoords='data',
            # bbox=dict(boxstyle="round", fc="0.8",alpha=0.4))
        ax.set_ylim(0, maxYim * (En + 2))
        if xEIS == "Y":
            ax.set_xlim(0, 0.07)
            xunit = "mS"
        elif xEIS == "Z":
            ax.set_xlim(0, maxYre * 1.1)
            xunit = "\Omega"
        elif xEIS == "-Zangle":
            xunit = "Hz"
            ax.set_xscale("log")

        ax.set_ylabel("$\mathrm{%s_{Im}\/ offset \//\/ %s}$" % (colIm, xunit))
        ax.set_xlabel("$\mathrm{%s_{Re} \//\/ %s}$" % (colRe, xunit))
        #        plt.show()
        #    *(ht/wd)
        #    fig_path = EIS_dest_dir.with_suffix('.png')
        #    plt.show()
        #%%
        plt.savefig(DestFile, bbox_inches="tight", dpi=200)
        plt.close()

    def EIS_plotting_combined(PlotCombinedEIS):
        #%% ===== MAKE SPECIAL STACKED PLOTS ======
        # PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        # postOVVout = PostEC.LoadPostOVV()
        PostDestDir = ""
        postOVVout = None
        # FIX ME
        postEIScom = postOVVout.loc[
            postOVVout["Type_Exp"] == "EIS_Combined"
        ].drop_duplicates()
        SampleCodes = "find_codes"

        #    SampleCodes  = pd.read_excel(PostDestDir.joinpath('SampleCodeLst.xlsx'))
        PostDestDirEIScom = PostDestDir.joinpath("EIS_Combined")
        for Elec, ElecGr in postEIScom.groupby(by="Electrolyte"):
            PDDirEIScom = PostDestDirEIScom.joinpath(Elec)
            for Sample, sGr in ElecGr.groupby(by="SampleID"):
                for sf, fGr in sGr.groupby(by="SourceFilename"):
                    AllData_E_file = pd.read_excel(sf)
                    EISovv = fGr
                    DestFile = PDDirEIScom.joinpath(
                        "_".join(
                            [
                                fGr.SampleID.values[0],
                                fGr.Gas.values[0],
                                fGr.Status.values[0],
                                fGr.EXP_date.values[0],
                            ]
                        )
                    ).with_suffix(".png")
                    AllData_E_file.to_excel(DestFile.with_suffix(".xlsx"))
                    SampleCode = SampleCodes.loc[
                        SampleCodes.SampleID == EISovv["SampleID"].unique()[0], :
                    ]
                    for i in ["Y", "Z"]:
                        PlotCombinedEIS(
                            AllData_E_file,
                            EISovv,
                            SampleCode,
                            DestFile.parent.joinpath("{0}_".format(i) + DestFile.stem),
                            i,
                        )

    def plot_par_with_Ev(
        mgrp,
        SampleCodes,
        OriginColor,
        DestFilePars,
        xlabel="E_AppV_RHE",
        force=False,
        ylabels=[],
    ):
        #        OriginColors = Characterization_TypeSetting.OriginColorList()
        EvRHE = "E_AppV_RHE"
        # mgrp.lmfit_var_names.unique()
        plst = list(
            set([a for i in mgrp.lmfit_var_names.unique() for a in i.split(", ")])
        )
        # list(modselect.parameters_guesses.keys())
        plst += [
            i for i in eisplot.extrapars if i in mgrp.columns and mgrp[i].mean() != 0
        ]
        if ylabels:
            plst += ylabels

        # TODO FIX ME
        print("TODO FIX ME pars")
        for par in plst:
            _pf_stem = Path(mgrp.PAR_file.unique()[0]).stem
            xstr = "".join(i for i in xlabel if i.isalnum() or i == "_")
            ystr = "".join(i for i in par if i.isalnum() or i == "_")
            _par_dest = Path(f"{DestFilePars}/{_pf_stem}_{xstr}_{par}.png")
            mgrp = mgrp.dropna(subset=[xlabel, par])
            if not _par_dest.is_file() or force or not mgrp.empty:
                _ylim, logy, logyscale = eisplot(par).pars_plot_lims()
                ylim = (_ylim[0], np.max([mgrp[par].max() * 1.2, _ylim[1]]))
                if par in ylabels:
                    ylim = (mgrp[par].min() * 0.5, mgrp[par].max() * 1.2)
                    logy = False
                ht, wd = 10, 10
                fig, ax = plt.subplots(figsize=(ht, wd))
                ax.set_ylim(ylim)

                ax.set_xlim(0.0 * mgrp[xlabel].min(), mgrp[xlabel].max() * 1.2)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(par)
                if logy:
                    ax.set_yscale(logyscale)

                for gas, ggrp in mgrp.groupby("Gas"):

                    msa = eisplot(gas).gas_plot()
                    msa.update(dict(ms=10))

                    for sID, sgr in ggrp.groupby("SampleID"):
                        sgr = sgr.sort_values(xlabel)
                        try:
                            _ogc = OriginColor.loc[
                                OriginColor.OriginCode
                                == SampleCodes.loc[
                                    SampleCodes.SampleID == sID, "Colorcode"
                                ].unique()[0],
                                "RGB_255",
                            ].iloc[0]
                        except IndexError:
                            _ogc = OriginColor["RGB_255"].iloc[0]

                        ogc = [*_ogc[0:3], 0.6]
                        msa.update(dict(color=ogc, label=f"{gas} {sID}"))
                        if par + "_stderr" in sgr.columns:
                            try:
                                _stderrs = [
                                    i[1] if i[1] / i[0] < 1 else 0
                                    for i in zip(
                                        sgr[par].values, sgr[par + "_stderr"].values
                                    )
                                ]
                            except:
                                _stderrs = []
                        if _stderrs and len(_stderrs) == len(sgr[par].values):
                            ax.errorbar(
                                sgr[xlabel].values,
                                sgr[par].values,
                                yerr=_stderrs,
                                **msa,
                            )
                        else:
                            ax.scatter(
                                sgr[xlabel].values,
                                sgr[par].values,
                                **dict(color=ogc, label=f"{gas} {sID}"),
                            )
                #                               markerstyle=msa[0], alpha=msa[1])
                plt.legend(ncol=2, loc=2)
                plt.savefig(_par_dest, dpi=200, bbox_inches="tight")
                plt.close(fig)
            plt.close()
