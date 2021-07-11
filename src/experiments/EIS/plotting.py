# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:44:59 2020

@author: User
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import cmath
# from cmath import phase
from pathlib import Path
from datetime import datetime
import matplotlib.font_manager

from file_py_helper.find_folders import FindExpFolder

# if __name__ == '__main__':
#    from models import EEC_models_color
# elif any(__name__ in i for i in  ['eis_run_ovv','fitting']):
#    from .models import EEC_models_color
# else:
#    print(__name__)
#    from models import EEC_models_color


globals()["EvRHE"] = "E_AppV_RHE"


from matplotlib.ticker import ScalarFormatter

plt.rc("text", usetex=False)
plt.rc("font", family="serif", size=14)
#%%
class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of
    magnitude"""

    def __init__(self, order_of_mag=0, useOffset=True, useMathText=True):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)

    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag


def plot_nyquist(ax, freq, Z, scale=1, units="Ohms", fmt=".-"):
    """Convenience function for plotting nyquist plots
    Parameters
    ----------
    ax: matplotlib.axes.Axes
        axes on which to plot the nyquist plot
    freq: np.array of floats
        frequencies
    Z: np.array of complex numbers
        impedance data
    scale: float
        the scale for the axes
    units: string
        units for :math:`Z(\\omega)`
    fmt: string
        format string passed to matplotlib (e.g. '.-' or 'o')
    Returns
    -------
    ax: matplotlib.axes.Axes
    """

    ax.plot(np.real(Z), -np.imag(Z), fmt, lw=3)

    # Make the axes square
    ax.axis("square")
    #    print(units)
    if units in ["Ohms", "\Omega"]:
        # Set the labels to -imaginary vs real Impedance
        xlbl = r"$Z^{\prime}(\omega)$ " + "$[{}]$".format(units)
        ylbl = r"$-Z^{\prime\prime}(\omega)$ " + "$[{}]$".format(units)

    elif "S" in units:
        # Set the labels to -imaginary vs real Admittance
        xlbl = r"$Y^{\prime}(\omega)$ " + "$[{}]$".format(units)
        ylbl = r"$-Y^{\prime\prime}(\omega)$ " + "$[{}]$".format(units)

    ax.set_xlabel(xlbl, fontsize=20)
    ax.set_ylabel(ylbl, fontsize=20)

    # Make the tick labels larger
    ax.tick_params(axis="both", which="major", labelsize=14)

    # Change the number of labels on each axis to five
    ax.locator_params(axis="x", nbins=5, tight=True)
    ax.locator_params(axis="y", nbins=5, tight=True)

    # Add a light grid
    ax.grid(b=True, which="major", axis="both", alpha=0.5)

    # Change axis units to 10**log10(scale) and resize the offset text
    ax.xaxis.set_major_formatter(FixedOrderFormatter(-np.log10(scale)))
    ax.yaxis.set_major_formatter(FixedOrderFormatter(-np.log10(scale)))
    y_offset = ax.yaxis.get_offset_text()
    y_offset.set_size(18)
    t = ax.xaxis.get_offset_text()
    t.set_size(18)
    return ax


def plot_linKK(EIS_data_KKvalid, EIS_data_raw, *args, **kwargs):
    #    kwargs = linKK_pars
    global EvRHE
    freq, Zdata, Z_linKK, res_real, res_imag = (
        EIS_data_KKvalid["Frequency(Hz)"].values,
        EIS_data_KKvalid.DATA_Z.values,
        EIS_data_KKvalid.linKK_Z.values,
        EIS_data_KKvalid.linKK_resRe.values,
        EIS_data_KKvalid.linKK_resIm.values,
    )
    #    EIS_invalid  = EIS_data.query('Valid == False')
    EIS_invalid = EIS_data_raw.loc[~EIS_data_raw["Frequency(Hz)"].isin(freq)]
    #    inv_freq,inv_Zdata, inv_Z_linKK, inv_res_real, res_imag =  EIS_invalid['Frequency(Hz)'].values,EIS_invalid.DATA_Z.values,EIS_invalid.linKK_Z.values,EIS_invalid.linKK_resRe.values,EIS_invalid.linKK_resIm.values
    linkKK_invalid_prefit = kwargs.get("linkKK_invalid_prefit", pd.DataFrame())
    res_scale = kwargs.get("res_scale", 1)
    res_ylim = kwargs.get("res_ylim", 0.25)

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 6)
    ax1 = fig.add_subplot(gs[:2, 0:2])
    ax1_zoom = fig.add_subplot(gs[:2, 2:4])
    axY = fig.add_subplot(gs[:2, 4:6])
    ax2 = fig.add_subplot(gs[2, :])
    #    for key, value in kwargs.items():
    #        print("{0}".format(key))
    sup1, sup2 = [], []
    E_Vrhe = kwargs[EvRHE]
    E_mOCP = float(kwargs["Measured_OCP"]) * 1e-3 + float(kwargs["RHE_OCP"])
    basenm = kwargs["basename"]
    sup1 = "{:s} at {:.2f} Vrhe with OCP {:.2f}".format(basenm, E_Vrhe, float(E_mOCP))
    if "linKK_M" in kwargs.keys():
        sup2 = f'Number of RC elements {kwargs["linKK_M"]} with c = {kwargs["linKK_mu"]:.3f}'
    if sup1 and sup2:
        time_now = datetime.now().strftime(format="%c")
        plt.suptitle(f"{sup1}\n{sup2} {time_now}")
    if "type_res" in kwargs.keys():
        set_res_type = kwargs["type_res"]
    else:
        set_res_type = "Z"
    #    print(kwargs.keys())
    # plot original data
    # plot_nyquist(ax1, freq, Zdata, fmt='s',label='data')
    ax1.scatter(
        Zdata.real, -1 * Zdata.imag, c="blue", marker="s", alpha=1, label="data"
    )
    # plot measurement model
    ax1.plot(
        Z_linKK.real,
        -1 * Z_linKK.imag,
        c="orangered",
        ls="-",
        alpha=1,
        label="model Lin-KK",
    )
    # plot_nyquist(ax1, freq, Z_linKK, fmt='-', scale=1, units='\Omega',label=)

    ax1.plot(
        Zdata.real,
        Zdata.real - Zdata.real.min(),
        c="grey",
        ls="--",
        alpha=0.5,
        label="45 angle",
    )
    ax1.scatter(
        EIS_invalid.DATA_Zre,
        -1 * EIS_invalid.DATA_Zim,
        c="grey",
        s=40,
        alpha=0.6,
        marker="x",
        label="removed points",
    )
    # _ax1legend = ['data', 'model Lin-KK','45 angle line','removed points']

    # ax1.plot(Zdata.real+4*(Zdata.real.min()),Zdata.real-Zdata.real.min(), c='grey', ls='--', alpha=0.5)

    # _x = np.linspace(0,Zdata.real.max()*1.5)
    # _y = [1*(i-(Zdata.real[4]-np.abs(Zdata.imag[4]))) for i in _x]
    # ax1.plot(_x, _y, c='grey', ls='--', alpha=0.5)
    # Zdata.real[-6:]-Zdata.imag[-6:]
    # plot Z zoom on HF part
    _lincols = [i for i in EIS_data_KKvalid.columns if "_lintangent" in i]
    _lincols_freq = set([i[0:-3] for i in _lincols])
    if _lincols:
        for _fnm in list(_lincols_freq):
            _re0 = EIS_data_KKvalid[f"{_fnm}_re"].min()
            _lbl = f"45 lintang {_fnm[-2:]}: {_re0:.1f}"
            EIS_data_KKvalid.plot(
                x=f"{_fnm}_re", y=f"{_fnm}_im", ax=ax1, ls="--", label=_lbl
            )
            EIS_data_KKvalid.plot(
                x=f"{_fnm}_re", y=f"{_fnm}_im", ax=ax1_zoom, ls="--", label=_lbl
            )
            _linY = (
                EIS_data_KKvalid[f"{_fnm}_re"] + 1j * EIS_data_KKvalid[f"{_fnm}_im"]
            ) ** -1
            axY.plot(
                np.real(_linY.values), -1 * np.imag(_linY.values), ls="--", label=_lbl
            )
            # EIS_data_KKvalid.plot(x = f'{_fnm}_re', y= f'{_fnm}_im',ax=ax1_zoom,

            # _ax1legend.append(_lbl)
    ax1.legend(loc="best", fontsize=12)
    _ZmaxScale = np.max([Zdata.real.max() + 10, np.max(-1 * Zdata.imag) + 10])
    ax1.set_xlim(0, _ZmaxScale)
    ax1.set_ylim(0, _ZmaxScale)

    Zmax = EIS_data_KKvalid[EIS_data_KKvalid["Frequency(Hz)"] > 1e1].DATA_Zre.max()

    ax1_zoom.scatter(
        Zdata.real, -1 * Zdata.imag, c="blue", marker="s", alpha=1, label="data"
    )
    # plot measurement model
    ax1_zoom.plot(
        Z_linKK.real,
        -1 * Z_linKK.imag,
        c="orangered",
        ls="-",
        alpha=1,
        label="model Lin-KK",
    )
    ax1_zoom.scatter(
        EIS_invalid.DATA_Zre,
        -1 * EIS_invalid.DATA_Zim,
        c="grey",
        s=40,
        alpha=0.6,
        marker="x",
    )
    # plot_nyquist(ax1_zoom, freq, Zdata, fmt='s')
    ax1_zoom.plot(
        Zdata.real - 2,
        Zdata.real - Zdata.real.min(),
        c="grey",
        ls="--",
        alpha=0.5,
        label="45 angle",
    )
    # plot measurement model
    # plot_nyquist(ax1_zoom, freq, Z_linKK, fmt='-', scale=1, units='\Omega')
    # ,c='orangered')

    ax1_zoom.set_xlim(Zdata.real.min() / 2, Zmax)
    ax1_zoom.set_ylim(0, Zmax)

    # xlbl = (r'$Z^{\prime}(\omega)$ ' + '$[{}]$'.format( ['Ohms','\Omega']))
    # ylbl = (r'$-Z^{\prime\prime}(\omega)$ ' + '$[{}]$'.format( ['Ohms','\Omega']))
    ax1_zoom.set_xlabel("Z_re / Ohm")
    ax1_zoom.set_ylabel("Z_im / Ohm")
    ax1.set_xlabel("Z_re / Ohm")
    ax1.set_ylabel("Z_im / Ohm")
    # plot original data
    # === Admittance ===
    axY.scatter(
        (Zdata ** -1).real,
        (Zdata ** -1).imag,
        c="blue",
        marker="s",
        alpha=1,
        label="data",
    )
    # plot measurement model
    axY.plot(
        (Z_linKK ** -1).real,
        (Z_linKK ** -1).imag,
        c="orangered",
        ls="-",
        alpha=1,
        label="model Lin-KK",
    )
    # plot_nyquist(axY, freq,  (Zdata**-1).real+1j*-1*(Zdata**-1).imag, fmt='s')
    # plot_nyquist(axY, freq, (Z_linKK**-1).real+1j*-1*(Z_linKK**-1).imag, fmt='-', scale=1, units='S')
    # ,c='orangered')
    axY.scatter(
        EIS_invalid.DATA_Yre,
        EIS_invalid.DATA_Yim,
        c="grey",
        s=40,
        alpha=0.6,
        marker="x",
    )

    Yxlbl = "Y_re / S"  # + '$[{}]$'.format('S'))
    Yylbl = "Y_im / S"  # (r'$-Y^{\prime\prime}(\omega)$ ' + '$[{}]$'.format('S'))
    axY.set_xlabel(Yxlbl)
    axY.set_ylabel(Yylbl)
    _Ymax = (Zdata ** -1).real.max()
    axY.set_xlim(0, _Ymax)
    axY.set_ylim(0, _Ymax)

    # Plot residuals
    res_Real_c, res_Imag_c = "orangered", "limegreen"
    ax2.scatter(
        freq, res_real * res_scale, label=r"$\Delta_{\mathrm{Real}}$", c=res_Real_c
    )
    ax2.scatter(
        freq, res_imag * res_scale, label=r"$\Delta_{\mathrm{Imag}}$", c=res_Imag_c
    )
    res_ylim_max = np.max([np.abs(res_real).max(), np.abs(res_real).max()])
    if not linkKK_invalid_prefit.empty:
        ax2.scatter(
            linkKK_invalid_prefit["Frequency(Hz)"],
            linkKK_invalid_prefit.linKK_resRe_prefit * res_scale,
            c=res_Real_c,
            s=40,
            alpha=0.6,
            marker="x",
        )
        ax2.scatter(
            linkKK_invalid_prefit["Frequency(Hz)"],
            linkKK_invalid_prefit.linKK_resIm_prefit * res_scale,
            c=res_Imag_c,
            s=40,
            alpha=0.6,
            marker="x",
        )
        invalid_max = (
            np.abs(linkKK_invalid_prefit[["linKK_resRe_prefit", "linKK_resIm_prefit"]])
            .max()
            .max()
        )
        res_ylim_max = np.max([res_ylim_max, invalid_max])
    #    ax2.scatter(freq, EIS_data_KKvalid.linKK_resRe, c='grey', s=40, alpha=0.6, marker='x')
    if "linKK_limit_Re" in kwargs.keys():
        ax2.plot(
            freq,
            len(freq) * [kwargs["linKK_limit_Re"] * res_scale],
            "--",
            alpha=0.7,
            c=res_Real_c,
        )
        ax2.plot(
            freq,
            len(freq) * [kwargs["linKK_limit_Re"] * -1 * res_scale],
            "--",
            alpha=0.7,
            c=res_Real_c,
        )

        ax2.plot(
            freq,
            len(freq) * [kwargs["linKK_limit_Im"] * res_scale],
            "--",
            alpha=0.7,
            c=res_Imag_c,
        )
        ax2.plot(
            freq,
            len(freq) * [kwargs["linKK_limit_Im"] * -1 * res_scale],
            "--",
            alpha=0.7,
            c=res_Imag_c,
        )
    #        ax2.plot(freq, kwargs['KK_valid_limit'][1]*res_scale, '--', label=r'$\Delta_{\mathrm{Real}}$',alpha=0.7)
    #        ax2.plot(freq, kwargs['KK_valid_limit'][1]*-1*res_scale, '--', label=r'$\Delta_{\mathrm{Real}}$', alpha=0.7)
    ax2.set_title(
        f'Lin-KK Model Error for {set_res_type} (factor: {kwargs["linKK_trimming_factor"]})',
        fontsize=14,
    )

    ax2.tick_params(axis="both", which="major", labelsize=12)
    ax2.set_ylabel("$\Delta$ $(\%)$", fontsize=14)
    ax2.set_xlabel("$f$ [Hz]", fontsize=14)
    ax2.set_xscale("log")

    if res_ylim_max < 0.05:
        res_ylim = 0.05
    elif res_ylim_max < 0.1:
        res_ylim = 0.1
    elif res_ylim_max < 0.2:
        res_ylim = 0.2

    ax2.set_ylim(-res_ylim * res_scale, res_ylim * res_scale)
    ax2.legend(loc="upper left", fontsize=14, ncol=2)
    #    vals = ax2.get_yticks()
    #    ax2.set_yticklabels(['{:.0%}'.format(x) for x in vals])
    plt.tight_layout()
    #    plt.show()
    if "save_target" in kwargs.keys():
        plt.savefig(Path(kwargs["save_target"]), dpi=100, bbox_inches="tight")
    plt.close()


#
# def old_EEC_models_color(modname=''):
#    model_color = {'Model(Singh2015_RQR)' : 'cyan', 'Model(Singh2015_RQRQ)' : 'darkcyan',
#      'Model(Singh2015_RQRQR)' : 'orangered', 'Model(Bandarenka_2011_RQRQR)' : 'gray',
#      'Model(Singh2015_RQRWR)' : 'purple', 'Model(Randles_RQRQ)' : 'gold',
#      'Model(Singh2015_R3RQ)' : 'green'}


def EEC_models_color(modname=""):
    model_color = {
        "Model(EEC_Randles_RWpCPE)": "gray",
        "Model(EEC_2CPE)": "cyan",
        "Model(EEC_2CPEpW)": "orangered",
        "Model(EEC_RQ_RC_RW)": "red",
        "Model(EEC_RQ_RC_RQ)": "fuchsia",
        "Model(Singh2015_RQRWR)": "purple",
        "Model(Randles_RQRQ)": "gold",
        "Model(EEC_Randles_RWpCPE_CPE)": "red",
        "Model(EEC_2CPE_W)": "green",
        "Model(EEC_2CPEpRW)": "gold",
    }
    if modname:
        return model_color.get(modname, "fuchsia")
    else:
        return model_color


def EIS_plotting_per_EV(
    EIS_fit_data,
    pars_models,
    EIS_outPath_target_png,
    plot_show=False,
    std_print_model="RL-TLM(Rct-Qad-W)",
):
    """This method plots each EIS fit spectrum from the LMfit_method output"""
    #%%
    # EIS_fit_data, pars_models,EIS_outPath_target_png =  self.EIS_fit_data,self.pars_models.query('FINAL_FIT == 1'), spectra_fit_outpath_png
    # plot_show = False
    # std_print_model = 'M_RC_CPE_W'

    global EvRHE
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(22, 22))
    ax, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    ax5, ax6 = axes[2, 0], axes[2, 1]

    sID = pars_models["SampleID"].iloc[0]
    # EIS_fit_data.index.get_level_values('SampleID').unique()[0]
    Elec, pH = pars_models["Electrolyte"].iloc[0], pars_models["pH"].iloc[0]
    # EIS_fit_data.index.get_level_values('Electrolyte').unique()[0], EIS_fit_data.index.get_level_values('pH').unique()[0]
    EappV = pars_models[EvRHE].iloc[0]
    # EIS_fit_data.index.get_level_values(EvRHE).unique()[0]
    time_now = datetime.now().strftime(format="%c")
    _sID_stem = f"{sID}, {EIS_outPath_target_png.stem}:" + "\n"
    _pH_elec = (
        "$\mathrm{E_{DC}\/ =\/}$"
        + f"{EappV:.2f}"
        + "$\mathrm{\/ V_{RHE} \/}$"
        + f"in {Elec} at pH = {pH:.0f}, {time_now}"
    )
    # _weight

    fig.suptitle(_sID_stem + _pH_elec)

    #    fig.suptitle(f'{sID}, {EIS_outPath_target_png.stem}: \n$\mathrm{E_{DC}\/ =\/ %.2f\/ V_{RHE} \/ in \/%s\/at\/pH\/=\/%.0f}$\n%s' %(sID, , EappV,Elec,pH,time_now))

    dataC, fitC, extraC, initC = "gray", "tab:red", "gold", "gray"
    #        outData.plot(x='FIT_Zre',y='FIT_-Zim',kind='line',c='orange',ax=ax)
    """axes: Data plots """
    #    EIS_fit_data.plot(x='DATA_Zre',y='DATA_-Zim',kind='scatter',c=dataC,s=100,ax=ax)
    ax.scatter(
        EIS_fit_data["DATA_Zre"], EIS_fit_data["DATA_-Zim"], c=dataC, s=100, alpha=0.5
    )
    ax.set(xlim=(0, ax.axis("equal")[1]), ylim=(0, ax.axis("equal")[-1]))
    # add 45 angle line TODO
    # ax.plot(EIS_fit_data['DATA_Zre']+EIS_fit_data['DATA_Zre'].min(),EIS_fit_data['DATA_Zre'], c='grey',alpha=0.5,ls='--')

    ax.grid(True)
    EIS_fit_data.plot(
        x="DATA_Yre", y="DATA_Yim", kind="scatter", c=dataC, s=100, ax=ax2, alpha=0.5
    )
    ax2.set(xlim=(0, ax2.axis("equal")[1]), ylim=(0, ax2.axis("equal")[-1]))
    ax2.grid(True)
    EIS_fit_data.plot(
        x="Frequency(Hz)",
        y="DATA_Zmod",
        kind="scatter",
        c=dataC,
        s=100,
        ax=ax3,
        logy=False,
        logx=True,
        alpha=0.5,
    )
    ax3.axis("tight")
    EIS_fit_data.plot(
        x="Frequency(Hz)",
        y="DATA_-Zangle",
        kind="scatter",
        c=dataC,
        s=100,
        ax=ax4,
        logy=False,
        logx=True,
        alpha=0.5,
    )
    """axes plot Fittings per model"""
    # pars_models.Model_EEC.nunique() == len(pars_models)
    pars_models = pars_models.sort_values("lmfit_aic", ascending=True)
    EIS_fit_data_grpmod = EIS_fit_data.groupby(["Model_EEC"])
    _ncount = 0
    _nmodels = 5
    len(pars_models)
    for n, r in pars_models.iterrows():

        modnm, _parsmodgrp = r["Model_EEC"], r
        # groupby(['Model_EEC']):
        """ax: Impedance (Nyquist) Plots"""
        modgr = EIS_fit_data_grpmod.get_group(modnm)
        # _parsmodgrp = pars_models.loc[pars_models.Model_EEC == modnm]
        mod_color = _parsmodgrp["plot_color"]
        _aic = _parsmodgrp["lmfit_aic"]
        mc = [float(i) for i in mod_color.split(", ")]

        _res_alpha = 1 - (1 / _nmodels) * _ncount
        _res_alpha = _res_alpha if _res_alpha >= 0.3 else 0.3
        _ncount += 1

        #        EEC_models_color(modnm)
        #     mod,mod2 = Model(EEC_ORR),Model(EEC_ORRpOx)
        #        EIS_fit_data.plot(x='FIT2_Zre',y='FIT2_Zim',kind='line',c=fitC,lw=2.5,ax=ax,label='FIT_2 (=EEC_ORRpOx)')
        #        modgr.plot(x='FIT_Zre',y='FIT_-Zim',kind='line',lw=2.5,ax=ax,label='FIT: {0}'.format(modnm))
        ax.plot(
            modgr["FIT_Zre"],
            modgr["FIT_-Zim"],
            lw=2.5,
            alpha=_res_alpha,
            label="FIT: {0}".format(modnm),
            c=mc,
        )

        #        ax.plot(EIS_fit_data['FIT2_Zre'].values,EIS_fit_data['FIT2_Zre'].values,ls=':',c='k',lw=1)
        #    outData.plot(x='FIT1_Zre',y='FIT1_Zre',kind='line',ls='.',c='k',lw=10,ax=ax)

        """ax2: Admittance Plots"""
        ax2.plot(modgr["FIT_Yre"], modgr["FIT_Yim"], lw=2.5, alpha=_res_alpha, c=mc)

        #        modgr.plot(x='FIT_Yre',y='FIT_Yim',kind='line',lw=2.5,ax=ax2,alpha=0.7)
        #        ax2.legend(ncol=2,fontsize=12,loc='upper right')
        #        outData.plot(x='FIT1_Yre',y='FIT1_Yim',kind='line',c=extraC,lw=2.5,ax=ax2)
        #        outData.plot(x='INIT2_Yre',y='INIT2_Yim',kind='line',c=initC,lw=2.5,ax=ax2,alpha=0.5)
        #        outData.plot(x='INIT1_Yre',y='INIT1_Yim',kind='line',c=initC,lw=2.5,ax=ax2,alpha=0.5)

        """ax3 & ax4: Frequency (Bode) Plots"""
        #        outData.plot(x='Frequency(Hz)',y='FIT2_Zmod',kind='line',c=fitC,lw=2.5,ax=ax3,logy=True,logx=True)
        modgr.plot(
            x="Frequency(Hz)",
            y="FIT_Zmod",
            kind="line",
            lw=2.5,
            alpha=_res_alpha,
            ax=ax3,
            logy=False,
            logx=True,
            c=mc,
        )
        ax3.set_xlim(0.5, 3e4)

        modgr.plot(
            x="Frequency(Hz)",
            y="FIT_-Zangle",
            kind="line",
            lw=2.5,
            alpha=_res_alpha,
            ax=ax4,
            logy=False,
            logx=True,
            c=mc,
        )
        ax4.set_xlim(0.5, 3e4)
        ax4.set_ylim(-40, 65)

        #    outData.plot(x='Frequency(Hz)',y='INIT2_Zphase',kind='line',c=initC,lw=1.5,ax=ax4,logy=False,logx=True)
        #        outData.plot(x='Frequency(Hz)',y='FIT1_Zphase',kind='line',c=extraC,lw=2.5,ax=ax4,logy=False,logx=True)

        """ax5 & ax5: Frequency Residual Plots"""

        modgr.plot(
            x="Frequency(Hz)",
            y="errRe",
            kind="line",
            ax=ax5,
            logy=False,
            logx=True,
            alpha=_res_alpha,
            label=f"errRe: {modnm}, {_aic:.2f}",
            c=mc,
        )

        ax5.set_ylim(-0.15, 0.15)
        ax6.set_ylabel("$\Delta$ Re $(\%)$")
        modgr.plot(
            x="Frequency(Hz)",
            y="errIm",
            kind="line",
            ax=ax6,
            logy=False,
            logx=True,
            alpha=_res_alpha,
            c=mc,
        )

        ax6.set_ylim(-0.15, 0.15)
        ax6.set_ylabel("$\Delta$ Im $(\%)$")

    ax.legend(ncol=3, fontsize=12, loc="upper left", bbox_to_anchor=(-0.05, 1.31))
    ax3.legend([])
    ax4.legend([])
    ax5.legend(ncol=1, fontsize=12)
    ax6.legend([])

    #        outData.plot(x='Frequency(Hz)',y='errRe1',kind='line',c=fitC,ax=ax5,logy=False,logx=True)
    #        outData.plot(x='Frequency(Hz)',y='errIm1',kind='line',c=extraC,ax=ax6,logy=False,logx=True)
    #    fig_path = EIS_dest_dir.joinpath(Path(str(Path(EISgr['File'].unique()[0]).stem)+'_%.0fmV.png' %(E_dc*1000)))
    #    if lmfit_out_lst:
    pars_std_print_model = pars_models.loc[
        pars_models.Model_EEC.str.contains(std_print_model)
    ]
    if not pars_std_print_model.empty:
        out1 = pars_std_print_model["lmfit_out"].to_list()[0]
        _second_pars_print_model = pars_models.loc[
            ~pars_models.index.isin(pars_std_print_model.index)
        ]
        if not _second_pars_print_model.empty:
            out2 = _second_pars_print_model.head(1)["lmfit_out"].to_list()[0]
        else:
            out2 = out1
    else:
        if len(pars_models.head(2)) == 2:
            out1, out2 = pars_models.head(2)["lmfit_out"].to_list()
        else:
            out1, out2 = (
                pars_models.iloc[0]["lmfit_out"],
                pars_models.iloc[0]["lmfit_out"],
            )
    #    out1,ou2 = lmfit_out_lst[best_2mods[0]]
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax4.text(
        2.25,
        1.3,
        out1.fit_report_wts,
        transform=ax.transAxes,
        fontsize=14,
        verticalalignment="top",
        bbox=props,
    )
    ax6.text(
        2.25,
        -0.6,
        out2.fit_report_wts,
        transform=ax.transAxes,
        fontsize=14,
        verticalalignment="top",
        bbox=props,
    )
    plt.savefig(EIS_outPath_target_png, bbox_inches="tight", dpi=100)
    plt.close()


#    else:
#        if not pars_models.sort_values('RedChisqr',ascending=True).head(1).empty:
#            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#            out1 = pars_models.sort_values('RedChisqr',ascending=True).head(1)
#            ax4.text(2.25, 1, out1.fit_report(min_correl=0.99), transform=ax.transAxes, fontsize=15,
#                 verticalalignment='top', bbox=props)
#%%
#    if plot_show:
#        plt.show()
#    xl_path = EIS_dest_dir.joinpath(Path(str(Path(EISgr_data_EV['File'].unique()[0]).stem)+'_%.0fmV.xlsx' %(E_dc*1000)))
#    out_cols = ['Frequency(Hz)','FIT2_Zre','FIT2_Zim','DATA_Zre','DATA_-Zim','FIT2_Yre','FIT2_Yim','DATA_Yre','DATA_Yim']
#    fig_path.parent.mkdir(parents=True,exist_ok=True)
#    outData[out_cols].to_excel(xl_path)
#    print('EIS fig saved: %s'  %fig_path)
#    plt.show()
#%%
def EIS_Trimming_plot(
    EISgr_data_EV, EISovv, EIS_dest_dir, outData, TrimmedOutData, E_dc
):
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 14))
    ax, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    time_now = datetime.now().strftime(format="%c")
    fig.suptitle(
        "%s, %s: \n$\mathrm{E_{DC}\/ =\/ %.2f\/ V_{RHE} \/ in \/%s\/at\/pH\/=\/%.0f}$ %s"
        % (
            EISgr_data_EV["SampleID"].unique()[0],
            Path(EISgr_data_EV["PAR_file"].unique()[0]).stem,
            E_dc,
            EISovv["Electrolyte"].unique()[0],
            EISovv["pH"].unique()[::],
            time_now,
        )
    )
    dataC, fitC, extraC, initC = "tab:blue", "tab:red", "gold", "gray"
    #        outData.plot(x='FIT_Zre',y='FIT_-Zim',kind='line',c='orange',ax=ax)
    """ax: Impedance (Nyquist) Plots"""
    #     mod,mod2 = Model(EEC_ORR),Model(EEC_ORRpOx)
    #    outData.plot(x='FIT2_Zre',y='FIT2_Zim',kind='line',c=fitC,lw=2.5,ax=ax,label='raw')
    #    TrimmedOutData.plot(x='FIT2_Zre',y='FIT2_Zim',kind='line',c=fitC,lw=2.5,ax=ax,label='trimmed')
    #    outData.plot(x='FIT1_Zre',y='FIT1_-Zim',kind='line',c=extraC,lw=2.5,ax=ax,label='FIT_1 (=EEC_ORR)')
    outData.plot(
        x="DATA_Zre",
        y="DATA_-Zim",
        kind="scatter",
        c=dataC,
        s=150,
        ax=ax,
        alpha=0.7,
        label="raw",
    )
    TrimmedOutData.plot(
        x="DATA_Zre",
        y="DATA_-Zim",
        kind="scatter",
        c=fitC,
        s=100,
        ax=ax,
        marker="^",
        label="trimmed",
    )
    #    ax.plot(outData['FIT2_Zre'].values,outData['FIT2_Zre'].values,ls=':',c='k',lw=1)
    #    outData.plot(x='FIT1_Zre',y='FIT1_Zre',kind='line',ls='.',c='k',lw=10,ax=ax)
    ax.axis("square")
    """ax2: Admittance Plots"""
    #    outData.plot(x='FIT2_Yre',y='FIT2_Yim',kind='line',c=fitC,lw=2.5,ax=ax2)
    outData.plot(
        x="DATA_Yre", y="DATA_Yim", kind="scatter", c=dataC, s=150, ax=ax2, alpha=0.7
    )
    TrimmedOutData.plot(
        x="DATA_Yre", y="DATA_Yim", kind="scatter", c=fitC, s=100, ax=ax2, marker="^"
    )
    #    outData.plot(x='FIT1_Yre',y='FIT1_Yim',kind='line',c=extraC,lw=2.5,ax=ax2)
    #    outData.plot(x='INIT2_Yre',y='INIT2_Yim',kind='line',c=initC,lw=2.5,ax=ax2,alpha=0.5)
    #    outData.plot(x='INIT1_Yre',y='INIT1_Yim',kind='line',c=initC,lw=2.5,ax=ax2,alpha=0.5)
    ax2.axis("square")
    """ax3 & ax4: Frequency (Bode) Plots"""
    #    outData.plot(x='Frequency(Hz)',y='FIT2_Zmod',kind='line',c=fitC,lw=2.5,ax=ax3,logy=True,logx=True)
    #    outData.plot(x='Frequency(Hz)',y='FIT1_Zmod',kind='line',c=extraC,lw=2.5,ax=ax3,logy=True,logx=True)
    ax3.axis("auto")
    outData.plot(
        x="Frequency(Hz)",
        y="DATA_Zmod",
        kind="scatter",
        c=dataC,
        s=150,
        ax=ax3,
        logx=True,
        alpha=0.7,
    )
    TrimmedOutData.plot(
        x="Frequency(Hz)",
        y="DATA_Zmod",
        kind="scatter",
        c=fitC,
        s=100,
        ax=ax3,
        logy=True,
        logx=True,
        marker="^",
    )

    #    outData.plot(x='Frequency(Hz)',y='FIT2_Zphase',kind='line',c=fitC,lw=2.5,ax=ax4,logy=False,logx=True)
    #    outData.plot(x='Frequency(Hz)',y='INIT2_Zphase',kind='line',c=initC,lw=1.5,ax=ax4,logy=False,logx=True)
    #    outData.plot(x='Frequency(Hz)',y='FIT1_Zphase',kind='line',c=extraC,lw=2.5,ax=ax4,logy=False,logx=True)
    ax4.axis("normal")
    outData.plot(
        x="Frequency(Hz)",
        y="DATA_Zim",
        kind="scatter",
        c=dataC,
        s=150,
        ax=ax4,
        logy=False,
        logx=True,
        alpha=0.7,
    )
    TrimmedOutData.plot(
        x="Frequency(Hz)",
        y="DATA_Zim",
        kind="scatter",
        c=fitC,
        s=100,
        ax=ax4,
        logy=False,
        logx=True,
        marker="^",
    )

    plt.show()


#    AllData_E_file,EISovv, Pars, E_data_combined_path_target = spectra_combined,EVgr, EVgr, PDD_rpm.joinpath('EIS_combined_{0}.jpg'.format('_'.join([str(i) for i in EV])))
# AllData_E_file,gr_EIS_ovv,Parsout,E_data_combined_path_target333
def EIS_plotting_EvRHE(AllData_E_file, EISovv, Pars, E_data_combined_path_target):

    #%% ===== MAKE SPECIAL STACKED PLOTS ======
    #    maxLim = (AllData_E_file[['DATA_Yim','DATA_Yre']].max()).max()
    # TODO: add impedance in plots
    global EvRHE
    maxYim = AllData_E_file["DATA_Yim"].max()
    maxYre = AllData_E_file["DATA_Yre"].max()

    maxZim = AllData_E_file["DATA_Yim"].max()
    maxZre = AllData_E_file["DATA_Yre"].max()

    combined_grouper = EvRHE
    AllData_E_file_undup = AllData_E_file.loc[:, ~AllData_E_file.columns.duplicated()]
    Lenrows = AllData_E_file_undup[combined_grouper].nunique()
    if Lenrows == 1 and AllData_E_file_undup["RPM_DAC"].nunique() > 1:
        combined_grouper = "RPM_DAC"
        Lenrows = AllData_E_file_undup[combined_grouper].nunique()
        if Lenrows == 1:
            combined_grouper = "Segment"
            Lenrows = AllData_E_file_undup[combined_grouper].nunique()

    #    fig,axes = plt.subplots(nrows=Lenrows ,sharex=True,sharey=True)
    ht, wd = 20, 15
    #    fig,ax = plt.subplots(figsize=(ht,wd))
    fig = plt.figure(constrained_layout=True, figsize=(ht, wd))
    gs = plt.GridSpec(4, 4, figure=fig)
    ax = fig.add_subplot(gs[:, :2])
    ax2 = fig.add_subplot(gs[0, 2:])
    ax2_twin = ax2.twinx()
    ax3 = fig.add_subplot(gs[1:2, 2:])
    ax3_twin = ax3.twinx()
    ax4 = fig.add_subplot(gs[2:3, 2:])
    ax4_twin = ax4.twinx()
    ax5 = fig.add_subplot(gs[3:, 2:])
    ax5_twin = ax5.twinx()

    fig.suptitle(
        "%s, %s: \n in  %s at pH = %.0f"
        % (
            EISovv["SampleID"].unique()[0],
            Path(EISovv["basename"].unique()[0]).stem,
            EISovv["Electrolyte"].unique()[0],
            EISovv["pH"].unique()[::],
        )
    )
    dataC, fitC, fit1C, extraC, initC = (
        "tab:blue",
        "tab:red",
        "tab:green",
        "gold",
        "gray",
    )

    #    ax.set_xlim(0,maxYre)
    #    ax.set_ylim(0,maxYre+maxYim*Lenrows)
    ax.grid(True)
    ax.axis("equal")
    #    for En,Ev in enumerate(AllData_E_file[EvRHE].unique()):
    #    AllData_E_file[EvRHE].max() - AllData_E_file[EvRHE].min()
    AllData_E_file_undup.sort_values(
        by=[combined_grouper, "Frequency(Hz)"], ascending=True, inplace=True
    )
    Lenrows = AllData_E_file_undup[combined_grouper].nunique()
    Emin, Emax = (
        AllData_E_file_undup[combined_grouper].min(),
        AllData_E_file_undup[combined_grouper].max(),
    )
    E_diff_minmax = Emax - Emin if not Emax == Emin else Emin
    E_lst_bestmod = []
    for Ev, Egr in AllData_E_file_undup.groupby(combined_grouper):

        Enum = Lenrows * (Ev / E_diff_minmax)
        vertical_data_shift = len(Egr) * [maxYim * Enum]
        ax.scatter(
            Egr["DATA_Yre"].values,
            Egr["DATA_Yim"].values + vertical_data_shift,
            c=dataC,
            s=100,
            alpha=0.6,
        )
        #        ax.annotate("$\mathrm{%.2f} \/ V_{RHE} $"%Ev,
        #                    xy=(maxYre*1.1, maxYim*Enum), xycoords='data')
        #        Pars.loc[np.isclose(Pars.loc[Pars[EvRHE] == Ev,'RedChisqr'],Pars.loc[Pars[EvRHE] == Ev].RedChisqr.min(),rtol=1E-2)]
        E_bestmods = (
            Pars.loc[Pars[combined_grouper] == Ev]
            .loc[
                np.isclose(
                    Pars.loc[Pars[combined_grouper] == Ev, "lmfit_redchi"],
                    Pars.loc[Pars[combined_grouper] == Ev].lmfit_redchi.min(),
                    rtol=5e-2,
                ),
                ["lmfit_redchi", "Model_EEC", "lmfit_chiqsr"],
            ]
            .sort_values(by="lmfit_redchi")
        )
        E_lst_bestmod.append((Ev, E_bestmods.head(3).Model_EEC.to_list()))
        #        E_bestmod_redchisq = Pars.loc[Pars.index == Pars.loc[Pars[EvRHE] == Ev].RedChisqr.idxmin(),'Model_EEC'].to_list()[0]
        #        E_bestmod_chisq = Pars.loc[Pars.index == Pars.loc[Pars[EvRHE] == Ev].Chisqr.idxmin(),'Model_EEC'].to_list()[0]
        if not E_bestmods.empty:
            E_bestmod_label = "\n".join(E_bestmods.Model_EEC.to_list())
            plot_fit_models = Egr.loc[
                Egr.Model_EEC.isin(E_bestmods.head(3).Model_EEC.tolist())
            ]
        else:
            E_bestmod_label = "\n".join(list(Pars.Model_EEC.unique()))
            plot_fit_models = Egr
        #        else:
        #            E_bestmod_label = '{}, \n RedChi: {}'.format(E_bestmod_chisq,E_bestmod_redchisq)
        ax.annotate(
            "{:.2f} {}, {} ".format(Ev, combined_grouper, E_bestmod_label),
            xy=(maxYre * 1.1, maxYim * Enum),
            xycoords="data",
        )
        E_count = 0
        for modnm, Emodgr in plot_fit_models.groupby("Model_EEC"):
            vertical_data_shift_mod = len(Emodgr) * [maxYim * Enum]
            ax.plot(
                Emodgr["FIT_Yre"].values,
                Emodgr["FIT_Yim"].values + vertical_data_shift_mod,
                lw=3,
                alpha=0.7,
                label="FIT: {0}".format(modnm),
                c=EEC_models_color(modnm),
            )
    #            ax.plot(Edata['FIT1_Yre'].values,Edata['FIT1_Yim'].values+len(Edata)*[maxYim*En],c=fit1C,ls='dotted',lw=3,alpha=0.7)
    #            if modcount == 0:
    #                ax.annotate("$\mathrm{%.2f} \/ V_{RHE} $"%Ev,
    #                    xy=(maxYre*1.1, maxYim*Enum), xycoords='data')
    #    ax.set_ylim(0,maxYim*(En+1))
    try:
        fast_checking_EEC_models = ["Model(RL-TLM(Rct-Qad-W))"]
        ["EEC_Randles_RWpCPE_CPE", "Model(EEC_2CPEpRWs)", "Model(EEC_2CPE)"]
        bestmodPars = Pars.loc[Pars.Model_EEC.isin(fast_checking_EEC_models)]
        for modnm, Pmodgr in bestmodPars.groupby("Model_EEC"):
            first_sct_title = ""
            if Pmodgr[combined_grouper].nunique() > 1:
                bestmod_xscatter = Pmodgr[combined_grouper].values
                first_sct_title = "{0} rpm".format(Pmodgr[combined_grouper].unique()[0])
            #        elif Pmodgr['RPM_DAC'].nunique() > 2:
            #            bestmod_xscatter = Pmodgr['RPM_DAC'].values
            #            first_sct_title = '{0} Vrhe'.format(Pmodgr['E_RHE'].unique()[0])
            else:
                bestmod_xscatter = Pmodgr[EvRHE].values
                first_sct_title = ""
            mc = EEC_models_color(modnm)
            sc2Rs = ax2.scatter(
                bestmod_xscatter,
                Pmodgr.Rs.values,
                label="Rs " + modnm,
                s=80,
                marker="s",
                alpha=0.7,
                c=mc,
            )
            ax2.set_title(first_sct_title)
            sc2Rct = ax2_twin.scatter(
                bestmod_xscatter,
                Pmodgr.Rct.values,
                label="Rct",
                marker="*",
                s=100,
                c=mc,
            )

            ax3.scatter(
                bestmod_xscatter,
                Pmodgr.Rct.values,
                s=100,
                label="Rct " + modnm,
                marker="*",
                c=mc,
            )
            ax3_twin.scatter(
                bestmod_xscatter,
                Pmodgr.Rorr.values,
                label="Rorr",
                marker="^",
                s=100,
                c=mc,
            )

            ax4.scatter(
                bestmod_xscatter,
                Pmodgr.Cdlp.values,
                s=100,
                marker="o",
                label="Cdpl " + modnm,
                c=mc,
            )
            ax4_twin.scatter(
                bestmod_xscatter,
                Pmodgr.nDL.values,
                marker="X",
                s=100,
                label="nDL",
                c=mc,
            )

            ax5.scatter(
                bestmod_xscatter,
                Pmodgr.Qad.values,
                s=100,
                marker="D",
                label="Qad " + modnm,
                c=mc,
            )
            ax5_twin.scatter(
                bestmod_xscatter,
                Pmodgr.nAd.values,
                marker="P",
                s=100,
                label="nAd",
                c=mc,
            )
        #    a = ax2.get_legend_handles_labels()
        #    b = ax2_twin.get_legend_handles_labels()
        ax2.legend(fontsize=10, loc="upper left")
        ax2_twin.legend(fontsize=10, loc="upper right")
        ax2.set_ylim(0, 60)
        ax2.set_ylabel("Rs")
        ax2_twin.set_yscale("log")
        ax2_twin.set_ylim(1, 5e3)
        ax2_twin.set_ylabel("R_ct")
        ax3.legend(fontsize=10)
        ax3.set_yscale("log")
        ax3.set_ylabel("R_ct")
        ax3_twin.set_yscale("log")
        ax3_twin.set_ylabel("R_orr")
        ax3.legend(fontsize=10, loc="upper left")
        ax3_twin.legend(fontsize=10, loc="upper right")
        ax4_twin.set_ylim(0.3, 1.2)
        ax4.set_ylabel("Cdlp")
        ax4.legend(fontsize=10, loc="upper left")
        ax4_twin.legend(fontsize=10, loc="upper right")
        if bestmodPars.Cdlp.mean() < 3e-3:
            ax4.set_ylim(0, 3e-3)
        else:
            ax4.set_ylim(0, bestmodPars.Cdlp.max())
        ax5_twin.set_ylim(0.3, 1.2)
        ax5.set_ylabel("Qad")
        ax5.legend(fontsize=10, loc="upper left")
        ax5_twin.legend(fontsize=10, loc="upper right")
        if bestmodPars.Qad.max() < 4e-3:
            ax5.set_ylim(0, 4e-3)
        else:
            ax5.set_ylim(0, bestmodPars.Qad.max())
        ax5.set_xlabel(f"{combined_grouper}")
    except Exception as e:
        print("ERROR EIS Plot combined in extra scatter plots", e)
    #        Pmodgr.plot(x=EvRHE,y='Rs',ax=ax2,label=modnm,kind='scatter')
    #    ax.set_xlim(0,maxYre)
    #    *(ht/wd)
    E_data_combined_path_png = E_data_combined_path_target.with_suffix(".png")
    #    EIS_dest_dir.joinpath(Path(str(Path(EISgr['PAR_file'].unique()[0]).stem)+'_Combined.png' ))
    #    fig_path.parent.mkdir(parents=True,exist_ok=True)
    #    print('EIS fig saved: %s'  %fig_path)
    #    plt.show()
    plt.savefig(E_data_combined_path_png, bbox_inches="tight", dpi=100)
    plt.close()


def plot_lin_Warburg(spec, _lin, lin_slopes, pars=pd.DataFrame(), dest_path=""):

    if not pars.empty:
        pars.to_pickle(dest_path.with_suffix(".pkl"))

    fig, (ax, ax1) = plt.subplots(nrows=2, figsize=(12, 8))
    spec.plot(x="DATA_Zre", y="DATA_-Zim", c="r", ax=ax, label="data")
    #    f'{yax}_slopes'
    if hasattr(spec, "Zre_lin_slopes"):
        spec.plot(
            x="DATA_Zre",
            y="Zre_lin_slopes",
            c="b",
            ax=ax,
            label=f'linear: {_lin["Zre_lin_slopes"]["slope"]:.3f}x + {_lin["Zre_lin_slopes"]["intercept"]:.3f} ',
        )
    for _slope, _color in lin_slopes:
        spec.plot(
            x="DATA_Zre",
            y=f"WB_lin_-Zim_a{_slope}",
            c=_color,
            ax=ax,
            label=f"1 to {1/_slope:.0f}"
            + f' {_lin[f"lin_slope_{_slope}"]["popt"]:.2f}',
            ls="--",
        )
    _lincols = [i for i in spec.columns if "_lintangent" in i]
    if _lincols:
        for freq in set([i[0:-3] for i in _lincols]):
            _re0 = spec[f"{freq}_re"].min()
            spec.plot(
                x=f"{freq}_re",
                y=f"{freq}_im",
                ax=ax,
                ls="--",
                label=f"lintang {freq[-2:]}: {_re0:.1f}",
            )

    #    spec.plot(x='DATA_Zre',y='Z_lin_1t1',c='lightgreen',ax=ax,label=f'1 to 1',ls='--' )
    #    spec.plot(x='DATA_Zre',y='Z_lin_1t2',c='grey',ax=ax,label=f'1 to 2',ls='-.' )
    #    spec.plot(x='DATA_Zre',y='Z_lin_1t4',c='orange',ax=ax,label=f'1 to 4',ls='-.' )
    ax.set_xlabel("Zre")
    ax.set_ylabel("-Zim")
    ax.set_title("Linear check ")
    _angW = "Ang_Warburg"
    #    fig,ax=plt.subplots(figsize=(12,8))
    spec.plot(x=_angW, y="DATA_Zre", c="r", ax=ax1, label="real", lw=4)
    spec.plot(x=_angW, y="DATA_-Zim", c="b", ax=ax1, label="-imag", lw=4)
    if hasattr(spec, "-Zim_lin_WB_angW_high"):
        spec.plot(
            x=_angW,
            y="-Zim_lin_WB_angW_high",
            c="b",
            ax=ax1,
            label=f'{_lin["-Zim_angW_lin_high"]["slope"]:.3f}x + {_lin["-Zim_angW_lin_high"]["intercept"]:.3f}',
            ls=":",
            alpha=0.8,
        )
    if hasattr(spec, "Zre_lin_WB_angW_high"):
        spec.plot(
            x=_angW,
            y="Zre_lin_WB_angW_high",
            c="r",
            ax=ax1,
            label=f'{_lin["Zre_angW_lin_high"]["slope"]:.3f}x + {_lin["Zre_angW_lin_high"]["intercept"]:.3f}',
            ls=":",
            alpha=0.8,
        )
    if hasattr(spec, "-Zim_lin_WB_angW_low"):
        spec.plot(
            x=_angW,
            y="-Zim_lin_WB_angW_low",
            c="b",
            ax=ax1,
            label=f'{_lin["-Zim_angW_lin_low"]["slope"]:.3f}x + {_lin["-Zim_angW_lin_low"]["intercept"]:.3f}',
            ls="--",
            alpha=0.7,
        )
    if hasattr(spec, "Zre_lin_WB_angW_low"):
        spec.plot(
            x=_angW,
            y="Zre_lin_WB_angW_low",
            c="r",
            ax=ax1,
            label=f'{_lin["Zre_angW_lin_low"]["slope"]:.3f}x + {_lin["Zre_angW_lin_low"]["intercept"]:.3f}',
            ls="--",
            alpha=0.7,
        )

    plt.savefig(dest_path, bbox_inches="tight", dpi=300)
    #        set_dest_dir(dest_path).joinpath(_key+'_check_linslope').with_suffix('.png')
    #    spec.plot(x='ang_Warburg',y='W_lin_Zre',c='r',ax=ax,label=f'{_lin["Zre"].slope:.3f}x + {_lin["Zre"].intercept:.3f} ')
    #    plt.show()
    plt.close()

    #%%


#        ax.text(-.05, 0, "$\mathrm{%.2f}$"%Ev, transform=ax.transAxes, ha="left", va="top")
#        ax4.text('%.2f'%Ev)

#        ax.scatter(Edata['DATA_Zre'].values,Edata['DATA_Zim'].values+len(Edata)*[maxYim*En],c=dataC,s=150)
#        ax.plot(Edata['FIT2_Zre'].values,Edata['FIT2_Zim'].values+len(Edata)*[maxYim*En],c=fitC,lw=2.5)
#        axes[En].set_xlim(maxLim)

#        axes[En].axis('equal')
#    fig.subplots_adjust(hspace=0)
#    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)


#    for Erhe,grErhe  in EISgr.groupby(by=EvRHE):
#
#    ax,ax2, ax3, ax4 =axes[0,0], axes[0,1], axes[1,0], axes[1,1]
#    fig.suptitle('%s, %s: \n$\mathrm{E_{DC}\/ =\/ %.2f\/ V_{RHE} \/ in \/%s\/at\/pH\/=\/%.0f}$' %(grErhe['SampleID'].unique()[0],
#    Path(grErhe['PAR_file'].unique()[0]).stem,grErhe[EvRHE].unique()[0],EISovv['Electrolyte'].unique()[0],EISovv['pH'].unique()[::]))
#    dataC,fitC,extraC,initC = 'tab:blue','tab:red','gold','gray'
##        outData.plot(x='FIT_Zre',y='FIT_-Zim',kind='line',c='orange',ax=ax)
#    '''ax: Impedance (Nyquist) Plots'''
##     mod,mod2 = Model(EEC_ORR),Model(EEC_ORRpOx)
#
#
#    outData.plot(x='FIT_2_Zre',y='FIT_2_-Zim',kind='line',c=fitC,lw=2.5,ax=ax,label='FIT_2 (=EEC_ORRpOx)')
#    outData.plot(x='FIT1_Zre',y='FIT1_-Zim',kind='line',c=extraC,lw=2.5,ax=ax,label='FIT_2 (=EEC_ORR)')
#    outData.plot(x='DATA_Zre',y='DATA_-Zim',kind='scatter',c=dataC,s=150,ax=ax)
#    ax.plot(outData['FIT_2_Zre'].values,outData['FIT_2_Zre'].values,ls=':',c='k',lw=1)
##    outData.plot(x='FIT1_Zre',y='FIT1_Zre',kind='line',ls='.',c='k',lw=10,ax=ax)
#    ax.axis('equal')
#    '''ax2: Admittance Plots'''
#    outData.plot(x='FIT_2_Yre',y='FIT_2_Yim',kind='line',c=fitC,lw=2.5,ax=ax2)
#    outData.plot(x='DATA_Yre',y='DATA_Yim',kind='scatter',c=dataC,s=150,ax=ax2)
#    outData.plot(x='FIT1_Yre',y='FIT1_Yim',kind='line',c=extraC,lw=2.5,ax=ax2)
##    outData.plot(x='INIT2_Yre',y='INIT2_Yim',kind='line',c=initC,lw=2.5,ax=ax2)
#    ax2.axis('equal')
#
#    '''ax3 & ax4: Frequency (Bode) Plots'''
#    outData.plot(x='Frequency(Hz)',y='FIT2_Zmod',kind='line',c=fitC,lw=2.5,ax=ax3,logy=True,logx=True)
#    outData.plot(x='Frequency(Hz)',y='FIT1_Zmod',kind='line',c=extraC,lw=2.5,ax=ax3,logy=True,logx=True)
#    outData.plot(x='Frequency(Hz)',y='DATA_Zmod',kind='scatter',c=dataC,s=150,ax=ax3,logy=True,logx=True)
#
#    outData.plot(x='Frequency(Hz)',y='FIT2_Zphase',kind='line',c=fitC,lw=2.5,ax=ax4,logy=False,logx=True)
##    outData.plot(x='Frequency(Hz)',y='INIT2_Zphase',kind='line',c=initC,lw=1.5,ax=ax4,logy=False,logx=True)
#    outData.plot(x='Frequency(Hz)',y='FIT1_Zphase',kind='line',c=extraC,lw=2.5,ax=ax4,logy=False,logx=True)
#    outData.plot(x='Frequency(Hz)',y='DATA_Zphase',kind='scatter',c=dataC,s=150,ax=ax4,logy=False,logx=True)
#
#    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#    ax3.text(2.25, 1, out1.fit_report(min_correl=0.99), transform=ax.transAxes, fontsize=15,
#                verticalalignment='top', bbox=props)
#    ax4.text(2.25, -0.2, out2.fit_report(min_correl=0.99), transform=ax.transAxes, fontsize=15,
#                verticalalignment='top', bbox=props)
#
#    fig_path = EIS_dest_dir.joinpath(Path(str(Path(EISgr['PAR_file'].unique()[0]).stem)+'_%.0fmV.png' %(EISgr[EvRHE].unique()[0]*1000)))
#    fig_path.parent.mkdir(parents=True,exist_ok=True)
##    print('EIS fig saved: %s'  %fig_path)
#    plt.savefig(fig_path,bbox_inches='tight',dpi=100)
##    plt.show()
#    plt.close()
#%%
def EIS_ParsPlotting(Pars2, dest_dir, SaveFigs=True):
    for yPar in ["Rct", "Rs", "Rorr"]:
        fig, ax = plt.subplots(1, 1)
        for fn, Pgr in Pars2.groupby("PAR_file"):
            ax.scatter(
                Pgr["E_AppV_RHE"].values, Pgr[yPar].values, label=Path(fn).stem, s=80
            )
            ax.set_ylabel(yPar)
            ax.set_xlabel("E / V v RHE")
            ax.grid(True)
            plt.legend(bbox_to_anchor=(0.4, 1.24), ncol=2, loc="upper center")
            #                        plt.legend()
            fig.title = yPar
        #        plt.show()
        if SaveFigs:
            plt.savefig(
                dest_dir.joinpath("EIS_%s.png" % (yPar)), dpi=300, bbox_inches="tight"
            )
        plt.close()
    for yPar in [("Cdlp", "nDL"), ("Qad", "nAd")]:
        fig1, ax1 = plt.subplots(1, 1)
        ax2 = ax1.twinx()
        minl, maxl = 0.5 * Pars2[yPar[0]].min(), 1.1 * Pars2[yPar[0]].max()
        for fn, Pgr in Pars2.groupby("PAR_file"):
            ax1.scatter(
                Pgr["E_AppV_RHE"].values,
                Pgr[yPar[0]].values,
                label=str(yPar[0] + ": " + Path(fn).stem),
                s=80,
            )
            ax2.scatter(
                Pgr["E_AppV_RHE"].values,
                Pgr[yPar[1]].values,
                label=str(yPar[1] + ": " + Path(fn).stem),
                s=80,
                marker="^",
                c="orange",
            )
            ax1.set_ylabel(yPar[0])
            ax2.set_ylabel(yPar[1])
            ax1.set_xlabel("E / V v RHE")
            ax1.set_ylim([minl, maxl])
            ax1.grid(True)
            ax1.legend(bbox_to_anchor=(0.4, 1.40), ncol=2, loc="upper center")
            ax2.legend(bbox_to_anchor=(0.4, 1.24), ncol=2, loc="upper center")

        if SaveFigs:
            plt.savefig(
                dest_dir.joinpath("EIS_%s.png" % (yPar[0])),
                dpi=300,
                bbox_inches="tight",
            )
        plt.close()


#                    Pgr.plot(x='E_AppV_RHE',y='Rct',kind='scatter',c= label=Path(fn).name,ax=ax)


def EIS_postprocessing():
    #%%
    Pars2 = pd.read_excel(FindExpFolder("VERSASTAT").DestDir.joinpath("EIS_Pars2.xlsx"))
    Pars1 = pd.read_excel(FindExpFolder("VERSASTAT").DestDir.joinpath("EIS_Pars1.xlsx"))
    EISfitPath = Path(FindExpFolder("VERSASTAT").DestDir / "EIS_FIT")
    EISfitPath.mkdir(parents=True, exist_ok=True)
    yPars = [
        "Rs",
        "Cdlp",
        "nDL",
        "Rct",
        "Qad",
        "nAd",
        "Rorr",
        "Rct_kin",
        "Chisqr",
        "RedChisqr1",
        "RedChisqr2",
    ]
    for y in yPars[::]:
        if y not in Pars2.columns:
            continue
        else:
            pass

        for pH, grph in Pars2.groupby(by="pH"):
            #        grph.plot(x=EvRHE,y='Rct',kind='scatter',title='pH = %s' %pH)
            grph = grph.sort_values(by=["PAR_file", EvRHE])
            #        grph.to_excel(FindExpFolder('VERSASTAT').DestDir.joinpath('EIS_Pars2_pH%s.xlsx' %pH))
            for Gas, grG in grph.groupby(by="Gas"):
                ymin = 0.5 * grG[y].min()

                if grG[y].mean() * 50 < grG[y].max():
                    ymax = grG[y].max()
                else:
                    ymax = grG[y].mean() * 5

                fig, ax = plt.subplots(1, figsize=(10, 8))
                ax.set_title("%s in pH = %s in %s gas" % (y, pH, Gas))
                FitDir = EISfitPath.joinpath(str(pH), Gas)
                FitDir.mkdir(parents=True, exist_ok=True)
                for Fn, Fgr in grG.groupby(by="File"):

                    sID, Elect = (
                        Fgr["SampleID"].unique()[0],
                        Fgr["Electrolyte"].unique()[0],
                    )
                    if Fgr[EvRHE].empty or Fgr[y].empty:

                        continue
                    #                    print(Fn,Fgr[y])
                    ax.plot(Fgr[EvRHE], Fgr[y], label="%s" % (sID), marker="o")
                    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
                    ax.set_xlabel(EvRHE)
                    ax.set_ylabel(y)
                    ax.set_ylim([ymin, ymax])
                    Fgr.to_excel(FitDir.joinpath("%s_%s.xlsx" % (sID, Path(Fn).stem)))
                plt.savefig(
                    EISfitPath.joinpath(str(pH), "%s_%s.png" % (y, Gas)),
                    dpi=300,
                    bbox_inches="tight",
                )
                #                plt.show()
                plt.close()
    #            Egr.plot(x=EvRHE,y='Rct',kind='scatter',title='pH = %s, E = %s' %(pH,Ev))
