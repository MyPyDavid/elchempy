#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:34:53 2020

@author: zmg
"""


import numpy as np
import os
from pathlib import Path
from collections import namedtuple
import matplotlib as mpl
import matplotlib.pyplot as plt
import random as rnd
import math
from math import sin, cos, pi
import pandas as pd

from scipy.optimize import minimize, curve_fit
from scipy.stats import linregress


import torch
import torch.nn.functional as F

plt.rc("text", usetex=False)
plt.rc("font", family="serif", size=12)
mpl.style.use("default")

if __name__ == "__main__":
    from impedance.models.circuits import CustomCircuit
    from impedance.visualization import (
        plot_nyquist,
        plot_bode,
        plot_residuals,
        plot_altair,
    )
    from models import Model_Collection


def reduce_Z_data(spec):
    try:
        spec = spec.sort_values("Frequency(Hz)", ascending=True)
        R_ohm = abs(spec.DATA_Z).min()
        w_min = spec["Angular"].min()
        Zim_min = spec.loc[spec["Angular"] == w_min, "DATA_Z"].values.imag
        C_sub = 1 / (w_min * Zim_min)

        #     (1j*fmin*1E-3)**-1
        spec["DATA_Z_reduce"] = spec.DATA_Z - R_ohm + (1j * spec.Angular * C_sub) ** -1
        spec["DATA_Z_reduce_real"] = spec["DATA_Z_reduce"].values.real
        spec["DATA_Z_reduce_imag"] = -1 * spec["DATA_Z_reduce"].values.imag
        spec["ang_Warburg"] = 1 / (np.sqrt(spec.Angular))
    except Exception as e:
        print(e)

    return spec


def read_xl(xlfile):
    df = pd.read_excel(xlfile, index_col=[0]).sort_values(
        "Frequency(Hz)", ascending=True
    )
    if "Model_EEC" in df.columns:
        mgrp = df.groupby("Model_EEC")

        getgrp = (
            "Model(Singh2015_RQRQR)"
            if "Model(Singh2015_RQRQR)" in mgrp.groups.keys()
            else list(mgrp.groups.keys())[0]
        )
        spec = mgrp.get_group(getgrp)
    #                              mgrp.groups
    else:
        spec = df
    complex_cols = [
        i for i in spec.columns if "+" and "j" in str(spec.head(1)[i].iloc[0])
    ]
    #    spec[complex_cols] =
    spec = spec.assign(
        **{col: spec[col].apply(lambda x: np.complex(x)) for col in complex_cols}
    )
    #    spec[complex_cols].applymap(lambda x: np.complex(x))
    return spec


def add_EIS_data(spec):
    plt.rc("text", usetex=False)
    #    spec['Zcm2'] = spec['DATA_Z']*0.238
    #    plt.plot(np.real(spec['Zcm2'] ), -np.imag(spec['Zcm2']), "o", markersize=10, color="black", label="synth exp")
    spec.plot(x="DATA_Zre" * 0.238, y="DATA_-Zim")
    spec.plot(x="DATA_Z_reduce_real", y="DATA_Z_reduce_imag")

    N_freqs = len(spec)
    Z_exp = spec.DATA_Z.values
    return N_freqs, Z_exp


def read_eis_excel():
    xl_files = list(Path.cwd().parent.rglob("testing_data/*spectrum*xlsx"))
    #    spec = pd.read_excel(xl_files[1],index_col=[0])
    #    converters={'DATA_Z': lambda s: np.complex(s.replace('i', 'j'))}
    all_data = {
        a.stem: {"Filepath": a, "spectrum": reduce_Z_data(read_xl(a))}
        for a in xl_files
        if not "_GP_" in a.name
    }
    specs = [i["spectrum"] for i in all_data.values()]
    return all_data


def check_Warburg(_key, spec):
    _lin = {}
    for yax in ["Zre", "-Zim"]:
        _lin.update(
            {
                yax: linregress(
                    spec.query("ang_Warburg > 0.3").ang_Warburg,
                    spec.query("ang_Warburg > 0.3")["DATA_" + yax],
                )
            }
        )
        spec["W_lin_" + yax] = _lin[yax].slope * spec.ang_Warburg + _lin[yax].intercept

    fig, ax = plt.subplots(figsize=(12, 8))
    spec.plot(x="ang_Warburg", y="DATA_Zre", c="r", ax=ax, label="real")
    spec.plot(x="ang_Warburg", y="DATA_-Zim", c="b", ax=ax, label="-imag")
    spec.plot(
        x="ang_Warburg",
        y="W_lin_-Zim",
        c="b",
        ax=ax,
        label=f'{_lin["-Zim"].slope:.3f}x + {_lin["-Zim"].intercept:.3f} ',
    )
    spec.plot(
        x="ang_Warburg",
        y="W_lin_Zre",
        c="r",
        ax=ax,
        label=f'{_lin["Zre"].slope:.3f}x + {_lin["Zre"].intercept:.3f} ',
    )
    plt.savefig(
        Path.cwd().joinpath(_key + "_check_Warburg").with_suffix(".png"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.show()
    plt.close()


def func_lin(a):
    def func(x, b):
        return a * x + b

    return func


def check_linZ(_key, spec, _lin_window_size=7):
    _lin = {}
    lin_slopes = [(0.25, "lightgreen"), (0.5, "grey"), (1, "orange")]
    #    zip([0.25, 0.5, 1], ['lightgreen','grey','orange'])
    for yax in ["-Zim"]:
        _lin.update(
            {
                yax: linregress(
                    spec.query("Angular < 30")["DATA_Zre"],
                    spec.query("Angular < 30")["DATA_" + yax],
                )
            }
        )
        spec["Z_lin_" + yax] = _lin[yax].slope * spec.DATA_Zre + _lin[yax].intercept

        for _slope, _ in lin_slopes:
            perr_set = 1000
            #            for _win_size in [7,10,15,25]:
            for win in spec.rolling(_lin_window_size):
                popt, pcov = curve_fit(
                    func_lin(_slope), win.DATA_Zre, win["DATA_" + yax]
                )
                perr = np.sqrt(np.diag(pcov))
                #                print(win.index,popt,pcov,perr)
                if perr < perr_set:
                    perr_set = perr
                    best = (_slope, win, popt, perr)
            #            popt, pcov = curve_fit(func_lin(_slope), spec.query('Angular > 30').DATA_Zre, spec.query('Angular > 30')['DATA_'+yax])
            spec[f"Z_lin_a{_slope}"] = func_lin(best[0])(spec.DATA_Zre, best[2][0])
            _lin.update(
                {
                    _slope: {
                        "popt": best[2][0],
                        "win_size": len(best[1]),
                        "perr": best[-1],
                    }
                }
            )

    #        spec['Z_lin_1t4'] = 0.25* spec.DATA_Zre + 12
    #        spec['Z_lin_1t2'] = 0.5* spec.DATA_Zre + -20
    #        spec['Z_lin_1t1'] = 1* spec.DATA_Zre + -95

    fig, ax = plt.subplots(figsize=(12, 8))

    spec.plot(x="DATA_Zre", y="DATA_-Zim", c="r", ax=ax, label="data")
    spec.plot(
        x="DATA_Zre",
        y="Z_lin_-Zim",
        c="b",
        ax=ax,
        label=f'{_lin["-Zim"].slope:.3f}x + {_lin["-Zim"].intercept:.3f} ',
    )
    for _slope, _color in lin_slopes:
        spec.plot(
            x="DATA_Zre",
            y=f"Z_lin_a{_slope}",
            c=_color,
            ax=ax,
            label=f'1 to {1/_slope:.0f} + {_lin[_slope]["popt"]:.2f}',
            ls="--",
        )
    #    spec.plot(x='DATA_Zre',y='Z_lin_1t1',c='lightgreen',ax=ax,label=f'1 to 1',ls='--' )
    #    spec.plot(x='DATA_Zre',y='Z_lin_1t2',c='grey',ax=ax,label=f'1 to 2',ls='-.' )
    #    spec.plot(x='DATA_Zre',y='Z_lin_1t4',c='orange',ax=ax,label=f'1 to 4',ls='-.' )
    ax.set_xlabel("Zre")
    ax.set_ylabel(yax)
    ax.set_title("Linear check")
    plt.savefig(
        Path.cwd().joinpath(_key + "_check_linslope").with_suffix(".png"),
        bbox_inches="tight",
        dpi=300,
    )
    #    spec.plot(x='ang_Warburg',y='W_lin_Zre',c='r',ax=ax,label=f'{_lin["Zre"].slope:.3f}x + {_lin["Zre"].intercept:.3f} ')
    plt.show()
    plt.close()


def compare_O2_N2():
    xl_files = list(Path.cwd().rglob("testing_data/*xlsx"))
    all_data = {
        a.stem: {"Filepath": a, "spectrum": (pd.read_excel(a, index_col=[0]))}
        for a in xl_files
        if "_GP_DRT" in a.name
    }

    _lst = []
    for k, val in all_data.items():
        _spec = val["spectrum"]
        _spec.columns = [k[0:2] + "_" + c for c in _spec.columns]
        _lst.append(_spec)
    DRT_compare = pd.concat(_lst, sort=False, axis=1)

    fig, ax = plt.subplots(figsize=(12, 8))
    DRT_compare.plot(
        x="N2_freq_vec_star", y="N2_gamma_vec_star", c="b", ax=ax, label="N2"
    )
    ax.fill_between(
        DRT_compare["N2_freq_vec_star"],
        DRT_compare["N2_gamma_vec_star"]
        - 3 * np.sqrt(abs(DRT_compare["N2_Sigma_gamma_vec_star"])),
        DRT_compare["N2_gamma_vec_star"]
        + 3 * np.sqrt(abs(DRT_compare["N2_Sigma_gamma_vec_star"])),
        color="0.4",
        alpha=0.25,
    )
    DRT_compare.plot(
        x="O2_freq_vec_star", y="O2_gamma_vec_star", c="r", ax=ax, label="O2"
    )
    DRT_compare.loc[np.isclose(DRT_compare["N2_freq_vec_star"], 0.5, atol=0.05)].plot(
        x="N2_freq_vec_star",
        y="N2_gamma_vec_star",
        c="b",
        ax=ax,
        label="lowest frequency measured",
        kind="scatter",
        s=80,
    )
    DRT_compare.loc[np.isclose(DRT_compare["O2_freq_vec_star"], 0.5, atol=0.05)].plot(
        x="O2_freq_vec_star", y="O2_gamma_vec_star", c="r", ax=ax, kind="scatter", s=80
    )
    ax.fill_between(
        DRT_compare["O2_freq_vec_star"],
        DRT_compare["O2_gamma_vec_star"]
        - 3 * np.sqrt(abs(DRT_compare["O2_Sigma_gamma_vec_star"])),
        DRT_compare["O2_gamma_vec_star"]
        + 3 * np.sqrt(abs(DRT_compare["O2_Sigma_gamma_vec_star"])),
        color="0.4",
        alpha=0.25,
    )
    ax.set_xscale("log")
    ax.set_ylim(-50, 500)
    ax.set_xlabel(r"$f/{\rm Hz}$", fontsize=20)
    ax.set_ylabel(r"$\gamma/\Omega$", fontsize=20)
    ax.set_title(f"{list(all_data.keys())[0]}\n\n")
    plt.savefig(
        Path.cwd()
        .joinpath("testing_data", "GP_DRT_comparison_large")
        .with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
    )


# all_test_data = read_eis_excel()


def choose_test(
    all_test_data,
    name="O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm",
    spec_type="spectrumraw",
    reduce=False,
):

    #    name = 'O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm_spectrumfit_v20'
    jos2 = [i for i in list(all_test_data.keys()) if name in i and spec_type in i]
    #    spec = all_test_data.get('O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm_spectrumfit_v20')['spectrum']
    _key = jos2[0]  # TODO FILE SELECTOR TODO
    spec = all_test_data.get(_key)["spectrum"]
    check_Warburg(_key, spec)
    check_linZ(_key, spec)
    N_freqs = len(spec)
    freq_vec = spec["Frequency(Hz)"].to_numpy()
    Z_exp = spec.DATA_Z.to_numpy()
    if reduce:
        Z_exp = spec.DATA_Z_reduce.values
    print(_key)
    return N_freqs, freq_vec, Z_exp, _key


all_test_data = read_eis_excel()
_test_name_select = "N2_EIS-range_1500rpm_JOS3_288_758mV_1500rpm_3"
N_freqs, freq_vec, Z_exp, _key = choose_test(all_test_data)

#%%
# freq_vec, Z_exp, _key = freq_KKv, Z_KKv, fit_run_arg.PAR_file.name

#  Z_KKv,ang_KKv = EIS_data_KKvalid.DATA_Z.to_numpy(), EIS_data_KKvalid.Angular.to_numpy()
#    freq_KKv = EIS_data_KKvalid['Frequency(Hz)'].to_numpy()

circuit = "R0-p(R1,C1)-p(R2-CPE1,C2)"

circ_C_W = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 90, 247, 1, 4e-4], circuit="R0-p(R1,C1)-p(R2-Ws1,C2)"
)

circ_CPE_W = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 0.7, 90, 247, 1, 4e-4],
    circuit="R0-p(R1,CPE1)-p(R2-Ws1,C2)",
)


circ_RC1_Ws = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 90, 247, 1, 4e-4, 1e-3],
    circuit="R0-p(R1,C1)-p(R2-Ws1,C2)-L",
)


circ_RC1_CPE_Ws = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 0.7, 90, 247, 1, 4e-4, 0.01],
    circuit="R0-p(R1,CPE1)-p(R2-Ws1,C2)-L",
)
best_mod = CustomCircuit(
    initial_guess=[25, 90, 4e-4, 0.7, 100, 10, 0.01, 1e-4],
    circuit="R0-p(R1,CPE1)-p(R2,Ws2)-L0",
)
best_mod_N2 = CustomCircuit(
    initial_guess=[25, 90, 4e-4, 0.7, 100, 10, 0.01, 1e-4],
    circuit="R0-p(R1,CPE1)-p(R2,Ws2)-L0",
)
best_mod_Wser = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 0.7, 90, 4e-4, 0.7, 341, 1, 1e-4],
    circuit="R0-p(R1,CPE1)-p(R2,CPE2)-Wo1-L",
)
best_mod3_RW = CustomCircuit(
    initial_guess=[25, 56, 1e-04, 0.7, 50, 1e-2, 0.9, 500, 1e-03, 1e-5],
    circuit="R0-p(R1,CPE1)-p(R2,CPE2)-p(R3,W3)-L0",
)


# TODO BUILT IN THESE BEST MODELS TO STANDARD FITTING
best_mod_RandlesW = CustomCircuit(
    initial_guess=[25, 100, 3e02, 0.7e-03, 0.7, 1e-4], circuit="R0-p(R1-W1,CPE1)-L0"
)
best_mod_Randles = CustomCircuit(
    initial_guess=[25, 100, 1e-4, 0.5, 0.7e-03, 1e-4], circuit="R0-p(R1-CPE2,C1)-L0"
)


best_mod2_RCPE = CustomCircuit(
    initial_guess=[25, 100, 1e-04, 0.7, 1000, 1e-3, 0.7, 1e-4],
    circuit="R0-p(R1,CPE1)-p(R2,CPE2)-L0",
)

best_mod2_RWpCPE = CustomCircuit(
    initial_guess=[25, 100, 1e-04, 0.7, 400, 4e2, 1e-3, 0.7, 1e-4],
    circuit="R0-p(R1,CPE1)-p(R2-W2,CPE2)-L0",
)

best_mod2_W_2CPE = CustomCircuit(
    initial_guess=[25, 4e2, 100, 1e-04, 0.7, 400, 1e-3, 0.7, 1e-4],
    circuit="R0-W1-p(R1,CPE1)-p(R2,CPE2)-L0",
)

best_UEEC = CustomCircuit(
    initial_guess=[30, 1e-5, 30, 1e-05, 1e-04, 0.7, 25, 1e-04, 0.7, 500, 1e-4],
    circuit="R4-L4-p(R0-L0,CPE0)-p(R1-CPE1,R2-C2)",
)

best_mod3_midC_W3 = CustomCircuit(
    initial_guess=[25, 56, 0.7e-04, 0.7, 50, 1e-2, 560, 2.7e02, 1e-5],
    circuit="R0-p(R1,CPE1)-p(R2,C2)-p(R3,W3)-L0",
)
best_mod3_midC_CPE3 = CustomCircuit(
    initial_guess=[25, 56, 0.7e-04, 0.7, 50, 1e-2, 560, 1.7e-03, 0.5, 1e-5],
    circuit="R0-p(R1,CPE1)-p(R2,C2)-p(R3,CPE3)-L0",
)

models = [
    best_mod_RandlesW,
    best_mod2_RCPE,
    best_mod2_RWpCPE,
    best_mod3_midC_W3,
    best_mod3_midC_CPE3,
][::]
# TODO ===================
circ_RC1_CPE_W = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 0.7, 90, 430, 2e-04, 2e-05],
    circuit="R0-p(R1,CPE1)-p(R2,Wo1)-L",
)
# circ_RC1_CPE_W = CustomCircuit(initial_guess=[25,100, 2E-04, 0.7, 90, 4E-4,0.9,2E-05 ],
#                              circuit='R0-p(R1,CPE1)-p(R2,CPE2)-L')
circ_C1W1_RWo_L = CustomCircuit(
    initial_guess=[25, 80, 2e-06, 90, 430, 2e-04, 2e-05],
    circuit="R0-p(R1,C1)-p(R2,Wo1)-L",
)


circ_RC1_CPE_Wo = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 0.7, 90, 247, 1, 4e-4, 0.9],
    circuit="R0-p(R1,CPE1)-p(R2-Wo1,CPE2)",
)


c_RCPE_RWC = CustomCircuit(
    initial_guess=[25, 50, 0.005, 0.5, 70, 2], circuit="R0-p(R1-Ws1,CPE2)"
)
#%%
# ====== NEW ML MODELS CHECK =======
type1 = CustomCircuit(
    initial_guess=[25, 100, 3e02, 0.7e-03, 0.7], circuit="R0-p(R1-W0,CPE1)", name="1"
)
type2 = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 0.7, 90, 4e-4, 0.9, 300, 1e-05],
    circuit="R0-p(R1,CPE1)-p(R2,CPE2)-W0-L0",
    name="2",
)
type3 = CustomCircuit(
    initial_guess=[25, 1e-05, 100, 2e-04, 0.7, 90, 300, 4e-4, 0.9],
    circuit="R0-L0-p(R1,CPE1)-p(R2-W0,CPE2)",
    name="3",
)
type3a = CustomCircuit(
    initial_guess=[25, 1e-05, 100, 2e-04, 0.7, 90, 4e-4, 0.9],
    circuit="R0-L0-p(R1,CPE1)-p(R2,CPE2)",
    name="3a",
)


type4 = CustomCircuit(
    initial_guess=[25, 100, 0.7e-03, 0.7, 3e02], circuit="R0-p(R1,CPE1)-W0", name="4"
)
type5 = CustomCircuit(
    initial_guess=[25, 100, 3e02, 0.7e-03, 0.7, 4e-4, 0.9],
    circuit="R0-p(R1-W0,CPE1)-CPE2",
    name="5",
)


initial_guess = [0.01, 0.01, 100, 1, 0.05, 100, 1]

# circuit1 = CustomCircuit(circ_C_W, initial_guess=initial_guess)
# circuit1 = CustomCircuit(circ_C_W, initial_guess=initial_guess)
frequencies, Z = freq_vec, Z_exp
f_pred = np.logspace(6, -2)

test_fits = {}
[1, 2, 3][1:2]
models = [
    best_mod_RandlesW,
    best_mod_Randles,
    best_mod2_RCPE,
    best_mod2_RWpCPE,
    best_mod3_midC_W3,
    best_mod3_midC_CPE3,
    best_mod2_W_2CPE,
][1:4]
ok_models = [best_mod, best_mod_Wser, best_mod_N2, best_mod3_RW]
bad_models = [
    circ_C_W,
    circ_RC1_Ws,
    circ_CPE_W,
    circ_RC1_CPE_Ws,
    circ_C1W1_RWo_L,
    circ_RC1_CPE_W,
]
models = [best_UEEC]
type_models = [type1, type2, type3, type3a, type4, type5]


def testing_models(models):
    for mod in models:
        test_fits.update(
            {
                f"{mod.name} {mod.circuit}": {
                    "fit": mod.fit(frequencies, Z),
                    "predict": mod.predict(f_pred),
                    "mod": mod,
                    "name": mod.name,
                    "res_real": (Z - mod.predict(frequencies)).real / np.abs(Z),
                    "res_imag": (Z - mod.predict(frequencies)).imag / np.abs(Z),
                    "MSE": (Z - mod.predict(frequencies)).imag ** 2
                    + (Z - mod.predict(frequencies)).imag ** 2,
                    "params": pd.DataFrame(
                        data=mod.parameters_.reshape(-1, len(mod.parameters_)),
                        columns=mod.get_param_names()[0],
                        index=[f"{mod.name} {mod.circuit}"],
                    ),
                    "params_err": pd.DataFrame(
                        data=mod.conf_.reshape(-1, len(mod.conf_)),
                        columns=[f"{i}_err" for i in mod.get_param_names()[0]],
                        index=[f"{mod.name} {mod.circuit}"],
                    ),
                }
            }
        )
    # circ_C_W.fit(frequencies, Z)
    # circ_CPE_W.fit(frequencies, Z)
    # circ_RC1_CPE_W.fit(frequencies,Z)
    # circ_C_W_fit = circ_C_W.predict(f_pred)
    # circ_CPE_W_fit = circ_CPE_W.predict(f_pred)
    # circ_RC1_CPE_W_fit = circ_RC1_CPE_W.predict(f_pred)
    # Z_fit = circuit.predict(frequencies)

    mod_legends = [
        f'{i} : {np.sum(test_fits[i]["MSE"]):.2f}' for i in list(test_fits.keys())
    ]

    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(4, 2)
    ax1 = fig.add_subplot(gs[:3, :])
    ax2 = fig.add_subplot(gs[3, :])
    # fig, ax = plt.subplots(figsize = (10,10))
    plot_nyquist(ax1, Z, fmt="o")
    for key, val in test_fits.items():
        plot_nyquist(ax1, val["predict"], fmt="-")

        print(val["mod"])
        plot_residuals(
            ax2,
            frequencies,
            val["res_real"],
            val["res_imag"],
            fmt="-",
            y_limits=(-10, 10),
            extra_label=val["name"],
        )
    #    ax2.text(1,6,f'MSE:{np.sum(val["MSE"]):.2f}')
    ax1.set_ylim((0, abs(Z.imag).max() * 2))
    ax1.set_xlim((0, abs(Z.imag).max() * 2))
    ax1.legend(["Data"] + mod_legends)
    ax2.legend(ncol=4)
    plt.savefig(
        Path.cwd().joinpath(_key + "_IMPY_nyquist").with_suffix(".png"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.show()
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 10), nrows=2)
    plot_bode(ax, frequencies, Z, fmt="o")
    for key, val in test_fits.items():
        plot_bode(ax, f_pred, val["predict"], fmt="-", label=key)
    #    print(val['mod'])
    plt.legend(["Data"] + mod_legends)
    # 'C', 'CPE1','RC1'])
    plt.savefig(
        Path.cwd().joinpath(_key + "_IMPY_bode").with_suffix(".png"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.show()
    plt.close()

    params_out = pd.concat(
        [
            pd.concat([val["params"], val["params_err"]], axis=1).assign(
                **{
                    "MSE": np.sum(val["MSE"]),
                    "n_params": len(val["mod"].parameters_),
                    "File": _key,
                }
            )
            for key, val in test_fits.items()
        ],
        sort=False,
    )
    params_out.to_excel(
        Path.cwd().joinpath(_key + "_IMPY_params_TYPES").with_suffix(".xlsx")
    )


# fig, ax = plt.subplots(figsize = (10,10),nrows=2)
# plot_residuals(ax, frequencies, Z, fmt='o')
# for key,val in test_fits.items():
#    plot_residuals(ax, f_pred, val['predict'],fmt='-',label=key)
##    print(val['mod'])
# plt.legend(['Data']+list(test_fits.keys()))
## 'C', 'CPE1','RC1'])
# plt.show()

# plot_nyquist(ax, circ_C_W_fit, fmt='-')
# plot_nyquist(ax, circ_CPE_W_fit, fmt='-')
# plot_nyquist(ax, circ_RC1_CPE_W_fit, fmt='-')

# plt.legend(['Data']+list(test_fits.keys()))
# 'C', 'CPE1','RC1'])
# plt.show()
# print(circ_C_W)
# print(circ_CPE_W)
# print(circ_RC1_CPE_W)
