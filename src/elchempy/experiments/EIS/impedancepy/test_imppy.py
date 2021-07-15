#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:34:53 2020

@author: zmg
"""


import numpy as np
import os
from pathlib import Path
from collections import namedtuple, OrderedDict
import matplotlib as mpl
import matplotlib.pyplot as plt
import random as rnd
import math
from math import sin, cos, pi
import pandas as pd
import random
from scipy.optimize import minimize, curve_fit
from scipy.stats import linregress
import sys
import pickle
import torch
import torch.nn.functional as F
import datetime as dt

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

    sys.path.append(str(Path(__file__).parent.parent))
    from models import Model_Collection

_skip_files = ["N2_EIS-range_0rpm_JOS2_272_650mV_0rpm_5"]
_test_name_select = [
    "N2_EIS-range_1500rpm_JOS2_288_758mV_1500rpm",
    "O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm",
    "O2_EIS-range_1500rpm_JOS4_270_650mV_1500rpm_5",
    "O2_EIS-range_1500rpm_JOS4_268_188mV_1500rpm_11",
    "O2_EIS-range_1500rpm_JOS13_211_614mV_1500rpm_5",
    "O2_EIS-range_1500rpm_JOS12_211_464mV_0rpm_8",
    "O2_EIS-rpm-range-400mV_JOS12_233_633mV_1500rpm_3",
    "O2_EIS-range_1500rpm_JOS4_low-load_263_700mV_1472rpm_4",
] + _skip_files
dest_path = Path.cwd().parent.joinpath("testing_data/impy_results")
dest_path_complex = Path.cwd().parent.joinpath("testing_data/impy_complex_plots")
dest_path_complex.mkdir(exist_ok=True)
dest_path_params = dest_path_complex.joinpath("params")
dest_path_params.mkdir(exist_ok=True)


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
        spec = spec.reset_index()
    except Exception as e:
        print(e)

    return spec


def set_dest_dir(dest_path):
    _dest_dir = Path.cwd()
    if dest_path:
        if Path(dest_path).is_dir():
            _dest_dir = Path(dest_path)
    return _dest_dir


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
    xl_files = list(Path.cwd().parent.rglob("testing_data/spectrumfits/*spectrum*xlsx"))
    #    spec = pd.read_excel(xl_files[1],index_col=[0])
    #    converters={'DATA_Z': lambda s: np.complex(s.replace('i', 'j'))}
    spec_files = [i for i in xl_files if not "_GP_" in i.name]
    # all_data = {a.stem : {'Filepath' : a, 'spectrum' : reduce_Z_data(read_xl(a))} for a in _spec_files}
    set_sfls = set([i.stem for i in spec_files])
    _pickle_path = Path.cwd().parent.joinpath("testing_data/spec_pickle.pkl")
    all_data = {}
    if _pickle_path.is_file():
        try:
            with open(_pickle_path, "rb") as handle:
                all_data = pickle.load(handle)
        except Exception as e:
            print("Load error", e)
        if not set(all_data.keys()) == set(set_sfls):
            all_data = {}

        # set(all_data.keys())
    if not all_data:
        all_data = {
            a.stem: {"Filepath": a, "spectrum": reduce_Z_data(read_xl(a))}
            for a in spec_files
        }
        # Store data (serialize)
        with open(_pickle_path, "wb") as handle:
            pickle.dump(all_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # specs = [i['spectrum'] for i in all_data.values()]
    return all_data


def check_Warburg(_key, spec, _lin_window_size=15, dest_path="", export_plot=False):
    _lin = {}
    for yax in ["Zre", "-Zim"]:
        _lin.update(
            {
                yax: linregress(
                    spec.query("ang_Warburg > 0.3").ang_Warburg,
                    spec.query("ang_Warburg > 0.3")["DATA_" + yax],
                ),
                yax
                + "_low": linregress(
                    spec.query("ang_Warburg < 0.045").ang_Warburg,
                    spec.query("ang_Warburg < 0.045")["DATA_" + yax],
                ),
            }
        )

        spec = spec.assign(
            **{
                "W_lin_"
                + yax: _lin[yax].slope * spec.ang_Warburg
                + _lin[yax].intercept,
                "W_lin_"
                + yax
                + "_low": _lin[yax + "_low"].slope * spec.ang_Warburg
                + _lin[yax + "_low"].intercept,
            }
        )
        # spec['W_lin_'+yax] = _lin[yax].slope * spec.ang_Warburg + _lin[yax].intercept
        # spec['W_lin_'+yax+'_low'] = _lin[yax+'_low'].slope * spec.ang_Warburg + _lin[yax+'_low'].intercept
    #    for win in spec.rolling(_lin_window_size):
    #        for yax in ['Zre','-Zim']:
    #            popt, pcov = curve_fit(func_lin(_slope), win.DATA_Zre, win['DATA_'+yax])
    #            perr = np.sqrt(np.diag(pcov))
    # TODO ADD extra linear fits
    if export_plot:
        fig, ax = plt.subplots(figsize=(12, 8))
        spec.plot(x="ang_Warburg", y="DATA_Zre", c="r", ax=ax, label="real", lw=4)
        spec.plot(x="ang_Warburg", y="DATA_-Zim", c="b", ax=ax, label="-imag", lw=4)
        spec.plot(
            x="ang_Warburg",
            y="W_lin_-Zim",
            c="b",
            ax=ax,
            label=f'{_lin["-Zim"].slope:.3f}x + {_lin["-Zim"].intercept:.3f}',
            ls=":",
            alpha=0.8,
        )
        spec.plot(
            x="ang_Warburg",
            y="W_lin_Zre",
            c="r",
            ax=ax,
            label=f'{_lin["Zre"].slope:.3f}x + {_lin["Zre"].intercept:.3f}',
            ls=":",
            alpha=0.8,
        )

        spec.plot(
            x="ang_Warburg",
            y="W_lin_-Zim_low",
            c="b",
            ax=ax,
            label=f'{_lin["-Zim_low"].slope:.3f}x + {_lin["-Zim_low"].intercept:.3f}',
            ls="--",
            alpha=0.7,
        )
        spec.plot(
            x="ang_Warburg",
            y="W_lin_Zre_low",
            c="r",
            ax=ax,
            label=f'{_lin["Zre_low"].slope:.3f}x + {_lin["Zre_low"].intercept:.3f}',
            ls="--",
            alpha=0.7,
        )
        plt.savefig(
            set_dest_dir(dest_path)
            .joinpath(_key + "_check_Warburg")
            .with_suffix(".png"),
            bbox_inches="tight",
            dpi=300,
        )
        #    plt.show()
        plt.close()

    spec = spec.assign(
        **{
            "Zre_Rs": spec.DATA_Zre / spec.DATA_Zre.min(),
            "-Zim_Rs": spec["DATA_-Zim"] / spec.DATA_Zre.min(),
            "Z_bode_phase_ang": np.tan(
                np.abs(spec["DATA_-Zim"]) / np.abs(spec.DATA_Zre)
            ),
        }
    )
    return spec


def func_lin(a):
    def func(x, b):
        return a * x + b

    return func


def check_linZ(_key, spec, _lin_window_size=7, dest_path="", export_plot=False):
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
            #            for win in spec.rolling(_lin_window_size):
            for i in range((len(spec) - _lin_window_size)):
                win = spec.iloc[i : i + _lin_window_size]
                #                print(win.index)
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
    if export_plot:
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
            set_dest_dir(dest_path)
            .joinpath(_key + "_check_linslope")
            .with_suffix(".png"),
            bbox_inches="tight",
            dpi=300,
        )
        #    spec.plot(x='ang_Warburg',y='W_lin_Zre',c='r',ax=ax,label=f'{_lin["Zre"].slope:.3f}x + {_lin["Zre"].intercept:.3f} ')
        plt.close()


def semicircle_func(x, a, b, aR, bR):
    y = np.sqrt(a - x ** 2) * aR + np.sqrt(b - x ** 2) * bR
    return y


def check_semicircle(_key, spec):
    _semi = []
    yax = ["-Zim"][0]
    popt, pcov = curve_fit(semicircle_func, spec.DATA_Zre, spec["DATA_" + yax])
    spec["Z_lin_a{_slope}"] = semicircle_func(spec.DATA_Zre)


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
    plt.close()


# all_test_data = read_eis_excel()


def choose_test(
    all_test_data,
    name="O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm",
    spec_type="spectrumraw",
    reduce=False,
    dest_path=dest_path,
    freqlim=1e5,
    export_plot=False,
):
    jos2 = [i for i in list(all_test_data.keys()) if name in i and spec_type in i]
    _key = jos2[0]  # TODO FILE SELECTOR TODO
    spec = all_test_data.get(_key)["spectrum"]
    spec = spec.iloc[1::]
    spec = check_Warburg(_key, spec, dest_path=dest_path, export_plot=export_plot)
    check_linZ(_key, spec, dest_path=dest_path, export_plot=export_plot)
    _spec_lim = spec.loc[spec["Frequency(Hz)"] < freqlim]
    N_freqs = len(_spec_lim)
    freq_vec = _spec_lim["Frequency(Hz)"].to_numpy()
    Z_exp = _spec_lim.DATA_Z.to_numpy()
    if reduce:
        Z_exp = _spec_lim.DATA_Z_reduce.values
    print(_key)
    return N_freqs, freq_vec, Z_exp, _key, spec


def Z_to_Y(Z):
    return np.array([np.real(j) + -1j * np.imag(j) for j in [i ** -1 for i in Z]])


all_test_data = read_eis_excel()


N_freqs, freq_vec, Z_exp, _key, spec = choose_test(
    all_test_data, name=_test_name_select[2], spec_type="spectrumfit"
)

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

best_mod_Randles_Rorr = CustomCircuit(
    initial_guess=[10, 1e-5, 0.8, 38, 150, 1e-03, 0.8],
    circuit="R0-p(CPE0,R1-p(R2,CPE1))",
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
# type1 = CustomCircuit(initial_guess=[25,5E-05,100, 300,2, 0.7E-03,0.7],
#                              circuit='R0-L0-p(R1-Ws0,CPE1)',name='1')
type1 = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 0.7e-03, 0.7, 3e-04, 0.7],
    circuit="R0-L0-p(R1-W0,CPE1)-CPE2",
    name="1",
)

type1C = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 0.7e-03, 0.7, 3e-04],
    circuit="R0-L0-p(R1-W0,CPE1)-C2",
    name="1+C",
)

# type1C = CustomCircuit(initial_guess=[25,5E-05,100,300,2, 0.7E-03,0.7,3E02],
#                               circuit='R0-L0-p(R1-Wo0,CPE1)-W0',name='1+C')

# 'R0-L0-p(R1-Wo0,CPE1)-W0'
# type1b = CustomCircuit(initial_guess=[25,5E-05,100,300, 0.7E-03,0.7,3E-04,0.7],
#                              circuit='R0-L0-p(R1-W0-CPE2,CPE1)',name='1b')
# type1b =
# CustomCircuit(initial_guess=[20,5E-05,30,300,0.5, 0.7E-03,0.7,50, 30, 0.5,3E-04,0.7],
#                              circuit='R0-L0-p(R1-Wo0-CPE2,R2-Ws0-CPE1)',name='1b') # slecht
typeRandlesC = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 0.7e-03, 3e02, 0.5],
    circuit="R0-L0-p(R1,C1)-Wo1",
    name="Randles+Wo",
)


typeRandlesCPE = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 0.7e-03, 0.7, 3e02, 0.5],
    circuit="R0-L0-p(R1,CPE1)-Wo1",
    name="RandlesCPE+Wo",
)
# 'R0-p(R1,C1)-Wo1'

type1RW_C = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 2, 0.7e-03, 0.7, 3e-04, 0.7],
    circuit="R0-L0-p(R1-Wo0,CPE1)-CPE2",
    name="1RW+CPE",
)
type1C_RC = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100, 3e-03],
    circuit="R0-L0-p(R1-Wo0,CPE1)-p(R2,C2)",
    name="1c+RC",
)
type1RWoCPE_C2 = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 0.2, 3e-04, 0.7, 3e-03],
    circuit="R0-L0-p(R1-Wo0,CPE1)-C2",
    name="1RWoCPE+C",
)
type1RWsCPE_C2 = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 0.2, 3e-04, 0.7, 3e-03],
    circuit="R0-L0-p(R1-Ws0,CPE1)-C2",
    name="1RWsCPE+C",
)


type1RCPE_C = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 3e-04, 0.7, 3e-03],
    circuit="R0-L0-p(R1,CPE1)-C2",
    name="1RCPE+C",
)


# type3C, type3CW, type3CWo]

type1C_RCWs = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100, 3e-03],
    circuit="R0-L0-p(R1-Ws0,CPE1)-p(R2,C2)",
    name="1c+RWs+C",
)

type1C_W = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100],
    circuit="R0-L0-p(R1-Wo0,CPE1)-W0",
    name="1c+W",
)
type1C_RCPE = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100, 3e-03, 0.7],
    circuit="R0-L0-p(R1-Wo0,CPE1)-p(R2,CPE2)",
    name="1c+RCPE",
)
type1C_CPE = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 3e-03, 0.7],
    circuit="R0-L0-p(R1-Wo0,CPE1)-CPE2",
    name="1c+CPE",
)


# type_models = [type1, type1C,type1C_RC,type1C_RCWs, type1C_W, type3C, type3CW, type3CWo]
# type1 = CustomCircuit(initial_guess=[25,5E-05,100,300,2,1E-3,0.7, 3E02,0.7, 1E-03,0.7],
#                              circuit='R0-L0-p(R1-Ws0,CPE1)-Wo0-CPE2',name='1 test')


type2 = CustomCircuit(
    initial_guess=[25, 100, 2e-04, 0.7, 90, 4e-4, 0.9, 1e-05, 300],
    circuit="R0-p(R1,CPE1)-p(R2,CPE2)-L0-W0",
    name="2",
)

type3 = CustomCircuit(
    initial_guess=[25, 5e-06, 14, 2e-04, 0.7, 90, 300, 4e-4, 0.9],
    circuit="R0-L0-p(R1,CPE1)-p(R2-W0,CPE2)",
    name="3",
)

type3b = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 2e-04, 0.7, 90, 300, 0.5, 4e-4, 0.9],
    circuit="R0-L0-p(R1,CPE1)-p(R2-Ws0,CPE2)",
    name="3+L",
)
type3b = CustomCircuit(
    initial_guess=[25, 5e-05, 50, 10, 0.05, 2e-04, 0.7, 500, 200, 0.67, 4e-4, 0.9],
    circuit="R0-L0-p(R1-Wo0,CPE1)-p(R2-Ws0,CPE2)",
    name="3-Ws",
)

type3C = CustomCircuit(
    initial_guess=[25, 5e-06, 14, 2e-04, 0.7, 90, 100, 4e-4],
    circuit="R0-L0-p(R1,CPE1)-p(R2-W0,C2)",
    name="3RWC",
)
type3CWo = CustomCircuit(
    initial_guess=[25, 5e-06, 14, 2e-04, 0.7, 90, 100, 0.5, 4e-4],
    circuit="R0-L0-p(R1,CPE1)-p(R2-Wo0,C2)",
    name="3RWoC",
)


type3CW = CustomCircuit(
    initial_guess=[25, 5e-06, 14, 2e-04, 0.7, 90, 4e-4, 300],
    circuit="R0-L0-p(R1,CPE1)-p(R2,C2)-W0",
    name="3C-W",
)


# type3b = CustomCircuit(initial_guess=[25, 100, 2E-04, 0.7, 90, 4E-4,0.9 ],
#                              circuit='R0-p(R1,CPE1)-p(R2,CPE2)',name='3-Ws')


type3a = CustomCircuit(
    initial_guess=[25, 1e-05, 100, 2e-04, 0.7, 90, 4e-4, 0.9],
    circuit="R0-L0-p(R1,CPE1)-p(R2,CPE2)",
    name="3a",
)


type4 = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 0.7e-03, 0.7, 3e02],
    circuit="R0-L0-p(R1,CPE1)-W0",
    name="4",
)
type5 = CustomCircuit(
    initial_guess=[25, 100, 3e02, 0.67, 0.7e-03, 0.7, 4e-4, 0.9, 5e-05],
    circuit="R0-p(R1-Ws0,CPE1)-CPE2-L0",
    name="5",
)

type1RWsCPE_C2 = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 300, 0.2, 3e-04, 0.7, 3e-03],
    circuit="R0-L0-p(R1-Ws0,CPE1)-C2",
    name="1RWsCPE+C",
)

type1RWsC = CustomCircuit(
    initial_guess=[25, 100, 300, 1, 3e-04], circuit="R0-p(R1-Ws0,C1)", name="1RWsC"
)


t1_RTCPE = CustomCircuit(
    initial_guess=[25, 5e-05, 20, 300, 0.2, 0.7, 3e-04, 0.7],
    circuit="R0-L0-p(T0,CPE1)",
    name="R0-L0-p(T0,CPE1)",
)
t1_RT_CPE = CustomCircuit(
    initial_guess=[25, 5e-05, 20, 300, 0.2, 0.7, 3e-04, 0.7],
    circuit="R0-L0-T0-CPE1",
    name="R0-L0-T0-CPE1",
)


t1_Gs = CustomCircuit(
    initial_guess=[25, 5e-05, 100, 0.2, 1], circuit="R0-L0-Gs0", name="Gs"
)

t1_RG0CPE = CustomCircuit(
    initial_guess=[25, 5e-05, 75, 1e3, 5, 6e-06, 0.7],
    circuit="R0-L0-p(R1-G0,CPE1)",
    name="R-G,CPE",
)

t1_RG0CPER2 = CustomCircuit(
    initial_guess=[25, 5e-05, 75, 1e3, 5, 6e-06, 0.7, 400],
    circuit="R0-L0-p(R1-G0,CPE1,R2)",
    name="R-G,CPE,R2",
)


t1_RG0C = CustomCircuit(
    initial_guess=[25, 5e-05, 75, 1e3, 5, 6e-06],
    circuit="R0-L0-p(R1-G0,C1)",
    name="R-G,C",
)
guesses = {
    "Rs": 20,
    "Ls": 5e-5,
    "Rct": 95,
    "Cdlp": 7e-05,
    "R_G": 2e3,
    "t_G": 20,
    "phi_G": 1,
}
t1_RGsC = CustomCircuit(
    initial_guess=[25, 5e-05, 95, 2e3, 20, 1, 3e-05],
    circuit="R0-L0-p(R1-Gs0,C1)",
    name="R-Gs,C",
)
t1_RGsCPE = CustomCircuit(
    initial_guess=[25, 5e-05, 95, 2e3, 20, 1, 3e-05, 0.7],
    circuit="R0-L0-p(R1-Gs0,CPE1)",
    name="R-Gs,CPE",
)


t1_RGsCPE_C = CustomCircuit(
    initial_guess=[25, 5e-05, 20, 100, 0.2, 1, 3e-05, 0.7, 1e-03],
    circuit="R0-L0-p(R1-Gs0,CPE1)-C2",
    name="R-Gs,CPE-C",
)

H_a = CustomCircuit(
    initial_guess=[5e-04, 30, 10, 5e-05, 100, 100, 1],
    circuit="p(C1,R1,R0-L0-p(R2,Ws0))",
    name="H_a",
)
H_b = CustomCircuit(
    initial_guess=[5e-04, 30, 200, 100, 1], circuit="p(C1,R1-p(R2,Ws0))", name="H_b"
)
H_c = CustomCircuit(
    initial_guess=[20, 5e-05, 5e-04, 0.7, 50, 100, 0.5],
    circuit="R0-L0-p(CPE1,R1,Ws0)",
    name="H_c, CPE",
)

mech = CustomCircuit(
    initial_guess=[25, 1e-05, 100, 2e-04, 20, 30, 3e-3, 50],
    circuit="R0-L0-W0-p(C0,p(R2,R1-L1)-W1)",
    name="mechFC",
)

macro = CustomCircuit(
    initial_guess=[0.5, 25, 70, 0.05, 1e-04], circuit="R0-T0", name="T,macro"
)

macro_Ls = CustomCircuit(
    initial_guess=[0.5, 5e-05, 25, 70, 0.05, 1e-04],
    circuit="R0-L0-T0",
    name="L-T,macro",
)
macro_CPE = CustomCircuit(
    initial_guess=[0.5, 5e-05, 25, 70, 0.05, 1e-04, 5e-04, 0.7],
    circuit="R0-L0-p(T0,CPE1)",
    name="L-T,CPE,macro",
)

FC_ice = CustomCircuit(
    initial_guess=[20, 5e-04, 50, 100, 0.5, 100, 5e-04],
    circuit="R0-p(C1,R2-Ws1,R3-C2)",
    name="FC_ice",
)
FC_iceWs_CPE = CustomCircuit(
    initial_guess=[20, 5e-05, 5e-06, 50, 100, 1, 100, 5e-04, 0.7],
    circuit="R0-L0-p(C1,R2-Ws1,R3-CPE2)",
    name="FC_ice Ws,CPE",
)
FC_iceW_CPE = CustomCircuit(
    initial_guess=[20, 5e-05, 5e-06, 50, 100, 100, 5e-04, 0.7],
    circuit="R0-L0-p(C1,R2-W1,R3-CPE2)",
    name="FC_ice W,CPE",
)
FC_iceWo_CPE = CustomCircuit(
    initial_guess=[20, 5e-05, 5e-06, 50, 100, 0.5, 100, 5e-04, 0.7],
    circuit="R0-L0-p(C1,R2-Wo1,R3-CPE2)",
    name="FC_ice Wo,CPE",
)

FC_ice_L = CustomCircuit(
    initial_guess=[20, 5e-04, 50, 100, 0.5, 100, 5e-04],
    circuit="R0-p(C1,R2-Ws1,R3-C2)",
    name="FC_ice",
)

FC_iceRandles = CustomCircuit(
    initial_guess=[20, 5e-05, 5e-04, 50, 100, 0.5, 5],
    circuit="R0-L0-p(C1,R2-Ws1,R2)",
    name="FC_iceRandles",
)

FC_iceG_L = CustomCircuit(
    initial_guess=[20, 5e-04, 50, 100, 0.5, 100, 5e-04],
    circuit="R0-p(C1,R2-G1,R3-C2)",
    name="FC_iceG",
)

initial_guess = [0.01, 0.01, 100, 1, 0.05, 100, 1]

# circuit1 = CustomCircuit(circ_C_W, initial_guess=initial_guess)
# circuit1 = CustomCircuit(circ_C_W, initial_guess=initial_guess)


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
type_models = [type1, type2, type3, type3b, best_mod3_midC_CPE3, type3a, type4, type5]
type_models = [type3]
type_models = [type1, type2]
type_models = [type1, type2, type3, type4, type3C, type3CW]
type_models = [
    type1,
    type1C,
    type1C_RC,
    type1C_RCWs,
    type1C_W,
    type3C,
    type3CW,
    type3CWo,
]
type_models = [
    type1,
    type1C,
    type1C_RC,
    type1RCPE_C,
    type4,
    typeRandlesC,
    typeRandlesCPE,
]
type_models = [type1C, type1RWoCPE_C2, type1RWsCPE_C2]
type_models = [t1_Gs, type1RWsCPE_C2]
type_models = [type1RWsCPE_C2, t1_RTCPE, t1_RT_CPE, t1_Gs, t1_RGsCPE_C]
type_models = [mech, t1_RT_CPE, t1_RG0CPE, t1_RG0C]
type_models = [
    macro,
    macro_Ls,
    macro_CPE,
    t1_RG0CPE,
    t1_RG0C,
    t1_RG0CPER2,
    FC_ice,
    FC_iceG_L,
    type1RWsC,
    FC_iceWs_CPE,
    FC_iceWs_CPE,
][-3:]
type_models = [FC_iceW_CPE, FC_iceWs_CPE, FC_iceWo_CPE, H_c]
""" bad models:
    [t1_RGsC, t1_RGsCPE, FC_iceRandles]
"""


def test_circuit_predict():
    circuit = CustomCircuit(
        initial_guess=[15, 1e-4, 0.7, 20, 1e-03, 0.8, 150],
        circuit="R0-p(CPE1,R1-p(CPE2,R2))",
    )
    Z = circuit.predict(np.logspace(4, -1))
    fig, ax = plt.subplots(figsize=(5, 5))
    plot_nyquist(ax, Z, fmt=".")

    bad_mod_Randles_Rorr = CustomCircuit(
        initial_guess=[1e1, 1e-5, 0.8, 10, 1e-05, 0.8, 10],
        circuit="R0-p(CPE1,R1-p(CPE2,R2))",
    )
    badZ = bad_mod_Randles_Rorr.predict(np.logspace(4, -1))

    circuit = CustomCircuit(
        initial_guess=[10, 1e-4, 0.7, 20, 1e-03, 0.8, 150],
        circuit="R0-p(CPE1,R1-p(CPE2,R2))",
    )
    Z_predict = circuit.predict(np.logspace(4, -1))
    fig, ax = plt.subplots(figsize=(5, 5))
    plot_nyquist(ax, Z_predict, fmt=".")

    circuit = CustomCircuit(
        initial_guess=[20, 1e-4, 0.7, 20, 1e-03, 0.8, 150],
        circuit="R0-p(CPE1,R1-p(CPE2,R2))",
    )
    Z_predict = circuit.predict(np.logspace(4, -1))
    fig, ax = plt.subplots(figsize=(5, 5))
    plot_nyquist(ax, Z_predict, fmt="-")

    circuit = CustomCircuit(
        initial_guess=[20, 1e-4, 0.7, 20, 1e-03, 0.8, 150],
        circuit="R0-p(CPE1,R1-p(CPE2,R2))",
    )
    Z_predict = circuit.predict(np.logspace(4, -1))
    fig, ax = plt.subplots(figsize=(5, 5))
    plot_nyquist(ax, Z_predict, fmt="-")

    fig, ax1 = plt.figure(figsize=(8, 8))
    plot_nyquist(ax1, Z, fmt="o")
    gs = fig.add_gridspec(6, 4)
    ax1 = fig.add_subplot(gs[:3, 0:2])
    axY = fig.add_subplot(gs[:3, 2:4])
    ax2 = fig.add_subplot(gs[3, :])
    # fig, ax = plt.subplots(figsize = (10,10))
    plot_nyquist(ax1, Z, fmt="o")
    plot_nyquist(axY, Z_to_Y(Z), fmt="o")
    # plot_nyquist(ax1, badZ, fmt='o')
    # plot_nyquist(axY, Z_to_Y(badZ), fmt='o')


def chi_sq(Z, Zcalc, nparams):
    #    Zcalc = mod.predict(frequencies)
    _chi_sq = np.sum(
        ((Z.real - Zcalc.real) / Zcalc.real) ** 2
        + ((Z.imag - Zcalc.imag) / Zcalc.imag) ** 2
    )
    _chi_sq_norm = _chi_sq / (2 * len(Z) - nparams)
    return {"chisqr": np.round(_chi_sq, 4), "redchi": np.round(_chi_sq_norm, 6)}


def testing_models(
    run_models, *_args, dest_path="", best_params=pd.DataFrame(), _tfile=""
):
    test_fits = OrderedDict()
    #    print(args)
    N_freqs, freq_vec, Z_exp, _key, spec = _args
    _dest_dir = set_dest_dir(dest_path)

    _dest_dir.joinpath("_model_fits_dict.pkl")
    _dest_file = _dest_dir.joinpath(_key + "_IMPY_params_TYPES").with_suffix(".xlsx")
    _dest_spectra_file = _dest_dir.joinpath(_key + "_IMPY_spectra_TYPES").with_suffix(
        ".xlsx"
    )

    if not best_params.empty:
        _preparams = best_params
    else:
        if _dest_file.is_file():
            _preparams = pd.read_excel(
                _dest_dir.joinpath(_key + "_IMPY_params_TYPES").with_suffix(".xlsx"),
                index_col=[0],
            )
        else:
            _preparams = pd.DataFrame()

    frequencies, Z = freq_vec, Z_exp
    _ang = frequencies * 2 * np.pi
    f_pred = np.logspace(5, -2)
    _angpred = f_pred * 2 * np.pi
    Rs_guess = np.min(Z_exp.real) if 0.5 < np.min(Z_exp.real) < 200 else 20

    # modc = Model_Collection()
    # + modc.lmfit_models

    for mod in run_models:

        if "impedance.models.circuits.circuits.CustomCircuit" in str(type(mod)):
            try:
                try:  # input initial guesses from previous fits
                    if not _preparams.loc[_preparams.mod_name == mod.name].empty:
                        _parguess = [
                            _preparams.loc[_preparams.mod_name == mod.name, i].iloc[0]
                            for i in mod.get_param_names()[0]
                        ]
                        mod.initial_guess = [
                            i * 1e-2 * random.randint(80, 120) for i in _parguess
                        ]  # add random noise
                except KeyError:
                    print("key not in params")
                except AttributeError:
                    print("attr not in params")

                if mod.get_param_names()[0][0] == "R0":
                    mod.initial_guess[0] = Rs_guess

                mod.fit(frequencies, Z)
                _predict = mod.predict(frequencies)
                test_fits.update(
                    {
                        f"{mod.name}: {mod.circuit}": {
                            "fit": mod,
                            "predict": mod.predict(f_pred),
                            "mod_fit": _predict,
                            "mod_data": Z,
                            "mod_data_freq": frequencies,
                            "mod_obj": mod,
                            "mod_name": mod.name,
                            "n_params": len(mod.parameters_),
                            "mod_circuit": mod.circuit,
                            "res_real": (Z - _predict).real / np.abs(Z),
                            "res_imag": (Z - _predict).imag / np.abs(Z),
                            "MSE": (Z - mod.predict(frequencies)).imag ** 2
                            + (Z - _predict).imag ** 2,
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
                            **chi_sq(Z, _predict, len(mod.parameters_)),
                        }
                    }
                )
            except RuntimeError:
                print(f"RuntimeError: {_key} with {mod.name}")
            except ValueError:
                print(f"ValueError: {_key} with {mod.name}")

        elif "models." in str(type(mod)):
            if "Rs" in mod.parameters_guesses.keys():
                mod.parameters_guesses["Rs"].set(value=Rs_guess)

            methods = ["least_squares", "differential_evolution", "ampgo"]
            # mod.initial_guess[0] = Rs_guess
            try:

                try:  # input initial guesses from previous fits
                    if not _preparams.loc[_preparams.mod_name == mod.name].empty:
                        _parguess = {
                            i: _preparams.loc[_preparams.mod_name == mod.name, i].iloc[
                                0
                            ]
                            for i in mod.parameters_guesses.keys()
                        }
                        _parguess = {
                            k: val * 1e-2 * random.randint(80, 120)
                            for k, val in _parguess.items()
                        }  # add random noise
                        for k, val in _parguess.items():
                            mod.parameters_guesses[k].set(value=val)
                except KeyError:
                    print("key not in params")
                except AttributeError:
                    print("attr not in params; {_key} with {mod.name}")

                DataWeights_modulus_Y = np.sqrt((Z.real ** 2 + Z.imag ** 2)) ** -1

                weights_store = {
                    "weight_Y_mod": DataWeights_modulus_Y,
                    "weight_Z_mod": DataWeights_modulus_Y ** -1,
                }

                weights_name = "weight_Z_mod"
                result = mod.model.fit(
                    Z,
                    mod.parameters_guesses,
                    ang=_ang,
                    method=methods[0],
                    weights=weights_store.get("weights_name"),
                )

                _eval = result.eval(ang=_ang)

                outP = result.best_values
                out_params_stderss = [
                    (i + "_stderr", result.params.__getitem__(i).stderr)
                    for i in result.params
                ]
                out_params_correl = [
                    (i + "_correl", result.params.__getitem__(i).correl)
                    for i in result.params
                ]
                par_errs = dict(
                    zip(
                        [i[0] for i in out_params_stderss],
                        [i[1] for i in out_params_stderss],
                    )
                )
                paramDF = pd.DataFrame(result.best_values, index=[f"{mod.name}"])
                paramDF_err = pd.DataFrame(par_errs, index=[f"{mod.name}"])
                # pd.DataFrame(data=mod.parameters_.reshape(-1,len(mod.parameters_)),columns=mod.get_param_names()[0],
                # index=[f'{mod.name}'])
                # 'params_err' : pd.DataFrame(data=mod.conf_.reshape(-1,len(mod.conf_)),columns=[f'{i}_err' for i in mod.get_param_names()[0]],
                # index=[f'{mod.name}'])

                test_fits.update(
                    {
                        f"{mod.name}": {
                            "fit": result,
                            "predict": result.eval(ang=_angpred),
                            "mod_fit": result.best_fit,
                            "mod_data": result.data,
                            "mod_data_freq": frequencies,
                            "mod_obj": mod,
                            "mod_name": mod.name,
                            "mod_circuit": mod.name,
                            "mod_color": [float(i) for i in mod.color.split(", ")],
                            "mod_weights_name": weights_name,
                            "n_params": len(mod.parameters_guesses),
                            "res_real": (Z - _eval).real / np.abs(Z),
                            "res_imag": (Z - _eval).imag / np.abs(Z),
                            "MSE": (Z - _eval).imag ** 2 + (Z - _eval).imag ** 2,
                            "lmfit_message": result.message,
                            "params": paramDF,
                            "params_err": paramDF_err,
                            **chi_sq(Z, _eval, result.nvarys),
                        }
                    }
                )
            except Exception as e:
                print(f"fit fail: : {_key} with {mod.name}", e)

    # circ_C_W.fit(frequencies, Z)
    # circ_CPE_W.fit(frequencies, Z)
    # circ_RC1_CPE_W.fit(frequencies,Z)
    # circ_C_W_fit = circ_C_W.predict(f_pred)
    # circ_CPE_W_fit = circ_CPE_W.predict(f_pred)
    # circ_RC1_CPE_W_fit = circ_RC1_CPE_W.predict(f_pred)
    # Z_fit = circuit.predict(frequencies)
    # mod_legends = [f'{i} : {np.sum(test_fits[i]["MSE"]):.2f}, {np.sum(test_fits[i]["MSE"])/test_fits[i]["n_params"]:.2Ef}, {test_fits[i]["redchi"]:.3G}' for i in list(test_fits.keys())]

    test_fits = OrderedDict(sorted(test_fits.items(), key=lambda x: x[1]["redchi"]))

    n_tf = len(test_fits.keys())

    test_fits = {
        key: {**val, **{"n_alpha": np.max([0.3, 1 - n / 10])}}
        for n, (key, val) in enumerate(test_fits.items())
    }

    mod_legends = [
        f'{i} : {test_fits[i]["chisqr"]:.2f}, {test_fits[i]["redchi"]:.3G}'
        for i in list(test_fits.keys())
    ]

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(10, 16)

    axZ_small = fig.add_subplot(gs[:7, :5])
    ax1 = fig.add_subplot(gs[:7, 5:11])
    axY = fig.add_subplot(gs[:7, 11:])
    ax_res = fig.add_subplot(gs[8:, :])

    # fig, ax = plt.subplots(figsize = (10,10))
    Z_small = Z[-20:]
    Z_45 = [i.real + (-1j * i.real + 1j * Z.real.min()) for i in Z]
    ax1.set_title(f"{_key} at {dt.datetime.now()}")
    plot_nyquist(ax1, Z, fmt="o")
    plot_nyquist(ax1, Z_45, fmt="--")

    plot_nyquist(axZ_small, Z_small, fmt="o")
    plot_nyquist(axZ_small, Z_45, fmt="--")

    plot_nyquist(axY, Z_to_Y(Z), fmt="o")
    plot_nyquist(axY, [i.real + (-1j * i.real * 0.5) for i in Z_to_Y(Z)], fmt="--")

    for key, val in test_fits.items():
        vc = val.get("mod_color", "")
        kw = {}
        if vc:
            kw = {"color": vc}
        plot_nyquist(ax1, val["predict"], fmt="-", alpha=val["n_alpha"], **kw)
        plot_nyquist(
            axZ_small, val["predict"][:25], fmt="-", alpha=val["n_alpha"], **kw
        )
        plot_nyquist(axY, Z_to_Y(val["predict"]), fmt="-", alpha=val["n_alpha"], **kw)
        print(val["mod_name"])
        plot_residuals(
            ax_res,
            frequencies,
            val["res_real"],
            val["res_imag"],
            fmt="-",
            y_limits=(-10, 10),
            extra_label=val["mod_name"],
            alpha=val["n_alpha"],
            **kw,
        )
    #    ax2.text(1,6,f'MSE:{np.sum(val["MSE"]):.2f}')
    ax1.set_ylim((0, abs(Z.imag).max() * 2))
    ax1.set_xlim((0, abs(Z.imag).max() * 2))
    axZ_small.set_ylim(0, Z_small.real.min() * 2)
    axZ_small.set_xlim(Z_small.real.min() * 0.5, Z_small.real.min() * 2)

    _lbls = ["Data"] + ["45 line"] + mod_legends
    ax1.legend(_lbls, ncol=4, fontsize=10, loc="upper left", bbox_to_anchor=(-1.3, 1.4))
    ax_res.legend(ncol=4)
    plt.savefig(
        _dest_dir.joinpath(_key + "_IMPY_nyquist").with_suffix(".png"),
        bbox_inches="tight",
        dpi=300,
    )
    #    plt.show()
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 10), nrows=2)
    ax[0].set_title(f"{_key} at {dt.datetime.now()}")
    plot_bode(ax, 2 * np.pi * frequencies, Z, fmt="o", mag_log=True)
    plot_bode(ax, 2 * np.pi * frequencies, Z_45, fmt="--", mag_log=True)
    for key, val in test_fits.items():
        vc = val.get("mod_color", "")
        kw = {}
        if vc:
            kw = {"color": vc}
        plot_bode(
            ax,
            2 * np.pi * f_pred,
            val["predict"],
            fmt="-",
            label=key,
            mag_log=True,
            alpha=val["n_alpha"],
            **kw,
        )
    #    print(val['mod'])
    # plt.legend(['Data']+['45 line']+mod_legends)
    plt.legend(
        _lbls, ncol=4, fontsize=10, loc="upper left", bbox_to_anchor=(-0.5, -0.30)
    )
    # 'C', 'CPE1','RC1'])
    plt.savefig(
        _dest_dir.joinpath(_key + "_IMPY_bode").with_suffix(".png"),
        bbox_inches="tight",
        dpi=300,
    )
    #    plt.show()
    plt.close()

    if test_fits:
        params_out = pd.concat(
            [
                pd.concat([val["params"], val["params_err"]], axis=1)
                .assign(
                    **{
                        "MSE": np.sum(val["MSE"]),
                        "n_params": val["n_params"],
                        "chisqr": val["chisqr"],
                        "redchi": val["redchi"],
                        "File": _key,
                        "mod_name": val["mod_name"],
                        "mod_circuit": val["mod_circuit"],
                    }
                )
                .set_index(["File", "mod_name"])
                for key, val in test_fits.items()
            ],
            sort=False,
        ).assign(**{"testfile": _tfile})
        params_out.to_excel(_dest_file)

        spectra_out = pd.concat(
            [
                pd.DataFrame(
                    {
                        "mod_data": val["mod_data"],
                        "mod_fit": val["mod_fit"],
                        "File": _key,
                        "mod_name": val["mod_name"],
                    }
                )
                for key, val in test_fits.items()
            ]
        ).assign(**{"testfile": _tfile})
        spectra_out.to_excel(_dest_spectra_file)
    else:
        params_out, spectra_out = pd.DataFrame(), pd.DataFrame()
    return params_out, spectra_out


def multi_run(all_test_data, run_models, dest_path, dest_path_complex, pre_params):
    _params, _spectra = [], []

    _filter = (
        '(pre_params.File.str.startswith("O2")) & (pre_params.File.str.contains("EIS")'
    )
    best_params = pre_params.loc[
        pre_params.groupby("mod_name")["chisqr"].transform("idxmin").unique()
    ]
    best_params = pd.DataFrame()
    _raw_test_data = {k: val for k, val in all_test_data.items() if "fit" in k}
    for _testfile, testval in _raw_test_data.items():
        # _testfile = _test_name_select[0]
        print("Loading: ", _testfile, "\n")
        _args = choose_test(
            all_test_data,
            name=_testfile,
            spec_type="spectrumfit",
            dest_path=dest_path,
            freqlim=15e3,
        )
        # _args = N_freqs, freq_vec, Z_exp, _key, spec
        params_out, spectra_out = testing_models(
            run_models,
            *_args,
            dest_path=dest_path_complex,
            best_params=best_params,
            _tfile=_testfile,
        )
        # .assign(**{'testfile' : _testfile})
        _params.append(params_out)
        _spectra.append(spectra_out)
    models_params = pd.concat(_params)
    models_spectra = pd.concat(_spectra)
    return models_params, models_spectra


def check_params():

    for par in ["R1", "R0", "C2", "CPE0", "CPE1"]:
        for _f, _fgrp in models_params.groupby("testfile"):
            _fgrp.dropna(subset=[par]).plot.bar(
                x="mod_circuit", y=par, yerr=f"{par}_err", rot=60, logy=0, title=_f
            )

    for par in ["R1", "R0", "G0_0", "G0_1", "C2", "CPE0", "CPE1"][2:4]:
        for _f, _fgrp in models_params.query("redchi < 2").groupby("mod_circuit"):
            if par in _fgrp.dropna(axis=1, how="all").columns:
                _fgrp.dropna(subset=[par]).plot.bar(
                    x="testfile", y=par, yerr=f"{par}_err", rot=60, logy=0, title=_f
                )

    _check_mod = [
        "R0-L0-p(C1,R2-Ws1,R3-CPE2)",
        "R0-p(C1,R2-Ws1,R3-C2)",
        "R0-L0-p(C1,R2-G1,R3-C2)",
        "R0-L0-p(R1-G0,CPE1)",
    ][0]
    _check_mod = best_mod_params.index[2]
    par_mod = models_params.query("mod_circuit == @_check_mod").reset_index()
    _pars = [
        i
        for i in par_mod.columns
        if i.split("_")[0] in _check_mod and not i.endswith("err")
    ]
    par_mod.plot.bar(x="File", y="redchi", rot=45, logy=True, title=_check_mod)

    for par in _pars:
        _parvals = par_mod.query(f"{par}_err < 1E3 & {par} < 1E7")[par]
        _ymax = (
            _parvals.max()
            if not _parvals.max() > _parvals.mean() * 5
            else _parvals.mean() * 4
        )
        par_mod.dropna(subset=[par]).plot.bar(
            x="File",
            y=par,
            yerr=f"{par}_err",
            rot=45,
            logy=0,
            title=_check_mod,
            ylim=(0, _ymax),
        )

    #    fig,ax = plt.subplots()
    models_params.groupby("File").plot.bar(x="mod_name", y="L0", rot=60)
    models_params.set_index(["File", "mod_name"]).plot.bar(
        y="CPE2_0", rot=40, logy=True
    )
    pass


def test_plot_models_collection():
    modc = Model_Collection()

    fig, ax = plt.subplots(figsize=(12, 12))
    for mod_inst in modc.lmfit_models:
        mc = [float(i) for i in mod_inst.color.split(", ")]
        ax.scatter(
            x=mod_inst.mod_eval.real,
            y=-1 * mod_inst.mod_eval.imag,
            color=mc,
            alpha=0.6,
            label=mod_inst.name,
        )

    ax.set_xlim(0, 160)
    ax.set_ylim(-10, 150)
    ax.grid(True)
    ax.legend(ncol=3, fontsize=12, loc="upper left", bbox_to_anchor=(-0.1, 1.2))

    plt.save
    # ax.legend(ncol=3)


def plot_params(models_params, models_spectra, dest_path_params):
    models_params
    for modn, mgrp in models_params.groupby("mod_name"):
        modn, mgrp
        # modn = 'L-TLM(Rct-Qad-W)'
        # mgrp = models_params.groupby('mod_name').get_group(modn)
        mgrp.to_excel(dest_path_params.joinpath(f"{modn}.xlsx"))
        varnames = [i.split("_stderr")[0] for i in mgrp.columns if "stderr" in i]
        mgrp = mgrp.dropna(how="all", axis=1)
        # mgrp = mgrp.dropna(subset=varnames , how='all',axis=1)
        # mgrp = mgrp.dropna(how='all',axis=1)
        varnames = [i for i in varnames if i in mgrp.columns] + ["chisqr", "redchi"]
        for var in varnames:
            if not mgrp[var].dropna().empty:
                ymax = mgrp[var].max()
                fig, ax = plt.subplots()
                mgrp.plot(x="testfile", y=var, ylim=[0, ymax], rot=90, title=modn)
                plt.savefig(
                    dest_path_params.joinpath(f"{modn}_{var}.png"),
                    bbox_inches="tight",
                    dpi=100,
                )
                plt.close()


def plot_mod_multi_spectra(modn, models_spectra, dest_path_params):

    for modn, mgrp in models_spectra.groupby("mod_name"):
        modn, mgrp
        fig, ax = plt.subplots(nrows=mgrp.testfile.nunique(), ncols=2, figsize=(20, 40))
        for n, (tf, tgrp) in enumerate(mgrp.groupby("testfile")):
            n, (tf, tgrp)
            axZ = ax[n][0]
            axY = ax[n][1]
            axZ.scatter(
                x=tgrp["mod_data"].to_numpy().real,
                y=-1 * tgrp["mod_data"].to_numpy().imag,
                c="r",
            )
            axZ.plot(
                tgrp["mod_fit"].to_numpy().real,
                -1 * tgrp["mod_fit"].to_numpy().imag,
                c="blue",
            )
            axZ.annotate(tf, xy=(-1, 0.5), xycoords="axes fraction")
            axZ.set_xlim([0, 100])
            axZ.set_ylim([0, 100])

            Y = Z_to_Y(tgrp["mod_data"].to_numpy())
            Y_fit = Z_to_Y(tgrp["mod_fit"].to_numpy())
            axY.scatter(x=Y.real, y=-1 * Y.imag, c="r")
            axY.plot(Y_fit.real, -1 * Y_fit.imag, c="blue")
            # axZ.annotate(tf, xy=(-1, 0.5), xycoords='axes fraction')
            x_max = 0.07 if Y.real.max() < 0.1 else 0.2
            y_max = 0.03 if np.abs(Y.imag).max() < 0.03 else 0.1
            axY.set_xlim([0, x_max])
            axY.set_ylim([0, y_max])

        plt.savefig(
            dest_path_params.joinpath(f"{modn}_spectra.png"),
            bbox_inches="tight",
            dpi=100,
        )
        plt.close()


if __name__ == "__main__":
    all_test_data
    pre_params = pd.concat(
        [
            pd.read_excel(a, index_col=[0, 1])
            for i in all_test_data.keys()
            for a in (dest_path.rglob(f"{i}*_IMPY_params_TYPES*"))
        ],
        sort=False,
        ignore_index=False,
    )
    # .set_index(['File','mod_name'])
    models_params = pre_params
    _type_models = type_models
    # mc = Model_Collection(_startswith='F_')
    modc = Model_Collection(_startswith="F_")
    run_models = modc.lmfit_models
    models_params, models_spectra = multi_run(
        all_test_data, run_models, dest_path, dest_path_complex, pre_params
    )
    _err = "chisqr"
    best_mod_params = (
        models_params.loc[models_params.redchi < 1e2]
        .groupby("mod_circuit")[[_err, "n_params"]]
        .agg(["sum", "mean", "count", "std"])
        .sort_values(by=[(_err, "mean")])
    )
    print(
        best_mod_params[
            [(_err, "sum"), (_err, "mean"), (_err, "count"), ("n_params", "mean")]
        ]
    )


# fig, ax = plt.subplots(figsize = (10,10),nrows=2)
# plot_residuals(ax, frequencies, Z, fmt='o')f
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
