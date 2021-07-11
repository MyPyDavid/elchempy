#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:00:48 2020

@author: zmg
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from pathlib import Path
import json

from scipy.optimize import minimize
from scipy.stats import linregress
import pandas as pd


# print(f'GP_DRT_FITTING: Name: {__name__} for file {__file__}')
from . import GP_DRT_Ciucci_Liu_2019 as GP_DRT


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
    xl_files = list(Path.cwd().rglob("testing_data/*xlsx"))
    #    spec = pd.read_excel(xl_files[1],index_col=[0])
    #    converters={'DATA_Z': lambda s: np.complex(s.replace('i', 'j'))}
    all_data = {
        a.stem: {"Filepath": a, "spectrum": reduce_Z_data(read_xl(a))}
        for a in xl_files
        if not "_GP_" in a.name
    }
    specs = [i["spectrum"] for i in all_data.values()]
    return all_data


def check_Warburg(spec):
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
# all_test_data = read_eis_excel()
def choose_test(
    all_test_data,
    name="O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm_spectrumfit_v20",  # pragma: allowlist secret
    reduce=False,
):

    jos2 = [
        i
        for i in list(all_test_data.keys())
        if "JOS2_288_758mV_1500rpm_3" in i and "raw" in i
    ]
    #    spec = all_test_data.get('O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm_spectrumfit_v20')['spectrum'] # pragma: allowlist secret
    _key = jos2[0]
    spec = all_test_data.get(_key)["spectrum"]
    check_Warburg(spec)
    N_freqs = len(spec)
    freq_vec = spec["Frequency(Hz)"].to_numpy()
    Z_exp = spec.DATA_Z.to_numpy()
    if reduce:
        Z_exp = spec.DATA_Z_reduce.values
    print(_key)
    return N_freqs, freq_vec, Z_exp, _key


def GP_DRT_analysis(DRT_data_KKvalid, fit_display=False):
    DRT_data_KKvalid = DRT_data_KKvalid.sort_values(by="Frequency(Hz)", ascending=True)
    #    N_freqs, freq_vec, Z_exp, _file =  choose_test(all_test_data,reduce=False)
    freq_vec, Z_exp = (
        DRT_data_KKvalid["Frequency(Hz)"].to_numpy(),
        DRT_data_KKvalid.DATA_Z.to_numpy(),
    )
    #    freq_vec, Z_exp = DRT_data_KKvalid['Frequency(Hz)'].to_numpy(), DRT_data_KKvalid.DATA_Z.to_numpy()
    #%%
    # define the frequency range
    N_freqs = len(freq_vec)
    xi_vec = np.log(freq_vec)
    tau = 1 / freq_vec

    # define the frequency range used for prediction, we choose a wider range to better display the DRT
    freq_vec_star = np.logspace(-4.0, 6.0, num=101, endpoint=True)
    xi_vec_star = np.log(freq_vec_star)

    # finer mesh for plotting only
    freq_vec_plot = np.logspace(-4.0, 6.0, num=1001, endpoint=True)

    # initial parameters parameter to maximize the marginal log-likelihood as shown in eq (31)
    sigma_n, sigma_f, ell = 1.0e-4, 1.0e-3, 1.0  # standard initial params
    sigma_n, sigma_f, ell = -1.4330427, -116.6128790, -3.4878355
    #    sigma_n = -7.3E-01
    #    sigma_f = -2.06E2
    #    ell = -3

    theta_0 = np.array([sigma_n, sigma_f, ell])
    seq_theta = np.copy(theta_0)

    def print_results(theta):
        global seq_theta
        seq_theta = np.vstack((seq_theta, theta))
        print("{0:.7f}  {1:.7f}  {2:.7f}".format(theta[0], theta[1], theta[2]))

    #    print('sigma_n,   sigma_f,   ell')
    #
    #    GP_DRT.NMLL_fct(theta_0, Z_exp, xi_vec)
    GP_DRT.grad_NMLL_fct(theta_0, Z_exp, xi_vec)

    # minimize the NMLL $L(\theta)$ w.r.t sigma_n, sigma_f, ell using the BFGS method as implemented in scipy
    res = minimize(
        GP_DRT.NMLL_fct,
        theta_0,
        args=(Z_exp, xi_vec),
        method="BFGS",
        jac=GP_DRT.grad_NMLL_fct,
        callback=None,
        options={"disp": fit_display},
    )

    # collect initial parameters
    res.init_sigma_n = theta_0[0]
    res.init_sigma_f = theta_0[1]
    res.init_ell = theta_0[2]

    # collect the optimized parameters
    sigma_n, sigma_f, ell = res.x
    res.sigma_n = sigma_n
    res.sigma_f = sigma_f
    res.ell = ell

    # calculate the matrices shown in eq (18)
    K = GP_DRT.matrix_K(xi_vec, xi_vec, sigma_f, ell)
    L_im_K = GP_DRT.matrix_L_im_K(xi_vec, xi_vec, sigma_f, ell)
    L2_im_K = GP_DRT.matrix_L2_im_K(xi_vec, xi_vec, sigma_f, ell)
    Sigma = (sigma_n ** 2) * np.eye(N_freqs)

    # the matrix $\mathcal L^2_{\rm im} \mathbf K + \sigma_n^2 \mathbf I$ whose inverse is needed
    K_im_full = L2_im_K + Sigma

    # Cholesky factorization, L is a lower-triangular matrix
    L = np.linalg.cholesky(K_im_full)

    # solve for alpha
    alpha = np.linalg.solve(L, Z_exp.imag)
    alpha = np.linalg.solve(L.T, alpha)

    # estimate the gamma of eq (21a), the minus sign, which is not included in L_im_K, refers to eq (65)
    gamma_fct_est = -np.dot(L_im_K.T, alpha)

    # covariance matrix
    inv_L = np.linalg.inv(L)
    inv_K_im_full = np.dot(inv_L.T, inv_L)

    # estimate the sigma of gamma for eq (21b)
    cov_gamma_fct_est = K - np.dot(L_im_K.T, np.dot(inv_K_im_full, L_im_K))
    sigma_gamma_fct_est = np.sqrt(np.diag(cov_gamma_fct_est))

    # initialize the imaginary part of impedance vector
    Z_im_vec_star = np.empty_like(xi_vec_star)
    Sigma_Z_im_vec_star = np.empty_like(xi_vec_star)

    gamma_vec_star = np.empty_like(xi_vec_star)
    Sigma_gamma_vec_star = np.empty_like(xi_vec_star)

    # calculate the imaginary part of impedance at each $\xi$ point for the plot
    for index, val in enumerate(xi_vec_star):
        xi_star = np.array([val])

        # compute matrices shown in eq (18), k_star corresponds to a new point
        k_star = GP_DRT.matrix_K(xi_vec, xi_star, sigma_f, ell)
        L_im_k_star = GP_DRT.matrix_L_im_K(xi_vec, xi_star, sigma_f, ell)
        L2_im_k_star = GP_DRT.matrix_L2_im_K(xi_vec, xi_star, sigma_f, ell)
        k_star_star = GP_DRT.matrix_K(xi_star, xi_star, sigma_f, ell)
        L_im_k_star_star = GP_DRT.matrix_L_im_K(xi_star, xi_star, sigma_f, ell)
        L2_im_k_star_star = GP_DRT.matrix_L2_im_K(xi_star, xi_star, sigma_f, ell)

        # compute Z_im_star mean and standard deviation using eq (26)
        Z_im_vec_star[index] = np.dot(L2_im_k_star.T, np.dot(inv_K_im_full, Z_exp.imag))
        Sigma_Z_im_vec_star[index] = L2_im_k_star_star - np.dot(
            L2_im_k_star.T, np.dot(inv_K_im_full, L2_im_k_star)
        )

        # compute Z_im_star mean and standard deviation
        gamma_vec_star[index] = -np.dot(
            L_im_k_star.T, np.dot(inv_K_im_full, Z_exp.imag)
        )
        Sigma_gamma_vec_star[index] = k_star_star - np.dot(
            L_im_k_star.T, np.dot(inv_K_im_full, L_im_k_star)
        )

    #    _Z_star.plot(x='freq_vec_star',y='gamma_vec_star',logx=True)
    # collect the optimized paramete
    _Z_exp = pd.DataFrame(data=zip(Z_exp, freq_vec), columns=["Z_exp", "freq_vec"])
    _Z_exp = _Z_exp.assign(
        **{
            "Z_exp_real": _Z_exp["Z_exp"].to_numpy().real,
            "Z_exp_imag": _Z_exp["Z_exp"].to_numpy().imag,
        }
    )
    _Z_star = pd.DataFrame(
        zip(freq_vec_star, gamma_vec_star, Sigma_gamma_vec_star, Z_im_vec_star),
        columns=[
            "freq_vec_star",
            "gamma_vec_star",
            "Sigma_gamma_vec_star",
            "Z_im_vec_star",
        ],
    )
    res_fit_params = {
        key: val for key, val in dict(res).items() if "array" not in str(type(val))
    }
    res_fit_arrays = {
        key: val for key, val in dict(res).items() if "array" in str(type(val))
    }

    #%%
    return _Z_exp, _Z_star, res_fit_params, res_fit_arrays


def run_GP_DRT_fit(DRT_data_KKvalid, **kwargs):
    _Z_exp, _Z_star, res_fit_params, res_fit_arrays = GP_DRT_analysis(DRT_data_KKvalid)
    savefig = kwargs.get("GP_DRT_savefig_path", False)
    if savefig:
        _ny_path = plot_nyquist(_Z_exp.Z_exp, _Z_exp.freq_vec, savefig=savefig)
        _GP_DRT_path = plot_GP_DRT(
            _Z_star.freq_vec_star,
            _Z_star.gamma_vec_star,
            _Z_star.Sigma_gamma_vec_star,
            savefig=savefig,
        )
        _GP_DRT_imag_path = plot_imag(
            _Z_exp.freq_vec,
            _Z_exp.Z_exp.to_numpy(),
            _Z_star.freq_vec_star,
            _Z_star.Z_im_vec_star,
            _Z_star.Sigma_gamma_vec_star,
            savefig=savefig,
        )
        res_fit_params.update(
            {
                "GP_DRT_plot_nyquist": str(_ny_path),
                "GP_DRT_plot_gamma": str(_GP_DRT_path),
                "GP_DRT_plot_fit": str(_GP_DRT_imag_path),
            }
        )
        #        if kwargs.get('GP_DRT_pickle_path',False):
        _pkl_path_Z_exp = savefig.with_name(savefig.name + "_GP_Z_exp").with_suffix(
            ".pkl"
        )
        _Z_exp.to_pickle(_pkl_path_Z_exp)
        _pkl_path_Z_star = savefig.with_name(
            savefig.name + "_GP_DRT_Z_star"
        ).with_suffix(".pkl")
        _Z_star.to_pickle(_pkl_path_Z_star)
        res_fit_params.update(
            {
                "GP_DRT_DF_Z_exp": str(_pkl_path_Z_exp),
                "GP_DRT_DF_Z_star": str(_pkl_path_Z_star),
            }
        )

        _json_path_res = savefig.with_name(
            savefig.name + "_GP_DRT_res_params"
        ).with_suffix(".json")
        #        res_fp_json = json.dumps(res_fit_params)
        with open(_json_path_res, "w") as outfile:
            json.dump(res_fit_params, outfile)
        res_fit_params.update({"GP_DRT_json_res": str(_json_path_res)})

    return _Z_exp, _Z_star, res_fit_params, res_fit_arrays


def plot_nyquist(Z_exp, freq_vec, savefig=""):
    N_freqs = len(freq_vec)
    # Nyquist plot of the EIS spectrum
    plt.plot(
        np.real(Z_exp),
        -np.imag(Z_exp),
        "o",
        markersize=10,
        fillstyle="none",
        color="red",
        label="experiment",
    )
    plt.plot(
        np.real(Z_exp[1:N_freqs:10]),
        -np.imag(Z_exp[1:N_freqs:10]),
        "o",
        markersize=10,
        color="black",
    )

    plt.rc("text")
    plt.rc("font", family="serif", size=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.legend(frameon=False, fontsize=15)
    plt.axis("scaled")

    # this depends on the data used - if you wish to use your own data you may need to modify this
    #    plt.xlim(1.42, 1.52)
    #    plt.ylim(-0.001, 0.051)
    #    plt.xticks(np.arange(1.42, 1.521, 0.02))
    #    plt.yticks(np.arange(0.00, 0.051, 0.01))
    plt.autoscale(True)
    #    plt.yscale('log')
    #    plt.xscale('log')
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel(r"$Z_{\rm re}/\Omega$", fontsize=20)
    plt.ylabel(r"$-Z_{\rm im}/\Omega$", fontsize=20)

    # label the frequencies - if you wish to use your own data you may need to modify this
    label_index = range(1, N_freqs, 10)
    move = [[-0.005, 0.008], [-0.005, 0.008], [-0.005, 0.008], [-0.005, 0.01]]
    for k, ind in enumerate(label_index):
        power = int(np.log10(freq_vec[ind]))
        num = freq_vec[ind] / (10 ** (power))
        plt.annotate(
            f"{num:.1f}E{power:}",
            xy=(np.real(Z_exp[ind]), -np.imag(Z_exp[ind])),
            xytext=(np.real(Z_exp[ind]) * 1.005, 1.15 * -np.imag(Z_exp[ind])),
            arrowprops=dict(arrowstyle="-", connectionstyle="arc"),
        )
    if savefig:
        _savepath = savefig.with_name(savefig.name + "_GP_nyquist").with_suffix(".png")
        plt.savefig(_savepath, dpi=200, bbox_inches="tight")
        plt.close()
        return _savepath

    else:
        plt.show()
        plt.close()
        return ""


def plot_GP_DRT(freq_vec_star, gamma_vec_star, Sigma_gamma_vec_star, savefig=""):
    # plot the DRT and its confidence region
    plt.semilogx(
        freq_vec_star, gamma_vec_star, linewidth=4, color="red", label="GP-DRT"
    )
    plt.fill_between(
        freq_vec_star,
        gamma_vec_star - 3 * np.sqrt(abs(Sigma_gamma_vec_star)),
        gamma_vec_star + 3 * np.sqrt(abs(Sigma_gamma_vec_star)),
        color="0.4",
        alpha=0.3,
    )
    plt.rc("text", usetex=False)
    plt.rc("font", family="serif", size=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    #    plt.axis([1E-4,1E6,-0.01,50])
    plt.autoscale(True)
    #    plt.yticks(np.arange(-0.01, 0.025, 0.01))
    plt.legend(frameon=False, fontsize=15)
    plt.xlabel(r"$f/{\rm Hz}$", fontsize=20)
    plt.ylabel(r"$\gamma/\Omega$", fontsize=20)

    if savefig:
        _savepath = savefig.with_name(savefig.name + "_GP_DRT").with_suffix(".png")
        plt.savefig(_savepath, dpi=200, bbox_inches="tight")
        plt.close()
        return _savepath
    else:
        plt.show()
        plt.close()
        return ""


def plot_imag(
    freq_vec, Z_exp, freq_vec_star, Z_im_vec_star, Sigma_Z_im_vec_star, savefig=""
):
    plt.semilogx(freq_vec, -Z_exp.imag, "o", markersize=10, color="black", label="exp")
    plt.semilogx(
        freq_vec_star, -Z_im_vec_star, linewidth=4, color="red", label="GP-DRT"
    )
    plt.fill_between(
        freq_vec_star,
        -Z_im_vec_star - 3 * np.sqrt(abs(Sigma_Z_im_vec_star)),
        -Z_im_vec_star + 3 * np.sqrt(abs(Sigma_Z_im_vec_star)),
        alpha=0.3,
    )
    plt.rc("text", usetex=False)
    plt.rc("font", family="serif", size=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.autoscale(True)
    #    plt.axis([1E-3,1E5,Z_im_vec_star.min(),Z_im_vec_star.max()])
    plt.legend(frameon=False, fontsize=15)
    plt.xlabel(r"$f/{\rm Hz}$", fontsize=20)
    plt.ylabel(r"$-Z_{\rm im}/\Omega$", fontsize=20)

    if savefig:
        _savepath = savefig.with_name(savefig.name + "_GP_imag").with_suffix(".png")
        plt.savefig(_savepath, dpi=200, bbox_inches="tight")
        plt.close()
        return _savepath
    else:
        plt.show()
        plt.close()
        return ""
