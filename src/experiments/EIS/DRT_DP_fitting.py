#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:53:16 2020

@author: zmg

In this tutorial we will reproduce Figure 2 in
Liu, J., & Ciucci, F. (2020). The Deep-Prior Distribution of Relaxation Times. Journal of The Electrochemical Society, 167(2), 026506
https://iopscience.iop.org/article/10.1149/1945-7111/ab631a/meta

The DP-DRT method is our next newly developed deep learning based approach to obtain the DRT from the EIS data.
The DP-DRT is trained on a single electrochemical impedance spectrum.
 A single random input is given to the neural network underlying the DP-DRT.
"""


import numpy as np
import os
from pathlib import Path
from collections import namedtuple
import matplotlib.pyplot as plt
import random as rnd
import math
from math import sin, cos, pi
import pandas as pd

from warnings import warn

try:

    import torch
    import torch.nn.functional as F

    # check the device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
except ModuleNotFoundError:
    warn("Torch not found")


if __name__ == "__main__":

    import compute_DP_DRT_Ciucci_2020 as compute_DRT
else:
    try:
        import compute_DP_DRT_Ciucci_2020 as compute_DRT
    except:
        from . import compute_DP_DRT_Ciucci_2020 as compute_DRT


# print('Using device:', device)
# if device.type == 'cuda':
#    print(torch.cuda.get_device_name(0))
#    print('Memory Usage:')
#    print('Allocated:', round(torch.cuda.memory_allocated(0)/1024**2,1), 'MB')
#    print('Cached:   ', round(torch.cuda.memory_cached(0)/1024**2,1), 'MB')
#
# we will assume you have a cpu
# if you want to use a GPU, you will need to use cuda


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
    except Exception as e:
        print(e)

    return spec


def read_xl(xlfile):
    df = pd.read_excel(xlfile, index_col=[0])

    if "Frequency(Hz)" in df.columns:
        df = df.sort_values("Frequency(Hz)", ascending=True)
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
        a.stem: {"Filepath": a, "spectrum": reduce_Z_data(read_xl(a))} for a in xl_files
    }
    specs = [i["spectrum"] for i in all_data.values()]
    return all_data


# all_test_data = read_eis_excel()
def choose_test(
    all_test_data,
    name="O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm_spectrumfit_v20",  # pragma: allowlist secret
    reduce=False,
):
    all_test_data = read_eis_excel()
    spec = all_test_data.get(
        "O2_EIS-range_1500rpm_JOS2_899_499mV_1500rpm_spectrumfit_v20"  # pragma: allowlist secret
    )["spectrum"]
    N_freqs = len(spec)
    Z_exp = spec.DATA_Z.values
    if reduce:
        Z_exp = spec.DATA_Z_reduce.values
    return N_freqs, Z_exp * 0.238


def standard_DRT(N_freqs=81):
    #    N_freqs = 81
    freq_vec = np.logspace(-3.0, 3.0, num=N_freqs, endpoint=True)
    tau_vec = 1.0 / freq_vec
    # define parameters for ZARC model and calculate the impedance and gamma following the above equations
    R_inf = 10
    R_ct = 10
    phi = 0.8
    tau_0 = 0.1
    C = tau_0 ** phi / R_ct
    tau_2, R_2 = 5, 20
    C_2 = tau_2 ** phi / R_2
    C_3 = 1e-06
    # exact Z and gamma
    Z = (
        R_inf
        + 1.0 / (1.0 / R_ct + C * (1j * 2.0 * pi * freq_vec) ** phi)
        + 1.0 / (1.0 / R_2 + C_2 * (1j * 2.0 * pi * freq_vec) ** 0.8)
    )
    gamma_exact = (
        (R_ct)
        / (2.0 * pi)
        * sin((1.0 - phi) * pi)
        / (np.cosh(phi * np.log(tau_vec / tau_0)) - cos((1.0 - phi) * pi))
    )

    # adding noise to the impedance data
    sigma_n_exp = 0.01
    Z_exp = Z + sigma_n_exp * (
        np.random.normal(0, 1, N_freqs) + 1j * np.random.normal(0, 1, N_freqs)
    )
    #    plt.plot(np.real(Z_exp), -np.imag(Z_exp), "o", markersize=10, color="black", label="synth exp")
    return Z_exp, gamma_exact


#    N_freqs, Z_exp =  choose_test(all_test_data,reduce=False)
def DP_DRT_analysis(Z_exp):
    #%%
    # set the seed for the random number generators
    rng = rnd.seed(214975)
    rng_np = np.random.seed(213912)
    torch.manual_seed(213912)

    # DW PY ADDED

    N_freqs = len(Z_exp)
    Z_exact, gamma_exact = standard_DRT(N_freqs)
    # define frequency range, from 1E-4 to 1E4 with 10 ppd
    freq_vec = np.logspace(-4, 4.0, num=N_freqs, endpoint=True)
    tau_vec = 1.0 / freq_vec
    omega_vec = 2.0 * pi * freq_vec
    #
    # define the matrices that calculate the impedace from DRT,
    #  i.e., Z_re = A_re * gamma,
    # Z_im = A_im * gamma
    A_re = compute_DRT.A_re(freq_vec)
    A_im = compute_DRT.A_im(freq_vec)

    # transform impedance variables to tensors
    Z_exp_re_torch = (
        torch.from_numpy(np.real(Z_exp)).type(torch.FloatTensor).reshape(1, N_freqs)
    )
    Z_exp_im_torch = (
        torch.from_numpy(np.imag(Z_exp)).type(torch.FloatTensor).reshape(1, N_freqs)
    )
    # tranform gamma
    gamma_exact_torch = torch.from_numpy(gamma_exact).type(torch.FloatTensor)

    freq_vec_torch = torch.from_numpy(freq_vec).type(torch.FloatTensor)
    #    omega_vec_torch= torch.from_numpy(2.*pi*np.logspace(-3,4., num=
    #                               compute_DRT.count_parameters(model), endpoint=True)).type(torch.FloatTensor)
    omega_vec_torch = torch.from_numpy(omega_vec).type(torch.FloatTensor)
    tau_vec_torch = torch.from_numpy(tau_vec).type(torch.FloatTensor)

    # transform these matrices into tensors
    A_re_torch = torch.from_numpy(A_re.T).type(torch.FloatTensor)
    A_im_torch = torch.from_numpy(A_im.T).type(torch.FloatTensor)

    # size of the arbitrary zeta input
    N_zeta = 1

    # define the neural network
    # N is batch size, D_in is input dimension, H is hidden dimension, D_out is output dimension.
    N = 1
    D_in = N_zeta
    H = max(N_freqs + 2, 10 * N_zeta)
    # the output also includes the R_inf, so it has dimension N_freq+1
    # note that
    # 1) there is no inductance (in this specific example - the DP-DRT can include inductive features, see article)
    # 2) R_inf is stored as the last item in the NN output

    D_out = N_freqs + 2
    # Construct the neural network structure
    class vanilla_model(torch.nn.Module):
        def __init__(self):
            super(vanilla_model, self).__init__()
            self.fct_1 = torch.nn.Linear(D_in, H)
            self.fct_2 = torch.nn.Linear(H, H)
            self.fct_3 = torch.nn.Linear(H, H)
            self.fct_4 = torch.nn.Linear(H, D_out)

            # initialize the weight parameters
            torch.nn.init.zeros_(self.fct_1.weight)
            torch.nn.init.zeros_(self.fct_2.weight)
            torch.nn.init.zeros_(self.fct_3.weight)
            torch.nn.init.zeros_(self.fct_4.weight)

        # forward
        def forward(self, zeta):
            h = F.elu(self.fct_1(zeta))
            h = F.elu(self.fct_2(h))
            h = F.elu(self.fct_3(h))
            gamma_pred = F.softplus(self.fct_4(h), beta=5)
            return gamma_pred

        # transform impedance variables to tensors

    Z_exp_re_torch = (
        torch.from_numpy(np.real(Z_exp)).type(torch.FloatTensor).reshape(1, N_freqs)
    )
    Z_exp_im_torch = (
        torch.from_numpy(np.imag(Z_exp)).type(torch.FloatTensor).reshape(1, N_freqs)
    )
    # tranform gamma
    gamma_exact_torch = torch.from_numpy(gamma_exact).type(torch.FloatTensor)
    # transform these matrices into tensors
    A_re_torch = torch.from_numpy(A_re.T).type(torch.FloatTensor)
    A_im_torch = torch.from_numpy(A_im.T).type(torch.FloatTensor)

    # def DP_DRT_model():
    def loss_fn(
        output,
        Z_exp_re_torch,
        Z_exp_im_torch,
        A_re_torch,
        A_im_torch,
        _lambda,
        omega_vec_torch,
    ):
        #    output = gamma
        # we assume no inductance and the R_inf is stored as the last item in the NN output
        _gamma_only = output[:, 0:-2]
        _R_inf, _L_inf, _C_ser = output[:, -2], output[:, -1], []

        MSE_re = torch.sum(
            (Z_exp_re_torch - torch.mm(_gamma_only, A_re_torch) - _R_inf) ** 2
        )
        MSE_im = torch.sum(
            (
                (
                    Z_exp_im_torch
                    - torch.mm(_gamma_only, A_re_torch)
                    - _L_inf * omega_vec_torch
                )
            )
            ** 2
        )

        #        MSE_re = torch.sum((output[:, -3] + torch.mm(output[:, 0:-3], A_re_torch) - Z_exp_re_torch)**2)
        #        MSE_im = torch.sum((((output[:, -2]*omega_vec_torch)**1)**-1 + (output[:, -1]*omega_vec_torch) + torch.mm(output[:, 0:-3], A_im_torch) - Z_exp_im_torch)**2)
        #    reg = torch.sum(_lambda*(torch.mm(output[:, 0:-1], A_im_torch) )/(torch.mm(output[:, 0:-1], A_im_torch) ))
        #    C_sub = 1/(w_min*Zim_min)
        #     (1j*EIS_data_KKvalid.Angular*C_sub)**-1

        #        plt.semilogy(g_diff.detach().reshape(-1).numpy() , linewidth=4, color="black")
        #        plt.semilogy(g_diff_diff.detach().reshape(-1).numpy() , linewidth=4, color="black")
        #        plt.semilogy(tau_vec, gamma_only.detach().reshape(-1).numpy() , linewidth=4, color="black")
        #    filter = torch.nn.Conv1d(in_channels=1, out_channels=1, kernel_size=2, stride=1, padding=1, groups=1, bias=False)
        #    MSE = 0.2*MSE_re + MSE_im +  _lambda*torch.sum(torch.from_numpy(1*(np.abs(np.diff(np.diff(output[:, 0:-1].detach().reshape(-1)))))))
        if _lambda:
            g_diff = _gamma_only - F.pad(_gamma_only, (1, 0))[:, :-1]
            #           g_diff = gamma_only[:,1:] - gamma_only[:,:-1]
            g_diff_diff = g_diff - F.pad(g_diff, (1, 0))[:, :-1]
            #        g_diff[:,1:] - g_diff[:,:-1]
            d1_tau = tau_vec_torch - F.pad(tau_vec_torch, (1, 0))[:-1]
            d2_tau = d1_tau - F.pad(d1_tau, (1, 0))[:-1]
            #        tau_vec_torch[1:] - tau_vec_torch[:-1]
            #        d2_tau = d_tau[1:] -  d_tau[:-1]
            #        torch.sum(_lambda*(g_diff_diff/d2_tau)**2*d1_tau.abs())
            reg0_term = torch.sum(_lambda * _gamma_only)
            reg1_term = torch.sum(_lambda * (g_diff / d1_tau) ** 2 * d1_tau.abs())
            reg2_term = torch.sum(_lambda * (g_diff_diff / d2_tau) ** 2 * d1_tau.abs())
            reg_term = reg0_term + reg1_term + reg2_term
            MSE = MSE_re + MSE_im + reg_term

        #            _lambda*torch.sum(g_diff_diff)**2
        #            torch.sum(torch.mm(output,_lambda * torch.ones(output.size()[1],output.size()[1])))
        #            torch.mm(_lambda * torch.ones(output.size()),output.t())**2
        else:
            MSE = MSE_re + MSE_im
            reg_term = np.array(_lambda)
        return MSE, reg_term

    model = vanilla_model()
    # initialize following variables
    zeta = torch.randn(N, N_zeta)
    loss_vec = np.array([])
    reg_vec = np.array([])
    distance_vec = np.array([])
    lambda_vec = np.array([])

    # optimize the neural network
    learning_rate = 1e-5
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    _lambda = 0
    # max iteration
    max_iters = int(100e3)
    gamma_NN_store = torch.zeros((max_iters, N_freqs))
    R_inf_NN_store = torch.zeros((max_iters, 1))
    C_const_NN_store = torch.zeros((max_iters, 1))
    L_ser_NN_store = torch.zeros((max_iters, 1))

    for t in range(max_iters):
        # Forward pass: compute predicted y by passing x to the model.
        gamma = model(zeta)

        # Compute the loss
        loss, reg = loss_fn(
            gamma,
            Z_exp_re_torch,
            Z_exp_im_torch,
            A_re_torch,
            A_im_torch,
            _lambda,
            omega_vec_torch,
        )
        # save it
        loss_vec = np.append(loss_vec, loss.item())
        reg_vec = np.append(reg_vec, reg.item())

        #         _gamma_only = output[:, 0:-2]
        #        _R_inf,_L_inf, _C_ser = output[:, -2], output[:, -1], []
        # store gamma
        gamma_NN = gamma[:, 0:-2].detach().reshape(-1)
        gamma_NN_store[t, :] = gamma[:, 0:-2].detach().reshape(-1)
        # store R_inf
        R_inf_NN_store[t, :] = gamma[:, -2].detach().reshape(-1)
        # store C_ser
        #        C_const_NN_store[t,:] = gamma[:, -2].detach().reshape(-1)
        # store L_ser
        L_ser_NN_store[t, :] = gamma[:, -1].detach().reshape(-1)

        # Compute the distance
        distance = math.sqrt(torch.sum((gamma_NN - gamma_exact_torch) ** 2).item())
        # save it
        distance_vec = np.append(distance_vec, distance)
        #        '; distance_NOT USED=', distance
        # and print it
        if not t % 10000:
            print("iter=", t, "; loss=", loss.item())

        # zero all gradients (purge any cache)
        optimizer.zero_grad()

        # compute the gradient of the loss with respect to model parameters
        loss.backward()

        # Update the optimizer
        optimizer.step()

    #    distance_vec,loss_vec,gamma_NN_store,R_inf_NN_store = train_model()

    index_opt = np.argmin(loss_vec)
    index_early_stop = np.flatnonzero(np.abs(np.diff(loss_vec)) < 1e-8)

    C_const_vec = C_const_NN_store.detach().numpy()
    R_inf_vec = R_inf_NN_store.detach().numpy()
    L_ser_vec = L_ser_NN_store.detach().numpy()

    gamma_DIP_torch_opt = gamma_NN_store[index_opt, :]
    R_inf_DIP_torch_opt = R_inf_NN_store[index_opt, :]
    C_const_DIP_torch_opt = C_const_NN_store[index_opt, :]
    L_ser_DIP_torch_opt = L_ser_NN_store[index_opt, :]

    gamma_DIP_opt = gamma_DIP_torch_opt.detach().numpy()
    R_DIP_opt = R_inf_DIP_torch_opt.detach().numpy()
    C_DIP_opt = C_const_DIP_torch_opt.detach().numpy()
    L_DIP_opt = L_ser_DIP_torch_opt.detach().numpy()
    #    C_const_NN_store = torch.zeros((max_iters, 1))

    if len(index_early_stop):
        gamma_DIP_torch_early_stop = gamma_NN_store[index_early_stop[0], :]
        gamma_DIP = gamma_DIP_torch_early_stop.detach().numpy()
        R_DIP_torch_early = R_inf_NN_store[index_early_stop[0], :]
        R_DIP = R_DIP_torch_early.detach().numpy()
        C_DIP_torch_early = C_const_NN_store[index_early_stop[0], :]
        C_DIP = C_DIP_torch_early.detach().numpy()
        L_DIP_torch_early = L_ser_NN_store[index_early_stop[0], :]
        L_DIP = L_DIP_torch_early.detach().numpy()
    #        print('distance_early_stop = ', distance_vec[index_early_stop[0]])
    else:
        gamma_DIP = gamma_DIP_opt
        R_DIP = R_DIP_opt

    #    Z_DIP = R_DIP +L_DIP + np.matmul(A_re, gamma_DIP) + 1j*np.matmul(A_im, gamma_DIP)

    Z_DIP = (
        R_DIP
        + np.matmul(A_re, gamma_DIP)
        + 1j * np.matmul(A_im, gamma_DIP)
        + np.divide(
            1,
            (1j * omega_vec * C_DIP),
            out=np.zeros_like((1j * omega_vec * C_DIP)),
            where=(1j * omega_vec * C_DIP) != 0,
        )
        + (1j * omega_vec * L_DIP)
    )

    print("total number parameters = ", compute_DRT.count_parameters(model))
    DIP_out_templ = namedtuple(
        "DP_DRT_info",
        "Z_exp Z_DIP gamma_DIP freq_vec_DIP tau_vec_DIP loss_vec_DIP N_freqs_DIP R_DIP  total_pars",
    )
    DIP_out = DIP_out_templ(
        Z_exp,
        Z_DIP,
        gamma_DIP,
        freq_vec,
        tau_vec,
        loss_vec,
        N_freqs,
        R_DIP,
        compute_DRT.count_parameters(model),
    )
    #%%
    return DIP_out


def DP_DRT_plots(Z_exp, DIP_out, savepath=""):
    #    plt.semilogy(L_ser_vec, linewidth=4, color="black")
    #    plt.semilogy(freg_vec, linewidth=4, color="black")
    #    plt.semilogy(C_const_vec, linewidth=4, color="black")
    #    plt.semilogy(R_inf_vec, linewidth=4, color="black")
    #    plot_loss(loss_vec)
    if savepath:
        _save = savepath
    else:
        _save = None
    plot_imp(Z_exp, DIP_out.Z_DIP)
    plot_DRT(DIP_out.tau_vec_DIP, DIP_out.gamma_DIP)


def plot_test(DP_DRT_out):
    Z_DIP, gamma_DIP, freq_vec, tau_vec, loss_vec, N_freqs, R_DIP = DP_DRT_out[:-1]


def plot_loss(loss_vec, save=None):
    plt.semilogy(loss_vec, linewidth=4, color="black")
    #    plt.semilogy(np.array([index_early_stop[0], index_early_stop[0]]), np.array([1E-3, 1E7]),
    #                  ':', linewidth=3, color="red")
    #    plt.semilogy(np.array([index_opt, index_opt]), np.array([1E-3, 1E7]),
    #                  ':', linewidth=3, color="blue")
    plt.text(
        30000,
        1e2,
        r"early stop",
        {
            "color": "red",
            "fontsize": 20,
            "ha": "center",
            "va": "center",
            "rotation": 90,
            "bbox": dict(boxstyle="round", fc="white", ec="red", pad=0.2),
        },
    )
    plt.text(
        0.93e5,
        1e2,
        r"optimal",
        {
            "color": "blue",
            "fontsize": 20,
            "ha": "center",
            "va": "center",
            "rotation": 90,
            "bbox": dict(boxstyle="round", fc="white", ec="blue", pad=0.2),
        },
    )

    plt.rc("font", family="serif", size=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.xlabel(r"iter", fontsize=20)
    plt.ylabel(r"loss", fontsize=20)
    plt.axis([0, 1.01e5, 9, 1.1e7])
    fig = plt.gcf()
    fig.set_size_inches(5, 4)
    if save:
        plt.savefig(save, dpi=100, bbox_inches="tight")
    plt.show()
    plt.close()


def plot_R_C_L(loss_vec, index_early_stop, index_opt, save=None):
    plt.semilogy(loss_vec, linewidth=4, color="black")
    plt.semilogy(
        np.array([index_early_stop[0], index_early_stop[0]]),
        np.array([1e-3, 1e7]),
        ":",
        linewidth=3,
        color="red",
    )
    plt.semilogy(
        np.array([index_opt, index_opt]),
        np.array([1e-3, 1e7]),
        ":",
        linewidth=3,
        color="blue",
    )
    plt.text(
        30000,
        1e2,
        r"early stop",
        {
            "color": "red",
            "fontsize": 20,
            "ha": "center",
            "va": "center",
            "rotation": 90,
            "bbox": dict(boxstyle="round", fc="white", ec="red", pad=0.2),
        },
    )
    plt.text(
        0.93e5,
        1e2,
        r"optimal",
        {
            "color": "blue",
            "fontsize": 20,
            "ha": "center",
            "va": "center",
            "rotation": 90,
            "bbox": dict(boxstyle="round", fc="white", ec="blue", pad=0.2),
        },
    )

    plt.rc("font", family="serif", size=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.xlabel(r"iter", fontsize=20)
    plt.ylabel(r"loss", fontsize=20)
    plt.axis([0, 1.01e5, 9, 1.1e7])
    fig = plt.gcf()
    fig.set_size_inches(5, 4)
    if save:
        plt.savefig(save, dpi=100, bbox_inches="tight")
    plt.show()
    plt.close()


def plot_err_iter(distance_vec, index_early_stop, index_opt, save=None):
    plt.semilogy(distance_vec, linewidth=4, color="black")
    plt.semilogy(
        np.array([index_early_stop[0], index_early_stop[0]]),
        np.array([1e-3, 1e7]),
        ":",
        linewidth=4,
        color="red",
    )
    plt.semilogy(
        np.array([index_opt, index_opt]),
        np.array([1e-3, 1e7]),
        ":",
        linewidth=4,
        color="blue",
    )
    plt.text(
        30000,
        2e1,
        r"early stop",
        {
            "color": "red",
            "fontsize": 20,
            "ha": "center",
            "va": "center",
            "rotation": 90,
            "bbox": dict(boxstyle="round", fc="white", ec="red", pad=0.2),
        },
    )
    plt.text(
        0.93e5,
        2e1,
        r"optimal",
        {
            "color": "blue",
            "fontsize": 20,
            "ha": "center",
            "va": "center",
            "rotation": 90,
            "bbox": dict(boxstyle="round", fc="white", ec="blue", pad=0.2),
        },
    )
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif", size=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.xlabel(r"iter", fontsize=20)
    plt.ylabel(r"error", fontsize=20)
    plt.axis([0, 1.01e5, 5, 2.1e2])
    fig = plt.gcf()
    fig.set_size_inches(5, 4)
    if save:
        plt.savefig(save, dpi=100, bbox_inches="tight")
    plt.show()
    plt.close()


def plot_imp(Z_exp, Z_DIP, save=None):
    plt.plot(
        np.real(Z_exp),
        -np.imag(Z_exp),
        "o",
        markersize=10,
        color="black",
        label="synth exp",
    )
    plt.plot(np.real(Z_DIP), -np.imag(Z_DIP), linewidth=4, color="red", label="DP-DRT")
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif", size=20)
    plt.annotate(
        r"$10^{-2}$",
        xy=(np.real(Z_exp[20]), -np.imag(Z_exp[20])),
        xytext=(np.real(Z_exp[20]) - 2, 10 - np.imag(Z_exp[20])),
        arrowprops=dict(arrowstyle="-", connectionstyle="arc"),
    )
    plt.annotate(
        r"$10^{-1}$",
        xy=(np.real(Z_exp[30]), -np.imag(Z_exp[30])),
        xytext=(np.real(Z_exp[30]) - 2, 6 - np.imag(Z_exp[30])),
        arrowprops=dict(arrowstyle="-", connectionstyle="arc"),
    )
    plt.annotate(
        r"$1$",
        xy=(np.real(Z_exp[40]), -np.imag(Z_exp[40])),
        xytext=(np.real(Z_exp[40]), 10 - np.imag(Z_exp[40])),
        arrowprops=dict(arrowstyle="-", connectionstyle="arc"),
    )
    plt.annotate(
        r"$10$",
        xy=(np.real(Z_exp[-1]), -np.imag(Z_exp[-1])),
        xytext=(np.real(Z_exp[-1]) - 1, 10 - np.imag(Z_exp[-1])),
        arrowprops=dict(arrowstyle="-", connectionstyle="arc"),
    )
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.legend(frameon=False, fontsize=15)
    plt.xlim(0, Z_exp.max().real * 1.2)
    plt.ylim(0, Z_exp.max().real * 1.2)

    #    plt.xlim(0, 0.05)
    #    plt.ylim(0, 0.05)
    #    plt.xticks(range(0, 150, 10))
    #    plt.yticks(range(0, 110, 10))
    plt.gca().set_aspect("equal", adjustable="box")
    plt.xlabel(r"$Z_{\rm re}/\Omega$", fontsize=20)
    plt.ylabel(r"$-Z_{\rm im}/\Omega$", fontsize=20)
    fig = plt.gcf()
    size = fig.get_size_inches()
    if save:
        plt.savefig(save, dpi=100, bbox_inches="tight")
    plt.show()
    plt.close()


def plot_DRT(tau_vec, gamma_DIP, save=None):
    #    plt.semilogx(tau_vec, gamma_exact, linewidth=4, color="black", label="exact")
    plt.semilogx(tau_vec, gamma_DIP, linewidth=4, color="red", label="early stop")
    #    plt.semilogx(tau_vec, gamma_DIP_opt, linewidth=2, color="blue", label="optimal") # linestyle='None', marker='o'
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif", size=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.axis([1e-4, 5e4, -0.4, 150])
    plt.legend(frameon=False, fontsize=15)
    plt.xlabel(r"$\tau/{\rm s}$", fontsize=20)
    plt.ylabel(r"$\gamma/\Omega$", fontsize=20)
    fig = plt.gcf()
    fig.set_size_inches(5, 4)
    if save:
        plt.savefig(save, dpi=100, bbox_inches="tight")
    plt.show()
    plt.close()
