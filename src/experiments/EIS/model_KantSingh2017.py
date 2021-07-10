#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 15:40:36 2020

@author: zmg
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi

# import GP_DRT_Ciucci_Liu_2019 as GP_DRT
from scipy.optimize import minimize
import pandas as pd
import pathlib
from pathlib import Path


def get_test_file():
    EIS_data = pd.read_excel(
        list(Path.cwd().rglob("*spectrumfit*xlsx"))[0], index_col=[0]
    )
    return EIS_data


EIS_data = get_test_file()


def local_addmittance_density(omega, H, K):
    pass


def adm_p(omega):
    pass


def dS_adj(dS, H, K, r_H):
    _t = dS * (1 - 2 * H * r_H + K * r_H ** 2)


def H_adj(H, K, r_H):
    _t = (H - K * r_H) / (1 - 2 * H * r_H + K * r_H ** 2)
    return _t


def K_adj(H, K, r_H):
    _t = K / (1 - 2 * H * r_H + K * r_H ** 2)
    return _t

    """
    By considering the boundary of Helmholtz layer to be at a
    distance r H and in the direction normal to the electrode surface,
    the adjusted area dS′, adjusted mean curvature H′ and adjusted
    Gaussian curvature K′ can be obtained by Steiner’s formula 30,44
    as
    adjusted area  dS ′ = dS (1 − 2 Hr H + Kr H 2 )
    adjusted mean curvature  H ′ = (H − Kr H) / (1 − 2 Hr H + Kr H**2)
    Gaussian curvature K ′ = K / (1 − 2 Hr H + Kr H 2)
    (34)
    If the compact layer is thin with small r H , we have dS′ ≈ dS, H′
    ≈ H, and K′ ≈ K. On substituting H′ and K′ in eq 32, we
    obtain the local admittance
    """


def Table1(r, r_H):
    """
    Table 1. List of Mean (H) and Gaussian (K) Curvatures of
    Idealized Geometries along with Corrected Mean (H′) and
    Gaussian (K′) Curvature.
    Note H and K are at the electrode surface and H′ and K′ are at the
    outer Helmholtz surface.
    """
    table1 = {
        "geometry": {
            "cylindrical_pore": {
                "H": 1 / (2 * r),
                "K": 0,
                "H_adj": 1 / (2 * (r - r_H)),
                "K_adj": 0,
            },
            "cylindrical_rod": {
                "H": -1 / (2 * r),
                "K": 0,
                "H_adj": -1 / (2 * (r + r_H)),
                "K_adj": 0,
            },
            "sphere": {
                "H": -1 / r,
                "K": 1 / r ** 2,
                "H_adj": -1 / (r + r_H),
                "K_adj": 1 / (r + r_H) ** 2,
            },
            "cavity": {
                "H": 1 / r,
                "K": 1 / r ** 2,
                "H_adj": 1 / (r - r_H),
                "K_adj": 1 / (r - r_H) ** 2,
            },
        }
    }
    return pd.DataFrame(table1["geometry"]).T


def R_curv(kappa_e, H_adj, K_adj):
    term1 = delta_c(omega) / (kappa_e * (kappa_e + delta_c(omega)))
    term2 = (delta_c(omega) * (3 * kappa_e + delta_c(omega))) / (
        2 * kappa_e ** 2 * (kappa_e + delta_c(omega)) ** 2
    )
    term3 = delta_c(omega) / (2 * kappa_e ** 2 * (kappa_e + delta_c(omega)))
    expression = 1 - term1 * H_adj - term2 * H_adj ** 2 + term3 * K_adj
    return expression


def gamma_star(omega):
    pass
