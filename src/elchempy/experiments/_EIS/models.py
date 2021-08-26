# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:43:07 2020

@author: User
"""

import numpy as np
from lmfit import Model, Parameters
from lmfit.models import ExpressionModel
import inspect
import sys
from pathlib import Path
import matplotlib.pyplot as plt
from collections import OrderedDict
import pandas as pd
import math


# from .impedancepy.impedance import validation as imppy_validation
from .impedancepy.impedance.models.circuits import CustomCircuit, buildCircuit
from .impedancepy.impedance.models.circuits.elements import R, C, CPE, W, Wo, Ws

EvRHE = "E_AppV_RHE"


def get_element_from_name(name):
    excluded_chars = "0123456789_"
    return "".join(char for char in name if char not in excluded_chars)


# inspect.getmembers(sys.modules[__name__], inspect.isclass)


class ParseCircuitExpression:

    test = "R0-p(R1,p(R3-CPE1)-W0)-p(R2,CPE2-W1)-L0"

    def __init__(self, circuit_str):
        self.circuit_str = circuit_str

    def _check_series(self):
        enumt = enumerate(self.circuit_str)
        dash = [(n, i) for n, i in enumerate(self.circuit_str) if i == "-"]
        _open = [(n, i) for n, i in enumerate(self.circuit_str) if i == "("]
        _close = [(n, i) for n, i in enumerate(self.circuit_str) if i == ")"]
        dash_out = [
            n
            for n, i in dash
            if not any(c[0][0] < n < c[1][0] for c in list(zip(_open, _close)))
        ]
        dash_in = [i[0] for i in dash if i[0] not in dash_out]
        t_series = "".join(
            [t if n not in dash_out else " + " for n, t in enumerate(self.circuit_str)]
        )
        self.series_split = t_series.split(" + ")

    def _check_elements(self):

        for ser_elem in self.series_split:
            ser_elem
            # if


def _check_series(circuit_str):
    enumt = enumerate(circuit_str)
    dash = [(n, i) for n, i in enumerate(circuit_str) if i == "-"]
    _open = [(n, i) for n, i in enumerate(circuit_str) if i == "("]
    _close = [(n, i) for n, i in enumerate(circuit_str) if i == ")"]
    dash_out = [
        n
        for n, i in dash
        if not any(c[0][0] < n < c[1][0] for c in list(zip(_open, _close)))
    ]
    dash_in = [i[0] for i in dash if i[0] not in dash_out]
    t_series = "".join(
        [t if n not in dash_out else " + " for n, t in enumerate(circuit_str)]
    )
    p_split = t_series.split(" + ")
    return p_split


# best_mod_RandlesW = CustomCircuit(initial_guess=[25,100, 3E02, 0.7E-03,0.7, 1E-4 ],
#                              circuit='R0-p(R1-W1,CPE1)-L0')


class BaseModel:
    """This BaseModel class adds some standard initializations to each (LMfit type) model class as a super class template.
    It has a collection of all parameters initial guesses and contains some mathematical function definitions."""

    _params_guesses_base = (
        ("Rs", 25, True, 0.5, 200, None, None),
        ("Rl", 20, True, 0.5, 200, None, None),
        ("R_ma", 20, True, 1e-6, 1e6, None, None),
        ("Rct", 10, True, 1e-03, 1e8, None, None),
        ("R_ion", 70, True, 1e-5, 1e3, None, None),
        ("Cdlp", 1e-05, True, 1e-9, 1e-1, None, None),
        ("nDL", 0.9, True, 0, 1, None, None),
        ("nAd", 0.9, True, 0, 1, None, None),
        ("Ls", 15e-06, True, 1e-15, 1e-1, None, None),
        ("Aw", 5, True, 1e-06, 1e8, None, None),
        ("tau", 3e-02, True, 1e-08, 1e2, None, None),
        ("Rorr", 1e3, True, 1e-03, 1e6, None, None),
        ("Qad", 1e-3, True, 1e-10, 0.5, None, None),
        ("Rad", 10, True, 1e-03, 1e8, None, None),
        ("Qorr", 1e-4, True, 1e-10, 1e-1, None, None),
        ("G0", 500, True, 1e-01, 1e8, None, None),
        ("R_G", 500, True, 1e-01, 1e8, None, None),
        ("t_G", 3e-02, True, 1e-08, 1e3, None, None),
        ("Ger_D", 10, True, 1e-8, 1e5, None, None),
        ("Ger_k", 1e-2, True, 1e-8, 1e2, None, None),
        ("L_pore", 1e-07, True, 1e-11, 1e1, None, None),
        ("D_pore", 20e-09, True, 1e-11, 1e1, None, None),
        ("N_pores", 20e10, True, 1e-2, 1e15, None, None),
        ("W_orr", 5, True, 1e-03, 1e8, None, None),
        ("t_orr", 1e-03, True, 1e-06, 1e3, None, None),
        ("ratio_S", 0.2, True, 1e-4, 1e3, None, None),
        ("alpha_ma", 20 * 1e-09, True, 1e-10, 1e-7, None, None),
        ("rho", 20e7, True, 1e-2, 1e11, None, None),
        ("n_ma", 20e10, True, 1e-2, 1e15, None, None),
        ("Xi_ma", 80e-9, True, 1e-10, 1e-5, None, None),
        ("ksi", 0.2, True, 0, 1, None, None),
    )
    # W_orr,t_orr)
    # ratio_S, alpha_ma, rho, n_ma

    _params_tuple_names = ["name", "value", "vary", "min", "max", "expr", "brute_step"]
    ang = np.logspace(4, -2)

    #    model = Model(self.func,name=self.name)
    def __init__(self):
        self.classname = self.__class__.__name__
        if hasattr(self, "func"):
            self.model = Model(self.func, name=self.name)
            self.guesses_model_from_base()

            self.add_params()
            self.eval_model()
        self.set_std_model()
        self.set_repr_str()

    #        setattr(self.model, 'class') = self.__name__
    #        self.model = Model(self.func,name=self.name)
    def set_repr_str(self):
        if hasattr(self, "func"):
            _modtxt = f"{self.name}, {len(self.model.param_names)}"
        else:
            _modtxt = f"empty {self.classname}"
        self.repr_str = _modtxt

    def set_std_model(self):
        if hasattr(self, "std_model"):
            _std = self.std_model
        else:
            _std = False
        self.std_model = _std

    def __repr__(self):
        return f"{self.repr_str}"

    def __str__(self):
        return f"{self.repr_str}"

    def guesses_model_from_base(self):
        self.guesses = tuple(
            (i for i in self._params_guesses_base if i[0] in self.model.param_names)
        )

    def add_params(self):

        for p in self.guesses:
            if type(p) == tuple and len(p) == 7:
                self.model.set_param_hint(
                    **dict(zip(self._params_tuple_names[:-1], p[:-1]))
                )
        self.parameters_guesses = self.model.make_params()
        # Parameters()
        # if type(self.guesses) == tuple:
        # if all(len(i) == 7 for i in self.guesses):
        # self.parameters_guesses.add_many(*self.guesses)

    def set_model(self):
        pass

    def eval_model(self):
        self.mod_eval = self.model.eval(ang=self.ang, params=self.parameters_guesses)

    @staticmethod
    def coth_func(a):
        #  coth(a) = cosh(a) / sinh(a)
        return [
            1
            if i > 200
            else -1
            if i < -200
            else ((np.exp(2 * i) + 1) / (np.exp(2 * i) - 1))
            for i in a
        ]

    @staticmethod
    def TLM_RionCdl(ang, R_ion, Cdlp, Zpore_branch):
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Zpore_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Zpore_branch))
        Zcoth = A * BaseModel.coth_func(B)
        return Zcoth

    @staticmethod
    def TLM_Levie(R_ion, Zpore_branch, L):
        A = np.sqrt(R_ion * Zpore_branch)
        B = np.sqrt(R_ion / Zpore_branch)
        Zcoth = A * BaseModel.coth_func(B * L)
        return Zcoth

    @staticmethod
    def Warburg_Z(ang, Z0, tau):
        """Journal of The Electrochemical Society, 166 (15) F1209-F1217 (2019)
        is a simplified Warburg element.
        defines a short (finite-length) Warburg element
        Z_W = (W0 / np.sqrt(1j*ang*tau)) * np.tanh(L* np.sqrt(1j*ang*tau) )
        """  # FIXME
        Z_W = (Z0 / np.sqrt(1j * ang * tau)) * np.tanh(np.sqrt(1j * ang * tau))
        return Z_W

    @staticmethod
    def Warburg_Z_Wo(ang, Z0, tau):
        """
        is a simplified Warburg element.
        defines an open (finite-space) Warburg element
        Z_W = (W0 / np.sqrt(1j*ang*tau)) * np.coth(L* np.sqrt(1j*ang*tau) )
        """
        Z_Wo = (Z0 / np.sqrt(1j * ang * tau)) * BaseModel.coth_func(
            np.sqrt(1j * ang * tau)
        )
        return Z_Wo

    @staticmethod
    def Warburg_Z_FLW(ang, Z0, tau, L):
        """Journal of The Electrochemical Society, 158 (8) B877-B884 (2011)
        is a finite lenght Warburg element.
        Z_W = (W0 / np.sqrt(1j*ang*tau)) * np.tanh(L* np.sqrt(1j*ang*tau) )
        """
        Z_FLW = (Z0 / np.sqrt(1j * ang * tau)) * np.tanh(L * np.sqrt((1j * ang) / tau))
        return Z_FLW

    @staticmethod
    def Gerischer_Z(ang, G0, Ger_k, Ger_D, L):
        """Journal of The Electrochemical Society, 158 (8) B877-B884 (2011)
           J. Phys. Energy 2 (2020) 042001
           General expression:
               Z_G = ( Z_0 / np.sqrt(1 + 1j*ang*tau))
        Finite length
        Z_G_FLW = ( Z_0 / np.sqrt(k + 1j*ang)*D)*tanh(L * np.sqrt(k + 1j*ang)/D)
        """
        Z_G_FLW = (G0 / np.sqrt((Ger_k + 1j * ang) * Ger_D)) * np.tanh(
            L * np.sqrt((Ger_k + 1j * ang) / Ger_D)
        )
        return Z_G_FLW


#    def parse_name(self):
#        split= self.name.split('-')
#        _expr = ''
#        for i in split:
#            if


"""

class EEC_(BaseModel):

    name =
    guesses =

    def func(self,
"""


class M_2CPE(BaseModel):
    """R0-p(R1,CPE1)-p(R2,CPE2)-L0,
    Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""

    name = "R0-L0-p(R1,CPE1)-p(R2,CPE2)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr):

        expr = (
            Rs
            + Ls * (1j * ang)
            + ((1 / (1 / (Cdlp * (1j * ang) ** nDL))) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
        )
        return expr

    def expr_func(self):
        """R0-p(R1,CPE1)-p(R2,CPE2)-L0,
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        _Exprr = "Rs + Ls*(1j*ang) + ((1/(1/(Cdlp*(1j*ang)**nDL))) + 1/Rct)**-1 + ((1/(1/(Qad**1*(1j*ang)**nAd))) + 1/Rorr)**-1"
        return ExpressionModel(_Exprr, independent_vars=["ang"])


class M_2CPE_W(BaseModel):
    """R0-p(R1,CPE1)-p(R2,CPE2)-L0-W0"""

    name = "R0-L0-W0-p(R1,CPE1)-p(R2,CPE2)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, Aw):

        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Zw
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr)) ** -1
        )


class M_CPE_WCPE(BaseModel):
    """R0-p(R1,CPE1)-p(R2-W2,CPE2)-L0,
    Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""

    name = "R0-L0-p(R1,CPE1)-p(R2-W2,CPE2)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, Aw):
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Zw)) ** -1
        )


class M_RC_CPE(BaseModel):

    name = "R0-L0-p(R1-C2,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad):
        """'R0-L0-p(R1-C2,CPE1)"""
        # Zw = Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + (
                (Cdlp * (1j * ang) ** nDL)
                + 1 / (Rct + (1 / (Qad ** 1 * (1j * ang) ** 1)))
            )
            ** -1
        )


class M_RCPE_CPE(BaseModel):

    name = "R0-L0-p(R1-CPE2,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd):
        """'R0-L0-p(R1-CPE2,CPE1)"""
        # Zw = Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + (
                (Cdlp * (1j * ang) ** nDL)
                + 1 / (Rct + (1 / (Qad ** 1 * (1j * ang) ** nAd)))
            )
            ** -1
        )


# === WARBURG VARIATIONS =====
class M_RW_CPE(BaseModel):

    name = "R0-L0-p(R1-W1,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Aw):
        """'R0-L0-p(R1-W1,CPE1)"""
        # Zw = Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Aw)) ** -1
        )


class M_RWC_C(BaseModel):

    name = "R0-L0-p(R1-W1,C1)-C2"

    def func(self, ang, Rs, Ls, Cdlp, Rct, Qad, Aw):
        """'R0-L0-p(R1-W1,C1)-C2"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** 1) + 1 / (Rct + Zw)) ** -1
            + (1 / (Qad ** 1 * (1j * ang) ** 1))
        )


class M_RWCPE_C(BaseModel):

    name = "R0-L0-p(R1-W1,CPE1)-C2"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, Aw):
        """'R0-L0-p(R1-W1,CPE1)-C2"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Zw)) ** -1
            + (1 / (Qad ** 1 * (1j * ang) ** 1))
        )


class M_WCPE_C(BaseModel):

    name = "R0-L0-p(W1,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, Aw):
        """'R0-L0-p(R1-W1,CPE1)-C2"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang) ** nDL) + 1 / (Zw)) ** -1


class M_RWoCPE_C(BaseModel):

    name = "R0-L0-p(R1-Wo1,CPE1)-C2"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, Aw, tau):
        """'R0-L0-p(R1-W1,CPE1)-C2"""
        # Z = Z0/(np.sqrt(1j*omega*tau)*np.tanh(np.sqrt(1j*omega*tau)))
        Z_Wo = Aw / (np.sqrt(1j * ang * tau) * np.tanh(np.sqrt(1j * ang * tau)))
        # Zw = Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Z_Wo)) ** -1
            + (1 / (Qad ** 1 * (1j * ang)))
        )


class M_RWoCPE(BaseModel):

    name = "R0-L0-p(R1-Wo1,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Aw, tau):
        """'R0-L0-p(R1-W1,CPE1)"""
        Z_Wo = Aw / (np.sqrt(1j * ang * tau) * np.tanh(np.sqrt(1j * ang * tau)))
        return (
            Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Z_Wo)) ** -1
        )


class M_WoCPE(BaseModel):

    name = "R0-L0-p(Wo1,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Aw, tau):
        """'R0-L0-p(R1-W1,CPE1)"""
        Z_Wo = Aw / (np.sqrt(1j * ang * tau) * np.tanh(np.sqrt(1j * ang * tau)))
        return Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang) ** nDL) + 1 / (Z_Wo)) ** -1


class M_RWsCPE_C(BaseModel):

    name = "R0-L0-p(R1-Ws1,CPE1)-C2"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Qad, Aw, tau):
        """'R0-L0-p(R1-W1,CPE1)-C2"""
        # Z = Z0/(np.sqrt(1j*omega*tau)*np.tanh(np.sqrt(1j*omega*tau)))
        # Z_Wo = Awo/(np.sqrt(1j*ang*tau)*np.tanh(np.sqrt(1j*ang*tau)))
        Z_Ws = Aw * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)
        # Zw = Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Z_Ws)) ** -1
            + (1 / (Qad ** 1 * (1j * ang)))
        )


class M_RWsCPE(BaseModel):

    name = "R0-L0-p(R1-Ws1,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Aw, tau):
        """'R0-L0-p(R1-W1,CPE1)"""
        Z_Ws = Aw * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)
        return (
            Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Z_Ws)) ** -1
        )


class S_pCRWs(BaseModel):

    name = "R0-L0-p(C1,R2-Ws1)"

    def func(self, ang, Rs, Ls, Cdlp, Rct, Aw, tau):
        """R0-L0-p(C1,R2-Ws1,CPE2)"""
        Z_Ws = Aw * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)

        return Rs + Ls * (1j * ang) + (Cdlp * (1j * ang) + 1 / (Rct + Z_Ws)) ** -1


class S_pCpRWs(BaseModel):

    name = "R0-L0-p(C1,p(R1,Ws1)"

    def func(self, ang, Rs, Ls, Cdlp, Rct, Aw, tau):
        """R0-L0-p(C1,p(R1,Ws1)"""
        Z_Ws = Aw * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)

        return (
            Rs
            + Ls * (1j * ang)
            + (Cdlp * (1j * ang) + 1 / (1 / Rct + 1 / Z_Ws) ** -1) ** -1
        )


class S_pCRWsCPE(BaseModel):

    name = "R0-L0-p(C1,R2-Ws1,CPE2)"

    def func(self, ang, Rs, Ls, Cdlp, Rct, Aw, tau, Qad, nAd):
        """R0-L0-p(C1,R2-Ws1,CPE2)"""
        Z_Ws = Aw * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)

        return (
            Rs
            + Ls * (1j * ang)
            + (Cdlp * (1j * ang) + 1 / (Rct + Z_Ws) + (Qad * (1j * ang)) ** nAd) ** -1
        )


class S_pCR_coth(BaseModel):

    name = "R0-L0-p(C1,cothRct)"

    def func(self, ang, Rs, Ls, Cdlp, Rct, Qad):
        """R0-L0-p(C1,R2-Ws1,CPE2)"""

        A = np.sqrt(Rct / (1j * ang * Qad))
        B = np.sqrt(Rct * (1j * ang * Qad))

        # def coth_func(i):
        #     return 1 if i > 200 else -1 if i <-200 else ((np.exp(2*i) + 1) / (np.exp(2*i) - 1))

        Zcoth = A * self.coth_func(B)
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + (Cdlp * (1j * ang) + 1 / (Zcoth)) ** 1


class S_pCRWsRCPE(BaseModel):

    name = "R0-L0-p(C1,R2-Ws1,R3-CPE2)"

    def func(self, ang, Rs, Ls, Cdlp, Rct, Aw, tau, Rorr, Qad, nAd):
        """R0-L0-p(C1,R2-Ws1,CPE2)"""
        Z_Ws = Aw * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)

        return (
            Rs
            + Ls * (1j * ang)
            + (
                Cdlp * (1j * ang)
                + 1 / (Rct + Z_Ws)
                + 1 / (Rorr + 1 / (Qad * (1j * ang) ** nAd))
            )
            ** -1
        )


class M_WsCPE(BaseModel):

    name = "R0-L0-p(Ws1,CPE1)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Aw, tau):
        """'R0-L0-p(R1-W1,CPE1)"""
        Z_Ws = Aw * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)
        return Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang) ** nDL) + 1 / (Z_Ws)) ** -1


class N_RWCPE_CPE(BaseModel):

    name = "R0-L0-p(R1-W1,CPE1)-CPE2"

    def func(self, ang, Rs, Ls, Rct, Aw, Cdlp, nDL, Qad, nAd):
        """'R0-L0-p(R1-W1,CPE1)-CPE2"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Zw)) ** -1
            + (1 / (Qad ** 1 * (1j * ang) ** nAd))
        )


class M_RC_CPE_W(BaseModel):

    name = "R0-L0-W0-p(R1,CPE1)-p(R2,C2)"

    def func(self, ang, Rs, Ls, Aw, Rct, Cdlp, nDL, Qad, Rorr):
        """R0-p(R1,CPE1)-p(R2,C2)-L0-W0"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Zw
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** 1))) + 1 / (Rorr)) ** -1
        )


class M_CPE_WsCPE(BaseModel):

    name = "R0-L0-p(R1,CPE1)-p(R2-Ws0,CPE2)"

    def func(self, ang, Rs, Ls, Rct, Cdlp, nDL, Rorr, Aws, tau, Qad, nAd):
        """R0-L0-p(R1,CPE1)-p(R2-Ws0,CPE2)',
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        #    Z_Ws = Z_0 / (np.sqrt(1j*ang))
        Z_Ws = Aws * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)
        #    Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Z_Ws)) ** -1
        )


class M_CPE_WoCPE(BaseModel):

    name = "R0-L0-p(R1,CPE1)-p(R2-Wo0,CPE2)"

    def func(self, ang, Rs, Ls, Rct, Cdlp, nDL, Rorr, Awo, tau, Qad, nAd):
        """R0-L0-p(R1,CPE1)-p(R2-Wo0,CPE2)',
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        #    Z_Ws = Z_0 / (np.sqrt(1j*ang))
        Z_Wo = Awo / (np.sqrt(1j * ang * tau) * np.tanh(np.sqrt(1j * ang * tau)))
        #        Z_Ws = Aws*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        #    Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Z_Wo)) ** -1
        )


class M_RC_WsCPE(BaseModel):

    name = "R0-L0-p(R1,C1)-p(R2-Ws0,CPE2)"

    def func(self, ang, Rs, Ls, Rct, Cdlp, Rorr, Aws, tau, Qad, nAd):
        """R0-L0-p(R1,C1)-p(R2-Ws0,CPE2)',
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        #    Z_Ws = Z_0 / (np.sqrt(1j*ang))
        Z_Ws = Aws * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)
        #    Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang)) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Z_Ws)) ** -1
        )


class M_RC_WoCPE(BaseModel):

    name = "R0-L0-p(R1,C1)-p(R2-Wo0,CPE2)"

    def func(self, ang, Rs, Ls, Rct, Cdlp, Rorr, Awo, tau, Qad, nAd):
        """R0-L0-p(R1,C1)-p(R2-Wo0,CPE2)',
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        #    Z_Ws = Z_0 / (np.sqrt(1j*ang))
        Z_Wo = Awo / (np.sqrt(1j * ang * tau) * np.tanh(np.sqrt(1j * ang * tau)))
        #        Z_Ws = Aws*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        #    Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang)) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Z_Wo)) ** -1
        )


class N_RWoCPE_W(BaseModel):

    name = "R0-L0-p(R1-Wo0,CPE1)-W0"

    def func(self, ang, Rs, Ls, Rct, Awo, tau, Cdlp, nDL, Aw):
        "R0-L0-p(R1,Wo0-CPE1)-W0"

        Z_Wo = Awo / (np.sqrt(1j * ang * tau) * np.tanh(np.sqrt(1j * ang * tau)))
        #        Z_Wo = Awo*np.cosh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Z_Wo)) ** -1
            + Zw
        )


class N_RWoC_W(BaseModel):
    """Not used anymore, used before for fitting but was worst of all models..."""

    name = "R0-L0-p(R1-Wo0,C1)-W0"

    def func(self, ang, Rs, Ls, Rct, Awo, tau, Cdlp, Aw):
        "R0-L0-p(R1,Wo0-C1)-W0"
        Z_Wo = Awo / (np.sqrt(1j * ang * tau) * np.tanh(np.sqrt(1j * ang * tau)))
        #        Z_Wo = Awo*np.cosh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang)) + 1 / (Rct + Z_Wo)) ** -1 + Zw
        )


# ================================ #

#  ==== LEVIE VARIATIONS ==== #
# def coth_func(a):
#     #  coth(a) = cosh(a) / sinh(a)
#     return [1 if i > 200 else -1 if i <-200 else ((np.exp(2*i) + 1) / (np.exp(2*i) - 1)) for i in a]

# if  a > 200:
#     return 1
# elif a < -200:
#     return -1
# else:
#     return (np.exp(2*a) + 1) / (np.exp(2*a) - 1)
# return np.cosh(a) / np.sinh(a)


class N_RandlesTLM(BaseModel):

    name = "R0-L0-Z_RandelsTLM"

    def func(self, ang, Rs, Ls, Cdlp, Rct):
        """R0-L0-TLM(Rion,Cdl)"""

        A = np.sqrt(Rct / (1j * ang * Cdlp))
        B = np.sqrt(Rct * (1j * ang * Cdlp))
        Zcoth = A * self.coth_func(B)
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + Zcoth


#%% Model Classes with Q


class Q_Randles_RctZw(BaseModel):

    name = "RL-p(Qdl,Rct-Zw)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Aw, tau):
        """R0-L0-p(Cdl,Rct-Zw)"""
        Zw = self.Warburg_Z(ang, Aw, tau)
        Zout = Rs + Ls * (1j * ang) + ((1 / (Rct + Zw)) + 1j * ang * Cdlp ** nDL) ** -1
        return Zout


class Q_Randles_RctZwQad(BaseModel):

    name = "RL-p(Qdl,Rct-Zw-Qad)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Aw, tau, Qad):
        """R0-L0-p(Cdl,Rct-Zw)"""
        Zw = self.Warburg_Z(ang, Aw, tau)
        Zout = (
            Rs
            + Ls * (1j * ang)
            + ((1 / (Rct + Zw + (1j * ang * Qad) ** -1)) + 1j * ang * Cdlp ** nDL) ** -1
        )
        return Zout


class Q_Randles_RctZw_pQad_Rorr(BaseModel):

    name = "RL-p(Qdl,Rct-Zw-p(Qad,Rorr)"

    def func(self, ang, Rs, Ls, Cdlp, nDL, Rct, Aw, tau, Qad, Rorr):
        """R0-L0-p(Cdl,Rct-Zw)"""
        Zw = self.Warburg_Z(ang, Aw, tau)
        Zout = (
            Rs
            + Ls * (1j * ang)
            + (
                (1 / (Rct + Zw + 1 / (1 / Rorr + 1j * ang * Qad)))
                + 1j * ang * Cdlp ** nDL
            )
            ** -1
        )
        return Zout


class Q_RandlesZwORR(BaseModel):

    name = "R0-L0-p(Cdl,Rct-p(Zw,Rorr))"

    def func(self, ang, Rs, Ls, Cdlp, Rct, Aw, tau, Rorr):
        """R0-L0-p(Cdl,Rct-Zw)"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Zw_ORR = (1 / Zw + 1 / Rorr) ** -1
        Zout = Rs + Ls * (1j * ang) + ((1 / (Rct + Zw_ORR)) + 1j * ang * Cdlp) ** -1
        return Zout


class Q_TLM_Rct(BaseModel):

    name = "RL-TLM(Rct)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct):
        """R0-L0-TLM(Rion,p(Cdl,Rct-Zw))"""

        Z_RctW_branch = 1 / Rct
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_RctW(BaseModel):

    name = "RL-TLM(Rct-W)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau):
        """R0-L0-TLM(Rion,p(Cdl,Rct-Zw))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + Zw)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_ser_RctQad(BaseModel):

    name = "RL-TLM(Rct-Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Qad):
        """R0-L0-Z_TLM-(Rct-Qad)"""

        Z_RctW_branch = 1 / (Rct + (1j * ang * Qad) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_ser_RctQadW(BaseModel):

    name = "RL-TLM(Rct-Qad-W)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-Z_TLM-(Rct-Qad-W)"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + (1j * ang * Qad) ** -1 + Zw)
        Z_TLM = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Rs + Ls * (1j * ang) + Z_TLM


class Q_TLM_RctQadWo_(BaseModel):

    name = "RL-TLM(Rct-Qad-Wo)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-Z_TLM-(Rct-Qad-W)"""

        Zw = self.Warburg_Z_Wo(ang, Aw, tau)
        Zpore_branch = 1 / (Rct + (1j * ang * Qad) ** -1 + Zw)
        Z_TLM = self.TLM_RionCdl(ang, R_ion, Cdlp, Zpore_branch)
        return Rs + Ls * (1j * ang) + Z_TLM


class Q_TLM_Qdl_RctQadW(BaseModel):

    name = "RL-TLM(Qdl-(Rct-W-Qad))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, nDL, Rct, Aw, tau, Qad):
        """R0-L0-Z_TLM-(Rct-Qad-W)"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        # Z_RctW_branch = (1 / (Rct + Zw))
        Z_RctW_branch = 1 / (Rct + Zw + (1j * ang * Qad) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp ** nDL) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp ** nDL) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_QadW(BaseModel):

    name = "RL-TLM(Qad-W)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Aw, tau, Qad):
        """R0-L0-Z_TLM-(Qad-W)"""
        Zw = self.Warburg_Z_Wo(ang, Aw, tau)
        Zpore_branch = 1 / ((1j * ang * Qad) ** -1 + Zw)
        Z_TLM = self.TLM_RionCdl(ang, R_ion, Cdlp, Zpore_branch)
        return Rs + Ls * (1j * ang) + Z_TLM


class Q_TLM_pQad_W(BaseModel):

    name = "RL-TLM(p(Qad,W))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Aw, tau, Qad):
        """R0-L0-Z_TLM-(Qad-W)"""
        Zw = self.Warburg_Z_Wo(ang, Aw, tau)
        Zpore_branch = 1 / (1 / ((1j * ang * Qad) + 1 / Zw))
        Z_TLM = self.TLM_RionCdl(ang, R_ion, Cdlp, Zpore_branch)
        return Rs + Ls * (1j * ang) + Z_TLM


class Q_TLM_pRctQad_W(BaseModel):

    name = "RL-TLM(p(Qad-Rct,W))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-Z_TLM-(Rct-Qad,W)"""
        Zw = self.Warburg_Z_Wo(ang, Aw, tau)
        Zpore_branch = 1 / (1 / (1 / ((1j * ang * Qad) ** -1 + Rct) + 1 / Zw))
        Z_TLM = self.TLM_RionCdl(ang, R_ion, Cdlp, Zpore_branch)
        return Rs + Ls * (1j * ang) + Z_TLM


class Q_TLM_Qdl_RctW_pQadRorr(BaseModel):

    name = "RL-TLM(Qdl-(Rct-W-p(Qad,Rorr))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, nDL, Rct, Aw, tau, Qad, Rorr):
        """R0-L0-Z_TLM-(Rct-Qad-W)"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + Zw + 1 / ((1j * ang * Qad + 1 / Rorr)))
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp ** nDL) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp ** nDL) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_W_pRct_Qad(BaseModel):
    name = "RL-TLM(W-p(Qad,Rct))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        Zw = self.Warburg_Z(ang, Aw, tau)
        Zpore_branch = 1 / (Zw + 1 / (1j * ang * Qad + 1 / Rct))
        Z_TLM = self.TLM_RionCdl(ang, R_ion, Cdlp, Zpore_branch)
        return Rs + Ls * (1j * ang) + Z_TLM


class Q_TLM_W_pRct_QadRorr(BaseModel):
    name = "RL-TLM(W-p(Qad-Rorr,Rct))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad, Rorr):
        Zw = self.Warburg_Z(ang, Aw, tau)
        Zpore_branch = 1 / (Zw + 1 / (1 / ((1j * ang * Qad) ** -1 + Rorr) + 1 / Rct))
        Z_TLM = self.TLM_RionCdl(ang, R_ion, Cdlp, Zpore_branch)
        return Rs + Ls * (1j * ang) + Z_TLM


class Q_TLM_Gerischer_ser_RctQadW(BaseModel):

    name = "RL-TLM(Rct-Qad-Ger)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Qad, G0, Ger_k, Ger_D, L_pore):
        """R0-L0-Z_TLM-(Rct-Qad-W)"""

        Z_Ger = self.Gerischer_Z(ang, G0, Ger_k, Ger_D, L_pore)
        # Z_RctW_branch = (1 / (Rct + Zw))
        Z_RctW_branch = 1 / (Rct + (1j * ang * Qad) ** -1 + Z_Ger)

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_RctCdl_TLM(BaseModel):

    name = "RL-p(Rion,Cdl)-TLM-p(Rct,Qad)"

    def func(self, ang, Rs, Ls, Cdlp, Rct, R_ion, Qad):
        """R0-L0-(Cdl,Rct)-Z_TLM))"""

        Z_Rct_Cdl = (1 / Rct + (1j * ang * Cdlp)) ** -1
        A = np.sqrt(R_ion / (1j * ang * Qad))
        B = np.sqrt(R_ion * (1j * ang * Qad))
        # def coth_func(i):
        #     return 1 if i > 200 else -1 if i <-200 else ((np.exp(2*i) + 1) / (np.exp(2*i) - 1))
        Zcoth = A * self.coth_func(B)
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + Z_Rct_Cdl + Zcoth


class Q_TLM_W_Rside(BaseModel):

    name = "RL-p(TLM-W-Rct,Rorr-Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-Zw))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + (1 / (Rct + Zw)))
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion / (Rct + Zw)))
        Zcoth = A * self.coth_func(B)
        Rside = Rorr + (1j * ang * Qad)
        return Rs + Ls * (1j * ang) + (1 / Zcoth + 1 / Rside) ** -1


class Q_TLM_WarburgRorr(BaseModel):

    name = "RL-TLM-Rct-p(W,Rorr)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + (1 / Zw + 1 / Rorr) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_Rct_pQadW(BaseModel):

    name = "RL-TLM-Rct-p(W,Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + (1 / Zw + (1j * ang * Qad)) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_W_pQadRct(BaseModel):

    name = "RL-TLM-W-p(Rct,Qad)"
    std_model = False

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Zw + (1 / Rct + (1j * ang * Qad)) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_RctQad_p_WRorr(BaseModel):

    name = "RL-TLM-Rct-Qad-p(Rorr,W)"
    std_model = True

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + (1j * ang * Qad) ** -1 + (1 / Rorr + 1 / Zw) ** -1)
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Rs + Ls * (1j * ang) + Zcoth
        # A = np.sqrt(R_ion)/ np.sqrt((1j*ang*Cdlp) + Z_RctW_branch)
        # B = np.sqrt(R_ion*(1j*ang*Cdlp) + (R_ion * Z_RctW_branch))
        # Zcoth = A* self.coth_func(B)


class Q_TLM_RctQad_p_Rorr_WQorr(BaseModel):

    name = "RL-TLM-Rct-Qad-p(Rorr,W-Qorr)"
    std_model = True

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad, Qorr):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (
            Rct
            + (1j * ang * Qad) ** -1
            + (1 / Rorr + 1 / (Zw + (1j * ang * Qorr) ** -1)) ** -1
        )
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Rs + Ls * (1j * ang) + Zcoth


# class Q_TLM_RctQad_p_WRorr(BaseModel):

#     name = 'RL-TLM-Rct-Qad-p(Rorr,W-Qorr)'
#     std_model = True

#     def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad, Qorr):
#         '''R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))'''

#         Zw = self.Warburg_Z(ang, Aw, tau)
#         Z_RctW_branch = (1 / (Rct + (1j*ang*Qad)**-1 + (1/Rorr + 1/ (Zw+(1j*ang*Qorr)**-1) )**-1))
#         Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
#         return Rs + Ls*(1j*ang) +  Zcoth


class Q_TLM_p_RctQadW_Rorr(BaseModel):

    name = "RL-TLM(p(Rct-Qad-W,Rorr))"
    std_model = False

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / ((1 / (Rct + Zw + (1j * ang * Qad) ** -1) + 1 / Rorr) ** -1)
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Rs + Ls * (1j * ang) + Zcoth

    # A = np.sqrt(R_ion)/ np.sqrt((1j*ang*Cdlp) + Z_RctW_branch)
    # B = np.sqrt(R_ion*(1j*ang*Cdlp) + (R_ion * Z_RctW_branch))
    # Zcoth = A* self.coth_func(B)
    # return Rs + Ls*(1j*ang) +  Zcoth


class Q_TLM_W_p_QadRorr(BaseModel):

    name = "RL-TLM-Rct-W-p(Rorr,Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        # Z_RctW_branch = (1 / (Rct + Zw))
        Z_RctW_branch = 1 / (Rct + Zw + (1 / Rorr + (1j * ang * Qad)) ** -1)

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_noW_p_Qad_QorrRorr(BaseModel):

    name = "RL-TLM-Rct-p(Rorr-Qorr,Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Rorr, Qad, Qorr):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        # Z_RctW_branch = (1 / (Rct + Zw))
        Z_RctW_branch = 1 / (
            Rct + (1 / (Rorr + (1j * ang * Qorr) ** -1) + (1j * ang * Qad)) ** -1
        )

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_W_p_RctQad_Rorr(BaseModel):

    name = "RL-TLM-W-p(Rorr,Rct-Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)

        # Z_RctW_branch = (1 / (Rct + Zw))
        Z_RctW_branch = 1 / (Zw + (1 / Rorr + 1 / (Rct + (1j * ang * Qad) ** -1)) ** -1)

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        # def coth_func(i):
        #     return 1 if i > 200 else -1 if i <-200 else ((np.exp(2*i) + 1) / (np.exp(2*i) - 1))
        Zcoth = A * self.coth_func(B)
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_p_QadRorr(BaseModel):

    name = "RL-TLM-p(Rorr,Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Qad,Rorr)))"""

        Z_RctW_branch = 1 / (Rct + (1 / Rorr + (1j * ang * Qad)) ** -1)

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


# class Q_TLM_Warburg_RctQad_pWRorr(BaseModel):

#     name = 'RL-TLM(Rct-Qad-p(W,Rorr))'
#     std_model=False
#     def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
#         '''R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))'''

#         Zw = self.Warburg_Z(ang, Aw, tau)

#         Z_RctW_branch = (1 / (Rct + (1j*ang*Qad)**-1 + (1/Zw + 1/Rorr)**-1 ) )

#         A = np.sqrt(R_ion)/ np.sqrt((1j*ang*Cdlp) + Z_RctW_branch)
#         B = np.sqrt(R_ion*(1j*ang*Cdlp) + (R_ion * Z_RctW_branch))
#         Zcoth = A* self.coth_func(B)
#         return Rs + Ls*(1j*ang) +  Zcoth


class Q_TLM_Warburg_Rct_pW_QadRorr(BaseModel):

    name = "RL-TLM(Rct-p(W,Qad-Rorr))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-Z_TLM-Rct-p(W,Qad-Rorr)"""

        Zw = self.Warburg_Z(ang, Aw, tau)

        Z_RctW_branch = 1 / (Rct + (1 / Zw + 1 / (Rorr + (1j * ang * Qad) ** -1)) ** -1)

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_Warburg_Rct_p_Qad_WRorr(BaseModel):

    name = "RL-TLM(Rct-p(Qad,W-Rorr))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-Z_TLM-Rct-p(W,Qad-Rorr)"""

        Zw = self.Warburg_Z(ang, Aw, tau)

        Z_RctW_branch = 1 / (Rct + ((1j * ang * Qad) + 1 / (Rorr + Zw)) ** -1)

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class F_noRs_TLM_Warburg_pQadWRorr(BaseModel):

    name = "L-TLM(Rct-p(Qad-W,Rorr))"

    def func(self, ang, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + (1 / (Zw + (1j * ang * Qad) ** -1) + 1 / Rorr) ** -1)
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Ls * (1j * ang) + Zcoth


class F_noRs_TLM_QadWRct(BaseModel):

    name = "L-TLM(Rct-Qad-W)"

    def func(self, ang, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""
        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + Zw + (1j * ang * Qad) ** -1)
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Ls * (1j * ang) + Zcoth


class Q_TLM_Warburg_pQadWRorr(BaseModel):

    name = "RL-TLM(Rct-p(Qad-W,Rorr))"
    std_model = True

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)

        Z_RctW_branch = 1 / (Rct + (1 / (Zw + (1j * ang * Qad) ** -1) + 1 / Rorr) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_pRctQadW_Rorr(BaseModel):

    name = "RL-TLM(p(Rct-Qad-W,Rorr))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / ((1 / (Rct + Zw + (1j * ang * Qad) ** -1) + 1 / Rorr) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_Warburg_pQadW_QorrRorr(BaseModel):

    name = "RL-TLM(Rct-p(Qad-W,Qorr-Rorr))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad, Qorr):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)

        Z_RctW_branch = 1 / (
            Rct
            + (1 / (Zw + (1j * ang * Qad) ** -1) + 1 / (Rorr + (1j * ang * Qorr) ** -1))
            ** -1
        )

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_2W_RctWorr_pQadW_Rorr(BaseModel):

    name = "RL-TLM(Rct-Worr-p(Qad-W,Rorr))"
    std_model = False

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad, W_orr, t_orr):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Zw_orr = self.Warburg_Z(ang, W_orr, t_orr)

        Z_RctW_branch = 1 / (
            Rct + Zw_orr + (1 / (Zw + (1j * ang * Qad) ** -1) + 1 / (Rorr)) ** -1
        )

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_2W_p_QadW_RorrWorr(BaseModel):

    name = "RL-TLM(p(Qad-W,Rorr-Worr))"
    std_model = True

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Aw, tau, Qad, W_orr, t_orr, Rorr):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""
        Zw = self.Warburg_Z(ang, Aw, tau)
        Zw_orr = self.Warburg_Z(ang, W_orr, t_orr)
        Z_RctW_branch = 1 / (
            (1 / (Zw + (1j * ang * Qad) ** -1) + 1 / (Zw_orr + Rorr)) ** -1
        )
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Rs + Ls * (1j * ang) + Zcoth
        # Rct + Zw_orr +
        # A = np.sqrt(R_ion)/ np.sqrt((1j*ang*Cdlp) + Z_RctW_branch)
        # B = np.sqrt(R_ion*(1j*ang*Cdlp) + (R_ion * Z_RctW_branch))
        # Zcoth = A* self.coth_func(B)

        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_2W_Rct_p_QadW_RorrWorr(BaseModel):
    name = "RL-TLM(Rct-p(Qad-W,Rorr-Worr))"
    std_model = True

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad, Rorr, W_orr, t_orr):
        Zw = self.Warburg_Z(ang, Aw, tau)
        Zw_orr = self.Warburg_Z(ang, W_orr, t_orr)
        Z_RctW_branch = 1 / (
            Rct + (1 / (Zw + (1j * ang * Qad) ** -1) + 1 / (Zw_orr + Rorr)) ** -1
        )
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_RctW_pQadR_Rorr(BaseModel):

    name = "RL-TLM(Rct-W-p(Qad-Rad,Rorr))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Rorr, Qad, Rad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (
            Rct + Zw + (1 / (Rad + (1j * ang * Qad) ** -1) + 1 / (Rorr)) ** -1
        )

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_Warburg_p_Qad(BaseModel):

    name = "RL-TLM(Rct-p(W,Qad))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-Z_TLM-p(W,Qad)"""

        Zw = self.Warburg_Z(ang, Aw, tau)

        Z_RctW_branch = 1 / (Rct + (1 / Zw + (1j * ang * Qad)) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class Q_TLM_W_p_Qad(BaseModel):
    # FIXME
    name = "RL-TLM(Rct-p(W-Rorr,Qad))"
    std_model = True

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad, Rorr):
        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + (1 / ((1j * ang * Qad) + 1 / (Zw + Rorr))))
        Zcoth = self.TLM_RionCdl(ang, R_ion, Cdlp, Z_RctW_branch)
        return Rs + Ls * (1j * ang) + Zcoth


#%% Model Classes with N
class N_TLM_Warburg_Rctp_Qad(BaseModel):

    name = "R0-L0-Z_TLM-p(Rct-W,Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-Z_TLM-p(Rct-W,Qad)"""

        Zw = self.Warburg_Z(ang, Aw, tau)

        Z_RctW_branch = 1 / (Rct + (1 / Zw + (1j * ang * Qad)) ** -1)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class N_TLM_Warburg_ser_Qad(BaseModel):

    name = "R0-L0-Z_TLM-(W-Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-Z_TLM-p(Warburg-Qad)"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        Z_RctW_branch = 1 / (Rct + Zw + (1j * ang * Qad) ** -1)

        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class N_TLM_Qad(BaseModel):
    """Model requires the extra Warburg Zw for best fitting"""

    name = "R0-L0-Z_TLM-Rct-Qad"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Qad):
        """R0-L0-Z_TLM-Rct-Qad"""

        Zw = (1j * ang * Qad) ** -1
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + (1 / (Rct + Zw)))
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion / (Rct + Zw)))
        Zcoth = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Zcoth


class N_noL_TLM_Warburg(BaseModel):
    """Model requires the extra HF inductor Ls for best fitting"""

    name = "R0-Z_TLM-Warburg"

    def func(self, ang, Rs, R_ion, Cdlp, Rct, Aw, tau):
        """R0-TLM(Rion,p(Cdl,Rct-Zw))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + (1 / (Rct + Zw)))
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion / (Rct + Zw)))
        Z_TLM_W = A * self.coth_func(B)
        return Rs + Z_TLM_W


class N_TLM_Warburg_Qad(BaseModel):

    name = "R0-L0-Z_TLM-(Rct-W-Qad)"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Aw, tau, Qad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-Zw-Qad))"""

        Zw = self.Warburg_Z(ang, Aw, tau)
        A = np.sqrt(R_ion) / np.sqrt(
            (1j * ang * Cdlp) + (1 / (Rct + Zw + (1j * ang * Qad) ** -1))
        )
        B = np.sqrt(
            R_ion * (1j * ang * Cdlp) + (R_ion / (Rct + Zw + (1j * ang * Qad) ** -1))
        )
        Z_TLM_W = A * self.coth_func(B)
        return Rs + Ls * (1j * ang) + Z_TLM_W


#%% Model Classes with Levie model


class N_Levie_pore(BaseModel):

    name = "Levie_pore"

    def func(self, ang, rho_el, l_pore, n_pore, r_pore, npores, Z_eq):
        # Rct, Awo, tau, Cdlp, Aw):

        R_0 = rho_el / (np.pi * r_pore ** 2)
        Z_0 = Z_eq / (2 * np.pi * r_pore)
        Z_deLevie = np.sqrt(R_0 * Z_0) * BaseModel.coth_func(
            l_pore * np.sqrt(R_0 * Z_0)
        )


class N_Levie_W(BaseModel):
    """Barcia et al., Electrochimica Acta 47 (2002) 2109 -2116"""

    name = "Levie"

    def func(self, ang, rho_el, l_pore, n_pore, r_pore, npores, Z_eq):
        # Rct, Awo, tau, Cdlp, Aw):

        R_0 = rho_el / (np.pi * r_pore ** 2)
        Z_0 = Z_eq / (2 * np.pi * r_pore)

        Z_deLevie = np.sqrt(R_0 * Z_0) * BaseModel.coth_func(
            l_pore * np.sqrt(R_0 * Z_0)
        )


class N_Levie_ESR_Cdl(BaseModel):
    """
    L.M. Da Silva, R. Cesar, C.M.R. Moreira, J.H.M. Santos, L.G. De Souza,
    B.M. Pires, R. Vicentini, W. Nunes, H. Zanin, Reviewing the fundamentals of supercapacitors and the
    difficulties involving the analysis of the electrochemical findings obtained for porous electrode materials,
    Energy Storage Materials, https://doi.org/10.1016/j.ensm.2019.12.015.
    [127,254,257]
    """

    name = "ESR_Levie"

    _guesses = (
        ("Rs", 20, True, 0.1, 200, None, None),
        ("Rct", 280, True, 1e-03, 1e8, None, None),
        ("c_dl", 1, True, 1e-9, 1e4, None, None),
        ("Qad", 3e-4, True, 1e-10, 200, None, None),
        ("rho_el", 5, True, 1e-9, 1e4, None, None),
        ("N_pores", 2e7, True, 1e5, 1e12, None, None),
        ("l_pore", 1e-6, True, 1e-10, 1e-5, None, None),
        ("r_pore", 20e-9, True, 1e-10, 1e-5, None, None),
    )
    # ('kappa_CL', 5E-2, True, 1E-9, 10, None, None))

    _extra = (
        ("r_pore_mean", 40, True, 1, 200, None, None),
        ("pore_psd_var", 7e-01, True, 1e-03, 1, None, None),
        ("Qad", 3e-4, True, 1e-10, 200, None, None),
        ("Rside", 550, True, 1e-8, 1e4, None, None),
    )

    guesses = _guesses
    # guesses = OrderedDict({'Rs':20, 'Rct':5, 'Cdlp':7E-06,'rho_el':10,
    #                        'l_pore':1E-06,'r_pore':1E-06,'n_pores':1E10})
    def func(self, ang, Rs, Rct, c_dl, rho_el, l_pore, r_pore, N_pores, Qad):
        #
        # Rct, Awo, tau, Cdlp, Aw):
        # R_0 = rho_el / (np.pi*r_pore**2)
        # Z_0 = Z_eq / (2*np.pi*r_pore)
        # Cdl_pore = Cdlp / (2.*np.pi*n_pores)
        # ((np.exp(2*i) + 1) / (np.exp(2*i) - 1))

        Cdl_star = c_dl / (2 * np.pi * N_pores * r_pore * l_pore)
        B = np.sqrt(2 * rho_el * l_pore ** 2 * Cdl_star / r_pore)
        # np.sqrt( (2. * rho_el * l_pore**2. * Cdl_pore) / r_pore )
        #  coth(a) = cosh(a) / sinh(a)
        # coth_fun = ( np.cosh(B * np.sqrt(1j*ang))/
        Z_pore = (B * np.sqrt(1j * ang)) * self.coth_func(B * np.sqrt(1j * ang))

        Zads = Rct + 1.0 / (1j * ang * Qad)
        Z_par = 1.0 / (1.0 / Z_pore + 1.0 / Zads)
        Z_out = Rs + Z_par
        # Z_deLevie = Rs + ( Rct / (B*np.sqrt(1j*ang))) * self.coth_func(B * np.sqrt(1j*ang))
        return Z_out


class P_Levie_Gassa(BaseModel):
    """
    Gassa, L. M., J. R. Vilche, M. Ebert, K. Jüttner, and W. J. Lorenz.
    “Electrochemical Impedance Spectroscopy on Porous Electrodes.”
    Journal of Applied Electrochemistry 20, no. 4 (July 1, 1990): 677–85. https://doi.org/10.1007/BF01008882.
    """

    name = "Levie_Gassa"
    guesses = OrderedDict(
        {"Ls": 5e-05, "Rs": 20, "Rel_p": 50, "Rct": 0.1, "Cdlp": 1, "nDL": 0.7}
    )

    def func(self, ang, Ls, Rs, Rel_p, Rct, Cdlp, nDL, Qad):
        """System parameters = Rct, Cdlp, r_pore,l_pore, and rho_el
        _alpha = r / 2*rho*l_pore**2
        A = _alpha * Rct
        B = Cdlp/_alpha
        Rct == R_el_pore = rho_el*l_pore / (n_pores*np.pi*r_pore**2)
        nDL :_alpha, Rct :  A
        """
        _alpha = nDL
        A = Rct
        B = Cdlp
        _lambda = (1.0 + (1j * ang * A * B) ** _alpha) / A
        return (
            Ls
            + Rs
            + (Rel_p * np.sqrt(1 / _lambda) * np.sqrt(BaseModel.coth_func(_lambda)))
            + (1 / (Qad ** 1 * (1j * ang) ** 1))
        )


def pore_PSD_simul(r_pore, r_pore_mean, pore_psd_var, N_pores):
    _preexp = 1.0 / (pore_psd_var * np.sqrt(2.0 * np.pi))
    _exp = np.exp(-1.0 * (r_pore - r_pore_mean) ** 2 / 2.0 * pore_psd_var ** 2)
    return N_pores * _preexp * _exp


#  test plt.plot(np.log10(np.linspace(10,1E3)),[guess_PSD(i,20,0.3E-01,1E9) for i in np.linspace(1,1000)])


def ink_recipe():
    """dropping 6 uL on 0.238 cm2"""
    _std_ink_20 = {"EtOH": 2.84, "H2O": 1.66, "Nafion_5wt%": 0.5}
    std_ink = {k: val / 20 for k, val in _std_ink_20.items()}
    sum_ink = np.sum(list(std_ink.values()))
    std_ink["sum"] = sum_ink
    std_ink.update({f"{k}_vol%": 100 * val / sum_ink for k, val in std_ink.items()})
    std_ink["Nafion/Catalyst"] = std_ink["Nafion_5wt%"] * 1 * 0.05 * 1000 / 5

    return std_ink


class T_Levie_Ober(BaseModel):
    """
    Obermaier, Michael, Aliaksandr S. Bandarenka, and Cyrill Lohri-Tymozhynsky.
    “A Comprehensive Physical Impedance Model of Polymer Electrolyte Fuel Cell Cathodes
    in Oxygen-Free Atmosphere.” Scientific Reports 8, no. 1 (March 21, 2018): 4933.
    https://doi.org/10.1038/s41598-018-23071-5.
    Include a Levie style PSD-TLM and possible side reactions.
    """

    name = "Levie_Ober"
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    # guesses = Parameters()

    _guesses = (
        ("Rs", 20, True, 0.1, 200, None, None),
        ("Rct", 280, True, 1e-03, 1e8, None, None),
        ("Qad", 3e-4, True, 1e-10, 200, None, None),
        ("Rside", 550, True, 1e-8, 1e4, None, None),
        ("c_dl", 795, True, 1e-9, 1e4, None, None),
        ("kappa_CL", 5e-2, True, 1e-9, 10, None, None),
    )

    _extra = (
        ("r_pore_mean", 40, True, 1, 200, None, None),
        ("pore_psd_var", 7e-01, True, 1e-03, 1, None, None),
        ("N_pores", 2e7, True, 1e5, 1e12, None, None),
        ("lenpore", 1e-6, True, 1e-10, 1e-5, None, None),
    )
    guesses = (*_guesses, *_extra)

    r_pore_mean, pore_psd_var, N_pores = 20, 0.05, 1e9

    # OrderedDict({'Rs' :20, 'Rct':25, 'Qad':1E-01, 'Rside':4E2, 'c_dl':0.5,
    #                        'kappa_CL':4, 'r_pore':30E-09, 'l_pore':1E-09})
    def func(
        self,
        ang,
        Rs,
        Rct,
        Qad,
        Rside,
        c_dl,
        kappa_CL,
        r_pore_mean,
        pore_psd_var,
        N_pores,
        lenpore,
    ):
        """System parameters =
        kappa_CL : ionomer cond
        c_dl : specific DL cap
        l_pore : pore length, r_pore
        VF = ratio of ionomer volume to total pore volume
        function Zout = PASR(x, freq, rpore, npore, lenpore, RelHawa, VF)
        % x = fitting parameter
        % rpore = pore radii in m
        % npore = number of pores of given radius
        % lenpore = length of pores in m = t_ccl ⋅ tortuosity == 7E-06
        % MemArea = active area in m2
        % RelHawa = electronic resistance of the cell
        % VF = ratio of ionomer volume to total pore volume, == 1

        Rmembrane=x(1) % Membrane resistance as 1. fitting parameter (in Ohm)
        cdl=x(2) % Specific double layer capcitance as 2. fitting parameter (in F/m^2)
        kappa=x(3) % CL's ionomer conductivity as 3. fitting parameter (in 1/(Ohm*m))
        Rads=x(4) % Adsorption resistance as 4. fitting parameter (in Ohm)
        Cads=x(5) % Adsorption capacitance 5. fitting parameter (in F)
        RSR=x(6);
        Rct, Cdlp, r_pore,l_pore, and rho_el
        _alpha = r / 2*rho*l_pore**2
        A = _alpha * Rct
        B = Cdlp/_alpha
        Rct == R_el_pore = rho_el*l_pore / (n_pores*np.pi*r_pore**2)
        nDL :_alpha, Rct :  A
        R0=zeros(length(rpore),1);                                  % Electrolyte resistance in a pore (in Ohm/m)
        Z0=ones(length(w),length(rpore))*(1+1j);                    % Inerfacial impedance in (Ohm*m)
        Zp=ones(length(w),length(rpore))*(1+1j);                    % Impedance of a single pore (in Ohm)
        Ztot=ones(length(w),1)*(1+1j);                              % Impedance of parallel connection of all pores (in Ohm)

        Zads=ones(length(w),1)*(1+1j);

        """
        # lenpore = 20E-09
        VF = 1
        pore_radii = np.linspace(1 * 1e-9, 1e3 * 1e-9)

        R0 = np.zeros((len(pore_radii), 1))
        Z0 = (np.ones((len(ang), len(pore_radii)))) * (1 + 1j)
        Zpore = (np.ones((len(ang), len(pore_radii)))) * (1 + 1j)
        Ztot = np.ones((len(ang), 1)) * (1 + 1j)

        Zads = np.ones((len(ang), 1)) * (1 + 1j)

        PSD_simul = [
            pore_PSD_simul(r_pore * 1e9, r_pore_mean, pore_psd_var, N_pores)
            for r_pore in pore_radii
        ]
        for _na, _a in enumerate(ang):
            for nrp, rpore in enumerate(pore_radii):

                _t = rpore * (1 - np.sqrt(1 - 1 / VF))
                R0[nrp] = 1 / (kappa_CL * np.pi * ((rpore ** 2) - (rpore - _t) ** 2))
                Z0[_na, nrp] = 1 / (2 * np.pi * rpore * 1j * c_dl * _a)
                Zpore[_na, nrp] = np.sqrt(R0[nrp] * Z0[_na, nrp]) * BaseModel.coth_func(
                    lenpore * np.sqrt(R0[nrp] / Z0[_na, nrp])
                )

            Ztottemp = 0
            for nrp, rpore in enumerate(pore_radii):
                Ztottemp = Ztottemp + PSD_simul[nrp] / Zpore[_na, nrp]
                Ztot[_na] = 1 / Ztottemp

            Zads[_na] = Rct + 1.0 / (1j * _a * Qad)

        Z_par = 1.0 / (1.0 / Ztot + 1.0 / Zads + 1.0 / Rside)

        Z_sideads = 1 / (1.0 / Zads + 1.0 / Rside)
        Z_out = Rs + Z_par
        Z_out = Z_out.flatten()
        return Z_out

        # Z_dL_0 = 1. / np.sqrt(1j * kappa_CL * r_pore**3 * ang * c_dl)
        # Z_dL_coth = coth_func(np.srqt( (2*1j*c_dl * l_pore**2) / (kappa_CL * r_pore) ))
        # Z_deLevie = Z_dL_0 + Z_dL_coth


def testplots(Z_out, Ztot, Z_sideads):
    plt.plot(Z_out.real, -1 * Z_out.imag)
    plt.plot(Ztot.real, -1 * Ztot.imag)
    plt.plot(Z_sideads.real, -1 * Z_sideads.imag)


#%% Model Classes PSD Song


class L_Song_PSDTLM(BaseModel):
    """
    Z Cathode = L(iω) + R e + [R C + Z W ] γ 1 coth ( γ 1 ) / (1 + Y ( i ω )^P [ R C + Z W ])

    Include a Levie style PSD-TLM and possible side reactions.
    """

    name = "Song_PSDTLM"
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    # guesses = Parameters()

    _guesses = (
        ("Rs", 20, True, 0.1, 200, None, None),
        ("Rct", 50, True, 1e-03, 1e8, None, None),
        ("Qad", 1e-3, True, 1e-10, 200, None, None),
        ("Rside", 150, True, 1e-8, 1e4, None, None),
        ("c_dl", 1e-4, True, 1e-9, 10, None, None),
        ("kappa_CL", 1e-2, True, 1e-9, 10, None, None),
    )

    _extra = (
        ("r_pore_mean", 20, True, 1, 200, None, None),
        ("pore_psd_var", 5e-02, True, 1e-03, 1, None, None),
        ("N_pores", 1e9, False, 1e6, 1e12, None, None),
    )
    guesses = (*_guesses, *_extra)

    r_pore_mean, pore_psd_var, N_pores = 20, 0.05, 1e9

    # OrderedDict({'Rs' :20, 'Rct':25, 'Qad':1E-01, 'Rside':4E2, 'c_dl':0.5,
    #                        'kappa_CL':4, 'r_pore':30E-09, 'l_pore':1E-09})
    def func(
        self,
        ang,
        Rs,
        Rct,
        Qad,
        Rside,
        c_dl,
        kappa_CL,
        r_pore_mean,
        pore_psd_var,
        N_pores,
        alpha,
        y,
    ):

        # R_mu, n_mu, lamda_mu, alpha_mu

        def Z_pore(alpha):
            # Zp = Zp_star / R_pore
            _ar = 1 / alpha
            _real = (np.sinh(_ar) - np.sin(_ar)) / (np.cosh(_ar) - np.cos(_ar))
            _imag = (np.sinh(_ar) + np.sin(_ar)) / (np.cosh(_ar) - np.cos(_ar))
            Zp_star = alpha * (_real - 1j * _imag)
            return Zp_star

        def PSD_simul(y):
            np.exp(-0.5 * y ** 2)

        Ztot_1 = np.sum(1 / Z_pore(alpha) * PSD_simul(y))

        # Z_out = L(iω) + R e + [ R C + Z W ] γ 1 coth ( γ 1 ) / (1 + Y ( i ω )^P [ R C + Z W ])

        return Ztot_1

    # Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang)) + 1 / (Rct + Z_G)) ** -1


class F_fractal_pore(BaseModel):
    name = "RL-TLM_ma(TLM_me(TLM_mi))"

    def func(self, ang, Ls, Cdlp, ratio_S, alpha_ma, rho, n_ma, Xi_ma, ksi):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        # Zw = self.Warburg_Z(ang, Aw, tau)

        # def Pore():
        # Rsol_ma = R_ma
        # Rsol_me = R_ion
        # Rsol_mi = Rct
        #
        Xi_me = ratio_S * Xi_ma
        alpha_me = ratio_S * alpha_ma

        Xi_mi = ratio_S ** 2 * Xi_ma
        alpha_mi = ratio_S ** 2 * alpha_ma
        n_mi = 1 / ((np.sqrt(3) / 4) * alpha_mi ** 2)

        Zss_pore_mi = 1 / (1j * ang * Cdlp)
        Zs_pore_mi = Zss_pore_mi / (3 * alpha_mi)
        Rsol_mi = rho / ((np.sqrt(3) / 4) * alpha_mi ** 2)

        Zmi_A = np.sqrt(Rsol_mi * Zs_pore_mi)
        Zmi_B = np.sqrt((Rsol_mi / Zs_pore_mi) * Xi_mi)
        Z_mi = Zmi_A * self.coth_func(Zmi_B)

        Zss_me = 1 / (1 / (Z_mi / (n_mi * (1 - ksi))) + 1j * ang * ksi * Cdlp)
        Zs_pore_me = Zss_me / (3 * alpha_me)
        Rsol_me = rho / ((np.sqrt(3) / 4) * alpha_me ** 2)
        n_me = 1 / ((np.sqrt(3) / 4) * alpha_mi ** 2)

        Zme_A = np.sqrt(Rsol_me * Zs_pore_me)
        Zme_B = np.sqrt((Rsol_me / Zs_pore_me) * Xi_me)
        Z_me = Zme_A * self.coth_func(Zme_B)

        Zss_ma = 1 / (1 / (Z_me / (n_me * (1 - ksi))) + 1j * ang * Cdlp)

        Zs_pore_ma = Zss_me / (3 * alpha_me)
        R_sol_ma = rho / ((np.sqrt(3) / 4) * alpha_ma ** 2)

        Zma_A = np.sqrt(R_sol_ma * Zs_pore_ma)
        Zma_B = np.sqrt((R_sol_ma / Zs_pore_ma) * Xi_ma)
        Z_ma = Zma_A * self.coth_func(Zma_B)
        Z_fractalTLM = Z_ma / n_ma

        # Z_RctW_branch = (1 / (( 1/(Rct + Zw + (1j*ang*Qad)**-1) + 1/Rorr)**-1 ) )
        # A = np.sqrt(R_ion)/ np.sqrt((1j*ang*Cdlp) + Z_RctW_branch)
        # B = np.sqrt(R_ion*(1j*ang*Cdlp) + (R_ion * Z_RctW_branch))
        # Zcoth = A* self.coth_func(B)
        return Ls * (1j * ang) + Z_fractalTLM


class F_fractal_TLMTLM(BaseModel):

    name = "RL-(Rion-TLM(Cdp,TLM(Rct)))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Rorr, Qad, nDL):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""
        # Zw = self.Warburg_Z(ang, Aw, tau)
        Z_branch_me = 1 / (Rorr)
        A_me = np.sqrt(Rct) / np.sqrt((1j * ang * Qad) + Z_branch_me)
        B_me = np.sqrt(Rct * (1j * ang * Qad) + (Rct * Z_branch_me))
        Zcoth_me = A_me * self.coth_func(B_me)
        Z_RctW_branch = nDL * (1 / (Zcoth_me))
        # Zw = self.Warburg_Z(ang, Aw, tau)
        # Z_RctW_branch = (1 / (( 1/(Zw+(1j*ang*Qad)**-1) + 1/Rorr)**-1 ) )
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)

        return Rs + Ls * (1j * ang) + Zcoth


class F_fractal_TLMTLM_noRs(BaseModel):

    name = "Rion-TLM(Cdp,TLM(Rct))"

    def func(self, ang, R_ion, Cdlp, Rct, Rorr, Qad, nDL):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""

        # Zw = self.Warburg_Z(ang, Aw, tau)
        Z_branch_me = 1 / (Rorr)

        A_me = np.sqrt(Rct) / np.sqrt((1j * ang * Qad) + Z_branch_me)
        B_me = np.sqrt(Rct * (1j * ang * Qad) + (Rct * Z_branch_me))
        Zcoth_me = A_me * self.coth_func(B_me)
        Z_RctW_branch = nDL * (1 / (Zcoth_me))
        # Zw = self.Warburg_Z(ang, Aw, tau)
        # Z_RctW_branch = (1 / (( 1/(Zw+(1j*ang*Qad)**-1) + 1/Rorr)**-1 ) )
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)

        return Zcoth


class F_TLMTLM(BaseModel):

    name = "RL-p(TLM(Rad),TLM(Rct))"

    def func(self, ang, Rs, Ls, R_ion, Cdlp, Rct, Rorr, Qad, Rad):
        """R0-L0-TLM(Rion,p(Cdl,Rct-p(Zw,Rorr)))"""
        # Zw = self.Warburg_Z(ang, Aw, tau)
        Z_branch_me = 1 / (Rorr)
        A_me = np.sqrt(Rct) / np.sqrt((1j * ang * Qad) + Z_branch_me)
        B_me = np.sqrt(Rct * (1j * ang * Qad) + (Rct * Z_branch_me))
        Zcoth_me = A_me * self.coth_func(B_me)
        Z_RctW_branch = 1 / (Rad)
        # Zw = self.Warburg_Z(ang, Aw, tau)
        # Z_RctW_branch = (1 / (( 1/(Zw+(1j*ang*Qad)**-1) + 1/Rorr)**-1 ) )
        A = np.sqrt(R_ion) / np.sqrt((1j * ang * Cdlp) + Z_RctW_branch)
        B = np.sqrt(R_ion * (1j * ang * Cdlp) + (R_ion * Z_RctW_branch))
        Zcoth = A * self.coth_func(B)

        return Rs + Ls * (1j * ang) + 1 / (1 / Zcoth_me + 1 / Zcoth)


class F_TLM_Levie_Rct(BaseModel):
    name = "Levie(Rct)"

    def func(self, ang, Cdlp, Rct, L_pore, D_pore, N_pores, rho):
        "Levie(Rct)"
        Z_Rct_branch = Rct
        Zpore_cm2 = 1 / (1 / (Z_Rct_branch + (1j * ang * Cdlp)))
        R_ion = rho / (np.pi * D_pore ** 2)
        Zpore_branch = Zpore_cm2 / (2 * np.pi * D_pore)
        Z_TLM = self.TLM_Levie(R_ion, Zpore_branch, L_pore)
        Ztot = Z_TLM / N_pores
        return Ztot


class F_TLM_Levie_RctQad(BaseModel):
    name = "Levie(Rct-Qad)"

    def func(self, ang, Cdlp, Rct, Qad, L_pore, D_pore, N_pores, rho):
        "Levie(Rct)"
        Z_Rct_branch = Rct + (1j * ang * Qad) ** -1
        Zpore_cm2 = 1 / (1 / (Z_Rct_branch + (1j * ang * Cdlp)))
        R_ion = rho / (np.pi * D_pore ** 2)
        Zpore_branch = Zpore_cm2 / (2 * np.pi * D_pore)
        Z_TLM = self.TLM_Levie(R_ion, Zpore_branch, L_pore)
        Ztot = Z_TLM / N_pores
        return Ztot


class S_Keiser_pore(BaseModel):
    """
    Z 0 = R 0 + 1 / (jwC1 +1/ Z_i)

    Z_i = R_i + 1 / (jwC_i+1 + 1/ Z_i+1)

    Z_i/R = Zstar_i = 1 /N*g_i**2 + 1/ (j * 1/2N * (1/lambda)**2 + g_i+1 + 1/Zstar_i)
    R = l / (k* pi* r**2)
    r: pore_radius and l total_length
    g_i = r_i/r shapefactor r_i radius of disc
    lambda = 0.5 * np.sqrt(k*r / (w*Cdl))
    k: ionic cond.

    """


#  ==== Gerischer VARIATIONS ==== #


class Ger_RGCPE(BaseModel):
    """R0-L0-p(R1-G0,CPE1)"""

    name = "R0-L0-p(R1-G0,CPE1)"
    guesses = OrderedDict(
        {
            "Rs": 20,
            "Ls": 5e-5,
            "Rct": 95,
            "Cdlp": 7e-05,
            "nDL": 0.65,
            "R_G": 2e3,
            "t_G": 20,
        }
    )

    def func(self, ang, Rs, Ls, Rct, Cdlp, nDL, R_G, t_G):
        Z_G = R_G / (np.sqrt(1 + 1j * ang * t_G))
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return (
            Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Z_G)) ** -1
        )


class Ger_RGCPE_C2(BaseModel):
    """R0-L0-p(R1-G0,CPE1)-C2"""

    name = "R0-L0-p(R1-G0,CPE1)-C2"
    guesses = OrderedDict(
        {
            "Rs": 20,
            "Ls": 5e-5,
            "Rct": 95,
            "Cdlp": 7e-05,
            "nDL": 0.65,
            "R_G": 2e3,
            "t_G": 20,
        }
    )

    def func(self, ang, Rs, Ls, Rct, Cdlp, nDL, R_G, t_G, Qad):
        Z_G = R_G / (np.sqrt(1 + 1j * ang * t_G))
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / (Rct + Z_G)) ** -1
            + (1 / (Qad ** 1 * (1j * ang) ** 1))
        )


class Ger_RGC(BaseModel):
    """R0-L0-p(R1-G0,C1)"""

    name = "R0-L0-p(R1-G0,C1)"
    guesses = OrderedDict(
        {"Rs": 20, "Ls": 5e-5, "Rct": 95, "Cdlp": 7e-05, "R_G": 2e3, "t_G": 20}
    )

    def func(self, ang, Rs, Ls, Rct, Cdlp, R_G, t_G):
        Z_G = R_G / (np.sqrt(1 + 1j * ang * t_G))
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang)) + 1 / (Rct + Z_G)) ** -1


class Ger_RGsC(BaseModel):
    """R0-L0-p(R1-G0,C1)"""

    name = "R0-L0-p(R1-Gs0,C1)"
    guesses = OrderedDict(
        {
            "Rs": 20,
            "Ls": 5e-5,
            "Rct": 95,
            "Cdlp": 7e-05,
            "R_G": 2e3,
            "t_G": 20,
            "phi_G": 1,
        }
    )

    def func(self, ang, Rs, Ls, Rct, Cdlp, R_G, t_G, phi_G):
        Z_Gs = R_G / (
            np.sqrt(1 + 1j * ang * t_G) * np.tanh(phi_G * np.sqrt(1 + 1j * ang * t_G))
        )
        # Z_Ws = Aw*np.tanh(np.sqrt(1j*ang*tau))/np.sqrt(1j*ang*tau)
        return Rs + Ls * (1j * ang) + ((Cdlp * (1j * ang)) + 1 / (Rct + Z_Gs)) ** -1

        # Z_Wo = Awo/(np.sqrt(1j*ang*tau)*np.tanh(np.sqrt(1j*ang*tau)))
        # Zw = Aw*(1-1j)/np.sqrt(ang)
        # return Rs + Ls*(1j*ang) + ( (Cdlp*(1j*ang)) + 1/(Rct + Z_Wo) )**-1 + Zw


# ====================         ===============#
#%% === Model Collection ===
class Model_Collection:
    """Iterable collection of all defined costum (lmfit type) models."""

    _bad_models = [
        "Q_RandlesZw",
        "Q_TLM_Warburg",
        "Q_TLM_Warburg_p_Qad",
        "Q_TLM_Warburg_ser_RctQadW",
        "Q_TLM_W_p_QadRorr",
        "Q_TLM_Warburg_Rct_p_Qad_WRorr",
        "Q_TLM_Warburg_Rct_pW_QadRorr",
        "Q_TLM_RctQadWo_",
        "Q_TLM_W_Rside",
        "Q_TLM_p_RctQadW_Rorr",
    ]
    _std_models = [
        "Q_TLM_Warburg_ser_RctQadW",
        "Q_TLM_Warburg_pQadWRorr",
        "Q_TLM_RctQad_p_WRorr",
        "Q_TLM_W_pRct_Qad",
        "Q_TLM_W_pRct_QadRorr",
        "Q_TLM_pRctQad_W",
        "Q_TLM_ser_RctQadW",
    ]
    _skip_models = [
        "Q_TLM_RctW",
        "Q_TLM_Rct_pQadW",
        "Q_TLM_QadW",
        "Q_TLM_pQad_W",
        "Q_TLM_W_pQadRct",
        "Q_TLM_p_RctQadW_Rorr",
    ]
    # 'Q_TLM_Rct','Q_TLM_ser_RctQad'
    # "Q_TLM_W_p_QadRorr"
    # 'Q_TLM_Qdl_RctQadW','Q_TLM_Qdl_RctW_pQadRorr']
    # 'Q_TLM_Qdl_WRctQadW']-*
    _Randles_models = [
        "Q_Randles_RctZw",
        "Q_Randles_RctZwQad",
        "Q_Randles_RctZw_pQad_Rorr",
    ]
    # 'Q_RandlesZw','Q_RandlesZwQad']
    # 'Q_TLM_RctW_pQadR_Rorr']
    # 'Q_TLM_Warburg_pQadW_W2orrRorr',,'Q_TLM_Warburg_pQadW_QorrRorr',
    # 'Q_TLM_pRctQadW_Rorr']

    def __init__(self, model_prefixes=["Q"], standard_models=True, _startswith=""):
        self.model_prefixes = model_prefixes
        self._startswith = _startswith
        # self._bad_models = _bad_models
        self._skipped_models = set(self._bad_models + self._skip_models)
        self.std_model_selection = bool(standard_models)
        # if any(m[0].startswith(_pre) for _pre in model_prefixes)
        self.validation_inspect_models()
        self.model_selection = [m for m in self.valid_models]
        if self.std_model_selection:
            self.model_selection = [
                m
                for m in self.valid_models
                if (
                    (
                        m().classname in self._std_models
                        and not m().classname in self._skipped_models
                    )
                )
            ]
        if self._startswith:
            self.model_selection = [
                m
                for m in self.valid_models
                if (
                    m().classname.startswith(self._startswith)
                    and not m().classname in self._skipped_models
                )
            ]
        # if len([i for i in m.name if i == 'R']) == 2

        self._bak_models = [
            M_2CPE,
            M_2CPE_W,
            M_CPE_WCPE,
            M_RC_CPE_W,
            M_RC_WsCPE,
            N_RWoCPE_W,
            M_RC_WsCPE,
            M_RC_WoCPE,
        ]

        self._nonused_models = [N_RWoC_W]
        self.assign_colors_to_mod_inst()
        self.add_standard_init_params()
        self.add_model_names_var_names()
        # self.lmfit_models = [(m(), ', '.join([str(i) for i in self.cmap_set(n)]))
        # for n,m in enumerate(self.model_selection)]

    def assign_colors_to_mod_inst(self):
        self.cmap_set = plt.get_cmap(
            "Dark2" if not len(self.model_selection) > 10 else "tab20"
        )
        _mod_inst = []
        for n, m in enumerate(self.model_selection):
            _m_inst = m()
            _m_inst.color = ", ".join([str(i) for i in self.cmap_set(n)])
            # _m_inst._funcname = str(m).split('__main__.')[-1][:-2]
            _m_inst._lenpars = len(_m_inst.model.param_names)
            _mod_inst.append(_m_inst)
        _mod_inst = sorted(_mod_inst, key=lambda x: x._lenpars)
        self.lmfit_models = _mod_inst

    def add_standard_init_params(self):
        self.standard_init_params = Parameters()
        self.standard_init_params.add_many(*BaseModel._params_guesses_base)

    def add_model_names_var_names(self):
        _modvars = {i.model.name: i.model.param_names for i in self.lmfit_models}
        self.modpars = _modvars

    def get_df_models_parameters(self):
        _EIS_models = pd.DataFrame(
            [
                (i.model.name, len(i.model.param_names), ", ".join(i.model.param_names))
                for i in self.lmfit_models
            ],
            columns=["Model_EEC", "model_lenpars", "model_parnames"],
        )
        return _EIS_models

    # def select(self,name):
    #     if

    def validation_inspect_models(self):

        self._inspect_models = [
            i
            for i in inspect.getmembers(sys.modules[__name__], inspect.isclass)
            if issubclass(i[1], BaseModel)
        ]
        _validmods, _badfuncs = [], []
        for _nm, m in self._inspect_models:
            # if issubclass(m,BaseModel):
            try:
                _inst = m()
                if hasattr(_inst, "model"):
                    _validmods.append(m)
                # if ;Basestr(type(_inst))
                # assert ha
            except Exception as e:
                _badfuncs.append(m)
                # print(_nm,e)
        self.valid_models = _validmods
        self.bad_funcs = _badfuncs

    def __iter__(self):
        for mod_inst in self.lmfit_models:
            yield mod_inst
        # self.params_set = set([a for m in self.lmfit_models for a in m[1].param_names])

    #    lmfit_models = [Model(i.func,name=i.name) for i in model_selection]
    def __repr__(self):
        return (
            f"Model Collections, {len(self.model_selection)} models from prefixes {self.model_prefixes}: "
            + "\n\t- "
            + "\n\t- ".join(
                [
                    f'{i.name} \t "{i.classname}", {i._lenpars}'
                    for i in self.lmfit_models
                ]
            )
        )


#    def colors(self, cmap = 'tab20'):
#        self.lmift_models_color = ((*i, plt.get_cmap('tab20')(n)) for n,i in enumerate(self.lmfit_models))
# {mod_inst : {'color'Angular' mod_inst.ang
# ====================   ^^^    ===============#
#%% older stuff
