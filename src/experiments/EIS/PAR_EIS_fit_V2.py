# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 11:17:23 2018

@author: User
"""
import lmfit
from lmfit.models import VoigtModel
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cmath
from cmath import phase

matplotlib.rcParams.update({"font.size": 14})
plt.rcParams["interactive"]  # or: plt.isinteractive()
import cycler

# from FolderOr import *
from lmfit import minimize, Parameters, Parameter, Model, report_errors

# import PyOrigin
# import FolderOps
# try:
#    from .FolderOrganizing import FolderOps
#    print('imported FolderOrganizing.FolderOps')
# except:
#    try:
#        import FolderOps
#    except:
#        print('PAR EIS FIT error imported Folderops',__name__)
#    pass
from collections import namedtuple
import logging

try:
    import impedance
except:
    print("no import impedancepy")
from scipy.optimize import least_squares

import FileHelper

#%%
EvRHE = "E_AppV_RHE"
# Gets or creates a logger
def rmse(a, b):
    """
    Code copied from : https://github.com/ECSHackWeek/impedance.py [2019-12-02]
    A function which calculates the root mean squared error
    between two vectors.
    Notes
    ---------
    .. math::
        RMSE = \\sqrt{\\frac{1}{n}(a-b)^2}
    """

    n = len(a)
    return np.linalg.norm(a - b) / np.sqrt(n)


def typeChecker(p, f, name, length):
    assert isinstance(p, list), "in {}, input must be of type list".format(name)
    for i in p:
        assert isinstance(
            i, (float, int, np.int32, np.float64)
        ), "in {}, value {} in {} is not a number".format(name, i, p)
    for i in f:
        assert isinstance(
            i, (float, int, np.int32, np.float64)
        ), "in {}, value {} in {} is not a number".format(name, i, f)
    assert len(p) == length, "in {}, input list must be length {}".format(name, length)
    return


def element_metadata(num_params, units):
    """decorator to store metadata for a circuit element
    Parameters
    ----------
    num_params : int
        number of parameters for an element
    units : list of str
        list of units for the element parameters
    """

    def decorator(func):
        def wrapper(p, f):
            typeChecker(p, f, func.__name__, num_params)
            return func(p, f)

        wrapper.num_params = num_params
        wrapper.units = units
        wrapper.__name__ = func.__name__
        wrapper.__doc__ = func.__doc__

        return wrapper

    return decorator


# def StartLogging():
#    log_fn = FolderOps.FindExpFolder('VERSASTAT').DestDir.joinpath('EC_EIS_logger.log')
#    logging.basicConfig(filename=log_fn, filemode='w', level=logging.DEBUG,format='%(asctime)s %(message)s')
#    logging.info('Starting PAR EIS fit log...')


def EEC_Bandarenka_Ads(ang, Rs, Cdlp, nDL, Rct, Qad, nAd):
    """Bondarenko_Langmuir_2011, Fig6b, adsorption of species"""
    #    return Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Qad*(1j*ang)**-0.5)))))
    return Rs + (
        1
        / (
            Cdlp ** 1 * (1j * ang) ** nDL
            + 1 / (Rct + 1 / ((Qad ** 1 * (1j * ang) ** nAd)))
        )
    )


def EEC_Randles_RQRQ(ang, Rs, Cdlp, nDL, Rct, Qad, nAd):
    """Bondarenko_Langmuir_2011, Fig6b, adsorption of species"""
    return (
        Rs
        + (
            (1 / (Cdlp * (1j * ang) ** nDL) ** -1)
            + (1 / (Rct + (Qad * (1j * ang) ** nAd) ** -1))
        )
        ** -1
    )


def model_ORR(params, ang, Zdata):
    """Model Bondarenko_Langmuir_2011, Fig6b"""
    Rs = params["Rs"].value
    Cdlp = params["Cdlp"].value
    nDL = params["nDL"].value
    Rct = params["Rct"].value
    Qad = params["Qad"].value
    Zfit = EEC_ORR(ang, Rs, Cdlp, n, Rct, Aw)
    #    (Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-n)+(1/(Ra+Aw*(1j*ang)**-0.5))))))
    diff = Zfit - Zdata
    return diff.view(np.float)


#
# def EEC_ORRpOx(ang,Rs,Cdlp,nDL,Rct,Aw,Rad,Qad,nAd,Rorr):
#    '''Other model with Faradaic and 2 CPEs'''
#    return Rs+ (1/(Cdlp*(1j*ang)**nDL+ 1/( Rct+ 1/( (Qad*(1j*ang)**nAd)*Rorr**-1))))
##Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Aw*(1j*ang)**-0.5+ Rad+ 1/(Qad**-1*(1j*ang)**-nAd+Rorr**-1))))))


def EEC_ORRpOx(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
    """Other model with Faradaic and 2 CPEs"""
    return Rs + (
        1
        / (
            Cdlp ** 1 * (1j * ang) ** nDL
            + 1 / (Rct + 1 / ((Qad ** 1 * (1j * ang) ** nAd) + Rorr ** -1))
        )
    )


# Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Aw*(1j*ang)**-0.5+ Rad+ 1/(Qad**-1*(1j*ang)**-nAd+Rorr**-1))))))

""" Singh et al. jES(2015) studied EIS of the ORR on a RDE, effect of ionomer content and carbon support.
 no Ionomer
 Pt black """


def EEC_Singh2015_3RQ(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr, R3, Q3, n3):
    """Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
    return (
        Rs
        + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / Rct) ** -1
        + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
        + ((1 / (1 / (Q3 ** 1 * (1j * ang) ** n3))) + 1 / R3) ** -1
    )


def EEC_Singh2015_RQRQR(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
    """Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
    return (
        Rs
        + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / Rct) ** -1
        + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
    )


def EEC_Singh2015_RQRWR(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
    """Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)
    with only Warburg as set nAd=0.5 (fixed) during fitting"""
    return (
        Rs
        + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / Rct) ** -1
        + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
    )


# (1/(Qad**1*(1j*ang)**nAd)+Rorr**-1)


def EEC_Singh2015_RQRQ(ang, Rs, Cdlp, nDL, Rct, Qad, nAd):
    """Other model with 2 (R/Q-CPEs)from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
    return (
        Rs
        + (Rct ** -1 + (Qad ** 1 * (1j * ang) ** nAd) ** -1) ** -1
        + (Cdlp ** 1 * (1j * ang) ** nDL) ** -1
    )


# Rs + (1/ (Cdlp**1*(1j*ang)**nDL + 1/Rct) + (1/(Qad**1*(1j*ang)**nAd))

# !!! Switch position of Rorr with Rct
# Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Aw*(1j*ang)**-0.5+ Rad+ 1/(Qad**-1*(1j*ang)**-nAd+Rorr**-1))))))
def EEC_Singh2015_RQR(ang, Rs, Cdlp, nDL, Rct):
    """Other model with 2 (R/Q-CPEs)from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
    return Rs + (Rct ** -1 + (Cdlp ** 1 * (1j * ang) ** nDL) ** -1) ** -1


# Rs + ((1/(Cdlp*(1j*ang)**nDL)**-1)+1/Rct)**-1
# Rs + ((1/(1/(Cdlp**1*(1j*ang)**nDL))) + 1/Rct)**-1
#           Rs + (1/Cdlp**1*(1j*ang)**nDL+ 1/Rct)


def EIS_ORR_fit(PAR_EIS_data):
    PAR_EIS_data["-Z Imag"] = -1 * PAR_EIS_data["Z Imag"]
    Zre, Zim = PAR_EIS_data["Z Real"], PAR_EIS_data["Z Imag"]
    gr.plot(x="Z Real", y="-Z Imag", kind="scatter")


def ORR_model_7_pars(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
    #    ang,Rs,Cdlp,nDL,Rct,Qad,nAd,Rorr = **args
    EEC_choices = {
        "EEC_ORRpOx": Rs
        + (
            1
            / (
                Cdlp ** 1 * (1j * ang) ** nDL
                + 1 / (Rct + 1 / ((Qad ** 1 * (1j * ang) ** nAd) + Rorr ** -1))
            )
        ),
        "EEC_Singh2015_RQRQR": Rs
        + (1 / Cdlp ** 1 * (1j * ang) ** nDL + 1 / Rct)
        + (1 / (Qad ** 1 * (1j * ang) ** nAd) + Rorr ** -1),
        "EEC_Singh2015_RQRQ": Rs
        + (1 / Cdlp ** 1 * (1j * ang) ** nDL + 1 / Rct)
        + (1 / (Qad ** 1 * (1j * ang) ** nAd)),
    }
    return EEC_choices


def model_ORRpOx(params, ang, Zdata):
    Rs = params["Rs"].value
    Cdlp = params["Cdlp"].value
    nDL = params["nDL"].value
    Rct = params["Rct"].value
    #    Rad, Aw = params['Aw'].value, params['Rad'].value
    Qad, nAd, Rorr = params["Qad"].value, params["nAd"].value, params["Rorr"].value
    Zfit = EEC_ORRpOx(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr)
    #    (Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-n)+(1/(Ra+Aw*(1j*ang)**-0.5))))))
    diff = Zfit - Zdata
    return diff.view(np.float)


def EIS_simpleRRC(ang, Rs, Cdlp, Rc, nDL):
    return Rs + (1 / (Cdlp ** 1 * (1j * ang) ** nDL + 1 / Rc))


def model_EEC(params, ang, Zdata):
    Rs = params["Rs"].value
    Cdlp = params["Cdlp"].value
    Rc = params["Rct"].value
    nDL = params["nDL"].value
    Zfit = EIS_simpleRRC(ang, Rs, Cdlp, Rc, nDL)
    #    (Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-n)+(1/(Ra+Aw*(1j*ang)**-0.5))))))
    diff = Zfit - Zdata
    return diff.view(np.float)


def EEC_models_index():
    model_index = [
        (1, Model(EEC_Singh2015_RQR, name="Singh2015_RQR")),
        (2, Model(EEC_Singh2015_RQRQ, name="Singh2015_RQRQ")),
        (3, Model(EEC_Singh2015_RQRQR, name="Singh2015_RQRQR")),
        (4, Model(EEC_ORRpOx, name="Bandarenka_2011_RQRQR")),
        (6, Model(EEC_Singh2015_RQRQR, name="Singh2015_RQRWR")),
        (7, Model(EEC_Randles_RQRQ, name="Randles_RQRQ")),
        (8, Model(EEC_Singh2015_3RQ, name="Singh2015_R3RQ")),
    ]
    return model_index


def create_example_series():
    DestDirTop = FindExpFolder("VERSASTAT").DestDir.joinpath("EIS_example")
    DestDirTop.mkdir(parents=True, exist_ok=True)
    #    default_pars = {'Rs' : 60, 'Cdlp' : 2E-04,'nDL' : 1, 'Rct' : 100,'Qad' : 1E-03,'nAd' : 0.5,'Rorr' : 7E5}
    default_pars = {
        "Rs": 40,
        "Cdlp": 2e-04,
        "nDL": 0.9,
        "Rct": 100,
        "Qad": 1e-03,
        "nAd": 0.7,
        "Rorr": 7e5,
    }
    ranges = {
        "Rct": [10, 50, 500, 1e3, 5e3],
        "Rs": [10, 20, 40, 80],
        "nDL": [0.5, 0.7, 0.8, 0.9, 1],
        "nAd": [0.5, 0.7, 0.8, 0.9, 1],
        "Qad": [1e-05, 5e-05, 1e-04, 5e-04, 1e-03],
        "Cdlp": [1e-06, 5e-06, 1e-05, 5e-05],
        "Rorr": [1, 1e1, 1e2, 1e3, 5e4],
    }

    DestDir = DestDirTop.joinpath(
        "test_{0}".format(np.abs(hash(default_pars.values())))
    )
    DestDir.mkdir(parents=True, exist_ok=True)

    for par in ranges.keys():
        outpar, outdata = [], []

        color = plt.cm.viridis(np.linspace(0, 1, len(ranges[par])))
        plt.rcParams["axes.prop_cycle"] = cycler.cycler("color", color)
        plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20c.colors)
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12))
        ax, ax2 = axes[0, 0], axes[0, 1]
        ax3, ax4 = axes[1, 0], axes[1, 1]
        plt.suptitle(
            ", ".join(
                "{0} = {1}".format(k, v)
                for k, v in default_pars.items()
                if k not in par
            )
        )
        #        matplotlib.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
        #        hexcolor = map(lambda rgb:'#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255),
        #               tuple(color[:,0:-1]))

        for val, cl in zip(ranges[par], color):
            dfile = DestDir.joinpath("{0}_{1}_data.xlsx".format(par, val))
            #            print(par,val)
            od, oprs, mod, mod2 = make_example_plot(
                par=par, val=val, default=default_pars, single_plot=0
            )
            od["color"] = val
            od.to_excel(dfile)
            #            matplotlib.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)
            outpar.append(oprs)
            outdata.append(od)
            #            od.plot(x='INIT2_Zre',y='INIT2_-Zim',kind='scatter',s=80, c=[cl]*len(od),ax=ax,label='{0} = {1}'.format(par,val))
            ax.plot(
                od["INIT2_Zre"].values,
                od["INIT2_-Zim"].values,
                lw=3,
                label="{0} = {1}".format(par, val),
            )
            ax.plot(
                od["INIT1_Zre"].values,
                od["INIT1_-Zim"].values,
                lw=3,
                label=" ",
                ls="dotted",
                alpha=0.85,
            )
            if val == ranges[par][-1]:
                ax2.plot(
                    od["INIT2_Yre"].values,
                    od["INIT2_Yim"].values,
                    lw=3,
                    label=mod2.name,
                )
                ax2.plot(
                    od["INIT1_Yre"].values,
                    od["INIT1_Yim"].values,
                    lw=3,
                    ls="dotted",
                    alpha=0.85,
                    label=mod.name,
                )
            else:
                ax2.plot(od["INIT2_Yre"].values, od["INIT2_Yim"].values, lw=3)
                ax2.plot(
                    od["INIT1_Yre"].values,
                    od["INIT1_Yim"].values,
                    lw=3,
                    ls="dotted",
                    alpha=0.85,
                )

            ax3.plot(
                od["Frequency(Hz)"].values,
                od["INIT2_Zabs"].values,
                lw=3,
                label="{0} = {1}".format(par, val),
            )
            ax3.plot(
                od["Frequency(Hz)"].values,
                od["INIT1_Zabs"].values,
                lw=3,
                ls="dotted",
                alpha=0.85,
            )

            ax4.plot(
                od["Frequency(Hz)"].values,
                od["INIT2_Zangle"].values,
                lw=3,
                label="{0} = {1}".format(par, val),
            )
            ax4.plot(
                od["Frequency(Hz)"].values,
                od["INIT1_Zangle"].values,
                lw=3,
                ls="dotted",
                alpha=0.85,
            )
        #            od.plot(x='INIT2_Yre',y='INIT2_Yim',kind='scatter',s=80, c=[cl]*len(od),ax=ax2,alpha=0.5)

        all_oprs = pd.concat(
            [pd.DataFrame(i, index=[0]) for i in outpar], ignore_index=True, sort=False
        )
        all_data = pd.concat([i for i in outdata], ignore_index=True, sort=False)
        """ax2: Admittance Plots"""
        ax.set_xlim(0, all_data["INIT2_Zre"].max())
        ax.set_xlabel("Z real")
        ax.set_ylim(0, all_data["INIT2_Zre"].max())
        ax.set_ylabel("Z imag")
        legend0 = ax.get_legend_handles_labels()

        ax.legend(loc="upper left", ncol=1, fontsize=12)
        ax.grid(True), ax2.grid(True)
        ax2.set_xlim(0, all_data["INIT2_Yre"].max())
        ax2.set_xlabel("Y real")
        ax2.set_ylabel("Y imag")
        ax2.set_ylim(0, all_data["INIT2_Yre"].max())
        ax2.set_label([ax2.get_legend_handles_labels()][-2:])
        ax2.legend(loc="upper right")

        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_ylabel("log(|Z| / Ohm)")
        ax3.set_xlabel("log(Frequency / Hz)")

        ax4.set_xscale("log")
        ax4.set_ylabel("log(phase agle / degree)")
        ax4.set_xlabel("log(Frequency / Hz)")
        ax4.grid(True), ax3.grid(True)
        plt.savefig(DestDir.joinpath("{0}.png".format(par)), dpi=300)
        plt.show()
        print("Figure saved to: {0}".format(DestDir))
        plt.close()
        #        all_data.plot(x='INIT2_Zre',y='INIT2_-Zim',kind='scatter',c='color',lw=2.5,label='{0} = {1}'.format(par,val))
        #        all_data.plot(x='INIT2_Yre',y='INIT2_Yim',kind='scatter',c='color',lw=2.5,label='{0} = {1}'.format(par,val))
        all_oprs.to_excel(DestDir.joinpath("{0}.xlsx".format(par)))


def mktest(**kwargs):
    print(kwargs)
    return kwargs


def make_example_plot(**kwargs):
    freq = np.logspace(-0.5, 4, num=50, dtype="float")
    ang = freq * 2 * np.pi

    mod, mod2 = Model(EEC_Singh2015_RQRQR, name="EEC_Singh2015_RQRQR"), Model(
        EEC_ORRpOx, name="EEC_ORRpOx"
    )
    #%%
    params = Parameters()
    params.add("Rs", value=30, min=0.1, max=1e3)
    params.add("Cdlp", value=1e-06, min=1e-08, max=1e-1)
    params.add("nDL", value=1, min=0, max=1)
    params.add("Rct", value=100, min=1e-3, max=1e6, vary=True)
    #    params.add('Aw', value= 3000, min=5.0, max=1E6)
    #    params.add('Rad', value= 5E5, min=0.01, max=1E7)
    params.add("Qad", value=1e-03, min=1e-07, max=1e-02)
    params.add("nAd", value=1, min=0, max=1)
    params.add("Rorr", value=1e3, min=0.01, max=1e8)

    #    kwargs = {'Rs' : 20, 'Rorr' : 3000}
    if kwargs:
        #        print(kwargs)
        keys = kwargs.keys()
        #        print('par' in  and 'val 'in
        if "default" in keys:
            for dk in kwargs["default"].keys():
                params[dk].value = kwargs["default"][dk]
        print(keys)
        if "par" in keys or "val " in keys:
            par = kwargs["par"]
            val = kwargs["val"]
            params[par].value = val
            print(kwargs, par, val)
        #        for i in kwargs:
        #            print(kwargs)
        #            if i in ['Rs','Cdlp','nDL','Rct','Qad','nAd','Rorr']:
        #                params[i].value = kwargs[i]

        if "single_plot" in kwargs:
            single_plot = kwargs["single_plot"]
        else:
            single_plot = []

    #    outData = outData.assign(**{'Angular' : ang, 'DATA_Zre' : Zre, 'DATA_Zim' : Zim, 'DATA_-Zim' : -1*Zim,
    #                                'DATA_Z' : Zdata,'DATA_Zmod' : abs(Zdata),
    #                                'DATA_Zphase' : [phase(i) for i in Zdata] ,
    #                                'DATA_Y' : 1/Zdata,'DATA_Yre' : (1/Zdata).real,'DATA_Yim' : (1/Zdata).imag})
    #
    #%%
    outData = pd.DataFrame({"Frequency(Hz)": freq})
    init = mod.eval(params, ang=ang)
    #    print(params.pretty_print())
    #    ax3, ax4  = axes[1,0], axes[1,1]
    #    ax5, ax6 = axes[2,0] ,axes[2,1]
    #    fig.suptitle('%s, %s: \n$\mathrm{E_{DC}\/ =\/ %.2f\/ V_{RHE} \/ in \/%s\/at\/pH\/=\/%.0f}$' %(EISgr['SampleID'].unique()[0],
    #    Path(EISgr['File'].unique()[0]).stem,E_dc ,EISovv['Electrolyte'].unique()[0],EISovv['pH'].unique()[::]))
    dataC, fitC, extraC, initC = "tab:blue", "tab:red", "gold", "gray"
    #        outData.plot(x='FIT_Zre',y='FIT_-Zim',kind='line',c='orange',ax=ax)
    """ax: Impedance (Nyquist) Plots"""
    #     mod,mod2 = Model(EEC_ORR),Model(EEC_ORRpOx) 'Init_Z1' : init1,'Init_Z2' : init2,
    add_dict = {
        "INIT2_Yim": (1 / init2).imag,
        "INIT2_Yre": (1 / init2).real,
        "INIT1_Zre": init1.real,
        "INIT1_-Zim": -1 * init1.imag,
        "INIT2_Zre": init2.real,
        "INIT2_-Zim": -1 * init2.imag,
        "INIT2_Zphase": np.array([phase(i) for i in init2]),
        "INIT1_Zphase": np.array([phase(i) for i in init1]),
        "INIT2_Zangle": np.angle(init2, deg=True),
        "INIT1_Zangle": np.angle(init1, deg=True),
        "INIT2_Zabs": np.abs(init2),
        "INIT1_Zabs": np.abs(init1),
        "INIT1_Yim": (1 / init1).imag,
        "INIT1_Yre": (1 / init1).real,
    }
    #    add_dict(list(add_dict.keys())[0])
    outData = outData.assign(**add_dict)
    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 8))
        ax, ax2 = axes[0], axes[1]
        outData.plot(
            x="INIT1_Zre",
            y="INIT1_-Zim",
            kind="scatter",
            c=extraC,
            lw=2.5,
            ax=ax,
            label="FIT: {0}".format(mod.name),
        )
        #        outData.plot(x='INIT2_Zre',y='INIT2_-Zim',kind='scatter',marker='*',c=fitC,lw=2.5,ax=ax,label='FIT : {0}'.format(mod2.name))
        outData.plot(
            x="INIT2_Zre",
            y="INIT2_-Zim",
            kind="line",
            c=fitC,
            lw=2.5,
            ax=ax,
            label="FIT : {0}".format(mod2.name),
        )

        #    outData.plot(x='DATA_Zre',y='DATA_-Zim',kind='scatter',c=dataC,s=150,ax=ax)
        #    ax.plot(outData['FIT2_Zre'].values,outData['FIT2_Zre'].values,ls=':',c='k',lw=1)
        #    outData.plot(x='FIT1_Zre',y='FIT1_Zre',kind='line',ls='.',c='k',lw=10,ax=ax)
        #    ax.axis('equal')

        ax.set_xlim(0, init2.real.max())
        ax.set_ylim(0, init2.real.max())
        ax.grid(True)
        """ax2: Admittance Plots"""
        #    outData.plot(x='FIT2_Yre',y='FIT2_Yim',kind='line',c=fitC,lw=2.5,ax=ax2)
        #    outData.plot(x='DATA_Yre',y='DATA_Yim',kind='scatter',c=dataC,s=150,ax=ax2)
        #    outData.plot(x='FIT1_Yre',y='FIT1_Yim',kind='line',c=extraC,lw=2.5,ax=ax2)
        outData.plot(
            x="INIT1_Yre",
            y="INIT1_Yim",
            kind="scatter",
            c=extraC,
            lw=2.5,
            ax=ax2,
            alpha=0.5,
        )
        #        outData.plot(x='INIT2_Yre',y='INIT2_Yim',kind='scatter',c=fitC,lw=2.5,ax=ax2,alpha=0.5)
        outData.plot(
            x="INIT2_Yre", y="INIT2_Yim", kind="line", c=fitC, lw=2.5, ax=ax2, alpha=0.5
        )

        ax2.set_xlim(0, (1 / init2).real.max())
        ax2.set_ylim(0, (1 / init2).real.max())
        #    ax2.autoscale(True)
        #    ax2.axis('equal')
        ax2.grid(True)

        props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

        #    if text_out:
        text_Zmin = """ Z_Re min = {0:.3} \n Z_Re max = {1:.3}\n\n Y_Re min = {2:.3} \n Y_Im max = {3:.3}""".format(
            init2.real.min(),
            init2.real.max(),
            (1 / init2).real.min(),
            (1 / init2).real.max(),
        )
        ax.text(
            0.2,
            -0.3,
            text_Zmin,
            transform=ax.transAxes,
            fontsize=15,
            verticalalignment="center",
            bbox=props,
        )
        text_pars = ""
        for a in params:
            text_pars += "\n{0} = {1:.4}".format(a, float(params[a].value))

        text_pars += "\n"

        ax2.text(
            1.2,
            -0.30,
            text_pars,
            transform=ax.transAxes,
            fontsize=15,
            verticalalignment="center",
            bbox=props,
        )
    #%%
    #    print(params.pretty_print())
    return outData, params.valuesdict(), mod, mod2


# ax.legend(True)


# freq,Zdata,Z_linKK,res_real,res_imag = EIS_data_valid['Frequency(Hz)'].values,EIS_data_valid.DATA_Z.values,EIS_data_valid.linKK_Z.values,EIS_data_valid.linKK_resRe.values,EIS_data_valid.linKK_resIm.values
def plot_linKK(EIS_data_valid, EIS_data, res_limit=0.06, res_scale=1, *args, **kwargs):
    EIS_invalid = EIS_data.query("Valid == False")
    freq, Zdata, Z_linKK, res_real, res_imag = (
        EIS_data_valid["Frequency(Hz)"].values,
        EIS_data_valid.DATA_Z.values,
        EIS_data_valid.linKK_Z.values,
        EIS_data_valid.linKK_resRe.values,
        EIS_data_valid.linKK_resIm.values,
    )
    fig = plt.figure(figsize=(10, 9))
    gs = fig.add_gridspec(3, 2)
    ax1 = fig.add_subplot(gs[:2, 0])
    ax3 = fig.add_subplot(gs[:2, 1])
    ax2 = fig.add_subplot(gs[2, :])
    #    for key, value in kwargs.items():
    #        print("{0}".format(key))
    sup1, sup2 = [], []
    if "meta" in kwargs.keys():
        meta_EIS = kwargs["meta"]
        E_Vrhe = meta_EIS["E_AppV_RHE"]
        E_mOCP = float(meta_EIS["Measured_OCP"]) * 1e-3 + float(meta_EIS["RHE_OCP"])
        basenm = meta_EIS["BaseName"]
        sup1 = "{:s} at {:.2f} Vrhe with OCP {:.2f}".format(
            basenm, E_Vrhe, float(E_mOCP)
        )
    if "linKK" in kwargs.keys():
        meta_linKK = kwargs["linKK"]
        sup2 = "Number of RC elements {} with c = {:.3f}".format(*meta_linKK)
    if sup1 and sup2:
        plt.suptitle(sup1 + "\n" + sup2)
    if "type_res" in kwargs.keys():
        set_res_type = kwargs["type_res"]
    else:
        set_res_type = "Z"

    #    print(kwargs.keys())
    # plot original data
    plot_nyquist(ax1, freq, Zdata, fmt="s")
    ax1.scatter(
        EIS_invalid.DATA_Zre,
        -1 * EIS_invalid.DATA_Zim,
        c="grey",
        s=40,
        alpha=0.6,
        marker="x",
    )
    # plot measurement model
    plot_nyquist(ax1, freq, Z_linKK, fmt="-", scale=1, units="\Omega")
    # plot original data
    plot_nyquist(ax3, freq, (Zdata ** -1).real + 1j * -1 * (Zdata ** -1).imag, fmt="s")
    # plot measurement model
    plot_nyquist(
        ax3,
        freq,
        (Z_linKK ** -1).real + 1j * -1 * (Z_linKK ** -1).imag,
        fmt="-",
        scale=1,
        units="S",
    )
    ax3.scatter(
        EIS_invalid.DATA_Yre,
        EIS_invalid.DATA_Yim,
        c="grey",
        s=40,
        alpha=0.6,
        marker="x",
    )
    ax1.legend(["Data", "Lin-KK model"], loc="best", fontsize=12)

    # Plot residuals
    ax2.plot(freq, res_real * res_scale, "-", label=r"$\Delta_{\mathrm{Real}}$")
    ax2.plot(freq, res_imag * res_scale, "-", label=r"$\Delta_{\mathrm{Imag}}$")
    if "KK_valid_limit" in kwargs.keys():
        ax2.plot(
            freq,
            len(freq) * [kwargs["KK_valid_limit"][0] * res_scale],
            "--",
            alpha=0.7,
            c="orange",
        )
        ax2.plot(
            freq,
            len(freq) * [kwargs["KK_valid_limit"][0] * -1 * res_scale],
            "--",
            alpha=0.7,
            c="orange",
        )

        ax2.plot(
            freq,
            len(freq) * [kwargs["KK_valid_limit"][1] * res_scale],
            "--",
            alpha=0.7,
            c="lightblue",
        )
        ax2.plot(
            freq,
            len(freq) * [kwargs["KK_valid_limit"][1] * -1 * res_scale],
            "--",
            alpha=0.7,
            c="lightblue",
        )

    #        ax2.plot(freq, kwargs['KK_valid_limit'][1]*res_scale, '--', label=r'$\Delta_{\mathrm{Real}}$',alpha=0.7)
    #        ax2.plot(freq, kwargs['KK_valid_limit'][1]*-1*res_scale, '--', label=r'$\Delta_{\mathrm{Real}}$', alpha=0.7)
    ax2.set_title("Lin-KK Model Error for {}".format(set_res_type), fontsize=14)

    ax2.tick_params(axis="both", which="major", labelsize=12)
    ax2.set_ylabel("$\Delta$ $(\%)$", fontsize=14)
    ax2.set_xlabel("$f$ [Hz]", fontsize=14)
    ax2.set_xscale("log")

    if np.abs(res_real).max() < 0.05 and np.abs(res_imag).max() < 0.05:
        res_ylim = 0.05
    elif np.abs(res_real).max() < 0.1 and np.abs(res_imag).max() < 0.1:
        res_ylim = 0.1
    elif np.abs(res_real).max() < 0.2 and np.abs(res_imag).max() < 0.2:
        res_ylim = 0.2

    ax2.set_ylim(-res_ylim * res_scale, res_ylim * res_scale)

    ax2.legend(loc="best", fontsize=14, ncol=2)
    #    vals = ax2.get_yticks()
    #    ax2.set_yticklabels(['{:.0%}'.format(x) for x in vals])
    plt.tight_layout()
    #    plt.show()
    if "meta" in kwargs.keys() and "save_target" in kwargs.keys():
        plt.savefig(Path(kwargs["save_target"]), dpi=100, bbox_inches="tight")
    plt.close()


#%%
# === DATA VALIDATION Kramers-Kronig ====
class KramersKronigValidation:
    def __init__(self):
        pass

    #        self.f = f
    #        self.Z = Z

    def assigncols(df_in, Z_linKK, res_real, res_imag, *args, **kwargs):

        linKK = pd.DataFrame(
            {
                "linKK_Z": Z_linKK,
                "linKK_Zreal": Z_linKK.real,
                "linKK_Zimag": Z_linKK.imag,
                "linKK_Y": Z_linKK ** -1,
                "linKK_Yreal": (Z_linKK ** -1).real,
                "linKK_Yimag": (Z_linKK ** -1).imag,
                "linKK_resRe": res_real,
                "linKK_resIm": res_imag,
            },
            index=df_in.index,
        )
        if "suffix" in kwargs.keys():
            suffix = kwargs["suffix"]
            linKK.rename(
                columns=dict(zip(linKK.columns, [i + suffix for i in linKK.columns])),
                inplace=True,
            )

        return pd.concat([df_in, linKK], axis=1)

    def linKK(f, Z, c=0.85, max_M=50, type_res="Z"):
        """Code copied from https://github.com/ECSHackWeek/impedance.py/blob/master/impedance/validation.py [2019-12-02]
        A method for implementing the Lin-KK test for validating linearity [1]
        Parameters
        ----------
        f: np.ndarray
            measured frequencies
        Z: np.ndarray of complex numbers
            measured impedances
        c: np.float
            cutoff for mu
        max_M: int
            the maximum number of RC elements
        Returns
        -------
        mu: np.float
            under- or over-fitting measure
        residuals: np.ndarray of complex numbers
            the residuals of the fit at input frequencies
        Z_fit: np.ndarray of complex numbers
            impedance of fit at input frequencies
        Notes
        -----
        The lin-KK method from Schönleber et al. [1] is a quick test for checking
        the
        validity of EIS data. The validity of an impedance spectrum is analyzed by
        its reproducibility by a Kramers-Kronig (KK) compliant equivalent circuit.
        In particular, the model used in the lin-KK test is an ohmic resistor,
        :math:`R_{Ohm}`, and :math:`M` RC elements.
        .. math::
            \\hat Z = R_{Ohm} + \\sum_{k=1}^{M} \\frac{R_k}{1 + j \\omega \\tau_k}
        The :math:`M` time constants, :math:`\\tau_k`, are distributed
        logarithmically,
        .. math::
            \\tau_1 = \\frac{1}{\\omega_{max}} ; \\tau_M = \\frac{1}{\\omega_{min}}
            ; \\tau_k = 10^{\\log{(\\tau_{min}) + \\frac{k-1}{M-1}\\log{{(
                \\frac{\\tau_{max}}{\\tau_{min}}}})}}
        and are not fit during the test (only :math:`R_{Ohm}` and :math:`R_{k}`
        are free parameters).
        In order to prevent under- or over-fitting, Schönleber et al. propose using
        the ratio of positive resistor mass to negative resistor mass as a metric
        for finding the optimal number of RC elements.
        .. math::
            \\mu = 1 - \\frac{\\sum_{R_k \\ge 0} |R_k|}{\\sum_{R_k < 0} |R_k|}
        The argument :code:`c` defines the cutoff value for :math:`\\mu`. The
        algorithm starts at :code:`M = 3` and iterates up to :code:`max_M` until a
        :math:`\\mu < c` is reached. The default of 0.85 is simply a heuristic
        value based off of the experience of Schönleber et al.
        If the argument :code:`c` is :code:`None`, then the automatic determination
        of RC elements is turned off and the solution is calculated for
        :code:`max_M` RC elements. This manual mode should be used with caution as
        under- and over-fitting should be avoided.
        [1] Schönleber, M. et al. A Method for Improving the Robustness of
        linear Kramers-Kronig Validity Tests. Electrochimica Acta 131, 20–27 (2014)
        `doi: 10.1016/j.electacta.2014.01.034
        <https://doi.org/10.1016/j.electacta.2014.01.034>`_.
        """

        def calc_mu(Rs):
            """Calculates mu for use in LinKK"""

            neg_sum = sum(abs(x) for x in Rs if x < 0)
            pos_sum = sum(abs(x) for x in Rs if x >= 0)
            return 1 - neg_sum / pos_sum

        def fitLinKK(f, ts, M, Z):
            """Fits the linKK model using scipy.optimize.least_squares"""
            initial_guess = np.append(
                min(np.real(Z)),
                np.ones(shape=(M,)) * ((max(np.real(Z)) - min(np.real(Z))) / M),
            )

            result = least_squares(
                residuals_linKK,
                initial_guess,
                method="lm",
                args=(ts, Z, f, "both", type_res),
                ftol=1e-13,
                gtol=1e-10,
            )
            #            print(result)
            p_values = result["x"]
            mu = calc_mu(p_values[1:])

            return p_values, mu

        def circuit_element_s(series):
            """sums elements in series
            Notes
            ---------
            .. math::
                Z = Z_1 + Z_2 + ... + Z_n
            """
            z = len(series[0]) * [0 + 0 * 1j]
            for elem in series:
                z += elem
            return z

        @element_metadata(num_params=1, units=["Ohm"])
        def circuit_element_R(p, f):
            """defines a resistor
            Notes
            ---------
            .. math::
                Z = R
            """
            return np.array(len(f) * [p[0]])

        @element_metadata(num_params=2, units=["Ohm", "sec"])
        def circuit_element_K(p, f):
            """An RC element for use in lin-KK model
            Notes
            -----
            .. math::
                Z = \\frac{R}{1 + j \\omega \\tau_k}
            """
            omega = np.array(f)
            return p[0] / (1 + 1j * omega * p[1])

        def eval_linKK(Rs, ts, f):
            """Builds a circuit of RC elements to be used in LinKK"""
            s, R, K = circuit_element_s, circuit_element_R, circuit_element_K
            circuit_string = "s([R({},{}),".format([Rs[0]], f.tolist())

            for i, (Rk, tk) in enumerate(zip(Rs[1:], ts)):
                circuit_string += "K({},{}),".format([Rk, tk], f.tolist())

            circuit_string = circuit_string.strip(",")
            circuit_string += "])"

            return eval(circuit_string)

        def residuals_linKK(Rs, ts, Z, f, residuals="both", type_residual="Z"):
            """Calculates the residual between the data and a LinKK fit"""

            err = Z - eval_linKK(Rs, ts, f)
            Y = Z ** -1
            errY = Y - 1 / eval_linKK(Rs, ts, f)

            if "Z" in type_residual:
                err_take = 0
            elif "Y" in type_residual:
                err_take = 1

            if residuals == "real":
                err_real = [err.real / (np.abs(Z)), errY.real / (np.abs(Y))]
                return err_real[err_take]
            elif residuals == "imag":
                err_imag = err.imag / (np.abs(Z)), errY.imag / (np.abs(Y))
                return err_imag[err_take]
            elif residuals == "both":
                z1d = np.zeros(Z.size * 2, dtype=np.float64)
                z1d[0 : z1d.size : 2] = err.real / (np.abs(Z))
                z1d[1 : z1d.size : 2] = err.imag / (np.abs(Z))

                y1d = np.zeros(Y.size * 2, dtype=np.float64)
                y1d[0 : y1d.size : 2] = errY.real / (np.abs(Y))
                y1d[1 : y1d.size : 2] = errY.imag / (np.abs(Y))
                err_both = [z1d, y1d]
                return err_both[err_take]

        def get_tc_distribution(f, M):
            """Returns the distribution of time constants for the linKK method"""

            t_max = 1 / np.min(f)
            t_min = 1 / np.max(f)

            ts = np.zeros(shape=(M,))
            ts[0] = t_min
            ts[-1] = t_max
            if M > 1:
                for k in range(2, M):
                    ts[k - 1] = 10 ** (
                        np.log10(t_min) + ((k - 1) / (M - 1)) * np.log10(t_max / t_min)
                    )

            ts *= 2 * np.pi
            return ts

        if c is not None:
            M = 0
            mu = 1
            while mu > c and M <= max_M:
                M += 1
                ts = get_tc_distribution(f, M)
                p_values, mu = fitLinKK(f, ts, M, Z)

                if M % 10 == 0:
                    pass  # print(M, mu, rmse(eval_linKK(p_values, ts, f), Z))
        else:
            M = max_M
            ts = get_tc_distribution(f, M)
            p_values, mu = fitLinKK(f, ts, M, Z)

        return (
            M,
            mu,
            eval_linKK(p_values, ts, f),
            residuals_linKK(p_values, ts, Z, f, residuals="real"),
            residuals_linKK(p_values, ts, Z, f, residuals="imag"),
        )


#%%
# def EEC_ORRpOx(ang,Rs,Cdlp,nDL,Rct,Aw,Rad,Qad,nAd,Rorr):
#    '''Other model with Faradaic '''
#    return Rs+ (1/(1j*ang*Cdlp+ 1/( Rct+ 1/( (Qad*1j*ang)*Rorr**-1))))
# Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Aw*(1j*ang)**-0.5+ Rad+ 1/(Qad**-1*(1j*ang)**-nAd+Rorr**-1))))))
# prefit,TrimData,FitOnlyTrimmedData,FreqLim, exportEIS = True, True, True, 9500, 'Text,Plot'
# def EIS_set_index_columns():
#    return [EvRHE,'RPM','Gas','SampleID','PAR_file','pH','Electrolyte','PAR_date','Model_EEC','Model_index']


def params_extra_setting(params, EISgr_data_EV):
    sID, gas, pH = (
        EISgr_data_EV["SampleID"].unique()[0],
        EISgr_data_EV.Gas.unique()[0],
        EISgr_data_EV.pH.unique()[0],
    )
    if any([i for i in sID if i in ["DW29", "DW21"]]):
        print("DW29")
        params["Rct"].max = 20e3
        params["Rct"].value = 800

    if pH > 7 and gas == "O2":
        params["Rct"].value = 20
        params["Rct"].max = 1e5
    elif pH < 7 and gas == "O2":
        params["Rct"].value = 100
        params["Rct"].max = 1e6

    if "O2" in EISgr_data_EV["Gas"].unique()[0]:
        params["Rorr"].value = 1000
        params["Rorr"].max = 1e8
    #        params['Rct'].value = 500
    #        params['Rct'].max = 1E4
    return params


# perform_prefit, TrimData, FitOnlyTrimmedData,FreqLim,exportEIS, export_raw_data  = True, False, False, 50000, 'Text,Plot',False
def EIS_fit_EEC(
    EISgr_data_EV,
    gr_EIS_ovv,
    InitP,
    EIS_dest_dir,
    perform_prefit=True,
    TrimData=False,
    FitOnlyTrimmedData=False,
    FreqLim=30000,
    exportEIS="Text,Plot",
    export_raw_data=False,
):
    #%%
    # if True:
    indexes_per_EV = []
    global EvRHE
    # FIXME preparing meta dicts
    EISgr_data_meta = (
        EISgr_data_EV[
            [i for i in EISgr_data_EV.columns if EISgr_data_EV[i].nunique() == 1]
        ]
        .iloc[0]
        .to_dict()
    )
    gr_EIS_ovv_meta = gr_EIS_ovv.iloc[0].to_dict()
    EISgr_meta_combined = {**gr_EIS_ovv_meta, **EISgr_data_meta}
    overlap_keys = [i for i in EISgr_data_meta.keys() if i in gr_EIS_ovv_meta.keys()]
    matching_keys = [
        i for i in overlap_keys if EISgr_data_meta[i] == gr_EIS_ovv_meta[i]
    ]
    unmatching_keys = [i for i in overlap_keys if i not in matching_keys]

    # Important variables defined for shortened
    Meta = namedtuple("meta", "PAR_file E_dc_RHE E_dc_RHE_mV RPM_DAC_file")
    #    PAR_file = E_dc_RHE E_dc_RHE_mV = E_dc_RHE * 1000 RPM_DAC_file =
    meta_info = Meta(
        Path(EISgr_meta_combined["PAR_file"]),
        float(EISgr_data_meta[EvRHE]),
        float(EISgr_data_meta[EvRHE]) * 1e3,
        float(EISgr_data_meta["RPM_DAC"]),
    )

    if len(unmatching_keys) == 0:
        pass
    elif len(unmatching_keys) == 1 and unmatching_keys[0] == "PAR_file":
        pass
        if (len(EISgr_data_meta) + len(gr_EIS_ovv_meta)) - (
            len(EISgr_meta_combined) + len(matching_keys) + len(unmatching_keys)
        ) != 0:
            logger.error(
                "EIS missing keys in EISgr_data_meta and ovv_meta dicts, {:s} at {:G}".format(
                    PAR_file.stem, meta_info.E_dc_RHE
                )
            )
    else:
        logger.error(
            "EIS non-matching keys in EISgr_data_meta and ovv_meta dicts, {:s} at {:G}".format(
                PAR_file.stem, meta_info.E_dc_RHE
            )
        )
    # PAR_file = Path(gr_EIS_ovv['PAR_file'].unique()[0]) # FALSE GROUP, contains different PAR files !!!!

    combined_meta_slice = [
        "postAST",
        "Electrolyte",
        "pH",
        "Gas",
        "RPM",
        "Loading_cm2",
        "Loading_name",
        EvRHE,
    ] + [
        "SampleID",
        "Measured_OCP",
        "Instrument",
        "Instrument_SN",
        "Electrode",
        "Scanrate",
        "RPM_DAC",
        "Type_action",
    ]

    EISgr_meta_add = {}
    if all([(i in EISgr_meta_combined) for i in combined_meta_slice]):
        # EISgr_meta_combined =  EISgr_data_meta[['pH','RPM','Electrolyte','SampleID','Gas','Measured_OCP','Instrument']].to_dict()
        EISgr_meta_add.update({k: EISgr_meta_combined[k] for k in combined_meta_slice})
    else:
        [i for i in combined_meta_slice if i not in EISgr_meta_combined]
        EISgr_meta_add.update(
            {
                k: EISgr_meta_combined[k]
                for k in combined_meta_slice
                if k in EISgr_meta_combined
            }
        )
        logger.error(
            "EIS outP missing columns in EISgr_data_meta, {:s} at {:G}  mV".format(
                meta_info.PAR_file.stem, meta_info.E_dc_RHE_mV
            )
        )

    EIS_outPath = EIS_dest_dir.joinpath(
        Path(
            str(Path(EISgr_data_meta["PAR_file"]).stem)
            + "_{:.0f}mV_{:.0f}rpm".format(
                meta_info.E_dc_RHE_mV, meta_info.RPM_DAC_file
            )
        )
    )
    EIS_outPath_OLD = EIS_dest_dir.joinpath(
        Path(
            str(Path(EISgr_data_meta["PAR_file"]).stem)
            + "_{:.0f}mV".format(meta_info.E_dc_RHE_mV)
        )
    )
    EIS_outPath.parent.mkdir(parents=True, exist_ok=True)
    EIS_outpath_linKK_target = FolderOps.FileOperations.CompareHashDFexport(
        pd.DataFrame(), Path(str(EIS_outPath) + "_linkK")
    ).with_suffix(".png")
    #    EISgr = EISgr.loc[EISgr["Frequency(Hz)"] < FreqLim,:]
    freq = EISgr_data_EV["Frequency(Hz)"].values
    ang = freq * 2 * np.pi
    Zre, Zim = EISgr_data_EV["Z Real"].values, EISgr_data_EV["Z Imag"].values
    Zdata = Zre + 1j * Zim
    Ydata = Zdata ** -1
    Yre, Yim = Ydata.real, Ydata.imag

    DataWeights_modulus = 1 / (Zre ** 2 + Zim ** 2)
    #    DataWeights_modulus_Y = (Yre**2+Yim**2)
    DataWeights_unitary = Zre / Zre

    EIS_data = pd.DataFrame(
        {
            "Frequency(Hz)": freq,
            "Angular": ang,
            "DATA_Z": Zdata,
            "DATA_Zre": Zre,
            "DATA_Zim": Zim,
            "DATA_-Zim": -1 * Zim,
            "DATA_Zphase": [phase(i) for i in Zdata],
            "DATA_Zmod": abs(Zdata),
            "DATA_Zangle": np.angle(Zdata, deg=True),
            "DATA_-Zangle": -1 * np.angle(Zdata, deg=True),
            "DATA_Y": Ydata,
            "DATA_Yre": Yre,
            "DATA_Yim": Yim,
            "Valid": True,
            EvRHE: meta_info.E_dc_RHE,
            "DATA_weightsmod_Z": DataWeights_modulus,
        },
        index=EISgr_data_EV.index,
    )
    EIS_data.loc[EIS_data["Frequency(Hz)"] >= FreqLim, "Valid"] = False
    EIS_data_freqlim = EIS_data.query("Valid == True")

    """ include here KK-test!! """

    restype_set = "Z"
    M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(
        EIS_data_freqlim["Frequency(Hz)"].values,
        EIS_data_freqlim["DATA_Z"].values,
        type_res=restype_set,
    )
    #    plot_linKK(EIS_data['Frequency(Hz)'].values,EIS_data.DATA_Z.values,Z_linKK,res_real,res_imag,meta=EISgr_data_meta,linKK=[M,mu],type_res=restype_set)
    # plot 1st KK validation for testing purposes
    prefit_suffix = "_prefit"
    EIS_data_linKK = KramersKronigValidation.assigncols(
        EIS_data_freqlim, Z_linKK, res_real, res_imag, suffix=prefit_suffix
    )
    #    linKK = pd.DataFrame({'linKK_Z' : Z_linKK, 'linKK_Zreal' : Z_linKK.real,'linKK_Zimag' : Z_linKK.imag,
    #                          'linKK_Y' : Z_linKK**-1, 'linKK_Yreal' : (Z_linKK**-1).real,'linKK_Yimag' : (Z_linKK**-1).imag,
    #                          'linKK_resRe' : res_real, 'linKK_resIm' : res_imag }, index=EIS_data_freqlim.index)
    #    EIS_data_linKK = pd.concat([EIS_data_freqlim,linKK],axis=1)
    #    linKK_limit = 0.0475
    linKK_limit_Re = (
        np.abs(EIS_data_linKK["linKK_resRe" + prefit_suffix]).mean()
        + 0.95 * np.abs(EIS_data_linKK["linKK_resRe" + prefit_suffix]).std()
    )
    linKK_limit_Im = (
        np.abs(EIS_data_linKK["linKK_resIm" + prefit_suffix]).mean()
        + 0.95 * np.abs(EIS_data_linKK["linKK_resIm" + prefit_suffix]).std()
    )
    KKvalid = EIS_data_linKK.loc[
        (
            (
                (
                    (
                        np.abs(EIS_data_linKK["linKK_resRe" + prefit_suffix])
                        < linKK_limit_Re
                    )
                    & (
                        np.abs(EIS_data_linKK["linKK_resIm" + prefit_suffix])
                        < linKK_limit_Im
                    )
                )
                & (EIS_data_linKK["Frequency(Hz)"] > 2)
            )
            | (EIS_data_linKK["Frequency(Hz)"] <= 2)
        )
        & (EIS_data_linKK["Frequency(Hz)"] < FreqLim)
    ]
    EIS_data.loc[~EIS_data.index.isin(KKvalid.index), "Valid"] = False
    M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(
        KKvalid["Frequency(Hz)"].values,
        KKvalid["DATA_Z"].values,
        type_res=restype_set,
        max_M=50,
    )

    EIS_data_valid = KramersKronigValidation.assigncols(
        KKvalid, Z_linKK, res_real, res_imag
    )
    #    EIS_data_valid_KK = KramersKronigValidation.assigncols(EIS_data,Z_linKK,res_real,res_imag)

    plot_linKK(
        EIS_data_valid,
        EIS_data,
        res_scale=1,
        meta=EISgr_data_meta,
        linKK=[M, mu],
        type_res=restype_set,
        save_target=EIS_outpath_linKK_target,
        plot_prefit=True,
        KK_valid_limit=(linKK_limit_Re, linKK_limit_Im),
    )

    #    outData.loc[outData['Frequency(Hz)'].isin([i for i in outData['Frequency(Hz)'].values if i not in KKvalid['Frequency(Hz)'].values])]
    EIS_data.loc[~EIS_data.index.isin(EIS_data_valid.index), "Valid"] = False
    if export_raw_data:
        EIS_outpath_linKK_target = FolderOps.FileOperations.CompareHashDFexport(
            EIS_data, Path(str(EIS_outPath) + "_raw-data")
        ).with_suffix(".xlsx")
        index_raw_output = {
            "PAR_file": meta_info.PAR_file,
            "Type_output": "EIS_raw_data",
            "Type_exp": "EIS",
            "DestFile": EIS_outpath_linKK_target,
            "E_V": meta_info.E_dc_RHE,
        }
        indexes_per_EV.append(index_raw_output)
    #        index_raw = {'PAR_file'  : nm2,'DestFile' : Parsout_path_target,'Type_output' : 'EIS_Pars', 'Type_exp' : 'EIS','nModels' : Pars.Model_EEC.nunique(), 'Models_used' : Pars.Model_EEC.unique()}
    linKK_pars = {
        "Unvalid_freqs": ", ".join(
            EIS_data.query("Valid == False")["Frequency(Hz)"].astype(str).to_list()
        ),
        "linKK_M": M,
        "linKK_mu": mu,
    }
    #    M, mu, Z_linKK, res_real, res_imag = KramersKronigValidation.linKK(KKvalid['Frequency(Hz)'].values,KKvalid['DATA_Z'].values,type_res=restype_set,max_M=2)
    #    outData = pd.DataFrame(data=EIS_data,index=EISgr.index)
    #    outData = pd.DataFrame(data=EIS_data,index=EISgr.index)
    #    'DATA_Zmod' : [cmath.phase(i) for i in Zdata],'DATA_Zmod' : [cmath.polar(i)[0] for i in Zdata][cmath.polar(i)[1]*180/np.pi for i in Zdata],
    #%%
    # if True:
    EIS_data_KKvalid = EIS_data.query("Valid == True")

    Z_KKv, ang_KKv = EIS_data_KKvalid.DATA_Z.values, EIS_data_KKvalid.Angular.values
    # create a set of Parameters
    Rs_guess = EIS_data_KKvalid.DATA_Zre.min()
    standard_init_params = Parameters()

    standard_init_params.add("Rs", value=Rs_guess, min=0.1, max=5e2)
    standard_init_params.add("Cdlp", value=5e-04, min=1e-08, max=1e-1)
    standard_init_params.add("nDL", value=0.90, min=0, max=1)
    standard_init_params.add(
        "Rct", value=3e2, min=1e-3, max=1e6, vary=True, brute_step=1000
    )
    #    params.add('Aw', value= 3000, min=5.0, max=1E6)
    #    params.add('Rad', value= 5E5, min=0.01, max=1E7)
    standard_init_params.add("Qad", value=5e-03, min=1e-07, max=1e-01)
    standard_init_params.add("nAd", value=0.9, min=0, max=1)
    standard_init_params.add("Rorr", value=3e3, min=0.01, max=1e10, brute_step=1000)
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    standard_init_params.add_many(
        ("R3", 100, True, 0.01, 1e8, None, None),
        ("Q3", 1e-04, True, 1e-09, 1, None, None),
        ("n3", 0.8, True, 0, 1, None, None),
    )

    standard_init_params = params_extra_setting(standard_init_params, EISgr_data_EV)
    # Old style:  mod,mod2 = Model(EEC_Singh2015_RQRQR,name='Singh2015_RQRQR'),Model(EEC_ORRpOx,name='Bandarenka_2011_RQRQR')
    # FIXME  skip(5,Model(EEC_Bandarenka_Ads,name='Bandarenka2011_RQRW'))
    mod_lst = EEC_models_index()
    # TODO check whats wrong with Bandarenka Fit model
    statiscial_weights = EIS_data_KKvalid.DATA_weightsmod_Z.values
    MTHDS = [
        "leastsq",
        "nelder",
        "lbfgsb",
        "anneal",
        "powell",
        "cobyla",
        "slsqp",
        "differential-evolution",
        "ampgo",
        "basinhopping",
        "dual_annealing",
        "brute",
    ]
    pars_lst, fit_data_valid_lst, chisqr_start, previous_params = [], [], 1, []
    for mod_indx, model_set in mod_lst[::]:
        modname = model_set.name
        params_model = model_set.make_params()
        lstsq_method = MTHDS[0]

        #        print(params_model['nAd'])
        if InitP != {}:
            pass
        #            print('use previous')
        #            params = InitP1.query('Model_EEC')
        else:
            for pn in params_model:
                if previous_params and pn in previous_params:
                    params_model.add(previous_params[pn])
                else:
                    params_model.add(standard_init_params[pn])

        if "nAd" in params_model:
            if modname in ["Model(Singh2015_RQRWR)", "Model(Bandarenka2011_RQRW)"]:
                params_model["nAd"].set(value=0.5, vary=False)
            else:
                params_model["nAd"].set(value=0.9, vary=True)

        # print('After', params_model['Rct'])
        # params_model.pretty_print()
        #    Weights = np.array([i/i for i in range(1,len(ang)+1)][::-1])
        #    result = minimize(model_ORR, params, args=(ang,Zdata),method=methods[0])
        #    result2 = minimize(model_ORRpOx, params, args=(ang,Zdata),method=methods[0])
        #    z == z.real + z.imag*1j
        #    if Initp1 is not {} or Initp2 and not {}:
        #        params = InitP1
        if perform_prefit == True:
            #        errDF = pd.DataFrame()
            #        pre_init1,pre_init2 =  mod.eval(pre_prefit1.params, ang=ang),mod2.eval(pre_prefit2.params, ang=ang)
            #        pre_prefit1 , pre_prefit2
            prefit_method = MTHDS[-5]
            #            print('Prefit on model "{:s}" with method {:s}' .format(modname, prefit_method))
            #            params_model.pretty_print()
            pre_prefit = model_set.fit(
                Z_KKv,
                params_model,
                ang=ang_KKv,
                weights=statiscial_weights,
                method=prefit_method,
            )
            best_trial = pre_prefit

            if "Rct" in params_model:
                for Rct_guess in [5, 50, 150, 500, 1e3, 3e3, 8e3, 15e3]:
                    best_trial.params["Rct"].set(value=Rct_guess, vary=True)
                    if "Rorr" in params_model:
                        for Rorr_guess in [100, 500, 1e3, 3e3, 1e4, 3e4, 1e5]:
                            best_trial.params["Rorr"].set(value=Rorr_guess)
                            trial_prefit_brute = model_set.fit(
                                Z_KKv,
                                best_trial.params,
                                ang=ang_KKv,
                                weights=statiscial_weights,
                                method=prefit_method,
                            )
                            if trial_prefit_brute.redchi < best_trial.redchi:
                                best_trial = trial_prefit_brute
                    #                                print('best trial updated: R_ct[{:.3G}] -> {:.3G}  & R_orr[{:.3G}] -> {:.3G}, improvement of {:.2%} ({:.2E})'.format(Rct_guess, trial_prefit_brute.params['Rct'].value,
                    #                                      Rorr_guess, trial_prefit_brute.params['Rorr'].value,  1E-2*(best_trial.chisqr-trial_prefit_brute.chisqr) /best_trial.chisqr, trial_prefit_brute.redchi))

                    #                                trial_prefit_brute.params.pretty_print()
                    else:
                        trial_prefit_brute = model_set.fit(
                            Z_KKv,
                            best_trial.params,
                            ang=ang_KKv,
                            weights=statiscial_weights,
                            method=prefit_method,
                        )
                        if trial_prefit_brute.redchi < best_trial.redchi:
                            #                            print('best trial updated: Rct[{:.3G}] -> {:.3G} improvement of {:.2%} ({:.2E})'.format(Rct_guess,trial_prefit_brute.params['Rct'].value, 1E-2*(best_trial.chisqr-trial_prefit_brute.chisqr) /best_trial.chisqr, trial_prefit_brute.redchi))
                            best_trial = trial_prefit_brute
            #                            trial_prefit_brute.params.pretty_print()
            #            print(pre_prefit.fit_report())
            prefit = model_set.fit(
                Z_KKv,
                best_trial.params,
                ang=ang_KKv,
                weights=statiscial_weights,
                method="leastsq",
            )
            init = model_set.eval(prefit.params, ang=ang_KKv)
            logger.info(
                "Prefit starting with method {:s} for {:s} at {:0.2G} on model {:s},ChiSqr = {:.3G}".format(
                    MTHDS[-5],
                    meta_info.PAR_file.stem,
                    meta_info.E_dc_RHE,
                    model_set.name,
                    prefit.chisqr,
                )
            )
            InitParam = prefit.params
            #        WeightsRes = ((DataRes2['index']*DataRes2[0]).values[1::2]**2 + (DataRes2['index']*DataRes2[0]).values[0::2]**2)
            #        WeightsRes =   abs(errRE_init+1j*errIM_init)
            #        Weights = WeightsRes
            #        DataRes2[(DataRes2 < Rmaxlim) & (DataRes2 > Rminlim)]
            #        DateResFiltered = DataRes2.loc[(DataRes2[0] < Rmaxlim) & (DataRes2[0] > Rminlim)]
            #            InitParam = prefit.params
            if TrimData == True:
                errRE_init, errIM_init = (init.real - Z_KKv.real) / abs(Z_KKv), (
                    init.imag - Z_KKv.imag
                ) / abs(Z_KKv)
                errRE_initSTD, errIM_initSTD = errRE_init.std(), errIM_init.std()
                errRE_initMEAN, errIM_initMEAN = errRE_init.mean(), errIM_init.mean()
                trimSTDmax = 2
                errRmaxlim, errRminlim = (
                    errRE_initMEAN + trimSTDmax * np.abs(errRE_initSTD)
                ), (errRE_initMEAN - trimSTDmax * np.abs(errRE_initSTD))
                errRmaxlimIm, errRminlimIm = (
                    errIM_initMEAN + trimSTDmax * np.abs(errIM_initSTD)
                ), (errIM_initMEAN - trimSTDmax * np.abs(errIM_initSTD))
                #        pd.DataFrame([freq,prefit2.residual])
                #        abs(errRE+1j*errIM)
                errDF = pd.DataFrame([errRE_init, errIM_init, freq]).T
                #        errDF.loc[~((errDF[0] < errRmaxlim) & (errDF[0] > errRminlim)),0] = errRE_initMEAN
                #        errDF.loc[~((errDF[1] < errRmaxlimIm) & (errDF[1] > errRminlimIm)),1] = errIM_initMEAN
                #        badFreq = errDF.loc[~((errDF[1] < errRmaxlimIm) & (errDF[1] > errRminlimIm)),2]
                badFreq = errDF.loc[
                    ~(
                        ((errDF[0] < errRmaxlim) & (errDF[0] > errRminlim))
                        & ((errDF[1] < errRmaxlimIm) & (errDF[1] > errRminlimIm))
                    ),
                    2,
                ]
                #        | ((errDF[0] < errRmaxlim) & (errDF[0] > errRminlim))),2]
                TrimmedOutData = EIS_data_KKvalid.loc[
                    ~EIS_data_KKvalid["Frequency(Hz)"].isin(badFreq.values)
                ]
                trimAng = TrimmedOutData["Frequency(Hz)"].values * 2 * np.pi
                TrimmedWeightsResOutlier = abs(
                    errDF.loc[~errDF[2].isin(badFreq.values), 0]
                    + errDF.loc[~errDF[2].isin(badFreq.values), 1] * 1j
                )
                WeightsResOutlier = abs(errDF[0] + errDF[1] * 1j)
                prefit_res_mean, prefit_res_std = (
                    prefit.residual.mean(),
                    prefit.residual.std(),
                )
                Rmaxlim, Rminlim = (
                    prefit_res_mean + trimSTDmax * np.abs(prefit_res_std)
                ), (prefit_res_mean - trimSTDmax * np.abs(prefit_res_std))
                prefit.residual[
                    (prefit.residual < Rmaxlim) & (prefit.residual > Rminlim)
                ]

                DataRes = pd.DataFrame(
                    prefit.residual, Zdata.view(np.float)
                ).reset_index()
                DataRes.loc[
                    ~((DataRes[0] < Rmaxlim) & (DataRes[0] > Rminlim)), 0
                ] = 1e-12

                trim_ang = TrimmedOutData["Frequency(Hz)"].values * 2 * np.pi
                #            Zre, Zim = EISgr_data_EV['Z Real'].values, EISgr_data_EV['Z Imag'].values
                Z_TrimmedData = TrimmedOutData["DATA_Z"].values

                #            Ydata = 1/Z_TrimmedData
                #            WeightsRes = TrimmedWeightsResOutlier
                EISgr_trimmed_out = EISgr_data_EV.loc[
                    ~EISgr_data_EV["Frequency(Hz)"].isin(
                        TrimmedOutData["Frequency(Hz)"].values
                    )
                ]
                EISgr_data_EV = EISgr_data_EV.loc[
                    EISgr_data_EV["Frequency(Hz)"].isin(
                        TrimmedOutData["Frequency(Hz)"].values
                    )
                ]

                prefit = model_set.fit(
                    Z_TrimmedData,
                    pre_prefit.params,
                    ang=trim_ang,
                    weights=1 / abs(Z_TrimmedData),
                    method=prefit_method,
                )
                #                prefit2 = mod2.fit(Z_TrimmedData,pre_prefit2.params,ang=trim_ang,weights= 1/abs(Z_TrimmedData), method='differential-evolution')
                init = model_set.eval(prefit.params, ang=trim_ang)
                InitParam = prefit.params
                EIS_Trimming_plot(
                    EISgr_data_EV,
                    gr_EIS_ovv,
                    EIS_dest_dir,
                    outData,
                    TrimmedOutData,
                    meta_info.E_dc_RHE,
                )
            #            WeightsRes =   abs(errRE_init+1j*errIM_init)
            #            1/(Zre**2+1/Zim**2)
            if FitOnlyTrimmedData == True:
                ang, Zdata = trim_ang, Z_TrimmedData
                outData = TrimmedOutData
        else:
            logger.warning(
                "No prefit for {:s} at {:0.2f} with model {:s}".format(
                    PAR_file.stem, meta_info.E_dc_RHE, modname
                )
            )
            init = model_set.eval(params_model, ang=ang_KKv)
            InitParam = params_model

        #        out1 = mod.fit(Z_TrimmedData,InitParam1,ang=trim_ang,weights= 1/Z_TrimmedData, method=MTHDS[0])
        #        out2 = mod2.fit(Z_TrimmedData,InitParam2,ang=trim_ang,weights=1/Z_TrimmedData, method=MTHDS[0])
        ##       out2_initW = mod2.fit(Zdata,InitParam2,ang=ang,weights=DataWeights, method=MTHDS[-2])
        #        fit1,fit2 = out1.eval(ang=trim_ang),out2.eval(ang=trim_ang)
        ##=== EIS REFIT USING RESIDUAL ===#
        #        errRE,errIM = (fit2.real-Z_TrimmedData.real)/abs(Z_TrimmedData), (fit2.imag-Z_TrimmedData.imag)/abs(Z_TrimmedData)
        #        errRE1,errIM1 = (fit1.real-Z_TrimmedData.real)/abs(Z_TrimmedData), (fit1.imag-Z_TrimmedData.imag)/abs(Z_TrimmedData)
        #        outData =
        errRE_init, errIM_init = (init.real - Z_KKv.real) / abs(Z_KKv), (
            init.imag - Z_KKv.imag
        ) / abs(Z_KKv)
        #        Weights = 1/abs(Z_KKv)
        if modname in ["Model(Singh2015_RQRWR)", "Model(Bandarenka2011_RQRW)"]:
            InitParam["nAd"].set(value=0.5, vary=False)
        out = model_set.fit(
            Z_KKv,
            InitParam,
            ang=ang_KKv,
            weights=statiscial_weights,
            method=lstsq_method,
        )
        logger.info(
            "Fit out with method {:s} for {:s} at {:0.2G} on model {:s},ChiSqr = {:.3G}, RedChiSqr = {:.3G}".format(
                lstsq_method,
                meta_info.PAR_file.stem,
                meta_info.E_dc_RHE,
                model_set.name,
                out.chisqr,
                out.redchi,
            )
        )
        ### == RETRY SECTION: when first fit is not good enough.. some initial values are re-set and fit is run again ... ####
        retry_attempt = "no"
        #        if out.chisqr > 1E-4:
        #            large_chisqr = out.chisqr
        #            retry_method = MTHDS[0]
        #            logger.warning('Restarting Fit >>> Chisqr larger than 0.07 (={:0.3G}) for model {:s}.\n  fit with method {:s} for {:s} at {:0.2G}'.format(large_chisqr, modname, retry_method, PAR_file.stem, E_dc_RHE))
        #            retry_params = out.params
        #            if 'O2' in gr_EIS_ovv['Gas'].unique()[0]:
        #                if 'Rct' in retry_params.valuesdict().keys():
        #                    retry_params['Rct'].value = 20
        #                if 'nAd' in retry_params.valuesdict().keys():
        #                    retry_params['nAd'].value = 0.77
        #                if 'Rorr' in retry_params.valuesdict().keys():
        #                    retry_params['Rorr'].value = 1000
        #            if 'N2' in gr_EIS_ovv['Gas'].unique()[0] and EISgr_data_EV.pH.unique()[0] < 7:
        #                if 'Rct' in retry_params.valuesdict().keys():
        #                    retry_params['Rct'].value = 1
        #                if 'nAd' in retry_params.valuesdict().keys():
        #                    retry_params['nAd'].value = 1
        #                if 'Rorr' in retry_params.valuesdict().keys():
        #                    retry_params['Rorr'].value = 5E5
        #            if 'nAd' in retry_params:
        #                if modname in ['Model(Singh2015_RQRWR)','Model(Bandarenka2011_RQRW)']:
        #                    retry_params['nAd'].set(value= 0.5, vary=False)
        #                else:
        #                    retry_params['nAd'].set(value= 0.8, vary=True)
        #
        #            retry_out = model_set.fit(Z_KKv,retry_params,ang=ang_KKv,weights= statiscial_weights , method = retry_method)
        #            retry_chisqr_diff = np.abs(large_chisqr - retry_out.chisqr)
        #            retry_attempt = 'tried'
        #            if retry_out.chisqr < large_chisqr:
        #                out = retry_out
        #                retry_attempt = 'tried and result'
        #                logger.warning('Results of Re-fit Chisqr {:0.3f}.\n used method {:s} for {:s} at {:G}'.format(retry_chisqr_diff,MTHDS[0],PAR_file.stem,E_dc_RHE))
        #            else:
        #                pass
        #        else:
        #            pass
        #        if out.chisqr < chisqr_start:
        #            previous_params = out.params
        #            chisqr_start = out.chisqr
        fit = out.eval(ang=ang_KKv)
        #    out2_initW = mod2.fit(Zdata,InitParam2,ang=ang,weights=DataWeights, method=MTHDS[-2])
        # === EIS REFIT USING RESIDUAL ===#
        errRE, errIM = (fit.real - Z_KKv.real) / abs(Z_KKv), (
            fit.imag - Z_KKv.imag
        ) / abs(Z_KKv)
        #    print('Prefit Weights: DataWeights)
        # 'Init_Z1' : init1,'Init_Z2' : init2, 'FIT1_Z' : fit1,
        EIS_data_KKvalid_fit = EIS_data_KKvalid.assign(
            **{
                "INIT_Yim": (1 / init).imag,
                "INIT_Yre": (1 / init).real,
                "INIT_Zphase": [phase(i) for i in init],
                "INIT_errRe": errRE_init,
                "INIT_errIm": errIM_init,
                "FIT_Zre": fit.real,
                "FIT_Zim": fit.imag,
                "FIT_-Zim": -1 * fit.imag,
                "FIT_Yre": (1 / fit).real,
                "FIT_Yim": (1 / fit).imag,
                "FIT_Zmod": abs(fit),
                "FIT_Zphase": [phase(i) for i in fit],
                "FIT_Zangle": np.angle(fit, deg=True),
                "FIT_-Zangle": -1 * np.angle(fit, deg=True),
                "errRe": errRE,
                "errIm": errIM,
                "Model_EEC": modname,
                "Model_index": mod_indx,
                EvRHE: meta_info.E_dc_RHE,
            }
        )
        # add the metadata to the EIS_data_KKvalid_fit DF !!
        EISgr_meta_add_to_fit = pd.DataFrame(
            [EISgr_meta_add] * len(EIS_data_KKvalid_fit),
            index=EIS_data_KKvalid_fit.index,
        )
        EIS_data_KKvalid_fit = pd.concat(
            [EIS_data_KKvalid_fit, EISgr_meta_add_to_fit], axis=1
        )

        # EIS_data_KKvalid_fit = EIS_data_KKvalid_fit.assign(**[{i[0] : [i[1]]*len(EIS_data_KKvalid_fit)} for i in EISgr_meta_add.items()])
        # EIS_data_KKvalid_fit = EIS_data_KKvalid_fit.set_index(EIS_set_index_columns())

        #    metadata to OutData
        #'FIT2_Zmod' : [cmath.polar(i)[0] for i in fit2], FIT2_Zphase'[cmath.polar(i)[1]*180/np.pi for i in fit2]
        #    print('%s\n'%out.message,out.fit_report(min_correl=0.50))
        #    print('%s\n'%out2.message,out2.fit_report(min_correl=0.90)) EvRHE : [EISgr[EvRHE].unique()[0]]*len(outData)
        # === Output organizing ===
        outP = out.best_values
        out_params_stderss = [
            (i + "_stderr", out.params.__getitem__(i).stderr) for i in out.params
        ]
        out_params_correl = [
            (i + "_correl", out.params.__getitem__(i).correl) for i in out.params
        ]
        outP.update(
            dict(
                zip(
                    [i[0] for i in out_params_stderss],
                    [i[1] for i in out_params_stderss],
                )
            )
        )

        outP.update(EISgr_meta_add)  # from beginning of function
        outP.update(linKK_pars)  # from beginning of function after linKK validation
        outP.update(
            {
                EvRHE: meta_info.E_dc_RHE,
                "PAR_file": meta_info.PAR_file,
                "PAR_date": EISgr_meta_combined["PAR_date"],
                "RPM_DAC_file": meta_info.RPM_DAC_file,
            }
        )
        #        outP.update(extraP)
        if "Qad" not in outP.keys():
            outP.update({"Qad": 0, "nAd": None})
        xtra = {
            "Rct_kin": outP["Rct"] ** -1,
            "Qad+Cdlp": outP["Qad"] + outP["Cdlp"],
            "Chisqr": out.chisqr,
            "RedChisqr": out.redchi,
            "Model_EEC": modname,
            "Model_index": mod_indx,
            "lmfit_method": out.method,
            "lmfit_message": out.message,
            "lmfit_out": out,
            "retry_attempt": retry_attempt,
            "lmfit_var_names": str(tuple(out.var_names)),
        }
        outP.update(xtra)
        pars_lst.append(pd.DataFrame(outP, index=[EIS_data_KKvalid_fit.index[0]]))
        fit_data_valid_lst.append(EIS_data_KKvalid_fit)
    #        lmfit_out_lst.append({'Model_EEC' : modname,'Model_index' : mod_indx,'lmfit_out' : out})
    pars_models = pd.concat(pars_lst, sort=False, ignore_index=True)
    EIS_fit_data = pd.concat(fit_data_valid_lst, sort=False, ignore_index=True)

    #%%
    #    EIS_fit_data.groupby('Model_index').plot(x=['DATA_Yre','FIT_Yre'],y=['DATA_Yim','FIT_Yim'],kind='scatter')
    # +++ PLOTTING of EIS spectrum +++

    EIS_outPath_target = FolderOps.FileOperations.CompareHashDFexport(
        EIS_fit_data, EIS_outPath.with_suffix(".xlsx")
    )
    index_data_output = {
        "PAR_file": meta_info.PAR_file,
        "Type_output": "EIS_fit_data",
        "Type_exp": "EIS",
        "DestFile": EIS_outPath_target,
        "E_V": meta_info.E_dc_RHE,
        "RPM_DAC": meta_info.RPM_DAC_file,
    }
    indexes_per_EV.append(index_data_output)

    if "Plot" in exportEIS:
        EIS_outPath_target_png = EIS_outPath_target.with_suffix(".png")
        index_dataplot_output = {
            "PAR_file": meta_info.PAR_file,
            "Type_output": "EIS_fit_data_png",
            "Type_exp": "EIS",
            "DestFile": EIS_outPath_target_png,
            "E_V": meta_info.E_dc_RHE,
            "RPM_DAC": meta_info.RPM_DAC_file,
        }
        indexes_per_EV.append(index_dataplot_output)
        EIS_plotting_per_EV(
            EIS_fit_data,
            pars_models,
            EISgr_meta_combined,
            EIS_outPath_target_png,
            plot_show=False,
        )
    #    elif 'Plot' and not 'Text' in exportEIS:
    #        EIS_plotting(EIS_fit_data, pars_models,)
    else:
        logger.warning(
            "EIS fitting only text output {:s}".format(str(EIS_outPath_target))
        )
    #%%
    if indexes_per_EV:
        indexes_per_EV_out = pd.DataFrame(indexes_per_EV)
    else:
        indexes_per_EV_out = pd.DataFrame()

    return EIS_fit_data, pars_models, indexes_per_EV_out


#          pd.DataFrame(outP1,index=[EISgr[EvRHE].unique()[0]]),pd.DataFrame(outP2,index=[EISgr[EvRHE].unique()[0]]),outData,out1,out2
def EIS_plotting_per_EV(
    EIS_fit_data,
    pars_models,
    EISgr_meta_combined,
    EIS_outPath_target_png,
    plot_show=False,
):
    #%%
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(16, 16))
    ax, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    ax5, ax6 = axes[2, 0], axes[2, 1]

    sID = EISgr_meta_combined["SampleID"]
    # EIS_fit_data.index.get_level_values('SampleID').unique()[0]
    Elec, pH = EISgr_meta_combined["Electrolyte"], EISgr_meta_combined["pH"]
    # EIS_fit_data.index.get_level_values('Electrolyte').unique()[0], EIS_fit_data.index.get_level_values('pH').unique()[0]
    EappV = EISgr_meta_combined[EvRHE]
    # EIS_fit_data.index.get_level_values(EvRHE).unique()[0]
    fig.suptitle(
        "%s, %s: \n$\mathrm{E_{DC}\/ =\/ %.2f\/ V_{RHE} \/ in \/%s\/at\/pH\/=\/%.0f}$"
        % (sID, EIS_outPath_target_png.stem, EappV, Elec, pH)
    )

    dataC, fitC, extraC, initC = "tab:blue", "tab:red", "gold", "gray"
    #        outData.plot(x='FIT_Zre',y='FIT_-Zim',kind='line',c='orange',ax=ax)
    """axes: Data plots """
    #    EIS_fit_data.plot(x='DATA_Zre',y='DATA_-Zim',kind='scatter',c=dataC,s=100,ax=ax)
    ax.scatter(
        EIS_fit_data["DATA_Zre"], EIS_fit_data["DATA_-Zim"], c=dataC, s=100, alpha=0.5
    )
    ax.set(xlim=(0, ax.axis("equal")[1]), ylim=(0, ax.axis("equal")[-1]))
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
    for modnm, modgr in EIS_fit_data.groupby("Model_EEC"):
        """ax: Impedance (Nyquist) Plots"""
        #     mod,mod2 = Model(EEC_ORR),Model(EEC_ORRpOx)
        #        EIS_fit_data.plot(x='FIT2_Zre',y='FIT2_Zim',kind='line',c=fitC,lw=2.5,ax=ax,label='FIT_2 (=EEC_ORRpOx)')
        #        modgr.plot(x='FIT_Zre',y='FIT_-Zim',kind='line',lw=2.5,ax=ax,label='FIT: {0}'.format(modnm))
        ax.plot(
            modgr["FIT_Zre"],
            modgr["FIT_-Zim"],
            lw=2.5,
            alpha=0.7,
            label="FIT: {0}".format(modnm),
        )
        ax.legend(ncol=2, fontsize=12, loc="upper left", bbox_to_anchor=(0.1, 1.3))
        #        ax.plot(EIS_fit_data['FIT2_Zre'].values,EIS_fit_data['FIT2_Zre'].values,ls=':',c='k',lw=1)
        #    outData.plot(x='FIT1_Zre',y='FIT1_Zre',kind='line',ls='.',c='k',lw=10,ax=ax)

        """ax2: Admittance Plots"""
        ax2.plot(modgr["FIT_Yre"], modgr["FIT_Yim"], lw=2.5, alpha=0.7)

        #        modgr.plot(x='FIT_Yre',y='FIT_Yim',kind='line',lw=2.5,ax=ax2,alpha=0.7)
        ax2.legend(ncol=2, fontsize=12, loc="upper right")
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
            ax=ax3,
            logy=False,
            logx=True,
            alpha=0.7,
        )
        ax3.set_xlim(0.5, 3e4)
        ax3.legend(ncol=2, fontsize=12)
        modgr.plot(
            x="Frequency(Hz)",
            y="FIT_-Zangle",
            kind="line",
            lw=2.5,
            ax=ax4,
            logy=False,
            logx=True,
            alpha=0.7,
        )
        ax4.legend(ncol=2, fontsize=12)
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
            alpha=0.6,
            label="errRe: {0}".format(modnm),
        )
        ax5.legend(ncol=1, fontsize=12)
        ax5.set_ylim(-0.15, 0.3)
        ax6.set_ylabel("$\Delta$ Re $(\%)$")
        modgr.plot(
            x="Frequency(Hz)",
            y="errIm",
            kind="line",
            ax=ax6,
            logy=False,
            logx=True,
            alpha=0.6,
        )
        ax6.legend([])
        ax6.set_ylim(-0.15, 0.3)
        ax6.set_ylabel("$\Delta$ Im $(\%)$")

    #        outData.plot(x='Frequency(Hz)',y='errRe1',kind='line',c=fitC,ax=ax5,logy=False,logx=True)
    #        outData.plot(x='Frequency(Hz)',y='errIm1',kind='line',c=extraC,ax=ax6,logy=False,logx=True)
    #    fig_path = EIS_dest_dir.joinpath(Path(str(Path(EISgr['File'].unique()[0]).stem)+'_%.0fmV.png' %(E_dc*1000)))
    #    if lmfit_out_lst:
    if len(pars_models.sort_values("RedChisqr", ascending=True).head(2)) == 2:
        out1, out2 = (
            pars_models.sort_values("RedChisqr", ascending=True)
            .head(2)["lmfit_out"]
            .to_list()
        )
        #    out1,ou2 = lmfit_out_lst[best_2mods[0]]
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
        ax4.text(
            2.25,
            1.1,
            out1.fit_report(min_correl=0.99),
            transform=ax.transAxes,
            fontsize=14,
            verticalalignment="top",
            bbox=props,
        )
        ax6.text(
            2.25,
            -0.7,
            out2.fit_report(min_correl=0.99),
            transform=ax.transAxes,
            fontsize=14,
            verticalalignment="top",
            bbox=props,
        )
    else:
        if not pars_models.sort_values("RedChisqr", ascending=True).head(1).empty():
            props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
            out1 = pars_models.sort_values("RedChisqr", ascending=True).head(1)
            ax4.text(
                2.25,
                1,
                out1.fit_report(min_correl=0.99),
                transform=ax.transAxes,
                fontsize=15,
                verticalalignment="top",
                bbox=props,
            )

    if plot_show:
        plt.show()
    plt.savefig(EIS_outPath_target_png, bbox_inches="tight", dpi=100)
    plt.close()


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
    fig.suptitle(
        "%s, %s: \n$\mathrm{E_{DC}\/ =\/ %.2f\/ V_{RHE} \/ in \/%s\/at\/pH\/=\/%.0f}$"
        % (
            EISgr_data_EV["SampleID"].unique()[0],
            Path(EISgr_data_EV["PAR_file"].unique()[0]).stem,
            E_dc,
            EISovv["Electrolyte"].unique()[0],
            EISovv["pH"].unique()[::],
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

    AllData_E_file, EISovv, Pars, E_data_combined_path_target = (
        spectra_combined,
        EVgr,
        EVgr,
        PDD_rpm.joinpath("EIS_combined_{0}.jpg".format("_".join([str(i) for i in EV]))),
    )


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
    AllData_E_file = AllData_E_file.loc[:, ~AllData_E_file.columns.duplicated()]
    Lenrows = AllData_E_file[combined_grouper].nunique()
    if Lenrows == 1 and AllData_E_file["RPM_DAC"].nunique() > 1:
        combined_grouper = "RPM_DAC"
    #    fig,axes = plt.subplots(nrows=Lenrows ,sharex=True,sharey=True)
    ht, wd = 15, 15
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
    AllData_E_file.sort_values(
        by=[combined_grouper, "Frequency(Hz)"], ascending=True, inplace=True
    )
    Lenrows = AllData_E_file[combined_grouper].nunique()
    Emin, Emax = (
        AllData_E_file[combined_grouper].min(),
        AllData_E_file[combined_grouper].max(),
    )
    E_diff_minmax = Emax - Emin if not Emax == Emin else Emin
    E_lst_bestmod = []
    for Ev, Egr in AllData_E_file.groupby(combined_grouper):

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
                    Pars.loc[Pars[combined_grouper] == Ev, "RedChisqr"],
                    Pars.loc[Pars[combined_grouper] == Ev].RedChisqr.min(),
                    rtol=5e-2,
                ),
                ["RedChisqr", "Model_EEC", "Chisqr"],
            ]
            .sort_values(by="RedChisqr")
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
            )
    #            ax.plot(Edata['FIT1_Yre'].values,Edata['FIT1_Yim'].values+len(Edata)*[maxYim*En],c=fit1C,ls='dotted',lw=3,alpha=0.7)
    #            if modcount == 0:
    #                ax.annotate("$\mathrm{%.2f} \/ V_{RHE} $"%Ev,
    #                    xy=(maxYre*1.1, maxYim*Enum), xycoords='data')
    #    ax.set_ylim(0,maxYim*(En+1))
    try:
        fast_checking_EEC_models = [
            "Model(Singh2015_RQRQR)",
            "Model(Singh2015_RQRWR)",
            "Model(Singh2015_R3RQ)",
            "Model(Bandarenka_2011_RQRQR)",
        ]
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

            sc2Rs = ax2.scatter(
                bestmod_xscatter,
                Pmodgr.Rs.values,
                label="Rs " + modnm,
                s=80,
                marker="s",
                alpha=0.7,
            )
            ax2.set_title(first_sct_title)
            sc2Rct = ax2_twin.scatter(
                bestmod_xscatter, Pmodgr.Rct.values, label="Rct", marker="*", s=100
            )

            ax3.scatter(
                bestmod_xscatter,
                Pmodgr.Rct.values,
                s=100,
                label="Rct " + modnm,
                marker="*",
            )
            ax3_twin.scatter(
                bestmod_xscatter, Pmodgr.Rorr.values, label="Rorr", marker="^", s=100
            )

            ax4.scatter(
                bestmod_xscatter,
                Pmodgr.Cdlp.values,
                s=100,
                marker="o",
                label="Cdpl " + modnm,
            )
            ax4_twin.scatter(
                bestmod_xscatter, Pmodgr.nDL.values, marker="X", s=100, label="nDL"
            )

            ax5.scatter(
                bestmod_xscatter,
                Pmodgr.Qad.values,
                s=100,
                marker="D",
                label="Qad " + modnm,
            )
            ax5_twin.scatter(
                bestmod_xscatter, Pmodgr.nAd.values, marker="P", s=100, label="nAd"
            )

        #    a = ax2.get_legend_handles_labels()
        #    b = ax2_twin.get_legend_handles_labels()
        #    a+b
        #    sc2Rs+sc2Rct
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


#%%
if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    # set log level
    logger.setLevel(logging.INFO)

#    try:
#        from FolderOps import *
#
#
#        # define file handler and set formatter
#        file_handler = logging.FileHandler(FindExpFolder('VERSASTAT').IndexDir.joinpath('PAR_EIS_fit_logfile.log'))
#
#    except Exception as e:
#        print('Expception EIS import {0}, from {1}'.format(e,__name__))
#        from .FolderOrganizing.FolderOps import *
#        print('PAR_EIS no import')
#        file_handler = logging.FileHandler(FolderOrganizing.FolderOps.FindExpFolder('VERSASTAT').IndexDir.joinpath('PAR_EIS_fit_logfile.log'))
#
#    formatter    = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s')
#    file_handler.setFormatter(formatter)
#
#    # add file handler to logger
#    logger.addHandler(file_handler)
#    logger.warning('Started logging PAR EIS fit...')
#
# else:
#    logger = logging.getLogger(__name__)
#        # set log level
#    logger.setLevel(logging.INFO)
##
#    try:
#        from . import FolderOps
#
#        # define file handler and set formatter
#        file_handler = logging.FileHandler(FolderOps.FindExpFolder('VERSASTAT').IndexDir.joinpath('PAR_EIS_fit_logfile.log'))
#
#    except Exception as e:
#        print('Expception EIS import {0}, from {1}'.format(e,__name__))
#        from FolderOrganizing.FolderOps import *
#        print('PAR_EIS no import')
##        file_handler = logging.FileHandler(FolderOrganizing.FolderOps.FindExpFolder('VERSASTAT').IndexDir.joinpath('PAR_EIS_fit_logfile.log'))
#
#    formatter    = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : [%(lineno)d] %(message)s')
#    file_handler.setFormatter(formatter)
#
#    # add file handler to logger
#    logger.addHandler(file_handler)
#    logger.warning('Started logging PAR EIS fit...')


#    StartLogging()
#    runEIS = input('Want to run EIS postprocessing y/n? ...')
#    if runEIS == 'y':
#        EIS_postprocessing()
#    else:
#        pass
