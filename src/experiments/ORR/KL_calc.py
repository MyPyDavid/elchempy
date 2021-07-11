import sys
from pathlib import Path
from collections import namedtuple
from datetime import datetime
import numpy as np

from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines

import os
import multiprocessing
from functools import partial
from itertools import repeat
import pandas as pd

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper.ExtraInfo import EC_Properties


if __name__ == "__main__":
    pass

import logging

logger = logging.getLogger(__name__)


class KL_operations:
    def __init__(self, KL_data, ORR_dest_dir_file, electrode_properties={}):
        self.KL_data = KL_data
        self.ORR_dest_dir_file = ORR_dest_dir_file
        self.electrode_properties = electrode_properties

    def prepare_data(self, KL_data):
        _KL_meta = (
            KL_data[[i for i in KL_data.columns if KL_data[i].nunique() == 1]]
            .iloc[0]
            .to_dict()
        )
        KL_coeff = KL_coefficients()

        KLcoeffpH = KL_coeff.loc[
            (KL_coeff.pH == _KL_meta.get("pH", 99))
            & (KL_coeff.Electrolyte == _KL_meta.get("Electrolyte", "Unknown")),
            :,
        ]
        F, D, nu, C = (
            96485,
            KLcoeffpH.D0.values[0],
            KLcoeffpH.kVis.values[0],
            KLcoeffpH.C0.values[0],
        )
        # FIXME
        Area = self.electrode_properties.get("area_cm2", 1)
        # TODO WE_SA_collection_eff('PINE')


def KL_plots(KL_data, ORR_dest_dir_file):
    #        PAR_file, KLout, ORR_dest_dir,EvRHE, plot_KL=True):
    EvRHE = "E_AppV_RHE"
    # %%
    #    O2_out_ovv,dest_dir,dest_file = O2_out_ovv,dest_dir,dest_file
    #    F,D,nu,C,Area = 96485, 7.6E-05, 1.1E-02, 1.1E-06, WE_SA_collection_eff('PINE')['Disk_cm2']
    _KL_meta = (
        KL_data[[i for i in KL_data.columns if KL_data[i].nunique() == 1]]
        .iloc[0]
        .to_dict()
    )
    KL_coeff = KL_coefficients()

    KLcoeffpH = KL_coeff.loc[
        (KL_coeff.pH == _KL_meta.get("pH", 99))
        & (KL_coeff.Electrolyte == _KL_meta.get("Electrolyte", "Unknown")),
        :,
    ]
    F, D, nu, C, Area = (
        96485,
        KLcoeffpH.D0.values[0],
        KLcoeffpH.kVis.values[0],
        KLcoeffpH.C0.values[0],
        EC_Properties.WE_SA_collection_eff("PINE")[
            "Disk_cm2"
        ],  # FIX ME JUST GET CONSTANT VALUE
    )
    #    print('KL pars used:',F,D,nu,C,Area)
    #    print('KLoutCol:',KLout.columns)
    """C0  is  the  bulk  concentration  of  O2,  (C0  =1.2×10-6mol·cm-3),ν  is  the  kinematic viscosity of the electrolyte (ν=0.01 cm2·s-1),
    D0 is the diffusion coefficient of O2 in 0.1 M KOH (1.9 × 10-5 cm2·s-1)"""

    """ D: diffusion coefficient uncertain: In the current experiment the diffusion coefficient for oxygen in electrolyte is 7.6 x 10–5 cm2 s–1 with the four-electron mechanism,
     and 2.2 x 10–5 cm2 s–1 when the two-electron pathway is dominant. The measured signal from the electrochemical signal reflects both these processes.
     This value is higher than that seen in literature (about 1.4 x 10–5 cm2 s–1). This may possibly be explained by the extreme sensitivity of this parameter
     to the oxygen concentration in the electrolyte, which introduces uncertainty into the measurement.https://www.azom.com/article.aspx?ArticleID=15450
     """
    fstem = Path(_KL_meta.get("PAR_file", "")).stem
    ORR_KL_dir = ORR_dest_dir_file.joinpath("KL")
    ORR_KL_dir.mkdir(parents=True, exist_ok=True)

    #    KLout.to_excel(ORR_KL_dir.joinpath('KL_ORR_out.xlsx'))
    #    index_info_ORR_KL = {'PAR_file': PAR_file, 'DestFile': ORR_KL_fn, 'Type_output': 'ORR_Jkin_calc_KL_data',
    #                         'Type_exp': 'ORR'}
    #    Starting KL fit and analysis, using KL out filtered columns
    #    KL_data.loc[:, [i for i in KL_data.columns if not '_x' in i and not '_y' in i]]
    #    for swp, swgrp in KL_data.groupby('Sweep_Type'):
    #        for i, Erow in swgrp.iterrows():
    #            RpmL = []
    #            for Elec in ['I_Disk', 'I_Ring']:
    #                #                RpmL = []
    #                try:
    #                    RpmL = [[swp, Erow[EvRHE], Elec, Erow[q], (Erow[q]) ** -1, float(q.split('_')[-1]),
    #                             (float(q.split('_')[-1]) ** 0.5)] for q in Erow.index if Elec in q]
    #                #                    print(Erow,RpmL)
    #                #                    RPM not in rad == 2*np.pi/60)
    #                except Exception as e:
    #                    logger.error('KL_Plots ERROR row RPM K-L for {0}, because {1}'.format(ORR_KL_fn, e))
    #                    RpmL = []
    #                #                KLrow = pd.DataFrame(RpmL,columns=['Sweep_Type',EvRHE,'Electrode','I','I**-1','RPM','RPM_sqrt(w)'])
    #
    #                RpmLout.append(RpmL)
    #    #            Ring = [i for i in b.index if 'I_Ring' in i]
    #    KLorg = pd.DataFrame(data=[y for x in RpmLout for y in x],
    #                         columns=['Sweep_Type', EvRHE, 'Electrode', 'I', 'I**-1', 'RPM', 'RPM_sqrt(w)'])
    KL_rpms = KL_data.query("RPM_DAC > 0")
    KL_rpms = KL_rpms.assign(
        **{
            "RPM**-1": KL_rpms["RPM_DAC"] ** -1.0,
            "RPM_sqrt(w)": np.sqrt(KL_rpms["RPM_DAC"]),
            "RPM_sqrt(w)**-1": np.sqrt(KL_rpms["RPM_DAC"]) ** -1.0,
        }
    )
    #        I_inv = KL_data[Elec]**-1.0
    KL_rpms = KL_rpms.assign(
        **{
            f"{Elec}_I**-1": KL_rpms[Elec] ** -1.0
            for Elec in ["KL_I_Disk", "KL_I_Ring"]
        }
    )

    if KL_rpms.empty:
        logger.error("KL_Plots KLparsOout empty KL_rpms")
    #    KLorg = KLorg.loc[KLorg['Sweep_Type'] != 'NA', :]

    #    cath = KLorg.query('RPM > 1').groupby(['Sweep_Type','Electrode']).get_group(('cathodic','I_Disk'))
    #    ano = KLorg.query('RPM > 1').groupby(['Sweep_Type','Electrode']).get_group(('anodic','I_Disk'))
    #    fig,ax = plt.subplots()
    #    for rcath,rcathgrp in cath.groupby('RPM'):
    #        rcathgrp.plot(x=EvRHE,y='I',ax=ax,label=rcath)
    #    fig,ax = plt.subplots()
    _KLpars = []

    model_KL0_2ne = 1 / ((2) * (0.2 * F * C * D ** (2 / 3) * nu ** (-1 / 6)) * Area)
    model_KL0_4ne = 1 / ((4) * (0.2 * F * C * D ** (2 / 3) * nu ** (-1 / 6)) * Area)
    #%%
    for swp, swgrp in KL_rpms.groupby(["Sweep_Type"]):
        ORR_KL_fn_base = ORR_KL_dir.joinpath(f"KL_data_{fstem}_{swp}.xlsx")
        ORR_KL_fn = FileOperations.CompareHashDFexport(swgrp, ORR_KL_fn_base)
        for Elec in ["KL_I_Disk", "KL_I_Ring"]:
            #        Agr.to_excel(ORR_KL_gr_fn)
            #        if FolderOps.FileOperations.CompareHashDFexport(Agr,ORR_KL_gr_fn)[0] == False:
            fig, ax = plt.subplots()
            for rpmA, rpmAgrp in swgrp.groupby("RPM_DAC"):
                rpmAgrp.plot(
                    x=f"KL_{EvRHE}",
                    y=Elec,
                    ax=ax,
                    label=rpmA,
                    title=f"KL_{fstem}_{swp}_{Elec}",
                )

            plt.legend()
            plt.savefig(
                ORR_KL_dir.joinpath(f"KL_RPMs_{fstem}_{swp}_{Elec}_rpms.png"),
                dpi=100,
                bbox_inches="tight",
            )
            plt.close()
            #            if plot_KL == True:
            fig = plt.figure(figsize=(8, 8))
            fig.suptitle(f"KL_{fstem}_{swp}_{Elec}")
            ax = plt.subplot()
            _KL1_mean = []
            for E, Egr in swgrp.groupby(f"KL_{EvRHE}"):
                #            print(E)
                _Egr_meta = (
                    Egr[[i for i in Egr.columns if Egr[i].nunique() == 1]]
                    .iloc[0]
                    .to_dict()
                )
                _KLpar_row = _Egr_meta
                x, y = (
                    Egr["RPM_sqrt(w)**-1"].values,
                    Egr[f"{Elec}_I**-1"].values,
                )  # I**-1 in mA^-1
                # SAVING KL EXP DATA ###
                KLfitting = linregress(x, y)

                y_fit_2e = model_KL0_2ne * x + KLfitting[1]
                y_fit_4e = model_KL0_4ne * x + KLfitting[1]

                if np.abs(KLfitting[2]) > 0.94:
                    y_fit = KLfitting[0] * x + KLfitting[1]
                    #                === KOUTECKY-LEVICH EQUATION ===
                    n_e = 1 / (
                        (KLfitting[0])
                        * (0.2 * F * C * D ** (2 / 3) * nu ** (-1 / 6))
                        * Area
                    )

                    # Factor 0.62 if  w in rad/s, Factor 0.2 if  w in rpm
                    #                0.2nFcO2(DO2)2/3v−1/6
                    #                   jD = 0.62*ne*F*D**(2/3)*nu**(-1/6)*C*rotr**0.5
                    #    Factor_jD = 0.62*F*Area*D**(2/3)*nu**(-1/6)*C
                    #                === KOUTECKY-LEVICH EQUATION ===
                    #                n_e = 1/((KLfitting[0])*Factor_jD) # w in rad/s
                    #                print(KLfitting)
                    ax.scatter(x, y)
                    #                           label='%.2f %.2f'%(E,KLfitting[2]))
                    ax.plot(
                        x,
                        y_fit,
                        label="%.2f Vrhe, %.1f; i = $%.4f*\omega^{1/2}\/+\/%.2f$"
                        % (E, n_e, KLfitting[0], KLfitting[1]),
                    )
                    _rmplst = ", ".join([str(i) for i in Egr["RPM_DAC"].values])
                    _KL1_mean.append(KLfitting[1])
                    _KLpar_row.update(
                        {
                            "Electrode": Elec,
                            EvRHE: E,
                            "nElectrons": n_e,
                            "KL_slope": KLfitting[0],
                            "KL_intercept": KLfitting[1],
                            "KL_fitR": KLfitting[2],
                            "RPM_list": _rmplst,
                            "KL_data_file": ORR_KL_fn,
                            "KL_data_x": ", ".join([str(i) for i in x]),
                            "KL_data_y": ", ".join([str(i) for i in y]),
                            "KL_fit_y": ", ".join([str(i) for i in y_fit]),
                            "KL_fit_y_2e": ", ".join([str(i) for i in y_fit_2e]),
                            "KL_fit_y_4e": ", ".join([str(i) for i in y_fit_4e]),
                        }
                    )
                    _KLpars.append(_KLpar_row)
                #                    _KLpars.append([swp, Elec, E, n_e, KLfitting[0], KLfitting[1], KLfitting[2]])
                #               ax.legend()
                else:
                    #                print('No KL fit plot for %.1f Vrhe' %E)
                    pass

            ax.plot(
                x,
                model_KL0_2ne * x + np.mean(_KL1_mean),
                label="2e",
                c="g",
                alpha=0.4,
                **{"ls": "-.", "lw": 10},
            )
            ax.plot(
                x,
                model_KL0_4ne * x + np.mean(_KL1_mean),
                label="4e",
                c="r",
                alpha=0.4,
                **{"ls": "-.", "lw": 10},
            )
            ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
            #        plt.legend(), plt.show()
            ORR_KL_sweep_base = ORR_KL_dir.joinpath(f"KL_{fstem}_{swp}_{Elec}.xlsx")
            ORR_KL_sweep = FileOperations.CompareHashDFexport(swgrp, ORR_KL_sweep_base)
            plt.savefig(ORR_KL_sweep.with_suffix(".png"), dpi=100, bbox_inches="tight")
            #                    ORR_KL_dir.joinpath('KL_%s_%s.png'%(A[0],A[1]))
            plt.close()

    KLparsOut = pd.DataFrame(_KLpars)

    #    KLparsOut = pd.DataFrame(KLpars,
    #                             columns=['Sweep_Type', 'Electrode', EvRHE, 'nElectrons', 'KL_slope', 'KL_intercept',
    #                                      'KL_fitR'])
    ORR_KL_pars_fn_base = ORR_KL_dir.parent.joinpath(f"KL_pars_{fstem}.xlsx")
    ORR_KL_pars_fn = FileOperations.CompareHashDFexport(KLparsOut, ORR_KL_pars_fn_base)
    #    index_info_ORR_KL_pars = {'PAR_file': PAR_file, 'DestFile': ORR_KL_pars_fn, 'Type_output': 'ORR_Jkin_calc_KL_pars',
    #                              'Type_exp': 'ORR'}
    #    KLparsOut.to_excel(ORR_KL_pars_fn)
    #        legend = ax.legend(loc='upper center', shadow=True)
    try:
        #        for a,gr in KLparsOut.query('Electrode == "I_Disk"').groupby('Sweep_Type'):
        #        KLparsOut.query('Electrode == "I_Disk"').plot(x=EvRHE,y='nElectrons',kind='scatter',title=ORR_dest_dir.parts[-2]+'_Disk')
        for elec, elgr in KLparsOut.groupby("Electrode"):
            fig, ax = plt.subplots()
            for swp, grpswp in elgr.groupby("Sweep_Type"):
                ymax = 8 if "Disk" in elec else elgr.nElectrons.max() * 1.5
                color = (
                    "blue"
                    if "cathodic" in swp
                    else "red"
                    if "anodic" in swp
                    else "black"
                )
                grpswp.plot(
                    x=EvRHE,
                    y="nElectrons",
                    kind="scatter",
                    title=ORR_dest_dir_file.parts[-2] + "_{0}".format(elec),
                    ylim=(0, ymax),
                    ax=ax,
                    c=color,
                    label=swp,
                )
            plt.savefig(
                ORR_KL_dir.parent.joinpath("KL_nElec_{0}_{1}.png".format(fstem, elec)),
                dpi=100,
                bbox_inches="tight",
            )
            plt.close()
    except Exception as e:
        logger.error("KL_Plots KLparsOout PLotting error Disk, {0}".format(e))
    return KLparsOut
    #        KLparsOut.query('Electrode == "I_Disk"').plot(x=EvRHE,y='nElectrons',kind='scatter',title=ORR_dest_dir.parts[-2]+'_Disk',ylim=(0,8))
    #        print(KLparsOut.query('Electrode == "I_Disk"'))

    #        try:
    #            KLparsOut.query('Electrode == "I_Ring"').plot(x=EvRHE,y='nElectrons',kind='scatter',title=ORR_dest_dir.parts[-2]+'_Ring')
    #            plt.savefig(ORR_KL_dir.joinpath('KL_nElec_Ring.png'),dpi=100,bbox_inches='tight')
    #            plt.close()
    #        except Exception as e:
    #            logger.error('KL_Plots KLparsOout PLotting error Ring {0}'.format(e))

    # %%    return print('KL plotted: %s' %ORR_dest_dir)


#    KLrow.plot(x='RPM_sqrt(w)',y='I**-1',kind='scatter',title='%s_%s at %2.f $V_{RHE}$'%(swp,Elec,Erow[EvRHE]))
#    plt.show()
#    plt.close()
def KL_coefficients():
    #    F,D,nu,C,Area = 96485, 7.6E-05, 1.1E-02, 1.1E-06, WE_SA_collection_eff('PINE')['Disk_cm2']
    KL_ORR_Coeff = pd.DataFrame(
        {
            "pH": [0.3, 13, 1],
            "Electrolyte": ["0.5MH2SO4", "0.1MKOH", "0.1MH2SO4"],
            "C0": [1.1e-06, 1.05e-06, 1.1e-06],
            "D0": [1.4e-05, 1.9e-05, 1.4e-05],
            "kVis": [0.01, 0.01, 0.01],
        }
    )
    """C0  is  the  bulk  concentration  of  O2,  (C0  =1.2×10-6mol·cm-3),ν  is  the  kinematic viscosity of the electrolyte (ν=0.01 cm2·s-1),
    D0 is the diffusion coefficient of O2 in 0.1 M KOH (1.9 × 10-5 cm2·s-1)"""

    """cO2 is the concentration of dissolved oxygen (1.1×10−6 mol cm−3) [38], DO2 its diffusion coefficient (1.4×10−5 cm2 s−1) [38],
    ν the kinematic viscosity (0.01 cm2 s−1 for sulfuric acid) [38], F the Faraday constant and n the apparent number of electrons transferred per molecule of O2 in the overall reaction.
    The constant 0.2 is used when ω is expressed in revolutions per minute [39].
    Ye, S. & Vijh, A.K. J Solid State Electrochem (2005) 9: 146. https://doi.org/10.1007/s10008-004-0567-0"""

    """ D: diffusion coefficient uncertain: In the current experiment the diffusion coefficient for oxygen in electrolyte is 7.6 x 10–5 cm2 s–1 with the four-electron mechanism,
     and 2.2 x 10–5 cm2 s–1 when the two-electron pathway is dominant. The measured signal from the electrochemical signal reflects both these processes.
     This value is higher than that seen in literature (about 1.4 x 10–5 cm2 s–1). This may possibly be explained by the extreme sensitivity of this parameter
     to the oxygen concentration in the electrolyte, which introduces uncertainty into the measurement.https://www.azom.com/article.aspx?ArticleID=15450
     """
    return KL_ORR_Coeff


#   i_Diff_lim = 0.62*ne*F*D**(2/3)*nu**(-1/6)*C*rotr**0.5
#       LEVICH (KOUTECKY) PLOTS
# I L = ( 0.620 ) n F A D 2 3 ω 1 2 v − 1 6 C {\displaystyle I_{L}=(0.620)nFAD^{\frac {2}{3}}\omega ^{\frac {1}{2}}v^{\frac {-1}{6}}C} {\displaystyle I_{L}=(0.620)nFAD^{\frac {2}{3}}\omega ^{\frac {1}{2}}v^{\frac {-1}{6}}C}
# IL is the Levich current (A)
# n is the number of moles of electrons transferred in the half reaction (number)
# F is the Faraday constant (C/mol)
# A is the electrode area (cm2)
# D is the diffusion coefficient (see Fick's law of diffusion) (cm2/s)
# ω is the angular rotation rate of the electrode (rad/s)
# v is the kinematic viscosity (cm2/s)
# C is the analyte concentration (mol/cm3
