#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:03:49 2020

@author: zmg
"""

from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
from scipy.stats import linregress
import seaborn as sns


from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper.FindSampleID import GetSampleID

from file_py_helper.PostChar import (
        SampleSelection,
        Characterization_TypeSetting,
        SampleCodesChar,
    )



from .collect_load import Load_from_Indexes
EC_index, SampleCodes = Load_from_Indexes.get_EC_index()

from ..experiments.ORR import ORR_analyze_scans

OriginColor = FindExpFolder().LoadOriginColor()

import logging
logger = logging.getLogger(__name__)

from . import EvRHE, EXP_FOLDER



def make_subdir(DD):
    DD.mkdir(parents=True, exist_ok=True)
    return DD


def N2export_to_Origin(postOVVout, ORR_pars):
    PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")

    SeriesIDs = [
        SampleSelection.Series_CB_paper,
        SampleSelection.Series_Porhp_SiO2,
        SampleSelection.Series_Co_PANI,
        SampleSelection.Series_ML_SiO2,
        {"name": "all_EIS"},
    ]
    SeriesID_set = SeriesIDs[1]
    PPDN2 = PostDestDir.joinpath("ORR_{0}".format(SeriesID_set["name"]))
    PPDN2.mkdir(parents=True, exist_ok=True)

    PPDCV = PPDN2.joinpath("O2_ORR")
    PPDCV.mkdir(parents=True, exist_ok=True)
    #    Cdl_pars
    ORR_parsSIO2 = ORR_pars.loc[ORR_pars.SampleID.isin(SeriesID_set["sIDs"])]
    #    N2_CV_index = postOVVout.groupby('Type_output').get_group('N2_CV')
    ORR_parsSIO2 = ORR_parsSIO2.loc[ORR_pars.SampleID.isin(SeriesID_set["sIDs"])]
    #    & ORR_pars.Source.str.contains('ExpDir')]
    Acid_slice = ORR_parsSIO2.query('pH < 7  & postAST == "no"').dropna(
        subset=["ECexp"]
    )
    jOS5 = Acid_slice.query('SampleID == "JOS5"')

    ORR_acid_no = Acid_slice.loc[Acid_slice.ECexp.str.contains("05-06")]

    "ORR_Jkin_calc_KL_data"
    "ORR_Jkin_calc_RRDE"
    "ORR_Jkin_calc_RRDE_Chrono"
    ORR_acid_no.PAR_file
    postOVVout.loc[
        postOVVout.PAR_file.isin(
            list(ORR_acid_no.PAR_file.unique())
            + [i.replace("Disc", "Ring") for i in ORR_acid_no.PAR_file.unique()]
        )
    ]

    postOVVout.loc[postOVVout.PAR_file.isin(list(ORR_acid_no.PAR_file.unique()))]

    postOVVout.loc[
        postOVVout.PAR_file.isin(list(ORR_acid_no.PAR_file.unique()))
    ].SampleID.unique()
    selpost = postOVVout.loc[
        postOVVout.PAR_file.isin(list(ORR_acid_no.PAR_file.unique()))
    ]
    #%%
    Jkin_PARS_data, Jkin_PARS_path = [], PPDCV.joinpath("PostEC_ORR_Jkin_pars.xlsx")
    for sID, gr in selpost.groupby(["SampleID"]):
        #        sID,gr
        for (exp, elec, PARfile), expgrp in gr.groupby(
            ["ECexp", "Electrolyte", "PAR_file"]
        ):
            exp, elec, expgrp

            PPDelec = make_subdir(PPDCV.joinpath(elec))

            sID, metalprec = expgrp.SampleID.unique()[0], expgrp.MetalPrec.unique()[0]
            metal = metalprec[0:2]
            xldate = exp.split("_")[-1]
            PARf_name = Path(PARfile).name
            #            extra = '_'.join(Path(filename).stem.split('_')[-2:])
            typegrp = expgrp.groupby("Type_output")
            KLdataout = [1]

            def export_Jkin_RRDE(typegrp):
                Jkin_RRDE_fls_cath = typegrp.get_group(
                    "ORR_Jkin_calc_RRDE_Chrono"
                ).DestFile.values[0]

                Jkin_RRDE_an = pd.read_excel(
                    Jkin_RRDE_fls_cath.replace("cathodic", "anodic"), index_col=[0]
                )
                #                parfKLpars = Path(typegrp.get_group('ORR_Jkin_calc_Tafel').PAR_file.values[0]).name
                #            print(f'{sID} {exp} {parf}')
                Jkin_RRDE_cath = pd.read_excel(Jkin_RRDE_fls_cath, index_col=[0])
                #            Jkin = Jkin_read[['E_AppV_RHE_disk', 'Jcorr','Jkin_min', 'Frac_H2O2', 'J_ring', 'n_ORR','J_N2_scan']].set_index(EvRHE+'_disk')

                ORR_Jkin_pars = []

                def output_file(Jkin_RRDE):
                    for n, gr in Jkin_RRDE.groupby(
                        ["SampleID", "RPM_DAC", "Sweep_Type_disk"]
                    ):
                        EvRHE = [
                            i
                            for i in gr.columns
                            if "E_APPV_RHE" in i.upper() and not "RING" in i.upper()
                        ][0]
                        AST_par = GetSampleID.determine_postAST_from_filename(
                            gr.BaseName.unique()[0]
                        )

                        JKIN = gr[
                            [
                                "E_AppV_RHE_disk",
                                "Jcorr",
                                "Jkin_min",
                                "Frac_H2O2",
                                "J_ring",
                                "n_ORR",
                                "J_N2_scan",
                            ]
                        ].set_index(EvRHE)
                        jkinf = f"{n[-1]}_{n[1]}_{sID}_{metal}_{xldate}_Jkin.xlsx"
                        JKIN.to_excel(PPDelec.joinpath(jkinf))

                        ORR_Jkinpars_row = ORR_analyze_scans.ORR_extract_pars(
                            gr, postAST=AST_par
                        )
                        ORR_Jkinpars_row.update(
                            {"RRDE_DataFile": PPDelec.joinpath(jkinf)}
                        )
                        return ORR_Jkinpars_row

                #                ORR_Jkin_pars.update( output_file(Jkin_RRDE_an))
                #                ORR_Jkin_pars.update( output_file(Jkin_RRDE_cath))
                ORR_Jkin_pars.append(output_file(Jkin_RRDE_an))
                ORR_Jkin_pars.append(output_file(Jkin_RRDE_cath))
                return ORR_Jkin_pars

            Jkin_PARS_data.append(export_Jkin_RRDE(typegrp))
            #                print(gr.File.unique()[0])

            def export_KLdata(typegrp, sID, metalprec, metal, PARf_name, xldate):
                KLdata_fls = typegrp.get_group("ORR_Jkin_calc_KL_data").DestFile.values[
                    0
                ]
                #                parf = Path(typegrp.get_group('ORR_Jkin_calc_KL_data').PAR_file.values[0]).name
                KLdata = pd.read_excel(KLdata_fls, index_col=[0])
                for swp, sgr in KLdata.groupby("Sweep_Type"):
                    #                    swp,sgr
                    #                    sID,metalprec = sgr.SampleID.unique()[0], sgr.MetalPrec.unique()[0]
                    print(f"{sID} {exp}\nPAR: {PARf_name}")
                    #                    metal = metalprec[0:2]
                    #                    xldate = exp.split('_')[-1]
                    #                    extra = '_'.join(Path(PARf_name).stem.split('_')[-2:])
                    ShortLabel = f"{sID}_{metal}_RRDE"
                    fpath = f"{swp}_{ShortLabel}_{xldate}_{exp}.png"
                    #                    sgr.plot(x=EvRHE,y=[i for i in sgr.columns if 'Ring' in i],title=f'{swp} {sID},\n {exp}')
                    sgr.plot(
                        x=EvRHE,
                        y=[i for i in sgr.columns if "Disk" in i],
                        title=f"{swp} {sID},\n {exp}",
                    )
                    plt.savefig(PPDelec.joinpath(fpath), bbox_inches="tight", dpi=100)
                    plt.close()
                    Elec_info = {i: sgr[i].unique()[0] for i in ["Electrolyte", "pH"]}
                    RRDE = sgr.drop(
                        columns=["Electrolyte", "pH", "Sweep_Type"]
                    ).set_index(EvRHE)
                    sortcols = sorted(
                        RRDE.columns, key=lambda x: float(x.split("_")[-1])
                    )
                    RRDE = RRDE[
                        [i for i in sortcols if "Disk" in i]
                        + [i for i in sortcols if "Ring" in i]
                    ]
                    xlpath_RRDE = f"{swp}_{ShortLabel}_{xldate}.xlsx"
                    RRDE.to_excel(PPDelec.joinpath(xlpath_RRDE))
                    xlpath_KL_pars = f"{swp}_{sID}_{metal}_KL_pars_{xldate}.xlsx"
                    xlpath_KL_data = f"{swp}_{sID}_{metal}_KL_data_disk_{xldate}.xlsx"
                    KL_pars, KL_data = KLcalc_data(sgr, Elec_info)
                    with pd.ExcelWriter(PPDelec.joinpath(xlpath_KL_pars)) as writer:
                        KL_pars.loc[KL_pars[EvRHE].isin(np.arange(0, 1, 0.1))].to_excel(
                            writer, sheet_name="short"
                        )
                        KL_pars.to_excel(writer, sheet_name="full")

                    with pd.ExcelWriter(PPDelec.joinpath(xlpath_KL_data)) as writer:
                        for E, Egr in KL_data.loc[
                            KL_data["E_RHE"].isin(np.arange(0, 1, 0.1))
                        ].groupby("E_RHE"):
                            Egr.to_excel(writer, sheet_name=str(np.round(E, 2)))

            #                        N2cdl.iloc[1000:].to_excel(writer, sheet_name='anodic')
            export_KLdata(typegrp, sID, metalprec, metal, PARf_name, xldate)
            KLparsout = []

            def export_KLpars(typegrp, sID, metalprec, metal, PARf_name, xldate):
                KLpars_fls = typegrp.get_group("ORR_Jkin_calc_KL_pars").DestFile.values[
                    0
                ]
                #                sID,metalprec = expgrp.SampleID.unique()[0], expgrp.MetalPrec.unique()[0]
                #                metal = metalprec[0:2]
                #            print(f'{sID} {exp} {parf}')
                KLpars = pd.read_excel(KLpars_fls, index_col=[0])
                for (swp, electrode), KLsgr in KLpars.groupby(
                    ["Sweep_Type", "Electrode"]
                ):
                    swp, KLsgr
                    slKLpars = f"{sID}_{metal}_KLcalc"
                    KL = KLsgr.drop(columns=["Sweep_Type", "Electrode"]).set_index(
                        EvRHE
                    )
                    xlpathKL = f"{swp}_{electrode}_{slKLpars}_{xldate}.xlsx"
                    KL.to_excel(PPDelec.joinpath(xlpathKL))

            export_KLpars(typegrp, sID, metalprec, metal, PARf_name, xldate)
    Jkin_PARS = pd.concat(
        [pd.DataFrame(i, index=[0]) for sl in Jkin_PARS_data for i in sl],
        sort=False,
        ignore_index=True,
    )
    Jkin_PARS.to_excel(Jkin_PARS_path)


#%%

#            Jkin.RPM_DAC.unique()[0]
def plotxy(x, y):
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    plt.show()


def KLcalc_data(sgr, Elec_info):
    pH = Elec_info.get("pH", 0)
    KLcoeffpH = KL_coefficients(pH)
    electrode_info = {
        "Electrode Type": "PINE",
        "CollEff": 0.38,
        "Disk_cm2": 0.2376,
        "Ring_cm2": 0.2356,
    }
    F, D, nu, C, Area = (
        96485,
        KLcoeffpH.D0,
        KLcoeffpH.kVis,
        KLcoeffpH.C0,
        electrode_info.get("Disk_cm2"),
    )

    KL_pars, KL_data = [], []
    #    sgr.loc[sgr[EvRHE] == 0.4]
    for E, Egr in sgr.groupby(EvRHE):
        Erow = Egr.iloc[0].to_dict()
        disk = {key: val for key, val in Erow.items() if "Disk" in key}
        disk_rpms = [
            (float(i.split("_")[-1]), val)
            for i, val in disk.items()
            if float(i.split("_")[-1]) > 0
        ]
        disk_rads = [(2 * np.pi * i[0] / 60, i[1]) for i in disk_rpms]

        #        _model_disk_rpms =
        #        ['Sweep_Type', EvRHE, 'Electrode', 'I', 'I**-1', 'RPM', 'RPM_sqrt(w)']
        #        RpmL = [[swp, Erow[EvRHE], Elec, Erow[q], (Erow[q]) ** -1, float(q.split('_')[-1]),
        #                             (float(q.split('_')[-1]) ** 0.5)] for q in Erow.index if Elec in q]
        #        print(E)
        #       ω_rad/s = 2π*RPM / 60
        x = np.array([i[0] ** -0.5 for i in disk_rads])  # 'RPM_sqrt(w)**-1' in rad/s
        y = np.array([i[1] ** -1 for i in disk_rads])  # ['I**-1'].
        #        x, y = Egr['RPM_sqrt(w)**-1'].values, Egr['I**-1'].values  # I**-1 in mA^-1
        #        plotxy(x,y)
        # SAVING KL EXP DATA ###
        KLfitting = linregress(x, y)
        #        if np.abs(KLfitting[2]) > 0.94:
        y_fit = KLfitting[0] * x + KLfitting[1]
        #                === KOUTECKY-LEVICH EQUATION ===
        n_e = 1 / (
            (KLfitting[0]) * (0.62 * F * C * D ** (2 / 3) * nu ** (-1 / 6)) * Area
        )

        KLfitting_2e = 1 / ((2) * (0.62 * F * C * D ** (2 / 3) * nu ** (-1 / 6)) * Area)
        KLfitting_4e = 1 / ((4) * (0.62 * F * C * D ** (2 / 3) * nu ** (-1 / 6)) * Area)

        _y_fitting_2e = KLfitting_2e * x + KLfitting[1]
        _y_fitting_4e = KLfitting_4e * x + KLfitting[1]

        # Factor 0.62 if  w in rad/s, Factor 0.2 if  w in rpm
        #                0.2nFcO2(DO2)2/3v−1/6
        #                   jD = 0.62*ne*F*D**(2/3)*nu**(-1/6)*C*rotr**0.5
        #    Factor_jD = 0.62*F*Area*D**(2/3)*nu**(-1/6)*C
        #                === KOUTECKY-LEVICH EQUATION ===
        #                n_e = 1/((KLfitting[0])*Factor_jD) # w in rad/s
        #                print(KLfitting)
        #            ax.scatter(x, y)
        #                           label='%.2f %.2f'%(E,KLfitting[2]))
        #            ax.plot(x, y_fit, label='%.2f Vrhe, %.1f; i = $%.4f*\omega^{1/2}\/+\/%.2f$' % (
        #            E, n_e, KLfitting[0], KLfitting[1]))
        KL_data_row = {
            "E_RHE": [E for i in disk_rpms],
            "RPM": [i[0] for i in disk_rpms],
            "rad_s**-0.5": x,
            "I_**-1": y,
            "I": [i[1] for i in disk_rpms],
            "KL_fit_I**-1": y_fit,
            "KL_2e_fit_I**-1": _y_fitting_2e,
            "KL_4e_fit_I**-1": _y_fitting_4e,
        }
        KL_pars_row = {
            EvRHE: E,
            "nE": n_e,
            "KL_slope": KLfitting[0],
            "KL_intercept": KLfitting[1],
            "KL_r": KLfitting[2],
            "KL_2e_slope": KLfitting_2e,
            "KL_4e_slope": KLfitting_4e,
        }
        KL_pars.append(KL_pars_row)
        KL_data.append(pd.DataFrame(KL_data_row))
    return pd.DataFrame(KL_pars), pd.concat(KL_data)
    #               ax.legend()


#        else:
#            #                print('No KL fit plot for %.1f Vrhe' %E)
#            pass


def KL_coefficients(pH):
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
    return KL_ORR_Coeff.loc[pH]


#
#            for swp,Jkgr in Jkin_read.groupby(['Sweep_Type']):
#                swp,Jkgr
#                slKLpars = f'{sID}_{metal}_KLcalc'
#                KL = KLsgr.drop(columns=['Sweep_Type','Electrode']).set_index(EvRHE)
#                xlpathKL = f'{swp}_{electrode}_{slKLpars}_{xldate}.xlsx'
#                KL.to_excel(PPDelec.joinpath(xlpathKL))
#
#            ORR_Jkin_calc_Tafel
#%%
##            pd.read_excel(KLdata.DestFile.values[0],index_col=[0])
#            if expgrp.Cdl_CV_data_files.values[0] == expgrp.Cdl_CV_data_files.values[-1]:
#                CVfiles = expgrp.Cdl_CV_data_files.values[0]
#            else:
#                CVfiles = expgrp.Cdl_CV_data_files.values[-1]
#            if len(CVfiles) == 12:
#                CV12 = CVfiles
#                CVfiles = [i for i in CV12 if 'standard' in i]
#                if len(CVfiles) != 6:
#                    CVprob = sorted(CVfiles)
#
#                    if len(CVprob) == 12:
#                        CVfiles = CVprob[0:6]
#                    else:
#                        print(f'12 to 6 problem  {exp}, len is {len(CVfiles)} ')
#
#            if len(CVfiles) != 6:
#                print(f'len problem {exp}, len is {len(CVfiles)}')
##                print(f'len problem {exp}, len is {len(CVfiles)}')
#                lentest = [i for i in  expgrp.Cdl_CV_data_files.values if len(i) == 6]
#                if lentest:
#                    if lentest[0] == lentest[-1]:
#                        CVfiles = lentest[0]
#                    else:
#                        CVfiles = lentest[-1]
#                else:
#                    CVfiles = []
#
#            if len(CVfiles) != 6:
#                print(f'DOUBLE len problem {exp}, len is {len(CVfiles)}')
##            def export_CVset(CVfiles,expgrp):
##            [(Path(i).stem.split('_')[-2],i) for i in  CVfiles if not 'first' in i]
#            if CVfiles:
#                cvsr = sorted([(Path(i).stem.split('_')[-2],Path(i)) for i in  CVfiles if not 'first' in i], reverse=True)
#                CVexp = FileHelper.FileOperations.ChangeRoot_DF(pd.DataFrame(cvsr,columns=['ScanRate','File']),[],coltype='string')
#                pdlst,pdmg = [], pd.DataFrame()
#                N2origin = pd.DataFrame()
#                for n,r in CVexp.iterrows():
#                    try:
#                        SR = pd.read_excel(r.File,index_col=[0]).rename(columns={'j A/cm2' : f'CV_{r.ScanRate}'}).reset_index()
##                        SR[f'CV_{r.ScanRate}'] = SR[f'CV_{r.ScanRate}']*1000 # TODO FORGET about mA
#                        pdlst.append(SR)
#                        if pdmg.empty:
#                            pdmg = SR
#                        else:
#                            pdmg = pd.merge(pdmg, SR[f'CV_{r.ScanRate}'], left_index=True, right_index=True)
#        #                        pd.merge(pdmg,SR, on = EvRHE)
#                    except:
#                        print(f'read problem {exp} ')
#                N2origin = pdmg.set_index(EvRHE).drop('index',axis=1)
#        #                [pd.read_excel(i[1]) for i in cvsr]
#                PPDelec = make_subdir(PPDCV.joinpath(elec))
#
#
#                sID,metalprec = expgrp.SampleID.unique()[0], expgrp.MetalPrec.unique()[0]
#                metal = metalprec[0:2]
#                xldate = exp.split('_')[-1]
#                extra = '_'.join(Path(filename).stem.split('_')[-2:])
#                ShortLabel = f'{sID}_{metal}_CVs'
#                fpath = f'{ShortLabel}_{xldate}_{extra}_{exp}.png'
#
#                N2origin.plot(title=f'{ShortLabel},{extra}\n{exp}')
#                plt.savefig(PPDelec.joinpath(fpath),bbox_inches='tight',dpi=100)
#                plt.close()
#
#                xlpath = f'{ShortLabel}_{xldate}.xlsx'
#                if PPDelec.joinpath(xlpath).is_file():
#                    xlpath = f'{ShortLabel}_{xldate}_{extra}.xlsx'
#                N2cdl = calc_cdl(N2origin)
#                N2cdl.to_excel(PPDelec.joinpath(xlpath))


class ExportfromORR:

    ORR_kin_cols = [
        "E_onset",
        "E_half",
        "J_diff_lim",
        "Jkin_075",
        "Jkin_080",
        "TSa_l",
        "TSb_l",
        "TSa_h",
        "TSb_h",
        "Jring_050",
        "FracH2O2_050",
    ]

    def __init__(self):
        pass

    def ORR_standard_plots(ORR_pars):
        fig, ax = plt.subplots()
        ORR_pars.query("(pH < 15) & (RPM > 900)").plot(
            y="J_diff_lim",
            x="Loading_cm2",
            kind="scatter",
            logy=0,
            c="BET_cat_agg",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars.query("(pH < 14) & (RPM > 900)").plot(
            x="TSb_l",
            y="E_onset",
            kind="scatter",
            xlim=(0.5, 1.5),
            logx=True,
            ylim=(0.5, 1),
            c="pH",
            colormap="rainbow_r",
            ax=ax,
        )
        #        plot(y='J_diff_lim',x='Loading_cm2',kind='scatter',logy=0,c='BET_cat_agg',colormap='viridis',ax=ax)
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars.query("(pH < 14) & (RPM > 900)").plot(
            x="TSa_l",
            y="E_onset",
            kind="scatter",
            xlim=(10, 200),
            ylim=(0.5, 1),
            c="pH",
            colormap="rainbow_r",
            ax=ax,
        )
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars.query("(pH < 7) & (RPM > 900)").plot(
            y="J_diff_lim",
            x="Loading_cm2",
            kind="scatter",
            logy=0,
            c="BET_cat_agg",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ORR_pars.query(
            '(pH < 15) & (RPM > 1100) & (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
        ).plot(
            y="E_onset",
            x="N_content",
            kind="scatter",
            logy=0,
            c="pH",
            colormap="viridis",
            ax=ax,
            xlim=(0, 10),
        )
        plt.show()
        plt.close()
        fig, ax = plt.subplots()
        ORR_pars.query("(pH < 15) & (RPM > 900)").plot(
            y="E_half",
            x="AD/AG",
            kind="scatter",
            logy=0,
            c="pH",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()
        fig, ax = plt.subplots()
        ORR_pars.query("(pH < 7) & (RPM > 900)").plot(
            y="Jkin_075",
            x="BET_cat_agg",
            kind="scatter",
            logy=0,
            c="N_content",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()
        fig, ax = plt.subplots()
        ORR_pars.query("(pH < 7) & (RPM > 1200)").plot(
            y="Jkin_075",
            x="D1_pop_Fe_wt",
            kind="scatter",
            logy=0,
            c="N_content",
            colormap="viridis",
            ax=ax,
        )
        plt.show()
        plt.close()

    def LoadingSeries(ORR_pars):
        ORR_Loadingpars = pd.concat(
            [
                i[1]
                for i in [
                    (n, gr, gr.Loading_cm2.unique())
                    for n, gr in ORR_pars.query("RPM > 1000").groupby(
                        ["SampleID"] + SampleSelection.EC_exp_cols[0:-1]
                    )
                    if gr.Loading_cm2.nunique() > 2
                ]
            ]
        )
        ORR_Loadingpars.SampleID.unique()
        corr_loading = []
        for nID, grID in ORR_Loadingpars.groupby("SampleID"):
            nID, grID
            _normal = grID.query("(Loading_cm2 < 0.5) & (Loading_cm2 > 0.3)")
            for col in ["TSa_l", "TSb_l", "TSa_h", "TSb_h"]:
                grID.plot(
                    x="Loading_cm2",
                    y=col,
                    kind="scatter",
                    title="{0}, {1}".format(nID, col),
                    logy=0,
                )
                _lin = scipy.stats.linregress(
                    grID["Loading_cm2"].values, grID[col].values
                )
                corr_loading.append(
                    {
                        "SampleID": nID,
                        "X_col": "Loading_cm2",
                        "Y_col": col,
                        "slope": _lin.slope,
                        "intercept": _lin.intercept,
                        "r_val": _lin.rvalue,
                    }
                )
                if "Jkin_" in col:
                    _lin = scipy.stats.linregress(
                        grID["Loading_cm2"].values, np.log10(grID[col].values)
                    )
                    corr_loading.append(
                        {
                            "SampleID": nID,
                            "X_col": "Loading_cm2",
                            "Y_col": col + "_log",
                            "slope": _lin.slope,
                            "intercept": _lin.intercept,
                            "r_val": _lin.rvalue,
                        }
                    )
        #           slope, intercept, r_value, p_value, std_err
        Loading_series_result = pd.DataFrame(corr_loading).sort_values("r_val")
        # TODO small interpretation of results:
        """ From the correlations with ORR measurements at several loadings the following can be interpreted.
        It seems that log(Jkin) follows upward trend with loading.  TSa, F"""

    #%%
    def ORR_paper(ORR_pars,ORR_CB_paper): # FIX ME
        """Parameters of all measurements are plotted together and a linear fit is done that show for ('Jkin' with 'E_onset') specifically
        figures are saved."""
        PDD_ORR = FindExpFolder("VERSASTAT").DestDir.joinpath(
            "PostEC/{0}/ORR".format(SampleSelection.Series_CB_paper["name"])
        )
        PDD_ORR.mkdir(exist_ok=True, parents=True)
        #        ORR_data_df = pd.concat([pd.read_excel(i)  for i in ORR_CB_paper.query('RPM > 1000').RRDE_DataFile.unique() if Path(i).is_file()] )

        for xtest in [i for i in SampleSelection.EC_ORR_kin_par_cols if "Jkin" in i]:
            fig, ax = plt.subplots(figsize=(10, 8))
            logy_set = False if not "Jkin" in xtest else True
            #            ORR_test = ORR_pars.query('(pH < 14) & (RPM > 1000)  & (E_onset > 0.71) & (postAST == "no") & ((Loading_cm2 < 0.4) & (Loading_cm2 > 0.35))').drop_duplicates(subset=[xtest,'E_onset'])
            #            ORR_test =
            ORR_test = ORR_pars.query(
                '(pH < 14) & (RPM > 1000)  & (E_onset > 0.71) & (postAST == "no")'
            ).drop_duplicates(subset=[xtest, "E_onset"])
            #            ORR_test = ORR_pars.query('(pH < 14) & (RPM > 1000)  & (E_onset > 0.2) & (postAST != "no") & ((Loading_cm2 < 0.4) & (Loading_cm2 > 0.35))').drop_duplicates(subset=[xtest,'E_onset'])
            ORR_test = ORR_test.dropna(subset=[xtest, "E_onset"])
            ORR_ignored = ORR_pars.loc[
                ~ORR_pars.index.isin(ORR_test.index)
                & (ORR_pars.SampleID.isin(SampleCodes.SampleID))
            ]
            for Sern, Sergr in ORR_test.groupby("SeriesID"):
                mkset = Characterization_TypeSetting.Series_Colors(Sern)[
                    "marker"
                ]
                #                Sergr.plot(y=xtest,x='E_onset',kind='scatter',logy=logy_set, ax=ax, marker=mkset, markersize = 100)
                ax.scatter(
                    Sergr["E_onset"], Sergr[xtest], marker=mkset, s=100, label=Sern
                )

            #                c='pH' colormap='rainbow_r'
            if logy_set:
                logj, Eonset_fit, log_exp = (
                    np.log10(ORR_test[xtest].values),
                    ORR_test["E_onset"].values * 1000,
                    10,
                )
                lin2 = scipy.stats.linregress(y=logj, x=ORR_test["E_onset"].values)
                y_fit = 10 ** (np.arange(0.1, 1, 0.01) * lin2.slope + lin2.intercept)
                ax.set_xlim(0.5, 0.95)
                ax.set_ylim(0.005, 500)
            else:
                logj, Eonset_fit, log_exp = (
                    ORR_test[xtest].values,
                    ORR_test["E_onset"].values,
                    1,
                )
                lin2 = scipy.stats.linregress(y=logj, x=Eonset_fit)
                y_fit = np.arange(0.1, 1, 0.01) * lin2.slope + lin2.intercept
            slope, intercept, rval, pval, stder = scipy.stats.linregress(
                y=Eonset_fit, x=logj
            )
            if logy_set:
                ax.set_yscale("log")
            ax.scatter(ORR_ignored["E_onset"], ORR_ignored[xtest], alpha=0.05)
            ax.plot(np.arange(0.1, 1, 0.01), y_fit, label=f"Slope: {slope:.1f} mV/dec")
            ax.legend(loc="best")
            #        ORR_test.plot(y='meta_logJkin75',x='E_onset',kind='line',ls='--', logy=True,xlim=(0.5,1),ylim=(0.1,40),ax=ax)
            plt.show()
            plt.savefig(
                PDD_ORR.joinpath(f"{xtest}_Eonset_ovv.png"),
                bbox_inches="tight",
                dpi=300,
            )
            plt.close()

        acid = (
            ORR_CB_paper.query("RPM > 1000")
            .groupby(SampleSelection.EC_exp_cols)
            .get_group(("no", "0.5MH2SO4", 0.3, "O2", 1500, 0.379))
        )

        for nm, gr in ORR_CB_paper.query("RPM > 1000").groupby(
            SampleSelection.EC_exp_cols
        ):
            PDD_ORRgr = FindExpFolder("VERSASTAT").DestDir.joinpath(
                "PostEC/{0}/ORR/{1}".format(
                    SampleSelection.Series_CB_paper["name"],
                    "_".join([str(i) for i in nm]),
                )
            )
            PDD_ORRgr.mkdir(exist_ok=True, parents=True)
            print(nm)
            print(gr.SampleID.unique())
            missing = []

            for IDnm, IDgr in gr.groupby("SampleID"):
                fig, ax = plt.subplots(figsize=(8, 8))
                ax2 = ax.twinx()
                for PFnm, PFgr in IDgr.groupby("PAR_file"):
                    datafiles = [
                        pd.read_excel(i)
                        for i in PFgr.RRDE_DataFile.to_list()
                        if Path(i).is_file()
                    ]
                    #                    datafiles = ORR_data_df.query('File_disk == PFgr.RRDE_DataFile.to_list()[0]')
                    if datafiles:
                        RRDE_PF = pd.concat(datafiles)
                        RRDE_PF.to_excel(
                            PDD_ORRgr.joinpath(
                                "{0}_{1}_{2}.xlsx".format(
                                    IDnm,
                                    Path(PFnm).stem,
                                    PFgr.EXP_date.to_list()[0].strftime("%Y-%m-%d"),
                                )
                            )
                        )
                    else:
                        missing.append(
                            [
                                nm,
                                IDnm,
                                PFnm,
                                [
                                    i
                                    for i in PFgr.RRDE_DataFile.to_list()
                                    if not Path(i).is_file()
                                ],
                            ]
                        )
                    #                RRDE_PF.plot(x='E(V)',y=['I(A)_disk','I(A)_ring'])
                    lbl = "{0} ({1}) {2}[{3}]".format(
                        PFgr.SampleCode_x.to_list()[0],
                        PFgr.SampleID.to_list()[0],
                        Path(PFnm).stem,
                        PFgr.EXP_date.to_list()[0].strftime("%Y-%m-%d"),
                    )
                    RRDE_PF.plot(x="E(V)", y="Jcorr", ax=ax, label=lbl)
                    RRDE_PF.plot(x="E(V)", y="Frac_H2O2", ax=ax2, ls="--", legend=False)
                ax2.set_ylim(0, 10)
                ax2.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))

    def ORR_plotting_parameters1(postOVVout): # FIXME
        #%% ===== MAKE PLOTS of Parameters versus E v RHE ======
        #        OnlyRecentMissingOVV = run_PAR_DW.ECRunOVV(load=1).index
        ORR_pars = Load_from_Indexes.ORR_pars_OVV(
            pd.DataFrame(), pd.DataFrame(), reload=False
        )  # EIS_Pars2
        PostDestDir = FindExpFolder("VERSASTAT").PostDir
        SeriesIDs = [
            SampleSelection.Series_CB_paper,
            SampleSelection.Series_Porhp_SiO2,
            SampleSelection.Series_Co_PANI,
            SampleSelection.Series_ML_SiO2,
            {"name": "all_EIS"},
        ]
        SeriesID_set = SeriesIDs[1]
        PDDORR = PostDestDir.joinpath("ORR_Pars_Char_{0}".format(SeriesID_set["name"]))
        PDDORR.mkdir(parents=True, exist_ok=True)

        ORR_series_pars_raw = ORR_pars.loc[
            ORR_pars.SampleID.isin(SeriesID_set["sIDslice"])
        ]
        if "all" in SeriesID_set["name"]:
            ORR_series_pars_raw = ORR_pars

        jos4 = ORR_pars.query(
            'SampleID == "JOS4" & pH == 1 & RPM > 900 & postAST == "no"'
        ).drop_duplicates(ExportfromORR.ORR_kin_cols)
        jos1 = ORR_pars.query(
            'SampleID == "JOS1" & pH == 1 & RPM > 900 & postAST == "no"'
        ).drop_duplicates(ExportfromORR.ORR_kin_cols)
        jos5 = ORR_pars.query(
            'SampleID == "JOS5" & pH == 1 & RPM > 900 & postAST == "no"'
        ).drop_duplicates(ExportfromORR.ORR_kin_cols)
        pta7 = postOVVout.query('basename == "O2_ORR_PTA7_0.1MH2SO4_#2_Ch1_Disk_PDX"')

        ORR_CB_paper_DW20 = ORR_pars.loc[
            ORR_pars.SampleID.isin(["DW20", "DW28", "DW29"])
        ].query(SampleSelection.acid1500)

    def plot_triangle_hm(
            grE,
            rcorr,
            Ev,
            target_dir,
            corr_method="pearson",
            corr_cutoff=0.5,
            plot_option=False,
        ):
            #            rcorr = dfcorr[corr_Cols].corr(method=corr_method)
            heatmap_title = "EIS_heatmap_{0}_{1}_{2}.png".format(
                corr_method, int(corr_cutoff * 100), Ev
            )
            print(heatmap_title)
            selcorr = rcorr[(rcorr > corr_cutoff) | (rcorr < corr_cutoff)]
            corr_triu = rcorr.where(np.tril(np.ones(rcorr.shape)).astype(np.bool))
            scorr_triu = corr_triu.stack()
            filter_corr = scorr_triu[
                (scorr_triu.abs() > corr_cutoff) & (scorr_triu.abs() < 0.97)
            ]
            fcorr = filter_corr.unstack()
            fcorr_sliced = fcorr.dropna(how="all")

            if not fcorr_sliced.empty and plot_option is "heatmap":
                plt.subplots(figsize=(40, 40))
                sns.heatmap(fcorr_sliced)
                #                heatmap_title = 'EIS_heatmap_{0}_{1}_{2}'.format(corr_method,int(corr_cutoff*100),EvRHE)
                plt.suptitle(heatmap_title)
                plt.savefig(
                    target_dir.joinpath(heatmap_title), dpi=300, bbox_inches="tight"
                )
                plt.grid(True)
                plt.close()
            else:
                fcorr_sliced.loc[
                    fcorr_sliced.index.isin(
                        SampleSelection.EC_EIS_par_cols + ["Qad+Cdlp"]
                    )
                ]
                drop_cols = ["Colorcode", "RedChisqr2", "RedChisqr1"]
                drop_cols_set = [i for i in drop_cols if i in fcorr_sliced.columns]
                promising = (
                    fcorr_sliced.loc[
                        fcorr_sliced.index.isin(
                            SampleSelection.EC_EIS_par_cols[0:-2] + ["Qad+Cdlp"]
                        )
                    ]
                    .drop(columns=drop_cols_set)
                    .unstack()
                    .dropna()
                    .sort_values()
                )
                #                if not promising.loc[[('Rct','Rct_kin'),('Rct_kin','Rct')]].dropna().empty:
                promising = promising.drop(
                    index=[
                        ("Rct", "Rct_kin"),
                        ("Rct_kin", "Rct"),
                        ("Cdlp", "Qad+Cdlp"),
                        ("Cdlp", "Qad+Cdlp"),
                        ("Qad", "Qad+Cdlp"),
                        ("Qad+Cdlp", "Qad"),
                    ],
                    errors="ignore",
                )
                toptail = pd.concat([promising.tail(5), promising.head(5)])
                PPDEIS_Ev = target_dir.joinpath(str(Ev))
                PPDEIS_Ev.mkdir(parents=True, exist_ok=True)
                def plot_pd_SampleIDs(*args):
                    print('fix me')# FIXME

                if plot_option == "corr":
                    for (xc, yc), corr_val in toptail.iteritems():
                        plot_pd_SampleIDs(grE, xc, yc, corr_val, PPDEIS_Ev)
                return promising

    def ORR_plotting_parameters2(plot_triangle_hm, EIS_O2_065_acid_no = pd.DataFrame()): # FIXME
        #%% ===== MAKE PLOTS of Parameters versus E v RHE ======
        PostDestDir = FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        PostDestDirEIScom = PostDestDir.joinpath("ORR_Pars_Char_CB_paper")
        PPDEIS = PostDestDirEIScom.joinpath("EIS_corrs")
        PPDORR = PostDestDirEIScom.joinpath("ORR_corr")
        PPDORR.mkdir(parents=True, exist_ok=True)
        SampleSelection.EC_exp_cols
        #        EIS_O2_065_acid_no.drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)

        corr_method, corr_cutoff = "pearson", 0.5  # or default is pearson spearman
        rcorr = EIS_O2_065_acid_no[
            SampleSelection.InterestingCols
            + SampleSelection.EC_EIS_par_cols
            + SampleSelection.RAMAN_cols_corr
            + ["Qad+Cdlp"]
        ].corr()
        plt.subplots(figsize=(20, 15))
        selcorr = rcorr[(rcorr > corr_cutoff) | (rcorr < corr_cutoff)]
        sns.heatmap(selcorr)
        plt.savefig(
            PostDestDirEIScom.joinpath(
                "EIS_heatmap_{0}_{1}.png".format(corr_method, corr_cutoff)
            ),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
        #        melt= pd.melt(rcorr.reset_index(), id_vars='index')
        #        melt.columns
        ##            FileHelper.PlotAnalysis.corrplot(rcorr)
        #        ===== Making plots of correlations of EIS Parameters versus E v RHE ======
        Gases, pHses = ["N2", "O2"], [0.3, 1, 13]
        #        EIS_O2_no = EIS_CB_paper.query('(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))').drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
        #                pH03, pH1 = EIS_O2_acid_no.query('(pH == 0.3)'), EIS_O2_acid_no.query('(pH == 1)')

        corr_Cols = (
            SampleSelection.InterestingCols
            + SampleSelection.EC_EIS_par_cols
            + SampleSelection.RAMAN_cols_corr
            + ["Qad+Cdlp"]
        )
        corr_method, corr_cutoff = "pearson", 0.5  # or default is pearson spearman
        #        rcorr = EIS_O2_065_acid_no[corr_Cols].corr(method=corr_method)
        EIS_CB_paper = pd.DataFrame() # FIXME


        out_topcorrs_lst = []
        for Gas_set in Gases:
            for pH_set in pHses:
                EIS_O2_no_query = EIS_CB_paper.query(
                    '(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
                ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
                ORREIS_O2_no_query = EIS_CB_paper.query(
                    '(Gas == @Gas_set) & (pH == @pH_set)& (postAST == "no") & ((Loading_cm2 < 0.5) & (Loading_cm2 > 0.3))'
                ).drop_duplicates(subset=SampleSelection.EC_EIS_par_cols)
                target_dir = PPDEIS.joinpath(
                    "EIS_corr_{0}_pH{1}".format(Gas_set, pH_set)
                )
                target_dir.mkdir(parents=True, exist_ok=True)
                for EvRHE, grE in EIS_O2_no_query.groupby("E_RHE"):
                    if len(grE) > 2:
                        EvRHE, grE
                        rcorr = grE[corr_Cols].corr(method=corr_method)
                        prom_corr = plot_triangle_hm(
                            grE,
                            rcorr,
                            np.round(EvRHE, 2),
                            target_dir,
                            plot_option=False,
                        )
                        out_topcorrs_lst.append([Gas_set, pH_set, EvRHE, prom_corr])

        pd.DataFrame(out_topcorrs_lst[0][-1], columns=[corr_method]).assign(
            **{
                "Gas": out_topcorrs_lst[0][0],
                "pH": out_topcorrs_lst[0][1],
                "E_RHE": out_topcorrs_lst[0][2],
            }
        )
        topcorrs = pd.concat(
            [
                pd.DataFrame(i[-1], columns=[corr_method]).assign(
                    **{"Gas": i[0], "pH": i[1], "E_RHE": i[2]}
                )
                for i in out_topcorrs_lst
            ]
        )

        topcorrs.groupby(level=[0, 1]).size()
        topcorrs.groupby(level=[0, 1]).sum()
        topcorrs.loc[np.abs(topcorrs[corr_method] > 0.5)].groupby(
            level=[0, 1]
        ).sum().sort_values(by=corr_method).iloc[-10::].plot.barh(x=corr_method)
        topcorrs[corr_method].groupby(level=[0, 1]).sum().sort_values(
            by=corr_method
        ).iloc[-10::].plot.barh(x=corr_method)
        top_score_sum = topcorrs[corr_method].groupby(level=[0, 1]).sum().sort_values()
        top_score_sum.plot.barh()
        top_score_sum[(12 < top_score_sum) | (-11 > top_score_sum)].plot.barh(
            figsize=(10, 16)
        )
        top_score_abssum = (
            np.abs(topcorrs[corr_method]).groupby(level=[0, 1]).sum().sort_values()
        )
        abs_best = top_score_abssum.loc[
            top_score_abssum > top_score_abssum.mean() + 1 * top_score_abssum.std()
        ]
        abs_best.plot.barh(figsize=(12, 12))

        topcorrs.query('(pH < 5) & (Gas == "O2")').index
        sumbest_acid_E = (
            topcorrs.query('(pH < 5) & (Gas == "O2")')[[corr_method, "E_RHE"]]
            .reset_index()
            .set_index(["level_0", "level_1", "E_RHE"])
            .groupby(level=[0, 1, 2])
            .sum()
        )
        sumbest_acid = (
            sumbest_acid_E[np.abs(sumbest_acid_E) > 0.7].dropna().reset_index("E_RHE")
        )
        sumbest_acid.plot(x="E_RHE", y=corr_method, kind="scatter")
        sumbest_acid.groupby(level=[0, 1])

        #                rcorr = gr_sID_desc_Codes[SampleSelection.InterestingCols+SampleSelection.RAMAN_cols_corr].corr(method=corr_method)

