# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 11:33:55 2020

@author: User
"""
import sys
from pathlib import Path
import functools

# import post_helper
# import plotting
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd
import numpy as np

import types
import matplotlib as mpl

mpl.rcParams["figure.dpi"] = 100
# from sklearn.cluster import KMeans


# print ('Name prepare input:', __name__ )
if __name__ == "__main__":
    #    print(f'Package: {__package__}, File: {__file__}')
    #    FH_path = Path(__file__).parent.parent.parent.joinpath('FileHelper')
    #    sys.path.append(str(FH_path))
    #    sys.path.append(str(Path(__file__).parent.parent.joinpath('indexer')))
    sys.path.append(str(Path(__file__).parent.parent.parent))
    #    sys.path.append("..")
    #    print(sys.path)
    #    import FileHelper
    from FileHelper.PostChar import Characterization_TypeSetting, SampleCodesChar
    from FileHelper.PostPlotting import *
    from FileHelper.FindSampleID import GetSampleID
    from FileHelper.FindFolders import FindExpFolder

    #    from FileHelper.FindExpFolder import FindExpFolder
    from plotting import eisplot
elif "prepare_input" in __name__:
    import RunEC_classifier
    from FileHelper.FindSampleID import FindSampleID

# from FileHelper.PostChar import SampleSelection, Characterization_TypeSetting
def mkfolder(folder):
    folder.mkdir(exist_ok=True, parents=True)
    return folder


OriginColors = Characterization_TypeSetting.OriginColorList()
Pfolder = FindExpFolder().TopDir.joinpath(
    Path("Preparation-Thesis/SiO2_projects/SiO2_Me_EC+Struc")
)
plotsfolder = mkfolder(Pfolder.joinpath("correlation_plots"))


def get_float_cols(df):
    return [key for key, val in df.dtypes.to_dict().items() if "float64" in str(val)]


# class PorphSamples():
#    def __init__(self):
#        self.template = PorphSamples.template()
def decorator(func):
    @functools.wraps(func)
    def wrapper_decorator(*args, **kwargs):
        # Do something before
        value = func(*args, **kwargs)
        # Do something after
        return value

    return wrapper_decorator


def PorphSiO2_template():
    #        'SerieIDs' : ('Porph_SiO2')*5,
    Series_Porph_SiO2 = {
        "SampleID": ("JOS1", "JOS2", "JOS3", "JOS4", "JOS5"),
        "Metal": ("Fe", "Co", "MnTPP", "FeTPP", "H2"),
        "color": (2, 4, 6, 15, 3),
    }

    Porphyrins = {
        "TMPP": {"Formula": "C48H38N4O4", "MW": 734.8382},
        "TMPP-Fe(III)Cl": {"Formula": "C48H36ClFeN4O4", "MW": 824.1204},
        "TMPP-Co(II)": {"Formula": "C48H36CoN4O4", "MW": 791.7556},
        "TTP-Mn(III)Cl": {"Formula": "C44H28ClMnN4", "MW": 703.1098},
        "TPP-Fe(III)Cl": {"Formula": "C44H28ClFeN4", "MW": 704.0168},
        "TPP": {"Formula": "C44H30N4", "MW": 614.7346},
    }

    Porph_template = pd.DataFrame(Series_Porph_SiO2)
    return Porph_template


class EC_PorphSiO2:

    folder = FindExpFolder("PorphSiO2").compare
    Porph_template = PorphSiO2_template()
    EIS_models = [
        "Model(EEC_Randles_RWpCPE)",
        "Model(EEC_2CPE)",
        "Model(EEC_2CPEpW)",
        "Model(EEC_RQ_RQ_RW)",
        "Model(EEC_RQ_RQ_RQ)",
        "Model(Randles_RQRQ)",
    ]
    #    ['Model(Singh2015_RQRQ)', 'Model(Singh2015_RQRQR)', 'Model(Bandarenka_2011_RQRQR)',
    #                  'Model(Singh2015_RQRWR)', 'Model(Randles_RQRQ)', 'Model(Singh2015_R3RQ)']
    #    model_select = EC_PorphSiO2.EIS_models[1]

    def __init__(self):
        self.pars = EC_PorphSiO2.mergedEC()

    #        self.par_export = EC_OHC.to_excel(self.folder.joinpath('EC_ORR_HPRR.xlsx'))

    def mergedEC():
        template = PorphSiO2_template()
        HPRR, N2CV = EC_PorphSiO2.HPRR(), EC_PorphSiO2.N2cv()
        ORR = EC_PorphSiO2.ORR_updated_pars()
        EIS = EC_PorphSiO2.EIS_pars()
        mcols = ["SampleID", "Sweep_Type"]
        ECmerged = pd.merge(ORR, pd.merge(N2CV, HPRR, on=mcols, how="left"), on=mcols)
        EC_EIS = pd.merge(ECmerged, EIS, on=mcols)

        EC_OHN_merged = pd.merge(template, EC_EIS, on="SampleID")
        EC_PorphSiO2.export_to_xls(EC_OHN_merged)
        return EC_OHN_merged

    #        EC_all_merged_lst.append(EC_OHN_merged)
    #        EC_all_merged = pd.concat(EC_all_merged_lst)
    #        ORR_cath = EC_PorphSiO2.ORR_updated_pars(sweep_type_select='cathodic')
    #        ORR_an = EC_PorphSiO2.ORR_updated_pars(sweep_type_select='anodic')
    #        EC_OHN2 = pd.merge(template, pd.merge(ORR_an,pd.merge(HPRR, N2CV),on='SampleID'), on='SampleID')
    #        EC_OHN2_cath = pd.merge(template, pd.merge(ORR,pd.merge(HPRR, N2CV),on='SampleID'), on='SampleID')
    #        EC_OHN2.to_excel(FindExpFolder('PorphSiO2').compare.joinpath('EC_ORR_HPRR_N2.xlsx'))
    def export_to_xls(EC_OHN_merged):
        export_path = FindExpFolder("PorphSiO2").compare.joinpath(f"EC_pars_all.xlsx")
        if "Sweep_Type" in EC_OHN_merged.columns:
            with pd.ExcelWriter(export_path) as writer:
                for swp, swpgr in EC_OHN_merged.groupby("Sweep_Type"):
                    swpgr.to_excel(writer, sheet_name=swp)
                    swpgr.to_excel(export_path.with_name(f"EC_pars_{swp}.xlsx"))
        else:
            export_path = FindExpFolder("PorphSiO2").compare.joinpath(
                "EC_pars_no-sweep.xlsx"
            )
            EC_OHN_merged.to_excel(export_path)
        print(f"EC pars saved to:\n{export_path}")
        return export_path

    def edit_columns(func, template=PorphSiO2_template()):
        def wrapper(*args, **kwargs):
            if kwargs:
                pars_out, suffx = func(**kwargs)
            else:
                pars_out, suffx = func()
            cols = [
                i
                for i in pars_out.columns
                if i not in template.columns and i != "Sweep_Type"
            ]
            pars_out = pars_out.rename(
                columns=dict(zip(cols, [f"{suffx}_" + i for i in cols]))
            )
            return pars_out

        return wrapper

    @edit_columns
    def HPRR(sweep_type_select=["anodic", "cathodic"]):
        hfs = []
        for swp in sweep_type_select:
            hprr_files = list(EC_PorphSiO2.folder.rglob(f"*{swp}*HPRR*disk*"))
            #            print(hprr_files)
            for hf in hprr_files:
                hprr_raw = pd.read_excel(hf)
                hprr_raw["file"] = hf.stem
                E_App_col = [i for i in hprr_raw.columns if "E_APP" in i.upper()][0]
                E_jmin = hprr_raw.iloc[np.abs(hprr_raw["jmAcm-2"]).idxmin()][E_App_col]
                sID = GetSampleID.try_find_sampleID(hf)[0]
                fit_lin_fit = linregress(hprr_raw[E_App_col], hprr_raw["HPRR_j0_Fit"])

                hfs.append(
                    {
                        "SampleID": sID,
                        "E_onset": E_jmin,
                        "dj/dE": fit_lin_fit[0],
                        "Sweep_Type": swp,
                    }
                )
        HPRR_pars_origin = pd.DataFrame(hfs)
        return HPRR_pars_origin, "HPRR"

    @edit_columns
    def ORR():
        orr_files, orrfs = list(EC_PorphSiO2.folder.rglob("*ORR_pars*")), []
        for of in orr_files:
            orr_raw = pd.read_excel(of, index_col=[0])
            orr_raw.query("RPM > 1400")
            orrfs.append(orr_raw.query("RPM > 1400"))
        ORR_pars_origin = pd.concat(orrfs, ignore_index=True).reset_index(drop=True)
        return ORR_pars_origin, "ORR"

    @edit_columns
    def ORR_updated_pars(sweep_type_select=["anodic", "cathodic"]):
        orr_files, orrfs = list(EC_PorphSiO2.folder.rglob("PostEC*ORR*Jkin*pars*")), []
        for of in orr_files:
            orr_raw = pd.read_excel(of, index_col=[0])
            orr_raw.query("RPM > 1400")
            orrfs.append(orr_raw.query("RPM > 1400"))
        ORR_pars_origin = pd.concat(orrfs, ignore_index=True).reset_index(drop=True)
        sweep_col = [
            i
            for i in ORR_pars_origin.columns
            if "SWEEP" in i.upper() and not "RING" in i.upper()
        ][0]
        ORR_pars_origin_swp = ORR_pars_origin.loc[
            ORR_pars_origin[sweep_col].isin(sweep_type_select)
        ]
        ORR_pars_origin_swp = ORR_pars_origin_swp.rename(
            columns={sweep_col: "Sweep_Type"}
        )
        #        ORR_pars_swp  = {n : gr for n,gr in ORR_pars.groupby('Sweep_Type_disk')}
        #        ORR_anod, ORR_cath = ORR_pars_swp.get('anodic'), ORR_pars_swp.get('cathodic')
        return ORR_pars_origin_swp, "ORR"

    @edit_columns
    def N2cv(sweep_type_select=["anodic", "cathodic"], unit="F"):
        N2_files, N2fs = list(EC_PorphSiO2.folder.rglob("*CVs*xlsx")), []
        if unit == "mF":
            unit_factor = 1
        elif unit == "F":
            unit_factor = 1e-3
        else:
            unit_factor = 1
        for n2f in N2_files:
            n2_raw = pd.read_excel(n2f, index_col=[0], sheet_name=None)
            sID = GetSampleID.try_find_sampleID(n2f)[0]
            for swp in sweep_type_select:
                anod = n2_raw.get(swp)
                evlst = {"SampleID": sID, "Sweep_Type": swp}
                for Ev in np.arange(0, 1000, 100):
                    anod_05 = anod.loc[np.isclose(anod.index, Ev * 1e-3, atol=0.001)]
                    mean_anod05 = np.abs(anod_05.mean(axis=0)["Cdl_mFcm-2"])
                    if not np.isnan(mean_anod05):
                        #                        N2fs.append({'SampleID' :  sID, f'Cdl_mFcm-2_{Ev}' : mean_anod05, 'Sweep_Type' : swp})
                        evlst.update(
                            {f"Cdl_{unit}cm-2_{Ev}": mean_anod05 * unit_factor}
                        )
                N2fs.append(pd.DataFrame(evlst, index=[(sID, swp)]))
        N2_orig = pd.concat(N2fs)
        #        N2_orig = pd.DataFrame(N2fs) #.set_index('SampleID','Sweep_Type')
        return N2_orig, "N2"

    @edit_columns
    def EIS_pars(model_select="Model(EEC_2CPE)"):
        #        sweep_type_select = ['anodic','cathodic'], model_select = 'Model(Singh2015_RQRQR)'
        eis_files, eisfs = (
            list(
                EC_PorphSiO2.folder.parent.joinpath(
                    f"EIS_Porph_SiO2\{model_select}"
                ).rglob("JOS*.xlsx")
            ),
            [],
        )
        if eis_files:
            for ef in eis_files:
                eis_raw = pd.read_excel(ef, index_col=[0])
                eisfs.append(eis_raw)
            EIS_pars_origin = pd.concat(eisfs, ignore_index=True).reset_index(drop=True)
        else:
            print("EIS pars file list empty!!")

        EPgrp = EIS_pars_origin.groupby("Gas")
        EP_N2, EP_O2 = EPgrp.get_group("N2").drop(columns="Gas"), EPgrp.get_group(
            "O2"
        ).drop(columns="Gas")
        EIS_N2O2 = pd.merge(EP_N2, EP_O2, suffixes=["_N2", "_O2"], on="SampleID")
        redchis = [i for i in eis_raw.columns if "RedChi" in i]
        Rsis = [i for i in EIS_N2O2.columns if "Rs" in i and not "stderr" in i]
        Rct_cols = [
            i
            for i in EIS_N2O2.columns
            if "Rct" in i and not any(c in i for c in ("stderr", "_kin_"))
        ]
        #        EIS_pars_origin[Rsis] = EIS_pars_origin[Rsis].mask(EIS_pars_origin[Rsis] < 1)
        EIS_N2O2[Rsis] = EIS_N2O2[Rsis].mask(EIS_N2O2[Rsis] < 1)
        print("EIS Rs mask applied")
        EIS_N2O2[Rct_cols] = EIS_N2O2[Rct_cols].mask(EIS_N2O2[Rct_cols] > 1e5)
        print("EIS Rct mask applied")
        EIS_N2O2 = EIS_N2O2.dropna(axis=1, how="all")
        #        RedChiSq_limit = ORReis_merge.query('Rs > 1').RedChisqr.mean()+ 1*ORReis_merge.query('Rs > 1').RedChisqr.std()
        #        ORReis_neat = ORReis_merge.query('RedChisqr < @RedChiSq_limit & Rs > 2 & Rct < 9E05')
        EIS_N2O2_an, EIS_N2O2_cat = EIS_N2O2.copy(), EIS_N2O2.copy()
        EIS_N2O2_an["Sweep_Type"] = "anodic"
        EIS_N2O2_cat["Sweep_Type"] = "cathodic"
        EIS_N2O2_new = pd.concat([EIS_N2O2_an, EIS_N2O2_cat], axis=0)
        #        EIS_pars_orig_mod = EIS_pars_origin.query('Model_EEC == @model_select')
        return EIS_N2O2_new, "EIS"

    def EIS_spectra_origin_prep(model_select="Model(EEC_2CPE)"):
        eis_metaf, _specs = (
            list(
                EC_PorphSiO2.folder.parent.rglob(
                    "EIS_Porph_SiO2\meta_data*EIS*origin.xlsx"
                )
            ),
            [],
        )
        EISmeta = pd.read_excel(eis_metaf[0], index_col=[0])
        EISmeta.columns
        for (sID, gas), pgrp in EISmeta.groupby(["SampleID", "Gas"]):  # 'PAR_file'
            PF, pgrp
            EIScombined = pd.read_excel(pgrp.SpectraFile.iloc[0], index_col=[0])
            EISspectra_mod = EIScombined.query("Model_EEC == @model_select")
            EISspectra_mod = make_uniform_EvRHE(EISspectra_mod)
            for Ev, Egrp in EISspectra_mod.groupby("E_RHE"):
                Egrp = Egrp.assign(**{"SampleID": sID, "Gas": gas})
                _specs.append(Egrp)
        #                _specs.update({(sID,gas,Ev) : Egrp})
        spectra = pd.concat(_specs)
        spectra.to_excel(eis_metaf[0].with_name("clean_spectra.xlsx"))

    def EIS_spectra_origin(model_select="Model(EEC_2CPE)"):
        eis_metaf, _specs = (
            list(
                EC_PorphSiO2.folder.parent.rglob(
                    f"EIS_Porph_SiO2\{model_select}\meta_data*EIS*origin.xlsx"
                )
            ),
            [],
        )

        specdir = mkfolder(eis_metaf[0].parent.joinpath("spectra"))
        spectra = pd.read_excel(
            eis_metaf[0].with_name("clean_spectra.xlsx"), index_col=[0]
        )
        spectra.columns
        for ax_type in [("Zre", "-Zim"), ("Yre", "Yim")]:
            cols = [i + a for i in ["DATA_", "FIT_"] for a in ax_type]
            for gas, Ggrp in spectra.groupby("Gas"):
                for sID, sgrp in Ggrp.groupby("SampleID"):
                    with pd.ExcelWriter(
                        specdir.joinpath(f"{ax_type[0][0]}_{gas}_{sID}.xlsx")
                    ) as writer:
                        for Ev, Egrp in sgrp.groupby("E_RHE"):
                            # sheet_name = Ev
                            EmV = f"{1E3*Ev:.0f}"
                            Egrp[["Frequency(Hz)"] + cols].to_excel(
                                writer, sheet_name=EmV
                            )
                            # === plotting
                            fig, ax = plt.subplots()
                            Egrp.plot(
                                x=cols[0],
                                y=cols[1],
                                kind="scatter",
                                ax=ax,
                                label=cols[1],
                            )
                            Egrp.plot(x=cols[2], y=cols[3], c="r", ax=ax, label=cols[3])
                            plt.legend()
                            ax.set_xlabel(ax_type[0])
                            ax.set_ylabel(ax_type[1])
                            ax.set_title(f"{gas} {sID} {EmV}")
                            ax.grid(True)
                            plt.savefig(
                                specdir.joinpath(f"{ax_type[0][0]}_{gas}_{sID}_{EmV}"),
                                bbox_inches="tight",
                            )
                            plt.close()
                            # ===

    def corr_plots():
        EC_OHC.query('SampleID != "JOS5"').corr()
        corrstk = EC_OHC.query('SampleID != "JOS5"').corr().stack()

        EC_OHC.plot(x="E_onset", y="HPRR_E_onset", kind="scatter")
        EC_OHC.plot(x="FracH2O2_050", y="HPRR_E_onset", kind="scatter")
        EC_OHC.plot(x="N2_Cdl_mFcm-2_0.5", y="HPRR_dj/dE", kind="scatter")
        EC_OHC.plot(x="N2_Cdl_mFcm-2_0.5", y="E_half", kind="scatter")

        EC_OHC.corr(method="pearson")


def EC_PorphSio():
    #    folder = Path('F:\EKTS_CloudStation\CloudStation\Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc\EC_Porph_SiO2_0.1MH2SO4\Compare_parameters')
    #    folder = Path('G:\CloudStation\Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc\EC_Porph_SiO2_0.1MH2SO4\Compare_parameters')
    #    HPRR = pd.concat([pd.read_excel(i)['file'] for i in hprr_files])
    EC_ORR_HPRR = pd.merge(ORR_pars_origin, HPRR_pars_origin)
    HPRR_pars_origin.join(N2_orig, on="SampleID")
    EC_OHC = pd.merge(
        ORR_pars_origin, pd.merge(HPRR_pars_origin, N2_orig), on="SampleID"
    )
    #        orr_raw.query('RPM > 1400')
    #        orrfs.append(orr_raw.query('RPM > 1400'))
    EC_OHC.to_excel(folder.joinpath("EC_ORR_HPRR.xlsx"))


class postChar:
    def __init__(self):
        self.folder = FindExpFolder("PorphSiO2").folder
        #        FindExpFolder().TopDir.joinpath(Path('Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc'))
        self.raman = postChar.postRaman
        self.bet = postChar.postBET
        self.ea = postChar.postEA
        self.prep = postChar.prec_ea_ratios(
            PorphSiO2_template(), self.ea, postChar.postPrep
        )
        self.merged = postChar.merged(self)

    def merged(self):
        cols = list(PorphSiO2_template().columns)
        merged_out = pd.merge(
            self.prep,
            pd.merge(self.ea, pd.merge(self.raman, self.bet, on=cols), on=cols),
            on=cols,
        )
        return merged_out

    def decorator(func):
        @functools.wraps(func)
        def wrapper_decorator(*args, **kwargs):
            # Do something before
            value = postChar.slice_out(*args, slice_lst=postChar.template().SampleID)
            # Do something after
            return value

        return wrapper_decorator

    def slice_out(func, template=PorphSiO2_template()):
        pars, suffx = func()
        try:
            slice_lst = template.SampleID
            pars_out = pars.loc[pars.SampleID.isin(slice_lst)]
            pars_out = pd.merge(template, pars_out, on="SampleID")
            if suffx != "":
                cols = [i for i in pars_out.columns if i not in template.columns]
                pars_out = pars_out.rename(
                    columns=dict(zip(cols, [f"{suffx}_" + i for i in cols]))
                )
        except:
            pars_out = pars
            print("No slice out for {func.__name__} ")
        return pars_out

    @slice_out
    def postRaman(peak_model="5peaks", plot_resid=False):
        ramandir = FindExpFolder("PorphSiO2").folder / "SiO2_Me_EC+Struc" / "Raman"
        raman_files = ramandir.rglob("*xlsx")
        fitpars_fls = [
            i
            for i in raman_files
            if all([a in i.stem for a in ["FitParameters_Model", "peaks"]])
        ]
        FitPars_raw = pd.concat(
            [
                pd.read_excel(i).assign(**{"Model": i.stem.split("Model_")[-1]})
                for i in fitpars_fls
            ],
            sort=False,
        )
        FitPars = FitPars_raw.loc[
            FitPars_raw.SampleID.isin(PorphSiO2_template().SampleID.values)
        ]

        if plot_resid:
            Model_redchi = pd.DataFrame(
                [
                    (n, np.abs(gr.redchi).mean(), np.abs(gr.redchi).sum())
                    for n, gr in FitPars.groupby("Model")
                ],
                columns=["Model", "redchi_mean", "redchi_sum"],
            ).set_index("Model")
            Model_chi = pd.DataFrame(
                [
                    (n, gr.chisqr.mean(), gr.chisqr.sum())
                    for n, gr in FitPars.groupby("Model")
                ],
                columns=["Model", "redchi_mean", "redchi_sum"],
            ).set_index("Model")
            Model_redchi.plot.bar()
            Model_chi.plot.bar()
        if peak_model:
            FPars_out_1st = FitPars.loc[FitPars.Model.isin([peak_model])]
        else:
            FPars_out_1st = FitPars

        t2nd_mod = "2ndOrder_4peaks"
        if t2nd_mod in FitPars_raw.Model.unique():
            FPars_out_2nd = FitPars.loc[FitPars.Model == t2nd_mod].dropna(axis=1)
            flt2nd = get_float_cols(FPars_out_2nd)
            FPars_out_2nd = FPars_out_2nd.rename(
                columns=dict(zip(flt2nd, [f"2nd_{i}" for i in flt2nd]))
            )
            FPars_out = pd.merge(
                FPars_out_1st.dropna(axis=1),
                FPars_out_2nd,
                on="SampleID",
                how="left",
                suffixes=["_1st", "_2nd"],
            )

        else:
            FPars_out = FPars_out_1st
        return FPars_out, "R"

    @slice_out
    def postBET():
        betdir = FindExpFolder("PorphSiO2").folder.joinpath("SiO2_Me_EC+Struc/BET")
        bet_files = betdir.rglob("Overview_*xlsx")
        BET_ovv = pd.concat(
            [pd.read_excel(i, index_col=[0], sort=False) for i in bet_files]
        )
        return BET_ovv, "B"

    @slice_out
    def postEA():
        EAcols = [
            "SampleID",
            "C/N_ratio",
            "N_content",
            "C_content",
            "H_content",
            "100-CHN",
        ]
        EA_results = pd.read_excel(list(FindExpFolder("EA").DestDir.rglob("*xlsx"))[0])[
            EAcols
        ]
        EA_results["C/H_ratio"] = EA_results["C_content"] / EA_results["H_content"]
        return EA_results, "EA"

    @slice_out
    def postPrep():
        Prep_Porph_SiO2 = {
            "SampleID": ("JOS1", "JOS2", "JOS3", "JOS4", "JOS5"),
            "WL_precursor": (31.0, 56.0, 53.0, 38.0, 41.0),
            "Precursor_type": (
                "FeTMPPCl",
                "CoTMPPCl",
                "MnTPPCl",
                "FeTPPCl",
                "H2TMPPCl",
            ),
            "MW_g-mol": (824.12, 791.67, 703.11, 704.02, 734.84),
            "Metal_element": (26.0, 27.0, 25.0, 26.0, 1.0),
        }
        Prep_Porph = pd.DataFrame(Prep_Porph_SiO2)
        Prep_Porph["prec_N"] = 100 * (4 * 14) / Prep_Porph["MW_g-mol"]
        Prep_Porph["prec_C"] = (
            100 * (12 * np.array([48, 48, 44, 44, 48])) / Prep_Porph["MW_g-mol"]
        )
        Prep_Porph["prec_H"] = (
            100 * (1 * np.array([36, 36, 28, 28, 38])) / Prep_Porph["MW_g-mol"]
        )
        Prep_Porph["prec_100-CHN"] = 100 - (
            Prep_Porph["prec_H"] + Prep_Porph["prec_N"] + Prep_Porph["prec_C"]
        )
        Prep_Porph["prec_C/N_ratio"] = Prep_Porph["prec_C"] / Prep_Porph["prec_N"]

        return Prep_Porph, "P"

    def prec_ea_ratios(templ, ea, prep):
        pcols = prep.columns
        _ms = []
        for eacol in [i for i in ea.columns if i not in templ.columns]:
            if eacol.endswith("_content"):
                _mcol = [
                    i for i in pcols if i.split("P_prec_")[-1] in eacol.split("EA_")[-1]
                ]
            else:
                _mcol = [
                    i for i in pcols if i.split("P_prec_")[-1] == eacol.split("EA_")[-1]
                ]

            if _mcol:
                _ms.append((eacol, _mcol[0]))
        _coldict = {}
        for mc in _ms:
            elem = mc[1].split("P_prec_")[-1]
            ratiocol = f"P_ratioEA_{elem}"
            prep = prep.assign(**{ratiocol: prep[mc[1]] / ea[mc[0]]})
            _coldict.update(
                {elem: {"ea_col": mc[0], "prep_col": mc[1], "ratio_col": ratiocol}}
            )
        print("Added Ratio cols to prep")
        return prep


#        EA_results.loc[EA_results.SampleID.isin(isin_lst)]
class postCorrs:
    def __init__(self):
        pass

    def prep_ea(prep, ea, plots=False):
        prep_ea = pd.merge(prep, ea)
        corr_prep_ea = allchars_corr(prep_ea)
        if plots == True:
            testfolder = mkfolder(plotsfolder.joinpath("prep_ea_test"))
            for i in corr_prep_ea.index.values:
                plot(prep_ea, x=i[0], y=i[1], savefig=testfolder)
        return corr_prep_ea


def collect_pars():
    ECpars = EC_PorphSiO2().pars
    pst = postChar()
    raman, bet, ea, prep = pst.raman, pst.bet, pst.ea, pst.prep
    struc = pst.merged
    Allchars = pd.merge(ECpars, struc)

    char = {"EC": ECpars, "raman": raman, "bet": bet, "ea": ea, "prep": prep}
    return ECpars, raman, bet, ea, prep, struc, char, Allchars


def export_chars(char, suffix=""):
    for nm, df in char.items():
        swpcol = [i for i in df.columns if "Sweep" in i]
        if swpcol:
            for swp, swgr in df.groupby(swpcol[0]):
                swgr.to_excel(plotsfolder.joinpath(f"char_{nm}_{swp}{suffix}.xlsx"))
        else:
            df.to_excel(plotsfolder.joinpath(f"char_{nm}{suffix}.xlsx"))


def make_corrs(pars, charname):
    lst = []
    #    ECpars[[key for key,val in ECpars.dtypes.to_dict().items() if 'float64' in str(val)]]
    only_floats = [
        key for key, val in pars.dtypes.to_dict().items() if "float64" in str(val)
    ]
    for n, r in pars.iterrows():
        without1 = pars.drop(n)[only_floats]
        wcorrstk = without1.corr(method="pearson").stack()
        lst.append(wcorrstk.loc[wcorrstk < 1].rename(f"{charname}_{r.SampleID}"))
    corrs_outlier = pd.concat(lst, axis=1)
    abs_outliers = (
        np.abs(corrs_outlier).describe().sort_values("mean", axis=1, ascending=False)
    )
    abs_print = ",".join(
        [
            f"{key} : {val:.2f}"
            for key, val in (abs_outliers.loc["mean"]).to_dict().items()
        ]
    )
    print("Correlation without samples in collumns, indicates outliers...\n", abs_print)
    return {"corrs_outliers": corrs_outlier, "describe": abs_outliers}


def corrloop(char):
    out = []
    for charname, pars in char.items():
        out.append(make_corrs(pars, charname))
    return out


def corrplot():
    corr = PlotAnalysis.corrplot(ECpars.corr())
    corr = corr.loc[
        ~(corr.y.isin(template.columns) | corr.x.isin(template.columns))
        & (corr.value < 1)
    ]
    tops = corr.loc[corr.y.str.contains("HPRR")].sort_values(
        by="value", ascending=False
    )
    for n, r in tops.tail(10).iterrows():
        plot(ECpars, x=r.x, y=r.y)


def plot_heatmap():
    fig, ax = plt.subplots(figsize=(20, 20))
    ax = sns.heatmap(allchflrtd[allcorrs.name].unstack())
    plt.savefig(FindExpFolder("PorphSiO2").compare.joinpath("heatmap.png"))
    plt.close()


def allchars_corr(Allchars):
    AllCorrs = {}
    if "Sweep_Type" not in Allchars.columns:
        Allchars["Sweep_Type"] = "none"

    for swp, swgrp in Allchars.groupby("Sweep_Type"):
        corvalnm = f"corr_val_{swp}"
        only_floats = [
            key for key, val in swgrp.dtypes.to_dict().items() if "float64" in str(val)
        ]
        filter_cols_lst = [
            "error",
            "chi",
            "Coefficient",
            "rvalue",
            "pvalue",
            "Isotherm_P/P0",
            "relPrange",
            "wrong",
            "B_DFT",
            "AMBIENT ",
            "R_bic",
            "R_aic",
            "BET_P/Po_",
            "BET_1/(W((Po/P)-1))",
            "B_Area_in_cell_m2",
        ]
        filter_cols_lst += [
            i for i in only_floats if i.startswith("B_") and "_RPT" in i
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("B_")
            and any(i.endswith(c) for c in ["_microArea", "_Vmicro", "_r", "mass_RAW"])
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("B_t") and any(i.endswith(c) for c in ["_slope"])
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("R_")
            and any(i.endswith(c) for c in ["_amplitude", "_sigma", "_bic"])
        ]
        filter_cols_lst += [
            i
            for i in only_floats
            if i.startswith("EIS_")
            and any(
                c in i.split("_")
                for c in ["Chisqr", "stderr", "RedChisqr", "AppV", "aic", "bic"]
            )
        ]
        #        filter_cols_lst += [i for i in only_floats if i.startswith('B_') and i.endswith()]
        #        filter_cols_lst += [i for i in only_floats if i.startswith('B_') and i.endswith()]
        filtred_cols = [
            i for i in only_floats if not any(t in i for t in filter_cols_lst)
        ]
        allcorrstk = swgrp[filtred_cols].corr(method="pearson").stack()
        allcorrs = allcorrstk[np.abs(allcorrstk) < 1]
        allcorrs.name = corvalnm
        allchrrs = pd.concat(
            [
                allcorrs,
                pd.Series(
                    [(i[0].split("_")[0], i[1].split("_")[0]) for i in allcorrs.index],
                    name="chars",
                    index=allcorrs.index,
                ),
            ],
            axis=1,
        )
        allchrrs = (
            allchrrs[np.abs(allchrrs[allcorrs.name]) < 0.991]
            .drop_duplicates(allcorrs.name)
            .sort_values(allcorrs.name)
        )

        allchflrtd = remove_corrdups(allchrrs)
        AllCorrs.update({f"{swp}": allchflrtd})

    if "none" in Allchars["Sweep_Type"].unique()[0]:
        return AllCorrs.get("none").sort_values(by="corr_val_none")
    else:
        return AllCorrs


def uniq(AllCan, combo_select, corvalnm, self_corr=False, Ev: int = 0):
    uniqchars = AllCan["chars"].unique()

    if "*" in combo_select:
        uniq_char = [i for i in uniqchars if any(c in i for c in combo_select)]
    else:
        uniq_char = [i for i in uniqchars if set(i) == set(combo_select)]
    if self_corr == True:
        uniq_char += [i for i in uniqchars if (combo_select[0], combo_select[0]) == i]
    chC = AllCan.loc[AllCan["chars"].isin(uniq_char)].sort_values(corvalnm)
    if Ev > 50:
        AllC_Ev = get_corrs_E(chC, Ev).sort_values(corvalnm)
        return AllC_Ev
    else:
        return chC


def select(chC, one_col: str, exclude="", endswith_set=""):
    if not any(one_col in i for i in ["", "*"]):
        idx_select = [i for i in chC.index.values if any(one_col in a for a in i)]
        if endswith_set:
            idx_select = [
                i for i in idx_select if any(a.endswith(endswith_set) for a in i)
            ]
        if exclude:
            idx_select = [i for i in idx_select if not any(exclude in a for a in i)]
        return chC.loc[idx_select].sort_values(
            by=[i for i in chC.columns if "corr_val" in i][0]
        )
    else:
        return chC


def take_same_Ev(AllCan, allowed_diff=200):
    _take, _not = [], []
    for n, r in AllCan.iterrows():
        _num_idx0, _num_idx1 = -1, -1
        if n[0].startswith("EIS"):
            _num_idx0 = -2
        if n[1].startswith("EIS"):
            _num_idx1 = -2
        _n0, _n1 = n[0].split("_")[_num_idx0], n[1].split("_")[_num_idx1]

        if all([_n0.isnumeric(), _n1.isnumeric()]):
            if np.abs(float(_n0) - float(_n1)) < allowed_diff:
                _take.append(n)
            else:
                _not.append(n)
        else:
            _take.append(n)
    print(f"Rows left out: {len(AllCan) - len(_take)}")
    return AllCan.loc[_take]


# ========= TDOD Main testing of correlations STARTS HERE --====---- ===========
def sort_corrs(AllCorrs):
    swp = "anodic"
    corvalnm = f"corr_val_{swp}"
    AllCan = AllCorrs.get(swp)
    AllChwoJOS5 = Allchars.query('SampleID != "JOS5"')
    AllCanwoJOS5 = allchars_corr(AllChwoJOS5).get(swp)
    Allcan_EV = take_same_Ev(AllCan, allowed_diff=200)
    AllCanwoJOS5_EV = take_same_Ev(AllCanwoJOS5, allowed_diff=200)

    Allcan_EV_dE50 = take_same_Ev(AllCan, allowed_diff=50)
    AllCanwoJOS5_EV_dE50 = take_same_Ev(AllCanwoJOS5, allowed_diff=50)

    chEA = uniq(AllCan, "EA", corvalnm)
    chBET_SA = select(AllCan, "P_", exclude="")
    chBET_SA = select(chEA, "P_", exclude="P_MW")
    char_cols = ["P", "EA", "B", "R"]
    #    for n,r in chBET_SA.head(10).iterrows():
    #        plot(Allchars,x=n[0],y=n[1])
    #    for elem in ['C','H','N']:
    #        plot(Allchars,x=f'P_prec_{elem}',y=f'EA_{elem}_content',equal=True,savefig=mkfolder(plotsfolder.joinpath('prep_ea_out')))
    #    plot(Allchars,x=f'P_prec_100-CNH',y=f'EA_100-CHN',equal=True,savefig=mkfolder(plotsfolder.joinpath('prep_ea_out')))
    #   main_ch,second_ch ='N2', 'B',
    #   main_ch,second_ch ='ORR', '*',
    #    AllCan, main_ch,second_ch, Allchars = Allcan_EV_dE50,'*', 'EIS', Allchars #, include = char_cols + eisplot.parlst,exclude_second = ['linKK'],sep_topdir=f'EIS_char_best',corrval_zero = 1,cutoff_set=0.6
    def combo_makers(
        AllCan,
        main_ch,
        second_ch,
        Allchars,
        exclude_second="",
        include="",
        endswith="",
        separate_dir="",
        corrval_zero=1,
        cutoff_set=0,
        sep_topdir="",
        logy_set=False,
        line=False,
        force_replot=True,
        max_plots=100,
    ):
        chEA = uniq(AllCan, (main_ch, second_ch), corvalnm)
        chEA_excl = select(chEA, second_ch, endswith_set=endswith)

        if type(exclude_second) == type(""):
            exclude_second = [exclude_second]
        if type(exclude_second) == type([]):
            for excl in exclude_second:
                chEA_excl = select(chEA_excl, second_ch, exclude=excl)

        if type(include) == type(""):
            include = [include]
        if type(include) == type([]):
            for incl in include:
                chEA_excl = select(chEA_excl, incl)
        #                print(len(chEA_excl))
        #        chBET_SA = select(chEA, second_ch, exclude = exclude_second)
        #        chBET_SA = select(chBET_SA , '',exclude = 'Volume')
        cutoff = 0.4 if len(chEA_excl) > 300 else 0.1
        if cutoff_set:
            cutoff = cutoff_set
        corr_cutoff = chEA_excl[chEA_excl[corvalnm].abs() > cutoff]
        if len(corr_cutoff) > max_plots:
            corr_cutoff = pd.concat(
                [chEA_excl.head(int(max_plots / 2)), chEA_excl.tail(int(max_plots / 2))]
            )
        print(
            f"Starting to plot ({main_ch}, {second_ch};+{include}+ ends({endswith}), -{exclude_second}): {len(corr_cutoff)}"
        )
        for n, r in corr_cutoff.iterrows():

            combo_dest_dir = mkfolder(
                plotsfolder.joinpath(f"corrs_{r.chars[0]}_{r.chars[1]}")
            )
            if separate_dir:
                combo_dest_dir = mkfolder(combo_dest_dir.joinpath(separate_dir))
            if sep_topdir:
                combo_dest_dir = mkfolder(plotsfolder.joinpath(f"corrs_{sep_topdir}"))
            #            if corrval_zero == 0:
            _corr_val_set = 0 if corrval_zero != 1 else r[corvalnm]
            plot(
                Allchars,
                x=n[0],
                y=n[1],
                corr_val=_corr_val_set,
                savefig=combo_dest_dir,
                logy=logy_set,
                add_line=line,
                force_replot=force_replot,
            )
            plt.clf()
            plt.close("all")
        return chEA_excl

    def all_combos_chars_EC(EC_name):
        _dict = {}
        for Charnm in ["P", "EA", "B", "R"]:
            N2P = combo_makers(AllCan, EC_name, Charnm, Allchars, exclude_second="")
            _dict.update({f"{EC_name}_{Charnm}": N2P})
        return _dict

    P_EA = combo_makers(AllCan, "P", "EA", Allchars, corrval_zero=0)
    P_P = combo_makers(AllCan, "P", "P", Allchars)

    R_EA = combo_makers(AllCan, "R", "EA", Allchars, exclude_second="2nd")
    R_P = combo_makers(
        AllCanwoJOS5,
        "P",
        "R",
        AllChwoJOS5,
        exclude_second="2nd",
        include="Metal",
        separate_dir="Metal2",
    )

    R_EA2 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="2nd",
        include="N_content",
        separate_dir="N_content",
    )

    R_EA4 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="D4D4",
        include="D4",
        separate_dir="D4_peak",
    )
    R_EA3 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="",
        include="D3",
        separate_dir="D3_peak",
    )
    R_EA3 = combo_makers(
        AllCanwoJOS5,
        "EA",
        "R",
        Allchars,
        exclude_second="",
        include="G",
        separate_dir="G_peak",
    )

    N2_combos = all_combos_chars_EC("N2")
    N2_B_Extra = combo_makers(
        AllCan,
        "N2",
        "B",
        Allchars,
        include=["N2_Cdl", "B_BET_RAW_calc_SA"],
        separate_dir="Cdl_extra",
        corrval_zero=0,
    )
    N2_B_Extra = combo_makers(
        AllCan,
        "N2",
        "B",
        Allchars,
        include=["N2_Cdl", "B_BET_RAW_calc_constant"],
        separate_dir="Cdl_extra_C",
        corrval_zero=0,
        logy_set=True,
    )

    ### EIS correlations
    #    EIS_combos = all_combos_chars_EC('EIS')
    N2_EIS_dE50 = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_nAd"],
        separate_dir="N2_eis_nAd_dE50",
        corrval_zero=0,
    )
    N2_EIS_dE50 = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_nDL"],
        separate_dir="N2_eis_nDL_dE50",
        corrval_zero=0,
    )

    N2_EIS_dE50 = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "Qad+Cdlp"],
        separate_dir="N2_eis_Qad+Cdlp_dE50",
        corrval_zero=0,
    )
    N2_EIS_dE50_Qad = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_Qad_"],
        separate_dir="N2_eis_Qad_dE50",
        corrval_zero=0,
    )
    N2_EIS_dE50_Qad = combo_makers(
        Allcan_EV_dE50,
        "N2",
        "EIS",
        Allchars,
        include=["N2_Cdl", "EIS_Cdlp_"],
        separate_dir="N2_eis_Cdl_dE50",
        corrval_zero=0,
    )

    best_EIS_char = combo_makers(
        Allcan_EV_dE50,
        "*",
        "EIS",
        Allchars,
        include="",
        exclude_second=["linKK"],
        sep_topdir=f"EIS_char_best",
        corrval_zero=1,
        cutoff_set=0.9,
    )

    def make_all_EIS_combos():
        for par in eisplot.parlst[:3]:
            not_par = [i for i in eisplot.parlst if i not in par]
            #        for C in ['Qad_','_Cdlp_','Qad+Cdlp']:
            plt.close("all")
            pdir = f"EIS_pars/{par}"
            pexcl = ["linKK"] + not_par

            for gas in ["O2", "N2"]:
                for orrC in [
                    "J_diff_lim",
                    "_Jkin_min_600",
                    "H2O2_500",
                    "n_ORR",
                    "J_ring",
                ][:]:
                    ORR_EIS_Extra = combo_makers(
                        AllCanwoJOS5_EV_dE50,
                        "ORR",
                        "EIS",
                        Allchars,
                        include=[orrC, par],
                        exclude_second=pexcl,
                        sep_topdir=f"{pdir}/{orrC}/{gas}",
                        corrval_zero=1,
                    )
                #            plt.close('all')
                #                N2P = combo_makers(Allcan_EV_dE50,'B', 'EIS', Allchars, include = ['B', par, 'B_BET_RAW_calc_SA m2/g'],exclude_second = pexcl, endswith = gas, sep_topdir=f'{pdir}/{Charnm}_BET/{gas}',corrval_zero = 1)
                #                N2P = combo_makers(Allcan_EV_dE50,'EA', 'EIS', Allchars, include = ['EA', par, 'EA_N_content'],exclude_second = pexcl, endswith = gas, sep_topdir=f'{pdir}/{Charnm}_EA_N_content/{gas}',corrval_zero = 1)
                #                N2P = combo_makers(Allcan_EV_dE50,'EA', 'EIS', Allchars, include = ['EA', par, 'EA_C_content'],exclude_second = pexcl, endswith = gas, sep_topdir=f'{pdir}/{Charnm}_EA_C_content/{gas}',corrval_zero = 1)
                N2P = combo_makers(
                    Allcan_EV_dE50,
                    "EA",
                    "EIS",
                    Allchars,
                    include=["EA", par, "EA_H_content"],
                    exclude_second=pexcl,
                    endswith=gas,
                    sep_topdir=f"{pdir}/{Charnm}_EA_H_content/{gas}",
                    corrval_zero=1,
                )
                for Charnm in char_cols:
                    N2P = combo_makers(
                        Allcan_EV_dE50,
                        Charnm,
                        "EIS",
                        Allchars,
                        include=[Charnm, par],
                        exclude_second=pexcl,
                        endswith=gas,
                        sep_topdir=f"{pdir}/{Charnm}/{gas}",
                        corrval_zero=1,
                    )
            plt.close("all")

            #                plt.close()
            #                    AllCan,EC_name, Charnm, Allchars, exclude_second = '')

            N2_EIS_dE50_Qad = combo_makers(
                Allcan_EV_dE50,
                "N2",
                "EIS",
                Allchars,
                include=["N2_Cdl", par],
                exclude_second=pexcl,
                sep_topdir=f"{pdir}",
                corrval_zero=0,
            )
            #            N2_EIS_dE50_Qad = combo_makers(Allcan_EV_dE50,'N2', 'EIS', Allchars, include = ['N2_Cdl','EIS_Qad_'],separate_dir='N2_eis_Qad_dE50',corrval_zero = 0)

            for HPRRc in ["HPRR_dj/dE", "HPRR_E_onset"]:
                HPRR_EIS_Extra = combo_makers(
                    AllCanwoJOS5_EV,
                    "HPRR",
                    "EIS",
                    Allchars,
                    include=[HPRRc, par],
                    exclude_second=pexcl,
                    sep_topdir=f'{pdir}/{HPRRc.replace("/","_")}',
                    corrval_zero=1,
                )
            #                plt.close()
            for orrC in ["J_diff_lim", "_Jkin_min_600", "H2O2_500", "n_ORR", "J_ring"][
                :
            ]:
                ORR_EIS_Extra = combo_makers(
                    AllCanwoJOS5_EV_dE50,
                    "ORR",
                    "EIS",
                    Allchars,
                    include=[orrC, par],
                    exclude_second=pexcl,
                    sep_topdir=f"{pdir}/{orrC}",
                    corrval_zero=1,
                )
            plt.close("all")

    make_all_EIS_combos()

    HPRR_Extra = combo_makers(
        AllCan,
        "HPRR",
        "ORR",
        Allchars,
        include=["HPRR_dj/dE"],
        separate_dir="djdE",
        corrval_zero=0,
    )
    HPRR_Extra = combo_makers(
        AllCan,
        "HPRR",
        "ORR",
        Allchars,
        include=["HPRR_dj/dE", "J_ring"],
        separate_dir="djdE",
        corrval_zero=0,
    )
    HPRR_orr = combo_makers(
        AllCan, "HPRR", "ORR", Allchars
    )  # , include = ['HPRR_dj/dE'],separate_dir='djdE',corrval_zero = 0)

    TS_corr = combo_makers(
        AllCan, "ORR", "EA", Allchars, include=["ORR_TSa_min"], sep_topdir="TSa"
    )
    TS_corr = combo_makers(
        AllCan, "ORR", "*", Allchars, include=["ORR_TSa_max"], sep_topdir="TSa_max"
    )

    TS_corr = combo_makers(
        AllCanwoJOS5, "ORR", "EA", Allchars, include=["ORR_TSa"], sep_topdir="TSa_EA"
    )

    HPRR_combos = all_combos_chars_EC("HPRR")

    HPRR_P = combo_makers(AllCan, "N2", "P", Allchars, exclude_second="")
    HPRR_EA = combo_makers(AllCan, "N2", "EA", Allchars, exclude_second="")
    HPRR_B = combo_makers(
        AllCan, "N2", "B_BET_RAW_calc_SA", Allchars, exclude_second=""
    )
    HPRR_R = combo_makers(AllCan, "N2", "R", Allchars, exclude_second="")

    HPPR_EIS = combo_makers(
        AllCanwoJOS5, "HPRR", "EIS", Allchars, exclude_second="linKK"
    )

    N2_EIS = combo_makers(Allcan_EV, "N2", "EIS", Allchars, exclude_second="linKK")

    ORR_EIS = combo_makers(
        AllCanwoJOS5_EV, "ORR", "EIS", Allchars, exclude_second="linKK"
    )

    plot(
        Allchars,
        x="EA_H_content",
        y="B_BET_RAW_calc_constant",
        equal=True,
        logy=True,
        savefig=mkfolder(plotsfolder.joinpath("prep_ea_out")),
    )

    def plot_only(
        Chars, Corrs, selection, number_plots=100, exclusion="", chars_excl=""
    ):
        if chars_excl:
            Corrs = Corrs.loc[Corrs.chars != ("R", "R")]
        chD2 = select(Corrs, selection, exclude=exclusion)
        for n, r in pd.concat([chD2.head(50), chD2.tail(50)]).iterrows():
            plot(
                Chars,
                x=n[0],
                y=n[1],
                corr_val=r[corvalnm],
                savefig=mkfolder(plotsfolder.joinpath(f"corrs_only_{selection}")),
            )

    plot_only(
        Allchars.query('SampleID != "JOS5"'),
        select(AllCanwoJOS5, "C/N_ratio"),
        "Metal_element",
    )
    plot_only(Allchars, select(AllCanwoJOS5, "C/N_ratio"), "C/N_ratio")

    plot_only(Allchars, AllCan, "2nd_", chars_excl=("R", "R"))

    def D2_only():
        chD2 = select(AllCan, "R_D2_height", exclude="_sigma")
        for n, r in pd.concat([chD2.head(50), chD2.tail(50)]).iterrows():
            plot(
                Allchars,
                x=n[0],
                y=n[1],
                corr_val=r[corvalnm],
                savefig=mkfolder(plotsfolder.joinpath(f"corrs_R_D2height_only")),
            )

    def HPRR_only():
        chD2 = select(AllCanwoJOS5, "HPRR_E_onset", exclude="_sigma")
        for n, r in chD2.iterrows():
            plot(
                Allchars,
                x=n[0],
                y=n[1],
                corr_val=r[corvalnm],
                savefig=mkfolder(plotsfolder.joinpath(f"corrs_HPRR_onset_only")),
            )

    for n, r in AllCan.head(10).iterrows():
        plot(swgrp, x=n[0], y=n[1])

    for ch in uniqchars:
        chc = AllCan.query("chars == @ch").sort_values(corvalnm)
        for n, r in chc.tail(1).iterrows():
            plot(swgrp, x=n[0], y=n[1])

    # == RAMAN corrs ==
    #    uniqR = [i for i in uniqchars if 'R' in i and i != ('R','R')]
    #    chR = AllCan.loc[AllCan['chars'].isin(uniqR)].sort_values(corvalnm)
    chR2 = uniq(AllCan, "R", corvalnm)
    for n, r in chR.tail(5).iterrows():
        plot(swgrp, x=n[0], y=n[1])
    for n, r in chR.head(5).iterrows():
        plot(swgrp, x=n[0], y=n[1])

    AllC_500 = get_corrs_E(AllCan, 600).sort_values(corvalnm)
    uniqEA = [i for i in uniqchars if "EA" in i]
    uniqEA = [("EA", "EA")]
    chEA = AllCan.loc[AllCan["chars"].isin(uniqEA)].sort_values(corvalnm)
    for n, r in chEA.tail(20).iterrows():
        plot(Allchars, x="WL_precursor", y=n[1])
    AllC_500

    AllC_500 = get_corrs_E(AllCan, 600).sort_values(corvalnm)

    uniqEA = [("B", "B")]
    chEA = AllCan.loc[AllCan["chars"].isin(uniqEA)].sort_values(corvalnm)
    for n, r in chEA.head(20).iterrows():
        plot(Allchars, x=n[0], y=n[1])


def get_corrs_E(AllCorrs_test, Ev):
    lst, sets = [], []
    for ix in AllCorrs_test.index:
        pt = []
        for n, pos in enumerate(ix):
            last = pos.split("_")[-1]
            match = re.match(r"^[0-9]{3}\b", last)
            if match:
                if float(match.group(0)) == Ev:
                    pt.append([n, pos])  # Not a Potential
                #                    lst += ix
                else:
                    pass  # Wrong Potential Ev
            else:
                pt.append([n, pos])
        if len(pt) == 2:
            idx = (pt[0][1], pt[1][1])
            _ridx = (pt[1][1], pt[0][1])
            #            sets.append(set(idx))
            if not _ridx in lst:
                lst.append(idx)
    return AllCorrs_test.loc[lst]


def make_uniform_EvRHE(df, rounding_set=2):
    lst = []
    for E in df["E_AppV_RHE"].values:
        match = 0
        for i in np.arange(-0.10, 1.2, 0.05):
            if np.isclose(i, E, atol=0.025):
                match = 1
                lst.append((E, i))
        if match == 0:
            if E < 0 and E > -0.04:
                lst.append((E, i))
            else:
                print(E, i)

    if len(df["E_AppV_RHE"].values) == len(lst):
        df = df.assign(**{"E_RHE": [np.round(float(i[1]), rounding_set) for i in lst]})
        print(
            'Len({0}) matches, new column: "E_RHE"'.format(len(df["E_AppV_RHE"].values))
        )
    else:
        print(
            "make_uniform_EvRHE lengths do not match LenEIS : {0}, len(lst) : {1}".format(
                len(df["E_AppV_RHE"].values), len(lst)
            )
        )
    return df


# def plot_chars_Ev(Allchars)
#         PlotAnalysis.corrplot(allchflrtd[allcorrs.name].unstack().dropna())


def remove_corrdups(allchrrs):
    sl = allchrrs.loc[allchrrs.chars == ("ORR", "ORR")]
    sl_non_idx = allchrrs.loc[~allchrrs.index.isin(sl.index)]
    non_dupl = [
        i
        for i in sl.index
        if len(set(["_".join(i[0].split("_")[1:-1]), "_".join(i[1].split("_")[1:-1])]))
        > 1
    ]
    rings = ["J_ring", "Frac_H2O2", "n_ORR"]
    ring_dups = [
        i
        for i in non_dupl
        if not (any(r in i[0] for r in rings) and any(r in i[1] for r in rings))
    ]
    jkin_dups = [i for i in ring_dups if not (("Jkin_" in i[0] and "Jkin_" in i[1]))]
    return allchrrs.loc[list(sl_non_idx.index) + jkin_dups].sort_index()


if __name__ == "__main__":

    EvRHE = "E_AppV_RHE"
    templ = PorphSiO2_template()
    ECpars, raman, bet, ea, prep, struc, char, Allchars = collect_pars()
    export_chars(char)
    export_chars({"all_chars": Allchars})

    list(struc.columns)
    list(ECpars.columns)
    corrs = corrloop(char)

    AllCorrs = allchars_corr(Allchars)
    corr_prep_ea = postCorrs.prep_ea(prep, ea)
    corr_prep_ea = postCorrs.prep_ea(ea, raman)

    makeplots = 0
    if not ECpars.empty and makeplots == True:

        for i in np.arange(100, 1000, 100):
            plot(ECpars, x=f"ORR_Frac_H2O2_{i}", y="HPRR_dj/dE", savefig=plotsfolder)

        plot(Allchars, x="EA_N_content", y="ORR_Frac_H2O2_400")
        plot(ECpars, x="ORR_Frac_H2O2_500", y="HPRR_dj/dE")

        plot(ECpars, x="ORR_Frac_H2O2_500", y="HPRR_E_onset")
        plot(ECpars, x="ORR_Frac_H2O2_500", y="ORR_Jkin_min_500")
        plot(ECpars, x="N2_Cdl_mFcm-2_700", y="ORR_Jkin_min_750")
        plot(ECpars, x="N2_Cdl_mFcm-2_700", y="ORR_Frac_H2O2_500")

        plot(ECpars, x="N2_Cdl_mFcm-2_500", y="ORR_n_ORR_500")
        #        === POWERPOINT slides
        ppslides = mkfolder(FindExpFolder("PorphSiO2").compare.joinpath("pp_slides"))
        plot(ECpars, x="HPRR_E_onset", y="HPRR_dj/dE", savefig=ppslides)
        plot(ECpars, x="ORR_E_onset", y="ORR_Jkin_min_750", logy=True, savefig=ppslides)

        plot(ECpars, x="ORR_E_onset", y="HPRR_E_onset", savefig=ppslides)
        plot(ECpars, x="HPRR_dj/dE", y="ORR_Jkin_min_750", savefig=ppslides)

        plot(ECpars, x="ORR_Frac_H2O2_500", y="ORR_Jkin_min_750", savefig=ppslides)
        plot(ECpars, x="ORR_Frac_H2O2_500", y="HPRR_E_onset", savefig=ppslides)

        plot(ECpars, x="N2_Cdl_mFcm-2_500", y="HPRR_dj/dE")  # , savefig=plotsfolder)
        ECpars["relHPRR"] = ECpars["HPRR_dj/dE"] / ECpars["N2_Cdl_mFcm-2_500"]
        plot(ECpars, x="ORR_Jkin_min_750", y="relHPRR")  # , savefig=plotsfolder)

        plot(ECpars, x="N2_Cdl_mFcm-2_500", y="ORR_J_diff_lim", savefig=plotsfolder)

        plot(
            Allchars,
            y="B_BET_RAW_calc_constant",
            x="P_MW_g-mol",
            savefig=plotsfolder,
            logy=True,
        )

        #        plot(Allchars,x='B_tBoer_microArea',y='P_MW_content',savefig=plotsfolder)
        plot(ECpars, x="HPRR_E_onset", y="HPRR_dj/dE", savefig=plotsfolder)
        plot(
            ECpars,
            x="ORR_E_onset",
            y="ORR_Jkin_min_750",
            logy=False,
            savefig=plotsfolder,
        )
        plot(
            Allchars,
            x="B_BET_RAW_calc_SA m2/g",
            y="HPRR_dj/dE",
            logy=False,
            savefig=plotsfolder,
        )
        plot(
            Allchars,
            x="B_BET_RAW_calc_SA m2/g",
            y="P_WL_precursor",
            logy=False,
            savefig=plotsfolder,
        )


#    ECpars = EC_PorphSiO2().pars
#    pst = postChar()
#    raman,bet,ea = pst.raman, pst.bet, pst.ea
#    struc = pst.merged


#
#    fig,ax= plt.subplots(figsize= (16,16))
#    scatter_matrix(ECpars[[key for key,val in ECpars.dtypes.to_dict().items() if 'float64' in str(val)]], alpha=0.5, ax=ax, diagonal='kde')
#    plt.savefig(FindExpFolder('PorphSiO2').compare.joinpath('EC_scmatrix.png'))
#    plt.close()
