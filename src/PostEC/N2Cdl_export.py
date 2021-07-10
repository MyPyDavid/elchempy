#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:04:01 2020

@author: zmg
"""

from scipy.stats import linregress


from FileHelper.PostChar import SampleSelection, Characterization_TypeSetting
import post_helper
import plotting


def N2export_to_Origin():
    PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")

    SeriesIDs = [
        SampleSelection.Series_CB_paper,
        SampleSelection.Series_Porhp_SiO2,
        SampleSelection.Series_Co_PANI,
        SampleSelection.Series_ML_SiO2,
        {"name": "all_EIS"},
    ]
    SeriesID_set = SeriesIDs[1]
    PPDN2 = PostDestDir.joinpath("N2CV_{0}".format(SeriesID_set["name"]))
    PPDN2.mkdir(parents=True, exist_ok=True)

    PPDCV = PPDN2.joinpath("N2_CV")
    PPDCV.mkdir(parents=True, exist_ok=True)
    Cdl_pars
    CdlSIO2 = Cdl_pars.loc[Cdl_pars.SampleID.isin(SeriesID_set["sIDs"])]
    N2_CV_index = postOVVout.groupby("Type_output").get_group("N2_CV")
    N2idx = N2_CV_index.loc[
        N2_CV_index.SampleID.isin(SeriesID_set["sIDs"])
        & N2_CV_index.Source.str.contains("ExpDir")
    ]
    Acid_slice = CdlSIO2.query(
        'pH < 7 & Electrode != "Pt_ring" & postAST == "no" & Loading_name == "standard"'
    )
    N2_acid_no = Acid_slice.loc[~Acid_slice.ECexp.str.contains("postAST")]
    #%%
    for sID, gr in N2_acid_no.groupby(["SampleID"]):
        #        sID,gr
        for (exp, elec, filename), expgrp in gr.groupby(
            ["ECexp", "Electrolyte", "Filename"]
        ):
            exp, elec, expgrp
            if (
                expgrp.Cdl_CV_data_files.values[0]
                == expgrp.Cdl_CV_data_files.values[-1]
            ):
                CVfiles = expgrp.Cdl_CV_data_files.values[0]
            else:
                CVfiles = expgrp.Cdl_CV_data_files.values[-1]

            if len(CVfiles) == 12:
                CV12 = CVfiles
                CVfiles = [i for i in CV12 if "standard" in i]
                if len(CVfiles) != 6:
                    CVprob = sorted(CVfiles)

                    if len(CVprob) == 12:
                        CVfiles = CVprob[0:6]
                    else:
                        print(f"12 to 6 problem  {exp}, len is {len(CVfiles)} ")

            if len(CVfiles) != 6:
                print(f"len problem {exp}, len is {len(CVfiles)}")
                #                print(f'len problem {exp}, len is {len(CVfiles)}')
                lentest = [i for i in expgrp.Cdl_CV_data_files.values if len(i) == 6]
                if lentest:
                    if lentest[0] == lentest[-1]:
                        CVfiles = lentest[0]
                    else:
                        CVfiles = lentest[-1]
                else:
                    CVfiles = []

            if len(CVfiles) != 6:
                print(f"DOUBLE len problem {exp}, len is {len(CVfiles)}")
            #            def export_CVset(CVfiles,expgrp):
            #            [(Path(i).stem.split('_')[-2],i) for i in  CVfiles if not 'first' in i]
            if CVfiles:
                cvsr = sorted(
                    [
                        (Path(i).stem.split("_")[-2], Path(i))
                        for i in CVfiles
                        if not "first" in i
                    ],
                    reverse=True,
                )
                CVexp = FileHelper.FileOperations.ChangeRoot_DF(
                    pd.DataFrame(cvsr, columns=["ScanRate", "File"]),
                    [],
                    coltype="string",
                )
                pdlst, pdmg = [], pd.DataFrame()
                N2origin = pd.DataFrame()
                for n, r in CVexp.iterrows():
                    try:
                        CVread_raw = pd.read_excel(r.File, index_col=[0])
                        j_mA_col = "j mA/cm2"
                        CVread_mA = convert_mAcm2(CVread_raw, j_mA_col)

                        SR = CVread_mA.rename(
                            columns={j_mA_col: f"CV_{r.ScanRate}"}
                        ).reset_index(drop=True)
                        #                        SR[f'CV_{r.ScanRate}'] = SR[f'CV_{r.ScanRate}']*1000 # TODO FORGET about mA
                        pdlst.append(SR)
                        if pdmg.empty:
                            pdmg = SR
                        else:
                            pdmg = pd.merge(
                                pdmg,
                                SR[f"CV_{r.ScanRate}"],
                                left_index=True,
                                right_index=True,
                            )
                    #                        pd.merge(pdmg,SR, on = EvRHE)
                    except:
                        print(f"read problem {exp} ")
                N2origin = pdmg.set_index(EvRHE)
                N2cdl = calc_cdl(N2origin)
                #                [pd.read_excel(i[1]) for i in cvsr]
                PPDelec = make_subdir(PPDCV.joinpath(elec))
                sID, metalprec = (
                    expgrp.SampleID.unique()[0],
                    expgrp.MetalPrec.unique()[0],
                )
                metal = metalprec[0:2]
                xldate = exp.split("_")[-1]
                extra = "_".join(Path(filename).stem.split("_")[-2:])
                ShortLabel = f"{sID}_{metal}_CVs"
                fpath = f"{ShortLabel}_{xldate}_{extra}_{exp}.png"
                N2origin.plot(title=f"{ShortLabel},{extra}\n{exp}")
                plt.savefig(
                    PPDelec.joinpath(fpath).with_suffix(".png"),
                    bbox_inches="tight",
                    dpi=100,
                )
                plt.close()

                xlpath = f"{ShortLabel}_{xldate}.xlsx"
                if PPDelec.joinpath(xlpath).is_file():
                    xlpath = f"{ShortLabel}_{xldate}_{extra}.xlsx"

                with pd.ExcelWriter(PPDelec.joinpath(xlpath)) as writer:
                    N2cdl.iloc[:1000].to_excel(writer, sheet_name="cathodic")
                    N2cdl.iloc[1000:].to_excel(writer, sheet_name="anodic")


#
#                N2cdl.to_excel(PPDelec.joinpath(xlpath))

#
#                cdl_fit = linregress(gr['ScanRate'], gr['j_anodic'])
#
#                for swp, swpgrp in expgrp.groupby('Sweep_Type_N2'):
#                    swp,swpgrp
#%%
def convert_mAcm2(CVread, j_mA_col):
    jAcm2 = "j A/cm2"
    if jAcm2 in CVread.columns:
        CVread[j_mA_col] = CVread[jAcm2] * 1000
        CVread = CVread.drop(columns=jAcm2)
    return CVread


def calc_cdl(N2origin):
    """Typical values for carbon supercapacitor should be 10-100 F/g , or 100 mF/cm2
    double-layer capacitance on eg. Pt 20 uF/cm2."""
    lst = []
    for n, r in N2origin.iterrows():
        n, r
        srs = [float(i.split("_")[-1]) for i in r.index]
        if np.min(srs) > 5:
            sr_unit = "mV/s"
            srs = [i * 1e-3 for i in srs]
        r.index = srs
        xy = r.sort_index()
        ln = linregress(xy.index, xy.values)  # *1E-3 TODO FORGET about mA
        lst.append((n, ln.slope))
    Cdl_origin = pd.DataFrame(lst, columns=[EvRHE, "Cdl_mFcm-2"]).set_index(EvRHE)
    N2cdl = pd.concat([N2origin, Cdl_origin], axis=1)
    return N2cdl


def make_subdir(DD):
    DD.mkdir(parents=True, exist_ok=True)
    return DD


class ExportECfromN2:
    def __init__(self, **kwargs):
        if "reload" in kwargs:
            self.postOVVout = PostEC.LoadPostOVV(kwargs["reload"])
            print(
                "Exp types found in overview: {0}".format(
                    ", ".join([str(i) for i in self.postOVVout.Type_Exp.unique()])
                )
            )
            pass

    def out_Cdl_pars():
        E_range = np.arange(0, 1, 0.1)
        LK_cdl = Cdl_pars.loc[Cdl_pars.SampleID.str.contains("LK")]
        LK_cdlE = pd.concat(
            [LK_cdl.loc[np.isclose(i, LK_cdl.E_AppV_RHE, atol=0.01)] for i in E_range],
            sort=False,
        )
        LK_dest_path = PostEC().DestDir.joinpath("LK_samples")
        LK_dest_path.mkdir(exist_ok=True, parents=True)
        for n, gr in LK_cdlE.groupby(["Sweep_Type_N2", "E_RHE"]):
            #            gr.plot.bar(x='SampleID',y=('Cdl') )
            fig, ax = plt.subplots(figsize=(6, 6))
            gr_describe = (
                gr[SampleSelection.EC_N2Cdl_cols + ["SampleID"]]
                .groupby("SampleID")
                .describe()
            )
            gr_describe.plot.bar(
                y=("Cdl", "mean"),
                yerr=gr_describe.loc[:, ("Cdl", "std")],
                ylim=(0, 0.015),
                title=f"{n[0]} at {n[1]} Vrhe",
                rot=60,
                ax=ax,
            )
            leg1 = ax.get_legend_handles_labels()[1]
            ax.legend(leg1)
            plt.savefig(
                LK_dest_path.joinpath(f"{n[0]}_{n[1]}.jpg"),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()

    @staticmethod
    def N2_scans_Origin():
        PostDestDir = FileHelper.FindExpFolder("VERSASTAT").DestDir.joinpath("PostEC")
        PPD_N2CV = PostDestDir.joinpath("N2_CV")
        N2_CV_ovv = postOVVout.query('Type_Exp == "N2_CV"')
        N2sfs = [
            "_".join(i.stem.split("_")[0:-1]) for i in N2_CV_ovv.SourceFilename.values
        ]
        N2_CV_ovv = N2_CV_ovv.assign(**{"N2_Sf": N2sfs})
        N2_CVs = N2_CV_ovv.loc[
            N2_CV_ovv.SourceFilename.isin(
                [
                    i
                    for i in N2_CV_ovv.SourceFilename.values
                    if not i.name.endswith("_first_v20.xlsx")
                ]
            ),
            :,
        ]
        #    N2sfs = set(['_'.join(i.split('_')[0:-1]) for i in N2_CV_ovv.SourceFilename.values if 'first' not in i.split('_')[-1]])

        for Grnm, N2gr in N2_CVs.groupby(by="N2_Sf"):
            N2lst = []
            PPD_N2_Elec = PPD_N2CV.joinpath(N2gr.Electrolyte.unique()[0])
            sID, Date, Status = (
                N2gr.SampleID.unique()[0],
                N2gr.EXP_date.unique()[0],
                N2gr.Status.unique()[0],
            )
            extra = (
                Path(Grnm).stem.split(N2gr.SampleID.unique()[0])[-1][1::] + "_OriginSRs"
            )
            for nm, N2f in N2gr.groupby(by="SourceFilename"):
                N2sr = pd.read_excel(nm)
                if "j A/cm2" in N2sr.columns:
                    N2sr["j mA/cm2"] = N2sr["j A/cm2"] * 1000
                    N2sr.drop(columns="j A/cm2", inplace=True)

                N2sr.rename(
                    columns=dict(
                        zip(
                            N2sr.columns,
                            [
                                i + "_%s" % int(N2f.Scanrate.values[0])
                                for i in N2sr.columns
                            ],
                        )
                    ),
                    inplace=True,
                )
                N2lst.append(N2sr.reset_index(drop=True))
                print(Grnm)
            N2sr_Out = pd.concat(N2lst, sort=False, axis=1)
            PPD_N2_out = PPD_N2_Elec.joinpath(
                "_".join([sID, Status, Date, extra])
            ).with_suffix(".xlsx")
            N2sr_Out.to_excel(PPD_N2_out)
