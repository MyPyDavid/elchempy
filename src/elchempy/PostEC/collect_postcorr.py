#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 13:59:55 2021

@author: zmg
"""


class postChar:

    suffixes = ["R", "P", "EA", "B"]

    def __init__(self):
        self.folder = FindExpFolder("PorphSiO2").folder
        #        FindExpFolder().TopDir.joinpath(Path('Preparation-Thesis\SiO2_projects\SiO2_Me_EC+Struc'))
        self.raman = self.postRaman
        self.bet = self.postBET
        self.ea = self.postEA
        self.prec_ea_ratios()
        self.merged = self.merged()

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
        betdir = FindExpFolder("BET").DestDir
        bet_files = betdir.rglob("BET_pars_index*pkl")
        BET_ovv = pd.concat([pd.read_pickle(i) for i in bet_files])
        BET_ovv_template = BET_ovv
        # .loc[BET_ovv.SampleID.isin(PorphSiO2_template().SampleID.unique())]
        BET_ovv_template = BET_ovv_template.loc[
            BET_ovv_template.fullPath.str.endswith("RAW")
            & ~BET_ovv_template.fullPath.str.contains("_FAIL")
        ]
        return BET_ovv_template, ""

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

    def prec_ea_ratios(self):
        templ = PorphSiO2_template()
        ea = self.ea
        prep = self.postPrep
        # self.prep = self.prec_ea_ratios(, self.ea, postChar.postPrep)
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
        self.prep = prep
        # return prep
