import numpy as np
import pandas as pd


class EC_Properties:
    """Some properties of electrochemical experiments collected and reference values"""

    def __init__(self):
        self.electorode_SA = self.WE_surface_area_cm2()

    def WE_surface_area_cm2(self, TYPE="PINE"):
        coll_eff = []
        if TYPE == "ALS":
            r1, r2, r3 = 0.1 * 4 * 0.5, 0.1 * 5 * 0.5, 0.1 * 7 * 0.5
            SA = np.pi * (r3 ** 2 - r2 ** 2)
            r1, r2, r3 = 0.1 * 4 * 0.5, 0.1 * 5 * 0.5, 0.1 * 7 * 0.5
            SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
        elif TYPE == "PINE":
            r1, r2, r3 = 0.1 * 5.5 * 0.5, 0.1 * 6.5 * 0.5, 0.1 * 8.5 * 0.5
            coll_eff = 0.38
            SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
        elif TYPE == "PINE-ring":
            r1, r2, r3 = 0.1 * 5.5 * 0.5, 0.1 * 6.5 * 0.5, 0.1 * 8.5 * 0.5
            coll_eff = 0.38
            SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
        #        SA = np.pi*(r3**2-r2**2)
        if coll_eff == []:
            a, b = (r2 / r1) ** 3 - 1, (r3 / r1) ** 3 - (r2 / r1) ** 3
            c = a / b
            coll_eff = (
                1
                - r3 ** 2
                + b ** (2 / 3)
                - coll_func(c)
                - b ** (2 / 3) * coll_func(a)
                + r3 ** 2 * coll_func(c * r3 ** 3)
            )
        "r1 = disk, r2 = ring ID, r3 = ring OD in mm"
        #    print('%s According to manufacturer: disk(dia:%.2f cm   %.4f cm2), ring (%.4f cm2)' %(TYPE,r1*2,SAdisk,SAring))
        return {
            "Typename": TYPE,
            "CollEff": coll_eff,
            "Disk_cm2": SAdisk,
            "Ring_cm2": SAring,
        }

    @staticmethod
    def loading_ref(PAR_date):
        WE_SA = EC_Properties.WE_surface_area_cm2()["Disk_cm2"]
        if PAR_date < pd.to_datetime("20190901") and PAR_date > pd.to_datetime(
            "20190827"
        ):
            loading_var = {
                "high": ((9.5 / 250) * 5) / WE_SA,
                "normal": ((5 / 250) * 5) / WE_SA,
                "low": ((2.5 / 250) * 5) / WE_SA,
                "standard": (5 / 500) * 9 / WE_SA,
                "half": 0.5 * (5 / 250) * 5 / WE_SA,
                "loading-ref-unknown": (5 / 500) * 9 / WE_SA,
            }
        elif PAR_date <= pd.to_datetime("20190911") and PAR_date >= pd.to_datetime(
            "20190909"
        ):
            loading_var = {
                "high": ((9.5 / 250) * 5) / WE_SA,
                "normal": ((5 / 250) * 5) / WE_SA,
                "low": ((2 / 250) * 5) / WE_SA,
                "standard": (5 / 250) * 5 / WE_SA,
                "half": 0.5 * (5 / 250) * 5 / WE_SA,
                "loading-ref-unknown": (5 / 500) * 9 / WE_SA,
            }
        elif PAR_date == pd.to_datetime("20190912"):
            loading_var = {
                "high": ((9.5 / 250) * 5) / WE_SA,
                "normal": ((5 / 250) * 5) / WE_SA,
                "low": ((2 / 250) * 5) / WE_SA,
                "standard": (5 / 250) * 5 / WE_SA,
                "half": 0.5 * (5 / 250) * 5 / WE_SA,
                "loading-ref-unknown": (5 / 500) * 9 / WE_SA,
            }
        else:
            loading_var = {
                "high": ((9.5 / 250) * 5) / WE_SA,
                "normal": ((5 / 250) * 5) / WE_SA,
                "low": ((2 / 250) * 5) / WE_SA,
                "standard": (5 / 500) * 9 / WE_SA,
                "half": 0.5 * (5 / 250) * 5 / WE_SA,
                "loading-ref-unknown": (5 / 500) * 9 / WE_SA,
            }
        return loading_var

    @staticmethod
    def guess_RHE_from_Electrolyte(Elec):
        if "KOH" in Elec:
            RHE_OCP = 900
        elif "0.1MH2SO4" in Elec:
            RHE_OCP = 298
        elif "0.5MH2SO4" in Elec:
            RHE_OCP = 256
        else:
            RHE_OCP = 0
        return RHE_OCP

    @staticmethod
    def EC_Exp_Cols():
        EC_exp_cols = ["postAST", "Electrolyte", "pH", "Gas", "RPM", "Loading_cm2"]
        EC_EIS_par_cols = [
            "Qad",
            "nAd",
            "Cdlp",
            "nDL",
            "Rct",
            "Rct_kin",
            "Rs",
            "Rorr",
            "RedChisqr2",
            "RedChisqr1",
        ]
        EC_ORR_kin_par_cols = [
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
        return EC_exp_cols + EC_EIS_par_cols + EC_ORR_kin_par_cols

    def PAR_init_text(self, Etype):
        # WE_surface_area = np.pi*(0.55/2)**2
        #'G:\\Cloudstation\\Experimental data\\Raw_data\\VERSASTAT\\Plots'
        # EXTS = ['.par']
        # ext = ['.par']
        # for a in ['ALS','PINE']:
        # coll = self.WE_surface_area_cm2(a)
        _SA_msg = "Surface Area %s: %.3f cm2 \n Collection efficiency %s : %.3f" % (
            Etype,
            coll["Disk_cm2"],
            Etype,
            coll["CollEff"],
        )
        _loading_msg = str(
            "9 ul (0.09 mg cat) of ink on PINE GCdisk (0.238 cm2), LoadingRRDE = 0.378 mgcat/cm2"
            "\nIn all (DW16-21 & DW24-29 experiments 5 µl H2O was added after the ink was dropped."
            "\nA/g = 1000*Jkin /(1000*0.378*0.01*Fe-loading(wt%))"
        )
        return _SA_msg + "\n" + _loading_msg


if __name__ == "__main__":
    ecp = EC_Properties()
