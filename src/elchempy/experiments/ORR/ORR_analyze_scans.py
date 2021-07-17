# import sys
from pathlib import Path
from collections import namedtuple
from datetime import datetime
import numpy as np
import logging

from scipy.stats import linregress

# from scipy.interpolate import UnivariateSpline
from scipy import constants

# from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

# import os
# import multiprocessing
# from functools import partial
# from itertools import repeat
import pandas as pd

print("File", __file__, "\nName;", __name__)
# from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations


if __name__ == "__main__":

    from elchempy.experiments.Loader.CreateCV import create_CVs

    # from elchempy.experiments.EC_conditions.electrode import WE_SA_collection_eff
    from O2_chrono import O2_Chrono
    from plotting import ORR_plot_ring_disk
    from KL_calc import KL_plots
#    from EC_logging_config import start_logging
#    pass
else:
    from ..EC_DataLoader.CreateCV import create_CVs as create_CVs

    # from .EC_conditions.electrode import WE_SA_collection_eff
    from .O2_chrono import O2_Chrono
    from .plotting import ORR_plot_ring_disk
    from .KL_calc import KL_plots
#    from ..Loader.CreateCV import create_CVs as create_CVs


_logger = logging.getLogger(__name__)
#
# All_ovv.loc[(All_ovv['Elapsed Time(s)'].isin(Te['Elapsed Time(s)'])) &  (All_ovv['Type'] =='Chronoamperometry') & (All_ovv['DATE'].isin(O2_act.loc[O2_act['hash'] == h1,'DATE'].unique()))]
#    n= 4*ID/ (ID +IR/N)


def ORR_extract_pars(ORR_N2_corrected, **kwargs):
    """
    Extract values of interest from an ORR
    """
    #    ORR_N2_corrected = Jkin_RRDE_an
    E_onset, Diff_lim = pd.DataFrame([]), pd.DataFrame([])
    sweep_col = [
        i
        for i in ORR_N2_corrected.columns
        if "SWEEP" in i.upper() and not "RING" in i.upper()
    ][0]
    EvRHE = [
        i
        for i in ORR_N2_corrected.columns
        if "E_APPV_RHE" in i.upper() and not "RING" in i.upper()
    ][0]
    RPM_col = [i for i in ORR_N2_corrected.columns if "RPM" in i.upper()][0]
    Segment_col = [
        i
        for i in ORR_N2_corrected.columns
        if "SEGMENT" in i.upper() and not "RING" in i.upper()
    ][0]
    if "PAR_file" in ORR_N2_corrected.columns:
        PAR_file = ORR_N2_corrected.PAR_file.unique()[0]
    elif "File_disk" in ORR_N2_corrected.columns:
        PAR_file = ORR_N2_corrected.File_disk.unique()[0]
    else:
        PAR_file = "TF_PAR_file_collumn_missing"

    #    ORR_N2_corrected
    output_pars = {}
    for (rpm, seg, swp), O2_join in ORR_N2_corrected.groupby(
        [RPM_col, Segment_col, sweep_col]
    ):

        next_row = {
            "SampleID": O2_join["SampleID"].unique()[0],
            "PAR_file": PAR_file,
            "DATE": O2_join["DATE"].unique()[0],
            "RPM": rpm,
            "Segment": seg,
            sweep_col: swp,
            "Electrolyte": O2_join.iloc[0]["Electrolyte"],
            "pH": O2_join.iloc[0]["pH"],
            "Analysis_date": datetime.datetime.now(),
        }
        if "postAST" in kwargs.keys():
            next_row.update({"postAST": kwargs.get("postAST", "error")})

        try:
            #                    E_onset = O2_join.loc[(O2_join['Jcorr'] > -0.119) & (O2_join['Jcorr'] < -0.0999) & (O2_join[EvRHE] <(J2nd_deriv_max[EvRHE].values[0]+0.2)) & (O2_join[EvRHE] > (J2nd_deriv_max[EvRHE].values[0]-0.2))& (O2_join[EvRHE] < 1.19),:]
            E_onset = (
                O2_join.loc[
                    (O2_join["Jcorr"] > -0.129)
                    & (O2_join["Jcorr"] < -0.0999)
                    & (O2_join[EvRHE] < 1.19),
                    :,
                ]
                .sort_values(by=EvRHE)
                .head(1)
            )
        except Exception as e:
            #            logger.warning('Jkin calculation ORR; Jkin calc, E_onset expansion: %s' % e)
            i = 0
            while E_onset.empty:
                E_onset = O2_join.loc[
                    (O2_join["Jcorr"] > -0.119 + i)
                    & (O2_join["Jcorr"] < -0.0999 + i)
                    & (O2_join[EvRHE] < 1),
                    :,
                ]
                i += 0.04
        #                Diff_lim = O2_join.loc[np.isclose(O2_join['J_O2_diff'],0,atol=1E-07) | (O2_join['J_O2_diff'] == O2_join['J_O2_diff'].min()) & (O2_join[EvRHE] <  E_onset[EvRHE].mean()),:]

        try:
            #                    E_onset[EvRHE].mean()-0.1
            Diff_lim = O2_join.loc[
                np.isclose(O2_join["J_O2_diff"], -0.003, rtol=50)
                & (O2_join[EvRHE] < 0.15)
                & (O2_join["J_O2_diff"] < 0)
                & (O2_join["Jcorr"] > -8),
                :,
            ]
        except Exception as e:
            #            logger.warning('Jkin calculation, Diff lim problem: %s' % e)
            Diff_lim = O2_join.loc[
                np.isclose(O2_join["J_O2_diff"], -0.003, rtol=1000)
                & (O2_join[EvRHE] < 0.8)
                & (O2_join["J_O2_diff"] < 0)
                & (O2_join["Jcorr"] > -8),
                :,
            ]
        #                    Diff_lim = O2_join.loc[np.isclose(O2_join['J_O2_diff'],-0.003,rtol=100) & (O2_join[EvRHE] <  0.8) & (O2_join['J_O2_diff'] <  0),:]
        E_lowdiff = 0.03
        if E_onset.empty or Diff_lim.empty:
            if E_onset.empty:
                #                logger.warning('Jkin calculation E_onset empty! expanding... %s' % PAR_file)
                i = 0
                while E_onset.empty:
                    E_onset = (
                        O2_join.loc[
                            (O2_join["Jcorr"] > -1 * (0.529 + i))
                            & (O2_join["Jcorr"] < 0.0)
                            & (O2_join[EvRHE] < 1.19),
                            :,
                        ]
                        .sort_values(by=EvRHE)
                        .head(1)
                    )
                    #                            E_onset = O2_join.loc[(O2_join['Jcorr'] > -0.119+i) & (O2_join['Jcorr'] < -0.0999+i) & (O2_join[EvRHE] < 1),:]
                    i += 0.04
            elif Diff_lim.empty:
                #                logger.warning(
                #                    'Jkin calculation Diff_lim empty! %s take Jdiff at %.2f Vrhe' % (PAR_file, E_lowdiff))
                Diff_lim = O2_join.loc[(O2_join[EvRHE] < E_lowdiff), :]
            #                        Diff_lim = O2_join.loc[np.isclose(O2_join[EvRHE],E_lowdiff,rtol=20),:]
        if (
            Diff_lim["Jcorr"].min()
            < O2_join.loc[(O2_join[EvRHE] < E_lowdiff), "Jcorr"].min()
            and Diff_lim[EvRHE].min() > E_lowdiff
        ):
            #            logger.warning(
            #                'Jkin calculation Diff_lim is smaller at higher potentials ! %s taking Jdiff at %.2f Vrhe' % (
            #                PAR_file, E_lowdiff))
            Diff_lim = O2_join.loc[(O2_join[EvRHE] < E_lowdiff), :]

        O2_join = O2_join.assign(
            **{
                "Jkin_max": (
                    -1
                    * (Diff_lim["Jcorr"].max() * O2_join["Jcorr"])
                    / (Diff_lim["Jcorr"].max() - O2_join["Jcorr"])
                ),
                "Jkin_min": (
                    -1
                    * (Diff_lim["Jcorr"].min() * O2_join["Jcorr"])
                    / (Diff_lim["Jcorr"].min() - O2_join["Jcorr"])
                ),
            }
        )

        E_half = O2_join.loc[
            (O2_join["Jcorr"] < Diff_lim["Jcorr"].min() * 0.4980)
            & ((O2_join["Jcorr"] > Diff_lim["Jcorr"].min() * 0.5540))
            & (O2_join[EvRHE] < E_onset[EvRHE].mean()),
            :,
        ]
        if E_half.empty:
            _logger.warning("Jkin calculation E_half empty! expanding... %s" % PAR_file)
            E_half = O2_join.loc[
                (np.isclose(O2_join["Jcorr"], Diff_lim["Jcorr"].min() * 0.5, rtol=1))
                & (O2_join["Sweep_Type"] == "cathodic")
                & (O2_join[EvRHE] < E_onset[EvRHE].mean()),
                :,
            ]
        #                    i=0 i += 1     while E_half.empty:
        #                        E_half = O2_join.loc[(O2_join['Jcorr'] < Diff_lim['Jcorr'].mean()*(0.4980-0.001*i)) & ((O2_join['Jcorr'] > Diff_lim['Jcorr'].mean()*(0.5540+0.001*i))) & (O2_join['Sweep_Type'] == "cathodic")& (O2_join[EvRHE] < E_onset[EvRHE].mean()),:]
        Jkin_075 = O2_join.loc[(np.isclose(O2_join[EvRHE], 0.75, rtol=0.003)), :]
        Jkin_080 = O2_join.loc[(np.isclose(O2_join[EvRHE], 0.8, rtol=0.003)), :]

        Jkin_pars = {
            "E_onset": E_onset[EvRHE].mean(),
            "E_half": E_half[EvRHE].mean(),
            "J_diff_lim": Diff_lim["Jcorr"].min(),
            "Jkin_750": np.abs(Jkin_075["Jkin_min"].mean()),
            "Jkin_800": np.abs(Jkin_080["Jkin_min"].mean()),
        }
        #        'Jring_050': out_Jring, 'FracH2O2_050': out_FracH2O2,
        Tafel_ovv, TFxy, Tafel_pars_out = ORR_Tafel.calculations(O2_join)

        ring_pars_cols, ringpars_out = ["Jkin_min", "J_ring", "Frac_H2O2", "n_ORR"], {}
        if all([i in O2_join.columns for i in ring_pars_cols]):
            for i in np.arange(0, 1, 0.05):
                for ring_par in ring_pars_cols:
                    ring_par_value = np.round(
                        O2_join.loc[
                            np.isclose(O2_join[EvRHE], i, rtol=0.004), ring_par
                        ].mean(),
                        3,
                    )
                    #                    out_FracH2O2 = O2_join.loc[np.isclose(O2_join[EvRHE], i, rtol=0.003), ring_par].mean()
                    if np.isnan(ring_par_value) == False:
                        ringpars_out.update(
                            {f"{ring_par}_{i*1000:.0f}": ring_par_value}
                        )
        next_row.update(Jkin_pars)
        next_row.update(Tafel_pars_out)
        next_row.update(ringpars_out)
    output_pars.update(next_row)
    return output_pars


def tesplot(O2_join):  # noqa: F821
    O2_join.loc[(O2_join[EvRHE] < 0.75) & (O2_join[EvRHE] > 0.4)].plot(
        x="Jkin_min",
        y="Frac_H2O2",
        kind="scatter",
        xlim=(0.1, 10),
        ylim=(10, 300),
        logy=False,
        logx=True,
        c=EvRHE,
    )
    O2_join.loc[(O2_join[EvRHE] < 0.95) & (O2_join[EvRHE] > 0.1)].plot(
        x=EvRHE, y="Jcorr", kind="scatter", xlim=(0.1, 1), logy=False, c=EvRHE
    )


class ORR_Tafel:
    """Class for calculating Tafel plots of the ORR scans"""

    def __init__():
        pass

    def export(
        ORR_dest_dir,
        rTFxy,
        i,
        TFll,
        dest_file,
        TFfit,
    ):  # noqa : F821
        TafelDir = ORR_dest_dir.joinpath("TAFEL")
        TafelDir.mkdir(parents=True, exist_ok=True)

        rTFxy.plot(
            x=EvRHE, y=["log10_Jkin_min", TFll], title="Tafel_%.3f" % TFfit[0] * -1000
        )
        plt.savefig(TafelDir.joinpath("%s_" % i + dest_file + ".png"))
        #                            plt.show()
        plt.close()
        TF_out_base = TafelDir.joinpath(dest_file + ".xlsx")
        # TF_out = ileOperations.CompareHashDFexport(TFxy, TF_out_base)
        rTFxy.to_excel(TF_out_base)
        # index_info_ORR_TAFEL = {
        #     "PAR_file": PAR_file,
        #     "DestFile": TF_out,
        #     "Type_output": "ORR_Jkin_calc_Tafel",
        #     "Type_exp": "ORR",
        # }

    def calculations(O2_join):
        sweep_col = [
            i
            for i in O2_join.columns
            if "SWEEP" in i.upper() and not "RING" in i.upper()
        ]
        EvRHE = [
            i
            for i in O2_join.columns
            if "E_APPV_RHE" in i.upper() and not "RING" in i.upper()
        ][0]
        ### === TAFEL EQUATION AND PLOTTING ###
        #                TFxy = O2_join.loc[(E_half[EvRHE].mean()-0.20 < O2_join[EvRHE]) & (O2_join[EvRHE] < E_half[EvRHE].mean()+0.30) & (O2_join['Sweep_Type'] == 'cathodic')]
        TFxy = O2_join.loc[
            (0.55 < O2_join[EvRHE])
            & (O2_join[EvRHE] < 0.81)
            & (O2_join["Jkin_min"] > 0)
        ]
        #        [EvRHE, 'Jkin_min', 'Jkin_max', 'Jcorr']
        #                TFxy.plot(x=EvRHE,y='Jkin_min',xlim=(0.55,0.8),logy=True)
        #                 'log10_Jkin_max' : np.log10(TFxy['Jkin_max'])
        Tafel_ovv, TFlst = pd.DataFrame([]), []
        TafelDir = []
        try:
            TFxy = TFxy[[EvRHE, "Jkin_min", "Jkin_max", "Jcorr"]].assign(
                **{
                    "log10_Jkin_min": np.log10(TFxy["Jkin_min"]),
                    "E_overp": 1.23 - TFxy[EvRHE],
                }
            )
            TFxy["log10_Jkin_min"].dropna(inplace=True)
            TFxy["log_Jkinm_2ndDiff"] = TFxy["Jkin_min"].diff().diff().values
            #                    F,D,nu,C,A,Rc,Trt = 96485, 1.5E-05, 0.010, 1.1E-06,SA_disk,8.314,298
            #                    Tafel_fit_min = linregress(TFxy['E_overp'].values,TFxy['log10_Jkin_min'].values)
            w = 60
            rTFxy = TFxy.rolling(15).mean()
            for i in range(0, len(rTFxy) - w, w):
                rTFxy.dropna(inplace=True)
                tflogJx, tfEy = (
                    rTFxy.iloc[i : i + w]["log10_Jkin_min"].values,
                    rTFxy.iloc[i : i + w][EvRHE].values,
                )
                TFfit = linregress(tflogJx, tfEy)
                #                        print(i,rTFxy.iloc[i][EvRHE],rTFxy.iloc[i+w][EvRHE],TFfit[2])
                if np.abs(TFfit[2]) > 0.95:

                    TFll = "log10_TAFEL_%.3fVrhe_%.0f" % (
                        rTFxy.iloc[i][EvRHE],
                        TFfit[0] * -1000,
                    )
                    rTFxy = rTFxy.assign(**{TFll: (rTFxy[EvRHE] - TFfit[1]) / TFfit[0]})

                    TFxy = TFxy.assign(
                        **{"log10_TAFEL_%s" % i: (TFxy[EvRHE] - TFfit[1]) / TFfit[0]}
                    )

                    #                            print(' TAFEL: ',rTFxy.iloc[i][EvRHE],rTFxy.iloc[i+w][EvRHE],TFfit[2])
                    TFlst.append(
                        {
                            "index": i,
                            "E_low": tfEy.min(),
                            "E_high": tfEy.max(),
                            "TS": TFfit[0] * -1000,
                            "TSb": TFfit[1],
                            "TSerr": TFfit[2],
                        }
                    )
                else:
                    TFlst.append(
                        {
                            "index": i,
                            "E_low": 0,
                            "E_high": 0,
                            "TS": 0,
                            "TSb": 0,
                            "TSerr": 0,
                        }
                    )
            Tafel_ovv = pd.DataFrame(TFlst)
        #                    TFxy.to_excel(TF_out)
        #                    print('TAFEL SLOPE: %.2f mV/dec,\n OUTPUT: %s'%(Tafel_ovv.iloc[Tafel_ovv['TF_E_high'].idxmax()]['TS'],TF_out))
        except Exception as e:
            print("Jkin calculation ERROR in TAFEL of ORR: {0}".format(e))
            index_info_ORR_TAFEL = {}
            TFxy = pd.DataFrame()
            Tafel_ovv = pd.concat(
                [
                    Tafel_ovv,
                    pd.DataFrame(
                        {
                            "index": 0,
                            "E_low": 0,
                            "E_high": 0,
                            "TS": 0,
                            "TSb": 0,
                            "TSerr": 0,
                        },
                        index=[0],
                    ),
                ]
            )
        ##                print(Tafel_fit[0])
        #                O2_join = O2_join.assign(**{'E_onset' : E_onset[EvRHE].min(),'Diff_lim' : Diff_lim['Jcorr'].min(),
        #                                        'E_half' : E_half[EvRHE].mean()})
        if Tafel_ovv.empty:
            Tafel_ovv = pd.concat(
                [
                    Tafel_ovv,
                    pd.DataFrame(
                        {
                            "index": 0,
                            "E_low": 0,
                            "E_high": 0,
                            "TS": 0,
                            "TSb": 0,
                            "TSerr": 0,
                        },
                        index=[0],
                    ),
                ]
            )

        TSmin, TSmax = Tafel_ovv.iloc[Tafel_ovv["TS"].idxmin()].round(
            3
        ), Tafel_ovv.iloc[Tafel_ovv["TS"].idxmax()].round(3)
        Tafel_pars_out = {
            "TSa_min": TSmin["TS"],
            "TSb_min": TSmin["TSb"],
            "TSa_max": TSmax["TS"],
            "TSb_max": TSmax["TSb"],
        }
        return Tafel_ovv, TFxy, Tafel_pars_out


class ORR_calculations:
    """
    Prepares the ORR data. Subtracts background N2 scan from ORR current:
    [ N2_jcorr, N2_Icorr ]
    """

    E_slices = {"mid": (0.25, 0.5), "high": (0.85, 0.95)}
    N2_corr_E_slice = {"mid": (0.35, 0.65)}
    mA = 1000

    N2_testcols = ["jmAcm-2", "N2_jcorr"]

    globals()["EvRHE"] = "E_AppV_RHE"

    def __init__(self, fit_run_arg):
        assert "ORR_scan_data" in str(type(fit_run_arg))

        self.ORR_data = fit_run_arg
        self.add_attributes()

        self.prepare_N2_BG_data()

        self.read_CV_disk_data()
        self.check_ORR_capacity()

        self.read_ring_data()

    def add_attributes(self):

        for att in self.ORR_data.__dict__.keys():
            setattr(self, att, getattr(self.ORR_data, att))

        self.DestDir = self.ovv_disk.ORR_ecexp_destdir.unique()[0]
        self._set_electrode_info()

    def _doublecheck_O2_act(self):
        assert self.O2_act.PAR_file.nunique() == 1
        assert self.O2_act.PAR_file.unique()[0] == self.PAR_file_disk
        # logger.warning(f'ORR difference in PAR_files from OVV and CreateCV:\n{PAR_file_ovv_file} vs {PAR_file}')

    def read_CV_disk_data(self):
        O2_CVs, O2_action = create_CVs(self.ovv_disk)
        self.O2_disk_raw_data, self.O2_disk_raw_actions = O2_CVs, O2_action

        if not O2_CVs.empty:
            _logger.info(f"Starting ORR for {self.DestDir}/{self.PAR_file_disk.stem}")
            O2_act = self.O2_disk_raw_data.query(
                '(Gas == "O2") & (Type_action == "Cyclic Voltammetry (Multiple Cycles)") & (Scanrate == 0.01) & (SampleID != "Pt_ring")'
            )
            if not O2_act.empty:
                self.O2_act = O2_act
                self.O2_act_seggrp = self.O2_act.groupby("Segment #")
                self._doublecheck_O2_act()
        else:
            _logger.warning(
                f"Not starting ORR -> ERROR empty for {self.DestDir}/{self.PAR_file_disk.stem}"
            )

    def add_rpm_list(self):
        # O2_segs = self.O2_act['Segment #'].unique()
        # O2_act_seggrp = self.O2_act.groupby('Segment #')
        _longsegs = [n for n, gr in self.O2_act_seggrp if len(gr) > 1000]
        _rpm_ref_lists = {
            5: [0, 200, 400, 900, 1500],
            3: [0, 900, 1500],
            2: [1500, 1500],
            4: [0, 200, 400, 900],
            6: [0, 200, 400, 900, 1500, 2500],
        }
        _rpm_list = _rpm_ref_lists.get(len(_longsegs), [])
        self.rpm_list = _rpm_list
        _logger.debug("Jkin calculation rpm list used: {0}".format(_rpm_list))

    def _set_electrode_info(self):
        # TODO move out and input as kwargs
        # electrode_name, collN, SA_disk, SA_ring
        for k, val in self.kwargs["WE_SA_collection_eff"]("PINE").items():
            setattr(self, k, val)

    def read_ring_data(self):
        O2_chrono_Ring, O2_chrono_actions = ORR_read_Ring_file(
            self.ovv_all, self.ovv_disk
        )
        self.O2_ring_raw_data, self.O2_ring_raw_actions = (
            O2_chrono_Ring,
            O2_chrono_actions,
        )

        if not self.O2_ring_raw_data.empty:
            self.PAR_file_ring = O2_chrono_Ring.PAR_file.unique()[0]
        if (
            "DAC Control" not in self.O2_disk_raw_actions.Type_action.unique()
            and "DAC Control" in O2_chrono_actions.Type_action.unique()
            and sum(self.O2_act.RPM_DAC.unique()) == 0
        ):
            O2_chrono_Ring[["Segment #", "RPM_DAC"]]
            _seg_rpm_cols = ["Segment #", "RPM_DAC"]
            _seg_rpms_chrono = set(
                [(i[0], i[1]) for i in O2_chrono_Ring[_seg_rpm_cols].values]
            )
            _seg_rpm_DF = pd.DataFrame(_seg_rpms_chrono, columns=_seg_rpm_cols)
            O2_act = pd.merge(
                self.O2_act.drop(columns="RPM_DAC"),
                _seg_rpm_DF,
                on=_seg_rpm_cols[0],
                how="left",
            )
            self.O2_act = O2_act
            _logger.warning(f"Jkin calculation, RPM_DAC used from Ring self.fstem")
            # TODO Match Segments with RPM_DACs of Ring and Merge with Disk

    def prepare_N2_BG_data(self):
        if not self.N2_BG_data.empty and not self.N2_option_ovv.empty:

            self.N2_option_ovv = self.N2_option_ovv.assign(
                **{
                    "N2_PAR_date_diff_seconds_abs": self.N2_option_ovv.N2_PAR_date_diff_seconds.abs()
                }
            )

            self.N2_option_ovv = self.N2_option_ovv.sort_values(
                "N2_PAR_date_diff_seconds_abs"
            )
            self.N2_option_ovv_grp = self.N2_option_ovv.groupby("PAR_file")

            if "Segm" in self.N2_BG_data.columns:
                self.N2_BG_data.Segm = self.N2_BG_data.Segm.fillna(1)
            else:
                self.N2_BG_data = self.N2_BG_data.assign(**{"Segm": 1})
            self.N2_BG_data_grp = self.N2_BG_data.groupby("PAR_file")
            self._correct_N2_BG_data()
            # self.N2_option_ovv.loc[self.N2_option_ovv.N2_PAR_date_diff_seconds.abs().idxmin()]
            # try:
            #     N2_bg_local = pd.read_excel(N2_bg_file_used, index_col=[0])
            # except Exception as e:
            #     logger.error('Jkin calculation N2 file can not be read make fake copy with 0 current: {0}'.format(e))
            #     N2_copy = O2_act.loc[(O2_act['Segment #'] == O2_segs[0])]
            #     N2_bg_local = N2_copy.assign(**{'j A/cm2' : 0, 'I(A)' : 0, 'jmAcm-2' : 0, 'Gas' : 'N2','N2_fake_0_current' : True})
        else:
            pass

    def _correct_N2_BG_data(self):
        _res = []
        _pars = []
        for n, (E_low, E_high) in self.N2_corr_E_slice.items():
            n, (E_low, E_high)
            E_name = f"{n}_{E_low}_{E_high}"
            for pf, pgrp in self.N2_BG_data_grp:
                pf, pgrp
                pgrp = pgrp.sort_values("Elapsed Time(s)")
                _anod = pgrp.loc[pgrp.Sweep_Type == "anodic"]
                _cath = pgrp.loc[pgrp.Sweep_Type == "cathodic"]
                _anod_Eslice = _anod.loc[
                    (_anod[EvRHE] > E_low) & (_anod[EvRHE] < E_high)
                ]

                _Emerge = pd.merge_asof(
                    _anod.sort_values(EvRHE),
                    _cath.sort_values(EvRHE),
                    on=EvRHE,
                    suffixes=["_anod", "_cath"],
                )
                _Emerge = _Emerge.assign(
                    **{
                        "jmean_swp": _Emerge[["jmAcm-2_anod", "jmAcm-2_cath"]].mean(
                            axis=1
                        )
                    }
                )
                _Eslice = pgrp.loc[(pgrp[EvRHE] > E_low) & (pgrp[EvRHE] < E_high)]

                _Emsl = _Emerge.loc[
                    (_Emerge[EvRHE] > E_low) & (_Emerge[EvRHE] < E_high)
                ]
                _Emsl = _Emsl.assign(
                    **{
                        "jmean_swp": _Emsl[["jmAcm-2_anod", "jmAcm-2_cath"]].mean(
                            axis=1
                        )
                    }
                )
                # _Emsl.plot(x=EvRHE,y=['jmAcm-2_anod','jmAcm-2_cath', 'jdiff_swp'])
                _linres = linregress(
                    _Emsl[EvRHE].to_numpy(), _Emsl.jmean_swp.to_numpy()
                )
                _linres_I = linregress(
                    _Emsl[EvRHE].to_numpy(),
                    _Emsl.jmean_swp.to_numpy() * 1e-3 * self.Disk_cm2,
                )
                # _linres_anod = linregress(_Emsl[EvRHE].to_numpy(), _Emsl['jmAcm-2_anod'].to_numpy())
                _linfit_b0 = pgrp[EvRHE] * _linres.slope + 0 * _linres.intercept
                pgrp = pgrp.assign(
                    **{
                        "N2_linfit_b0": _linfit_b0,
                        "N2_linfit_b0_jcorr": pgrp["jmAcm-2"] - _linfit_b0,
                    }
                )
                _linfit_b_E_slice = {"high": (0.85, 0.95)}
                _pgrp_highE_slice = pgrp.loc[
                    (pgrp[EvRHE] < np.max(_linfit_b_E_slice.get("high")))
                    & (pgrp[EvRHE] > np.min(_linfit_b_E_slice.get("high")))
                ]
                _linfit_b_val = _pgrp_highE_slice["N2_linfit_b0_jcorr"].mean()
                _linfit = pgrp[EvRHE] * _linres.slope + _linfit_b_val
                # pgrp.plot(x=EvRHE,y=['jmAcm-2','N2_linfit_b0_jcorr'])
                # (_Emerge[EvRHE]*_linres.slope  +_linres.intercept)
                _title = f"{pf.stem}, corr: {_linres.slope:.2f} + {_linfit_b_val:.2f}"
                # _Emerge = _Emerge.assign(**{'jlin' :  _linfit,
                #                             'jcorr_anod' : _Emerge['jmAcm-2_anod']-_linfit,
                #                             'jcorr_cath' : _Emerge['jmAcm-2_cath']-_linfit})
                pgrp = pgrp.assign(
                    **{
                        "N2_jcorr": pgrp["jmAcm-2"] - _linfit,
                        "N2_lincorr": _linfit,
                        "N2_Icorr": pgrp["jmAcm-2"] * 1e-3 * self.Disk_cm2
                        - (
                            pgrp[EvRHE] * _linres_I.slope
                            + _linfit_b_val * 1e-3 * self.Disk_cm2
                        ),
                    }
                )
                # pgrp.plot(x=EvRHE,y=['jmAcm-2','N2_jcorr','N2_Icorr'])
                _N2_EC_index = pgrp.N2_EC_index.unique()[0]
                _N2pars = dict(
                    zip([f"N2_BG_lin_{i}" for i in _linres._fields], _linres)
                )
                _N2pars.update(
                    {
                        "N2_BG_Eslice_name": E_name,
                        "N2_BG_Par_file": pf,
                        "N2_BG_EC_index": _N2_EC_index,
                        "N2_BG_lin_B_value": _linfit_b_val,
                    }
                )
                _pars.append(_N2pars)
                _res.append(pgrp)
                # _Emerge.plot(x=EvRHE,y=['jmAcm-2_anod','jmAcm-2_cath', 'jmean_swp','jcorr_anod',  'jcorr_cath'], title= _title)
                j_high = _anod_Eslice["jmAcm-2"].tail(20).mean()
                j_low = _anod_Eslice["jmAcm-2"].head(20).mean()
                _dJdE = (j_high - j_low) / (E_high - E_low)
        N2_BG_corr = pd.concat(_res)
        if not N2_BG_corr.empty:
            self.N2_BG_data = N2_BG_corr
            self.N2_BG_data_grp = self.N2_BG_data.groupby("PAR_file")
        if _pars:
            N2_BG_pars = pd.concat(
                [pd.DataFrame(i, index=[0]) for i in _pars], ignore_index=True
            )
            self.N2_BG_pars = N2_BG_pars
            # _anod = _anod.assign(**{'djdE' : _anod['jmAcm-2'].diff()/_anod[EvRHE].diff(),
            #                         'djdE2' : (_anod['jmAcm-2'].diff()/_anod[EvRHE].diff())/_anod[EvRHE].diff(),
            #                         'diff2j' : _anod['jmAcm-2'].diff().diff(),
            #                         'lin' : _dJdE*_anod[EvRHE],
            #                         'jcorr' :  _anod['jmAcm-2']-_dJdE*_anod[EvRHE]})
            # _cath = _cath.assign(**{'jcorr' :  _cath['jmAcm-2']-_dJdE*_cath[EvRHE]})
            # pgrp.loc[pgrp.Sweep_Type == 'anodic']
            # fig,ax = plt.subplots()
            # _anod.plot(x=EvRHE,y=['jmAcm-2','lin','jcorr'],ax=ax)
            # _cath.plot(x=EvRHE,y=['jmAcm-2','jcorr'],ax=ax)
            # ax.set_title(_title))
            # .plot(x=EvRHE,y='jmAcm-2')
            # _anod.djdE.rolling(60).mean().plot(x=EvRHE)
            # _anod.plot(x=EvRHE,y='diff2j')
            # _anod.loc[(_anod.djdE > -3) & (_anod.djdE < 3)].plot(x=EvRHE,y='djdE')

    def _check_capacity_sweeps(self, dataDF):

        all_res = []
        E_res = {}
        for Ename, Eval in self.E_slices.items():
            Elow, Ehigh = Eval
            _ORRcheck = dataDF.loc[(dataDF[EvRHE] < Ehigh) & (dataDF[EvRHE] > Elow)]
            _name = f"{Ename}_{Elow}_{Ehigh}"

            if not _ORRcheck.empty:
                _segres = {}
                # _segr = []
                for seg, sgrp in _ORRcheck.groupby("Segm"):
                    _res = {}
                    for swp, swpgrp in sgrp.groupby("Sweep_Type"):
                        jmean = swpgrp["jmAcm-2"].mean()
                        _res = {**_res, **{swp: jmean}}
                    if "anodic" and "cathodic" in _res.keys():
                        _diff = _res["anodic"] - _res["cathodic"]
                        _res = {**_res, **{"diff": _diff}}
                    # _segr.append({seg : _res})
                    _segres = {**_segres, **{seg: _res}}
                E_res = {**E_res, **{_name: _segres}}
                # all_res.append({_name : _segr})

        _E_results = [
            {"name": f"{k}", "seg": f"{k1}_{k3}", "val": v2["diff"]}
            for k, val in E_res.items()
            for k1, v2 in val.items()
            for k3, v3 in v2.items()
            if "diff" in k3
        ]
        _check_cap = pd.DataFrame(_E_results)
        return _check_cap

    def check_ORR_capacity(self):
        ORR_check_cap = self._check_capacity_sweeps(
            self.O2_act.loc[(self.O2_act.Segm > 1)]
        )
        self.ORR_check_cap = ORR_check_cap

    def check_N2_capacity(self):
        _res = []
        for n, N2r in self.N2_option_ovv.iterrows():
            print(N2r.PAR_file)
            _N2_data = self.N2_BG_data_grp.get_group(N2r.PAR_file)
            _N2_check_cap = self._check_capacity_sweeps(_N2_data)
            _N2_check_cap = _N2_check_cap.assign(**{"PAR_file": N2r.PAR_file})
            _res.append(_N2_check_cap)
        N2_check_cap = pd.concat(_res)
        self.N2_check_cap = N2_check_cap

    def compare_cap(self):

        if hasattr(self, "ORR_check_cap") and hasattr(self, "N2_check_cap"):
            _ORR_cap_agg = self.ORR_check_cap.groupby("name").agg("val").mean()
            _N2_diff_lst = [
                gr.assign(
                    **{
                        "val_diff": gr.val - _ORR_cap_agg.loc[n],
                        "O2_val": _ORR_cap_agg.loc[n],
                    }
                )
                for n, gr in self.N2_check_cap.groupby("name")
            ]
            N2_cap_diff = pd.concat(_N2_diff_lst, ignore_index=True).sort_values(
                ["name", "val_diff"]
            )

    def __iter__(self):

        for N2_cols in self.N2_testcols:
            for N2_pf, N2_BG in self.N2_BG_data_grp:
                N2_cols
                N2_pf, N2_BG
                # O2_act_slice = O2_act.loc[(O2_act['Segment #'] == seg)]
                for seg, O2_disk_seg in self.O2_act_seggrp:
                    seg, O2_disk_seg

                    yield self, O2_disk_seg, N2_BG, N2_cols


class ORR_KL_loop:
    """
    This class takes the ORR_calculation class instance and run the loop for K-L calculations over it.
    That is loop over each N2_BG file, each type of N2 current ['raw' or 'corrected']
    and finnaly over the ORR scan segments
    """

    def __init__(self, ORR_calc, run_loop=True):
        self.ORR_calc = ORR_calc
        self.fstem = self.ORR_calc.PAR_file_disk.stem
        if run_loop:
            self.calc_loop()

    def calc_loop(self):

        for N2_pf, N2_BG in self.ORR_calc.N2_BG_data_grp:
            N2_pf, N2_BG
            for N2_cols in self.ORR_calc.N2_testcols:
                N2_cols
                ORR_ops_col = []
                if hasattr(self, "ORR_ops_col"):
                    del self.ORR_ops_col
                for seg, O2_disk_seg in self.ORR_calc.O2_act_seggrp:
                    seg, O2_disk_seg
                    try:
                        ORR_ops = ORR_operations(
                            self.ORR_calc, O2_disk_seg, N2_BG, N2_cols
                        )
                        self.ORR_dest_dir_file = ORR_ops.ORR_dest_dir_file
                    except Exception as e:
                        _logger.error("ORR_KL_loop error :======= {e}")
                    if hasattr(ORR_ops, "ORR_KL_data_rpm"):
                        ORR_ops_col.append(ORR_ops)
                if ORR_ops_col:
                    self.ORR_ops_col = ORR_ops_col
                    self.KL_calc()

    # def collect_rpm_data(self):
    #      ORR_KL_data_rpm = ORR_select_KL_data(O2_CalcOut)
    #      _KL_data_all_rpms.append(ORR_KL_data_rpm)

    def KL_calc(self):
        _dt_now = datetime.now()
        # ORR_dest_dir_file = self.ORR_ops_col[-1].ORR_dest_dir_file
        try:
            # _ops_col_KL = []
            # for i in self.ORR_ops_col:
            # if hasattr(i,'ORR_KL_data_rpm'):
            KL_data = pd.concat([i.ORR_KL_data_rpm for i in self.ORR_ops_col])
            KL_data = KL_data.assign(**{"ORR_KL_Analysis_date_dt": _dt_now})
        except Exception as e:
            _logger.error("ORR_KL_loop :======= {0}".format(e))
        #    KLout = KLout.loc[(KLout['Sweep_Type'].str.contains("cathodic|anodic")), :]
        try:
            #            KLdest = 'KL_'+os.path.basename(O2_act.loc[O2_act['hash'] == h1,'File'].unique()[0]).split('.')[0]
            KL_fit = KL_plots(KL_data, self.ORR_dest_dir_file)
        #        KL_plots(KL_data, ORR_dest_dir_file)
        except Exception as e:
            _logger.error("ORR_KL_loop No KL FITTING plot:======= {0}".format(e))
            KL_fit = ["NA", "N"]

        ORR_pars_all_rpms = pd.concat([i.ORR_pars_out for i in self.ORR_ops_col])
        ORR_pars_all_rpms = ORR_pars_all_rpms.assign(
            **{"ORR_Analysis_date_dt": _dt_now}
        )

        ORR_pars_target_base = self.ORR_dest_dir_file.joinpath(f"ORR_pars_{self.fstem}")
        ORR_pars_target = FileOperations.CompareHashDFexport(
            ORR_pars_all_rpms, ORR_pars_target_base
        )
        _logger.info("Jkin calculation ORR Pars EXCEL DEST:{ORR_pars_target}")


class ORR_operations:

    # usecols_N2_correction = 'jmAcm-2'
    sweep_col = "Sweep_Type"
    E_col = "E_AppV_RHE"
    # set_jcorr_col = 'Jcorr_minus_factor'
    set_jcorr_col = ["Jcorr_minus_factor", "Jcorr_minus_lin", "Jcorr_raw"][-1]
    KL_ring_col = "I(A)_ring"
    KL_disk_col = ["Icorr_raw", "Icorr_raw_orig"]
    KL_select_data_cols = {"KL_I_Disk": "Icorr_raw_orig", "KL_I_Ring": "I(A)_ring"}

    def __init__(self, _ORR_obj, O2_disk_seg, N2_BG_scan, N2_col="jmAcm-2"):
        self.ORR_calc = _ORR_obj
        self.O2_disk_seg = O2_disk_seg
        self.N2_BG_scan = N2_BG_scan

        self.usecols_N2_correction = N2_col
        try:
            self.N2_BG_extract_meta_pars()
            self.merge_O2_with_N2_BG()
            self.add_mean_Jcorr_col()
            self.add_destfile_dirs()
            self.ORR_disk_calc_pars()
            self.read_ring_chrono()
            self.plot_RRDE()
            self.export_data()
        except Exception as e:
            _logger.warning(f"{repr(self)} error: {e}")
            pass

    def __repr__(self):
        _parent = self.ORR_calc.PAR_file_disk.parent.name
        _fstem = self.ORR_calc.PAR_file_disk.stem
        _RPM_DAC = int(self.O2_disk_seg["RPM_DAC"].unique()[0])
        _segment = int(self.O2_disk_seg["Segm"].unique()[0])
        _N2_pfstem = self.N2_BG_scan.PAR_file.unique()[0].stem
        _repr = f"{_parent}/{_fstem} seg: {_segment}, rpm: {_RPM_DAC}, N2_BG, {_N2_pfstem}, {self.usecols_N2_correction}"
        return _repr

    def add_destfile_dirs(self):
        N2_EC_index = self.N2_BG_scan.N2_EC_index.unique()[0]
        _N2_corr_col = self.usecols_N2_correction.split("-")[0]
        ORR_ecexp_destdir = self.ORR_calc.ovv_disk.ORR_ecexp_destdir.unique()[0]
        ORR_dest_dir_file = ORR_ecexp_destdir.parent.joinpath(
            ORR_ecexp_destdir.name + f"_{N2_EC_index}_{_N2_corr_col}"
        )
        ORR_dest_dir_file.mkdir(exist_ok=True, parents=True)
        self.ORR_dest_dir_file = ORR_dest_dir_file
        self.fstem = self.ORR_calc.PAR_file_disk.stem
        self.RPM_DAC = int(self.O2_disk_seg["RPM_DAC"].unique()[0])
        self.segment = int(self.O2_disk_seg["Segm"].unique()[0])
        # self.O2_disk_seg['RPM_n']
        self.dest_file = f"{self.fstem}_{self.segment}_{self.RPM_DAC}"

    def N2_BG_extract_meta_pars(self):
        N2_bg_local = self.N2_BG_scan
        _uniqcols = [i for i in N2_bg_local.columns if N2_bg_local[i].nunique() == 1]
        _N2_meta = {
            i
            if i.startswith("N2_BG_")
            else f"N2_BG_{i[3:]}"
            if i.startswith("N2_")
            else f"N2_BG_{i}": N2_bg_local[i].unique()[0]
            for i in _uniqcols
        }
        self.N2_meta = _N2_meta

    def merge_O2_with_N2_BG(self):
        O2_act_slice = self.O2_disk_seg
        N2_bg_local = self.N2_BG_scan

        mA = self.ORR_calc.mA
        O2_act_slice = O2_act_slice.assign(
            **{
                "orr_extra_E_AppV_RHE_full_prec": O2_act_slice[EvRHE],
                "Analysis_date": datetime.now(),
            }
        )
        O2_act_slice[EvRHE] = O2_act_slice[EvRHE].round(7)
        N2_bg_local = N2_bg_local.assign(
            **{"orr_extra_E_AppV_RHE_full_prec": N2_bg_local[EvRHE]}
        )
        N2_bg_local[EvRHE] = N2_bg_local[EvRHE].round(7)
        _N2grps = N2_bg_local.groupby("Sweep_Type")
        _O2grps = O2_act_slice.groupby("Sweep_Type")
        if not all([len(gr) == 1000 for n, gr in _N2grps]):
            _logger.warning("Jkin calculation merge_O2_with_N2_BG not len(1000) for N2")
        if not all([len(gr) == 1000 for n, gr in _O2grps]):
            _logger.warning("Jkin calculation merge_O2_with_N2_BG not len(1000) for O2")

        _N2cathan = pd.merge_asof(
            _N2grps.get_group("cathodic").sort_values(EvRHE),
            _N2grps.get_group("anodic").sort_values(EvRHE),
            on=EvRHE,
            direction="nearest",
            suffixes=["_cath", "_an"],
        )
        _N2cathan = _N2cathan.assign(
            **{"jmA_mean": _N2cathan[["jmAcm-2_an", "jmAcm-2_cath"]].mean(axis=1)}
        )
        _N2cathan = _N2cathan.assign(
            **{
                "jmA_mean_diff": (
                    _N2cathan["jmA_mean"].rolling(31).mean().diff() * 1e3
                ),
                "jmA_mean_diff_diff": (
                    _N2cathan["jmA_mean"]
                    .rolling(31)
                    .mean()
                    .diff()
                    .rolling(31)
                    .mean()
                    .diff()
                    .rolling(31)
                    .mean()
                )
                * 1e5,
            }
        )
        #    _N2cathan.plot(x=EvRHE,y=['jmAcm-2_an','jmAcm-2_cath','jmA_mean'])
        #    _N2cathan.plot(x=EvRHE,y=['jmA_mean_diff','jmA_mean_diff_diff','jmA_mean'],ylim=[-3,3])
        _N2cathan_slc = _N2cathan.query("(E_AppV_RHE > 0.3) & (E_AppV_RHE < 0.7)")
        _linp = linregress(
            _N2cathan_slc[EvRHE].to_numpy(), _N2cathan_slc["jmA_mean"].to_numpy()
        )
        _N2cathan = _N2cathan.assign(
            **{
                "lin_jmA_mean": linear(_N2cathan[EvRHE].to_numpy(), *_linp),
                "jmA_mean_lincorr": _N2cathan.jmA_mean
                - linear(_N2cathan[EvRHE].to_numpy(), *_linp),
            }
        )
        N2_bg_local = N2_bg_local.assign(
            **{
                "lin_jmA_linear_N2": linear(N2_bg_local[EvRHE].to_numpy(), *_linp),
                "lin_jmA_lincorr_N2": N2_bg_local["jmAcm-2"]
                - linear(N2_bg_local[EvRHE].to_numpy(), *_linp),
            }
        )
        #    _N2cathan.plot(x=EvRHE,y=['jmA_mean_diff','jmA_mean_diff_diff','jmA_mean','lin_jmA_mean'],ylim=[-3,3])
        #    N2_bg_local.plot(x=EvRHE,y=['jmAcm-2','lin_jmA_linear','lin_jmA_lincorr'])
        _N2grps = N2_bg_local.groupby("Sweep_Type")
        _N2_merge_cols = set(
            [
                EvRHE,
                "jmAcm-2",
                "I(A)",
                "lin_jmA_lincorr_N2",
                "N2_Icorr",
                "PAR_file",
                "Scan Rate (V/s)",
            ]
            + [self.usecols_N2_correction]
        )
        if not self.usecols_N2_correction in _N2_merge_cols:
            self.usecols_N2_correction = "jmAcm-2"
        _N2pars = dict(zip([f"N2_BG_lin_{i}" for i in _linp._fields], _linp))
        _N2pars.update(
            {
                "N2_BG_usecols_for_BG_correction": self.usecols_N2_correction,
                "N2_BG_Par_file": N2_bg_local.PAR_file.unique()[0],
            }
        )
        _O2N2_out = []
        for swp, _O2gr in _O2grps:
            _N2gr = _N2grps.get_group(swp)
            _N2min, _N2max = _N2gr[EvRHE].min(), _N2gr[EvRHE].max()
            _O2min, _O2max = _O2gr[EvRHE].min(), _O2gr[EvRHE].max()
            _min, _max = np.max([_N2min, _O2min]), np.min([_N2max, _O2max])

            _O2N2 = pd.merge_asof(
                _O2gr.query("(E_AppV_RHE > @_min) & (E_AppV_RHE < @_max)").sort_values(
                    EvRHE
                ),
                _N2gr.query("(E_AppV_RHE > @_min) & (E_AppV_RHE < @_max)")[
                    _N2_merge_cols
                ].sort_values(EvRHE),
                on=EvRHE,
                direction="nearest",
                suffixes=["", "_N2"],
            )
            _O2N2.columns.duplicated()
            if self.usecols_N2_correction in _O2gr:
                self.usecols_N2_correction_col = f"{self.usecols_N2_correction}_N2"
            else:
                self.usecols_N2_correction_col = self.usecols_N2_correction
            _N2pars.update(
                {"N2_BG_usecols_for_BG_correction": self.usecols_N2_correction_col}
            )

            _O2N2 = _O2N2.assign(
                **{
                    "Jcorr_raw": _O2N2["jmAcm-2"]
                    - _O2N2[f"{self.usecols_N2_correction_col}"],
                    "Jcorr_N2lincorr": _O2N2["jmAcm-2"] - _O2N2["lin_jmA_lincorr_N2"],
                    "Icorr_raw": _O2N2["I(A)"] - _O2N2["N2_Icorr"],
                    # (_O2N2['jmAcm-2'] - _O2N2[f'{self.usecols_N2_correction_col}'])*1E3*0.2376,
                    "Icorr_raw_orig": _O2N2["I(A)"] - _O2N2["I(A)_N2"],
                }
            )
            # FIXME original Icorr_raw used for KL I data, replaced with j now

            #                    _Jcorr_spl = UnivariateSpline(_O2N2[EvRHE],_O2N2['Jcorr_raw'])
            _O2N2 = _O2N2.assign(
                **{
                    "Jcorr_O2diff": (
                        _O2N2["Jcorr_raw"].rolling(31).mean().diff() * 1e3
                    ),
                    "Jcorr_O2diffdiff": (
                        _O2N2["Jcorr_raw"]
                        .rolling(31)
                        .mean()
                        .diff()
                        .rolling(31)
                        .mean()
                        .diff()
                        .rolling(31)
                        .mean()
                    )
                    * 1e5,
                }
            )
            _O2N2_out.append(_O2N2)
        O2_join = pd.concat(_O2N2_out)

        self.O2_disk_N2BG = O2_join
        self.N2_calc_pars = {**self.N2_meta, **_N2pars}

    def add_mean_Jcorr_col(self):
        #    mean_cols = ['Jcorr_raw','Icorr_raw','jmAcm-2','I(A)_N2','I(A)','j A/cm2' ]):
        #    O2DF = O2_join_raw
        sweep_col = self.sweep_col
        E_col = self.E_col
        O2DF_res = self.O2_disk_N2BG.reset_index()
        cath = O2DF_res.groupby(sweep_col).get_group("cathodic").sort_values(E_col)
        anod = O2DF_res.groupby(sweep_col).get_group("anodic").sort_values(E_col)
        _n1cols = [i for i in O2DF_res.columns if O2DF_res[i].nunique() <= 2]
        _nXcols = [i for i in O2DF_res.columns if O2DF_res[i].nunique() > 2]
        _mean_cols_used = [
            i
            for i in O2DF_res.columns
            if (i.upper().startswith("J") or i.upper().startswith("I"))
            and i not in _n1cols
        ]

        swp_merge = pd.merge_asof(
            cath,
            anod[[i for i in anod.columns if i not in _n1cols]],
            on=E_col,
            suffixes=["_cath", "_anod"],
        )
        _mean_dct = {
            i: swp_merge[[f"{i}_cath", f"{i}_anod"]].mean(axis=1)
            for i in _mean_cols_used
        }
        swp_merge = swp_merge.assign(**{**{"Sweep_Type": "mean"}, **_mean_dct})
        swp_merge = swp_merge.assign(
            **{
                "Jcorr_O2diff": (
                    swp_merge["Jcorr_raw"].rolling(31).mean().diff() * 1e3
                ),
                "Jcorr_O2diffdiff": (
                    swp_merge["Jcorr_raw"]
                    .rolling(31)
                    .mean()
                    .diff()
                    .rolling(31)
                    .mean()
                    .diff()
                    .rolling(31)
                    .mean()
                )
                * 1e5,
            }
        )
        O2DF_mean = pd.concat([O2DF_res, swp_merge], ignore_index=True)
        self.O2_disk_N2BG_mean = O2DF_mean
        # return O2DF_mean
        # === TEST PLOTTING
        #                t1.plot(x='E_AppV_RHE',y=['J_O2_diff','Jcorr','J_N2_scan'],xlim=(0.9,1.05),ylim=(-1E0,1E0))
        #                fig,ax = plt.subplots()
        #                O2_join.plot(x='E_AppV_RHE',y=['jmAcm-2','Jcorr','J_N2_scan'],xlim=(0,1.1),ylim=(-6E0,5E0),ax=ax)
        #                plt.savefig(ORR_dest_dir.joinpath(dest_file+'_test.png'),dpi=100,bbox_inches='tight')
        #                plt.close()
        #                Jderiv_Diffs = O2_join.loc[(O2_join['J_O2_diff'] == O2_join['J_O2_diff'].max()) | (O2_join['J_O2_diff'] == O2_join['J_O2_diff'].min()) & (O2_join[EvRHE] < 1.0),:]
        #                N2_scan.plot(x='E(V)',y='I(A)',xlim=(-1,1))
        ### ==== Calculation of Jdiff, then Jkin IMPORTANT PARAMETERS ====== ###

    #            'Jcorr_minus_factor', 'Jcorr_minus_lin'
    def ORR_disk_calc_pars(self):
        # Calculation of Jdiff, then Jkin IMPORTANT PARAMETERS !! #
        # PAR_file, ORR_dest_dir_file, dest_file, O2_join, **kwargs):
        #    _meta_O2_join = O2_join[[i for i in O2_join.columns if O2_join[i].nunique() == 1]].iloc[0].to_dict()
        PAR_file = self.ORR_calc.PAR_file_disk
        try:
            _ORR_disk_calc, _ORR_disk_Tafel = [], []
            _swpgr_Jkin_join = []
            for swp, swgrp in self.O2_disk_N2BG_mean.groupby("Sweep_Type"):
                swp, swgrp
                #            swp,swgrp = 'anodic', O2_join.groupby('Sweep_Type').get_group('anodic')
                _pars = (
                    swgrp[[i for i in swgrp.columns if swgrp[i].nunique() == 1]]
                    .iloc[0]
                    .to_dict()
                )
                # set_jcorr_col = kwargs.get('set_jcorr_col', 'Jcorr_minus_lin')
                swgrp, _lincorr_pars = determine_Jcorr_zero_correction(
                    swgrp, set_jcorr_col=self.set_jcorr_col
                )
                #            swgrp.plot(x=EvRHE,y=['Jcorr_minus_lin','Jcorr_minus_factor','Jcorr_raw'])
                #            swgrp,_horizontal_res = find_horizontal_stretch(swgrp)
                # _O2N2 = _O2N2.assign(**{'Jcorr' : _O2N2.Jcorr_raw - _Jcorr_factor})
                _pars.update(_lincorr_pars)
                #            _O2corr_pars_out.append(_O2corr_pars)
                #            _pars = {'PAR_file': PAR_file, 'RPM_n': rpm_n, 'Segment': seg,'Sweep_Type' : swp}
                #            _pars.update(**_meta_O2_join)
                E_onset = ORR_determine_E_onset(swgrp, PAR_file)
                E_limit_lowdiff = 0.03
                Diff_lim = ORR_determine_Jdiff_lim(swgrp, PAR_file, E_limit_lowdiff)
                RHE_OCP_0 = swgrp["RHE_OCP"].unique()[0] * 1000
                #            E_lowdiff = 0.03
                swgrp = swgrp.assign(
                    **{
                        "Jkin_max": (
                            -1
                            * (Diff_lim.max * swgrp["Jcorr"])
                            / (Diff_lim.max - swgrp["Jcorr"])
                        ),
                        "Jkin_min": (
                            -1
                            * (Diff_lim.min * swgrp["Jcorr"])
                            / (Diff_lim.min - swgrp["Jcorr"])
                        ),
                    }
                )
                E_half = swgrp.loc[
                    np.isclose(swgrp["Jcorr"], Diff_lim.min / 2, rtol=0.005)
                    & (swgrp[EvRHE] < E_onset[EvRHE].mean()),
                    :,
                ]
                #                E_half = swgrp.loc[(swgrp['Jcorr'] < Diff_lim['Jcorr'].min() * 0.4980) & (
                #                (swgrp['Jcorr'] > Diff_lim['Jcorr'].min() * 0.5540)) & (swgrp[EvRHE] < E_onset[EvRHE].mean()), :]
                if E_half.empty:
                    #                _logger.warning(f'Jkin calculation E_half empty! expanding...{PAR_file},{swp}')
                    E_half = swgrp.loc[
                        np.isclose(swgrp["Jcorr"], Diff_lim.min / 2, rtol=0.05)
                        & (swgrp[EvRHE] < E_onset[EvRHE].mean()),
                        :,
                    ]

                _pars.update(
                    {
                        "ORR_E_onset": E_onset[EvRHE].mean(),
                        "ORR_E_half": E_half[EvRHE].mean(),
                        "ORR_J_diff_lim": Diff_lim.min,
                        "ORR_J_diff_lim_max": Diff_lim.min,
                        "RHE_OCP_mV": RHE_OCP_0,
                    }
                )
                #                    i=0 i += 1     while E_half.empty:
                #                        E_half = O2_join.loc[(O2_join['Jcorr'] < Diff_lim['Jcorr'].mean()*(0.4980-0.001*i)) & ((O2_join['Jcorr'] > Diff_lim['Jcorr'].mean()*(0.5540+0.001*i))) & (O2_join['Sweep_Type'] == "cathodic")& (O2_join[EvRHE] < E_onset[EvRHE].mean()),:]
                for E_mV in range(500, 850, 50):
                    _EV = E_mV * 1e-3
                    for _parnm in ["Jkin_min", "Jkin_max"]:
                        _parv = swgrp.loc[
                            np.isclose(swgrp[EvRHE], _EV, atol=0.005), _parnm
                        ].mean()
                        _pars.update({f"ORR_{_parnm}_{E_mV}": _parv})

                #            Jkin_075 = swgrp.loc[(np.isclose(swgrp[EvRHE], 0.75, rtol=0.003)), :]
                #            Jkin_080 = swgrp.loc[(np.isclose(swgrp[EvRHE], 0.8, rtol=0.003)), :]
                #                Jkin_075 = O2_join.query('Sweep_Type == "cathodic"').loc[np.isclose(O2_join[EvRHE],0.75,rtol=0.003),:]
                #                & (O2_join['Jcorr'] < -0.085/mA) & (O2_join['Sweep_Type'] == "cathodic")]
                #            _pars.update(
                #             'ORR_Jkin_075': np.abs(Jkin_075['Jkin_min'].mean()), 'ORR_Jkin_080': np.abs(Jkin_080['Jkin_min'].mean())})
                _ORR_disk_calc.append(_pars)
                ### Run Tafel analysis ###
                Tafel_ovv = ORR_get_Tafel(swgrp, self.ORR_dest_dir_file, self.dest_file)
                _Tafel_pars = Tafel_ovv.assign(**_pars)
                _ORR_disk_Tafel.append(_Tafel_pars)
                _swpgr_Jkin_join.append(swgrp)

            ORR_disk_pars_Tafel = pd.concat(_ORR_disk_Tafel)
            ORR_disk_pars_Tafel = ORR_disk_pars_Tafel.assign(
                **self.N2_calc_pars
            ).drop_duplicates()
            _O2_join_Jkin = pd.concat(_swpgr_Jkin_join).sort_index()
            ORR_disk_pars = pd.DataFrame(_ORR_disk_calc)
            ORR_disk_pars = ORR_disk_pars.assign(**self.N2_calc_pars).drop_duplicates()
        except Exception as e:
            _logger.error(
                f"Jkin calculation ORR disk pars: {e}\n{self}, dest:{self.dest_file}"
            )
            ORR_disk_pars_Tafel = pd.DataFrame()
            ORR_disk_pars = pd.DataFrame()
        self.O2_join_Jkin = _O2_join_Jkin
        self.O2_disk_N2BG_mean_Jkin = self.O2_join_Jkin
        self.ORR_disk_pars_Tafel = ORR_disk_pars_Tafel
        self.ORR_disk_pars = ORR_disk_pars
        # return _O2_join_Jkin, ORR_disk_pars_Tafel

    def read_ring_chrono(self):
        try:
            Iring_Chrono, _ORR_RRDE_pars = O2_Chrono(
                self.segment,
                self.RPM_DAC,
                self.O2_disk_N2BG_mean_Jkin,
                self.ORR_calc.O2_ring_raw_data,
                self.ORR_dest_dir_file,
                plot_BipotChrono=True,
                export_excel=True,
                ring_filter=True,
            )
        #                   'Jring_050': out_Jring, 'FracH2O2_050': out_FracH2O2})
        except Exception as e:
            _logger.error(f"Jkin calculation ERROR: Iring_Chrono {e}, {self.dest_file}")
            Iring_Chrono, _ORR_RRDE_pars = pd.DataFrame(), pd.DataFrame()
        self.Iring_Chrono = Iring_Chrono
        self.ORR_RRDE_pars = _ORR_RRDE_pars

    def plot_RRDE(self):
        try:
            ORR_plot_ring_disk(
                self.O2_disk_N2BG_mean_Jkin,
                self.ORR_disk_pars,
                self.Iring_Chrono,
                self.ORR_RRDE_pars,
                self.ORR_dest_dir_file,
                self.dest_file,
            )
        except Exception as e:
            _logger.error(f"Jkin calculation ERROR: plot_RRDE {e}, {self.dest_file}")

    def export_data(self):
        # self.O2_disk_N2BG_mean_Jkin
        O2_CalcOut = pd.merge(
            self.O2_disk_N2BG_mean_Jkin,
            self.Iring_Chrono,
            on=[EvRHE, "Sweep_Type"],
            how="left",
        )
        #            O2_join.plot(x=EvRHE, y=['Jcorr']),  Iring_Chrono.plot(x=EvRHE, y=['J_ring'])
        #            O2_CalcOut.groupby('Sweep_Type').plot(x=EvRHE, y=['Jcorr','J_ring'])
        yOut = [
            "E_AppV_RHE",
            "jmAcm-2",
            "jmAcm-2_N2",
            "Jcorr",
            "Jkin_min",
            "Jkin_max",
            "J_ring",
            "Frac_H2O2",
            "n_ORR",
            "Elapsed Time(s)",
            "Sweep_Type",
            "RPM_DAC",
        ]
        #                O2_CalcOut[yOut].to_excel(ORR_dest_dir.joinpath(dest_file+'_RRDE.xlsx'))
        ORR_Excel_dest_base = self.ORR_dest_dir_file.joinpath(
            f"{self.dest_file}_RRDE.xlsx"
        )
        ORR_Excel_dest = FileOperations.CompareHashDFexport(
            O2_CalcOut[yOut], ORR_Excel_dest_base
        )
        self.ORR_data_out = O2_CalcOut
        self.ORR_data_out = self.ORR_data_out.assign(**self.N2_meta)
        O2_ParsOut = (
            pd.merge(self.ORR_disk_pars, self.ORR_RRDE_pars)
            .assign(**{"RRDE_merged_data_file": ORR_Excel_dest})
            .drop_duplicates()
        )
        self.ORR_pars_out = O2_ParsOut
        # _ORR_parslst.append(O2_ParsOut )
        self.ORR_select_KL_data()
        _logger.info(f"Jkin calculation ORR Jkin calc succes: {self.ORR_dest_dir_file}")
        # _KL_data_all_rpms.append(ORR_KL_data_rpm)

    def ORR_select_KL_data(self):
        # O2_CalcOut =
        KLrpm = pd.DataFrame([])
        try:
            KLoutRow = []
            for swp, swgrp in self.ORR_data_out.groupby("Sweep_Type"):
                _d_meta = (
                    swgrp[[i for i in swgrp.columns if swgrp[i].nunique() == 1]]
                    .iloc[0]
                    .to_dict()
                )
                for E in np.linspace(0, 1, 41):
                    _row = {}
                    _row.update(_d_meta)

                    swgrp_E_slice = swgrp.loc[(np.isclose(swgrp[EvRHE], E, atol=0.010))]
                    _row.update({f"KL_{EvRHE}": E})
                    for KLnm, KLcol in self.KL_select_data_cols.items():
                        KLcolmean = np.abs(swgrp_E_slice[KLcol]).mean()
                        _row.update({KLnm: KLcolmean})

                    # KLoutDisk = np.abs(swgrp_E_slice[self.KL_disk_col]).mean()
                    # #                        KLoutDisk = np.abs(SwpGrp.loc[(np.isclose(SwpGrp[EvRHE],E,atol=0.010)) & (SwpGrp['Sweep_Type'] == swp)]['I(A)_disk']).mean()
                    # KLoutRing = np.abs(swgrp_E_slice[self.KL_ring_col]).mean()
                    # _row.update({f'KL_{EvRHE}' : E, 'KL_I_Disk' : KLoutDisk, 'KL_I_Ring' : KLoutRing})
                    KLoutRow.append(_row)
            #                KLoutRow.append([swp, E, KLoutDisk, KLoutRing])
            #                        KLoutRow ={'SweepType' : swp,EvRHE : E,'I_Disk' : KLoutDisk,'I_Ring' : KLoutRing}
            #                        KLoutRow.update(KLoutRow)
            #                        +'_%s'%rpm_list[rpm], +'_%s'%rpm_list[rpm]
            KLrpm = pd.DataFrame(KLoutRow)
        #        _logger.info('ORR KL prep succes {0}'.format(PAR_file))
        except Exception as e:
            _logger.warning("ORR ERROR KL prep {e}")
        self.ORR_KL_data_rpm = KLrpm

        def _testplot():
            self.ORR_KL_data_rpm.groupby("Sweep_Type").plot(
                x=f"KL_{EvRHE}", y="KL_I_Disk"
            )

        # return KLrpm


#      def _test():
# #         ORR_disk_pars = ORR_disk_pars.assign(**_N2BG_pars).drop_duplicates()
# # #            (PAR_file,rpm_n, seg, ORR_dest_dir_file, dest_file, O2_join)
# # #                                'Analysis_date': datetime.now()}

# #             ### +++++++++-   ++++++++++++++++++++++    -+++++++++++ ###
# #     ### +++++++++- PLOTTING -+++++++++++ ###
# #             # https://stackoverflow.com/questions/37737538/merge-matplotlib-subplots-with-shared-x-axis
# #             #                fig, (axRing,axJ) = plt.subplots(2, sharex=True)
# #             #                SegGr.loc[:,[x]+y].to_csv(ORR_dest_dir.joinpath(Path(f).stem+'_%s'%rpm+'.csv'))
# #             #                                gr.plot(x=x,y=y,xlim=(0,1.2),ylim=(-6,0.5),title=Path(f).stem)
# #             #                                plt.savefig(ORR_dest_dir.joinath(Path(f).stem)+'_O2.png',bbox_inches='tight')
# #             #                                plt.close()
# #         ORR_plot_ring_disk( O2_join, ORR_disk_pars , Iring_Chrono, _ORR_RRDE_pars,
# #                            ORR_dest_dir_file, dest_file)
# # ### +++++++++- FINALIZING -+++++++++++ ###
# #         #                O2_out_ovv,O2_CalcOut = pd.DataFrame([]), pd.DataFrame([])

#         O2_CalcOut = pd.merge(O2_join, Iring_Chrono, on=[EvRHE,'Sweep_Type'],how='left')
# #            O2_join.plot(x=EvRHE, y=['Jcorr']),  Iring_Chrono.plot(x=EvRHE, y=['J_ring'])
# #            O2_CalcOut.groupby('Sweep_Type').plot(x=EvRHE, y=['Jcorr','J_ring'])
#         yOut = ['E_AppV_RHE', 'jmAcm-2', 'jmAcm-2_N2', 'Jcorr', 'Jkin_min', 'Jkin_max', 'J_ring',  'Frac_H2O2', 'n_ORR','Elapsed Time(s)','Sweep_Type','RPM_DAC']
#         #                O2_CalcOut[yOut].to_excel(ORR_dest_dir.joinpath(dest_file+'_RRDE.xlsx'))

#         ORR_Excel_dest_base = ORR_dest_dir_file.joinpath(f'{dest_file}_RRDE.xlsx')
#         ORR_Excel_dest = FileOperations.CompareHashDFexport(O2_CalcOut[yOut], ORR_Excel_dest_base)

#         O2_ParsOut = pd.merge(ORR_disk_pars,_ORR_RRDE_pars).assign(**{'RRDE_merged_data_file' : ORR_Excel_dest}).drop_duplicates()
#         _ORR_parslst.append(O2_ParsOut )
#         ORR_KL_data_rpm = ORR_select_KL_data(O2_CalcOut)
#         _KL_data_all_rpms.append(ORR_KL_data_rpm)


def ORR_read_Ring_file(gr_ovv, ORR_ovv_file):
    #                       PAR_file, fstem, ORR_file_PAR_date):
    #                pd.concat([E_onset,Diff_lim,E_half])
    ### +++++++++- Choosing Ring File and Adding I_ring from Chrono -+++++++++++ ###
    #            ('Disc','Ring'),('disc','ring'), ('Ch1_Disk','Ch2_Ring'), ('disk','ring'), ('Disc','Ring')
    Ring_ovv = pd.DataFrame()
    PAR_file, fstem = (
        ORR_ovv_file.PAR_file.values[0],
        ORR_ovv_file.basename.values[0],
    )

    ring_date_ovv = gr_ovv.loc[
        (
            gr_ovv.PAR_date.isin(ORR_ovv_file.PAR_date)
            | gr_ovv.PAR_date_min.isin(ORR_ovv_file.PAR_date_min)
        )
        & gr_ovv.EXP_dir.isin(ORR_ovv_file.EXP_dir)
        & (gr_ovv.Electrode == "Pt_ring")
    ]
    #    ring_date_ovv = gr_ovv.query('(PAR_date == @ORR_file_PAR_date) & (Electrode == "Pt_ring")')
    #    if ring_date_ovv.empty:
    #        ORR_f_par_date_min = f'{pd.to_datetime(ORR_file_PAR_date):%Y-%m-%dT%H:%M}'
    #        ring_date_ovv = gr_ovv.query('(PAR_date_min == @ORR_f_par_date_min) & (Electrode == "Pt_ring")')

    if ring_date_ovv.empty:
        ORR_file_PAR_date = ORR_ovv_file.PAR_date.values[0]
        gr_ovv["PAR_exp_delta_disk"] = [
            pd.to_datetime(i) - pd.to_datetime(ORR_file_PAR_date)
            for i in gr_ovv.PAR_date.values
        ]
        ring_date_ovv = gr_ovv.loc[
            (gr_ovv["PAR_exp_delta_disk"] > pd.to_timedelta(0))
            & (gr_ovv["PAR_exp_delta_disk"] < pd.to_timedelta(10, unit="s"))
        ]
    #                    ring_date_ovv = ovv.loc[ovv.Comment.str.contains('Pt|chrono') & np.isclose(ovv.PAR_date,ORR_file_PAR_date)]
    if len(ring_date_ovv) == 1:
        Ring_ovv = ring_date_ovv
        _logger.info(
            "Jkin calculation O2_chrono_Ring match found: {0} for {1}".format(
                ring_date_ovv.basename.iloc[0], fstem
            )
        )
    else:
        _logger.warning(
            "Jkin calculation O2_chrono_Ring matches error: {0} for {1}".format(
                ring_date_ovv.basename.values, fstem
            )
        )
        O2_PAR_file_upper = fstem.upper()
        O2_ring_fn_options_disC = [
            ("DISC", "RING"),
            ("CH1_DISC", "CH2_RING"),
            ("V3F_DISC", "V3_RING"),
            ("_V3.PAR", "_V3F.PAR"),
            ("_V3.PAR", "_V3F.PAR"),
        ]

        O2_ring_fn_options_disK = [
            (i[0].replace("DISC", "DISK"), i[1])
            for i in O2_ring_fn_options_disC
            if "DISC" in i[0]
        ]
        O2_ring_fn_options = O2_ring_fn_options_disC + O2_ring_fn_options_disK

        O2_file_ring_opts = []
        for _disc_opt, _ring_opt in O2_ring_fn_options:
            _ring_fn_test = [
                i
                for i in gr_ovv.PAR_file.unique()
                if Path(i).name.upper()
                == O2_PAR_file_upper.replace(_disc_opt, _ring_opt)
            ]
            if _ring_fn_test:
                O2_file_ring_opts.append(_ring_fn_test)

        O2_file_ring = list(
            set(
                [
                    i
                    for i in O2_file_ring_opts
                    if O2_PAR_file_upper not in str(i).upper()
                ]
            )
        )
        if O2_file_ring:
            #                    O2_chrono_CV = O2_CVs.loc[(O2_CVs.File.isin(O2_file_ring) | (O2_CVs.File == PAR_file))].query('SampleID == @O2_act_SampleID')
            Ring_ovv = gr_ovv.query("PAR_file == @O2_file_ring")
        else:
            _logger.warning(
                "Jkin calculation O2_chrono_Ring file options empty: {0}".format(
                    ring_date_ovv.basename.values
                )
            )
            Ring_ovv = pd.DataFrame()
    #        O2_fr_1 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('DISC', 'RING')]
    #        O2_fr_2 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('DISK', 'RING')]
    #        O2_fr_3 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('CH1_DISK', 'CH2_RING')]
    #        O2_fr_4 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('CH1_DISC', 'CH2_RING')]
    #        O2_fr_5 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('V3F_DISK', 'V3_RING')]
    #        O2_fr_6 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('_V3.PAR', '_V3F.PAR')]
    #        O2_file_ring_opts = O2_fr_1 + O2_fr_2 + O2_fr_3 + O2_fr_4 + O2_fr_5 + O2_fr_6
    #        O2_file_ring = list(set([i for i in O2_file_ring_opts if O2_PAR_file_upper not in str(i).upper()]))[0]
    #                print('Ring files {0} for PAR file: {1}'.format(O2_file_ring,PAR_file))
    if not Ring_ovv.empty:
        if Ring_ovv.basename.nunique() != 1:
            _logger.warning(
                "Jkin calculation create CV from O2_Chrono Ring multiple Ring Files in ovv!!: {0}".format(
                    Ring_ovv.basename.values
                )
            )
        _logger.info(
            "Jkin calculation used; DISK: {0},\n RING: {1}".format(
                Path(PAR_file).name, Ring_ovv.basename.iloc[0]
            )
        )
        O2_chrono_Ring, O2_chrono_actions = create_CVs(Ring_ovv)
        if O2_chrono_Ring.empty:
            O2_chrono_Ring = pd.DataFrame()
            _logger.warning(
                "Jkin calculation create CV from O2_Chrono Ring failed: {0}".format(
                    Ring_ovv.basename.iloc[0]
                )
            )
    else:
        _logger.warning(
            "Jkin calculation O2_chrono_Ring MATCH OVV empty: {0}".format(
                ring_date_ovv.basename.values
            )
        )
    return O2_chrono_Ring, O2_chrono_actions


def linear(x, a, b, *args):
    return a * x + b


def _get_merge_cols(ORR_disk_pars, _ORR_RRDE_pars):
    _l, _r = ORR_disk_pars, _ORR_RRDE_pars
    _l1 = [(i, _l[i].unique()[0]) for i in _l.columns if _l[i].nunique() <= 1]
    _r1 = [(i, _r[i].unique()[0]) for i in _r.columns if _r[i].nunique() <= 1]
    _mcls = set([i[0] for i in _l1]).intersection([i[0] for i in _r1])
    pd.merge(_l, _r, on=list(_mcls) + ["Sweep_Type"])


# def add_mean_Jcorr_col(_DF):
#    cath = _DF.groupby('Sweep_Type').get_group('cathodic')
#    anod = _DF.groupby('Sweep_Type').get_group('anodic')
#    _n1cols = [i for i in _DF.columns if _DF[i].nunique() ==1]
#    swp_merge = pd.merge_asof(cath, anod[[i for i in anod.columns if i not in _n1cols]],on=EvRHE,suffixes=['_cath','_anod'])
#    swp_merge = swp_merge.assign(**{'Sweep_Type' : 'mean','Jcorr_raw': swp_merge[['Jcorr_raw_cath', 'Jcorr_raw_anod']].mean(axis=1),
#                                    'Jcorr_raw': swp_merge[['Jcorr_raw_cath', 'Jcorr_raw_anod']].mean(axis=1)})
##                _n1cols_an = [i for i in swp_merge if swp_merge[i].nunique() ==1 if i.endswith('_anod')]
##                _n1cols_cath = [i for i in swp_merge if swp_merge[i].nunique() ==1 if i.endswith('_cath')]
##                [i for i in _n1cols_an if swp_merge[i].unique()[0] == swp_merge[i.replace('_anod','_cath')].unique()[0] ]
##                [swp_merge[i].unique() for i in _n1cols_an]
##                swp_merge.plot(x='E_AppV_RHE',y=['Jcorr_raw_cath', 'Jcorr_raw_anod','Jcorr_raw'],xlim=(0,1.1),ylim=(-3E0,1E0))
#    O2_join_Jmean = pd.concat([_DF, swp_merge],ignore_index=True)
#    return O2_join_Jmean
#


class ORR_determine_Jcorr_zero:
    def __init__(self, swgrp, set_jcorr_col="Jcorr_minus_factor"):
        self.swgrp = swgrp
        self.set_jcorr_col = set_jcorr_col

    def _find_hor_stretch(self):
        swgrp, _hor_mean, _horizontal_res, _overall_lin_pars = find_horizontal_stretch(
            self.swgrp
        )


def determine_Jcorr_zero_correction(swgrp, set_jcorr_col="Jcorr_minus_factor"):
    #    _O2N2 = swgrp
    swgrp, _hor_mean, _horizontal_res, _overall_lin_pars = find_horizontal_stretch(
        swgrp
    )
    #    swgrp.plot(EvRHE,y=['Jcorr_raw', 'Jcorr_hor'])
    #    _horizontal_res.loc[_horizontal_res.hor_lin_R < 0.4].mean()
    _lin_pars = {f"ORR_lin_{key}": val for key, val in _hor_mean.to_dict().items()}
    _lin_pars = {**_lin_pars, **_overall_lin_pars}
    #    _Jdiff_slc = swgrp.query('(E_AppV_RHE > 0.3) & (E_AppV_RHE < 0.95)')
    _J_highE_plateau = swgrp.query(
        "(E_AppV_RHE > @_hor_mean.hor_E_lower) & (E_AppV_RHE < @_hor_mean.hor_E_upper)"
    )
    #    _J_highE_check_max = swgrp.query('(E_AppV_RHE > 0.75) & (E_AppV_RHE < 1)')
    _linp = linregress(
        _J_highE_plateau[EvRHE].to_numpy(), _J_highE_plateau["Jcorr_raw"].to_numpy()
    )
    _J_highE_plateau = _J_highE_plateau.assign(
        **{"plateau_Jcorr_raw": linear(_J_highE_plateau[EvRHE].to_numpy(), *_linp)}
    )
    #    _J_highE_check_max = _J_highE_check_max.assign(**{'plateau_Jcorr_raw' : linear(_J_highE_check_max[EvRHE].to_numpy(),*_linp),
    #                                                   'hor_Jcorr_raw_diff' :  _J_highE_check_max.Jcorr_raw - _hor_popt})
    #    _J_highE_check_max.plot(EvRHE,y=['Jcorr_raw','plateau_Jcorr_raw', 'hor_Jcorr_raw_diff'])
    #    _J_highE_plateau.plot(EvRHE,y=['Jcorr_raw','plateau_Jcorr_raw'])
    swgrp = swgrp.assign(
        **{"Jcorr_minus_lin": swgrp.Jcorr_raw - linear(swgrp[EvRHE].to_numpy(), *_linp)}
    )
    _lin_pars = {
        **_lin_pars,
        **dict(zip([f"ORR_lin_{i}" for i in _linp._fields], _linp)),
    }
    _Jfactor_mean = _J_highE_plateau.Jcorr_raw.mean()

    if (
        _overall_lin_pars["hor_overall_slope"] > 3
        and _overall_lin_pars["hor_overall_rvalue"] > 0.85
    ):
        #        J_highE_slope = swgrp.query('(E_AppV_RHE > @_overall_lin_pars.get("hor_overall_E_lower")) & (E_AppV_RHE < @_overall_lin_pars.get("hor_overall_E_upper"))')
        #        J_highE_slope.plot(EvRHE,y=['Jcorr_raw'])
        _Jfactor_mean = _overall_lin_pars.get("hor_overall_J_mean", 0)

    if _hor_mean.hor_E_lower < 0.77 and _overall_lin_pars["hor_overall_J_mean"] > 0.1:
        #        or (_hor_mean.hor_Ewindow < 0.06)) :
        _Jfactor_mean = 0

    #    _J_highE_check_max = _J_highE_check_max.assign(**{'Jcorr_minus_lin' : _J_highE_check_max.Jcorr_raw - linear(_J_highE_check_max[EvRHE].to_numpy(),*_linp),
    #                           'Jcorr_minus_factor' : _J_highE_check_max.Jcorr_raw - _Jfactor_mean })
    #    if 'Jcorr_minus_factor' in set_jcorr_col:
    #        _Jfactor_mean  = _Jfactor_mean
    swgrp = swgrp.assign(
        **{
            "Jcorr_minus_lin": swgrp.Jcorr_raw
            - linear(swgrp[EvRHE].to_numpy(), *_linp),
            "Jcorr_minus_factor": swgrp.Jcorr_raw - _Jfactor_mean,
        }
    )
    if set_jcorr_col in swgrp.columns:
        swgrp = swgrp.assign(**{"Jcorr": swgrp[set_jcorr_col]})
    _lin_pars.update(
        {"ORR_lin_Jcorr_factormean": _Jfactor_mean, "ORR_set_Jcorr_col": set_jcorr_col}
    )
    return swgrp, _lin_pars


#    swgrp.plot(EvRHE,y=['Jcorr_raw', 'Jcorr_minus_lin','Jcorr'])
#     _Jdiff_slc = swgrp.query('(E_AppV_RHE > 0.3) & (E_AppV_RHE < 0.95)')
#    _E_jdiffmax = _Jdiff_slc.loc[_Jdiff_slc.Jcorr_O2diff.idxmax()][EvRHE]
#    _E_jdiffdiffmin = _Jdiff_slc.loc[_Jdiff_slc.Jcorr_O2diffdiff.idxmin()][EvRHE]
#    _Ewindow = (_E_jdiffdiffmin,_E_jdiffmax)
#    _Ew_min, _Ew_max = np.min(_Ewindow), np.max(_Ewindow)
#    _Jcorr_factor = swgrp.query('(E_AppV_RHE > @_Ew_min) & (E_AppV_RHE < @_Ew_max)').Jcorr_raw.max()


def find_horizontal_stretch(swgrp):
    def hor_y(x, b):
        a = 0
        return a * x + b

    def lin_y(x, a, b):
        return a * x + b

    _J_highE_check_max = swgrp.query("(E_AppV_RHE > 0.5)")
    #    _J_highE_check_max.plot(x=EvRHE,y='Jcorr_raw')
    _res = []
    for c in range(70, 1, -1):
        _i = c / 100
        #         _J_highE_check_max = swgrp.query('(E_AppV_RHE > @_i0) & (E_AppV_RHE < @_i1)')
        _J_close_to_0 = _J_highE_check_max.loc[
            (_J_highE_check_max["Jcorr_raw"] < c / 100)
            & (_J_highE_check_max["Jcorr_raw"] > -c / 100)
        ]
        #        if not _J_close_to_0.empty:

        if not len(_J_close_to_0) < 5:
            _i0, _i1 = (
                _J_close_to_0["E_AppV_RHE"].min(),
                _J_close_to_0["E_AppV_RHE"].max(),
            )
            dE = _i1 - _i0
            _J_0_window = swgrp.query("(E_AppV_RHE > @_i0) & (E_AppV_RHE < @_i1)")
            _hor_popt, _hor_pcov = curve_fit(
                hor_y,
                _J_0_window[EvRHE].to_numpy(),
                _J_0_window["Jcorr_raw"].to_numpy(),
            )
            _R = sum((_J_0_window.Jcorr_raw - _hor_popt) ** 2)
            _res.append(
                {
                    "hor_J0slice_j": c,
                    "hor_lin_B": _hor_popt[0],
                    "hor_Ewindow": dE,
                    "hor_lin_R": _R,
                    "hor_J0slice_len": len(_J_close_to_0),
                    "hor_E_lower": _i0,
                    "hor_E_upper": _i1,
                    "hor_Ewin_R_ratio": dE / _R,
                }
            )

    #        _J_highE_check_max = _J_highE_check_max.assign(**{'hor_Jcorr_raw_diff' :  _J_highE_check_max.Jcorr_raw - _hor_popt})
    hor_res = pd.DataFrame(_res)
    _hres_Rmin = hor_res.loc[
        (hor_res["hor_Ewindow"] > 0.021) & (hor_res["hor_E_lower"] > 0.5)
    ]["hor_lin_R"].min()

    _res_slice = hor_res.loc[
        (hor_res["hor_lin_R"] < 0.3)
        & (hor_res["hor_Ewindow"] > 0.03)
        & (hor_res["hor_E_lower"] > 0.75)
    ]
    if _res_slice.empty:
        _res_slice = hor_res.loc[
            (hor_res["hor_lin_R"] < 0.5)
            & (hor_res["hor_Ewindow"] > 0.03)
            & (hor_res["hor_E_lower"] > 0.5)
        ]
    if _res_slice.empty:
        _res_slice = hor_res.loc[
            (hor_res["hor_lin_R"] < 0.5)
            & (hor_res["hor_Ewindow"] > 0.021)
            & (hor_res["hor_E_lower"] > 0.5)
        ]

    if _res_slice.empty:
        _res_slice = hor_res.loc[
            (hor_res["hor_lin_R"] < _hres_Rmin * 1.1)
            & (hor_res["hor_Ewindow"] > 0.021)
            & (hor_res["hor_E_lower"] > 0.5)
        ]

    _res_mean = _res_slice.mean()
    _J_highE_check_max = swgrp.query(
        "(E_AppV_RHE > @_res_mean.hor_E_lower) & (E_AppV_RHE < @_res_mean.hor_E_upper)"
    )
    if len(_J_highE_check_max) > 10:
        _hor_popt, _hor_pcov = curve_fit(
            hor_y,
            _J_highE_check_max[EvRHE].to_numpy(),
            _J_highE_check_max["Jcorr_raw"].to_numpy(),
        )
    else:
        _hor_popt = _res_mean.hor_lin_B
    hor_overall_E_lower, hor_overall_E_upper = 0.8, 0.95
    _J_highE_overall = swgrp.query(
        "(E_AppV_RHE > @hor_overall_E_lower) & (E_AppV_RHE < @hor_overall_E_upper)"
    )
    #     _overall_popt, _overall_pcov = curve_fit(lin_y, _J_highE_overall[EvRHE].to_numpy(), _J_highE_overall['Jcorr_raw'].to_numpy())
    _overall_lin = linregress(
        _J_highE_overall[EvRHE].to_numpy(), _J_highE_overall["Jcorr_raw"].to_numpy()
    )
    _overall_lin_pars = dict(
        zip([f"hor_overall_{i}" for i in _overall_lin._fields], _overall_lin)
    )
    _overall_lin_pars.update(
        {
            "hor_overall_E_upper": hor_overall_E_upper,
            "hor_overall_E_lower": hor_overall_E_lower,
            "hor_overall_J_mean": _J_highE_overall["Jcorr_raw"].mean(),
        }
    )
    swgrp = swgrp.assign(**{"Jcorr_hor": swgrp.Jcorr_raw - _hor_popt})

    return swgrp, _res_mean, hor_res, _overall_lin_pars


# _J_highE_check_max.plot(EvRHE,y=['Jcorr_raw', 'hor_Jcorr_raw_diff'])
# _J_highE_overall.plot(EvRHE,y=['Jcorr_raw'])
# swgrp.plot(x=EvRHE,y=['Jcorr', 'Jcorr_hor','Jcorr_raw'])
#         _J_close_to_0.plot(EvRHE,y=['Jcorr_raw', 'hor_Jcorr_raw_diff'])
#    _J_highE_check_max = swgrp.query('(E_AppV_RHE > 0.75) & (E_AppV_RHE < 1)')
#    hor_y(swgrp[EvRHE].to_numpy(),*_hor_popt)
#    'hor_Jcorr_raw_diff' :  swgrp.Jcorr_raw - _hor_popt


def modeling_ORR_swp(_O2N2):
    def BV_func(x, i0, C_O, C_O_bulk, alpha, C_R, C_R_bulk):
        _f_RT = (constants.e * constants.Avogadro) / (constants.R * 293)
        return i0 * (
            (C_O / (1.3e-03)) * np.exp(-0.5 * _f_RT * (x))
            - (C_R / C_R_bulk) * np.exp((1 - 0.5) * _f_RT * (x))
        )

    def BV_simple0_no_CD(x, i0, alpha, E):
        _f_RT = (constants.e * constants.Avogadro) / (constants.R * 298)
        return -1 * i0 * (np.exp(-alpha * _f_RT * (x - 1.23))) / 0.237

    def BV_simple(x, k_rj_0i, alpha, Cp_i, E, E2, Cp_2, k_rj_02, c, d):
        _f_RT = (constants.e * constants.Avogadro) / (constants.R * 298)
        j1 = -1 * k_rj_0i * (np.exp(-alpha * _f_RT * (x - E)) * Cp_i) * 0.237
        # j2 = -1 * k_rj_02 * (np.exp(-alpha * _f_RT * (x - E2)) * Cp_2) * 0.237
        return j1 + (c * x + d)

    def fsigmoid(x, a, b, c, d):
        return 1.0 / (1 + np.exp(-a * (x - b))) + (c * x + d)

    _sig_popt, _sig_pcov = curve_fit(
        fsigmoid,
        _O2N2[EvRHE].to_numpy(),
        _O2N2["Jcorr_raw"].to_numpy(),
        method="dogbox",
    )
    _BV_popt, _BV_pcov = curve_fit(
        BV_func, _O2N2[EvRHE].to_numpy(), _O2N2["Jcorr_raw"].to_numpy()
    )
    _BVsimple_popt, _BVsimple_pcov = curve_fit(
        BV_simple, _O2N2[EvRHE].to_numpy(), _O2N2["Jcorr_raw"].to_numpy()
    )
    _O2N2 = _O2N2.assign(
        **{
            "sigmoid_Jcorr_raw": fsigmoid(_O2N2[EvRHE].to_numpy(), *_sig_popt),
            "BVcomp_Jcorr_raw": BV_func(_O2N2[EvRHE].to_numpy(), *_BV_popt),
            "BVsimpl_Jcorr_raw": BV_simple(_O2N2[EvRHE].to_numpy(), *_BVsimple_popt),
        }
    )
    _O2N2.plot(
        EvRHE, y=["Jcorr_raw", "sigmoid_Jcorr_raw", "BVsimpl_Jcorr_raw"], ylim=(-6, 1)
    )

    def fsigmoid(x, a, b, c, d):
        return 1.0 / (1 + np.exp(-a * (x - b))) + (c * x + d)

    def BV_simple0(x, i0, alpha, E, c, d):
        _f_RT = (constants.e * constants.Avogadro) / (constants.R * 298)
        return -1 * i0 * (np.exp(-alpha * _f_RT * (x - E))) * 0.237 + (c * x + d)

    _sig_popt, _sig_pcov = curve_fit(
        BV_simple0, _O2N2[EvRHE].to_numpy(), _O2N2["Jcorr_raw"].to_numpy()
    )
    _sig_pars = dict(
        zip([f"BVsimple_{i}" for i in ["i0", "alpha", "E", "c", "d"]], _sig_popt)
    )

    _O2N2 = _O2N2.assign(
        **{"sigmoid_Jcorr_raw": BV_simple0(_O2N2[EvRHE].to_numpy(), *_sig_popt)}
    )

    _O2N2.plot(EvRHE, y=["Jcorr_raw", "sigmoid_Jcorr_raw"], ylim=(-6, 1))


def ORR_apply_Jcorr_zero_correction(O2_join, apply=False):
    if apply:
        t1 = O2_join.loc[
            (np.isclose(O2_join["J_O2_diff"], -0.03, rtol=50))
            & (np.isclose(O2_join[EvRHE], 0.95, atol=100e-3)),
            :,
        ]
        #            O2_join.loc[(O2_join['Sweep_Type'] == sweep), 'Jcorr']
        fig, ax = plt.subplots()
        O2_join.groupby("Sweep_Type").plot(x=EvRHE, y="Jcorr", ax=ax)

        for sweep in t1["Sweep_Type"].unique():
            #                    t1.loc[(O2_join['Sweep_Type'] == sweep),:]
            #                    Zero_fit = linregress(t1.loc[(O2_join['Sweep_Type'] == sweep),EvRHE],t1.loc[(O2_join['Sweep_Type'] == sweep),'Jcorr'])
            if not "NA" in sweep:
                Jzero_mean = t1.loc[(O2_join["Sweep_Type"] == sweep), "Jcorr"].mean()
                O2_join.loc[(O2_join["Sweep_Type"] == sweep), "Jcorr"] = (
                    O2_join.loc[(O2_join["Sweep_Type"] == sweep), "Jcorr"] - Jzero_mean
                )
                _logger.info(
                    "Jkin calculation Jcorr minus %.6f for 0-correction" % Jzero_mean
                )
    return O2_join


def ORR_determine_E_onset(O2_join, PAR_file):
    global EvRHE
    E_onset = pd.DataFrame([])
    #    for swp, swgrp in O2_join.groupby('Sweep_Type')
    try:
        #                    E_onset = O2_join.loc[(O2_join['Jcorr'] > -0.119) & (O2_join['Jcorr'] < -0.0999) & (O2_join[EvRHE] <(J2nd_deriv_max[EvRHE].values[0]+0.2)) & (O2_join[EvRHE] > (J2nd_deriv_max[EvRHE].values[0]-0.2))& (O2_join[EvRHE] < 1.19),:]
        E_onset = (
            O2_join.loc[
                (O2_join["Jcorr"] > -0.129)
                & (O2_join["Jcorr"] < -0.0999)
                & (O2_join[EvRHE] < 0.99),
                :,
            ]
            .sort_values(by=EvRHE)
            .head(1)
        )
    except Exception as e:
        _logger.warning(f"Jkin calculation ORR; Jkin calc, E_onset expansion: {e}")
        i = 0
        while E_onset.empty:
            E_onset = O2_join.loc[
                (O2_join["Jcorr"] > -0.119 + i)
                & (O2_join["Jcorr"] < -0.0999 + i)
                & (O2_join[EvRHE] < 1),
                :,
            ]
            i += 0.04

    if E_onset.empty:
        _logger.warning(f"Jkin calculation E_onset empty! expanding... {PAR_file}")
        i = 0
        while E_onset.empty:
            E_onset = (
                O2_join.loc[
                    (O2_join["Jcorr"] > -1 * (0.529 + i))
                    & (O2_join["Jcorr"] < 0.0)
                    & (O2_join[EvRHE] < 1.19),
                    :,
                ]
                .sort_values(by=EvRHE)
                .head(1)
            )
            #                            E_onset = O2_join.loc[(O2_join['Jcorr'] > -0.119+i) & (O2_join['Jcorr'] < -0.0999+i) & (O2_join[EvRHE] < 1),:]
            i += 0.04
    return E_onset
    #                Diff_lim = O2_join.loc[np.isclose(O2_join['J_O2_diff'],0,atol=1E-07) | (O2_join['J_O2_diff'] == O2_join['J_O2_diff'].min()) & (O2_join[EvRHE] <  E_onset[EvRHE].mean()),:]


#                E_lowdiff = 0.03


def ORR_determine_Jdiff_lim(swgrp, PAR_file, E_limit_lowdiff):
    global EvRHE
    Diff_lim = pd.DataFrame([])
    #    swp,swgrp
    #    swgrp.plot(EvRHE,y=['Jcorr_raw', 'Jcorr_O2diff', 'Jcorr_O2diffdiff'])
    _diff_lim_templ = namedtuple("Diff_lim", "df min max")
    try:
        #                    E_onset[EvRHE].mean()-0.1
        _E_jdiffmax = swgrp.loc[
            swgrp.loc[
                np.isclose(swgrp.Jcorr_O2diffdiff, 0, atol=1e-1)
            ].Jcorr_O2diff.idxmax()
        ][EvRHE]
        #        .plot(x=EvRHE,y='Jcorr_O2diff',kind='scatter')
        #        _E_jdiffmax = swgrp.loc[swgrp.Jcorr_O2diff.idxmin()][EvRHE]
        _E_Diff_lim = swgrp.loc[
            swgrp.loc[swgrp[EvRHE] < _E_jdiffmax].Jcorr_O2diff.idxmin()
        ][EvRHE]
        Diff_lim = swgrp.loc[np.isclose(swgrp[EvRHE], _E_Diff_lim, atol=0.005)]
    #        Diff_lim = swgrp.loc[
    #                   np.isclose(swgrp['Jcorr_O2diff'], -0.003, rtol=50) & (swgrp[EvRHE] < 0.5) & (
    #                               swgrp['Jcorr_O2diff'] < 0) & (swgrp['Jcorr'] > -8), :]
    except Exception as e:
        _logger.warning(f"Jkin calculation, Diff lim problem: {e}")
        swgrp.plot(EvRHE, y=["Jcorr_raw", "Jcorr_O2diff", "Jcorr_O2diffdiff"])
        Diff_lim = swgrp.loc[
            np.isclose(swgrp["Jcorr_O2diff"], -0.003, rtol=1000)
            & (swgrp[EvRHE] < 0.8)
            & (swgrp["Jcorr_O2diff"] < 0)
            & (swgrp["Jcorr"] > -8),
            :,
        ]
    #                    Diff_lim = O2_join.loc[np.isclose(O2_join['J_O2_diff'],-0.003,rtol=100) & (O2_join[EvRHE] <  0.8) & (O2_join['J_O2_diff'] <  0),:]
    if Diff_lim.empty:
        swgrp.plot(EvRHE, y=["Jcorr_raw", "Jcorr_O2diff", "Jcorr_O2diffdiff"])
        _logger.warning(
            f"Jkin calculation Diff_lim empty! {PAR_file} take Jdiff at {E_limit_lowdiff:.2f} Vrhe"
        )
        Diff_lim = swgrp.loc[(swgrp[EvRHE] < E_limit_lowdiff), :]

    if (
        Diff_lim["Jcorr"].min()
        < swgrp.loc[(swgrp[EvRHE] < E_limit_lowdiff), "Jcorr"].min()
        and Diff_lim[EvRHE].min() > E_limit_lowdiff
    ):
        _logger.warning(
            f"Jkin calculation Diff_lim is smaller at higher potentials ! {PAR_file} taking Jdiff at {E_limit_lowdiff:.2f} Vrhe"
        )
        Diff_lim = swgrp.loc[(swgrp[EvRHE] < E_limit_lowdiff), :]

    diff_lim = _diff_lim_templ(
        Diff_lim, Diff_lim["Jcorr"].min(), Diff_lim["Jcorr"].max()
    )
    return diff_lim


def ORR_get_Tafel(swgrp, ORR_dest_dir_file, dest_file, **_pars):
    ### === TAFEL EQUATION AND PLOTTING ###
    #                TFxy = O2_join.loc[(E_half[EvRHE].mean()-0.20 < O2_join[EvRHE]) & (O2_join[EvRHE] < E_half[EvRHE].mean()+0.30) & (O2_join['Sweep_Type'] == 'cathodic')]
    _swgpr_meta = (
        swgrp[[i for i in swgrp.columns if swgrp[i].nunique() == 1]].iloc[0].to_dict()
    )
    TFxy = swgrp.loc[
        (0.50 < swgrp[EvRHE]) & (swgrp[EvRHE] < 0.85) & (swgrp["Jkin_min"] > 0),
        [EvRHE, "Jkin_min", "Jkin_max", "Jcorr", "Sweep_Type"],
    ]
    #                TFxy.plot(x=EvRHE,y='Jkin_min',xlim=(0.55,0.8),logy=True)
    #                 'log10_Jkin_max' : np.log10(TFxy['Jkin_max'])
    Tafel_ovv = pd.DataFrame([])
    TFlst, TafelDir = [], ORR_dest_dir_file.joinpath("TAFEL")
    TafelDir.mkdir(parents=True, exist_ok=True)

    try:
        TFxy = TFxy.assign(
            **{
                "log10_Jkin_min": np.log10(TFxy["Jkin_min"]),
                "E_overp": 1.23 - TFxy[EvRHE],
            }
        )
        TFxy["log10_Jkin_min"].dropna(inplace=True)
        TFxy["log_Jkinm_2ndDiff"] = TFxy["Jkin_min"].diff().diff().values
        #                    F,D,nu,C,A,Rc,Trt = 96485, 1.5E-05, 0.010, 1.1E-06,SA_disk,8.314,298
        #                    Tafel_fit_min = linregress(TFxy['E_overp'].values,TFxy['log10_Jkin_min'].values)
        w = 60
        rTFxy = TFxy.rolling(10).mean()
        swp = _swgpr_meta.get("Sweep_Type", "swp_unkown")
        for i in range(0, len(rTFxy) - w, w):
            _TFipars = _swgpr_meta
            rTFxy.dropna(inplace=True)
            tflogJx, tfEy = (
                rTFxy.iloc[i : i + w]["log10_Jkin_min"].values,
                rTFxy.iloc[i : i + w][EvRHE].values,
            )
            tfxy_Eoverp = rTFxy.iloc[i : i + w]["E_overp"].values
            TFfit = linregress(tflogJx, tfEy)
            #                        print(i,rTFxy.iloc[i][EvRHE],rTFxy.iloc[i+w][EvRHE],TFfit[2])
            _TF_collabel = f"log10_TAFEL_{i}"
            if np.abs(TFfit[2]) > 0.95:

                TFll = "log10_TAFEL_%.3fVrhe_%.0f" % (
                    rTFxy.iloc[i][EvRHE],
                    TFfit[0] * -1000,
                )
                TFxy = TFxy.assign(
                    **{_TF_collabel: (TFxy[EvRHE] - TFfit[1]) / TFfit[0]}
                )
                rTFxy = rTFxy.assign(**{TFll: (rTFxy[EvRHE] - TFfit[1]) / TFfit[0]})
                TF_E_mean = rTFxy.loc[
                    ((rTFxy.log10_Jkin_min - rTFxy[TFll]).abs() < 5e-03), EvRHE
                ].mean()

                rTFxy.plot(
                    x=EvRHE,
                    y=["log10_Jkin_min", TFll],
                    title=f"TS: {TFfit[0]*-1000:.0f} mV/dec\n {dest_file}_{swp}_{i} ",
                )
                plt.ylabel("log10_Jkin_min")
                plt.savefig(
                    TafelDir.joinpath(f"{dest_file}_{swp}_{i}.png"), bbox_inches="tight"
                )
                #                            plt.show()
                plt.close()
                #                            print(' TAFEL: ',rTFxy.iloc[i][EvRHE],rTFxy.iloc[i+w][EvRHE],TFfit[2])
                _TFipars.update(
                    {
                        "TF_index": i,
                        "TF_E_low": tfEy.min(),
                        "TF_E_high": tfEy.max(),
                        "TF_TS": TFfit[0] * -1000,
                        "TF_TSb": TFfit[1],
                        "TF_TSerr": TFfit[2],
                        "TF_E_mean": TF_E_mean,
                        "TF_label": _TF_collabel,
                        "TF_Eoverp_min": tfxy_Eoverp.min(),
                        "TF_Eoverp_max": tfxy_Eoverp.max(),
                    }
                )
            else:
                _TFipars.update(
                    {
                        "TF_index": i,
                        "TF_E_low": 0,
                        "TF_E_high": 0,
                        "TF_TS": 0,
                        "TF_TSb": 0,
                        "TF_TSerr": 0,
                        "TF_label": _TF_collabel,
                    }
                )
            TFlst.append(_TFipars)

        TF_out_base = TafelDir.joinpath(f"{dest_file}_{swp}.xlsx")
        TF_out = FileOperations.CompareHashDFexport(TFxy, TF_out_base)
        TFlst = [{**i, **{"TF_data_file": TF_out}} for i in TFlst]
        Tafel_ovv = pd.DataFrame(TFlst)
    #            'TF_data_file'
    #            TFxy.to_excel(TF_out_base)
    #                    index_info_ORR_TAFEL = {'PAR_file': PAR_file, 'DestFile': TF_out,
    #                                            'Type_output': 'ORR_Jkin_calc_Tafel', 'Type_exp': 'ORR'}
    #                    TFxy.to_excel(TF_out)
    #                    print('TAFEL SLOPE: %.2f mV/dec,\n OUTPUT: %s'%(Tafel_ovv.iloc[Tafel_ovv['TF_E_high'].idxmax()]['TS'],TF_out))
    except Exception as e:
        _logger.error("Jkin calculation ERROR in TAFEL of ORR: {0}".format(e))
        #            index_info_ORR_TAFEL = {}
        Tafel_ovv = pd.DataFrame(
            {
                **_swgpr_meta,
                **{
                    "TF_index": 0,
                    "TF_E_low": 0,
                    "TF_E_high": 0,
                    "TF_TS": 0,
                    "TF_TSb": 0,
                    "TF_TSerr": 0,
                },
            },
            index=[0],
        )
    return Tafel_ovv


# if Tafel_ovv.empty:
#             Tafel_ovv = pd.DataFrame({**_swgpr_meta,**{'TF_index': 0, 'TF_E_low': 0, 'TF_E_high': 0,

#        pd.DataFrame(data=KLoutRow, columns=['Sweep_Type', EvRHE, 'I_Disk_%s' % rpm_list[rpm],
#                                                     'I_Ring_%s' % rpm_list[rpm]])
#                pd.merge(left,right,on='SampleID',how='outer').query('SweepType == "cathodic"')
#                KLout2 = pd.concat([KLout2,KLrpm],axis=1)
#        KLout = pd.merge(KLout, KLrpm, on=[f'KL_{EvRHE}', 'Sweep_Type'], how='outer')
#        KLout = KLout.assign(
#            **{'Electrolyte': ORR_gr_ovv.iloc[0]['Electrolyte'], 'pH': ORR_gr_ovv.iloc[0]['pH']})

#                fig,ax = plt.subplots()
#                for rcath,rcathgrp in KLout.groupby('Sweep_Type'):
#                    rcathgrp.dropna(axis=0,subset=['I_Disk_1500_y'],inplace=True)
#                    rcathgrp.plot(x=EvRHE,y='I_Disk_1500_y',ax=ax)

#                print('ERROR KL prep',e )
# %%
#            out_row.update(next_row)
