"""
Created on Fri Nov 12 12:27:27 2021

@author: DW
"""


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
