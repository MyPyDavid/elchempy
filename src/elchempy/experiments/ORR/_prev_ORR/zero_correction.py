"""
Created on Fri Nov 12 12:32:51 2021

@author: DW
"""


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
