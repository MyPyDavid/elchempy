"""
Created on Fri Nov 12 12:26:43 2021

@author: DW
"""


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
