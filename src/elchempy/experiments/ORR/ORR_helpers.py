"""
Created on Sun Nov 14 11:59:29 2021

@author: DW
"""


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
