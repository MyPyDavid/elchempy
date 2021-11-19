"""
This is an UNUSED OLD VERSION
of the jkin calculation
"""


def _debug():
    fit_run_arg = orr_run_args[1]


def Jkin_calc_multi(fit_run_arg, **ORR_kwargs):
    try:
        #                fit_run_arg = orr_run_args[1]
        Jkin_calculations(fit_run_arg, **ORR_kwargs)
    #                fit_export_EV_all.append([fit_export_EV])
    except Exception as e:
        _logger.warning(
            f"ORR Jkin multi attempt failed for {fit_run_arg[0]} because {e}"
        )


def _rest():
    ORR_operations(self, O2_disk_seg, N2_BG)

    ORR_operations(self, O2_disk_seg, N2_BG)

    O2_join_raw, _N2BG_pars = ORR_operations.merge_O2_with_N2_BG(
        O2_disk_seg, N2_BG, self.mA
    )
    O2_join = ORR_operations.add_mean_Jcorr_col(O2_join_raw)

    O2_join, ORR_disk_pars = ORR_operations.ORR_disk_calc_pars(
        PAR_file,
        ORR_dest_dir_file,
        dest_file,
        O2_join,
        set_jcorr_col="Jcorr_minus_factor",
    )
    # TEST PLOTTING #
    #            O2_join.groupby('Sweep_Type').plot(x=EvRHE,y=['Jcorr_raw','Jcorr_hor'])
    ORR_disk_pars = ORR_disk_pars.assign(**_N2BG_pars).drop_duplicates()
    #            (PAR_file,rpm_n, seg, ORR_dest_dir_file, dest_file, O2_join)
    #                                'Analysis_date': datetime.now()}
    try:
        Iring_Chrono, _ORR_RRDE_pars = O2_Chrono(
            seg,
            rpm_n,
            O2_join,
            O2_chrono_Ring,
            ORR_dest_dir_file,
            plot_BipotChrono=True,
            export_excel=True,
            ring_filter=True,
        )
    #                   'Jring_050': out_Jring, 'FracH2O2_050': out_FracH2O2})
    except Exception as e:
        _logger.error(
            "Jkin calculation ERROR: Iring_Chrono {0}, {1}".format(e, dest_file)
        )
        Iring_Chrono, _ORR_RRDE_pars = pd.DataFrame(), pd.DataFrame()
    ### +++++++++-   ++++++++++++++++++++++    -+++++++++++ ###
    ### +++++++++- PLOTTING -+++++++++++ ###
    # https://stackoverflow.com/questions/37737538/merge-matplotlib-subplots-with-shared-x-axis
    #                fig, (axRing,axJ) = plt.subplots(2, sharex=True)
    #                SegGr.loc[:,[x]+y].to_csv(ORR_dest_dir.joinpath(Path(f).stem+'_%s'%rpm+'.csv'))
    #                                gr.plot(x=x,y=y,xlim=(0,1.2),ylim=(-6,0.5),title=Path(f).stem)
    #                                plt.savefig(ORR_dest_dir.joinath(Path(f).stem)+'_O2.png',bbox_inches='tight')
    #                                plt.close()
    ORR_plot_ring_disk(
        O2_join,
        ORR_disk_pars,
        Iring_Chrono,
        _ORR_RRDE_pars,
        ORR_dest_dir_file,
        dest_file,
    )


#             usecols_N2_correction = 'lin_jmA_lincorr'
#                     _O2N2.plot(x='E_AppV_RHE',y=['J_O2_diff','Jcorr','J_N2_scan'],xlim=(0.9,1.05),ylim=(-1E0,1E0))
#                    Jz(EvRHE,y=['Jcorr_raw'])
#                    _O2N2ero_mean = t1.loc[(O2_join['Sweep_Type'] == sweep), 'Jcorr'].mean()
#                O2_join.loc[(O2_join['Sweep_Type'] == sweep), 'Jcorr'] =
#                    O2_join.loc[(O2_join['Sweep_Type'] == sweep), 'Jcorr'] - Jzero_mean
#            Jcorr_mA = (O2_act_slice['j A/cm2'].values - N2_bg_local['j A/cm2'].values) * mA
#            I_O2_N2_diff = O2_act_slice['I(A)'].values - N2_bg_local['I(A)'].values
#            Newcols = {'RPM_n': int(rpm_n),'Segment': seg, 'Jcorr': Jcorr_mA,
#                        'Jcorr_raw' : Jcorr_mA, 'Icorr': I_O2_N2_diff,
#                       'J_N2_scan': N2_bg_local['j A/cm2'].values * mA, 'N2_PAR_file' : N2_bg_file_used,
#                       'J_O2_diff': (O2_act_slice['j A/cm2'].diff().rolling(30).mean().values) * 1E05,
#                       'Analysis_date': datetime.now()}
#            #                'J_O2_diff_diff' : (O2_act.loc[O2_segslice,'j A/cm2'].diff().rolling(10).mean().diff().rolling(10).mean().values)
#            O2_join = O2_act_slice.assign(**Newcols)
#            #                O2_join['J_2nd_diff'] = O2_join['J_O2_diff'].diff().rolling(50).mean().values*1E2
#            O2_join = ORR_apply_Jcorr_zero_correction(O2_join,apply=False)
# _ORR_obj, O2_disk_seg, N2_BG_scan, N2_col = er
# === TESTING CALCULATIONS ===
def _testing_function():
    calc = ORR_calculations(fit_run_arg)
    self = calc
    et = iter(calc)
    er = next(et)
    for er in et:
        print(er[-1])
        orr = ORR_operations(*er)
    self = orr
    kl = ORR_KL_loop(calc, run_loop=False)
    kl = ORR_KL_loop(calc, run_loop=True)
    self = kl
    self = ORR_ops


def Jkin_calculations(fit_run_arg, rpm_list, **ORR_kwargs):

    ### Input is the read-in data from 1 ORR PAR file !! ###

    # %%
    #    N2_scan, O2_act  = N2_background, grORR
    #    O2_act,O2_CVs,ovv,ORR_dest_dir,N2_bg_local = grORR,O2_CVs,ovv,ORR_dest_dir,N2_bg_local
    EvRHE = "E_AppV_RHE"
    globals()["EvRHE"] = "E_AppV_RHE"

    PAR_file_ovv_file, ORR_ovv_file, gr_ovv, N2_bg_file_used, N2_bg_local = fit_run_arg

    ORR_dest_dir_file = ORR_ovv_file.ORR_ecexp_destdir.unique()[0]
    fstem = ORR_ovv_file.PAR_fstem.unique()[0]
    #    PAR_file_ovv_file = Path(ORR_ovv_file.PAR_file.unique()[0])

    O2_CVs, O2_action = create_CVs(ORR_ovv_file)

    if O2_CVs.empty:
        _logger.warning(
            f"Not starting ORR -> ERROR empty for {ORR_dest_dir_file}/{fstem}"
        )

    else:
        _logger.info("Starting ORR for {0}".format(ORR_ovv_file))
        O2_act = O2_CVs.query(
            '(Gas == "O2") & (Type_action == "Cyclic Voltammetry (Multiple Cycles)") & (Scanrate == 0.01) & (SampleID != "Pt_ring")'
        )

    PAR_file = O2_act["PAR_file"].unique()[0]

    if not Path(PAR_file_ovv_file) == PAR_file:
        _logger.warning(
            f"ORR difference in PAR_files from OVV and CreateCV:\n{PAR_file_ovv_file} vs {PAR_file}"
        )
    else:
        ORR_file_PAR_date = ORR_ovv_file.PAR_date.unique()[0]

    #    O2_out_ovv = pd.DataFrame([])
    mA = 1000
    #    out_row, indexes_out, KIN_PARS   = pd.DataFrame([]), [], pd.DataFrame([])
    #    O2_act.assign(Sweep_Type=np.where(O2_act[EvRHE].diff() > 0, 'anodic', np.where(O2_act[EvRHE].diff() < 0, 'cathodic', 'NA')))
    O2_join = pd.DataFrame([])
    outdf, KLout = pd.DataFrame([]), pd.DataFrame(
        data=[], columns=["Sweep_Type", EvRHE]
    )

    electrode_name, collN, SA_disk, SA_ring = WE_SA_collection_eff("PINE").values()
    #        N2_scan = preN2_scan.assign(**{'j A/cm2' : preN2_scan.loc[:,'j A/cm2']/N2_factor})
    #    print('!! N2 background wrong scan rate: %s, so j divided by %.0f '%(ScanRates.min(),N2_factor))
    #    KLout2 = pd.DataFrame([])
    #%%
    #    for h1 in O2_act['hash'].unique():
    #        h1 = O2_act['hash'].unique()[0]
    O2_segs = O2_act["Segment #"].unique()
    #        O2_act_SampleID = O2_act['SampleID'].unique()[0]

    if O2_act["PAR_file"].nunique() > 1:
        _logger.error(
            f'Jkin calculation multiple PAR_files in input DF {", ".join([str(i) for i in O2_act.PAR_file.unique()])}'
        )

    #        print(PAR_file,':',O2_segs)
    ###### ====== CHECKING FOR N2 scans from OVV index ========== #######
    #        ORR_ovv = ovv.loc[ovv['PAR_file'] == str(PAR_file)]
    O2_chrono_Ring, O2_chrono_actions = ORR_read_Ring_file(gr_ovv, ORR_ovv_file)
    #                                                           PAR_file, fstem, ORR_file_PAR_date)
    #    ORR_read_Ring_file(gr_ovv, PAR_file, ORR_file_PAR_date, fstem)
    if (
        "DAC Control" not in O2_action.Type_action.unique()
        and "DAC Control" in O2_chrono_actions.Type_action.unique()
        and sum(O2_act.RPM_DAC.unique()) == 0
    ):
        O2_chrono_Ring[["Segment #", "RPM_DAC"]]
        _seg_rpm_cols = ["Segment #", "RPM_DAC"]
        _seg_rpms_chrono = set(
            [(i[0], i[1]) for i in O2_chrono_Ring[_seg_rpm_cols].values]
        )
        _seg_rpm_DF = pd.DataFrame(_seg_rpms_chrono, columns=_seg_rpm_cols)
        O2_act = pd.merge(
            O2_act.drop(columns="RPM_DAC"), _seg_rpm_DF, on=_seg_rpm_cols[0], how="left"
        )
        _logger.warning(f"Jkin calculation, RPM_DAC used from Ring {fstem}")
        # TODO Match Segments with RPM_DACs of Ring and Merge with Disk
    if N2_bg_local.empty and N2_bg_file_used:
        try:
            N2_bg_local = pd.read_excel(N2_bg_file_used, index_col=[0])
        except Exception as e:
            _logger.error(
                "Jkin calculation N2 file can not be read make fake copy with 0 current: {0}".format(
                    e
                )
            )
            N2_copy = O2_act.loc[(O2_act["Segment #"] == O2_segs[0])]
            N2_bg_local = N2_copy.assign(
                **{
                    "j A/cm2": 0,
                    "I(A)": 0,
                    "jmAcm-2": 0,
                    "Gas": "N2",
                    "N2_fake_0_current": True,
                }
            )
    else:
        pass
    #        N2_scan = N2_bg_local
    #        N2_bg_file_used = N2_bg_file_used
    #        if N2_scan.empty:
    #            N2_bg_PARf = ORR_gr_ovv.ORR_act_N2_bg.values[0]
    ##        N2_bg_PARf = ovv.loc[ovv['PAR_file'] == str(PAR_file)].ORR_act_N2_bg.values[0]
    #        else:
    #            N2_scan = N2_bg_local
    ###### ====== CHECKING FOR O2 scans, segments, rotation
    # ========== #######
    if O2_act.PAR_file.str.contains("Pt-ring").any():
        _logger.error(f"Jkin calculation PAR_file disk containts Pt-ring {PAR_file}")
        print("Pt-ring")
    newsegs = []
    for nums1, s1 in enumerate(O2_segs):
        _O2_act_seg = O2_act.loc[(O2_act["Segment #"] == s1)]
        seglen = len(_O2_act_seg)
        if seglen > 1000:
            _O2_act_seg.RPM_DAC.unique()
            newsegs.append(s1)
    #%%
    _ORR_parslst = []
    _KL_data_all_rpms = []
    for rpm, seg in enumerate(newsegs):
        rpm
        #            rpm,seg = 0,newsegs[0]
        #%%
        try:
            dest_file = fstem + "_" + str(int(rpm_list[rpm]))
            _logger.info("Jkin calculation DEST ORR: {0}".format(dest_file))

            rpm_n = rpm_list[rpm]
            O2_act_slice = O2_act.loc[(O2_act["Segment #"] == seg)]

            O2_join_raw, _N2BG_pars = merge_O2_with_N2_BG(O2_act_slice, N2_bg_local, mA)
            O2_join = add_mean_Jcorr_col(O2_join_raw)
            #             usecols_N2_correction = 'lin_jmA_lincorr'
            #                     _O2N2.plot(x='E_AppV_RHE',y=['J_O2_diff','Jcorr','J_N2_scan'],xlim=(0.9,1.05),ylim=(-1E0,1E0))
            #                    Jz(EvRHE,y=['Jcorr_raw'])
            #                    _O2N2ero_mean = t1.loc[(O2_join['Sweep_Type'] == sweep), 'Jcorr'].mean()
            #                O2_join.loc[(O2_join['Sweep_Type'] == sweep), 'Jcorr'] =
            #                    O2_join.loc[(O2_join['Sweep_Type'] == sweep), 'Jcorr'] - Jzero_mean
            #            Jcorr_mA = (O2_act_slice['j A/cm2'].values - N2_bg_local['j A/cm2'].values) * mA
            #            I_O2_N2_diff = O2_act_slice['I(A)'].values - N2_bg_local['I(A)'].values
            #            Newcols = {'RPM_n': int(rpm_n),'Segment': seg, 'Jcorr': Jcorr_mA,
            #                        'Jcorr_raw' : Jcorr_mA, 'Icorr': I_O2_N2_diff,
            #                       'J_N2_scan': N2_bg_local['j A/cm2'].values * mA, 'N2_PAR_file' : N2_bg_file_used,
            #                       'J_O2_diff': (O2_act_slice['j A/cm2'].diff().rolling(30).mean().values) * 1E05,
            #                       'Analysis_date': datetime.now()}
            #            #                'J_O2_diff_diff' : (O2_act.loc[O2_segslice,'j A/cm2'].diff().rolling(10).mean().diff().rolling(10).mean().values)
            #            O2_join = O2_act_slice.assign(**Newcols)
            #            #                O2_join['J_2nd_diff'] = O2_join['J_O2_diff'].diff().rolling(50).mean().values*1E2
            #            O2_join = ORR_apply_Jcorr_zero_correction(O2_join,apply=False)
            # === TEST PLOTTING
            #                t1.plot(x='E_AppV_RHE',y=['J_O2_diff','Jcorr','J_N2_scan'],xlim=(0.9,1.05),ylim=(-1E0,1E0))
            #                fig,ax = plt.subplots()
            #                O2_join.plot(x='E_AppV_RHE',y=['jmAcm-2','Jcorr','J_N2_scan'],xlim=(0,1.1),ylim=(-6E0,5E0),ax=ax)
            #                plt.savefig(ORR_dest_dir.joinpath(dest_file+'_test.png'),dpi=100,bbox_inches='tight')
            #                plt.close()
            #                Jderiv_Diffs = O2_join.loc[(O2_join['J_O2_diff'] == O2_join['J_O2_diff'].max()) | (O2_join['J_O2_diff'] == O2_join['J_O2_diff'].min()) & (O2_join[EvRHE] < 1.0),:]
            #                N2_scan.plot(x='E(V)',y='I(A)',xlim=(-1,1))
            ### ==== Calculation of Jdiff, then Jkin IMPORTANT PARAMETERS ====== ###
            ### Calculation of Jdiff, then Jkin IMPORTANT PARAMETERS !! ###
            #            'Jcorr_minus_factor', 'Jcorr_minus_lin'
            O2_join, ORR_disk_pars = ORR_disk_calc_pars(
                PAR_file,
                ORR_dest_dir_file,
                dest_file,
                O2_join,
                set_jcorr_col="Jcorr_minus_factor",
            )
            # TEST PLOTTING #
            #            O2_join.groupby('Sweep_Type').plot(x=EvRHE,y=['Jcorr_raw','Jcorr_hor'])

            ORR_disk_pars = ORR_disk_pars.assign(**_N2BG_pars).drop_duplicates()
            #            (PAR_file,rpm_n, seg, ORR_dest_dir_file, dest_file, O2_join)
            #                                'Analysis_date': datetime.now()}
            try:
                Iring_Chrono, _ORR_RRDE_pars = O2_Chrono(
                    seg,
                    rpm_n,
                    O2_join,
                    O2_chrono_Ring,
                    ORR_dest_dir_file,
                    plot_BipotChrono=True,
                    export_excel=True,
                    ring_filter=True,
                )
            #                   'Jring_050': out_Jring, 'FracH2O2_050': out_FracH2O2})
            except Exception as e:
                _logger.error(
                    "Jkin calculation ERROR: Iring_Chrono {0}, {1}".format(e, dest_file)
                )
                Iring_Chrono, _ORR_RRDE_pars = pd.DataFrame(), pd.DataFrame()
            ### +++++++++-   ++++++++++++++++++++++    -+++++++++++ ###
            ### +++++++++- PLOTTING -+++++++++++ ###
            # https://stackoverflow.com/questions/37737538/merge-matplotlib-subplots-with-shared-x-axis
            #                fig, (axRing,axJ) = plt.subplots(2, sharex=True)
            #                SegGr.loc[:,[x]+y].to_csv(ORR_dest_dir.joinpath(Path(f).stem+'_%s'%rpm+'.csv'))
            #                                gr.plot(x=x,y=y,xlim=(0,1.2),ylim=(-6,0.5),title=Path(f).stem)
            #                                plt.savefig(ORR_dest_dir.joinath(Path(f).stem)+'_O2.png',bbox_inches='tight')
            #                                plt.close()
            ORR_plot_ring_disk(
                O2_join,
                ORR_disk_pars,
                Iring_Chrono,
                _ORR_RRDE_pars,
                ORR_dest_dir_file,
                dest_file,
            )
            ### +++++++++- FINALIZING -+++++++++++ ###
            #                O2_out_ovv,O2_CalcOut = pd.DataFrame([]), pd.DataFrame([])

            O2_CalcOut = pd.merge(
                O2_join, Iring_Chrono, on=[EvRHE, "Sweep_Type"], how="left"
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

            ORR_Excel_dest_base = ORR_dest_dir_file.joinpath(f"{dest_file}_RRDE.xlsx")
            ORR_Excel_dest = FileOperations.CompareHashDFexport(
                O2_CalcOut[yOut], ORR_Excel_dest_base
            )

            O2_ParsOut = (
                pd.merge(ORR_disk_pars, _ORR_RRDE_pars)
                .assign(**{"RRDE_merged_data_file": ORR_Excel_dest})
                .drop_duplicates()
            )
            _ORR_parslst.append(O2_ParsOut)
            ORR_KL_data_rpm = ORR_select_KL_data(O2_CalcOut)
            _KL_data_all_rpms.append(ORR_KL_data_rpm)
        #            _logger.info('Jkin calculation ORR EXCEL DEST:{0}'.format(ORR_Excel_dest))
        #                for sweepname, O2swpgrp in O2_CalcOut.groupby('Sweep_Type'):
        #                    ORR_Excel_dest_base_swp = ORR_dest_dir.joinpath(dest_file+'RRDE_{0}.xlsx'.format(sweepname))
        #                    ORR_Excel_dest_swp = FolderOps.FileOperations.CompareHashDFexport(O2swpgrp[yOut],ORR_Excel_dest_base_swp)
        #            index_info_ORR = (
        #            {'PAR_file': PAR_file, 'DestFile': ORR_Excel_dest, 'Type_output': 'ORR_Jkin_calc_RRDE',
        #             'Type_exp': 'ORR', 'ORR_N2_PAR_file': Path(N2_bg_file_used).name,
        #             'ORR_RING_PAR_file': Ring_ovv.PAR_file.values[0]})
        #                if O2_out_ovv.empty:
        #                    O2_out_ovv = O2_CalcOut
        #                else: # O2_out_ovv.empty == False
        #                    O2_out_ovv = O2_out_ovv.append(O2_CalcOut,ignore_index=True)
        #                O2_out_file = ORR_dest_dir.joinpath(dest_file+'_RING_' +str(rpm_list[rpm])+'.csv')
        #                print('RING OUTPUT',O2_out_file)
        #                O2_out_ovv.loc[O2_out_ovv['Segment #'] == seg].to_csv(O2_out_file)
        except Exception as e3:
            _logger.error(
                "Jkin calculation %s: no ORR at %.f because:\n %s"
                % (O2_act["SampleID"].unique()[0], int(rpm_list[rpm]), e3)
            )
            _ORR_parslst.append(
                pd.DataFrame(
                    O2_join[[i for i in O2_join.columns if O2_join[i].nunique() == 1]]
                    .iloc[0]
                    .to_dict(),
                    index=[0],
                )
            )
            _logger.error(
                "Jkin calculation ORR EXCEL DEST:{0}".format(
                    ORR_dest_dir_file.joinpath(f"{dest_file}_RRDE.xlsx")
                )
            )
    #        #%%
    #%%
    #        outdf = pd.concat([outdf, nextdf], ignore_index=True, sort=False)
    #    outdf = pd.concat([outdf,N2_Cdl_ovv],ignore_index=True)
    #    plt.close()
    _logger.info("Jkin calculation ORR Jkin calc succes: %s" % ORR_dest_dir_file)
    _dt_now = datetime.now()
    KL_data = pd.concat(_KL_data_all_rpms)
    KL_data = KL_data.assign(**{"ORR_KL_Analysis_date_dt": _dt_now})

    #    KLout = KLout.loc[(KLout['Sweep_Type'].str.contains("cathodic|anodic")), :]
    try:
        #            KLdest = 'KL_'+os.path.basename(O2_act.loc[O2_act['hash'] == h1,'File'].unique()[0]).split('.')[0]
        KL_fit = KL_plots(KL_data, ORR_dest_dir_file)
    #        KL_plots(KL_data, ORR_dest_dir_file)
    except Exception as e:
        _logger.error("Jkin calculation No KL FITTING plot:======= {0}".format(e))
        KL_fit = ["NA", "N"]

    ORR_pars_all_rpms = pd.concat(_ORR_parslst)
    ORR_pars_all_rpms = ORR_pars_all_rpms.assign(**{"ORR_Analysis_date_dt": _dt_now})

    ORR_pars_target_base = ORR_dest_dir_file.joinpath(f"ORR_pars_{fstem}")
    ORR_pars_target = FileOperations.CompareHashDFexport(
        ORR_pars_all_rpms, ORR_pars_target_base
    )
    _logger.info("Jkin calculation ORR Pars EXCEL DEST:{0}".format(ORR_pars_target))


#    index_info_ORR_pars = {'PAR_file': PAR_file, 'DestFile': ORR_pars_target, 'Type_output': 'ORR_Jkin_calc_Pars',
#                           'Type_exp': 'ORR'}
#    ORR_indexes = [index_info_ORR, index_info_ORR_O2Chrono, index_info_ORR_TAFEL, index_info_ORR_KL,
#                   index_info_ORR_KL_pars, index_info_ORR_pars]
#    return ORR_indexes
# try:
#            next_row = {'SampleID': O2_join['SampleID'].unique()[0], 'PAR_file': O2_join['PAR_file'].unique()[0],
#                        'DATE': O2_join['DATE'].unique()[0],
#                        'RPM': rpm_list[rpm], 'Segment': seg, 'E_onset': E_onset[EvRHE].mean(),
#                        'E_half': E_half[EvRHE].mean(), 'J_diff_lim': Diff_lim['Jcorr'].min(),
#                        'Jkin_075': np.abs(Jkin_075['Jkin_min'].mean()),
#                        'Jkin_080': np.abs(Jkin_080['Jkin_min'].mean()),
#                        'TSa_l': Tafel_ovv.iloc[Tafel_ovv['E_low'].idxmax()]['TS'],
#                        'TSb_l': Tafel_ovv.iloc[Tafel_ovv['E_low'].idxmax()]['TSb'],
#                        'TSa_h': Tafel_ovv.iloc[Tafel_ovv['E_high'].idxmin()]['TS'],
#                        'TSb_h': Tafel_ovv.iloc[Tafel_ovv['E_high'].idxmin()]['TSb'],
#                        'Jring_050': out_Jring, 'FracH2O2_050': out_FracH2O2,
#                        'Electrolyte': ORR_gr_ovv.iloc[0]['Electrolyte'],
#                        'pH': ORR_gr_ovv.iloc[0]['pH'], 'postAST': ORR_gr_ovv.iloc[0]['postAST'],
#                        'Analysis_date': datetime.now()}
#        except Exception as e:
#            _logger.error('Jkin calculation ERROR in ORR next row', e
