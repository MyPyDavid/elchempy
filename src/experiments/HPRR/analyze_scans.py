def HPRR_scan(All_HPRR, HPRR_ovv_file, dest_dir):
    #    All_HPRR,dest_dir = Samples_ovv, Path(HPRR_ovv_file.Dest_dir.iloc[0])
    #    EvRHE = 'E_AppV_RHE'
    # %%    Eapp = 'E_Applied_VRHE'
    EvRHE = "E_AppV_RHE"
    HPRR_dest_dir = dest_dir.joinpath("HPRR_scans")
    HPRR_dest_dir.mkdir(parents=True, exist_ok=True)
    #    make_sure_path_exists(HPRR_dest_dir)
    SampleID = All_HPRR["SampleID"].unique()[0]
    All_HPRR = All_HPRR.assign(
        **{
            "jmAcm-2": All_HPRR["j A/cm2"] * 1000,
            "Abs_jmAcm-2": np.abs(All_HPRR["j A/cm2"] * 1000),
            "log_Abs_jmAcm-2": np.log10(np.abs(All_HPRR["j A/cm2"] * 1000)),
        }
    )
    HPRR_CV = All_HPRR.query(
        'EXP == "HPRR" & ScanRate_calc < 0.02 & SampleID != "Pt_ring" & Type_action == "Cyclic Voltammetry (Multiple Cycles)" '
    )
    #    HPRR_fn = Path(HPRR_ovv['PAR_file'].unique()[0]).stem
    HPRR_PAR_fn = Path(HPRR_ovv_file.PAR_file.iloc[0])
    HPRR_fn = HPRR_PAR_fn.stem
    HPRR_out_lst = []
    if HPRR_ovv_file.empty:
        #           ovv[~ovv['SampleID'].str.contains('Pt_ring')].loc[:,['PAR_exp' == 'N2']].empty:
        logger.warning("!! Critical HPRR empty: {0}!!".format(dest_dir))
    try:
        grA = HPRR_CV.groupby(
            by=["Gas", "Type_action", "EXP", "Scanrate", "RPM_DAC", "Segment #"]
        )
        #        grB = HPRR_CV.groupby(by=['Gas','Type','EXP'])
        #        for scan in grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','N2_act')):
        #            print(scan)
        #                grA.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','0.1'))
        #        grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','HPRR')).to_csv(HPRR_dest_dir.joinpath('%s.csv' %HPRR_fn))
        #        hp_data = grB.get_group(('N2','Cyclic Voltammetry (Multiple Cycles)','HPRR'))
        out, HPRR_out = [], []
        for nm, gr in grA:
            for swnm, sweep in gr.groupby(by="Sweep_Type"):
                if swnm == "NA":
                    continue

                swp_target_file = HPRR_dest_dir.joinpath(
                    "{0}_{1}_{2}.xlsx".format(swnm, nm[4], HPRR_fn)
                )
                try:
                    old_target = HPRR_dest_dir.joinpath(
                        "HPRR_Tafel_{0}_{1}.xlsx".format(swnm, nm[4])
                    )
                    if old_target.is_file():
                        old_target.unlink()
                        logger.warning(
                            "HPRR output deleted old target: {0}".format(old_target)
                        )
                except:
                    logger.warning(
                        "HPRR output delete old target fail: {0}".format(old_target)
                    )

                swp = sweep.loc[
                    :,
                    [
                        EvRHE,
                        "jmAcm-2",
                        "Abs_jmAcm-2",
                        "log_Abs_jmAcm-2",
                        "Sweep_Type",
                        "RPM_DAC",
                    ],
                ]
                j_use = "jmAcm-2_fltr"
                #                .rolling(window=5).mean()
                swp[j_use] = scipy.signal.savgol_filter(swp["jmAcm-2"], 21, 3)
                swp = swp.assign(
                    **{
                        "Abs_jmAcm-2_fltr": np.abs(swp["jmAcm-2_fltr"]),
                        "log_Abs_jmAcm-2_fltr": np.log10(np.abs(swp["jmAcm-2_fltr"])),
                        "j/E": swp[j_use] / swp[EvRHE],
                        "dJ": swp[j_use].diff(),
                        "d/d2": swp[j_use].diff().diff(),
                        "dE": swp[EvRHE].diff(),
                    }
                )
                swp["dj/dE"] = swp["dJ"] / swp["dE"]
                #                swp.plot(x=EvRHE,y=['log_Abs_jmAcm-2_fltr','jmAcm-2_fltr','dj/dE'])
                #                rw[EvRHE]
                #                rw['Jabs'] = np.log10(np.abs(rw['jmAcm-2']))
                #                swp['j/E'] =swp[j_use]/swp[EvRHE]
                #                swp['dJ'] = swp[j_use].diff()
                #                swp['dJ/d2'] = swp[j_use].diff().diff()
                #                swp['dE'] = swp[EvRHE].diff()
                #                swp['dj/dE'] = swp['dJ']/swp['dE']
                ###### ======Analyzing HPRR CV and extracting kinetic parameters ========== #######
                HPOR = swp.loc[
                    (np.isclose(swp[EvRHE], swp[EvRHE].max() - 0.002, atol=0.010))
                ][j_use].mean()
                HPRR_08 = swp.loc[(np.isclose(swp[EvRHE], 0.8, atol=0.010))].head(1)
                HPRR_02 = swp.loc[(np.isclose(swp[EvRHE], 0.14, atol=0.010))].head(1)
                HPRR_onset = (
                    swp.loc[
                        (swp[j_use] > -0.129)
                        & (swp[j_use] < -0.0999)
                        & (swp[EvRHE] < 0.85),
                        :,
                    ]
                    .sort_values(by=EvRHE)
                    .head(1)
                )

                swp_08_TF = swp.loc[
                    (swp[EvRHE] <= HPRR_onset[EvRHE].values[0] + 0.020)
                    & (swp[EvRHE] >= HPRR_onset[EvRHE].values[0] - 0.020),
                    :,
                ]
                TF08fit = linregress(swp_08_TF[EvRHE].values, swp_08_TF[j_use].values)
                TF08_out = [
                    "E_onset",
                    HPRR_onset[EvRHE].iloc[0],
                    TF08fit[0],
                    TF08fit[1],
                    TF08fit[2],
                ]

                swp_02_TF = swp.loc[
                    (swp[EvRHE] <= HPRR_02[EvRHE].iloc[0] + 0.020)
                    & (swp[EvRHE] >= HPRR_02[EvRHE].iloc[0] - 0.020),
                    :,
                ]
                TF02fit = linregress(swp_02_TF[EvRHE].values, swp_02_TF[j_use].values)
                TF02_out = [
                    "E_0.2",
                    HPRR_02[EvRHE].iloc[0],
                    TF02fit[0],
                    TF02fit[1],
                    TF02fit[2],
                ]

                swp_11_TF = swp.loc[
                    (swp[EvRHE] <= swp[EvRHE].max() - 0.010)
                    & (swp[EvRHE] >= swp[EvRHE].max() - 0.050),
                    :,
                ]
                TF11fit = linregress(swp_11_TF[EvRHE].values, swp_11_TF[j_use].values)
                TF11_out = [
                    "E_max",
                    swp[EvRHE].max(),
                    TF11fit[0],
                    TF11fit[1],
                    TF11fit[2],
                ]

                E_j0 = swp[EvRHE].loc[swp["log_Abs_%s" % j_use].idxmin()]
                swp_j0_TF = swp.loc[
                    (swp[EvRHE] <= E_j0 + 0.050) & (swp[EvRHE] >= E_j0 - 0.050), :
                ]
                TFj0fit = linregress(swp_j0_TF[EvRHE].values, swp_j0_TF[j_use].values)
                TFj0_out = ["E_j0", E_j0, TF11fit[0], TF11fit[1], TF11fit[2]]

                swp_Tafel = swp.loc[
                    (swp[EvRHE] <= E_j0 + 0.15) & (swp[EvRHE] >= E_j0 - 0.15), :
                ]
                swp_Tafel_red, swp_Tafel_ox = (
                    swp_Tafel.loc[swp_Tafel[EvRHE] < E_j0 - 0.040, :],
                    swp_Tafel.loc[swp_Tafel[EvRHE] > E_j0 + 0.040, :],
                )
                Tafel_red = linregress(
                    swp_Tafel_red["log_Abs_%s" % j_use].values,
                    swp_Tafel_red[j_use].values,
                )
                Tafel_ox = linregress(
                    swp_Tafel_ox["log_Abs_%s" % j_use].values,
                    swp_Tafel_ox[j_use].values,
                )
                #                swp_Tafel_red.plot(x=EvRHE,y=['log_Abs_jmAcm-2_fltr','jmAcm-2_fltr'])
                #                swp_Tafel_ox.plot(x=EvRHE,y=['log_Abs_jmAcm-2_fltr','jmAcm-2_fltr'])
                Tafel_red_out = [
                    "Tafel_red",
                    E_j0,
                    np.abs(Tafel_red[0]) * 100,
                    Tafel_red[1],
                    Tafel_red[2],
                ]
                Tafel_ox_out = [
                    "Tafel_ox",
                    E_j0,
                    np.abs(Tafel_ox[0]) * 100,
                    Tafel_ox[1],
                    Tafel_ox[2],
                ]
                ###### ======Saving all HPRR CV kinetic parameters to file and index ========== #######
                TF_lst = [
                    TF08_out,
                    TF02_out,
                    TF11_out,
                    TFj0_out,
                    Tafel_red_out,
                    Tafel_ox_out,
                ]
                TF_out = pd.DataFrame(
                    TF_lst,
                    columns=[
                        "E_name",
                        "E_fit_HPRR",
                        "fit_slope_HPRR",
                        "fit_intercept_HPRR",
                        "fit_r_HPRR",
                    ],
                )
                TF_index = pd.DataFrame(
                    {
                        "SampleID": SampleID,
                        "Sweep_Type_HPRR": swnm,
                        "RPM_HPRR": nm[4],
                        "Gas": nm[0],
                        "Type_action": nm[1],
                        "EXP": nm[2],
                        "Scanrate": nm[3],
                        "Analysis_date": datetime.now(),
                        "DataFile": swp_target_file,
                        "PAR_file": HPRR_PAR_fn,
                        "E_name": TF_out.E_name,
                    }
                )
                HPRR_out_swp = pd.merge(TF_out, TF_index, on="E_name")
                #                rw = rw.assign(**{'HPRR_TF_Fit' : (rw[EvRHE]-TF08fit[1])/TF08fit[0], 'HPRR_0.2_Fit' : (rw[EvRHE]-TF02fit[1])/TF02fit[0],
                #                                  'HPRR_1.1_Fit' : (rw[EvRHE]-TF11fit[1])/TF11fit[0],'HPRR_j0_Fit' : (rw[EvRHE]-TFj0fit[1])/TFj0fit[0]})
                #                rwTF = rwTF.assign(**{'HPRR_TF_Fit' : (rwTF[EvRHE]-TF08fit[1])/TF08fit[0]})
                swp = swp.assign(
                    **{
                        "HPRR_TF_Fit": (swp[EvRHE] * TF08fit[0]) + TF08fit[1],
                        "HPRR_0.2_Fit": (swp[EvRHE] * TF02fit[0]) + TF02fit[1],
                        "HPRR_1.1_Fit": (swp[EvRHE] * TF11fit[0]) + TF11fit[1],
                        "HPRR_j0_Fit": (swp[EvRHE] * TFj0fit[0]) + TFj0fit[1],
                        "HPRR_Tafel_red_Fit": (swp[EvRHE] * Tafel_red[0])
                        + Tafel_red[1],
                        "HPRR_Tafel_ox_Fit": (swp[EvRHE] * Tafel_ox[0]) + Tafel_ox[1],
                    }
                )
                swp.to_excel(swp_target_file)
                logger.info("HPRR output because to: {0}".format(swp_target_file))
                #                rwTF = rwTF.assign(**{'HPRR_TF_Fit' : (rwTF[EvRHE]-TF08fit[1])/TF08fit[0]})

                #                        print(TFfit)
                #                        print(i,rTFxy.iloc[i][EvRHE],rTFxy.iloc[i+w][EvRHE],TFfit[2])
                HPRR_out_lst.append(HPRR_out_swp)
                #                out.append({'SampleID' : SampleID,'Sweep_Type_HPRR' : swnm,'RPM_HPRR' : nm[4], 'Groupnm' : nm,'j_HPOR' : HPOR,
                #                            'HPRR_onset' : HPRR_onset[EvRHE],'Lin_fit_Onset' : TF08fit,'Lin_fit_0.2' : TF02fit,'Lin_fit_1.1' : TF11fit,
                #                            'Lin_fit_j0' : TFj0fit,'HPRR_08' : HPRR_08,'HPRR_02' : HPRR_02, 'Analysis_date' : datetime.now(),'DataFile' :  swp_target_file})
                #            rwd = rw[EvRHE].diff()
                #            rwd['dJ'] = rw['j A/cm2'].diff()
                #            gr.rolling(window=5,on=[EvRHE,'j A/cm2']).mean()
                #                rw.plot(x=EvRHE,y='j A/cm2',kind='scatter',ylim=(-0.001,0.001),label='%s_%s' %(swp,nm))
                #                rw.plot(x=Eapp,y=dj,ylim=(-10,10),label='%s_%s' %(swp,nm))
                #                rw.dropna(axis=0).plot(x=EvRHE,y=['dJ','dJ/d2'],xlim=(0.5,1))
                #                fig,ax = plt.subplots()
                #                plt.title('%s scan of %s at %s' %(sweep,SampleID,nm[4]))
                #                rw.dropna(axis=0).plot(x=EvRHE,y=['jmAcm-2','log_Abs_jmAcm-2'],ylim=(-5,5),xlim=(0,1.2),label='%s_%s' %(sweep,nm[4]),ax=ax)
                #                rw.dropna(axis=0).plot(x=EvRHE,y='log_Abs_jmAcm-2',ylim=(-5,5),xlim=(0,1.2),label='%s_%s' %(sweep,nm[4]),ax=ax)
                #                plt.savefig(HPRR_dest_dir+'\\HPRR_%s_%s.png' %(sweep,nm[4]),dpi=300,bbox_inches='tight')
                #                plt.close()
                if np.abs(TF_out["fit_r_HPRR"].mean()) > 0.11:
                    swp.plot(
                        x=EvRHE,
                        y=[
                            j_use,
                            "HPRR_TF_Fit",
                            "HPRR_0.2_Fit",
                            "HPRR_1.1_Fit",
                            "HPRR_j0_Fit",
                        ],
                        xlim=(0, 1.2),
                        label=[j_use, "TF", "0.2", "1", "j0"],
                    )
                    plt.legend(ncol=3)
                    plt.grid(True)
                    #                    xlim=(0,1.2),label=['%s_%s' %(swnm,nm),'TF','0.2','1','j0'])
                    plt.title("%s scan of %s at %s" % (swnm, SampleID, nm[4]))

                    swp_target_png = HPRR_dest_dir.joinpath(
                        "{0}_{1}_{2}.png".format(swnm, nm[4], HPRR_fn)
                    )
                    try:
                        old_target_png = HPRR_dest_dir.joinpath(
                            "HPRR_Tafel_{0}_{1}.png".format(swnm, nm[4])
                        )
                        if old_target_png.is_file():
                            old_target_png.unlink()
                            logger.warning(
                                "HPRR output deleted old target: {0}".format(
                                    old_target_png
                                )
                            )
                    except:
                        logger.warning(
                            "HPRR output delete old target fail: {0}".format(
                                old_target_png
                            )
                        )
                    plt.savefig(swp_target_png, dpi=100, bbox_inches="tight")
                    plt.close()
                else:
                    logger.warning(
                        "HPRR no plot output because TF < 0.11 ({0})".format(
                            swp_target_file
                        )
                    )

        #                rw.plot(x=EvRHE,y='j/E',ylim=(-1,1),label='%s_%s' %(swp,nm))
        HPRR_out = pd.concat([i for i in HPRR_out_lst], sort=False)
        logger.info("Done HPRR analysis of %s: %s" % (SampleID, HPRR_dest_dir))

    except Exception as e:
        print("No successfull HPRR: {0}".format(e))
        logger.error("No successfull HPRR: {0}".format(e))
        HPRR_out = pd.DataFrame([])
    # %%
    return HPRR_out
