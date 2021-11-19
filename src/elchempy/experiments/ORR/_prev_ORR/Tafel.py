"""
Created on Fri Nov 12 12:28:01 2021

@author: DW
"""


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
