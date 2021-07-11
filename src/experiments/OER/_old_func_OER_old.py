def old_OER():
    #    print('OER RPM',1500)
    ###### ====== Exporting all Segments of Scans of .PAR files to OER subfolder ========== #######
    for nm, gr in All_OER.groupby(["PAR_file", "Segment #", "Sweep_Type"]):
        OER_OutFile = OER_dest_dir.joinpath(
            nm[0].stem + "_" + str(int(nm[1])) + "_" + nm[2]
        )
        #        os.path.join(OER_dest_dir,os.path.basename(os.path.splitext(nm[0])[0])+'_'+str(int(nm[1]))+'_'+nm[2]+'.xlsx')
        if not os.path.isfile(OER_OutFile):
            print(nm)
            gr.to_excel(OER_OutFile.with_suffix(".xlsx"))
    ###### ====== Merging CV with Chronoe for Selectivity Exp ========== #######
    for nmBase, grBase in All_OER.groupby(["BaseName"]):
        grF = grBase.groupby(["File"])
        if (
            len(grF) == 2
            and "Chronoamperometry" in grBase.Type_action.unique()
            and "Cyclic Voltammetry (Multiple Cycles)" in grBase.Type_action.unique()
        ):
            print(nmBase)
            OER_CV = grBase.query(
                'Type_action == "Cyclic Voltammetry (Multiple Cycles)"'
            )
            OER_CV_SegMin = int(OER_CV["Segment #"].unique().min())
            OER_Chrono = grBase.query('Type_action == "Chronoamperometry"')
            OER_Chrono = OER_Chrono.loc[OER_Chrono["Segment #"] >= OER_CV_SegMin]
            OER1, OER2 = OER_CV.set_index("Elapsed Time(s)"), OER_Chrono.set_index(
                "Elapsed Time(s)"
            )
            OER_Both = OER1.join(OER2, lsuffix="_disk", rsuffix="_ring", how="outer")
            temp = OER_Both.loc[
                :, ["Segment #_disk", EvRHE + "_disk", "I(A)_disk", "I(A)_ring"]
            ].interpolate()
            OER_Both = temp.rename(
                columns=dict(
                    zip(
                        ["I(A)_disk", "I(A)_ring"],
                        ["I(A)_diskInterP", "I(A)_ringInterP"],
                    )
                )
            )
            OER_Both = OER_Both.assign(
                **{EvRHE + "_disk_diff": OER_Both[EvRHE + "_disk"].diff()}
            )
            #            , 'I(A)_diskInter' : OER_Both['I(A)_disk'].interpolate()})
            OER_Both = OER_Both.assign(
                **{
                    "SweepType": np.where(
                        OER_Both[EvRHE + "_disk_diff"] > 0,
                        "anodic",
                        np.where(OER_Both[EvRHE + "_disk_diff"] < 0, "cathodic", "NA"),
                    ),
                    "FarEff": OER_Both["I(A)_ringInterP"]
                    / (OER_Both["I(A)_ringInterP"] * CollEff),
                }
            )
            for nmBoth, grBoth in OER_Both.groupby(["Segment #_disk", "SweepType"]):
                if nmBoth[1] == "NA":
                    continue
                BothOutFile = OER_dest_dir.joinpath(
                    nmBase + "_Bipot_" + str(int(nmBoth[0])) + "_%s" % nmBoth[1]
                )
                print(BothOutFile)
                #                if not os.path.isfile()
                grBoth.loc[
                    :, [EvRHE + "_disk", "jmAcm-2_ring", "jmAcm-2_disk"]
                ].to_excel(BothOutFile.with_suffix(".xlsx"))
                fig, axOER = plt.subplots(1, figsize=(8, 6))
                axOER.set_title(nmBase + "(%s,%s)" % (nmBoth[0], nmBoth[1]))
                grBoth.plot(x=EvRHE + "_disk", y=["jmAcm-2_disk"], ax=axOER)
                try:
                    grBoth[[EvRHE + "_disk", "jmAcm-2_ring"]].dropna().plot(
                        x=EvRHE + "_disk", y=["jmAcm-2_ring"], ax=axOER
                    )
                except Exception as e:
                    print(e)
                #                plt.show()
                plt.savefig(
                    OER_OutFile.with_suffix(".png"), dpi=300, bbox_inches="tight"
                )
                plt.close()
            OER_Both.loc[:, [EvRHE + "_disk", "jmAcm-2_ring", "jmAcm-2_disk"]].plot(
                x=EvRHE + "_disk", y=["jmAcm-2_ring", "jmAcm-2_disk"]
            )
    return "OER"


# All_OER = All_ovv.query('EXP == "OER"')
# All_HER, HER_ovv_file, dest_dir = (Samples_ovv_cv,HER_ovv_file, Path(HER_ovv_file.Dest_dir.iloc[0]))
