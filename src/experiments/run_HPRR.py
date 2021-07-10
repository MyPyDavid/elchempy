def HPRR(exp, gr_ovv, ovv):
    ###### === Analyze the HPRR experiments ==== #######
    #        if 'HPRR' in Experiments:
    #        exp,ExpTypes_gr.get_group('HPRR'),ovv
    gr_ovv = ovv.groupby(by="PAR_exp").get_group("HPRR")
    DestDir = Path(gr_ovv.Dest_dir.unique()[0])
    #        gr_ovv = ExpTypes_gr
    HPRR_pars_lst, HPRR_pars = [], []
    for HPRR_file, HPRR_ovv_file in gr_ovv.groupby(by="PAR_file"):
        try:
            HPRR_CVs, HPRR_actions = create_CVs(HPRR_ovv_file)
            Samples_ovv = HPRR_CVs[
                (~HPRR_CVs["SampleID"].str.contains("Pt_ring|Pt-ring"))
                & (~HPRR_CVs["Comment"].str.contains("Pt_ring|Pt-ring"))
            ]
            if len(HPRR_CVs) - len(Samples_ovv) > 0:
                logger.info(
                    "HPRR filtered samples {0}".format(len(HPRR_CVs) - len(Samples_ovv))
                )
            HPRR_file_pars = HPRR_scan(
                Samples_ovv, HPRR_ovv_file, Path(HPRR_ovv_file.Dest_dir.iloc[0])
            )
            HPRR_pars_lst.append(HPRR_file_pars)
        except Exception as e:
            logger.error(
                "HPRR === ERROR in HPRR_scan === {0}\n because {1}".format(HPRR_file, e)
            )
    HPRR_pars = pd.concat([i for i in HPRR_pars_lst])
    HPRR_pars.to_excel(DestDir.joinpath("HPPR_Pars.xlsx"))
    index = pd.DataFrame(
        [
            ["HPRR", DestDir.joinpath("HPPR_Pars.xlsx"), i]
            for i in HPRR_pars.PAR_file.unique()
        ],
        columns=["Type_output", "DestFile", "PAR_file"],
    )
    return index


#        else:
#            print('=== no HPRR in Experimental folder ===')
#            sHPRR = 0
