"""
Created on Sun Jul 11 11:26:41 2021

@author: DW
"""


def old_HER_run(ovv_exp_grp, **HER_kwargs):
    ####### === Analyze the HER experiments ==== #######
    #        DestDir = Path(gr_ovv.Dest_dir.unique()[0])
    #        print('OER skipped')
    #        return pd.DataFrame()
    gr_ovv = ovv.groupby(by="PAR_exp").get_group("HER")
    DestDir = Path(gr_ovv.Dest_dir.unique()[0])
    #        gr_ovv = ExpTypes_gr
    HER_index_lst, HER_pars, faillst = [], [], []
    for HER_file, HER_ovv_file in gr_ovv.groupby(by="PAR_file"):
        try:
            HER_CVs, HER_actions = create_CVs(HER_ovv_file)
            if not HER_CVs.empty:
                Samples_ovv_cv = HER_CVs[
                    (~HER_CVs["SampleID"].str.contains("Pt_ring|Pt-ring"))
                ]
                if len(HER_CVs) - len(Samples_ovv_cv) > 0:
                    logger.info(
                        "HER filtered samples from Pt-ring files {0}".format(
                            len(HER_CVs) - len(Samples_ovv)
                        )
                    )
                if not Samples_ovv_cv.empty:
                    HER_file_index = HER_scan(
                        Samples_ovv_cv,
                        HER_ovv_file,
                        Path(HER_ovv_file.Dest_dir.iloc[0]),
                    )
                    HER_index_lst.append(HER_file_index)
                else:
                    logger.error(
                        "HER === ERROR HER_scan Samples ovv empty === {0}\n because {1}".format(
                            HER_file
                        )
                    )
                    faillst.append([HER_file, "Samples_ovv_cv.empty"])
            else:
                logger.error(
                    "HER === ERROR HER_CVs empty === {0}\n because {1}".format(HER_file)
                )
                faillst.append([OER_file, "OER_CV.empty"])
        except Exception as e:
            logger.error(
                "HER === ERROR in HER_scan === {0}\n because {1}".format(HER_file, e)
            )
    HER_index = pd.concat([i for i in HER_index_lst], ignore_index=True)
    HER_pars_target = FolderOps.FileOperations.CompareHashDFexport(
        HER_index, DestDir.joinpath("HER_Pars_index.xlsx")
    )
    return HER_index


#    @staticmethod
#    def HER(exp,gr_ovv,ovv):

##        if 'HER' in Experiments:
##            print('HER')
#        HER_CVs,HER_info = pd.DataFrame(), pd.DataFrame()
#        HER_ovv = ovv.loc[((ovv['PAR_exp'] == "HER") & (ovv.basename.str.contains('HER'))),:]
#        if HER_ovv.empty:
#            HER_ovv = ovv.loc[((ovv['PAR_exp'] == "EIS") & (ovv.basename.str.contains('HER'))),:]
##                    HER_CVs,HER_info = create_CVs(HER_ovv,PathDB,True)
#        HER_dest_dir = dest_dir.joinpath('HER')
#        HER_dest_dir.mkdir(parents=True,exist_ok=True)
#        HER_out_fn = HER_dest_dir.joinpath('HER_Pars.xlsx')
#        HER_pars = HER_calc(HER_ovv,HER_out_fn,PathDB)
#                if not HER_CVs.empty:
#                    HER_pars = H$R_scan(HER_CVs,dest_dir)
#                else:
#                    print('OER failed: %s'%exp_dir)
#                OER_pars,OER_action = OER_scan(create_CVs(ovv.query('PAR_exp == "OER"'),PathDB,False),dest_dir)

#    else:
#        plt.close()
#        print('No ORR')
#        sORR = 0
