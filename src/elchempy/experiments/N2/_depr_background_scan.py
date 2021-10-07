def _depr_background_scan():
    sr_min = N2_CVs.scanrate.min()
    srgrpby = N2_CVs.groupby("scanrate")

    try:
        # lenpreN2 = len(preN2_scan)
        # N2_scan = preN2_scan
        if len(BG_scan) > 2001:
            logger.warning("!! N2 scan more than 2000 data points.. ")
            """
                select scan with multiple steps: # TODO
                take only scans from length 2000, these are "good" measurements
                take the last from multiple segments
                analyze the current and potential values
            """

            #            N2_act_BG_over2000 = Path(N2_dest_dir.parent).joinpath(N2_fn+'_test2000.xlsx')# test
            #            preN2_scan.to_excel(N2_act_BG_over2000)
            #            logger.warning('!! N2 scan more than 2000 data points.. file saved {0} !! len({1})'.format(N2_fn,lenpreN2))
            # preN2_scan = preN2_scan.loc[
            #     preN2_scan.PAR_file == preN2_scan["PAR_file"].unique()[-1], :
            # ]
            sr_grp_min_2000 = pd.concat(
                [gr for n, gr in sr_grp_min.groupby("Segment #") if len(gr) == 2000]
            )
            if sr_grp_min["Segment #"].nunique() > 1:

                sr_grp_min_2000 = pd.concat(
                    [gr for n, gr in sr_grp_min.groupby("Segment #") if len(gr) == 2000]
                )
                # for i in preN2_scan["Segment #"].unique():
                # if len(preN2_scan.loc[preN2_scan["Segment #"] == i]) == 2000:
                # N2_scan = preN2_scan.loc[preN2_scan["Segment #"] == i]
            # elif len(preN2_scan["Segment #"].unique()) == 1:
            # N2_scan = preN2_scan.drop_duplicates()
        elif len(N2_scan) != 2000:
            logger.warning(
                "!! N2 scan not 2000 data points.. {0} !! len({1})".format(
                    N2_fn, len(N2_scan)
                )
            )
            #            print('!! N2 scan less than 2000 data points.. !! %s' %lenpreN2)
            N2_scan = preN2_scan.loc[
                preN2_scan.PAR_file == preN2_scan["PAR_file"].unique()[0], :
            ]
            if len(preN2_scan["Segment #"].unique()) == 1:
                for i in preN2_scan["Segment #"].unique():
                    if len(preN2_scan.loc[preN2_scan["Segment #"] == i]) == 2000:
                        N2_scan = preN2_scan.loc[preN2_scan["Segment #"] == i]

        elif lenpreN2 < 2000:
            logger.warning(
                "!! N2 scan less than 2000 data points..  {0} !! len({1})".format(
                    N2_fn, lenpreN2
                )
            )
            N2_scan = preN2_scan.loc[
                preN2_scan.PAR_file == preN2_scan["PAR_file"].unique()[0], :
            ]
            if len(preN2_scan["Segment #"].unique()) == 1:
                for i in preN2_scan["Segment #"].unique():
                    if len(preN2_scan.loc[preN2_scan["Segment #"] == i]) == 2000:
                        N2_scan = preN2_scan.loc[preN2_scan["Segment #"] == i]
        #                            if len(preN2_scan.loc[preN2_scan['Segment #'] == i]) != 2000:
        #                            N2_scan = pd.DataFrame([])
        #                            print('!! Wrong size for N2 background!!')
        else:
            logger.info(
                "!! N2 scan has 2000 data points..success  {0} !! len({1})".format(
                    N2_fn, lenpreN2
                )
            )
            N2_scan = preN2_scan
    except Exception as e:
        logger.error("N2 scan length problems %s", e)
        N2_scan = preN2_scan
    if N2_scan.empty:
        logger.warning("!! N2 background is Empty!! {0}".format(N2_FileName))
