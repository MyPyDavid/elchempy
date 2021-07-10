# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:31:52 2020

@author: User
"""


def makeplot(df):
    grJ_ovv = overall_rows.groupby(by=["SampleID", "DATE", "RPM"])
    Jkin1500_ovv = overall_rows.query('RPM == "1500"')
    plt_lst = [
        "E_half",
        "E_onset",
        "FracH2O2_050",
        "J_diff_lim",
        "Jkin_075",
        "Jkin_080",
        "Jring_050",
    ]
    for l in plt_lst:
        Jkin1500_ovv.plot(x="SampleID", y=l, label=l, kind="bar")
        #        plt.show()
        plt.close()
    t1 = str(
        "F:\\EKTS_CloudStation\\CloudStation\\Experimental data\\Raw_data\\VERSASTAT\\2019-10-Okt\\04.10.2019_0.1MKOH_cell2\\N2_EIS-range_1500rpm_905_DW21_cell2.par"
    )
    FileHelper.FindSampleID.try_find_sampleID(t1)
