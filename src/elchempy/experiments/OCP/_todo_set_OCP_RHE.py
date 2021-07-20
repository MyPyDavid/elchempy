''' #TODO rewrite this module '''

import pandas as pd
import numpy as np
import re

from scipy import stats
# import bs4 as BeautifulSoup
from pathlib import Path


from .read_file import (
    get_DF_soup_part,
    collect_action_data,
    collect_data_pars_segment,
    take_hash,
)
# from .CreateCV import check_lin

from file_py_helper.file_functions import FileOperations

import logging

logger = logging.getLogger(__name__)

# set log level

# import statsmodels.api as smf
# def check_lin(colX,colY):
#    data = pd.DataFrame([])
##    colX,colY = OCP_DATA['Elapsed Time(s)'],OCP_DATA['E(V)']
#    data['x'] = colX
#    data['y'] = colY
#    res = smf.ols( 'data.y ~ data.x', data=data).fit()
##    X = sm.add_constant(x), mod = sm.OLS(y,X),res = mod.fit()
#    intercept, slope,rsquared = res.params[0],res.params[1],res.rsquared
#    if np.abs(slope) < 1E-5:
#        check = True
#    else:
#        check = False
#    return intercept, slope,rsquared,check

def check_lin(colX, colY):
    #    data = pd.DataFrame([])
    #    colX,colY = OCP_DATA['Elapsed Time(s)'],OCP_DATA['E(V)']
    #    data['x'] = colX
    #    data['y'] = colY
    #    res = sm.ols( 'data.y ~ data.x', data=data).fit()
    #    X = sm.add_constant(x), mod = sm.OLS(y,X),res = mod.fit()
    #    intercept, slope,rsquared = res.params[0],res.params[1],res.rsquared
    slope, intercept, r_value, p_value, std_err = stats.linregress(colX, colY)
    rsquared = r_value ** 2
    check = True if np.abs(slope) < 1e-5 else False
    return intercept, slope, rsquared, check



def get_RHE_OCP(ovv):
    """takes in the ovv pd.DataFrame and check for the RHE value assignments of each row"""
    #    results_OCP = pd.DataFrame({'OCP Par' : [], 'Mean OCP' : []
    # %%
    results_OCP, action_out = [], []
    hash_f, out = [], []
    RHE_OCPs, RHE_ovv = ovv.query('PAR_exp == "RHE"')["PAR_file"], pd.DataFrame([])
    ovv = ovv.assign(**{"RHE_fn": 0})
    #    if RHE_OCPs.empty:
    try:
        #        RHE_suffixMatch = np.unique([float(i[-1]) for i in ovv.loc[(ovv['EXP_PAR_date_match'] == 'yes') & (ovv['PAR_exp'] == 'N2_act'),'basename'].str.split('_')])
        #                for ovv.loc[:,'basename'].str.split('_'):
        RHE_ovv_matches = []
        for n, rw in ovv.iterrows():

            M_RHE_fn_match = 0
            #                    rw = ovv.iloc[3]
            RHEmatches = 0
            try:
                MatchRHEformat = [
                    re.match(r"([0]{0,1}[0-9]{3}$)", i) for i in rw.basename.split("_")
                ]
            except:
                MatchRHEformat = [False]
                pass
            #            print(MatchRHEformat,rw.basename)
            if any(MatchRHEformat):
                Matches = [i for i in MatchRHEformat if i is not None]
                try:
                    RHEmatches = [
                        float(i.group(0))
                        for i in Matches
                        if float(i.group(0)) < 1000 and float(i.group(0)) > 50
                    ]
                    if "N2" in rw.basename:
                        RHEmatches = [
                            i for i in RHEmatches if i != 300 and i != 10 and i != 100
                        ]
                    elif rw.pH > 7:
                        RHEmatches = [
                            i
                            for i in RHEmatches
                            if i != 300 and i != 10 and i != 100 and i > 500
                        ]
                    elif rw.pH < 4:
                        RHEmatches = [
                            i
                            for i in RHEmatches
                            if i != 300 and i != 10 and i != 100 and i < 400
                        ]
                    else:
                        RHEmatches = [
                            i for i in RHEmatches if i != 300 and i != 10 and i != 100
                        ]
                #                    print(rw.basename,RHEmatches)
                except Exception as e:
                    RHEmatches = [i.group(0) for i in Matches]
                    logger.warning(
                        "RHE OCP filename float casting error (%s), %s \n %s"
                        % (RHEmatches, rw.basename, e)
                    )
                if len(RHEmatches) == 1:
                    try:
                        M_RHE_fn_match = float(RHEmatches[0])
                    except Exception as e:
                        logger.warning(
                            "RHE OCP filename float error ({0}), {1}\n because {2}"
                            % (RHEmatches[0], rw.basename, e)
                        )
                        # maybe change something at string ...
                        M_RHE_fn_match = RHEmatches[0]
                elif len(RHEmatches) > 1:
                    logger.warning(
                        "RHE OCP: more than 1 option for RHE value in filename %s. %s"
                        % (rw.basename, Matches)
                    )
                    if rw.pH < 5:
                        try:
                            M_RHE_fn_match = [
                                float(i)
                                for i in RHEmatches
                                if float(i) < 350 and float(i) > 50
                            ][-1]
                        except:
                            M_RHE_fn_match = [float(i) for i in RHEmatches][0]
                            logger.warning(
                                "RHE OCP: acidic (pH %s) error in filename %s. %s \n First value taken: %s mV"
                                % (rw.pH, rw.basename, RHEmatches, M_RHE_fn_match)
                            )
                    elif rw.pH > 5:
                        try:
                            M_RHE_fn_match = [
                                float(i)
                                for i in RHEmatches
                                if float(i) < 1000 and float(i) > 350
                            ][-1]
                        except:
                            M_RHE_fn_match = [float(i) for i in RHEmatches][0]
                            logger.warning(
                                "RHE OCP: alkaline (pH %s) error in filename %s. %s \n First value taken: %s mV"
                                % (rw.pH, rw.basename, RHEmatches, M_RHE_fn_match)
                            )
            else:
                #                print('No matching RHE number in filename: %s'%rw.basename)
                M_RHE_fn_match = 0
            RHE_ovv_matches.append(M_RHE_fn_match)

            if RHE_ovv_matches == []:
                RHE_mean = [0]
            elif np.std(RHE_ovv_matches) == 0:
                RHE_mean = [np.mean(RHE_ovv_matches)]
            else:
                if len(np.unique(RHE_ovv_matches)) == 1:
                    RHE_mean = [np.mean(RHE_ovv_matches)]
                elif len(np.unique(RHE_ovv_matches)) != 1:
                    #                    print('Multiple from RHE from files %s' %[np.unique(RHE_ovv_matches,return_counts=True)])
                    RHE_mean = np.unique(RHE_ovv_matches, return_counts=True)[0]

        try:
            #            print('Error assinging RHE values to OVV in exp dir: %s \n %s'%(rw.EXP_dir,e))
            #            ovv = ovv.assign(**{'RHE_fn' : RHE_mean[0]})
            if len(RHE_ovv_matches) == len(ovv):
                ovv = ovv.assign(**{"RHE_fn": RHE_ovv_matches})
            else:
                logger.warning(
                    " RHE OCP Error assinging RHE values to OVV in exp dir: %s \n, took RHE_mean"
                    % (rw.EXP_dir)
                )
                ovv = ovv.assign(**{"RHE_fn": RHE_mean[0]})
        #            assign(**{'RHE_fn' : })
        #                print('Assinging RHE values (%.0f) to %s \n OVV in exp dir: %s'%(RHE_mean[0],rw.basename,rw.EXP_dir))
        except Exception as e:
            logger.warning(
                " RHE OCP Error assinging RHE values to OVV in exp dir: %s \n %s"
                % (rw.EXP_dir, e)
            )
            pass
    #                if len(RHE_suffixMatch) == 1:
    #                    RHE_mean = [RHE_suffixMatch[0]*1E-3]
    #                else:
    #                    print('RHE suffixes of N2 scans do not match! %s'%RHE_suffixMatch)
    #                    RHE_mean = [(', ').join([str(int(i)) for i in RHE_suffixMatch])]
    #                except Exception as e:
    #                    print('RHE matching error: %s' %e)
    except Exception as e:
        if not ovv.empty:
            logger.warning(
                "Overall dir errror assinging RHE values to OVV in exp dir: %s \n %s"
                % (ovv.iloc[0].EXP_dir, e)
            )

        elif ovv.empty:
            logger.warning("Overall dir errror OVV = EMPTY in exp dir: \n %s" % (e))
        #        ovv['RHE_fn'] = 0
        RHE_mean = [0]

def OCP_RHE_measurement()
    if not RHE_OCPs.empty:

        try:
            RHE_elec = []
            for RHE_OCP in RHE_OCPs:
                RHE_PF = Path(RHE_OCP)
                with open(RHE_PF) as OCP_file:
                    OCP_read = OCP_file.read()
                OCP_soup_par = BeautifulSoup(OCP_read, "lxml")
                hash_f.append(
                    [FileOperations.take_hash(RHE_OCP)]
                )  # FIX ME define local hash func
                act0_DF = get_DF_soup_part(OCP_soup_par.action0, RHE_PF.name)
                RHE_elec = [i for i in Path(RHE_OCP).name.split("_") if "RHE" in i]

                Ag_elec = [
                    i for i in Path(RHE_OCP).name.split("_") if "Ag" in i or "Hg" in i
                ]
                try:
                    act1_dt = get_DF_soup_part(OCP_soup_par.action1, RHE_PF.name)
                    if act1_dt["Name"].values[0] == "Open Circuit":
                        OCP_DATA = collect_data_pars_segment(OCP_soup_par.segment1)
                        if OCP_DATA.empty:
                            continue
                        startlin = int(len(OCP_DATA) * 0.2)
                        check_lin_res = check_lin(
                            OCP_DATA.iloc[startlin::]["Elapsed Time(s)"],
                            OCP_DATA.iloc[startlin::]["E(V)"],
                        )
                        if np.abs(check_lin_res[1]) < 5e-05:
                            msrmnt = "stable"
                        #                            print('RHE_OCP: %.3f' %OCP_DATA['E(V)'][10::].mean())
                        elif np.abs(check_lin_res[1]) > 5e-05:
                            #                            print('unstable')
                            msrmnt = "unstable"
                        OCP_DATA.ActionId.unique()
                        results_OCP.append(
                            {
                                "RHE_file": RHE_OCP,
                                "RHE_file_hash": hash_f,
                                "Erhe_mean": OCP_DATA["E(V)"][startlin::].mean(),
                                "Lin_fit": check_lin_res,
                                "Stability": msrmnt,
                                "RHE_elec": RHE_elec[0],
                                "Ref_elec": Ag_elec[0],
                                "ActionID": OCP_DATA.ActionId.unique()[0],
                            }
                        )
                        #    'Comment' : act0_DF['Comment'].values[0]
                    else:
                        #                        print('Not Open Circuit %s'%(RHE_OCP))
                        logger.warning("OCP RHE Not Open Circuit %s" % (RHE_OCP))
                #                continue
                except Exception as e:
                    #                    print('RHE fail:%s' %e,RHE_OCP)
                    logger.error("OCP RHE fail: {0},{1}".format(RHE_OCP, e))
                    continue
            RHE_ovv = pd.DataFrame(data=results_OCP)
            if RHE_ovv.empty:
                logger.warning("OCP RHE Skipped because of RHE_ovv is empty")

            #       continue
            elif len(RHE_ovv.loc[RHE_ovv["Stability"] == "stable", "Erhe_mean"]) == 1:
                RHE_mean = float(
                    RHE_ovv.loc[RHE_ovv["Stability"] == "stable", "Erhe_mean"].values[0]
                )
                #                print('Used %.0f mV as RHE potential' %(RHE_mean*1000))
                logger.info("Used %.0f mV as RHE potential" % (RHE_mean * 1000))
            elif len(RHE_ovv.loc[RHE_ovv["Stability"] == "stable", "Erhe_mean"]) > 1:
                grRHE = RHE_ovv.loc[
                    RHE_ovv["Stability"] == "stable",
                    ["Erhe_mean", "RHE_elec", "Ref_elec"],
                ].groupby(by="RHE_elec")
                if len(grRHE.groups) > 1:
                    grRHEmax = grRHE.aggregate([np.mean, len])["Erhe_mean"][
                        "len"
                    ].idxmax()
                    RHE_mean = grRHE.get_group(grRHEmax)["Erhe_mean"].mean()
                else:
                    RHE_mean = RHE_ovv.loc[
                        RHE_ovv["Stability"] == "stable", "Erhe_mean"
                    ].mean()
                MultipleRHE_OCPSs = [
                    (Path(i[1].RHE_file).name, i[1].Erhe_mean)
                    for i in RHE_ovv.loc[RHE_ovv["Stability"] == "stable"].iterrows()
                ]
                #                print('Muiltpe OCPs available: %s \n Used %.0f mV as RHE potential' %(MultipleRHE_OCPSs ,RHE_mean*1000))
                logger.warning(
                    "OCP RHE: Multiple OCPs available: %s \n Used %.0f mV as RHE potential"
                    % (MultipleRHE_OCPSs, RHE_mean * 1000)
                )
            elif len(RHE_ovv.loc[RHE_ovv["Stability"] == "stable", "Erhe_mean"]) == 0:
                #                print('CHECK RHE VALUE!! Mean taken from unstable OCP measurement')
                logger.warning(
                    "OCP RHE: CHECK RHE VALUE!! Mean taken from unstable OCP measurement"
                )
                RHE_mean = RHE_ovv["Erhe_mean"].mean()

        except Exception as e:
            logger.error(
                "OCP RHE Error RHE OCP files in exp dir: %s \n %s"
                % (ovv.iloc[0].EXP_dir, e)
            )
    else:
        logger.warning("No RHE_OCP files in Ovv==empty!")
        RHE_mean = [0]
    if RHE_ovv.empty and len(RHE_mean) >= 1:
        logger.warning(f"set_OCP_RHE End problem RHE_mean:{RHE_mean}")
        #        ovv['RHE_mean'] = RHE_mean[0]
        ovv = ovv.assign(**{"RHE_mean": RHE_mean[0]})
    else:
        #        ovv['RHE_mean'] = RHE_mean
        ovv = ovv.assign(**{"RHE_mean": RHE_mean})

    # %%
    return ovv, RHE_ovv, RHE_mean
