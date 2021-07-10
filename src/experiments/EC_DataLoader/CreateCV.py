from pathlib import Path
from itertools import repeat
from collections import namedtuple
import multiprocessing
import os
import pandas as pd
import numpy as np

# import statsmodels.api as sm
from scipy import stats
from bs4 import BeautifulSoup
import re
import matplotlib.pyplot as plt

from .RHE_assignment import RHE_potential_assignment
from .read_file import (
    get_DF_soup_part,
    collect_data_pars_segment,
    take_hash,
    DAC_V_to_RPM,
    add_RRDE,
)

if __name__ == "__main__":
    from RHE_assignment import RHE_potential_assignment
    from read_file import (
        get_DF_soup_part,
        collect_data_pars_segment,
        take_hash,
        DAC_V_to_RPM,
        add_RRDE,
    )


#    from ..EC_conditions.electrode import WE_SA_collection_eff


# SORT OF GLOBAL CONSTANTS
EvRHE = "E_AppV_RHE"

CV_output_row = namedtuple("CV_data_from_pf", "data actions")


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


# def get_eis_data(ovv,**kwargs):
#    logger


def WE_SA_collection_eff(TYPE="PINE"):
    coll_eff = []
    if TYPE == "ALS":
        r1, r2, r3, coll_eff = 0.1 * 4 * 0.5, 0.1 * 5 * 0.5, 0.1 * 7 * 0.5, []
        SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
    if TYPE == "PINE":
        r1, r2, r3 = 0.1 * 5.5 * 0.5, 0.1 * 6.5 * 0.5, 0.1 * 8.5 * 0.5
        coll_eff = 0.38
        SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
    if TYPE == "PINE-ring":
        r1, r2, r3 = 0.1 * 5.5 * 0.5, 0.1 * 6.5 * 0.5, 0.1 * 8.5 * 0.5
        coll_eff = 0.38
        SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
    #        SA = np.pi*(r3**2-r2**2)
    if coll_eff == []:
        a, b = (r2 / r1) ** 3 - 1, (r3 / r1) ** 3 - (r2 / r1) ** 3
        c = a / b
        coll_eff = (
            1
            - r3 ** 2
            + b ** (2 / 3)
            - coll_func(c)
            - b ** (2 / 3) * coll_func(a)
            + r3 ** 2 * coll_func(c * r3 ** 3)
        )
    "r1 = disk, r2 = ring ID, r3 = ring OD in mm"
    #    print('%s According to manufacturer: disk(dia:%.2f cm   %.4f cm2), ring (%.4f cm2)' %(TYPE,r1*2,SAdisk,SAring))
    return {
        "Electrode Type": TYPE,
        "CollEff": coll_eff,
        "Disk_cm2": np.round(SAdisk, 4),
        "Ring_cm2": np.round(SAring, 4),
    }


def _drop_cols_from_data_segment():
    # Other optional columns: ['E Imag', 'ADC Sync Input(V)','I Imag', 'I Real']
    _drop_cols = ["E2 Imag", "E2 Real", "E2 Status", "E2(V)", "Z2 Imag", "Z2 Real"]
    return _drop_cols


def read_data_CV_row(ovv_row, **kwargs):
    """Takes a row from the EC_index and reads in the the data from the PAR_file segments
    returns a namedtuple: r.data and r.actions"""
    CVrow = ovv_row.iloc[0]
    n = CVrow.name
    global EvRHE
    OutPutRaw = kwargs.get("OutPutRaw", False)
    # if 'Series' in str(type(CVrow)):
    #     # RHE_potential = RHE_potential_assignment(CVrow.to_frame().T)
    #     CVrow_frame = CVrow_to
    RHE_potential = RHE_potential_assignment(CVrow)
    if CVrow.Electrode == "Pt_ring":
        #            print('Pt-ring electrode: %s'%CVrow.basename)
        WE_surface_area = WE_SA_collection_eff("PINE")["Ring_cm2"]
    else:
        #            print('Disk electrode: %s'%CVrow.basename)
        WE_surface_area = WE_SA_collection_eff("PINE")["Disk_cm2"]
    #        WE_surface_area = np.pi*(0.55/2)**2
    N2_file = Path(CVrow["PAR_file"])
    if N2_file.is_file() == False:
        logger.error(f"Create CV skipped file because {str(N2_file)} not exists")
    #            #continue
    DestDir = Path(CVrow["Dest_dir"])
    #        print(N2_file)
    #        N2_BaseFile = '_'.join([i for i in  os.path.basename(N2_file).split('_')[0:-2]])
    N2_BaseFile = N2_file.stem
    DestFile = DestDir.joinpath(N2_BaseFile)

    #            continue
    #               out_DATA_CVs=pd.DataFrame([])
    with open(N2_file) as N2f:
        N2_read = N2f.read()
        N2_soup_par = BeautifulSoup(N2_read, "html.parser")
    #    hash_f = take_hash(N2_file)
    N2f_bn = N2_file.name
    #    act0_DF = get_DF_soup_part(N2_soup_par.action0 ,hash_f)
    #    instr_DF = get_DF_soup_part(N2_soup_par.instrument ,hash_f)
    #    date_df = get_DF_soup_part(N2_soup_par.experiment ,hash_f)
    #    names_soup = [tag.name for tag in N2_soup_par.find_all(True)]
    #    try:
    #        date_df['DTC'] = pd.to_datetime(date_df['DateAcquired' ] +' ' +date_df['TimeAcquired'])
    #    except Exception as e:
    #        logger.warning('Create CV DTC failed: {0}'.format(e))
    #                continue
    _new_info = {"DestFile": DestFile, "RHE_OCP": np.abs(RHE_potential)}

    segment1 = collect_data_pars_segment(N2_soup_par.segment1)
    segment1 = segment1.assign(**{**CVrow.to_dict(), **_new_info})

    segment1.drop(_drop_cols_from_data_segment(), axis=1, inplace=True)

    if segment1.empty:
        logger.warning(
            "Create CV should skip because segment1 is empty: %s"
            % N2_file.parent.joinpath(N2_BaseFile)
        )
    #            continue
    elif len(segment1) < 2000 or len(segment1["Segment #"].unique()) > 1000:
        if any([i for i in segment1.ActionId.unique() if i == 21]):
            pass
        #                print('EIS measurements')
        else:
            try:
                fig, ax = plt.subplots()
                segment1.plot(x="Elapsed Time(s)", y="I(A)", ax=ax)
                DestDirSegs = DestDir.joinpath("Raw_PAR_data_errors")
                DestDirSegs.mkdir(parents=True, exist_ok=True)
                plt.savefig(DestDirSegs.joinpath(N2_BaseFile).with_suffix(".png"))
                plt.close()
                logger.warning(
                    "Create CV Skipped and plotted for too short/long (%s) data: %s"
                    % (len(segment1), N2_file)
                )
            except Exception as E:
                logger.warning(
                    "Create CV Skipped and NOT plotted for too short/long (%s) data: %s.\n %s"
                    % (len(segment1), N2_file, E)
                )
    #                continue
    #    segment1['RHE_OCP'] = np.abs(RHE_potential)
    if "Type" in segment1.columns:
        segment1 = segment1.rename(columns={"Type": "Type_action"})

    EC_important_cols = {
        EvRHE: segment1["E(V)"] + segment1["RHE_OCP"],
        "E_programmed_modulation": segment1["E Applied(V)"] + segment1["RHE_OCP"],
        "j A/cm2": segment1["I(A)"] / WE_surface_area,
        "jmAcm-2": 1000 * segment1["I(A)"] / WE_surface_area,
        "ScanRate_calc": np.round(
            np.abs(segment1["E(V)"].diff() / segment1["Elapsed Time(s)"].diff()), 3
        ),
        "ScanRate_calc_progr": np.round(
            (segment1["E Applied(V)"].diff() / segment1["Elapsed Time(s)"].diff()), 3
        ),
    }
    segment1 = segment1.assign(**EC_important_cols)
    #    _meta_info_columns = {'Gas': CVrow['Gas'], 'EXP': CVrow['PAR_exp'], 'Electrode': CVrow.Electrode,
    #                  'Comment': actComment, 'Measured_OCP': float(act0_DF['Measured Open Circuit'][0].split()[0]),
    #                  'pH': CVrow.pH, 'Electrolyte': CVrow.Electrolyte,
    #                  'Loading_name' : CVrow.Loading_name, 'Loading_cm2' : CVrow.Loading_cm2}

    #                        ,'Analysis_date' : datetime.now())}
    #    segment1 = segment1.assign(**extra_cols)
    if "E Applied(V)" in segment1.columns:
        segment1 = segment1.drop(columns="E Applied(V)")
    #    seg_colls = {'ScanRate_calc': np.abs(segment1['E(V)'].diff() / segment1['Elapsed Time(s)'].diff()),
    #                 'ScanRate_calc_progr': (segment1['E_programmed_modulation'].diff() / segment1['Elapsed Time(s)'].diff()),
    #                }
    #    segment1 = segment1.assign(**seg_colls)

    # TEST checking Segement data
    #        segment1.plot(x = 'E(V)', y = 'Elapsed Time(s)')

    #        segment1 = bb
    segment1 = segment1.assign(
        **{
            "Sweep_Type": np.where(
                segment1.ScanRate_calc_progr > 0,
                "anodic",
                np.where(
                    segment1.ScanRate_calc_progr < 0,
                    "cathodic",
                    np.where(segment1.ScanRate_calc_progr == 0, "chrono", None),
                ),
            )
        }
    )
    segment1.Sweep_Type = segment1.Sweep_Type.fillna(method="backfill")
    #    _seg_add_sweep = []
    for segn, seggrp in segment1.groupby(["Segment #"]):
        segn, seggrp
        if len(seggrp) == 2000:
            for _r0, _r1 in [(0, 1000), (1000, 2000)]:
                _r_swp_max = (
                    seggrp.loc[seggrp.iloc[_r0:_r1].index, "Sweep_Type"]
                    .value_counts()
                    .idxmax()
                )
                segment1.loc[seggrp.iloc[_r0:_r1].index, "Sweep_Type"] = _r_swp_max

    #                seggrp.iloc[1000:].loc[:,'Sweep_Type'] == 'anodic'
    #            _seg_add_sweep.append(seggrp)
    #        segment1 = pd.concat(_seg_add_sweep,sort=True)

    #        segment1 = segment1.assign(**{'Sweep_Type': np.where(segment1.ScanRate_calc_progr > 0, 'anodic',
    #        segment1.Sweep_Type = segment1.Sweep_Type.fillna(method='backfill')
    actions = [
        i.name
        for i in N2_soup_par.find_all(name=re.compile("action."))
        if len(i.name) < 9
    ][1::]

    SRunq = segment1["ScanRate_calc"].round(3).unique()
    seggies = segment1["Segment #"].unique()
    Seggr = segment1.groupby(by="Segment #")
    #        ===== Assing Stuff to each Seg of Segment1 ! ===
    ActionId_Reference = {
        "Cyclic Voltammetry (Multiple Cycles)": 38,
        "Chronoamperometry": 3,
        "Unknown": 0,
        "Potentiostatic EIS": 21,
        "Cyclic Voltammetry": 14,
    }

    for SegN, SegGr in Seggr:
        mean_sr = round(
            np.abs(
                segment1.loc[segment1["Segment #"] == SegN, "ScanRate_calc_progr"]
            ).mean(),
            2,
        )
        try:
            segment1.loc[segment1["Segment #"] == SegN, "Scanrate"] = mean_sr
        #                print('Seg %s, sr %s'%(SegN,mean_sr))
        #                segment1.loc[segment1['Segment #'] == SegN,'RPM'] = RPM
        except:
            segment1.loc[segment1["Segment #"] == SegN, "Scanrate"] = 0
        try:
            segAcId = int(
                segment1.loc[segment1["Segment #"] == SegN, "ActionId"].unique()[0]
            )
        except Exception as e:
            segAcId = 0
            logger.warning("Create CV assignt Type unkown %s %s" % (SegN, segAcId))
        try:
            segment1.loc[segment1["Segment #"] == SegN, "Type_action"] = [
                i for i in ActionId_Reference if ActionId_Reference[i] == segAcId
            ][0]
        except Exception as e:
            logger.warning(
                "Create CV assignt Type unknown %s %s in %s" % (SegN, segAcId, N2f_bn)
            )
    #        ===== Assing Stuff to each Seg of Segment1 ! ===
    #            O2_act = O2_act.assign(Sweep_Type=np.where(O2_act[EvRHE].diff() > 0, 'anodic', np.where(O2_act[EvRHE].diff() < 0, 'cathodic', 'NA')))
    action1_name = get_DF_soup_part(N2_soup_par.action1, N2_soup_par.action1.name)[
        "Name"
    ][0]
    #            actions = [add_RRDE(action1_name,i.name,N2_file) for i in N2_soup_par.find_all(name=re.compile('action.'))][1::]
    #            actions = [i.name for i in N2_soup_par.find_all(name=re.compile('action.')) if len (i.name) < 9][1::]
    action_out = pd.DataFrame([])
    #            gr_segs = segment1.groupby(by='Segment #')
    #        act_list = [[0]], tot_seg, fts\ sts = 0, 0, 0\ lst = []\ tot_points = 0\ActieList = []
    #            for actie in actions[::]:
    DAC_V = 0
    #        print(lst)
    tot_actie = 0
    #        ActionOvv = [(i,get_DF_soup_part(N2_soup_par.find_all(name=i)[0],i)['Name'].values[0],
    #          get_DF_soup_part(N2_soup_par.find_all(name=i)[0],i)['Segments'].values[0]) for i in actions]
    all_action_dfs = pd.concat(
        [get_DF_soup_part(N2_soup_par.find_all(name=i)[0], i) for i in actions],
        sort=False,
    )
    List_action_out = []
    for num, actie in enumerate(actions):
        actie_df = get_DF_soup_part(N2_soup_par.find_all(name=actie)[0], actie)
        if num == 0:
            RPM = 0
        try:
            actie_df.Segments = actie_df.Segments.astype(float)
        except:
            #                print('No segmetns:',actie_df['Name'].values[0] )
            actie_df["Segments"] = 0

        actie_name = actie_df["Name"].values[0]
        actie_df["Type_action"] = actie_name
        #            print(actie_name,int(actie_df.Segments.values[0]))
        if actie_name == "Chronoamperometry":
            actie_df["Cycles"] = int(actie_df.Segments.values[0])
            actie_df["Scan Rate (V/s)"] = 0.0

        elif actie_name == "Cyclic Voltammetry":
            actie_df["Cycles"] = int(actie_df.Segments.values[0])

        elif actie_name == "Cyclic Voltammetry (Multiple Cycles)":
            actie_df["Cycles"] = int(actie_df.Segments.values[0])

        elif actie_name == "RRDE":
            add_actie = add_RRDE(action1_name, actie, N2_file)
            add_actie_df = get_DF_soup_part(
                N2_soup_par.find_all(name=add_actie)[0], actie
            )
            actie_df = actie_df.join(add_actie_df, rsuffix="_add")
            #                    actie_df = pd.concat([actie_df,add_actie_df],axis=1)
            DAC_V = float(actie_df["RotationRate"][0])
            actie_df["Scan Rate (V/s)"] = 0.01
            actie_df["Cycles"] = 1
        elif actie_name == "Potentiostatic EIS":
            actie_df["Cycles"] = 1
            if not "DAC Control" in all_action_dfs.Name.to_list():
                if (
                    "rpm" in N2_BaseFile or str(1500) in N2_BaseFile
                ) and not "rpm-range" in N2_BaseFile:
                    try:
                        RPM = N2_BaseFile.split("rpm")[0].split("_")[-1]
                    except:
                        RPM = 1500.0
                else:
                    RPM = 0.0
        #                    continue
        #                actie_totpoints = int(actie_df['Total Points'][0] )
        elif actie_name == "Time Delay":
            actie_segs = 0
            actie_df["Cycles"] = 0
            actie_df["Segments"] = 0

        ### CONVERTS DAC control to RPM_DAC
        if actie_name == "DAC Control":
            actie_segs = 0
            actie_df["Cycles"] = 0
            actie_df["Segments"] = 0
            try:
                DAC_V = float(actie_df["Potential (V),Value"][0])
                RPM = DAC_V_to_RPM(DAC_V, test_PAR_file=N2_file)
            #                    print(f'DAC_V : {DAC_V}, RPM: {RPM}') #DEBUG
            except Exception as e:
                RPM = 0
                logger.error(" Create CV RPM fail {0}".format(e))

        #                print(DAC_V)
        #                    print('%s: %s V'%(actie_name,DAC_V))
        #                continue
        elif actie_name == "Open Circuit":
            actie_df["Cycles"] = 1

        elif actie_name == "Loop #1":
            actie_df["Cycles"] = 0
        actie_df["RPM"] = RPM
        #            else:
        #                print('ACTIE NAME UNKNOWN: %s for %s'%(actie_name,N2_BaseFile))
        try:
            sr = float(actie_df["Scan Rate (V/s)"][0])
        except Exception as e:
            #                    print(e,actie_df.columns)
            sr = "NAN"
        try:
            actie_segs = int(actie_df["Cycles"][0])
        except:
            actie_segs = 0
        tot_actie = actie_segs + tot_actie
        if actie_segs > 1:
            extra = 0
            if tot_actie == 1:
                extra = 1
                seglink = list(range(0 + extra, tot_actie + 1))
        #                if lst == []:
        #                    extra = 0
        #                    ListSegs = list(range(ActieList[-1][1]+extra,tot_actie+1))
        #                    seglink = list(range(0+extra,tot_actie+1))
        #                elif lst[-1][0] == 0:
        #                    extra = 1
        #                    ListSegs = list(range(ActieList[-1][1]+extra,tot_actie+1))
        #                    seglink = list(range(lst[-1][0]+extra,tot_actie+1))
        if actie_segs == 1 and tot_actie == 0:
            ListSegs = [0]
            seglink = [0]
        elif actie_segs == 1 and tot_actie != 0:
            ListSegs = [tot_actie]
            seglink = [tot_actie]
            last = 0
        #            print(actie_name,actie,num,seglink)
        # +++ TEST +++
        #                print(actie_name,actie_segs,tot_actie,seglink,last,sr)
        # +++ TEST +++
        #            actie_df['SegmentsNew'] = [seglink]
        #            ActionId_Reference = {'Cyclic Voltammetry (Multiple Cycles)' : 38,'Chronoamperometry' : 3}
        #            try:
        #                actie_df['ActionId'] = ActionId_Reference[actie_name]
        #            except:
        #                actie_df['ActionId'] = 0
        #            lst.append([num,actie,actie_name,actie_segs,tot_actie,seglink,DAC_V,sr])
        #            print([num,actie,actie_name,actie_segs,tot_actie,seglink,DAC_V,sr])
        #            set(action_out.columns) - set(actie_df.columns)
        #            action_out = action_out.append(actie_df)
        #            if not action_out.empty:
        List_action_out.append(actie_df)

    action_out = pd.concat([i for i in List_action_out], sort=False)
    #        action_out = pd.concat([action_out,actie_df],axis=0,sort=False)
    #            else:
    #                action_out = pd.concat([action_out,actie_df],axis=0)
    action_out["Cycles"] = action_out["Cycles"].astype(float)
    action_out["PAR_file"] = N2_file

    actionsDFsegs = action_out.loc[action_out.Segments > 0]
    act_seggies_diff = len(seggies) - actionsDFsegs.Segments.sum()
    out_join = []
    if act_seggies_diff == 0:
        seg_counter = -1
        for n, act in actionsDFsegs.iterrows():
            for s in np.arange(act.Segments):
                seg_counter = seg_counter + 1
                out_join.append([n, seg_counter, act.RPM, act["Scan Rate (V/s)"]])
        #                    print(n,s,seg_counter,act.RPM,act['Scan Rate (V/s)'])
        act_seg_match = pd.DataFrame(
            out_join, columns=["Name", "Segment #", "RPM_DAC", "Scan Rate (V/s)"]
        )
        segment1 = pd.merge(segment1, act_seg_match, on="Segment #", how="left")
        logger.info("Segments Actions matched and updated: {0}".format(N2_BaseFile))
    #            for sg in segment1['Segment #'].unique():
    #                out.append(act,sg)
    else:
        #            print('Segments Data and Actions do not match: {0}'.format(N2_BaseFile))
        logger.warning(
            "Segments Data and Actions do not match: {0}, diff {1}".format(
                N2_BaseFile, act_seggies_diff
            )
        )
        segment1 = segment1.assign(
            **{"Name": "NaN", "RPM_DAC": 0, "Scan Rate (V/s)": 0}
        )
    #        seg0 = 0
    segment1["Scan Rate (V/s)"] = segment1["Scan Rate (V/s)"].astype(float)
    segment1["RPM_DAC"] = segment1["RPM_DAC"].astype(float)
    segment1 = segment1.assign(**{"Segm": segment1["Segment #"]})

    if len(segment1["Segment #"].unique()) > 1000 and OutPutRaw == True:
        logger.warning(
            "OutPutRaw is cancelled because of large number of Segments (%s) in: %s"
            % (len(segment1["Segment #"].unique()), Path(CVrow.PAR_file).name)
        )
        OutPutRaw == False
    if OutPutRaw == True:
        DestDirSegs = DestDir.joinpath("Raw_PAR_data_errors")
        DestDirSegs.mkdir(parents=True, exist_ok=True)

        logger.warning("OutPutRaw is true in: %s" % (Path(CVrow.PAR_file).name))
        #            dlo.to_hdf(PathDB,DestFile.joinpath('ActionOut').as_posix(),table=True,mode='a')
        #            segment1.to_hdf(PathDB,DestFile.joinpath('Data').as_posix(),table=False,mode='a')
        FileActSegs = pd.DataFrame()
        FileInfoOut = pd.DataFrame(columns=segment1.columns)
        for seg, Sg in segment1.groupby("Segment #"):
            if "Series" in str(type(Sg)):
                logger.warning("Create CV Skipped series segment,string")
            #                    continue
            if len(Sg) < 2:
                logger.warning("Create CV Skipped series segment, length %s" % seg)
                continue
            try:
                Sg = Sg.set_index("Point #")
                Sg["Segment #"] = Sg["Segment #"].astype(float)
            except Exception as e:
                continue
                logger.warning("Create CV Segment set index col, problem")
            Sg_info = [
                (i, Sg[i].unique()) for i in Sg.columns if len(Sg[i].unique()) < 30
            ]
            Sg_info_Out = pd.DataFrame(columns=[i[0] for i in Sg_info])
            Sg_info_Out.loc[seg, [i[0] for i in Sg_info]] = [i[1] for i in Sg_info]
            Sg_info_Out["Seg"] = float(seg)
            #                dlo.loc[dlo.Seg == seg,:]
            FileInfoOut = pd.concat([Sg_info_Out, FileInfoOut], axis=0, sort=False)
            SgOut = Sg.drop(
                columns=[i[0] for i in Sg_info] + ["j A/cm2", "E_Applied_VRHE"]
            )
            SgOut.to_excel(
                DestDirSegs.joinpath("%s_%.0f" % (N2_BaseFile, int(seg))).with_suffix(
                    ".xlsx"
                )
            )
        action_out.to_excel(
            DestDirSegs.joinpath("%s_ActionsOut" % (N2_BaseFile)).with_suffix(".xlsx")
        )
    #            FileActSegs = pd.merge(dlo,FileInfoOut,on=['Seg'],how='outer')
    #            FileActSegs.to_csv(DestDirSegs.joinpath('%s_ActionsOut' %(N2_BaseFile)).with_suffix('.csv'))
    #                gr.plot(x='Z Real',y='Z Imag',kind='scatter')
    # %%
    #            dlo.to_csv(os.path.join(CVrow['dest_dir'],os.path.splitext(os.path.basename(N2_file))[0]+'_ACTIONS.csv'))
    #    DATA_list.append([segment1])
    #    Action_List.append(action_out)
    return (segment1, action_out)


def PF_fit_starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(PF_fit_apply_args_and_kwargs, args_for_starmap)


def PF_fit_apply_args_and_kwargs(fn, args, kwargs):
    return fn(args, **kwargs)


def create_CVs(ovv, **kwargs):
    # %%
    #    N2grF
    #    ovv, PathDB = HPRR_ovv_file,ECRunOVV.PathDB
    #    ovv, PathDB = N2_ovv_file, ECRunOVV.PathDB
    #     ovv, PathDB = ORR_gr_ovv,ECRunOVV.PathDB
    #    n,CVrow = 1,OER_ovv_file
    #    ovv = ORR_ovv_file
    # ovv = Ring_ovv
    #    EvR1HE = 'E_AppV_RHE'
    #    actions = []

    if not ovv.PAR_file.nunique() == len(ovv):
        _dups = ovv.PAR_file.value_counts().loc[(ovv.PAR_file.value_counts() > 1)].index
        logger.warning('Create CVs received duplicate PAR_files: {", ".join(_dups)}')

    multi_run = kwargs.get("multi_run", False)
    # === OVV FILTERS ====
    #    for N2_file in ovv.query('PAR_exp == "N2" & filesize < 13123896')['PAR_file']:
    #    .query('(PAR_exp == "N2_act"|PAR_exp == "HPRR"|PAR_exp == "OER"|Gas == "O2") & filesize < 17123896').iterrows():
    DATA_list, Action_List = [], []
    multi_failed = False
    if multi_run == True:
        run_args = [gr for n, gr in ovv.groupby("PAR_file")]
        run_kwargs_iter = repeat(kwargs)
        try:
            logger.info(f"{__name__} START multiprocessing for len{len(run_args)}")
            pool_size = os.cpu_count() - 2
            #            if 'linux' in sys.platform:
            #                os.system("taskset -p 0xff %d" % os.getpid())
            with multiprocessing.Pool(pool_size) as pool:
                _read_CV_data = PF_fit_starmap_with_kwargs(
                    pool, read_data_CV_row, run_args, run_kwargs_iter
                )
        except Exception as e2:
            #        print('EIS eis_run_group_ovv multiprocessing error: {0}'.format(e2))
            logger.error(
                f"{__name__} START multiprocessing for len{len(run_args), {e2}}"
            )
            #            logger.error(f'EIS eis_fit_PAR_file  multiprocessing erroe: {e2}, len out({len(fit_export_EV_all)})')
            multi_failed == True

    if multi_run == False or multi_failed == True:
        _read_CV_data = []
        for pf, ovv_row in ovv.groupby("PAR_file"):
            pf, ovv_row
            #        if CVrow.empty:
            #            return pd.DataFrame()
            _pf_read = read_data_CV_row(ovv_row)
            _read_CV_data.append(_pf_read)

    DATA_list = [i[0] for i in _read_CV_data]
    Action_List = [i[1] for i in _read_CV_data]

    try:
        ActionOut = pd.concat([i for i in Action_List], axis=0, sort=False)
    except Exception as e:
        logger.warning(
            "Create CV Create CVs empty:%s, %s \n %s"
            % (len(Action_List), e, [Path(i).name for i in ovv.PAR_file.values])
        )
        ActionOut = pd.DataFrame()
    try:
        DATA_Out = pd.concat([i for i in DATA_list], axis=0, sort=False)
    except Exception as e:
        logger.warning(
            "Create CV Create CVs empty:%s, %s \n %s"
            % (len(DATA_list), e, [Path(i).name for i in ovv.PAR_file.values])
        )
        DATA_Out = pd.DataFrame()
    return DATA_Out, ActionOut


# pd.concat(DATA_list,sort=False)


def CreateCV_multi(ovv, OutPutRaw=False):

    try:
        logger.info(
            f"EIS eis_run_group_ovv START multiprocessing {multi_par_fit} for len{len(fit_run_args)}"
        )
        pool_size = os.cpu_count() - 2
        #            if 'linux' in sys.platform:
        #                os.system("taskset -p 0xff %d" % os.getpid())
        with multiprocessing.Pool(pool_size) as pool:
            fit_export_EV_all_chunck = PF_fit_starmap_with_kwargs(
                pool, fit_EEC, fit_run_args, fit_kwargs_iter
            )
            fit_export_EV_all.append(fit_export_EV_all_chunck)
    except Exception as e2:
        #        print('EIS eis_run_group_ovv multiprocessing error: {0}'.format(e2))
        logger.error(
            f"EIS eis_fit_PAR_file  multiprocessing erroe: {e2}, len out({len(fit_export_EV_all)})"
        )
