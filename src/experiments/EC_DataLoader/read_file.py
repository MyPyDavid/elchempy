import pandas as pd
import hashlib
import numpy as np
from pathlib import Path


def take_hash(file):
    BLOCKSIZE = 65536
    hasher = hashlib.md5()
    with open(file, "rb") as afile:
        buf = afile.read(BLOCKSIZE)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(BLOCKSIZE)
    hash_out = hasher.hexdigest()
    return hash_out


def get_DF_soup_part(Soup_segment, hash_f):
    #    Soup_segment,hash_f = N2_soup_par.find_all(name=actie)[0],actie
    #    = N2_soup_par.instrument.string
    try:
        splt = [i.split("=") for i in Soup_segment.text.split("\n")[1:-1]]
        out = pd.DataFrame(dict([(i[0], i[1]) for i in splt]), index=[hash_f])
    except:
        splt = [i.split("=") for i in Soup_segment.text.split("\n")[1:-1]][0:20]
        out = pd.DataFrame(dict([(i[0], i[1]) for i in splt]), index=[hash_f])
    return out


def common(action_DATA):
    try:
        if (
            action_DATA["Comment"]
            .str.contains("Enter comments about the experiment here")
            .loc[0]
            == True
        ):
            print("Standard comment")
        if (
            action_DATA["Comment"]
            .str.contains("Enter comments about the experiment here")
            .loc[0]
            == False
        ):
            #       action_DATA['Comment']
            print("Special comment")
    except Exception as e:
        print("ERROR in common, %s" % e)
    return action_DATA.loc[:, ["Measured Open Circuit", "Comment", "Segments"]]


def Cyclic_Voltammetry(action_DATA):
    return action_DATA.loc[
        :,
        [
            "Vertex 1 Potential (V),Value",
            "Vertex 2 Potential (V),Value",
            "Scan Rate (V/s)",
            "Total Points",
            "Points Per Cycle",
            "Cycles",
        ],
    ]


def Open_Circuit(action_DATA):
    return action_DATA.loc[:, ["Duration (s)", "Time Per Point (s)"]]


def Chronoamperometry(action_DATA):
    return action_DATA.loc[:, ["Potential (V),Value", "Duration (s)"]]


def collect_action_data(action_name, actionsegment, instrument_df, date_df, hash_f):
    #    action_name,actionsegment = action.name,action
    #    a_read = pd.read_table((io.StringIO(str(actionsegment.get_text))))
    action_DATA = get_DF_soup_part(actionsegment, hash_f)
    Technique = action_DATA["Name"]
    action_out = pd.DataFrame(
        {
            "Action name": action_name,
            "Technique": Technique,
            "action#": int(action_name.strip("action")),
        },
        index=[0],
    )
    action_out = pd.concat(
        [action_out, instrument_df["Type_action"], date_df["DTC"]], axis=1, sort=False
    )
    if Technique.str.contains("Cyclic Voltammetry").loc[0] == True:
        action_out = pd.concat(
            [action_out, Cyclic_Voltammetry(action_DATA)], axis=1, sort=False
        )

    if Technique.str.contains("Common").loc[0] == True:
        action_out = pd.concat([action_out, common(action_DATA)], axis=1, sort=False)
    if Technique.str.contains("Open Circuit").loc[0] == True:
        action_out = pd.concat(
            [action_out, Open_Circuit(action_DATA)], axis=1, sort=False
        )
    #        print(Open_Circuit(action_DATA))
    if Technique.str.contains("Chronoamperometry").loc[0] == True:
        action_out = pd.concat(
            [action_out, Chronoamperometry(action_DATA)], axis=1, sort=False
        )
    #    print(action_out.T)
    return action_out.T


def collect_data_pars_segment(pars_segment):
    seg_cols = pars_segment.string.split("\n")[3].strip("Definition=").split(", ")[0:-1]
    seg_data = [i.split(",") for i in pars_segment.string.split("\n")[4:-1]]
    return pd.DataFrame(data=seg_data, columns=seg_cols).astype(float)


#        seg_cols = N2_soup_par.segment1.string.split('\n')[3].strip('Definition=').split(', ')[0:-1]
#            seg_data = [i.split(',') for i in N2_soup_par.segment1.string.split('\n')[4:-1]]
#            segment1 = pd.DataFrame(data=seg_data,columns=seg_cols).astype(float)
def add_RRDE(action1_name, action, N2_file):
    if action1_name == "RRDE":
        if "Ring|ring" in Path(N2_file).stem:
            # set Ring settings
            add_action = "_selectedringtechnique"
            add_action = ""
        elif "Ring|ring" not in Path(N2_file).stem:
            add_action = "_selecteddisktechnique"

    #        if os.path.basename(N2_file).split('.par')[0].split('_')[-1] != 'Ring':
    #            add_action = '_selecteddisktechnique'
    if action1_name != "RRDE":
        add_action = ""
    return action + add_action


def DAC_V_to_RPM(DAC_V, *args, **kwargs):
    cell_number = 0
    if "test_PAR_file" in kwargs:
        PAR_file_str = str(kwargs["test_PAR_file"])
        if "cell" in PAR_file_str:
            try:
                cell_number = int(PAR_file_str.split("cell")[-1][0])
            except:
                pass
    use_conv = 1
    rpm_list = [0, 200, 400, 900, 1500]
    conv2, conv1 = 1e-03, 1 / 3635
    v3l, pdx = [0, 0.055, 0.11, 0.2475, 0.245, 0.405, 0.410, 0.675], [
        0.0,
        0.2,
        0.4,
        0.9,
        1.5,
    ]
    match = []
    if cell_number > 0:
        if cell_number == 1 or cell_number == 3:
            use_conv = conv1
        elif cell_number == 2:
            use_conv = conv2

    if DAC_V == 0:
        RPM = 0
    #    elif use_conv != 1:  # converting the RPM number based on the voltage and cell number
    #        RPM = DAC_V / use_conv
    else:  # guessing the RPM number based on the voltage
        for n, r in enumerate(v3l + pdx):
            if np.isclose(r, DAC_V, rtol=0.01):
                match.append(r)
        if len(match) == 1:
            if match[0] in v3l:
                RPM = match[0] / conv1
            elif match[0] in pdx:
                RPM = match[0] / conv2
            else:
                RPM = 0

        elif len(match) == 0:
            RPM = DAC_V / conv1
        #            logger.warning(
        #                'DAC_V_to_RPM no match found in predict list DAC {0} used conv {1} for RPM {2}'.format(DAC_V, conv1,
        #                                                                                                       RPM))
        else:
            RPM = 0
    #            logger.warning('DAC_V_to_RPM two RPMS {0} for DAC: {1} set 0 RPM'.format(match, DAC_V))
    try:
        RPM = int(RPM)
    except Exception as e:
        RPM = f"error_int_{RPM}"
    #        logger.warning('RPM DAC no integer {0}'.format(e))
    return RPM
