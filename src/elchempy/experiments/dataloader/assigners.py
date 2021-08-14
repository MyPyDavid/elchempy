"""

These functions assign new columns to an existing DataFrame that contains e.g. the raw data

Created on Sat Aug 14 13:27:30 2021

@author: DW


"""

from typing import Tuple, Dict

import pandas as pd
import numpy as np

import logging
logger = logging.getLogger(__name__)

class AssignerError(Exception):
    pass

class AssignECDataColumns:

    __slots__ = ['data']

    def __init__(self, data):
        if data:
            if isinstance(data, pd.DataFrame):
                self.data = data
                self.assign_all()
            else:
                raise TypeError('argument data is not of type DataFrame')


def assign_all(data : pd.DataFrame,
               actions=pd.DataFrame(),
               raw_current_key: str='I(A)',
               **kwargs) -> pd.DataFrame:
    '''
    Checks first if input data is a DataFrame, assigns columns, then checks if actions DataFrame is also there

    Parameters
    ----------
    data : pd.DataFrame
        DESCRIPTION.
    actions : TYPE, optional
        DESCRIPTION. The default is pd.DataFrame().
    **kwargs : TYPE
        DESCRIPTION.

    Raises
    ------
    AssignerError
        DESCRIPTION.
    TypeError
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    '''
    NEW_COLUMNS = ['E_AppV_RHE', 'j A/cm2', 'jmAcm-2',
                    'scanrate_calc', 'scanrate_prog_calc',
                    'scanrate',
                    'SweepType'
                    ]

    if not isinstance(data, pd.DataFrame) :
        raise AssignerError('argument data or actions is not of type DataFrame')
    if not isinstance(actions, pd.DataFrame):
        raise AssignerError('argument data or actions is not of type DataFrame')

    if data.empty:
        return data

    _columns_init = set(data.columns.copy(deep=True))

    try:
        # j_mA  = get_current_density(data[raw_current_key].to_numpy(), **kwargs)
        j_mA, j_mA_info = get_current_density(data[raw_current_key].to_numpy(),convert_to_mA=True, **kwargs)
        j_A, j_A_info = get_current_density(data[raw_current_key].to_numpy(),convert_to_mA=False, **kwargs)
        data = data.assign(**{**j_A, **j_mA})

        data = assign_scanrate_calc(data)
        data = assign_sweep_type(data, scanrate_key = 'scanrate_calc')
        data = add_scanrate_per_segment(data)

        data = assign_E_vs_RE(data, **kwargs)
        # data = assign_electrochemical_data_columns(data, **kwargs)

    except Exception as exc:
        logger.warning(f"Could not assign all electrochemical columns to data.\n{exc}")

    if actions.empty:
        return data

    try:
        data = assign_action_number_to_data_from_actions_table(data , actions, **kwargs)
        data = assign_action_type_from_reference_table(data, **kwargs)
    except Exception as exp:
        logger.warning("Could not assign all electrochemical columns to data from actions.\n{exc}")

    _columns_end = set(data.columns)
    diff_cols = _columns_end - _columns_init
    # print(f'Diff columns: {", ".join(map(str, diff_cols))} ')
    return data, diff_cols


        #     else:
        #         raise TypeError('argument data is not of type DataFrame')

        # if not data.empty:
        #         print(f"{self.__class__.__qualname__} error {e} for\n {}")

def get_current_density(raw_current_I: np.array,
                           unit_I_raw: ['A', 'mA'] = 'A',
                           unit_SA: ['cm2', 'm2'] = 'cm2',
                           convert_to_mA: bool = False,
                           geometric_SA: float = 20) -> Tuple[Dict, Dict]:
    '''
    Normalized the current for the surface area of the electrode.
    Optionially, units can set for conversions.

    Parameters
    ----------
    raw_current_I : np.array
        DESCRIPTION.
    unit_I : ['A', 'mA'], optional
        DESCRIPTION. The default is 'A'.
    unit_j : ['Acm-2', 'mAcm-2'], optional
        DESCRIPTION. The default is 'mAcm-2'.
    unit_SA : ['cm2', 'm2'], optional
        DESCRIPTION. The default is 'cm2'.
    geometric_SA_cm2 : float, optional
        DESCRIPTION. The default is 20.

    Returns
    -------
    current_density_columns : dict
        {"j_unit" : np.array}

    '''


    if unit_SA == 'cm2' and geometric_SA > 100:
        logger.warning('Is this the Electrode Surface Area in {unit_SA}? {geometric_SA_cm2}.\nThis value seems too large..')

    # unit_j: ['Acm-2', 'mAcm-2'] = 'mAcm-2',
    if convert_to_mA and 'm' in unit_I_raw:
        conversion_mA, add_M  = 1, False
    elif convert_to_mA and not 'm' in unit_I_raw:
        conversion_mA, add_M = 1E3, True
    elif not convert_to_mA and 'm' in unit_I_raw:
        conversion_mA, add_M = 1E-3, False
    elif not convert_to_mA and not 'm' in unit_I_raw:
        conversion_mA, add_M = 1, False

    if add_M:
        current_density_key = f'j_m{unit_I_raw}_{unit_SA}'
    else:
        current_density_key = f'j_{unit_I_raw}_{unit_SA}'

    current_density_column = {
        current_density_key : conversion_mA * raw_current_I / geometric_SA
        }

    info = {'unit_I_raw' :unit_I_raw,
            # 'unit_j' :unit_j,
            'unit_SA' : unit_SA,
            'geometric_SA' : geometric_SA,
            'current_density_key' : current_density_key,
            'conversion_mA' : conversion_mA,
            'convert_to_mA' : convert_to_mA}
    return current_density_column, info



def assign_E_vs_RE(data: pd.DataFrame,
                   RE_potential_V: float = 0,
                   RE_name: ['RHE','AgAgCl', str]= 'RHE'):
    '''
    EappV_RHE: the applied potential versus the RHE potential.
        requires: RHE OCP potential value in Volt or mV
        notes: RHE_OCP potential value is read/guesses from the filename or taken from
        another experimental file or left at 0

    Parameters
    ----------
    data : pd.DataFrame
        DESCRIPTION.
    reference_potential_V : float, optional
        DESCRIPTION. The default is 0.
    reference_electrode : str, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    data : pd.DataFrame
    '''

    E_RHE_columns = {
        f'E_vs_{RE_name}': data['E(V)'] + RE_potential_V,
        f'E_programmed_modulation_vs_{RE_name}': data['E Applied(V)'] + RE_potential_V
        }

    data = data.assign(**E_RHE_columns)

    return data


def unit_ratio(A, B):
    pass



def assign_scanrate_calc(data: pd.DataFrame ):
    """
    Takes a row from the EC_index and reads in the the data from the PAR_file segments
    returns a namedtuple: r.data and r.actions

    Parameters
    ----------
    action : TYPE

    Returns
    -------
    data: same object as input, with assigned columns

    Comments
    -------
    Assigns extra colums that are important to the electrochemical data processing.


    j A/cm2: current density, current divided by geometric electrode surface area
        requires: geometric surface area of electrode in cm2

    scanrate_calc: the applied scanrate, ratio of dE/dt
        requires: absolute diff(E_column) and
        notes: defined as positive number so absolute of column is taken

    E_programmed_modulation: the programmed potential modulation
        notes: do not use this column for final data
    E_progr_scanrate_calc: the applied scanrate, ratio of dE_programmed_modulation/dt
        notes: this column is used for Sweep Type assignment


        notes: takes certain string values
    """
    # DR = DataReader(filepath)
    # segment1 = DR.data
    # actions = DR.actions

    # TODO split in separate functions !!!
    # E_RHE_columns = self.assign_E_vs_RE(data, 1)

    # current_density_columns = {
    #     'j A/cm2': data['I(A)'] / geometric_SA,
    #     'jmAcm-2': 1000 * data['I(A)'] / geometric_SA
    #     }

    scanrate_calc_columns = {
        'scanrate_calc': np.round(
            (data['E(V)'].diff() / data['Elapsed Time(s)'].diff()), 3
        ),
        'scanrate_prog_calc': np.round(
            (data['E Applied(V)'].diff() / data['Elapsed Time(s)'].diff()), 3
        ),
    }

    # TODO split in separate functions !!!
    # EC_important_columns = {
                            # **current_density_columns,
                            # **scanrate_calc_columns}
    data = data.assign(**scanrate_calc_columns)
    # data.scanrate_calc = data.scanrate_calc.fillna(method="backfill")
    # data.scanrate_prog_calc= data.scanrate_prog_calc.fillna(method="backfill")
    _fillna_cols = ['scanrate_calc', 'scanrate_prog_calc']
    data[_fillna_cols] = data[_fillna_cols].fillna(method='backfill')
    return data


def add_scanrate_per_segment(data,
                             segment_key: str = 'Segment #',
                             sr_round: int = 3,
                             scanrate_key: str = 'scanrate_calc'):
    # add simple scanrate per segment
    # segkey = 'Segment #'
    # sr_round = 3
    try:
        # check and add values per segment
        segment_sr_mapping_lst = []
        for n,gr in data.groupby(segment_key):
            sr_calc_mean = 0
            if len(gr) > 3:
                # sr_prog_mean = np.round(np.abs(gr.iloc[1:-1].scanrate_prog_calc).mean(), sr_round)
                sr_calc_mean = np.round(np.abs(gr.iloc[1:-1][scanrate_key]).mean(), sr_round)
                # sr_value = sr_prog_mean if np.isclose(sr_prog_mean,sr_calc_mean,atol=0.02) else sr_calc_mean
                # gr.ActionId.unique()
                segment_sr_mapping_lst.append((n, sr_calc_mean))
        data['scanrate'] = data[segment_key].map(dict(segment_sr_mapping_lst))
    except Exception as exc:
        logger.warning(f'error in adding simple sr values per segment.\n{exc}')
    return data

def assign_sweep_type(data, scanrate_key = 'scanrate_calc'):
    '''
    SweepType: the direction of the scan, ['anodic', 'cathodic', 'chrono', None]
        requires: scanrate with positive and negative numbers
    '''

    sr_arr = data[scanrate_key]

    sweeptype = {'SweepType': np.where(
        sr_arr > 0,
        'anodic',
        np.where(
            sr_arr < 0,
            "cathodic",
            np.where(sr_arr == 0, 'chrono', None),
        ),
    )}

    data = data.assign(**sweeptype )
    data.SweepType = data.SweepType.fillna(method="backfill")
    return data

def assign_action_type_from_reference_table(data : pd.DataFrame,
                                            type_action_key = 'action_type',
                                            actionID_key: str = 'ActionID',
                                            ActionId_to_Type_Reference: dict = {}):
    '''
    Goal:
        assign a type of action to each segment number
        in the data table,
        The action table contains action types, however,
        the data table does not.
        Assign the RPM values derived from actions to
        the corresponding data segments.

    How:
        - merge the actions with segments in the data depending on the number of segments
          found in both the actions and data tables

        - matches the actionID column with a reference dictionary
          containing the action types

        - assign type depending on shape of data in segment


    Parameters
    ----------

    data : pd.DataFrame
        contains the data table

    Returns
    -------
    data : pd.DataFrame
        contains the data table + new "action_type" column
    '''

    # SRunq = segment1["ScanRate_calc"].round(3).unique()
    data[type_action_key] = data[actionID_key].map({val : key for key,val in ActionId_to_Type_Reference.items()})
    return data

def assign_action_number_to_data_from_actions_table(
                                                    data: pd.DataFrame,
                                                    actions: pd.DataFrame,
                                                    actionID_key: str= 'ActionID'
                                                                            ):
    '''
    Parameters
    ----------
    actions : pd.DataFrame
        contains the actions table

    data : pd.DataFrame
        contains the data table

    Returns
    -------
    data : pd.DataFrame
        contains the data table + new "action_numer" column
    '''

    data_segs_sum = len(data['Segment #'].unique())
    action_segs_sum = actions.Segments.sum()

    matching_segments = data_segs_sum == action_segs_sum
    data_seg_grp = data.groupby('Segment #')

    _segcounter = 0
    action_seglist = []
    for actname, actrow in actions.iterrows():
        if actrow.Segments > 0 and _segcounter <= data_segs_sum:
            seglist = list(range(_segcounter, _segcounter + actrow.Segments))
            np_seglist = np.linspace(_segcounter, _segcounter + actrow.Segments, actrow.Segments)
            _segcounter += actrow.Segments
            action_seglist.append((actname, actrow.Name , seglist,_segcounter))
                                   # , [len(data_seg_grp.get_group(segn)) for segn in seglist]))
    if action_seglist:
        action_seg_dict = dict([(i,a[0]) for a in action_seglist for i in a[2]])
        data[actionID_key] = data['Segment #'].map(action_seg_dict)
    return data
