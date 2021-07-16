"""
Fetches data from a file and constructs the electrochemical data columns
"""

import pandas as pd
import numpy as np

from .reader import DataReader

from .converters import get_current_density, get_potential_vs_RE, get_RPM_from_DAC_V

class ElchemData:
    '''
    This class contains all functions
    which add several collumns to an
    existing DataFrame with information
    that is important for the analsysis
    of electrochemimcal experiments.

    '''

    def __init__(self, filepath):

        self.filepath = filepath

        self.DR = DataReader(filepath)
        self.raw_actions = self.DR.actions
        self.raw_data = self.DR.data
        self.data = self.raw_data.copy(deep=True)
        self.actions = self.raw_actions.copy(deep=True)

        self.data = self.assign_E_vs_RE(self.data)
        self.data = self.assign_electrochemical_data_columns(self.data)
        self.data = self.assign_action_type_from_reference_table(self.data )
        self.data = self.assign_action_number_to_data_from_actions_table(self.data , self.actions)


    @staticmethod
    def assign_E_vs_RE(data: pd.DataFrame,
                       reference_potential_V: float = 0,
                       reference_electrode: str= ''):
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
            'E_AppV_RHE': data['E(V)'] + reference_potential_V,
            'E_programmed_modulation_RHE': data['E Applied(V)'] + reference_potential_V
            }

        data = data.assign(**E_RHE_columns)

        return data

    def assign_electrochemical_data_columns(self,
                                            data,
                                            RHE_potential = 2,
                                            geometric_SA = 20,
                                            electrode_type = ''
                                            ):
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

        SweepType: the direction of the scan, ['anodic', 'cathodic', 'chrono', None]
            requires: scanrate with positive and negative numbers
            notes: takes certain string values
        """
        # DR = DataReader(filepath)
        # segment1 = DR.data
        # actions = DR.actions
        NEW_COLUMNS = ['E_AppV_RHE', 'j A/cm2', 'jmAcm-2',
                        'scanrate_calc', 'scanrate_prog_calc',
                        'scanrate',
                        'SweepType'
                        ]

        # TODO split in separate functions
        # E_RHE_columns = self.assign_E_vs_RE(data, 1)

        current_density_columns = {
            'j A/cm2': data['I(A)'] / geometric_SA,
            'jmAcm-2': 1000 * data['I(A)'] / geometric_SA
            }
        scanrate_calc_columns = {
            'scanrate_calc': np.round(
                (data['E(V)'].diff() / data['Elapsed Time(s)'].diff()), 3
            ),
            'scanrate_prog_calc': np.round(
                (data['E Applied(V)'].diff() / data['Elapsed Time(s)'].diff()), 3
            ),
        }

        EC_important_columns = {
                                **current_density_columns,
                                **scanrate_calc_columns}
        data = data.assign(**EC_important_columns)
        # data.scanrate_calc = data.scanrate_calc.fillna(method="backfill")
        # data.scanrate_prog_calc= data.scanrate_prog_calc.fillna(method="backfill")
        _fillna_cols = ['scanrate_calc', 'scanrate_prog_calc']
        data[_fillna_cols] = data[_fillna_cols].fillna(method='backfill')

        # add simple scanrate per segment
        segkey = 'Segment #'
        sr_round = 3
        try:
            # check and add values per segment
            segment_sr_mapping_lst = []

            for n,gr in data.groupby(segkey):
                sr_value = 0
                if len(gr) > 3:
                    sr_prog_mean = np.round(np.abs(gr.iloc[1:-1].scanrate_prog_calc).mean(), sr_round)
                    sr_calc_mean = np.round(np.abs(gr.iloc[1:-1].scanrate_calc).mean(), sr_round)
                    sr_value = sr_prog_mean if np.isclose(sr_prog_mean,sr_calc_mean,atol=0.02) else sr_calc_mean
                    # gr.ActionId.unique()
                    segment_sr_mapping_lst.append((n, sr_value))
            data['scanrate'] = data[segkey].map(dict(segment_sr_mapping_lst))
        except Exception as e:
            print('error in adding simple sr values per segment')

        data = data.assign(
            **{
                'SweepType': np.where(
                    data.scanrate_calc > 0,
                    'anodic',
                    np.where(
                        data.scanrate_calc < 0,
                        "cathodic",
                        np.where(data.scanrate_calc == 0, 'chrono', None),
                    ),
                )
            }
        )
        data.SweepType = data.SweepType.fillna(method="backfill")
        return data

    @staticmethod
    def assign_action_type_from_reference_table(data : pd.DataFrame):
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
        ActionId_to_Type_Reference = {
            "Cyclic Voltammetry (Multiple Cycles)": 38,
            "Chronoamperometry": 3,
            "Unknown": 0,
            "Potentiostatic EIS": 21,
            "Cyclic Voltammetry": 14,
        }
        type_action_key = 'action_type'

        # SRunq = segment1["ScanRate_calc"].round(3).unique()
        data[type_action_key] = data['ActionId'].map({val : key for key,val in ActionId_to_Type_Reference.items()})
        return data

    @staticmethod
    def assign_action_number_to_data_from_actions_table(
                                                        data: pd.DataFrame,
                                                        actions: pd.DataFrame
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
            data['action_number'] = data['Segment #'].map(action_seg_dict)
        return data

    def __repr__(self):
        _name = self.filepath.name
        _txt = f'actions = {len(self.actions)}, data = {len(self.data)}'
        return f'{self.__class__.__qualname__}: {_name}, {_txt}'

