

import datetime
from pathlib import Path
# from html.parser import HTMLParser

from typing import Tuple

# 3rd party
import datefinder

import pandas as pd

#%%


def get_starttime_from_parser(parser, source = False) -> Tuple[datetime.datetime, str]:
    ''' takes the parser and returns a datetime object, start time of experiment read out from metadata'''

    start_time = None
    source = ''

    if 'VersaStudioParser' in parser.__class__.__qualname__:
        # isinstance(parser, VersaStudioParser):
        # cast metadata from parser into DataFrame
        time = parser.metadata.get('experiment' , {}).get('TimeAcquired','')
        date = parser.metadata.get('experiment' , {}).get('DateAcquired','')
        date_time = time + ' ' + date
        dates = list(datefinder.find_dates(date_time, source=True))
        if len(dates) == 1:
            start_time,  source = dates[0]

    return start_time, source

def cast_parser_to_dataframe(parser):
    ''' takes the parser and returns DataFrames from the metadata, actions and data of the files'''
    metadata, actions, data = None, None, None
    pm_dict, pa_dict, data_dict = {},{},{}

    if 'VersaStudioParser' in parser.__class__.__qualname__:
        # isinstance(parser, VersaStudioParser) or

        # cast metadata from parser into DataFrame
        pm_dict = parser.metadata.copy()
        metadata = pd.DataFrame(pm_dict).T

        # cast actions from parser into DataFrame
        pa_dict = parser.actions.copy()
        actions = pd.DataFrame(pa_dict).T

        p_db_keys = parser.data_body.keys()
        data_body_key_default = ["segment1"]

        if set(p_db_keys) > set(data_body_key_default):
            _extra_keys = set(p_db_keys) - set(data_body_key_default)
            logger.warning(
                f'Unexpected extra keys found in parser,{", ".join(map(str,_extra_keys))}'
            )
        # cast data from parser into DataFrame
        _data_lst = []
        data_dict = parser.data_body.copy()
        for dbkey in p_db_keys:
            # if there are multiple data segments found in the parser
            # is most likely not the case and contains segment1 only

            if parser.data_body[dbkey]:
                pdb_dict = parser.data_body[dbkey].copy()
                # data_dict = {**data_dict, **pdb_dict}
                df = pd.DataFrame(
                    data=pdb_dict, columns=parser.data_keys
                )
                if len(p_db_keys) > 1:
                    df = df.assign(**{"parser_data_body_key": dbkey})
                _data_lst.append(df)
        if _data_lst:
            data = pd.concat(_data_lst)
        else:
            data = pd.DataFrame()
    _dicts =  {'dicts' : {'metadata': pm_dict, 'actions': pa_dict, 'data': data_dict}}
    _frames =  {'frames' : {'metadata': metadata, 'actions': actions, 'data': data}}
    parser_dict = { **_dicts, **_frames}
    # _frames = (metadata, actions, data)
    return parser_dict


def cast_metadata_to_flat_dict(metadata_dict, prefix='meta_'):
    # dct = ecd.DR.parser_dict['dicts']['metadata']
    _meta = {f'{prefix}{k0}_{k1}': v1 for k0,v0 in metadata_dict.items() if v0 for k1,v1 in v0.items()}# else for k1,v1 in {k1: 'None'}.items()}
    _meta_none = {f'{prefix}{k0}': None for k0,v0 in metadata_dict.items() if not v0}
    flat_dict = {**_meta, **_meta_none}
    flat_dict = {k.replace(':','_'): val for k,val in flat_dict.items()}

    return flat_dict

