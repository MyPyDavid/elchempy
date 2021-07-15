'''
empty placeholder
'''

from pathlib import Path

import pandas as pd
import logging
logger = logging.getLogger(__name__)

from parsers import read_PAR_file

def _drop_cols_from_data_segment():
    ''' Columns keys to possibly drop from data DataFrame'''
    # Other optional columns: ['E Imag', 'ADC Sync Input(V)','I Imag', 'I Real']
    _drop_cols = ["E2 Imag", "E2 Real", "E2 Status", "E2(V)", "Z2 Imag", "Z2 Real"]
    return _drop_cols

# read_par_file(filepath)

class DataReader:
    '''

    Wrapper for parsers in a class
    and handles the

    Class could be extended for reading out data
    from multiple source files.

    '''

    supported_filetypes = [".par"]

    def __init__(self, filepath: Path, max_bytesize=1*10**10):
         if not isinstance(filepath, Path):
            if isinstance(filepath, str):
                filepath = Path(filepath)
            else:
                raise TypeError("Argument given is not Path nor str")

         self.filepath = filepath

         self.parser = None
         self.actions = pd.DataFrame()
         self.data = pd.DataFrame()

         if filepath.exists():
            filesize = filepath.stat().st_size
            if filesize < max_bytesize and filepath.suffix in self.supported_filetypes:

                parser, actions, data = self.read_file(self.filepath)

                self.parser= parser
                self.actions = actions
                self.data = data

                # self.double_check_data(self.data, expected_values=self.expected_keys_values)
            else:
                logger.warn(f"File too large {filesize} > {max_bytesize}")

         else:
            logger.warn("File does not exist:\n{filepath}")

    def read_file(self, filepath, data_body_key_default = ['segment1']):
        '''

        Parameters
        ----------
        filepath : Path
            chooses the parsing method depending on the file suffix

        Returns
        -------
        parser_instance : object
            instance of parser
        actions : pd.DataFrame
            actions table
        data : pd.DataFrame
            data table
        '''

        actions = pd.DataFrame()
        data = pd.DataFrame()
        parser = None

        suffix = filepath.suffix
        if suffix in self.supported_filetypes:

            if '.par' in suffix:
                parser = read_PAR_file(filepath)
                # this parser has dictionary attributes which need to cast in to DataFrames

                actions = pd.DataFrame(parser.actions).T

                p_db_keys = parser.data_body.keys()

                if set(p_db_keys) > set(data_body_key_default):
                    _extra_keys = set(p_db_keys) -set(data_body_key_default)
                    logger.warn('Unexpected extra keys found in parser,{", ".join(map(str,_extra_keys))}')

                _data = []
                for dbkey in p_db_keys:
                    # if there are multiple data segments found in the parser
                    # is most likely not the case and contains segment1 only
                    df = pd.DataFrame(data=parser.data_body[dbkey], columns=parser.data_keys)
                    if len(p_db_keys) > 1:
                        df =  df.assign(**{'parser_data_boy_key' : dbkey})
                    _data.append(df)
                data = pd.concat(_data)

            elif '.other' in suffx:
                pass

        else:
            logger.warn('Filetype is not supported, not in {", ".join(map(str,self.supported_filetypes))} ')

        return parser, actions, data

    def __repr__(self):
        _name = self.filepath.name
        _txt = f'actions = {len(self.actions)}, data = {len(self.data)}'
        return f'{self.__class__.__qualname__} on {_name}, {_txt}'

