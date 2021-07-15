'''
empty placeholder
'''

from pathlib import Path

import pandas as pd
import logging
logger = logging.getLogger(__name__)


from parsers import read_PAR_file

def _dev():
    _n2files = Path.cwd().parent.parent.parent.parent.joinpath('data/raw').rglob('*N2*par')
    return _n2files

def _test_read():
    files = _dev()
    results = []
    for filepath in files:
        DR = DataReader(str(filepath))
        results.append(DR)
        actions = DR.actions
        data = DR.data
    # if True:
        data.plot(x='E(V)',y='I(A)', title=filepath.name)
        if any('EIS' in name for name in actions.Name.unique()):
            data.plot(x='Z Real',y='Z Imag', title=filepath.name)


def _drop_cols_from_data_segment():
    ''' Columns keys to possibly drop from data DataFrame'''
    # Other optional columns: ['E Imag', 'ADC Sync Input(V)','I Imag', 'I Real']
    _drop_cols = ["E2 Imag", "E2 Real", "E2 Status", "E2(V)", "Z2 Imag", "Z2 Real"]
    return _drop_cols

# read_par_file(filepath)

class DataReader:
    '''
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

         self.parser_instance = None
         self.actions = pd.DataFrame()
         self.data = pd.DataFrame()

         if filepath.exists():
            filesize = filepath.stat().st_size
            if filesize < max_bytesize and filepath.suffix in self.supported_filetypes:

                parser_instance, actions, data = self.read_file(self.filepath)

                self.parser_instance= parser_instance
                self.actions = actions
                self.data = data

                # self.double_check_data(self.data, expected_values=self.expected_keys_values)
            else:
                logger.warn(f"File too large {filesize} > {max_bytesize}")

         else:
            logger.warn("File does not exist:\n{filepath}")

    def read_file(self, filepath):
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
        parser_instance = None

        suffix = filepath.suffix
        if suffix in self.supported_filetypes:

            if '.par' in suffix:
                parser_instance = read_PAR_file(filepath)
                actions = pd.DataFrame(parser_instance.actions).T
                data = pd.DataFrame(data=parser_instance.data_body['segment1'], columns=parser_instance.data_keys)
            elif '.other' in suffx:
                pass

        else:
            logger.warn('Filetype is not supported, not in {", ".join(map(str,self.supported_filetypes))} ')

        return parser_instance, actions, data

    def __repr__(self):
        _name = self.filepath.name
        _txt = f'actions = {len(self.actions)}, data = {len(self.data)}'
        return f'{self.__class__.__qualname__} on {_name}, {_txt}'

