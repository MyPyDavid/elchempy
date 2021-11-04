"""
Created on Mon Nov  1 12:01:50 2021

@author: DW
"""



def get_file_from_secondary_electrode(disk_filepath):
    
    similar_files_options = set()
    
    files_in_folder = [i for i in list(disk_filepath.parent.rglob(f'*{disk_filepath.suffix}')) if i != disk_filepath]
    
    files_with_same_length = {f for f in files_in_folder if len(f.stem) == len(disk_filepath.stem)}
    similar_files_options.update(files_with_same_length)
    
    return similar_files_options
    
    

def _dev_test():
    # for developing
    from elchempy.experiments._dev_datafiles._dev_fetcher import (
        get_files,
        _dev_test_read,
    )
    
    files = get_files()
    filepath = files[8]
    similar_files_options = set()
    
    files_in_folder = [i for i in list(filepath.parent.rglob(f'*{filepath.suffix}')) if i != filepath]
    
    files_with_same_length = {f for f in files_in_folder if len(f.stem) == len(filepath.stem)}
    similar_files_options.update(files_with_same_length )
      
    files_with_same_first_letter = [i for i in files_with_same_length if i.stem.startswith(filepath.stem[0])]
    
    
    
    for i in f.stem:
        pass
        for f2 in files_with_same_first_letter:
            pass
            
    

def finds_file_with_similar_filename_in_folder(filepath):
    '''
    In case of ORR experiment.
    
    gets file from secondary electrode in same folder as filepath
    with a similar basename plus suffix
    
    '''
    
    Path(filepath).stem.split('_')
    return filepath