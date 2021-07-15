"""
Created on Wed Jul  7 13:31:14 2021

@author: DW
"""


def Post_EIS_pars(**kwargs):
    try:
        EIS_pars = get_EIS_pars(kwargs)
        EIS_all_check_redchi(EIS_pars)
        export_to_Origin(EIS_pars)
    except Exception as e:
        _logger.warning(f"Post EC failed because:: {e}")
