"""
Created on Fri Nov 12 12:35:19 2021

@author: DW
"""

###### === Calculations of Selectivity H2O2 yield, n electron from I(Ring and I(Disk) ==== #######
def RRDE_calc_H2O2(_DF, CollEff, SA_ring, mA):
    """Calculations of Selectivity H2O2 yield, n electron from I(Ring and I(Disk)"""
    _I_disk, _I_ring = np.abs(_DF["I(A)_disk"]).values, np.abs(_DF["I(A)_ring"].values)
    _FracH2O2 = 200 * (_I_ring / CollEff) / (_I_disk + _I_ring / CollEff)
    _J_ring = mA * _I_ring / (CollEff * SA_ring)
    n_ORR = 4 * _I_disk / (_I_disk + _I_ring / CollEff)
    _RRDE = {"Frac_H2O2": _FracH2O2, "J_ring": _J_ring, "n_ORR": n_ORR}
    _DF = _DF.assign(**_RRDE)
    return _DF
