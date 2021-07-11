"""
Created on Sun Jul 11 11:00:58 2021

@author: DW
"""


globals()["EvRHE"] = "E_AppV_RHE"


def _func():
    def test_plot_N2_BG_data(self):
        fig, ax = plt.subplots(figsize=(12, 12))
        for pf, gr in self.N2_BG_data.groupby("PAR_file"):
            gr.plot(x=EvRHE, y="jmAcm-2", label=pf.parent.name + "/" + pf.name, ax=ax)


def _testing_class(test):
    tt = ORR_run_loop(test, testing_mode=True)
    self = tt

    # fit_run_arg =
    er = iter(self.ORR_collection)
    fit_run_arg = next(er)
    self = fit_run_arg
    calc = ORR_calculations(fit_run_arg)

    cll = ORR_collection(tt)
    self = cll
    _dict = {}
    for pf, grp in self.gr_ovv_disk.groupby(by="PAR_file"):
        pf, grp
        _orr = ORR_scan_data(pf, grp, self.ovv_all, self.EC_index, self.run_kwargs)
        _dict.update({pf: _orr})
    self = _orr
