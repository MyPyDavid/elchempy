"""
Created on Fri Nov 12 12:31:56 2021

@author: DW
"""


class ORR_KL_loop:
    """
    This class takes the ORR_calculation class instance and run the loop for K-L calculations over it.
    That is loop over each N2_BG file, each type of N2 current ['raw' or 'corrected']
    and finnaly over the ORR scan segments
    """

    def __init__(self, ORR_calc, run_loop=True):
        self.ORR_calc = ORR_calc
        self.fstem = self.ORR_calc.PAR_file_disk.stem
        if run_loop:
            self.calc_loop()

    def calc_loop(self):

        for N2_pf, N2_BG in self.ORR_calc.N2_BG_data_grp:
            N2_pf, N2_BG
            for N2_cols in self.ORR_calc.N2_testcols:
                N2_cols
                ORR_ops_col = []
                if hasattr(self, "ORR_ops_col"):
                    del self.ORR_ops_col
                for seg, O2_disk_seg in self.ORR_calc.O2_act_seggrp:
                    seg, O2_disk_seg
                    try:
                        ORR_ops = ORR_operations(
                            self.ORR_calc, O2_disk_seg, N2_BG, N2_cols
                        )
                        self.ORR_dest_dir_file = ORR_ops.ORR_dest_dir_file
                    except Exception as e:
                        _logger.error("ORR_KL_loop error :======= {e}")
                    if hasattr(ORR_ops, "ORR_KL_data_rpm"):
                        ORR_ops_col.append(ORR_ops)
                if ORR_ops_col:
                    self.ORR_ops_col = ORR_ops_col
                    self.KL_calc()

    # def collect_rpm_data(self):
    #      ORR_KL_data_rpm = ORR_select_KL_data(O2_CalcOut)
    #      _KL_data_all_rpms.append(ORR_KL_data_rpm)

    def KL_calc(self):
        _dt_now = datetime.now()
        # ORR_dest_dir_file = self.ORR_ops_col[-1].ORR_dest_dir_file
        try:
            # _ops_col_KL = []
            # for i in self.ORR_ops_col:
            # if hasattr(i,'ORR_KL_data_rpm'):
            KL_data = pd.concat([i.ORR_KL_data_rpm for i in self.ORR_ops_col])
            KL_data = KL_data.assign(**{"ORR_KL_Analysis_date_dt": _dt_now})
        except Exception as e:
            _logger.error("ORR_KL_loop :======= {0}".format(e))
        #    KLout = KLout.loc[(KLout['Sweep_Type'].str.contains("cathodic|anodic")), :]
        try:
            #            KLdest = 'KL_'+os.path.basename(O2_act.loc[O2_act['hash'] == h1,'File'].unique()[0]).split('.')[0]
            KL_fit = KL_plots(KL_data, self.ORR_dest_dir_file)
        #        KL_plots(KL_data, ORR_dest_dir_file)
        except Exception as e:
            _logger.error("ORR_KL_loop No KL FITTING plot:======= {0}".format(e))
            KL_fit = ["NA", "N"]

        ORR_pars_all_rpms = pd.concat([i.ORR_pars_out for i in self.ORR_ops_col])
        ORR_pars_all_rpms = ORR_pars_all_rpms.assign(
            **{"ORR_Analysis_date_dt": _dt_now}
        )

        ORR_pars_target_base = self.ORR_dest_dir_file.joinpath(f"ORR_pars_{self.fstem}")
        ORR_pars_target = FileOperations.CompareHashDFexport(
            ORR_pars_all_rpms, ORR_pars_target_base
        )
        _logger.info("Jkin calculation ORR Pars EXCEL DEST:{ORR_pars_target}")
