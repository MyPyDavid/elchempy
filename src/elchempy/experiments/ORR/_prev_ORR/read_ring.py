"""
Created on Fri Nov 12 12:31:03 2021

@author: DW
"""


def ORR_read_Ring_file(gr_ovv, ORR_ovv_file):
    #                       PAR_file, fstem, ORR_file_PAR_date):
    #                pd.concat([E_onset,Diff_lim,E_half])
    ### +++++++++- Choosing Ring File and Adding I_ring from Chrono -+++++++++++ ###
    #            ('Disc','Ring'),('disc','ring'), ('Ch1_Disk','Ch2_Ring'), ('disk','ring'), ('Disc','Ring')
    Ring_ovv = pd.DataFrame()
    PAR_file, fstem = (
        ORR_ovv_file.PAR_file.values[0],
        ORR_ovv_file.basename.values[0],
    )

    ring_date_ovv = gr_ovv.loc[
        (
            gr_ovv.PAR_date.isin(ORR_ovv_file.PAR_date)
            | gr_ovv.PAR_date_min.isin(ORR_ovv_file.PAR_date_min)
        )
        & gr_ovv.EXP_dir.isin(ORR_ovv_file.EXP_dir)
        & (gr_ovv.Electrode == "Pt_ring")
    ]
    #    ring_date_ovv = gr_ovv.query('(PAR_date == @ORR_file_PAR_date) & (Electrode == "Pt_ring")')
    #    if ring_date_ovv.empty:
    #        ORR_f_par_date_min = f'{pd.to_datetime(ORR_file_PAR_date):%Y-%m-%dT%H:%M}'
    #        ring_date_ovv = gr_ovv.query('(PAR_date_min == @ORR_f_par_date_min) & (Electrode == "Pt_ring")')

    if ring_date_ovv.empty:
        ORR_file_PAR_date = ORR_ovv_file.PAR_date.values[0]
        gr_ovv["PAR_exp_delta_disk"] = [
            pd.to_datetime(i) - pd.to_datetime(ORR_file_PAR_date)
            for i in gr_ovv.PAR_date.values
        ]
        ring_date_ovv = gr_ovv.loc[
            (gr_ovv["PAR_exp_delta_disk"] > pd.to_timedelta(0))
            & (gr_ovv["PAR_exp_delta_disk"] < pd.to_timedelta(10, unit="s"))
        ]
    #                    ring_date_ovv = ovv.loc[ovv.Comment.str.contains('Pt|chrono') & np.isclose(ovv.PAR_date,ORR_file_PAR_date)]
    if len(ring_date_ovv) == 1:
        Ring_ovv = ring_date_ovv
        _logger.info(
            "Jkin calculation O2_chrono_Ring match found: {0} for {1}".format(
                ring_date_ovv.basename.iloc[0], fstem
            )
        )
    else:
        _logger.warning(
            "Jkin calculation O2_chrono_Ring matches error: {0} for {1}".format(
                ring_date_ovv.basename.values, fstem
            )
        )
        O2_PAR_file_upper = fstem.upper()
        O2_ring_fn_options_disC = [
            ("DISC", "RING"),
            ("CH1_DISC", "CH2_RING"),
            ("V3F_DISC", "V3_RING"),
            ("_V3.PAR", "_V3F.PAR"),
            ("_V3.PAR", "_V3F.PAR"),
        ]

        O2_ring_fn_options_disK = [
            (i[0].replace("DISC", "DISK"), i[1])
            for i in O2_ring_fn_options_disC
            if "DISC" in i[0]
        ]
        O2_ring_fn_options = O2_ring_fn_options_disC + O2_ring_fn_options_disK

        O2_file_ring_opts = []
        for _disc_opt, _ring_opt in O2_ring_fn_options:
            _ring_fn_test = [
                i
                for i in gr_ovv.PAR_file.unique()
                if Path(i).name.upper()
                == O2_PAR_file_upper.replace(_disc_opt, _ring_opt)
            ]
            if _ring_fn_test:
                O2_file_ring_opts.append(_ring_fn_test)

        O2_file_ring = list(
            set(
                [
                    i
                    for i in O2_file_ring_opts
                    if O2_PAR_file_upper not in str(i).upper()
                ]
            )
        )
        if O2_file_ring:
            #                    O2_chrono_CV = O2_CVs.loc[(O2_CVs.File.isin(O2_file_ring) | (O2_CVs.File == PAR_file))].query('SampleID == @O2_act_SampleID')
            Ring_ovv = gr_ovv.query("PAR_file == @O2_file_ring")
        else:
            _logger.warning(
                "Jkin calculation O2_chrono_Ring file options empty: {0}".format(
                    ring_date_ovv.basename.values
                )
            )
            Ring_ovv = pd.DataFrame()
    #        O2_fr_1 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('DISC', 'RING')]
    #        O2_fr_2 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('DISK', 'RING')]
    #        O2_fr_3 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('CH1_DISK', 'CH2_RING')]
    #        O2_fr_4 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('CH1_DISC', 'CH2_RING')]
    #        O2_fr_5 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('V3F_DISK', 'V3_RING')]
    #        O2_fr_6 = [i for i in gr_ovv.PAR_file.unique() if
    #                   Path(i).name.upper() == O2_PAR_file_upper.replace('_V3.PAR', '_V3F.PAR')]
    #        O2_file_ring_opts = O2_fr_1 + O2_fr_2 + O2_fr_3 + O2_fr_4 + O2_fr_5 + O2_fr_6
    #        O2_file_ring = list(set([i for i in O2_file_ring_opts if O2_PAR_file_upper not in str(i).upper()]))[0]
    #                print('Ring files {0} for PAR file: {1}'.format(O2_file_ring,PAR_file))
    if not Ring_ovv.empty:
        if Ring_ovv.basename.nunique() != 1:
            _logger.warning(
                "Jkin calculation create CV from O2_Chrono Ring multiple Ring Files in ovv!!: {0}".format(
                    Ring_ovv.basename.values
                )
            )
        _logger.info(
            "Jkin calculation used; DISK: {0},\n RING: {1}".format(
                Path(PAR_file).name, Ring_ovv.basename.iloc[0]
            )
        )
        O2_chrono_Ring, O2_chrono_actions = create_CVs(Ring_ovv)
        if O2_chrono_Ring.empty:
            O2_chrono_Ring = pd.DataFrame()
            _logger.warning(
                "Jkin calculation create CV from O2_Chrono Ring failed: {0}".format(
                    Ring_ovv.basename.iloc[0]
                )
            )
    else:
        _logger.warning(
            "Jkin calculation O2_chrono_Ring MATCH OVV empty: {0}".format(
                ring_date_ovv.basename.values
            )
        )
    return O2_chrono_Ring, O2_chrono_actions


def linear(x, a, b, *args):
    return a * x + b


def _get_merge_cols(ORR_disk_pars, _ORR_RRDE_pars):
    _l, _r = ORR_disk_pars, _ORR_RRDE_pars
    _l1 = [(i, _l[i].unique()[0]) for i in _l.columns if _l[i].nunique() <= 1]
    _r1 = [(i, _r[i].unique()[0]) for i in _r.columns if _r[i].nunique() <= 1]
    _mcls = set([i[0] for i in _l1]).intersection([i[0] for i in _r1])
    pd.merge(_l, _r, on=list(_mcls) + ["Sweep_Type"])
