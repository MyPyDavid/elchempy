'''
# TODO split functions
keep logic for calculations seperate from file handling/ conditions/
make calculations in simple functions with most basic args, kwargs only

'''

# import sys

from pathlib import Path

# from collections import namedtuple
from datetime import datetime
import numpy as np

from scipy.stats import linregress, zscore
import matplotlib.pyplot as plt

import pandas as pd


# TODO phase out file_py_helper dependencies
from file_py_helper.file_functions import FileOperations


# print("File", __file__, "\nName;", __name__)
if __name__ == "__main__":
    import plotting as n2plotting
else:
    import plotting as n2plotting

n2plotting.N2_plot_Cdl_sweeptype_scatter
n2plotting.N2_plot_Cdl_scans_scanrate
# pass

import logging

logger = logging.getLogger(__name__)

#%% notes from analyses


def _old_N2_analysis(ecdata: ElChemData):
    """
    Performs the steps in the N2 analysis and add the
    results to the ElchemData instance.

    Parameters
    ----------
    ecdata : ElchemData
        contains the raw data.

    Returns
    -------
    ecdata : ElchemData
        contains the results in an added method as attribute "N2".
    """
    # ElChemData
    # FIXME Select only CV types from Data segment
    # Select the data for N2 Cyclic Voltammograms
    N2_CVs = ecdata.data.loc[ecdata.data.ActionId == 38]

    N2_CVs = N2_CVs.dropna(subset=["scanrate"]).loc[N2_CVs.scanrate_calc != 0]

    Cdl_pars, Cdl_data = pd.DataFrame(), pd.DataFrame()
    # scanrates = N2_CVs.scanrate.unique()
    if N2_CVs.scanrate.nunique() > 2:
        # if multiple scanrates are present than calculations can start
        Cdl_pars, Cdl_data = CDL(N2_CVs, EvRHE=EvRHE)

    # Check if possible background (BG) scan is in the data
    BG_present = False
    N2_BG = pd.DataFrame()
    if N2_CVs.scanrate.min() < 0.015:
        # check presence of slow scanrates in data
        BG_present = True
    if BG_present:
        N2_BG = get_N2_background_data()

    N2_ecdata = N2_method(N2_CVs, Cdl_pars, Cdl_data, N2_BG)

    ecdata.add_analysis_method(N2_ecdata)

    return ecdata




#%%










def N2_scans(fit_run_arg, **N2_kwargs):
    N2_ovv_file = fit_run_arg.file_ovv
    try:
        N2_scan, Cdl_fit, Cdl_PARS = N2_analyze_scan(N2_ovv_file, **N2_kwargs)
    except Exception as err:
        logger.error(
            f"N2_scans analyze error : {err}, for {N2_ovv_file.PAR_file.values[0]}"
        )
        N2_scan, Cdl_fit, Cdl_PARS = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    return N2_scan, Cdl_fit, Cdl_PARS


def N2_prepare_file_info(N2_ovv_file, **N2_kwargs):
    pass


def N2_analyze_scan(N2_ovv_file, N2_CVs, N2_actions, N2_FileName, **N2_kwargs):
    """
    Contains the logic for a N2 or nitrogen scans experiment

    Parameters
    ----------
    # TODO change to take input as data and ignore the file
    N2_ovv_file : TYPE
        DESCRIPTION.
    **N2_kwargs : TYPE
        DESCRIPTION.

    Returns
    -------

    N2_scan : TYPE
        DESCRIPTION.
    Cdl_fit : TYPE
        DESCRIPTION.
    Cdl_PARS : TYPE
        DESCRIPTION.

    """
    # %%
    #    N2gr = N2_CVs
    #    gr_N2_ovv = N2_ovv_file
    #    gr_N2_ovv = N2act_ovv.loc[683]
    #    gr_N2_ovv = N2act_ovv.groupby(by='PAR_file').get_group(Path('G:/Cloudstation/Experimental data/Raw_data/VERSASTAT/2018-01-Jan/29.01.2018_DW24_bipot_0.1MKOH_RRDE22775/N2_20cls_300_100_10_DW24_900.par'))
    EvRHE = "E_AppV_RHE"
    #    RRDEloading = 0.09  # 0.09 mg_cat on PINE_disk (0.238 cm2)
    #    GeoArea = WE_SA_collection_eff('PINE')['Disk_cm2']

def _dev_N2_file_info():
    ''' prepares destination file names and folder'''
    SampleID = N2_ovv_file["SampleID"].values[0]
    Loading_name = N2_ovv_file["Loading_name"].values[0]
    Electrolyte = N2_ovv_file["Electrolyte"].values[0]
    postAST = N2_ovv_file["postAST"].values[0]
    N2_scans_folder_version = f"N2_scans_v{FileOperations.version}"
    N2_exp_folder = "_".join([Electrolyte, SampleID, Loading_name, postAST])
    N2_dest_dir = Path(N2_ovv_file.Dest_dir.iloc[0]).joinpath(
        f"{N2_scans_folder_version}/{N2_exp_folder}"
    )
    N2_dest_dir.mkdir(parents=True, exist_ok=True)
    #    FolderOps.FileOperations.make_path(N2_dest_dir)

    #    N2_fn = os.path.splitext(os.path.basename(N2gr['File'].unique()[0]))[0]
    #    gr_N2_ovv.PAR_file.iloc[0] = Path(gr_N2_ovv.PAR_file.iloc[0])
    N2_FileName = Path(N2_ovv_file["PAR_file"].unique()[0])
    N2_fn = N2_FileName.stem
    N2_act_BG = Path(N2_dest_dir.parent).joinpath(N2_fn + "_BG.pkl")

    # Read or load in experimental data file, PAR_file...
    if "N2_preload_CV" in N2_kwargs.keys():
        N2_CVs, N2_actions = N2_kwargs.get["N2_preload_CV"]
    else:
        pass
        # N2_CVs, N2_actions = create_CVs(N2_ovv_file)


    def calculate_Cdl():
        ###### ====== Analyze the Capacity (Cdl) of the N2 scan ========== #######
        Cdl_out = pd.DataFrame({EvRHE: []})
        lst_dct, index_out = [], []
        #        N2_Cdl_ovv = pd.DataFrame(data=ScanRates, index=ScanRates, columns=['ScanRate'])
        #        index_info = {'PAR_file': N2_FileName
        _Cdl_scans_plot = []
        if len(ScanRates) > 1:
            for SR in ScanRates:
                SR_mVs = int(SR * 1000)
                #               TODO export_cols = [EvRHE, 'j A/cm2','Scan Rate (V/s)']
                srA = grA.get_group(
                    ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", SR)
                )
                nuniq_segs = srA["Segment #"].nunique()
                uniq_segs = srA["Segment #"].unique()
                N2sr_fn_out_base = N2_dest_dir.joinpath(
                    "CV_%s_%s.xlsx" % (N2_fn, SR_mVs)
                )
                N2sr_DF_out = srA.loc[srA["Segment #"] == uniq_segs[-1], :]

                def check_if_last_seg_isnormal(nuniq_segs, uniq_segs, srA):
                    if nuniq_segs > 1:
                        _N2sr_DF_out = pd.DataFrame()
                        _idx = -1
                        while _N2sr_DF_out.empty and abs(_idx) < nuniq_segs:
                            _lastsegs = uniq_segs[
                                -(nuniq_segs - int(nuniq_segs * 0.8)) :
                            ]
                            _N2sr_segs = srA.loc[
                                srA["Segment #"].isin(_lastsegs),
                            ]
                            _minmean = (
                                _N2sr_segs.groupby("Segment #")
                                .agg("jmAcm-2")
                                .min()
                                .mean()
                            )

                            _N2sr_test = srA.loc[
                                srA["Segment #"] == uniq_segs[_idx],
                            ]
                            _last_mmn = _N2sr_test["jmAcm-2"].min().mean()
                            _dev_perc = 100 * (_minmean - _last_mmn) / _minmean
                            #                            print(f'{_dev_perc}')
                            if abs(_dev_perc) < 25:
                                _N2sr_DF_out = _N2sr_test
                            else:
                                _idx -= 1
                    else:
                        _N2sr_DF_out = srA.loc[
                            srA["Segment #"] == uniq_segs[0],
                        ]
                    return _N2sr_DF_out

                N2sr_DF_out = check_if_last_seg_isnormal(nuniq_segs, uniq_segs, srA)

                N2sr_DF_out = N2sr_DF_out.assign(
                    **{"PAR_file": N2_FileName, "ScanRate_mVs": SR_mVs}
                )
                N2sr_fn_out = FileOperations.CompareHashDFexport(
                    N2sr_DF_out, N2sr_fn_out_base, silent=True
                )
                N2sr_DF_out = N2sr_DF_out.assign(**{"N2_CV_datafile": N2sr_fn_out})

                #                , 'ScanRate': SR, 'DestFile': N2sr_fn_out, 'Segment': last_seg,
                #                              'Type_output': 'N2_CV'}
                #                index_out.append(index_info)
                #                    N2sr_DF_out.to_excel(N2sr_fn_out)
                if nuniq_segs > 1:
                    #                           srA.loc[srA['Segment #'] == srA['Segment #'].unique()[0],[EvRHE,'j A/cm2']].to_excel(N2_dest_dir.joinpath('%s_%s_first.xlsx' %(N2_fn,int(SR*1000))))
                    N2_sr_DF = srA.loc[srA["Segment #"] == uniq_segs[0], :]
                    N2_sr_fn_base = N2_dest_dir.joinpath(
                        "CV_%s_%s_first.xlsx" % (N2_fn, SR_mVs)
                    )
                    N2_sr_DF = N2_sr_DF.assign(
                        **{"PAR_file": N2_FileName, "ScanRate_mVs": SR_mVs}
                    )
                    N2_sr_fn = FileOperations.CompareHashDFexport(
                        N2_sr_DF, N2_sr_fn_base
                    )
                #                    index_info = {'PAR_file': N2_FileName, 'ScanRate': SR, 'DestFile': N2_sr_fn, 'Segment': first_seg,
                #                                  'Type_output': 'N2_CV'}
                #                    index_out.append(index_info)
                Cdl_scan = N2sr_DF_out
                _Cdl_scans_plot.append(Cdl_scan)
                #                srA.loc[srA['Segment #'] == srA['Segment #'].unique()[-1]]
                # Loop over Sweep Types
                for sweep in Cdl_scan["Sweep_Type"].unique():
                    if sweep == "NA":
                        continue
                    # Loop over potentials E This mean for example for E in [0.1, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.9]:
                    for E in np.linspace(0.1, 1, 50):

                        j_Cdl = np.abs(
                            Cdl_scan.loc[
                                (np.isclose(Cdl_scan[EvRHE], E, atol=0.010))
                                & (Cdl_scan["Sweep_Type"] == sweep)
                            ]["j A/cm2"]
                        ).mean()
                        #                    Cdl_std = np.abs(Cdl_scan.loc[(np.isclose(Cdl_scan[EvRHE],E,atol=0.010)) & (Cdl_scan['Sweep_Type'] == sweep)]['j A/cm2']).std()/SR
                        lst_dct.append(
                            {
                                "ScanRate": SR,
                                f"j_{sweep}": j_Cdl,
                                "E_V": E,
                                "N2_CV_datafile": N2sr_fn_out,
                            }
                        )
            #                    N2_Cdl_ovv.assign(**{'ScanRate': SR, f'j_{sweep}' : j_Cdl})
            Cdl_scans = pd.concat(_Cdl_scans_plot)
        def exportingstuff():
            N2_plot_Cdl_scans_scanrate(
                Cdl_scans, N2_dest_dir.joinpath(f"N2_ScanRates_{N2_fn}.png")
            )

            Cdl_data = pd.DataFrame(lst_dct)
            Cdl_fn_out_base = N2_dest_dir.joinpath("Cdl_data_%s.xlsx" % N2_fn)
            Cdl_data = Cdl_data.assign(
                **{"PAR_file": N2_FileName, "DestFile_N2_BG": N2_act_BG}
            )
            Cdl_fn_out = FileOperations.CompareHashDFexport(Cdl_data, Cdl_fn_out_base)
            #            index_info = {'PAR_file': N2_FileName, 'DestFile': Cdl_fn_out, 'Type_output': 'N2_Cdl_data'}
            #            index_out.append(index_info)
            #                Cdl_pars.to_csv(Cdl_fn_out)
        def cdl_lin():
            N2_CV_datafilenames = ", ".join(
                [Path(i).name for i in Cdl_data["N2_CV_datafile"].unique()]
            )
            Cdl_cath = Cdl_data.loc[:, ["j_cathodic", "E_V", "ScanRate"]].dropna(
                axis=0, how="any"
            )
            Cdl_anod = Cdl_data.loc[:, ["j_anodic", "E_V", "ScanRate"]].dropna(
                axis=0, how="any"
            )



            fitl = []

            _meta_info = {
                "SampleID": SampleID,
                "Filename": N2_FileName,
                "Cdl_datafile": Cdl_fn_out,
                "N2_CV_datafilenames": N2_CV_datafilenames,
            }

            ### Take the linear fits of the ScanRate groups and CHECK or TEST for consistency,
            #            if there was one certain scanrate that was not good ###
            def linregress_residual(gr, xcol="ScanRate", ycol="j_anodic"):
                _x, _y = gr[xcol], gr[ycol]
                _lin = linregress(_x, _y)
                _ymod = _x * _lin.slope + _lin.intercept
                _Res = abs(_y) - abs(_ymod)
                gr = gr.assign(
                    **{"lin_mod": _ymod, "lin_Res": _Res, "zscore_y": zscore(_y)}
                )
                return gr

            Cdl_anod_test = pd.concat(
                [
                    linregress_residual(gr, xcol="ScanRate", ycol="j_anodic")
                    for nm, gr in Cdl_anod.groupby(by="E_V")
                ]
            )
            _anodResmean = Cdl_anod_test.groupby("ScanRate").mean()
            _anodTest = _anodResmean.query(
                "lin_Res < 1E-3 & zscore_y < -0.5 & ScanRate > 0.02"
            )
            Cdl_cath_test = pd.concat(
                [
                    linregress_residual(gr, xcol="ScanRate", ycol="j_cathodic")
                    for nm, gr in Cdl_cath.groupby(by="E_V")
                ]
            )
            _cathResmean = Cdl_cath_test.groupby("ScanRate").mean()
            _cathTest = _cathResmean.query(
                "lin_Res < 1E-3 & zscore_y < -0.5 & ScanRate > 0.02"
            )

            if not _anodTest.empty and not _cathTest.empty:
                _wrong_SRs = set(_anodTest.index).union(_anodTest.index)
            elif not _anodTest.empty and _cathTest.empty:
                _wrong_SRs = set(_anodTest.index)
            elif _anodTest.empty and not _cathTest.empty:
                _wrong_SRs = set(_cathTest.index)
            else:
                _wrong_SRs = set()

            if _wrong_SRs:
                Cdl_anod = Cdl_anod.loc[~Cdl_anod.ScanRate.isin(_wrong_SRs)]
                Cdl_cath = Cdl_cath.loc[~Cdl_cath.ScanRate.isin(_wrong_SRs)]

            #            Cdl_scans.query('ScanRate_mVs == 300').plot(x=EvRHE,y='jmAcm-2',kind='scatter')
            ### Take the linear fits of the ScanRate groups and collect results for final Cdl pars ###
            def lin_fit(nm, _swpgrpE, _meta_info, _yfit="j_anodic"):
                _cdl_fit = linregress(_swpgrpE["ScanRate"], _swpgrpE[_yfit])
                #                gr.plot(x='ScanRate',y='j_anodic',label=nm)
                _linfitout = {
                    **_meta_info,
                    **{
                        "Sweep_Type_N2": _yfit.split("j_")[-1],
                        "ScanRates": ", ".join(
                            [str(i) for i in _swpgrpE["ScanRate"].values]
                        ),
                        "j_Acm2_fit": ", ".join(
                            [str(i) for i in _swpgrpE[_yfit].values]
                        ),
                        EvRHE: nm,
                        "Cdl": _cdl_fit.slope,
                        "Cdl_R": _cdl_fit.rvalue,
                        "Cdl_fit": _cdl_fit,
                        "EASA_m2": _cdl_fit.slope,
                    },
                }
                return _linfitout

            _fit_An = [
                lin_fit(nm, gr, _meta_info, _yfit="j_anodic")
                for nm, gr in Cdl_anod.groupby(by="E_V")
            ]
            _fit_Cath = [
                lin_fit(nm, gr, _meta_info, _yfit="j_cathodic")
                for nm, gr in Cdl_cath.groupby(by="E_V")
            ]
            Cdl_fit = pd.DataFrame(_fit_An + _fit_Cath)
            #
            #            for nm, gr in Cdl_anod.groupby(by='E_V'):
            #                lin_fit(nm, gr, _meta_info, _yfit = 'j_anodic')
            #                #            gr.plot(x='ScanRate',y='Cdl_anodic',label=nm,kind='scatter',ax=ax,ylim=(0,2E-3))
            #                cdl_fit = linregress(gr['ScanRate'], gr['j_anodic'])
            #                gr.plot(x='ScanRate',y='j_anodic',label=nm)
            #                fitl.append({**_meta_info, **{'Sweep_Type_N2': 'anodic',
            #                             'ScanRates': gr['ScanRate'].values, 'j_Acm2_fit': gr['j_anodic'].values, EvRHE: nm,
            #                             'Cdl': cdl_fit[0], 'Cdl_R': cdl_fit[2], 'Cdl_fit': cdl_fit,
            #                             'EASA_m2': cdl_fit[0]}})
            #            #                         cdl_fit[0],cdl_fit[1],cdl_fit[2]])
            #            for nm, gr in Cdl_cath.groupby(by='E_V'):
            #                #            gr.plot(x='ScanRate',y='Cdl_cathodic',label=nm,kind='scatter',ax=ax,color='r',ylim=(0,2E-3))
            #                cdl_fit = linregress(gr['ScanRate'], gr['j_cathodic'])
            #                gr.plot(x='ScanRate',y='j_cathodic',label=nm)
            #                fitl.append({**_meta_info, **{'Sweep_Type_N2': 'cathodic',
            #                             'ScanRates': gr['ScanRate'].values, 'j_Acm2_fit': gr['j_cathodic'].values, EvRHE: nm,
            #                             'Cdl': cdl_fit[0], 'Cdl_R': cdl_fit[2], 'Cdl_fit': cdl_fit,
            #                             'EASA_m2': cdl_fit[0]}})
            #            fitl.append([SampleID,'j_cathodic',nm,cdl_fit[0],cdl_fit[1],cdl_fit[2]])
            #            Cdl_fit = pd.DataFrame(fitl)
            Cdl_cath_slice = Cdl_fit.loc[
                (Cdl_fit["Sweep_Type_N2"].str.contains("cathodic") == True)
                & (Cdl_fit["Cdl_R"] > 0.85),
                :,
            ]
            Cdl_cath_corr = linregress(Cdl_cath_slice[EvRHE], Cdl_cath_slice["Cdl"])
            Cdl_cath_slice = Cdl_cath_slice.assign(
                **{
                    "Cdl_corr": Cdl_cath_slice["Cdl"]
                    - (
                        Cdl_cath_slice[EvRHE].values * Cdl_cath_corr[0]
                        + Cdl_cath_corr[1]
                    )
                    + Cdl_cath_slice["Cdl"].mean()
                }
            )
            #                               columns=['SampleID','Sweep',EvRHE,'Cdl','Cdl_std','Cdl_R'])
            Cdl_an_slice = Cdl_fit.loc[
                (Cdl_fit["Sweep_Type_N2"].str.contains("anodic") == True)
                & (Cdl_fit["Cdl_R"] > 0.85),
                :,
            ]
            #           ====# UGLY STUFF FOR EMPTY Anodic Cdl scan ===
            if not Cdl_an_slice.empty:
                Cdl_an_corr = linregress(Cdl_an_slice[EvRHE], Cdl_an_slice["Cdl"])
                Cdl_an_slice = Cdl_an_slice.assign(
                    **{
                        "Cdl_corr": Cdl_an_slice["Cdl"]
                        - (Cdl_an_slice[EvRHE].values * Cdl_an_corr[0] + Cdl_an_corr[1])
                        + Cdl_an_slice["Cdl"].mean()
                    }
                )
            elif Cdl_an_slice.empty:
                Cdl_an_corr = linregress(Cdl_cath_slice[EvRHE], Cdl_cath_slice["Cdl"])
                Cdl_an_slice = Cdl_cath_slice
                Cdl_an_slice = Cdl_an_slice.assign(**{"Cdl_corr": 0})
            #                               columns=['SampleID','Sweep',EvRHE,'Cdl','Cdl_std','Cdl_R'])
            Cdlprs = {
                "Cdl_cath_max": Cdl_cath_slice.loc[
                    Cdl_cath_slice["Cdl_corr"].idxmax(), :
                ],
                "Cdl_cath_mean": Cdl_cath_slice.loc[
                    :, ["Cdl_corr", "Cdl", "Cdl_R"]
                ].mean(),
                "Cdl_an_max": Cdl_an_slice.loc[Cdl_an_slice["Cdl_corr"].idxmax(), :],
                "Cdl_an_mean": Cdl_an_slice.loc[:, ["Cdl_corr", "Cdl", "Cdl_R"]].mean(),
                "Analysis_date": datetime.now(),
                "SampleID": SampleID,
                "PAR_file": N2_FileName,
            }
            Cdl_Eout = []
            for Et in np.linspace(0, 1, 21):
                try:
                    Cdl_Eout.append(
                        Cdl_cath_slice.loc[
                            (np.isclose(Cdl_cath_slice[EvRHE], Et, atol=50e-3) == True),
                            :,
                        ].head(1)
                    )
                    Cdl_Eout.append(
                        Cdl_an_slice.loc[
                            (np.isclose(Cdl_an_slice[EvRHE], Et, atol=50e-3) == True), :
                        ].head(1)
                    )
                except Exception as e:
                    logger.error("N2 Scan_Cdl Error at {0}: {1}".format(Et, e))
            Cdl_PARS = pd.concat([i for i in Cdl_Eout], axis=0)
            Cdl_PARS_ovv = pd.merge(Cdl_PARS, N2_ovv_file)

            Cdl_fn_pars_base = N2_dest_dir.joinpath("Cdl_pars_%s.xlsx" % N2_fn)
            Cdl_fn_pars_out = FileOperations.CompareHashDFexport(
                Cdl_PARS_ovv, Cdl_fn_pars_base
            )
            #            index_info = {'PAR_file': N2_FileName, 'DestFile': Cdl_fn_pars_out, 'Type_output': 'N2_Cdl_pars'}
            #            index_out.append(index_info)
            N2_plot_Cdl_sweeptype_scatter(
                SampleID,
                Cdl_cath_slice.ScanRates.unique()[0],
                Cdl_fit,
                Cdl_cath_slice,
                Cdl_an_slice,
                N2_dest_dir,
                N2_fn,
            )

        elif len(ScanRates) <= 1:
            Cdl_PARS = pd.DataFrame([])
            Cdl_fit = pd.DataFrame([])
            logger.warning("N2_scans few ScanRates {0},{1}".format(ScanRates, N2_fn))

        # FIRST LAST N2 SCANS#
        preN2_scan = grA.get_group(
            ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.min())
        )
        cls_fastN2_scan = grA.get_group(
            ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.max())
        )["Segment #"].unique()
        fastN2_scan = grA.get_group(
            ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.max())
        )
        gr_cls_fastN2_scan = fastN2_scan.groupby(by="Segment #")
        first_fastN2_scan, last_fastN2_scan = gr_cls_fastN2_scan.get_group(
            cls_fastN2_scan[0]
        ), gr_cls_fastN2_scan.get_group(
            cls_fastN2_scan[int(len(cls_fastN2_scan) * 0.85) + 1]
        )
    #        while len(first_fastN2_scan) > 2000 and len(last_fastN2_scan) > 2000:
    #            first_fastN2_scan, last_fastN2_scan = first_fastN2_scan.drop_duplicates(), last_fastN2_scan.drop_duplicates()
    #            first_fastN2_scan, last_fastN2_scan = first_fastN2_scan.drop_duplicates(), last_fastN2_scan.drop_duplicates()
    #            first_fastN2_scan, last_fastN2_scan = first_fastN2_scan.loc[
    #                                                  first_fastN2_scan.PAR_file == first_fastN2_scan['PAR_file'].unique()[
    #                                                      -1], :], last_fastN2_scan.loc[last_fastN2_scan.PAR_file ==
    #                                                                                    last_fastN2_scan[
    #                                                                                        'PAR_file'].unique()[-1], :]
    #                Last_fastN2_scan.get_group()
    def _except():
        logger.warning("N2 scan failure: %s" % e)
        preN2_scan = grA.get_group(
            ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.min())
        )
        cls_fastN2_scan = grA.get_group(
            ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.max())
        )["Segment #"].unique()
        fastN2_scan = grA.get_group(
            ("N2", "Cyclic Voltammetry (Multiple Cycles)", "N2_act", ScanRates.max())
        )
        gr_cls_fastN2_scan = fastN2_scan.groupby(by="Segment #")
        first_fastN2_scan, last_fastN2_scan = gr_cls_fastN2_scan.get_group(
            cls_fastN2_scan[0]
        ), gr_cls_fastN2_scan.get_group(
            cls_fastN2_scan[int(len(cls_fastN2_scan) * 0.85) + 1]
        )
    #        while len(first_fastN2_scan) > 2000 and len(last_fastN2_scan) > 2000:
    #            first_fastN2_scan, last_fastN2_scan = first_fastN2_scan.drop_duplicates(), last_fastN2_scan.drop_duplicates()
    #            first_fastN2_scan, last_fastN2_scan = first_fastN2_scan.drop_duplicates(), last_fastN2_scan.drop_duplicates()
    #            first_fastN2_scan, last_fastN2_scan = first_fastN2_scan.loc[
    #                                                  first_fastN2_scan.PAR_file == first_fastN2_scan['PAR_file'].unique()[
    #                                                      -1], :], last_fastN2_scan.loc[last_fastN2_scan.PAR_file ==
    #                                                                                    last_fastN2_scan[
    #                                                                                        'PAR_file'].unique()[-1], :]
    #                Last_fastN2_scan.get_group()


def plot_N2_background():
    try:
        fig, ax = plt.subplots()
        first_fastN2_scan.plot(
            x="E_AppV_RHE",
            y="jmAcm-2",
            xlim=[-0.2, 1.2],
            ylim=[-7, 7],
            ax=ax,
            label="1st N2 scan",
        )
        last_fastN2_scan.plot(
            x="E_AppV_RHE",
            y="jmAcm-2",
            xlim=[-0.2, 1.2],
            ax=ax,
            label="last N2 scan",
            title=(N2_scan.PAR_file.unique()[0]),
        )
        if len(N2_scan) < 2010 or len(N2_scan) > 1990:
            N2_scan.plot(
                x="E_AppV_RHE",
                y="jmAcm-2",
                xlim=[-0.2, 1.2],
                ylim=[-10, 10],
                ax=ax,
                label="N2 scan",
            )
        #        first_fastN2_scan.to_csv(N2_dest_dir.joinpath('N2_first.csv')) , last_fastN2_scan.to_csv(N2_dest_dir.joinpath('N2_last.csv'))
        #              N2_scan.plot(x='E_AppV_RHE',y='j A/cm2', xlim=[-0.2,1.2], ylim=[-0.02,0.02],title=os.path.basename(N2_scan.PAR_file.unique()[0]),ax=ax)
        plt.savefig(
            N2_dest_dir.joinpath("%s.png" % N2_fn), dpi=100, bbox_inches="tight"
        )
        plt.close()
    except Exception as e:
        logger.warning("No N2 plot: %s" % e)

        # %%
def prepare_data_for_export():
    seg = N2_scan["Segment #"].unique()[0]
    Sg_info = [
        (i, N2_scan[i].unique())
        for i in N2_scan.columns
        if len(N2_scan[i].unique()) < 30
    ]
    Sg_info_Out = pd.DataFrame(columns=[i[0] for i in Sg_info])
    Sg_info_Out.loc[seg, [i[0] for i in Sg_info]] = [i[1] for i in Sg_info]
    Sg_info_Out["Seg"] = float(seg)
    #                dlo.loc[dlo.Seg == seg,:]
    #    print(N2_scan.columns)
    N2_cols_drop_option = [
        "Segment #",
        "Point #",
        "Current Range",
        "Status",
        "Frequency(Hz)",
        "Z Real",
        "Z Imag",
        "ActionId",
        "AC Amplitude",
        "jcorr",
        "Gas",
        "PAR_exp",
        "j_ring",
        "RPM",
        "Comment",
        "Measured_OCP",
        "Electrode",
    ]
    N2_cols_drop_option_set = [i for i in N2_cols_drop_option if i in N2_scan.columns]
    Other_cols = [
        "E(V)",
        "I(A)",
        "Elapsed Time(s)",
        "RHE_OCP",
        "E_AppV_RHE",
        "E_Applied_VRHE",
        "j A/cm2",
        "jmAcm-2",
        "Gas",
        "PAR_exp",
        "j_ring",
        "RPM",
        "Comment",
        "Measured_OCP",
        "pH",
        "Electrolyte",
        "ScanRate_calc",
        "SampleID",
        "PAR_file",
        "BaseName",
        "hash",
        "Instrument",
        "DATE",
        "EvRHE_diff",
        "DestFile",
        "Type_action",
        "Sweep_Type",
        "Scanrate",
    ]
    #    drop [i[0] for i in Sg_info if not 'Sweep_Type' in i[0]]
    try:
        SgOut = N2_scan.drop(columns=N2_cols_drop_option_set)
    except Exception as e:
        logger.warning(
            "N2 scan background problem with dropping columns: {0},\n columns:{1}".format(
                e, N2_scan.columns
            )
        )
        SgOut = N2_scan
    #    N2_dest_dir = Path(gr_N2_ovv.Dest_dir.iloc[0]).joinpath('N2_scans/{0}'.format('_'.join([Electrolyte,SampleID,Loading_name])))
    #    N2_dest_dir.mkdir(parents=True,exist_ok=True)

    logger.info("N2 scan background exported to: %s" % N2_act_BG)
    SgOut.to_pickle(N2_act_BG)
    #    Cdl_PARS = Cdl_PARS.assign(**{})

    #    N2_act_BG = Path(N2_dest_dir.parent).joinpath(N2_fn + '_BG.xlsx')
    #    logger.info('N2 scan background exported to: %s' % N2_act_BG)
    #    SgOut.to_excel(N2_act_BG)
    #    index_info = {'PAR_file': N2_FileName, 'DestFile': N2_act_BG, 'Type_output': 'N2_bg'}
    #    index_out.append(index_info)

    # %%
    return N2_scan, Cdl_fit, Cdl_PARS
