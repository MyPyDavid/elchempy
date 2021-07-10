#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 17:26:54 2020

@author: DW

"""
# import pathlib
import datetime as dt
from pathlib import Path
import pandas as pd
from bs4 import BeautifulSoup

# from functools import partial
# import multiprocessing

# import os
# from os import cpu_count
# import hashlib
# import lxml
# import logging
# import numpy as np
# from FindExpFolder import *
# FindExpFolder('VERSASTAT').PARfiles
# import concurrent.futures
# from multiprocessing import Pool

from file_py_helper.find_folders import FindExpFolder
from file_py_helper.file_functions import FileOperations
from file_py_helper.FindSampleID import GetSampleID

# import lxml
import sys

import logging

logger = logging.getLogger(__name__)
# tt = PAR_file_parser(cc.par_files_run[42])

# t1 = '/mnt/DATA/EKTS_CloudStation/CloudStation/Experimental data/Raw_data/VERSASTAT/2018-01-Jan/19.01.2018_DW18_bipot_0.5MH2SO4_RRDE22434/OCP_RHE.par'
# tt2 = PAR_file_parser(t1)


class PAR_file_parser:
    """Input is a filename
    output self.parse_result is a dict"""

    def __init__(self, PARf, destdir=""):
        self.destdir = destdir
        self._check_destdir()

        self.PARf = PARf
        self._check_path()
        if not self.parse_result.get("_parser_fail", False):
            self._get_stat_results()
            self._file_or_dir()
            self._get_SampleID()
            self._get_Sample_date_from_filename()
            self._get_Gas_from_filename()
            self._get_postAST_from_filename()
            self._get_pH_from_filename()
            self._get_Loading_from_filename()
            self._get_Electrode_from_filename()
            self._get_PAR_introspector()

    #            self.EC_PAR_file_check()

    def _check_destdir(self):
        if not self.destdir:
            try:
                self.DestDir = FindExpFolder("VERSASTAT").DestDir
            except Exception as e:
                errm = (
                    f'PF inspector: FindExpFolder("VERSASTAT").DestDir \n because: {e}'
                )
        else:
            if self.destdir.is_dir():
                self.DestDir = self.destdir
            else:
                errm = f"PF inspector error: {self.destdir} is not a valid dir"
                self.parse_result = {
                    "PAR_file": self.PARf,
                    "_parser_fail": True,
                    "_parser_message": errm,
                }

    def _check_path(self):
        self.parse_result = {}
        if not "Path" in str(type(self.PARf)):
            try:
                self.PARf = Path(self.PARf)

            except Exception as e:
                errm = f"not valid filepath: {self.PARf}\n because: {e}"
                self.parse_result.update(
                    {"PAR_file": self.PARf, "_parser_fail": True, "_parser_message": errm}
                )
        self.parse_result.update({"PAR_file": self.PARf, "basename": self.PARf.stem})

    def _file_or_dir(self):
        if self.PARf.is_file():
            exp_dir = self.PARf.parent.parts[-1]
            exp_dir_path = self.PARf.parent
            self.file_type = "file"

        elif self.PARf.is_dir():
            exp_dir = self.PARf.name
            exp_dir_path = self.PARf.parent
            self.file_type = "dir"
        else:
            logger.error(
                "EC PAR file check. PAR file: {0} is neither file nor dir?!".format(
                    self.PARf
                )
            )
            exp_dir, exp_dir_path = self.PARf.parent.parts[-1], self.PARf.parent
        self.exp_dir, self.exp_dir_path = exp_dir, exp_dir_path

        try:
            dest_prts, exp_prts = self.DestDir.parts, self.exp_dir_path.parts
            diffeldl = len(exp_prts) - len(dest_prts)
            if diffeldl > 0:
                Dest_dir = Path(*(dest_prts + exp_prts[-diffeldl::]))
            else:
                Dest_dir = self.DestDir / "EC_plots" / diffeldl[-1]
        except Exception as e:
            print(e, exp_dir)
            Dest_dir = self.DestDir / "EC_plots" / exp_dir.parts[-1]
        self.Exp_Dest_dir = Dest_dir

        self.parse_result.update(
            {
                "EXP_dir": self.exp_dir_path,
                "Dest_dir": self.Exp_Dest_dir,
                "EXP_dir_name": self.exp_dir,
            }
        )

    def _get_stat_results(self):
        _stat = self.PARf.stat()
        fields = "mode ino dev nlink uid gid size atime mtime ctime blocks blksize rdev flags gen birthtime ftype attrs obtype"
        d = dict(
            (f"stat_{field}", getattr(_stat, "st_" + field, None))
            for field in fields.split()
        )
        for i in ["stat_mtime", "stat_ctime", "stat_atime"]:
            d.update({i: dt.datetime.fromtimestamp(d[i])})

        d = {key: val for key, val in d.items() if val}
        self.parse_result_stat = d

    def _get_stem_file(self):
        self.stem = self.PARf.stem
        self.split = self.stem.split("_")

    def _get_SampleID(self):
        try:
            method = GetSampleID(self.PARf)
            sID_file, sID_dict = method.SampleID, method.SampleID_dict
        except Exception as e:
            logger.error(
                "EC PAR file check. PAR file: {0} FindSampleID error {1}".format(
                    self.PARf, e
                )
            )
            method = GetSampleID(self.PARf.parent)
            sID_file, sID_dict = method.SampleID, method.SampleID_dict
        self.parse_result.update(
            {
                key: val
                for key, val in sID_dict.items()
                if any(
                    [
                        i == key
                        for i in [
                            "SampleID",
                            "sID_Folder",
                            "sID_Match_file_dir",
                            "EXP_PAR_folder_match",
                        ]
                    ]
                )
            }
        )

    #        self.sID,self.sID_dict = sID_file, sID_dict

    def _get_Sample_date_from_filename(self):
        #    def getting_date_from_PARf(PARf,exp_dir):
        exp_date = ""
        up_dir = -2
        while exp_date == "":
            #                print(up_dir)
            try:
                exp_date = dt.datetime.strptime(self.exp_dir.split("_")[0], "%d.%m.%Y")
            #                dt.date.fromisoformat(self.exp_dir.split('_')[0],format='%d.%m.%Y',errors='raise')
            #            exp_date = pd.to_datetime(exp_dir.split('_')[0],errors='ignore')
            except Exception as e1:
                #                faillst1.append((PARf,exp_dir,exp_date))
                try:
                    exp_date = pd.to_datetime(
                        "_".join(self.exp_dir.split("_")[0:3]),
                        format="%Y_%m_%d",
                        errors="raise",
                    )
                except Exception as e2:
                    #                exp_date = pd.to_datetime(exp_dir.split('_')[0],format='%d.%m.%Y',errors='raise')
                    exp_dir = self.PARf.parts[-up_dir]
                    up_dir -= 1
                    logger.warning(
                        "EXP date classifier error: {0} for {1}".format(
                            str(e1) + str(e2), exp_dir
                        )
                    )
            #                        exp_date = pd.to_datetime(PARf.stat().st_ctime,unit='s')
            #                    faillst2.append((PARf,exp_dir,exp_date))
            if up_dir == -(len(self.PARf.parts) - 1):
                exp_date = self.parse_result_stat["stat_ctime"]
                logger.warning(
                    "EXP date classifier error setting stat ctime: for {0}".format(
                        exp_dir
                    )
                )
        #                 faillst3.append((PARf,exp_dir,exp_date))

        self.parse_result.update(
            {
                "EXP_date": exp_date,
                "EXP_date_day_dt": dt.datetime.fromisoformat(str(exp_date)).date(),
            }
        )

    #        try:
    #            self.date = getting_date_from_PARf(self.PARf,self.exp_dir)
    #        except:
    #            self.date = 0

    def _get_Gas_from_filename(self):
        try:
            Exp_gas_type_name = GetSampleID.determine_Gas_Exp_from_filename(
                self.PARf.stem, "", False
            )
        except Exception as e:
            Exp_gas_type_name = {"Gas": "gas_error", "PAR_exp": "type_error"}
        self.parse_result.update(Exp_gas_type_name)

    #        self.Exp_gas,self.Exp_type_name = Exp_gas,Exp_type_name
    def _get_postAST_from_filename(self):
        try:
            postAST = GetSampleID.determine_postAST_from_filename(
                self.PARf, verbose=False
            )
        except Exception as e:
            postAST = {"postAST": "postAST_error"}
        self.parse_result.update(postAST)

    def _get_pH_from_filename(self):
        #        pf.split('\\')[-1]
        try:
            pH = GetSampleID.determine_pH_from_filename(self.PARf, False)
            if pH["Electrolyte"][0] == "FALSE":
                pH = GetSampleID.determine_pH_from_filename(self.exp_dir, False)
        except Exception as e:
            pH = {"pH": "pH_error", "Electrolyte": 0}
        #        {'pH' : pH['pH'][0],'Electrolyte' : pH['Electrolyte'][0]}
        self.parse_result.update(pH)

    def _get_Loading_from_filename(self):
        _loading_date = self.parse_result.get("EXP_date")
        try:
            _loading_result = GetSampleID.ink_loading_from_filename(
                self.PARf, _loading_date
            )
        except Exception as e:
            loading_name, loading_cm2 = f"loading_error, {e}", 0
            _loading_result = {
                "Loading_name": loading_name,
                "Loading_cm2": loading_cm2,
                "Loading_date": _loading_date,
            }
        self.parse_result.update(_loading_result)

    #        if 'ring' in basepf or 'Ring' in basepf:
    #            electrode = 'Pt_ring'
    def _get_Electrode_from_filename(self):
        if any([n in self.PARf.stem for n in ["Pt-ring", "ring", "Ring"]]):
            electrode = "Pt_ring"
        elif any([n in self.PARf.stem for n in ["disk", "Disk"]]):
            electrode = "disk"
        else:
            electrode = "Electrode_unknown"
        self.electrode = electrode
        self.parse_result.update({"Electrode": electrode})

    #        exps = ['RHE|HER_OCP|RHE_OCP|OCP_RHE','HPRR','N2','O2','_V3',' _V3F','_disk','_ring']
    #        exp_names = ['RHE','HPRR','N2','O2',' V3',' V3F','disk','ring']

    def _get_PAR_introspector(self):
        pp = PAR_file_introspector(self.PARf)

        try:
            if self.parse_result.get("EXP_date_day_dt", 1) == pp.read_result.get(
                "PAR_date_day_dt", 0
            ):
                exp_DATE_match = "yes"
            else:
                exp_DATE_match = "no"
        #            'EXP_PAR_date_match' : exp_DATE_match
        except:
            exp_DATE_match = "error"
            pass
        pp.read_result.update({"EXP_PAR_date_match": exp_DATE_match})
        self.parse_result.update(pp.read_result)


class PAR_file_introspector:
    def __init__(self, PARf):
        self.PARf = PARf
        self._get_hash()
        self.read_PAR_file()
        if self.PAR_soup != "error":
            self.read_segments()

    def _get_match_PAR_EXP_date(self, exp_date,PAR_date):
        try:
            if exp_date.date() == PAR_date.date():
                exp_DATE_match = "yes"
            else:
                exp_DATE_match = "no"
        except Exception as e:
            logger.warning(
                "EC_PAR_file_check {0} match date error: {1}".format(self.PARf, e)
            )
            #            print(e)
            exp_DATE_match = "NaN"

    def _get_hash(self):
        hash_f = FileOperations.take_hash(self.PARf)
        self.hash = hash_f

    def read_PAR_file(self):
        try:
            PAR_soup = BeautifulSoup(self.PARf.read_text(), features="html.parser")
        except:
            PAR_soup = "error"
            self.read_result = {"read_PAR_file": "error", "PAR_hash": self.hash}
        self.PAR_soup = PAR_soup

    def read_segments(self):
        _read_result = {"PAR_hash": self.hash}
        _read_result.update(
            self._get_DF_soup_part(self.PAR_soup.action0, prefix="_act0_")
        )

        _str_OCP = _read_result.get("_act0_Measured Open Circuit", "0 mV")
        try:
            _fl_OCP = float(_str_OCP.split()[0])
        except:
            _fl_OCP = 0
        _read_result.update({"Measured_OCP": _fl_OCP})

        _read_result.update(
            self._get_DF_soup_part(self.PAR_soup.instrument, prefix="_instr_")
        )

        _read_result.update(
            self._get_DF_soup_part(self.PAR_soup.application, prefix="_appl_")
        )
        #        if action0:
        #            _read_result.update({f'_act0_{key}' : val for key,val in action0.items()})
        _read_result.update(
            self._get_DF_soup_part(self.PAR_soup.experiment, prefix="_exp_")
        )
        try:
            _read_result.get("_exp_DateAcquired")
            PAR_date = pd.to_datetime(
                _read_result.get("_exp_DateAcquired")
                + " "
                + _read_result.get("_exp_TimeAcquired")
            )
            _PAR_date_src = "_exp"
        except Exception as e:
            #            logger.warning('EC_PAR_file_check {0} error: {1}'.format(self.PARf,e))
            PAR_date = pd.to_datetime(self.PARf.stat().st_ctime, unit="s")
            _PAR_date_src = "_stat_ctime"

        _PAR_date = {"PAR_date": PAR_date, "PAR_date_src": _PAR_date_src}
        try:
            _PAR_date.update(
                {
                    "PAR_date_day": dt.datetime.strftime(PAR_date, format="%Y-%m-%d"),
                    "PAR_date_min": dt.datetime.strftime(
                        PAR_date, format="%Y-%m-%dT%H:%M"
                    ),
                    "PAR_date_day_dt": dt.datetime.fromisoformat(str(PAR_date)).date(),
                }
            )
        except:
            pass

        _read_result.update(_PAR_date)

        self.read_result = _read_result

    @staticmethod
    def _get_DF_soup_part(soup_segment, prefix=""):
        #    Soup_segment,hash_f = N2_soup_par.find_all(name=actie)[0],actie#    = N2_soup_par.instrument.string
        out = {}
        if soup_segment:
            splt = []
            for ch in soup_segment.children:
                try:
                    #                splt = [i.split('=') for i in soup_segment.text.split('\n')[1:-1]]
                    splt += PAR_file_introspector._clean_up_comment(ch)
                #                    out = {f'{prefix}{i[0]}' : i[1]  for i in splt}
                #                out = dict([(i[0],i[1]) for i in splt])
                except TypeError:
                    if prefix == "_instr_":
                        splt += [i.split("=") for i in ch.text.split("\n")[1:-1]]
                        splt += ["".join(list(ch.attrs)).split(":"), ch.name.split(":")]
                except Exception as e:
                    splt += [["soup_error", e]]
        #                    splt = [i.split('=') for i in soup_segment.text.split('\n')[1:-1]]

        out = {f"{prefix}{i[0]}": i[1] for i in splt}
        return out

    @staticmethod
    def _clean_up_comment(ch):

        splt_clean = [i for i in ch.split("\n")[1:-1] if i]
        _comment = [
            (n, i) for n, i in enumerate(splt_clean) if i.startswith("Comment=")
        ]
        if _comment:
            _new_comment = ", ".join(str(i) for i in splt_clean[_comment[0][0] :])
            _new_splt = splt_clean[: _comment[0][0]] + [_new_comment]
            return [i.split("=") for i in _new_splt]
        else:
            return [i.split("=") for i in splt_clean]


#       soup_segment=
#       for ch in self.PAR_soup.instrument.children:
#           print(ch)
#           splt = PAR_file_introspector._clean_up_comment(soup_segment)
#           out = {f'{prefix}{i[0]}' : i[1]  for i in splt}
#        PARf_text = self.PARf.read_text()
#        p = etree.XMLParser(target = etree.TreeBuilder())
#        p = etree.XMLParser(huge_tree=True, remove_blank_text = True,ns_clean = True)
#        tree = etree.parse(str(self.PARf))
#        etree.tostring(tree)
#        root = tree.getroot()
#        etree.tostring(tree.getroot())
