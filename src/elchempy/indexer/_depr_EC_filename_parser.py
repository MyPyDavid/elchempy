"""
Created on Thu Aug 26 16:14:57 2021

@author: DW
"""


class PAR_file_parser:
    """
    Parses a filename

    Parameters
    ----------
    PARf : Path
        a Path object to a .par file

    Returns
    ---------
        parse_result : dict
    """

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
                    {
                        "PAR_file": self.PARf,
                        "_parser_fail": True,
                        "_parser_message": errm,
                    }
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
            if up_dir == -(len(self.PARf.parts) - 1):
                exp_date = self.parse_result_stat["stat_ctime"]
                logger.warning(
                    "EXP date classifier error setting stat ctime: for {0}".format(
                        exp_dir
                    )
                )
        self.parse_result.update(
            {
                "EXP_date": exp_date,
                "EXP_date_day_dt": dt.datetime.fromisoformat(str(exp_date)).date(),
            }
        )

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

    def _get_match_PAR_EXP_date(self, exp_date, PAR_date):
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


def _depr_collect_parse_results(
    self, read_data=False, store_data=False, **kwargs
) -> Dict:
    """performs all the steps for parsing the filepath"""
    parse_res_collect = {}

    self.stats_ = self.stat()

    _fnID = self.make_dict_from_keys(
        index_file_primary_keys, (self.get_rfID_from_path(self),)
    )
    _filepath = self.make_dict_from_keys(index_file_path_keys, (self.stem, self, self))
    # _sample = self.parse_sample_with_checks()
    _filestats = self.parse_filestats(self.stats_)
    if read_data == True:
        data = self.read_data(self)

        parse_res_collect = {**_fnID, **_filepath, **_filestats}
    else:
        logger.warning(f"{self._qcnm} {self} is not a file => skipped")
    # else:
    # logger.warning(f"{self._qcnm} {self} does not exist => skipped")
    return parse_res_collect


#%%


def _index_structure():
    dr = {
        "EXP_dir": exp_dir_path,
        "Dest_dir": Dest_dir,
        "EXP_date": exp_date,
        "PAR_date": PAR_date,
        "EXP_PAR_date_match": exp_DATE_match,
        "EXP_PAR_folder_match": sID_dict["sID_Match_file_dir"],
        "PAR_file": PARf,
        "PAR_hash": hash_f,
        "Gas": gas,
        "PAR_exp": tpnm,
        "postAST": postAST,
        "filesize": PARf.stat().st_size,
        "basename": basepf,
        "SampleID": sID_file,
        "SampleID_folder": sID_dict["sID_Folder"],
        "Creation_date": cr_date,
        "LastMod_date": mod_date,
        "Delta_PAR_LastMod": PAR_cr_Delta,
        "Electrode": electrode,
        "pH": pH["pH"][0],
        "Electrolyte": pH["Electrolyte"][0],
        "Comment": comment_act0_DF,
        "Loading_name": loading_name,
        "Loading_cm2": loading_cm2,
    }


index_file_primary_keys = {"fID": "string"}

index_file_path_keys = {"FileStem": "string", "FilePath": "Path", "PAR_file": "Path"}

index_folder_path_keys = {
    "DIR_name": "Path",
    "DIR_dest_name": "Path",
    "DIR_date": "datetime.date",
}

index_file_date_keys = {
    "PAR_date": "datetime.date",
    "PAR_introspec_data": "datetime.date",
}

index_file_sample_keys = {
    "SampleID": "string",
    "SampleGroup": "string",
}

index_file_read_text_keys = {"FileHash": "string", "FileText": "string"}

index_dtypes_collection = {
    **index_file_path_keys,
    **index_file_sample_keys,
    **index_file_read_text_keys,
}

# Extra name to sID mapper, if keys is in filename
# _extra_sID_name_mapper = {
#     "David": "DW",
#     "stephen": "SP",
#     "Alish": "AS",
#     "Aish": "AS"}
# Extra name to sID mapper, if key is in filepath parts
# _extra_sgrpID_name_mapper = {"Raman Data for fitting David": "SH"}


def _def_parse():
    for f in LOCAL_FILES:
        fp = FilePathParser(f)


def _depr_parse_filestats(self, fstat) -> Dict:
    """get status metadata from a file"""

    filestats = get_fstats(fstat)
    return self.make_dict_from_keys(index_file_stat_keys, filestats)


def _depr_make_dict_from_keys(self, _keys_attr: Dict, _result: tuple) -> Dict:
    """returns dict from tuples of keys and results"""
    if not isinstance(_result, tuple):
        logger.warning(
            f"{self._qcnm} input value is not a tuple, {_result}. Try to cast into tuple"
        )
        _result = (_result,)

    _keys = _keys_attr.keys()

    if not len(_result) == len(_keys) and not isinstance(_keys, str):
        # if len not matches make stand in numbered keys
        _keys = [f"{_keys_attr}_{n}" for n, i in enumerate(_result)]
    return dict(zip(_keys, _result))


def _depr_parse_sample_with_checks(self):
    """parse the sID and sgrpID from stem"""

    _parse_res = filestem_to_sid_and_pos(self.stem)

    if len(_parse_res) == 2:
        sID, position = _parse_res

        try:
            sID = _extra_overwrite_sID_from_mapper(sID)
        except Exception as exc:
            logger.info(
                f"{self._qcnm} {self} _extra_overwrite_sID_from_mapper failed => skipped.\n{exc}"
            )

        sgrpID = sID_to_sgrpID(sID)

        try:
            sgrpID = _extra_overwrite_sgrpID_from_parts(self.parts, sgrpID)
        except Exception as exc:
            logger.info(
                f"{self._qcnm} {self} _extra_overwrite_sgrpID_from_parts failed => skipped.\n{exc}"
            )

        _parse_res = sID, position, sgrpID
    else:
        logger.warning(
            f"{self._qcnm} {self} failed to parse filename to sID and position."
        )
    return self.make_dict_from_keys(index_file_sample_keys, _parse_res)


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
