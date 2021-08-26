"""Collection of helper methods for parsing a filename"""
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import Counter

import datetime
from typing import Tuple, List, Dict, Union

import re

import logging
logger = logging.getLogger(__name__)



from elchempy.indexer.filepath_parser import FilePathParser
from elchempy.indexer.extra_EC_info import loading_ref, WE_surface_area_cm2

### for Developing
from elchempy.config import LOCAL_FILES

# class ParserMethods

def _dev():
    f = LOCAL_FILES[-1]
    sid = InterpretFilePath(f)
    self = sid

class InterpretFilePath(FilePathParser):
    """Find or Guess the SampleID that is assiociated with the filename"""

    name_separators = ["_", "-"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.split_info = {}

        self.name_split = self.get_most_common_split(self.stem)
        self.parent_split = self.get_most_common_split(self.parent.name)

        self.EC_info = {}
        stem_info  = all_functions('stem', self.stem, self.name_split)
        parent_info = all_functions('parent',self.parent.name, self.parent_split)
        self.EC_info = {**stem_info, **parent_info}

    def get_most_common_split(self, name):
        sepcounter = Counter([i for i in name if i in self.name_separators])
        self.split_info.update(**{ name: { 'counter' : sepcounter}})
        if sepcounter:
            sep, _nsep = sepcounter.most_common()[0]
            return name.split(sep)
        else:
            return []

    def recognize_element(self, elem):
        sid = self.try_find_sampleID(elem)


    def try_find_sampleID(self, file):
        """Input string or Path and output is SampleID[0] and dict"""
        pf = Path(file)
        try:
            pf_Match = self.match_SampleID(pf)
        except:
            pf_Match = "NoneNoneNone"
        if pf.is_file():
            parent_folder = pf.parent.parts[-1]

            try:
                pf_StemMatch = self.match_SampleID(pf.stem)
            except:
                pf_StemMatch = "None"
            try:
                pf_ParentMatch = self.match_SampleID(parent_folder)
            except:
                pf_ParentMatch = "None"
            if pf_StemMatch == pf_ParentMatch:
                SampleID, match = pf_StemMatch, "yes"
            else:
                SampleID, match = pf_StemMatch, "no"
        elif pf.is_dir():
            parent_folder = pf.name
            SampleID, match = pf_Match, "dir"
            pf_StemMatch, pf_ParentMatch = "None", "None"
        else:
            parent_folder = pf.name
            SampleID, match = pf_Match, "other"
            pf_StemMatch, pf_ParentMatch = self.match_SampleID(pf.stem), "None"
        #             pf_StemMatch,pf_ParentMatch = 'None', 'None'
        #        if len(pf_Match) < 10:
        #            SampleID = pf_Match
        sID_MatchOut = {
            "SampleID": SampleID,
            "sID_Filename": pf_StemMatch,
            "Tested_Folder": parent_folder,
            "sID_Folder": pf_ParentMatch,
            "sID_Match_file_dir": match,
            "sID_FullPath": pf_Match,
            "Tested_File": file,
        }
        return SampleID, sID_MatchOut


skipped_Electrolytes = ["H2SO4", "HClO4", "MeOH", "KOH", "NaOH", "H2O2"]
skipped_Gas_Exp_parts = [
    "N2",
    "O2",
    "ORR",
    "RHE",
    "CV",
    "AST",
    "HPRR",
    "OER",
    "HER",
    "EIS",
    "LC",
    "SST",
    "sHA",
    "20CLS",
]

def all_functions(name, fname, fsplit) -> Dict:

    EC_info = {}

    elec = determine_electrode_name_from_stem(fname)
    # for nameparent in ['name', 'parent']:
    # nameparent
    # fname = getattr(self, nameparent )
    if not isinstance(fname, str):
        fname = fname.name
    # fsplit = getattr(self,f'{nameparent}_split')
    date = determine_date_from_filename(fsplit)
    gas_exp = determine_Gas_Exp_from_filename(fname)
    pH = determine_pH_from_filename(fname)
    postAST = determine_postAST_from_filename(fname)

    EC_info.update({name : {**{'fname' : fname}, **date, **gas_exp, **pH, **postAST } })

    return EC_info
        # determine_date_from_filename(self.parent_split)




# @staticmethod
def match_SampleID(file, message=False, include_extra=False):
    """finds a possible match in the filename with a known sample id"""
    #        message = True
    formatID_matches, sampleID = [], []

    # if '-' in file.name and not '_' in file.name
    # if not '_' in file and '-' in file:
    # filename = file.name.replace('-','_')
    PF0 = Path(file)

    if PF0.is_dir():
        FileParts_raw = PF0.parts[-1].replace(" ", "_").split("_")
    #            FileN = PF0.parts[-1]
    elif PF0.is_file():
        FileParts_raw = PF0.name.replace(" ", "_").split("_")
    #            FileN = PF0.stem
    else:
        FileParts_raw = PF0.name.replace(" ", "_").split("_")
        if len(FileParts_raw) == 1:
            FileParts_raw = PF0.stem.split()
    #            FileN = PF0.stem
    skip_parts = skipped_Electrolytes + skipped_Gas_Exp_parts
    FileParts = [
        i.upper()
        for i in FileParts_raw
        if i.upper() not in skip_parts and re.search("[a-zA-Z]", i)
    ]

    while sampleID == []:
        Pt_ring_match = [
            re.search("(?<=(pt)){0,1}.{0,1}ring{0,1}", i, re.IGNORECASE)
            for i in FileParts
        ]
        if any(Pt_ring_match):
            sampleID = "Pt_ring"
            formatID_matches = ["Pt-ring"]
        #                if not 'RING' in FileParts:
        #                    formatID_matches = [i for i in FileParts if re.match('^PT',i,re.IGNORECASE)]
        #                elif 'RING' in FileParts:
        #                    formatID_matches = [i for i in FileParts if re.match('(Pt-ring)',i,re.IGNORECASE)]
        #                    sampleID = 'Pt-ring'
        elif any(
            [i for i in FileParts if re.match("(Pt-ring|pt-ring)", i, re.IGNORECASE)]
        ):
            formatID_matches = ["Pt-ring"]
            sampleID = "Pt_ring"
        else:
            if (
                any([i for i in FileParts if re.match("(DWIM)", i, re.IGNORECASE)])
                and not "Pt-ring" in formatID_matches
            ):
                formatID_matches = [i for i in FileParts if re.match("(DWIM)", i)]
            elif any([i for i in FileParts if re.match("(SD)", i)]):
                formatID_matches = [
                    i for i in FileParts if re.match("(SD)", i, re.IGNORECASE)
                ]
            #        elif len(fileSplit) > 1 and not re.search('(KOH|NaOH|H2SO4|HClO4)',fileSplit[-1]):
            #            fs2 = fileSplit[0].split('.')
            #            if len(fs2) > 1:
            #                formatID_matches = [re.match('[A-Z]{2}[0-9]{2}_[0-9]\.[0-9]',file).group()]
            #            else:
            #                formatID_matches = [i for i in file.split('_') if re.match('([A-Z]{2}[0-9]{2,3}[a-z]{0,1})',i)]
            elif "Sandeep" in FileParts:
                formatID_matches = ["SD00"]
            elif any(
                [
                    i
                    for i in FileParts
                    if re.match("(?i)(xx|xxx|DWxx|JOS1xx)", i, re.IGNORECASE)
                ]
            ):
                formatID_matches = ["template"]
            elif "NW" in FileParts:
                formatID_matches = [
                    i for i in FileParts if re.match(".?NW_M.?", i, re.IGNORECASE)
                ]
            elif any("JOS" in s for s in FileParts):
                if any([i for i in FileParts if re.match("(ring)", i)]):
                    formatID_matches = ["Pt-ring"]
                elif any([i for i in FileParts if re.search("(?i)JOS[0-9]{0,3}XX", i)]):
                    formatID_matches = ["template"]
                else:
                    formatID_matches = [
                        i.upper()
                        for i in FileParts
                        if re.match(".?JOS.?", i, re.IGNORECASE)
                    ]
            #            elif any('jos' in s for s in FileParts):
            #                formatID_matches = [i for i in FileParts if re.match('.?jos.?',i)]
            elif any("PT" in s for s in FileParts):
                formatID_matches = [
                    i
                    for i in FileParts
                    if re.match(".?PT(A|B)([0-9]{1,3}).?", i, re.IGNORECASE)
                ]

            elif any(
                [
                    i
                    for i in FileParts
                    if re.match("([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i, re.IGNORECASE)
                ]
            ):
                formatID_matches = [
                    i
                    for i in FileParts
                    if re.match("([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i, re.IGNORECASE)
                ]

            elif any("CELL" in s for s in FileParts):
                formatID_matches = [
                    i
                    for i in FileParts
                    if re.match(".?CELL([0-9]{0,3}).?", i, re.IGNORECASE)
                ]

            else:
                formatID_matches = [
                    i
                    for i in FileParts
                    if re.match("([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i, re.IGNORECASE)
                ]
            #            print(FileParts,formatID_matches)
            try:
                fmtID_split = formatID_matches[0].split("-")
                if len(fmtID_split) > 1:
                    fmtID_split_matches = [
                        i
                        for i in fmtID_split
                        if re.match("(?i)([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i)
                    ]
                    if len(fmtID_split_matches) == 1:
                        formatID_matches = fmtID_split_matches
                    else:
                        formatID_matches = [fmtID_split[0]]
            except:
                pass
        if len(formatID_matches) == 1:
            sampleID = formatID_matches[0]
        elif len(formatID_matches) > 1:
            sampleID = formatID_matches[0]
        elif len(formatID_matches) == 0:
            for t in [
                "P-200",
                "P200",
                "P-250",
                "P-XC",
                "P-XP",
                "VXC72",
                "BP2K",
                "BP2000",
            ]:
                if t in FileParts:
                    sampleID = t
                else:
                    if len(FileParts) > 2:
                        sampleID = FileParts[0]
                    #                        sampleID = file[0:3]
                    elif len(FileParts) == 2:
                        sampleID = file[0:2] + "00"
                    else:
                        sampleID = file
        else:
            sampleID = "_".join(FileParts)

    if include_extra:
        extra_parts = [
            i.upper() for i in FileParts_raw if i.upper() not in sampleID.upper()
        ]
        extra_part_matches = [
            i
            for i in extra_parts
            if re.match("([0-9]|HT|AL|Ir|N2|repeat|[a-zA-Z]{0,1})", i, re.IGNORECASE)
        ]
        if any(extra_part_matches):
            sampleID = "_".join([sampleID, *extra_part_matches])

    if message:
        print('SampleID: "%s" used for %s' % (sampleID, file))
    return str(sampleID.upper())

def determine_pH_from_filename(filename: str) -> Dict:
    ''' determine pH from filename'''

    PAR_file_test = filename
    pH = {}
    while pH == {}:
        if re.search("H2SO4", PAR_file_test) and not re.search(
            "HClO4", PAR_file_test
        ):
            if re.search("MeOH", PAR_file_test):
                if re.search("1M.?MeOH", PAR_file_test):
                    pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+1M-MeOH"}
                elif re.search("0\.1.?MeOH", PAR_file_test):
                    pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+0.1M-MeOH"}
                else:
                    pH = {"pH": 97, "Electrolyte": "xxMH2SO4+MeOH"}
            elif re.search("H2O2", PAR_file_test):
                pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+10mMH2O2"}
            else:
                if re.search("(?<=(0\.1))[a-zA-Z]{0,2}(H2SO4)", PAR_file_test):
                    #                        pH = {'pH' : 1, 'Electrolyte' : '0.1MH2SO4'}
                    if re.search(
                        "(?<=(0\.1))[a-zA-Z]{0,2}(H2SO4)(?!(.{0,5}(MeOH|H2O2)))",
                        PAR_file_test,
                    ):
                        pH = {"pH": 1, "Electrolyte": "0.1MH2SO4"}
                    elif re.search(
                        "(?<=(0\.1))[a-zA-Z]{0,2}(H2SO4)(?=(.{0,5}(H2O2)))",
                        PAR_file_test,
                    ):
                        pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+xxH2O2"}

                elif re.search("(?<=(0\.5))[a-zA-Z]{0,2}(H2SO4)", PAR_file_test):
                    pH = {"pH": 0.3, "Electrolyte": "0.5MH2SO4"}
                else:
                    pH = {"pH": 98, "Electrolyte": "xxMH2SO4"}
        #                elif re.search('H2SO4(?!(.{0,5}(MeOH|H2O2)))',PAR_file_test):
        #                else:
        #                    pH = {'pH' : 96, 'Electrolyte' : 'xxMH2SO4'}
        #                    elif [a for a in pftest_split if '0.1MH2SO4+0.1M-MeOH' in a]:
        #                        pH = {'pH' : 1, 'Electrolyte' : '0.1MH2SO4+0.1M-MeOH'}
        elif re.search("HClO4", PAR_file_test):
            if re.search("(H2O2)", PAR_file_test):
                if re.search("0.1MHClO4(?=(.?10mMH2O2))", PAR_file_test):
                    pH = {"pH": 1, "Electrolyte": "0.1MHClO4+10mMH2O2"}
                else:
                    pH = {"pH": 1, "Electrolyte": "0.1MHClO4+xxH2O2"}

            elif re.search("(?<=(0\.1.))HClO4(?!(.{0,5}H2O2))", PAR_file_test):
                pH = {"pH": 1, "Electrolyte": "0.1MHClO4"}
            else:
                pH = {"pH": 1, "Electrolyte": "xxMHClO4"}
        elif re.search("KOH", PAR_file_test):
            if re.search(
                "(?<=(0\.1))[a-zA-Z]{0,2}(KOH)(?!(.{0,5}H2O2))", PAR_file_test
            ):
                pH = {"pH": 13, "Electrolyte": "0.1MKOH"}
            else:
                pH = {"pH": 13, "Electrolyte": "xxKOH"}
        #                pH = {'pH' : 13, 'Electrolyte' : 'xxKOH'}
        elif re.search("NaOH", PAR_file_test):
            if re.search(
                "(?<=(0\.1))[a-zA-Z]{0,2}(NaOH)(?!(.{0,5}H2O2))", PAR_file_test
            ):
                pH = {"pH": 13, "Electrolyte": "0.NaOH"}
        elif re.search("Acid", PAR_file_test):
            pH = {"pH": 1, "Electrolyte": "0.1MH2SO4_acid"}
        else:
            pH = {"pH": None, "Electrolyte": f"Other_{PAR_file_test}"}

    return pH

def determine_Gas_Exp_from_filename(filename) -> Dict:
    ''' returns a dict with Gas and Exp type from the filename'''
    basepf = filename

    if ('OCP' or 'RHE') in basepf:
        gas, tpnm, exp_match = "N2", "RHE", "yes"
    elif "N2" in basepf:
        gas = "N2"
        if "cls" in basepf and not any(
            [i in basepf for i in ["HER", "OER", "EIS", "ORR", "_AST_"]]
        ):
            tpnm = "N2_act"
        elif all([i in basepf for i in ["cls", "_post"]]):
            tpnm = "N2_act"
        elif "CV" in basepf and not "AST" in basepf:
            tpnm = "N2_act"
        elif (
            "CV" in basepf
            and "AST" in basepf
            and not ("post" in basepf or "pAST" in basepf)
        ):
            tpnm = "AST-SST"

        elif "HPRR" in basepf:
            tpnm = "HPRR"
        elif "OER" in basepf:
            tpnm = "OER"
        elif "EIS" in basepf:
            if "CV" in basepf:
                tpnm = "EIS+CV"
            elif "HER" in basepf:
                tpnm = "EIS_HER"
            else:
                tpnm = "EIS"
        elif "HER" in basepf:
            if "shA" in basepf and "AST" in basepf:
                tpnm = "AST-HER"
            else:
                tpnm = "HER"
        elif "short-cycling" in basepf:
            tpnm = "N2_short_cycling"
        elif "_AST_" in basepf and not "_post" in basepf:
            tpnm = "AST-SST"
        else:
            tpnm = "N2_act"
    elif "O2" in basepf:
        gas = "O2"
        if "ORR" in basepf and not "EIS" in basepf:
            tpnm = "ORR"
        elif "CV" in basepf:
            tpnm = "ORR-CV"
        elif "EIS" in basepf and not "ORR" in basepf:
            tpnm = "EIS"
        elif "LC" in basepf and not "ORR" in basepf:
            tpnm = "AST-LC"
        elif "SST" in basepf and not "ORR" in basepf:
            tpnm = "AST-SST"
        else:
            tpnm = "ORR"
    else:
        if ("OR" or "ORR") in basepf:
            gas, tpnm = "O2", "ORR"

        elif "HER" in basepf and not "EIS" in basepf:
            gas, tpnm = "N2", "HER"
        elif "OER" in basepf:
            gas, tpnm = "N2", "OER"
        elif "sHA" in basepf:
            gas, tpnm = "N2", "AST-HER"

        elif "HPRR" in basepf:
            gas, tpnm = "N2", "HPRR"
        elif "O2_" in basepf and not "H2O2" in basepf:
            gas, tpnm = "O2", "ORR_"
        elif "EIS" in basepf and not any(g in basepf for g in ["N2", "O2"]):
            gas, tpnm = "N2_guess", "EIS"
        elif "CV" in basepf and not any(g in basepf for g in ["N2", "O2"]):
            gas, tpnm = "N2_guess", "CV"
        else:
            gas, tpnm = None, None
    # if "AST SST" in comment_act0_DF:
        # gas, tpnm = "N2", "AST-SST"
    return {"Gas": gas, "PAR_exp": tpnm}


def determine_postAST_from_filename(filenamesplit) -> Dict:
    ''' determine postAST from filenamesplit'''

    basepf_split, postAST = filenamesplit, ""
    if any(
        (s in basepf_split for s in ["postAST", "pAST"])
        and (["post" and "AST" in basepf_split])
    ):
        if any(s in basepf_split for s in ["LC"]):
            postAST = "postAST_LC"
        elif any(s in basepf_split for s in ["pAST"]):
            postAST = "postAST_sHA"
        else:
            postAST = "postAST"
    elif any(s in basepf_split for s in ["postAST-LC"]):
        postAST = "postAST_LC"
    elif any(s in basepf_split for s in ["pAST-sHA"]):
        postAST = "postAST_sHA"
    elif any("AFTER" in i.upper() for i in filenamesplit) or any(s in basepf_split for s in ["postORR"]):
        postAST = "postORR"
    elif any(s in basepf_split for s in ["postAST"]):
        postAST = "postAST"
    else:
        postAST = "no"
    return {"postAST": postAST}

def determine_date_from_filename(
    filenamesplit: List, verbose=False, testing_formats=["%Y-%m-%d", "%d.%m.%Y"]
) -> Union[None, datetime.datetime]:
    '''determine date from splitted filename'''
    Dt_options = []
    _verbose = []
    for i in filenamesplit:
        for f in testing_formats:
            try:
                Dt_options.append(datetime.datetime.strptime(i, f))
            except ValueError as ve:
                _verbose.append((i, f, ve))
            except Exception as e:
                print(f"Datetime file {file}, other error {e} for {i}")

    if Dt_options:
        if len(Dt_options) == 1:
            outDt = Dt_options[0]
        else:
            logger.warning(f"Datetime file {'_'.join(filenamesplit)} multiple options {Dt_options}")
            outDt = Dt_options[0]
    else:
        logger.warning("Datetime file {'_'.join(filenamesplit)} has no formatted date options {Dt_options}")
        outDt = None
    return {'Date' : outDt}


def determine_ink_loading_from_filename(filename, reference_date=None) -> Dict:
    ''' guesses the ink loading from name'''
    if not reference_date:
        try:
            reference_date = determine_date_from_filename(filename)
        except Exception as exc:
            reference_date = datetime.datetime(2000, 1, 1)

    loading_variation_ref = loading_ref(reference_date)
    #        basepf_split = acid.query('Loading_name == "loading-unknown"').basename.iloc[0].split('_')
    guesses = ["load", "loading", "low-load", "high-load", "normal-load", "half-load"]
    matches = []
    if any(s in basepf for s in guesses):
        matches = [s for s in guesses if s in basepf]
        if any("low" in match for match in matches):
            loading_name = "low"
        #                    print('low', basepf)
        elif any("half" in match for match in matches):
            loading_name = "half"
        elif any("high" in match for match in matches):
            loading_name = "high"
        elif any("normal" in match for match in matches):
            loading_name = "normal"
        else:
            loading_name = "loading-ref-unknown"
        loading_cm2 = loading_variation_ref[loading_name]
    else:
        loading_name = "standard"
        loading_cm2 = loading_variation_ref[loading_name]
    return {
        "Loading_name": loading_name,
        "Loading_cm2": loading_cm2,
        "Loading_date": reference_date,
    }


def determine_electrode_name_from_stem(stem: str):
    """electrode name"""
    # add V3F, etc...
    if "Pt-ring" in stem or "ring" in stem or "Ring" in stem:
        electrode = "Pt_ring"
    elif "Disk" in stem or "disk" in stem:
        electrode = "disk"
    else:
        electrode = "Unknown"
    return {'electrode_type': electrode}


def _depr_filestem_to_sid_and_pos(stem: str, seps=("_", " ", "-")) -> Tuple[str, str]:
    """
    Parser for the filenames -> finds SampleID and sample position

    Parameters
    ----------
    # ramanfile_stem : str
    #    The filepath which the is parsed
    seps : tuple of str default
        ordered collection of seperators tried for split
        default : ('_', ' ', '-')

    Returns
    -------
    tuple of strings
        Collection of strings which contains the parsed elements.
    """

    split = None
    first_sep_match_index = min(
        [n for n, i in enumerate(seps) if i in stem], default=None
    )
    first_sep_match = (
        seps[first_sep_match_index] if first_sep_match_index is not None else None
    )
    split = stem.split(first_sep_match)
    _lensplit = len(split)

    if _lensplit == 0:
        sID, position = split[0], 0
    elif len(split) == 1:
        sID, position = split[0], 0
    elif len(split) == 2:
        sID = split[0]
        _pos_strnum = "".join(i for i in split[1] if i.isnumeric())
        if _pos_strnum:
            position = int(_pos_strnum)
        else:
            position = split[1]
    elif len(split) >= 3:
        sID = "_".join(split[0:-1])
        position = int("".join(filter(str.isdigit, split[-1])))
    #                split =[split_Nr0] + [position]
    return (sID, position)


def sID_to_sgrpID(sID: str, max_len=4) -> str:
    """adding the extra sample Group key from sample ID"""

    _len = len(sID)
    _maxalphakey = min(
        [n for n, i in enumerate(sID) if not str(i).isalpha()], default=_len
    )
    _maxkey = min((_len, _maxalphakey, max_len))
    sgrpID = "".join([i for i in sID[0:_maxkey] if i.isalpha()])
    return sgrpID


def format_filename(s):
    """Take a string and return a valid filename constructed from the string.
    Uses a whitelist approach: any characters not present in valid_chars are
    removed. Also spaces are replaced with underscores.

    Note: this method may produce invalid filenames such as ``, `.` or `..`
    When I use this method I prepend a date string like '2009_01_15_19_46_32_'
    and append a file extension like '.txt', so I avoid the potential of using
    an invalid filename.
    """
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    filename = "".join(c for c in s if c in valid_chars)
    filename = filename.replace(" ", "_")  # I don't like spaces in filenames.
    return filename

