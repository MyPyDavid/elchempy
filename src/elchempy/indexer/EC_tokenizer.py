"""Class and methods for parsing electrochemically relevant tokens from a filename"""
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import Counter
from functools import wraps

import datetime
from typing import Tuple, List, Dict, Union

import re
import copy

import logging

logger = logging.getLogger(__name__)

from elchempy.indexer.filepath_parser import FilePathParser
from elchempy.indexer.extra_EC_info import loading_ref, WE_surface_area_cm2

### for Developing
# from elchempy.config import LOCAL_FILES

### 3rd Party imports

import datefinder

import dateutil
from dateutil.parser import ParserError


#%%
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


def _dev_fname():
    fname = "O2_ORR_JOS12_223_#1_Disc_Versastat"
    fname = "N2_HPRR_JOS14_0.1MH2SO4+10mMH2O2_disk_221 2. Versuch"
    tokenize_name_into_remainder(fname)


if 0:
    date_patterns = [
        "^([1-9]|0[1-9]|1[0-9]|2[0-9]|3[0-1])(\.|-|/)([1-9]|0[1-9]|1[0-2])(\.|-|/)([0-9][0-9]|19[0-9][0-9]|20[0-9][0-9])$|^([0-9][0-9]|19[0-9][0-9]|20[0-9][0-9])(\.|-|/)([1-9]|0[1-9]|1[0-2])(\.|-|/)([1-9]|0[1-9]|1[0-9]|2[0-9]|3[0-1])$"
    ]
    date_res = search_pattern_and_cutout(name, patterns=date_patterns)
    date = determine_date_from_filename(cutsplit)


def get_most_common_split(name, name_separators=["_", "-"]):
    sepcounter = Counter([i for i in name if i in name_separators])
    # split_info.update(**{ name: { 'counter' : sepcounter}})
    if sepcounter:
        sep, _nsep = sepcounter.most_common(1)[0]
        split = name.split(sep)
    else:
        split, sep = [], None
    return split, sepcounter, sep


def tokenize_name_into_remainder(
    fname: str, name_separators=["_", "-", " "], date_from_parent=None
) -> Dict:
    """
    this method tries to find meaningful tokens, for the electrochemical experiment,
    in the filename and cuts those sequentially out of the name, leaves a remainder.
    """

    if not isinstance(fname, str):
        # in case of given path.parent
        if not isinstance(fname, Path):
            raise TypeError("Given name to tokenizer is not a string nor Path")
        fname = fname.name

    EC_info = {"fname": fname}
    sepstr = "".join(name_separators)

    split, sepcounter, sep = get_most_common_split(fname)
    # popsplit = split.copy()
    remainder_name = fname

    tok_dec_kwargs = dict(stripchars=sepstr, separator=sep)
    # remainder_name = sep.join(list(filter(bool, remainder_name.split(sep))))
    # patterns=electrode_patterns
    # sepstr
    remainder_name, EC_info = determine_electrode_name(
        remainder_name, EC_info, patterns=[], **tok_dec_kwargs
    )

    # TODO Continue with tokenizer
    # testing_formats = ["%Y-%m-%d", "%d.%m.%Y"]
    # day_date_finder = datefinder.DateFinder(first='day')
    # year_date_finder = datefinder.DateFinder(first='year')
    date_matches = list(datefinder.find_dates(remainder_name, source=True))
    if not date_matches:
        date_matches = []
        digitlensplit = [
            sp for sp in split if len("".join([i for i in sp if i.isdigit()])) >= 5
        ]
        for sp in digitlensplit:
            try:
                dt_parse = dateutil.parser.parse(sp)
                date_matches.append((dt_parse, sp))
            except ParserError:
                pass
    # date_matches = list(year_date_finder.find_dates(remainder_name, source=True))
    relevant_date_matches = list(
        filter(lambda x: 2015 < x[0].year < 2025, date_matches)
    )

    if len(relevant_date_matches) == 1:
        remainder_name = replace_and_strip(
            remainder_name, relevant_date_matches[0][1], sepstr
        )
        dt = relevant_date_matches[0][0]
        date = datetime.date(dt.year, dt.month, dt.day)
        date_info = {"date_dt": date, "date_source": relevant_date_matches[0][1]}
    elif len(relevant_date_matches) > 1:
        logger.warning(
            f"datefinder from {remainder_name}, {','.join(map(str,relevant_date_matches))}"
        )
    else:
        date_info = {"date_dt": None}
        # string_clean_end_character()
    EC_info.update(**date_info)

    remainder_name, EC_info = determine_Gas_Exp_from_filename(
        remainder_name, EC_info, **tok_dec_kwargs
    )

    remainder_name, EC_info = determine_pH_from_filename(
        remainder_name, EC_info, **tok_dec_kwargs
    )

    # @tokenizer_decorator
    remainder_name, EC_info = determine_postAST_from_filename(
        remainder_name, EC_info, **tok_dec_kwargs
    )

    if date_from_parent and not date_info.get("date_dt"):
        reference_date = date_from_parent
    elif date_from_parent and date_info.get("date_dt"):
        reference_date = date_info["date_dt"]
    else:
        reference_date = None

    remainder_name, EC_info = determine_ink_loading_from_filename(
        remainder_name, EC_info, reference_date=reference_date, **tok_dec_kwargs
    )

    # remainder_name, EC_info = determine_instrument_from_filename(remainder_name, EC_info, **tok_dec_kwargs )

    instrument_patterns = {"instrument": ["bipot"]}
    cell_number_patterns = {
        "cell_number": ["((C|c)ell[1-9]{1})|(internal-|)dummy-cell"]
    }
    rpm_patterns = {
        "RPM": [
            "([1-9]{1}(r|)pm(s|))|[0-9]{3,4}(rpm|rmp)|(0rpm)|(rpm-range-[0-9]{1,3}mV)"
        ]
    }
    template_patterns = {
        "template_file": [
            "template|TEMPLATE|([A-Z]{2,3}[1]{0,1}[x]{1,3})|((_|)[x]{2,3}_)|(([A-Z]{2,3}[1-9]{1,2}|)templ([A-Z]{2,3}[1-9]{1,2}|))",
            "settings|SETTINGS",
        ]
    }
    n2_cycles_patterns = {
        "N2_cycles": [
            "[0-9]{1,2}(cls|cl)( |_|-)([0-9]{2,3}(|( |_)[0-9]{2,3}(|( |_)[0-9]{2,3}|)))"
        ]
    }
    AST_cycles_patterns = {
        "AST_cycles": ["((20k\+CVs)(\+EIS|))|(20\.000cycles)|([0-9]{4,5}c)"]
    }
    trial_number_pattern = {"trial_number": ["\A([1-9]{1})(.|)", "test[1-9]{1}"]}
    fail_pattern = {"failure": ["(fail|fail_dup|noisy|wrong|badscan)"]}
    ink_status_pattern = {
        "ink_status": [
            "ink((-| )([0-3]{1}[0-9]{1}\.[0-1]{1}[1-9]{1}|)|)",
            "(new|old|bad)(_|-)ink",
            "[1-9]{1,2}(m|)yl",
            "[1-9]{1,2}ul",
        ]
    }
    file_status = {"file_status": ["((( \- |)(Kopie|Copy)){1,3}|Conflict|duplicate)"]}
    rhe_potential_mV_patterns = {"RHE_potential_mV": ["\A([0-9]{3})|([0-9]{3}\Z)"]}
    rhe_potential_mV_extra_patterns = {
        "RHE_potential_mV_extra": ["(| |-|_)([0-9]{3})(| |-|_)"]
    }
    persons_name = {"person_name": ["Nata|Sandeep|Nicole"]}

    all_extra_patterns = {
        **instrument_patterns,
        **cell_number_patterns,
        **rpm_patterns,
        **template_patterns,
        **n2_cycles_patterns,
        **AST_cycles_patterns,
        **fail_pattern,
        **ink_status_pattern,
        **file_status,
        **rhe_potential_mV_patterns,
        **persons_name,
    }
    # **trial_number_pattern,

    for pattname, patterns in all_extra_patterns.items():

        remainder_name, EC_info = tokenize_name_from_patterns(
            remainder_name,
            EC_info,
            patterns=patterns,
            token_name=pattname,
            **tok_dec_kwargs,
        )

    # template_file_found = EC_info.get("template_file_0", '')
    remainder_name, EC_info = match_SampleID(
        remainder_name,
        EC_info,
        template_file_found=EC_info.get("template_file_0", ""),
        **tok_dec_kwargs,
    )
    if remainder_name:
        # one more time try to find the rhe_potential pattern in remainder
        # rhe_pot_name = list(rhe_potential_mV_patterns.keys())[0]
        # print(remainder_name)
        remainder_name, EC_info = tokenize_name_from_patterns(
            remainder_name,
            EC_info,
            patterns=list(rhe_potential_mV_extra_patterns.values())[0],
            token_name=list(rhe_potential_mV_patterns.keys())[0],
            **tok_dec_kwargs,
        )
    if remainder_name:
        remainder_name, EC_info = tokenize_name_from_patterns(
            remainder_name,
            EC_info,
            patterns=list(trial_number_pattern.values())[0],
            token_name=list(trial_number_pattern.keys())[0],
            **tok_dec_kwargs,
        )

    EC_info.update(**{"token_remainder": remainder_name})

    return EC_info


#%%
def _test_pattern(name, patterns: dict):
    for pn, pl in patterns.items():
        for p in pl:
            m = re.search(p, name)
            print(f"{pn} {p}: {m}")


class TokenizerError(Exception):
    """raises errors from Tokenizer"""


def tokenizer_decorator(func, **kwargs):
    """
    This decorator wraps around a tokenizer function.
    It adds the token to the info dict and removes the found token from the given name.
    """
    if not callable(func):
        raise TypeError(f"func {func} not callable")

    @wraps(func)
    def wrapper(name, info, **kwargs):

        try:
            if ("patterns" and "token_name") in kwargs:
                token = func(name, **kwargs)
            elif "reference_date" in kwargs:
                token = func(name, reference_date=kwargs.get("reference_date", None))
            elif "template_file_found" in kwargs:
                token = func(
                    name, template_file_found=kwargs.get("template_file_found", None)
                )

            else:
                token = func(name)
        except TypeError as ex:
            logger.error(f"func: {func.__name__}, name: {name}\n{kwargs}")
            raise TokenizerError(ex) from ex
            # return name, info
        except Exception as ex:
            logger.error(f"func: {func.__name__}, name: {name}\n{kwargs}")
            raise TokenizerError(ex) from ex
            # return name, info
        if not token:
            # logger.warning(f'Wrapper no token found for {func}, {name}')
            return name, info
        str_token_values = [i for i in token.values() if isinstance(i, str)]
        str_token_values_in_name = [i for i in str_token_values if i in name]
        if str_token_values:
            for val in str_token_values_in_name:
                val_is_subset = [
                    i
                    for i in str_token_values_in_name
                    if val in i and len(i) > len(val)
                ]
                if not val_is_subset:
                    name = replace_and_strip(name, val, **kwargs)
        info.update(**token)
        # print("wrapper token:",info,'\nname',name)
        return name, info

    return wrapper


@tokenizer_decorator
def tokenize_name_from_patterns(name, patterns=None, token_name="", **kwargs):
    """tokenize a string from certain patterns"""

    pattern_res = search_pattern_and_cutout(name, patterns=patterns)
    info_res = {}
    for n, res in enumerate([i for i in pattern_res if i[-1]]):
        # name = replace_and_strip(name, res[-1], stripchars)
        info_res = {**info_res, **{f"{token_name}_{n}": res[-1]}}
    logger.debug(
        f"Info found for patterns{patterns} on {name}, {token_name}, {info_res}"
    )
    return info_res


def search_pattern_and_cutout(stem: str, patterns=None, **kwargs) -> List:
    """electrode name"""
    # add V3F, etc...
    # search_and_cut_out(stem, '', flags=re.IGNORECASE)
    # new_stem = copy.deepcopy(stem)
    res = []
    if not patterns:
        return res

    for n, pattern in enumerate(patterns):
        try:
            new_stem, match = search_and_cut_out(stem, pattern, **kwargs)
        except Exception as ex:
            logger.error(f"Pattern Error {ex}, stem {stem}, pattern {pattern}")
            raise ex from ex
        res.append((n, pattern, stem, new_stem, match))
    return res


def search_and_cut_out(stem: str, pattern: str, **kwargs):
    # pattern='RRDE[0-9]{5}'
    search = re.search(pattern, stem, **kwargs)
    if not search:
        return stem, search
    match = search.group()
    st, end = search.span()
    cutout_stem = stem[0:st] + stem[end:-1]
    cutout_stem = string_clean_end_character(cutout_stem)
    return cutout_stem, match


def replace_and_strip(
    name: str,
    value: str,
    stripchars: str = "_-",
    replace_with: str = "",
    separator=None,
    **kwargs,
):
    if not name or not value:
        return name

    if not isinstance(value, str):
        return name.strip(stripchars)
    logger.debug(f"replace and strip: {name} with {value}")
    if value in name:
        name = name.replace(value, replace_with)
        if separator:
            name = separator.join(list(filter(bool, name.split(separator))))
    logger.debug(f"replace and strip new name: {name}")
    return name.strip(stripchars)


def determine_instrument_from_filename(
    name: str, info: dict, **kwargs
) -> Tuple[str, dict]:
    """instrument from filename"""
    bipot = re.search("bipot|V3F", name)
    if bipot:
        source = bipot.group(0)
        name = name[0 : bipot.span()[0]] + name[bipot.span()[1] : :]
    else:
        source = None
    # v3f = re.search('V3F',remainder_name)
    info.update(**{"Instrument": source})
    return name, info


@tokenizer_decorator
def determine_electrode_name(name, **kwargs):
    """electrode names"""
    electrode_patterns = [
        "(Pt|pt)[-]{0,1}[\W]ring",
        "(PDX_Ch1)((_|)|dis(c|k))",
        "(PDX_Ch2)((_|)|ring)",
        "(PDX P[1-9]{1}K)(-|)(Ch[1-3]{1})",
        "(#[1-3]{1}_|)((D|d)is(c|k)|(R|r)ing)",
        "(#[1-3]{1}_|)Dis(c|k)(_Versastat|_Parstat)",
        "(#[1-3]{1}_|)Ring(_Versastat|_Parstat)",
        "(#[1-3]{1}_|)Ch[1-9]{1}( |_)(Disk|Ring)( |_)( |PDX)",
        "V3(F|)(_|)((D|d)is(c|k)|(R|r)ing|)",
        "AgAgCl[1-9]{0,1}",
        "Kalomel",
        "(Hg|)HgO[1-9]{0,1}",
        "Hg[1-9]{1}",
        "RRDE[1-4]{0,1}([0-9]{4}|)",
        "RD(E|e)[1-2]{0,1}",
        "RHE[1-9]{0,1}",
        "Pt-wire",
        "N-doped-graphene",
        "(Fe|Ni|)(-|)foam(-N-doped|-graphene|)",
        "RHE-vs-RE[1-9]{0,2}",
    ]

    elec_res = search_pattern_and_cutout(name, patterns=electrode_patterns)
    info_res = {}
    for n, res in enumerate([i for i in elec_res if i[-1]]):
        # name = replace_and_strip(name, res[-1], stripchars)
        info_res = {**info_res, **{f"electrode_{n}": res[-1]}}
    return info_res


def match_sampleid_Pt_ring():
    Pt_ring_match = [
        re.search("(?<=(pt)){0,1}.{0,1}ring{0,1}", i, re.IGNORECASE) for i in FileParts
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
        [i for i in FileParts if re.search("(Pt-ring|pt-ring)", i, re.IGNORECASE)]
    ):
        formatID_matches = ["Pt-ring"]
        sampleID = "Pt_ring"


@tokenizer_decorator
def match_SampleID(
    file, message=False, include_extra=False, template_file_found=None
) -> Dict:
    """finds a possible match in the filename with a known sample id"""
    #        message = True
    if not file:
        return None

    formatID_matches, sampleID = [], None
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
        i
        for i in FileParts_raw
        if i.upper() not in skip_parts and re.search("[a-zA-Z]", i)
    ]

    while not sampleID:

        if template_file_found:
            sampleID = template_file_found

        if "DW" in FileParts:
            formatID_matches = [
                i
                for i in FileParts
                if re.search("(DW[0-9]{1,3}[a-zA-Z]{0,1})", i, re.IGNORECASE)
            ]

        elif (
            any([i for i in FileParts if re.match("(DWIM)", i, re.IGNORECASE)])
            and not "Pt-ring" in formatID_matches
        ):
            formatID_matches = [i for i in FileParts if re.match("(DWIM)", i)]

        elif any([i for i in FileParts if re.search("(SD)", i)]):
            formatID_matches = [
                i for i in FileParts if re.search("(SD)", i, re.IGNORECASE)
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
                if re.search("(?i)(xx|xxx|DWxx|JOS1xx)", i, re.IGNORECASE)
            ]
        ):
            formatID_matches = ["template"]
        elif any(
            [
                i
                for i in FileParts
                if re.search("([A-Z]{2,3}xx)|(templ|template)\Z", i, re.IGNORECASE)
            ]
        ):
            formatID_matches = ["template"]

        elif "NW" in FileParts:
            formatID_matches = [
                i for i in FileParts if re.search(".?NW_M.?", i, re.IGNORECASE)
            ]
        elif any("JOS" in s for s in FileParts):
            if any([i for i in FileParts if re.search("(ring)", i)]):
                formatID_matches = ["Pt-ring"]
            elif any([i for i in FileParts if re.search("(?i)JOS[0-9]{0,3}XX", i)]):
                formatID_matches = ["template"]
            else:
                formatID_matches = [
                    i.upper()
                    for i in FileParts
                    if re.search(".?JOS.?", i, re.IGNORECASE)
                ]
        #            elif any('jos' in s for s in FileParts):
        #                formatID_matches = [i for i in FileParts if re.match('.?jos.?',i)]
        elif any("PT" in s for s in FileParts):
            formatID_matches = [
                i
                for i in FileParts
                if re.search("PT(A|B)([0-9]{1,3})(,\d|)", i, re.IGNORECASE)
            ]

        elif any(
            [
                i
                for i in FileParts
                if re.search("([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i, re.IGNORECASE)
            ]
        ):
            formatID_matches = [
                i
                for i in FileParts
                if re.search("([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i, re.IGNORECASE)
            ]

        elif any("CELL" in s for s in FileParts):
            formatID_matches = [
                i
                for i in FileParts
                if re.search(".?CELL([0-9]{0,3}).?", i, re.IGNORECASE)
            ]

        else:
            formatID_matches = [
                i
                for i in FileParts
                if re.search("([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i, re.IGNORECASE)
            ]
        #            print(FileParts,formatID_matches)
        try:
            fmtID_split = formatID_matches[0].split("-")
            if len(fmtID_split) > 1:
                fmtID_split_matches = [
                    i
                    for i in fmtID_split
                    if re.search("(?i)([A-Z]{2}[0-9]{2,3}[a-z]{0,1})", i)
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
        extra_parts = [i for i in FileParts_raw if i.upper() not in sampleID.upper()]
        extra_part_matches = [
            i
            for i in extra_parts
            if re.search("([0-9]|HT|AL|Ir|N2|repeat|[a-zA-Z]{0,1})", i, re.IGNORECASE)
        ]
        if any(extra_part_matches):
            sampleID = "_".join([sampleID, *extra_part_matches])

    if message:
        logger.warning('SampleID: "%s" used for %s' % (sampleID, file))
    return {"SampleID": str(sampleID)}


@tokenizer_decorator
def determine_pH_from_filename(filename: str) -> Dict:
    """determine pH from filename"""

    PAR_file_test = filename
    pH = {}
    # while not pH: # FIXME remove while loop ... possible infinite and no need
    if re.search("H2SO4", PAR_file_test) and not re.search("HClO4", PAR_file_test):
        if re.search("MeOH", PAR_file_test):
            if re.search("1M.?MeOH", PAR_file_test):
                pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+1M-MeOH"}
            elif re.search("0\.1.?MeOH", PAR_file_test):
                pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+0.1M-MeOH"}
            else:
                pH = {"pH": 97, "Electrolyte": "?MH2SO4+MeOH"}
        elif re.search("H2O2", PAR_file_test):
            pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+10mMH2O2"}
        else:
            if re.search("(?<=(0\.1))[a-zA-Z]{0,2}(|-)(H2SO4)", PAR_file_test):
                #                        pH = {'pH' : 1, 'Electrolyte' : '0.1MH2SO4'}
                if re.search(
                    "(?<=(0\.1))[a-zA-Z]{0,2}(H2SO4)(?!(.{0,5}(MeOH|H2O2)))",
                    PAR_file_test,
                ):
                    pH = {"pH": 1, "Electrolyte": "0.1MH2SO4"}
                elif re.search(
                    "(?<=(0\.1))[a-zA-Z]{0,2}-(H2SO4)(?!(.{0,5}(MeOH|H2O2)))",
                    PAR_file_test,
                ):
                    pH = {"pH": 1, "Electrolyte": "0.1M-H2SO4"}

                elif re.search(
                    "(?<=(0\.1))[a-zA-Z]{0,2}(H2SO4)(?=(.{0,5}(H2O2)))",
                    PAR_file_test,
                ):
                    pH = {"pH": 1, "Electrolyte": "0.1MH2SO4+?H2O2"}
                else:
                    pH = {"pH": 1, "Electrolyte": "?MH2SO4"}

            elif re.search("(?<=(0\.5))[a-zA-Z]{0,2}(H2SO4)", PAR_file_test):
                pH = {"pH": 0.3, "Electrolyte": "0.5MH2SO4"}
            else:
                pH = {"pH": 98, "Electrolyte": "?MH2SO4"}
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
                pH = {"pH": 1, "Electrolyte": "0.1MHClO4+?MH2O2"}

        elif re.search("(?<=(0\.1.))HClO4(?!(.{0,5}H2O2))", PAR_file_test):
            pH = {"pH": 1, "Electrolyte": "0.1MHClO4"}
        else:
            pH = {"pH": 1, "Electrolyte": "?MHClO4"}
    elif re.search("KOH|koh", PAR_file_test):
        if re.search("(?<=(0\.1))[a-zA-Z]{0,2}(KOH)(?!(.{0,5}H2O2))", PAR_file_test):
            pH = {"pH": 13, "Electrolyte": "0.1MKOH"}
        else:
            if "KOH" in PAR_file_test:
                token = "KOH"
            elif "koh" in PAR_file_test:
                token = "koh"
            else:
                token = "?koh?"
            pH = {"pH": 13, "Electrolyte": "?MKOH", "pH_token": token}
    #                pH = {'pH' : 13, 'Electrolyte' : 'xxKOH'}
    elif re.search("NaOH", PAR_file_test):
        if re.search("(?<=(0\.1))[a-zA-Z]{0,2}(NaOH)(?!(.{0,5}H2O2))", PAR_file_test):
            pH = {"pH": 13, "Electrolyte": "0.1MNaOH"}
        else:
            pH = {"pH": 13, "Electrolyte": "?MNaOH", "pH_token": "NaOH"}
    elif re.search("Acid", PAR_file_test):
        pH = {"pH": 1, "Electrolyte": "0.1MH2SO4_acid", "pH_token": "Acid"}
    else:
        pH = {"pH": None, "Electrolyte": None}
        # "Other_{PAR_file_test}"}
    return pH


@tokenizer_decorator
def determine_Gas_Exp_from_filename(filename) -> Dict:
    """returns a dict with Gas and Exp type from the filename"""
    basepf = filename

    if "OCP" in basepf:
        # or "RHE"
        return {"Gas": "N2", "PAR_exp": "OCP"}

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
            tpnm = "AST"

        elif "HPRR" in basepf:
            tpnm = "HPRR"
        elif "OER" in basepf:
            tpnm = "OER"
        elif "EIS" in basepf:
            if "CV" in basepf:
                tpnm = "EIS+CV"
            elif "EIS-range_HER" in basepf:
                tpnm = "EIS-range_HER"
            elif "EIS-range" in basepf and not "HER" in basepf:
                tpnm = "EIS-range"
            else:
                tpnm = "EIS"
        elif "HER" in basepf:
            if "shA" in basepf and "AST" in basepf and not "pAST" in basepf:
                if "sHA_short_HER-AST" in basepf:
                    tpnm = "sHA_short_HER-AST"
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
            if "ORR-act-sel" in basepf:
                tpnm = "ORR-act-sel"
            elif "ORR_act-sel" in basepf:
                tpnm = "ORR_act-sel"
            else:
                tpnm = "ORR"
        elif "CV" in basepf:
            tpnm = "ORR-CV"
        elif "EIS" in basepf and not "ORR" in basepf:
            if "EIS-range" in basepf:
                tpnm = "EIS-range"
            else:
                tpnm = "EIS"

        elif "LC" in basepf and not "ORR" in basepf:
            tpnm = "LC"
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
        # elif basepf.startswith("AST_"):
        # elif 'AST_20k+CVs+EIS' in basepf:
        #     tpnm = 'AST_20k+CVs+EIS'
        # elif 'AST_20k+CVs' in basepf:
        #     tpnm = 'AST_20k+CVs'
        else:
            gas, tpnm = None, None

    # if "AST SST" in comment_act0_DF:
    # gas, tpnm = "N2", "AST-SST"
    return {"Gas": gas, "PAR_exp": tpnm}


@tokenizer_decorator
def determine_postAST_from_filename(filenamesplit) -> Dict:
    """determine postAST from filenamesplit"""

    basepf_split, postAST = filenamesplit, ""
    if all([i in basepf_split for i in ("post", "AST")]):
        if any(s in basepf_split for s in ["LC"]):
            if "postAST_LC" in basepf_split:
                postAST = "postAST_LC"
            elif "postAST-LC" in basepf_split:
                postAST = "postAST-LC"
            else:
                postAST = "postAST?LC"
        elif "postAST_sHA" in basepf_split:
            postAST = "postAST_sHA"
        elif "post_AST" in basepf_split:
            postAST = "post_AST"
        elif "postAST" in basepf_split:
            postAST = "postAST"
        else:
            postAST = "postAST?"
    elif "pAST-sHA" in basepf_split:
        postAST = "pAST-sHA"
    elif any("AFTER" in i.upper() for i in filenamesplit) or any(
        s in basepf_split for s in ["postORR"]
    ):
        postAST = "postORR"
    elif "POST-EXP" in basepf_split:
        postAST = "POST-EXP"
    else:
        postAST = None
    return {"postAST": postAST}


def _determine_date_from_filename(
    filenamesplit: List, verbose=False, testing_formats=["%Y-%m-%d", "%d.%m.%Y"]
) -> Union[None, datetime.datetime]:
    """determine date from splitted filename"""
    date_options = []
    _verbose = []
    for i in filenamesplit:
        for f in testing_formats:
            try:
                date_options.append(datetime.datetime.strptime(i, f))
            except ValueError as ve:
                _verbose.append((i, f, ve))
            except Exception as e:
                logger.warning(f"Datetime file {file}, other error {e} for {i}")
    if date_options:
        if len(date_options) == 1:
            outDt = date_options[0]
        else:
            logger.warning(
                f"Datetime file {'_'.join(filenamesplit)} multiple options {date_options}"
            )
            outDt = date_options[0]
    else:
        logger.warning(
            f"Datetime file {'_'.join(filenamesplit)} has no formatted date options {date_options}"
        )
        outDt = None
    return {"Date": outDt}


@tokenizer_decorator
def determine_ink_loading_from_filename(filename, reference_date=None) -> Dict:
    """guesses the ink loading from name"""
    # print("ref date in ", reference_date)
    # if not reference_date:
    # reference_date = kwargs.get('date_dt', None)

    if not reference_date:
        reference_date = datetime.date(2000, 1, 1)
        # try:
        #     reference_date = determine_date_from_filename(filename).get(
        #         "Date", _default
        #     )
        # except Exception as exc:
        #     reference_date = _default

    loading_variation_ref = loading_ref(PAR_date=reference_date)
    guesses = ["load", "loading", "low-load", "high-load", "normal-load", "half-load"]
    matches = []
    match = None
    if any(s in filename for s in guesses):
        matches = [s for s in guesses if s in filename]
        if matches:
            match = max(matches, key=lambda x: len(x))

        if any("low" in match for match in matches):
            loading_name = "low"
        #                    print('low', filename)
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
        "Loading_mgcm2": round(loading_cm2, 3),
        "Loading_date": reference_date,
        "Loading_source": match,
    }


def string_clean_end_character(string):
    endsw = re.search("[^a-zA-Z0-9]\Z", string)
    if endsw:
        string = string[0 : endsw.start()]
    startsw = re.match("\A[^a-zA-Z0-9]", string)
    if startsw:
        string = string[startsw.end() : :]
    return string


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


# @tokenizer_decorator
def _depr_try_find_sampleID(file, template_file_found=None):
    """---unused-----
    #Input string or Path and output is SampleID[0] and dict
    """
    pf = Path(file)
    try:
        pf_Match = match_SampleID(pf)
    except:
        pf_Match = "NoneNoneNone"
    if pf.is_file():
        parent_folder = pf.parent.parts[-1]

        try:
            pf_StemMatch = match_SampleID(pf.stem)
        except:
            pf_StemMatch = "None"
        try:
            pf_ParentMatch = match_SampleID(parent_folder)
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
        pf_StemMatch, pf_ParentMatch = match_SampleID(pf.stem), "None"
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
