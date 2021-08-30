"""
Created on Thu Jul 15 10:25:51 2021

@author: DW
"""

from pathlib import Path
from html.parser import HTMLParser


class ParserError(ValueError):
    """unable to parse this file"""


def read_PAR_file(filepath):
    """
    Special parser for Versatstudio ".par" files

    opens the file, cleans the lines
    initializes the Parser and feeds the data
    returns the Parser Object

    """
    try:
        with open(filepath) as fp:
            fp_readlines = fp.readlines()
        fp_read = find_replace_line_endings(fp_readlines)
    except OSError:
        raise ParserError(
            "Can not open or read this file" "File: {filepath} is invalid."
        )

    VSP = VersaStudioParser()
    VSP.feed(fp_read)
    VSP.close()
    # VSP._interesting_data
    return VSP


def find_replace_line_endings(fp_readlines):
    """
    special find and replace function
    to clean up line endings in the file
    from end or starttag characters
    """
    clean = []
    for line in fp_readlines:
        if line.endswith("=<\n"):
            line = line.replace("<\n", "lt\n")
        clean.append(line)
    return "".join(clean)


class VersaStudioParser(HTMLParser):
    """
    Main VersaStudio .par file parser.
    It seperates the read-in data already
    in an actions part and data part,
    following the structure of the .par file.

    Usage is similar to HTMLParser:
        VSP = VersaStudioParser()
        VSP.feed(data)
        VSP.close

    Data can be found in attributes of VSP:
        actions
        data_body
        data_keys

    """

    _VSP_VERSION = "0.1.0"

    _skipped_tags = ("dockinglayout", "dockpanel", "graph1")

    _skipped_data = (">", "<", "\n\n")
    # _meta_tags = ('application', 'instrument','experiment')
    _data_name = "segment"
    _action_name = "action"

    actions = {}
    data_body = {}
    data_keys = []
    data_version = {}
    # wronglines = {}
    _tags = []
    _all_raw_data = {}

    def handle_starttag(self, tag, attrs):
        self.tag = ""
        # self.endtag = False
        if tag not in self._skipped_tags:
            # and not any(tag.startswith(i) for i in self._skipped_tags_startswith):
            # if tag in self.meta_tags:
            # self._tags.append(tag)
            self.tag = tag
        else:
            pass
            # print("skipped start tag :", tag)
        # if not hasattr(self, tag):
        # setattr(self, tag, tag)

    def handle_endtag(self, tag):
        if self.tag:
            pass

    def handle_data(self, data, max_len=None):
        """
        handles data depending on tag name
        different way of handling 'action' blocks
        and the data in a 'segment' block
        """

        if max_len:
            if len(data) > max_len:
                data = data[0:max_len]

        if data not in self._skipped_data:

            self._all_raw_data.update({self.tag: data})

            if self.tag.startswith(self._data_name):
                # parse main data
                self.data_version.update(self.parse_text(data, self.tag))

                for segkey, segval in self.data_version.items():
                    data_definition = ""
                    if isinstance(segval, dict):
                        data_definition = segval.get("Definition", "")

                    if data_definition:
                        data_keys = [i for i in data_definition.split(", ") if i != "0"]
                        self.data_keys = data_keys
                        self.data_body.update(
                            self.parse_data_body(data, segkey, data_keys)
                        )

            elif self.tag.startswith(self._action_name):
                # print("starting data handling:", self.tag)
                # print('datahandler:', self.tag)
                # parse actions and other metadata
                pass
                self.actions.update(self.parse_text(data, self.tag))
            else:
                pass
            # self._all_data.update({self.tag : data.strip()})

    def parse_data_body(self, text, segment_name: str, data_keys: list):

        data_body = []
        # data_body.append(data_keys)
        lenkeys = len(data_keys)

        for line in text.splitlines():
            splt = line.split(",")
            if len(splt) == lenkeys:
                try:
                    splt = [int(i) if not "." in i else float(i) for i in splt]
                except Exception as e:
                    pass

                data_body.append(splt)
        return {segment_name: data_body}

    def parse_text(self, text, current_tag):
        text = text.strip()
        textdict = {current_tag: None}
        try:
            # splt = [i.split("=") for i in text[1:-1].split("\n")]
            splitted_lines = [
                line.split(sep="=")
                for line in text.splitlines()
                if (line != "" and "=" in line)
            ]
            lines_notlen2 = [i for i in splitted_lines if len(i) != 2]

            if lines_notlen2:
                # FIXME report errors
                # print(f'parsemeta error for {current_tag} missing len2 lines {len(lines_notlen2)}')
                splitted_lines = [i for i in splitted_lines if len(i) != 2]

            # Cast values to numeric types
            splitted_lines_dict = cast_elements_to_numeric(splitted_lines)
            if splitted_lines_dict:
                textdict.update({current_tag: splitted_lines_dict})
            else:
                raise ValueError
        except Exception as e:
            pass
            # print(f'parsemeta error for {current_tag} {e}')
            # self.wronglines.update({current_tag: f'error {e},{text[0:100]}'})
        # print('PARSE',text,frame)
        return textdict


def cast_elements_to_numeric(
    splitted_lines: list, float_sep=(",", "."), minus_sign=("-")
) -> dict:
    """helper function for casting str to numeric in a list of lists"""
    result = {}
    for n, elem in enumerate(splitted_lines):

        if isinstance(elem, list) and len(elem) == 2:
            key, value = elem
            _isnum = "".join(
                [i for i in value if i.isnumeric() or i in (*float_sep, *minus_sign)]
            )
            if _isnum:
                # if numeric str characters are found
                if len(_isnum) == len(value):
                    # check if there is a decimal separator character in the numeric str
                    # cast in float or int
                    if any(sep in _isnum for sep in float_sep):
                        try:
                            value = float(value)
                        except Exception as e:
                            # type casting to float error
                            pass
                    else:
                        try:
                            value = int(value)
                        except Exception as e:
                            # type casting to int error
                            pass
                else:
                    # if both numeric and other characters are found in str
                    pass
            else:
                pass

        else:
            key = "error_cast_elem_{n}"
            value = elem

        result.update({key: value})
    return result
