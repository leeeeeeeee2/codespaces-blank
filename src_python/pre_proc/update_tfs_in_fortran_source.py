#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : check_transfer_funcs.py
Project Name: MPR
Description : analyzes mpr.nml file for transfer functions (TFs), compares it with source code and adds TFs to code
Author      : Stephan Thober (stephan.thober@ufz.de) and Robert Schweppe (robert.schweppe@ufz.de)
Created     : 2019-09-05 11:46
"""

# IMPORTS
import f90nml
import pathlib
from copy import copy
import re
from collections import OrderedDict
import argparse
import string
from shutil import copyfile
from src_python.pre_proc.mpr_interface import OPTIONS

# GLOBAL VARIABLES
FORTRAN_TF_SOURCEFILE = pathlib.Path('mo_mpr_transfer_func.F90')
FORTRAN_DA_SOURCEFILE = pathlib.Path('mo_mpr_data_array.F90')
MODIFIED_SUFFIX = '.mod'
BACKUP_SUFFIX = '.bak'
WORD_BOUNDARIES = string.ascii_letters + string.digits + '_'
# maximum line length for Fortran code
MAX_LINE_LENGTH = 100
EXPERIMENTAL = False
DEFAULT_CONFIG_FILE = pathlib.Path('mpr.nml')
DEFAULT_SOURCE_FOLDER = pathlib.Path('src')
FORTRAN_INDENT = '  '

TRANSLATE_DICT = OrderedDict([
    ('+', 'pl'),
    ('-', 'mi'),
    ('**', 'po'),
    ('*', 'ti'),
    ('/', 'di'),
    ('(', 'bs'),
    (')', 'be'),
    ('exp', 'ex'),
    ('log10', 'l1'),
    ('log', 'l2'),
    ('else', 'el'),
    ('if', 'if'),
    ('then', 'th'),
    ('end', 'en'),
    ('where', 'wh'),
    ('<=', 'le'),
    ('<', 'lt'),
    ('>=', 'ge'),
    ('>', 'gt'),
    ('==', 'eq'),
    ('.and.', 'ad'),
    ('.or.', 'or'),
    ('.not.', 'no'),
    ('asin', 'as'),
    ('acos', 'ac'),
    ('atan2', 'au'),
    ('atan', 'at'),
    ('sinh', 'sh'),
    ('cosh', 'ch'),
    ('tanh', 'tx'),
    ('sin', 'si'),
    ('cos', 'co'),
    ('tan', 'ta'),
    ('abs', 'ab'),
    ('max', 'ma'),
    ('min', 'mi'),
    ('sqrt', 'sq'),
    ('^', 'po')
])

# this dict is important as it stores characters that need to be ignored when inserting whitespaces around
# operators as this would break the meaning e.g. do not do "1<=2" -> "1 < =2" but "1 <= 2" for operator "<"
WORD_CHARS_DICT = {
    '<': '=',
    '>': '=',
    '*': '*',
    'log': '1',
    'sin': 'ah',
    'cos': 'ah',
    'tan': 'ah2',
    'atan': '2',
}


# FUNCTIONS
def insert_line_break_in_math_string(start_string, insert_string, end_string, *args, index_addon=0):
    trailing_whitespaces = len(start_string) - len(start_string.rstrip(' '))
    leading_whitespaces = len(end_string) - len(end_string.lstrip(' '))
    return_string = start_string.rstrip(' ') + insert_string + end_string.lstrip(' ')
    if args:
        return_string = return_string + ''.join(args)
    index_addon += len(insert_string) - trailing_whitespaces - leading_whitespaces
    return return_string, index_addon


def break_line(line, string_only=False):
    """break a long string (tf name) into multiple parts to be not longer than MAX_LINE_LENGTH"""
    parts = line.splitlines(True)
    for i_part, part in enumerate(parts):
        for break_index in range((len(part) - 1) // MAX_LINE_LENGTH, 0, -1):
            # if at position MAX_LINE_WIDTH, there is a string, split the string
            break_pos = break_index * MAX_LINE_LENGTH
            if part[:break_pos].count('\"') % 2 or string_only:
                part, _ = insert_line_break_in_math_string(
                    part[:break_pos], '"// &\n{}"'.format(FORTRAN_INDENT * 3), part[break_pos:])
            elif part[:break_pos].count('\'') % 2:
                part, _ = insert_line_break_in_math_string(
                    part[:break_pos] + "'// &\n{}'".format(FORTRAN_INDENT * 3) + part[break_pos:])
            # else insert part break at last occurrence of whitespace
            else:
                break_pos = part[:break_index * MAX_LINE_LENGTH].rfind(' ')
                part = part[:break_pos + 1] + '&\n{}'.format(FORTRAN_INDENT * 3) + part[break_pos + 1:]
        parts[i_part] = part
    return ''.join(parts)


def does_contain_pattern(pattern, search_string):
    """check whether if or where clauses as contained in raw tf string"""
    return re.search(r'\b{}\b'.format(pattern), search_string)


def get_index_in_string(long_string, part, word_chars=WORD_BOUNDARIES):
    """function from Fortran to check if a character sequence text is in a string s without being"""
    n_p = len(part)
    n_s = len(long_string)
    out_index = -1
    i = 0
    while True:
        # find index of the first character of part in string that has not been scanned so far
        if part[0] in long_string[i:n_s]:
            i_add = long_string[i:n_s].index(part[0])
        else:
            i_add = -1
        i = i + i_add
        if i_add == -1 or i + n_p > n_s:
            # the part cannot be in string as the first char is not even contained or
            # the part cannot be in string starting at i as it would be too long
            return out_index
        elif long_string[i:i + n_p] == part:
            # character matches the part
            is_begin_not_word = True
            # at beginning of string
            if i > 0:
                # is the part preceded by a alphanumeric character?
                if long_string[i - 1] in word_chars:
                    is_begin_not_word = False
                # hack so positive number is not found in negative number
                if long_string[i - 1] == '-' and part.startswith(tuple(string.digits)):
                    is_begin_not_word = False
            if is_begin_not_word:
                # is the part succeeded by a alphanumeric character?
                if i + n_p < n_s and long_string[i + n_p] in word_chars:
                    # part boundary end is violated, continue
                    i += 1
                else:
                    # index is found and part boundaries are checked
                    return i
            else:
                # part boundary start is violated, continue
                i += 1
        else:
            # word does not match, continue
            i += 1


def replace_in_string(search_string, string_pattern, replacement, *args, **kwargs):
    """
    replaces strings in a search string considering word boundaries, e.g.:
    _replace_in_string('s ss sss', 'ss', 's') -> 's s sss'"""
    index = get_index_in_string(search_string, string_pattern, *args, **kwargs)
    while index >= 0:
        search_string = search_string[:index] + replacement + search_string[index + len(string_pattern):]
        # continue the search from after the replacement
        covered_length = index+len(replacement)
        index = get_index_in_string(search_string[covered_length:], string_pattern, *args, **kwargs)
        if index >= 0:
            index += covered_length
    return search_string


# CLASSES
class TF(object):
    """class handling all string operations on a transfer function"""
    MATH_PREFIX_STRING = 'func_result(:) = '


    def __init__(self, raw_tf_string, predictors, is_test_mode=True, global_params=None, tfs=None):
        """init the TF class by providing the raw string information from the namelist
         with transfer_func, from_data_arrays and also a flag, the list of the global parameters and tfs"""
        # set args
        self.tfs = tfs or {}
        self.global_params = global_params or {}
        self.raw_tf_string = raw_tf_string
        self.predictors = predictors
        self.is_test_mode = is_test_mode

        # set main properties
        self.index_name = ''
        self.processed_tf_string = copy(self.raw_tf_string)
        self.translated_name = ''
        self.math_string = ''

    @property
    def is_contained(self):
        """is the translated tf name already contained in source code"""
        return '"{}"'.format(self.translated_name) in self.tfs

    @property
    def n_params(self):
        """get the number of parameters needed for tf"""
        return len(set(re.findall(r'p[0-9]+', self.translated_name)))

    @property
    def _max_index(self):
        """what is the next free index number for the transfer_function name"""

        return max([int(key.split('_')[2]) for key in self.tfs.values()] + [0]) + 1

    def translate_func(self):
        """
        translate the raw transfer function name (from namelist) to
        Fortran code, unique key and running index name
        """
        # initialize
        tf_string = copy(self.raw_tf_string)

        for ii, predictor in enumerate(self.predictors):
            # replace all occurrences of predictor as a whole word
            tf_string = replace_in_string(tf_string, predictor, 'x{}'.format(ii + 1))

        indices = []
        params = []
        # check for existence of parameters as word, get their index of first occurrence
        for param in self.global_params.keys():
            index = get_index_in_string(tf_string, param)
            if index >= 0:
                indices.append(index)
                params.append(param)
        # sort params according to indices
        params = [params[i] for i in sorted(range(len(indices)), key=lambda k: indices[k])]

        # replace parameters in right order, in our working copy as well as in the raw string
        for i_param, param in enumerate(params, 1):
            tf_string = replace_in_string(tf_string, param, 'p{}'.format(i_param))

        # temporary set the name so that math string creation works
        self.translated_name = tf_string

        if self.is_contained:
            self.math_string = ''
            self.index_name = self.tfs['"{}"'.format(self.translated_name)]
        else:
            self.math_string = self._set_math_string(tf_string)

            # work further to get tf key string:
            # ... replace operators
            for key, val in TRANSLATE_DICT.items():
                tf_string = tf_string.replace(key, '_{}_'.format(val))

            # ... eliminate blanks
            tf_string = tf_string.replace(' ', '')

            # ... remove underscores at beginning and end and multiple underscores
            self.translated_name = re.sub(r'[_]{2,}', r'_', tf_string).strip('_')

            # set the index name
            self.index_name = 'transfer_function_{}'.format(self._max_index)

    def _set_math_string(self, prepared_string):
        """
        create mathematical expression for that string, e.g.:
        "func_result(:) = param(1) + param(2) * x(1)%data_p(:) + param(3) * x(2)%data_p(:)"

        Returns
        -------
        """
        math_string = copy(prepared_string)
        remainder = copy(prepared_string)
        for key in TRANSLATE_DICT.keys():
            math_string = replace_in_string(math_string,
                                            key,
                                            ' {} '.format(key),
                                            word_chars=WORD_CHARS_DICT.get(key, ''))
            remainder = replace_in_string(remainder, key, '', word_chars=WORD_CHARS_DICT.get(key, ''))
        # replace the parameters and data arrays by the Fortran syntax
        for i_param in range(self.n_params, 0, -1):
            math_string = math_string.replace('p{}'.format(i_param), 'param({})'.format(i_param))
            remainder = remainder.replace('p{}'.format(i_param), '')
        for i_pred in range(len(self.predictors), 0, -1):
            math_string = math_string.replace('x{}'.format(i_pred), 'x({})%data_p(:)'.format(i_pred))
            remainder = remainder.replace('x{}'.format(i_pred), '')

        # polish some user-defined whitespace hoipolloi
        math_string = re.sub(r'[ ]{2,}', r' ', math_string).strip(' ')
        remainder = re.sub(r'[ ]{2,}', r' ', remainder).strip(' ')
        if remainder:
            raise Exception('Could not successfully parse the following characters : ' +
                            ', '.join(['"{}"'.format(_) for _ in remainder.split(' ')]) +
                            ' in transfer function "{}".'.format(self.raw_tf_string))
        # handle the much more complex where or if clauses
        math_string = self._handle_if_or_where_clause(math_string)

        # insert line breaks at appropriate places in the function if became very long
        if (len(math_string) - 1) > MAX_LINE_LENGTH:
            math_string = break_line(math_string)
        return math_string

    def _handle_if_or_where_clause(self, math_string):
        """create the math string for if and where clauses"""
        triggers = ['if', 'where']
        if re.match(r'|'.join([r'\b{}\b'.format(trigger) for trigger in triggers]), math_string) is None:
            # we prepend the fun_result = string
            math_string = self.MATH_PREFIX_STRING + math_string
        else:
            # now work on special commands like if-clauses and where-clauses
            for trigger in triggers:
                math_string = self._format_trigger_in_math_string(math_string, trigger)

        return math_string

    def _format_trigger_in_math_string(self, math_string, trigger):
        char1 = '('
        char2 = ')'
        # if (...) insert_strig_dict1 ... else ...
        # where (...) insert_strig_dict1 ... else where ...
        insert_strig_dict1 = {'if': 'then', 'where': ''}
        # if (...) then ... else insert_strig_dict2 ...
        # where (...) ... else insert_strig_dict2 ...
        insert_strig_dict2 = {'if': '', 'where': ' where'}
        # special case for where statements without else - set it to no_data in all cases
        insert_strig_dict3 = {'if': '', 'where': '{}{}{}{}{}{}'.format(
            '\n',
            FORTRAN_INDENT * 2,
            'else where',
            '\n',
            FORTRAN_INDENT * 3,
            self.MATH_PREFIX_STRING + 'nodata_dp')}
        end_string = '\n{}end {}'.format(FORTRAN_INDENT * 2, trigger)

        index = 0
        index_addon = 0
        contains_pattern = does_contain_pattern(trigger, math_string)
        if contains_pattern is not None:
            index = contains_pattern.end()
        while contains_pattern is not None:
            # open_parenthesis_counter of levels of nestedness ( char1 increases it, char2 decreases it)
            open_parenthesis_counter = 0
            # whether a pair of parenthesis is contained
            contained_parenthesis = False
            # loop over each character from trigger to end of math_string
            for index in range(index, len(math_string) + 1):
                # continue statements are omitted in each if clause, only execute something in if-clauses
                if math_string[index + index_addon] == ' ':
                    continue
                elif open_parenthesis_counter == 0 and contained_parenthesis:
                    if trigger == 'if' and not math_string[index + index_addon:].startswith('then'):
                        raise Exception('The transfer function is not valid, a "then" needs to follow the',
                                        'logical expression of an if-statement for transfer function:',
                                        self.raw_tf_string)
                    # this code gets inserted, nice linebreak and indentation
                    insert_string = '{}{}{}'.format(
                        ' {}\n'.format(insert_strig_dict1[trigger]),
                        FORTRAN_INDENT * 3,
                        self.MATH_PREFIX_STRING)
                    math_string, index_addon = insert_line_break_in_math_string(
                        math_string[:index + index_addon],
                        insert_string,
                        math_string[index + index_addon + len(insert_strig_dict1[trigger]):],
                        index_addon=index_addon,
                    )

                    # continue looking for else patterns
                    contains_pattern = does_contain_pattern(r'else[ ]+' + trigger, math_string[index + index_addon:])
                    contains_else = does_contain_pattern(r'else', math_string[index + index_addon:])

                    if contains_pattern:
                        insert_string = '{}{}{}'.format(
                            '\n',
                            FORTRAN_INDENT * 2,
                            'else {} '.format(trigger),
                        )
                        math_string, index_addon = insert_line_break_in_math_string(
                            math_string[:index + index_addon + contains_pattern.start()],
                            insert_string,
                            math_string[index + index_addon + contains_pattern.end():],
                            index_addon=index_addon + contains_pattern.start(),
                        )
                        contained_parenthesis = False
                    elif contains_else:
                        insert_string = '{}{}{}{}{}{}'.format(
                            '\n',
                            FORTRAN_INDENT * 2,
                            'else{}'.format(insert_strig_dict2[trigger]),
                            '\n',
                            FORTRAN_INDENT * 3,
                            self.MATH_PREFIX_STRING
                        )
                        math_string, index_addon = insert_line_break_in_math_string(
                            math_string[:index + index_addon + contains_else.start()],
                            insert_string,
                            math_string[index + index_addon + contains_else.end():],
                            end_string,
                            index_addon=index_addon + contains_else.start(),
                        )
                        contains_pattern = None
                        break
                    else:
                        # currently we only support one if or where clause per TF, so if no else is found:
                        math_string = '{}{}{}'.format(math_string, insert_strig_dict3[trigger], end_string)
                        # index_addon += len(insert_string)
                        contains_pattern = None
                        break
                elif math_string[index + index_addon] == char2:
                    open_parenthesis_counter -= 1
                elif math_string[index + index_addon] == char1:
                    open_parenthesis_counter += 1
                    contained_parenthesis = True
                elif not contained_parenthesis:
                    raise Exception('The transfer function is not valid, a "' + char1 +
                                    '" needs to follow a "' + trigger +
                                    '"-statement for transfer function: ' + self.raw_tf_string,
                                    '\nException occurred at character ', str(index + index_addon + 1))

        return math_string

    def insert_index(self, *args):
        """
        helper function to format the values set for indices property
        in TransferFunctionTable in mp_mpr_transfer_func.f90
        """
        return '{}_i4'.format(self.index_name.split('_')[-1], *args)

    def insert_name(self, *args):
        """
        helper function to format the values set for names property
        in TransferFunctionTable in mp_mpr_transfer_func.f90
        """
        return '"{}"'.format(break_line(self.translated_name, string_only=True), *args)


class SourceCode(object):
    """parent class for custom Fortran source code parts"""

    def __init__(self, filepath):
        self.source = self.read_fortran_tf_source(filepath)
        self.tfs = []

    def _retrieve_values(self, key, chars_to_delete=None):
        """retrieve the values of a 1d array parameter property"""
        chars_to_delete = chars_to_delete or []
        # set the pattern we look for in the string
        start_pattern = '{} = ['.format(key)
        # get index of key
        start_index = self.source.find(start_pattern)
        end_index = self.source[start_index + len(start_pattern):].find(']')
        # select the part of the source we are actually interested in
        part = self.source[start_index + len(start_pattern):start_index + len(start_pattern) + end_index]
        # remove the type specification
        if '::' in part:
            part = part.split('::')[-1]
        # remove the unneeded Fortran syntax
        for chars in chars_to_delete:
            part = part.replace(chars, '')

        values = [item.strip(' ') for item in part.split(',')]
        return values

    @staticmethod
    def read_fortran_tf_source(filepath):
        """read the Fortran source file and parse it to a string"""
        with open(filepath) as file:
            source = file.read()
        return source

    def add_tf(self, tf):
        """add a TF to the repo"""
        self.tfs.append(tf)

    @staticmethod
    def _paste_lines(t_lines, pos, p_lines, leading_blanks=''):
        """add a special line marking the inserted code"""
        for paste in p_lines:
            t_lines.insert(pos + 1, paste)
        t_lines.insert(pos + 1, '{}! >>> inserted automatically by Python script'.format(leading_blanks))
        return t_lines


class TFSource(SourceCode):
    """parent class for custom Fortran source code parts"""

    def __init__(self, *args, **kwargs):
        super(TFSource, self).__init__(*args, **kwargs)
        # retrieve the list of tf names as dict
        self.tf_names = self.get_tf_names()

    def get_tf_names(self):
        """retrieve names and indices of transfer_function_names and put them in a dict"""
        # modify the lookup table
        names = self._retrieve_values('names', ['&', '\n', ' ', '//'])
        indices = self._retrieve_values('indices', ['&', '\n', ' ', '//', '_i4'])

        return dict(zip(names, ['transfer_function_{}'.format(item) for item in indices]))

    def get_source(self):
        """update the source code and return it as a list of lines (string)"""
        tf_source = self.source.splitlines()
        # modify the lookup table
        for key in ['indices', 'names']:
            # start with the dimension length
            ii = [tf_source.index(ll) for ll in tf_source if '{} = ['.format(key) in ll][0]
            old_index = re.search(r'dimension\(([0-9]+)\)', tf_source[ii]).group(1)
            tf_source[ii] = tf_source[ii].replace('dimension({})'.format(old_index),
                                                  'dimension({})'.format(int(old_index) + len(self.tfs)))

            # now add the new values to the line
            while True:
                if ']' in tf_source[ii]:
                    # we are at last line of values, we need to insert the values here
                    for tf in self.tfs:
                        format_value = {'indices': tf.insert_index, 'names': tf.insert_name}[key]
                        # if new item would increase line length over max, insert line break
                        if len(tf_source[ii]) + len(format_value()) > MAX_LINE_LENGTH:
                            tf_source[ii] = tf_source[ii].replace(']', ', &')
                            tf_source.insert(ii + 1, ' ' * 6 + format_value() + ']')
                            ii += 1
                        # else simply append
                        else:
                            tf_source[ii] = tf_source[ii].replace(']', ', ' + format_value() + ']')
                    break
                else:
                    ii += 1

        # add function in mo_mpr_transfer_function
        # insert function declaration
        ii = [tf_source.index(ll) for ll in tf_source if 'private' in ll][0]
        for tf in self.tfs[::-1]:
            tf_source.insert(ii + 2, '  public :: {}'.format(tf.index_name))

        # insert function at the end
        for tf in self.tfs[::-1]:
            ii = [tf_source.index(ll) for ll in tf_source if 'end module mo_mpr_transfer_func' in ll][0] - 1
            # insert upside down
            func_string_list = [
                '  end function {}'.format(tf.index_name),
                '    {}'.format(tf.math_string),
                '    allocate(func_result(n))\n',
                '    n = size(x(1)%data_p, kind=i8)',
                '    ',
                '    integer(i8) :: n',
                '    real(dp), dimension(:), intent(in) :: param',
                '    type(InputFieldContainer), dimension(:), intent(in) :: x',
                '    real(dp), dimension(:), allocatable :: func_result',
                '  function {}(x, param) result(func_result)'.format(tf.index_name),
                '  ! ----------------------------------------------------------------------------------------',
            ]
            tf_source = self._paste_lines(tf_source, ii, func_string_list, '  ')

        return tf_source


class DASource(SourceCode):
    def __init__(self, *args, **kwargs):
        super(DASource, self).__init__(*args, **kwargs)

    def get_source(self):
        """update the source code and return it as a list of lines (string)"""
        da_source = self.source.splitlines()

        start_scan = False
        for ii in range(len(da_source)):
            line = da_source[ii].replace(' ', '').lower()[:-1]
            if 'subroutinecall_transfer_func' in line and not start_scan:
                start_scan = True
            if start_scan and 'case' in line:
                for tf in self.tfs[::-1]:
                    set_func_list = [
                        '      data = {}(inputFieldContainers, self%globalParameters)'.format(tf.index_name),
                        '      call self%check_transfer_func_args({}_i4, {}_i4, size(inputFieldContainers))'.format(
                            len(tf.predictors), tf.n_params),
                        '    case(\'{}\')'.format(tf.index_name)]
                    da_source = self._paste_lines(da_source, ii, set_func_list, '    ')
                break

        # insert use statement
        use_lines = [da_source.index(ll) for ll in da_source if 'implicit none' in ll][0]
        for tf in self.tfs[::-1]:
            da_source.insert(use_lines - 1, '  use mo_mpr_transfer_func, only: {}'.format(tf.index_name))

        return da_source


class TFConverter(object):
    """wrapper class for handling the modification of Fortran source files to incorporate new transfer functions"""

    def __init__(self):
        # all the attributes
        self.commandLineArgs = None
        self.do_create_backup = True

        # the (parsed) Fortran source
        self.tf_source = None
        self.da_source = None

        self.mpr_nml = {}
        self.mpr_global_parameter_nml = {}
        # final dict storing the Parameters
        self.global_params = {}

    def parse_args(self):
        """parse command line arguments"""
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description='''Preprocessor script for MPR.

            author: Stephan Thober, Robert Schweppe
            created: Mar 2018''')
        parser.add_argument('-c', '--config_file', action='store', type=pathlib.Path,
                            default=DEFAULT_CONFIG_FILE, dest='config_file', metavar='config_file',
                            help="path to config file for MPR (Default: {})".format(DEFAULT_CONFIG_FILE))
        parser.add_argument('-p', '--parameter_file', action='store', type=pathlib.Path,
                            dest='param_file',
                            help="path to config file with extra parameters for MPR")
        parser.add_argument('-s', '--source_folder', action='store', type=pathlib.Path,
                            default=DEFAULT_SOURCE_FOLDER, dest='src_folder',
                            help="path to source code for MPR (Default: {})".format(DEFAULT_SOURCE_FOLDER))
        parser.add_argument('-t', '--test_mode', action='store_true',
                            default=EXPERIMENTAL, dest='is_test_mode',
                            help="whether to write to temporary files with '{}' suffix (Default: {})".format(
                                MODIFIED_SUFFIX, EXPERIMENTAL))
        parser.add_argument('--clean', action='store_true',
                            dest='from_bak',
                            help="base Fortran source on '{}' files".format(
                                BACKUP_SUFFIX))

        self.commandLineArgs = parser.parse_args()

    @staticmethod
    def check_for_backup(files_to_check_against=None, suffix=''):
        """checks whether a set of file exists and performs some logic"""
        if files_to_check_against is None:
            return True
        else:
            # all files do not exist -> True
            # at least one file exists -> False
            return all([not file.with_suffix(file.suffix + suffix).exists() for file in files_to_check_against])

    def read_source_files(self):
        """reads all Fortran source files to be modified and also all configuration files"""

        da_file = pathlib.Path(self.commandLineArgs.src_folder, FORTRAN_DA_SOURCEFILE)
        transfer_func_file = pathlib.Path(self.commandLineArgs.src_folder, FORTRAN_TF_SOURCEFILE)

        # check if .bak files exist, if not create them
        do_create_backup = self.check_for_backup([da_file, transfer_func_file], suffix=BACKUP_SUFFIX)
        if do_create_backup:
            print('creating backup files at {}'.format(pathlib.Path(da_file.parent, '*' + BACKUP_SUFFIX)))
            copyfile(da_file, da_file.with_suffix(da_file.suffix + BACKUP_SUFFIX))
            copyfile(transfer_func_file, transfer_func_file.with_suffix(transfer_func_file.suffix + BACKUP_SUFFIX))

        # read Fortran source files
        for filepath, target, target_type in zip(
                [da_file, transfer_func_file],
                ['da_source', 'tf_source'],
                [DASource, TFSource],
        ):
            mod_path = filepath.with_suffix(filepath.suffix + MODIFIED_SUFFIX)
            if self.commandLineArgs.is_test_mode and mod_path.exists():
                print('reading Fortran source file from {}'.format(mod_path))
                setattr(self, target, target_type(mod_path))
            elif not do_create_backup and self.commandLineArgs.from_bak:
                bak_path = filepath.with_suffix(filepath.suffix + BACKUP_SUFFIX)
                print('reading Fortran source file from {}'.format(bak_path))
                setattr(self, target, target_type(bak_path))
            else:
                print('reading Fortran source file from {}'.format(filepath))
                setattr(self, target, target_type(filepath))

        # read namelists
        targets = ['mpr_nml', 'mpr_global_parameter_nml']
        for filepath, target in zip([self.commandLineArgs.config_file, self.commandLineArgs.param_file], targets):
            if filepath is not None:
                setattr(self, target, self._read_namelist_source(filepath))

        self.global_params = self._get_parameters()

    @staticmethod
    def _read_namelist_source(filepath):
        """read content from *.nml files and return a dict-like Namelist instance"""

        if not filepath.exists():
            return {}

        parser = f90nml.Parser()
        parser.global_start_index = 1
        return parser.read(filepath)

    def _get_parameters(self):
        """join the parameter names and values from the mpr config dicts to a global dict of parameters"""
        global_params = {}
        for _dict in self.mpr_nml, self.mpr_global_parameter_nml:
            if _dict and list(OPTIONS.keys())[18][0] in _dict:
                # create a dict with {parameter_names: parameter_values}
                global_params.update({str(k): v for k, v in zip(_dict[list(OPTIONS.keys())[18]],
                                                                _dict[list(OPTIONS.keys())[19]]) if k is not None})
        for key in global_params.keys():
            if re.match(r'[xp]{1}[\d]+', key):
                raise Exception('Please do not use parameter names starting with "x" or "p" and followed by a number. '
                                'You provided: ' + key)
        return global_params

    def parse_tfs(self):
        """parse the TFs from the Fortran source code and modify the source code according to configuration"""
        tfs = self.tf_source.tf_names
        # loop over effective params
        transfer_funcs = self.mpr_nml[list(OPTIONS.keys())[23]]
        predictors_key = list(OPTIONS.keys())[22]
        for ii, name in enumerate(self.mpr_nml[list(OPTIONS.keys())[20]]):
            if ii + 1 > len(transfer_funcs):
                transfer_func = None
            else:
                transfer_func = transfer_funcs[ii]

            if name is None or transfer_func is None:
                # no transfer function defined for this effective parameter
                continue

            if predictors_key[-1] in self.mpr_nml[predictors_key[0]]:
                predictors = [_ for _ in self.mpr_nml[predictors_key][ii] if _ is not None]
            else:
                predictors = []

            # handle the case, where there is a transfer function, but no from_data_arrays, then use the name ("self")
            if not predictors:
                predictors = [name]

            # initialize transfer function object
            tf = TF(
                transfer_func,
                predictors,
                self.commandLineArgs.is_test_mode,
                self.global_params,
                tfs
            )
            # generate the unique tf name and replace predictors and parameters in raw string
            tf.translate_func()

            # check if we already have the tf registered (in self.tfs)
            if not tf.is_contained:
                # add the TF to the source code
                self.tf_source.add_tf(tf)
                self.da_source.add_tf(tf)
                # register the new tf in the dict so future tfs know about the existing ones
                tfs['"{}"'.format(tf.translated_name)] = tf.index_name

        print('added {} new transfer functions to code:'.format(len(self.tf_source.tfs)))
        for tf in self.tf_source.tfs:
            print(f'-> "{tf.raw_tf_string}"')

    def write(self):
        """write the Fortran source files"""
        for filename, sourcecode in zip([FORTRAN_DA_SOURCEFILE, FORTRAN_TF_SOURCEFILE],
                                        [self.da_source, self.tf_source]):
            filepath = pathlib.Path(self.commandLineArgs.src_folder, filename)
            if self.commandLineArgs.is_test_mode:
                filepath = filepath.with_suffix(filepath.suffix + MODIFIED_SUFFIX)
            print('writing Fortran source file to {}'.format(filepath))
            with open(filepath, 'w') as file:
                for line in sourcecode.get_source():
                    file.write('{}\n'.format(line))


if __name__ == '__main__':
    converter = TFConverter()
    converter.parse_args()
    converter.read_source_files()
    converter.parse_tfs()
    converter.write()

    print('Done!')
