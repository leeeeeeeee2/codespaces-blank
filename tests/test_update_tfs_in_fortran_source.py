# encoding: utf-8

"""
File Name   : test_create_graphviz_plot.py
Project Name: MPR
Description : test suite for update_tfs_in_fortran_source.py script
Author      : Robert Schweppe (robert.schweppe@ufz.de)
Created     : 2019-03-30 13:30
"""

# IMPORTS
import pytest
from src_python.pre_proc.update_tfs_in_fortran_source import TFConverter, TF, DASource, TFSource, \
    replace_in_string, get_index_in_string, WORD_BOUNDARIES
import copy
from textwrap import dedent

class TestTFConverter:
    '''
    @staticmethod
    def check_fortran_code_for_completion(transfer_func_file, print_status=False):
        pass
        # read mo_mpr_transfer_funcs and check results
        # TODO: this was taken out, because the existing code does not work with more complex where and if clauses
        # also, the checking of the proper formulation of if and where clauses (if (...) then (...) else (...) end if)
        # was moved to the part where the math_string is created, so that should not be necessary anymore
        # maybe this part should be moved to a unit test anyway (pytest)

        if self.is_test_mode:
            fi = open(transfer_func_file + MODIFIED_SUFFIX, 'r')
        else:
            fi = open(transfer_func_file, 'r')
        lines = fi.readlines()
        fi.close()
        print('\n  start checking transfer functions: ')
        current_transfer_func = []

        for ll, line in enumerate(lines):
            if not re.search('func_result', line.replace(' ', '')) is None:
                # check whether line before contains if statement
                if not re.search('if.*then', lines[ll - 1][:-1].replace(' ', '')) is None:
                    start_if = True
                    store_func = lines[ll - 1][:-1].replace(' ', '')
                elif not re.search('endif', lines[ll + 1][:-1].replace(' ', '')) is None:
                    # if clause already processed
                    continue
                else:
                    start_if = False
                    store_func = ''
                # get the entire function
                store_func, ll, line = self._evaluate_func(store_func, ll, line, lines)
                #
                if start_if:
                    # check whether line contains else
                    if not re.search('else', line[:-1].replace(' ', '')) is None:
                        store_func += line[:-1].replace(' ', '')
                        # get the second part of the function function
                        store_func, ll, line = self._evaluate_func(store_func, ll + 1, lines[ll + 1], lines)
                    if not re.search('endif', line[:-1].replace(' ', '')) is None:
                        store_func += line[:-1].replace(' ', '')
                    else:
                        raise ValueError('Incomplete if statement in mo_mpr_transfer_func')

                if re.search('if.*if', store_func) and store_func[0] != '!':
                    current_transfer_func.append(store_func)
                elif re.search('func_result.*=', store_func) and store_func[0] != '!':
                    current_transfer_func.append(store_func)

        # remove func_result(i)= pattern
        current_transfer_func = [re.sub('func_result...=', '', func) for func in current_transfer_func]
        inverted_transfer_func = [self._invert(func) for func in current_transfer_func]

        if self.raw_tf_string.replace(' ', '') in inverted_transfer_func:
            print('')
            print('  Transfer function: ')
            print('  {}'.format(self.raw_tf_string))
            if print_status:
                print('  successfully added to Fortran code!')
            else:
                print('  already contained in Fortran code!')
            print('')
            print('  Please proceed with re-compiling mpr!')
            print('')
        else:
            raise ValueError('***ERROR: Transfer function {} has not been added correctly!'.format(self.raw_tf_string))

    @staticmethod
    def _evaluate_func(store_func, ll, line, lines):
        store_func += line[:-1].replace(' ', '')

        cont_line = True
        while cont_line:
            if line[-1] == '&':
                cont_line = True
                store_func += lines[ll + 1][:-1].replace(' ', '')
                ll += 1
            else:
                cont_line = False
        ll += 1
        return store_func, ll, lines[ll]

    def _invert(self, func):
        inv_func = copy(func)

        # replace parameters
        for pp in range(self.n_params + 1):
            inv_func = inv_func.replace('param({})'.format(pp), 'p{}'.format(pp))

        # replace predictors
        for ii in range(len(self.predictors)):
            inv_func = inv_func.replace('x({})%data_p(i)'.format(ii + 1), self.predictors[ii])

        return inv_func
    '''


    def test_one(self):
        assert True


class TestTF:
    def test_translation(self):
        raw_strings = (
            'p1 + p2 * d1 + p3 * d2',
            "thetas1+thetas2*CLYPPT_M+thetas3*BLDFIE_M+thetas4*SLTPPT_M**(2.0)+thetas5*ORCDRC_M**(2.0)+"
            "thetas6*CLYPPT_M**(-1.0)+thetas7*SLTPPT_M**(-1.0)+thetas8*log(SLTPPT_M)+thetas9*ORCDRC_M*CLYPPT_M-"
            "thetas10*BLDFIE_M*CLYPPT_M-thetas11*BLDFIE_M*ORCDRC_M+thetas12*topsoil*SLTPPT_M",
            'p1+exp((p2+d1)*p3)',
            'p1+exp(p3*(p2+d1))',
            'PTF_Ks_curveSlope * exp ((PTF_Ks_constant + PTF_Ks_sand * SAND + PTF_Ks_clay * CLAY) * log (10))',
            'where((z.lower_bound + (z.upper_bound - z.lower_bound) / 2.0) <= topsoil_boundary) 1.0 else 0.0',
            'where (x1 > l1 .and. x1 < l2) l5 else where (x1 < l3) l6 else where (x1 < l4) l7',
            'p1 + asin(d1)',
            'p1 + atan2(d1)',
            'p1 + tanh(d1)',
        )
        predictors_sets = (
            ['d1', 'd2'],
            ['CLYPPT_M', 'SLTPPT_M', 'BLDFIE_M', 'ORCDRC_M', 'topsoil'],
            ['d1'],
            ['d1'],
            ['SAND', 'CLAY'],
            ['SLTPPT_M', 'z.upper_bound', 'z.lower_bound'],
            ['x1'],
            ['d1'],
            ['d1'],
            ['d1'],
        )
        global_params_set = (
            None,
            {'thetas0': 0,
             'thetas1': 1,
             'thetas2': 2,
             'thetas3': 3,
             'thetas4': 4,
             'thetas5': 5,
             'thetas6': 6,
             'thetas7': 7,
             'thetas8': 8,
             'thetas9': 9,
             'thetas10': 10,
             'thetas11': 11,
             'thetas12': 12,
             '-1.0': -1,
             '2.0': 2},
            None,
            None,
            {'10': 10.0,
             'PTF_Ks_curveSlope': 0.000007055555,
             'PTF_Ks_constant': -0.6,
             'PTF_Ks_sand': 0.0126,
             'PTF_Ks_clay': -0.0064,
             },
            {'topsoil_boundary': 0.0501, '0.0': 0.0, '1.0': 1.0, '2.0': 2.0},
            {
                'l2': 2,
                'l6': 6,
                'l1': 1,
                'l7': 7,
                'l3': 3,
                'l5': 5,
                'l4': 4,
            },
            None,
            None,
            None,
        )
        translated_names = (
            'p1_pl_p2_ti_x1_pl_p3_ti_x2',
            'p1_pl_p2_ti_x1_pl_p3_ti_x3_pl_p4_ti_x2_po_bs_p5_be_pl_p6_ti_x4_po_bs_p5_be_pl_p7_ti_x1_po'
            '_bs_p8_be_pl_p9_ti_x2_po_bs_p8_be_pl_p10_ti_l2_bs_x2_be_pl_p11_ti_x4_ti_x1_mi_p12_ti_x3'
            '_ti_x1_mi_p13_ti_x3_ti_x4_pl_p14_ti_x5_ti_x2',
            'p1_pl_ex_bs_bs_p2_pl_x1_be_ti_p3_be',
            'p1_pl_ex_bs_p3_ti_bs_p2_pl_x1_be_be',
            'p1_ti_ex_bs_bs_p2_pl_p3_ti_x1_pl_p4_ti_x2_be_ti_l2_bs_p5_be_be',
            #'wh_bs_bs_x3_pl_bs_x2_mi_x3_be_di_p1_be_le_p2_be_p3_el_p4',
            'wh_bs_bs_x3_pl_bs_x2_mi_x3_be_di_p1_be_le_p2_be_p3_el_p4',
            'wh_bs_x1_gt_p1_ad_x1_lt_p2_be_p3_el_wh_bs_x1_lt_p4_be_p5_el_wh_bs_x1_lt_p6_be_p7',
            'p1_pl_as_bs_x1_be',
            'p1_pl_au_bs_x1_be',
            'p1_pl_tx_bs_x1_be',
        )
        math_strings = (
            'func_result(:) = param(1) + param(2) * x(1)%data_p(:) + param(3) * x(2)%data_p(:)',
            dedent('''\
            func_result(:) = param(1) + param(2) * x(1)%data_p(:) + param(3) * x(3)%data_p(:) + param(4) * &
                  x(2)%data_p(:) ** ( param(5) ) + param(6) * x(4)%data_p(:) ** ( param(5) ) + param(7) * x(1)%data_p(:) &
                  ** ( param(8) ) + param(9) * x(2)%data_p(:) ** ( param(8) ) + param(10) * log ( x(2)%data_p(:) ) + &
                  param(11) * x(4)%data_p(:) * x(1)%data_p(:) - param(12) * x(3)%data_p(:) * x(1)%data_p(:) - param(13) &
                  * x(3)%data_p(:) * x(4)%data_p(:) + param(14) * x(5)%data_p(:) * x(2)%data_p(:)'''),
            'func_result(:) = param(1) + exp ( ( param(2) + x(1)%data_p(:) ) * param(3) )',
            'func_result(:) = param(1) + exp ( param(3) * ( param(2) + x(1)%data_p(:) ) )',
            dedent('''\
            func_result(:) = param(1) * exp ( ( param(2) + param(3) * x(1)%data_p(:) + param(4) * &
                  x(2)%data_p(:) ) * log ( param(5) ) )'''),
            dedent('''\
            where ( ( x(3)%data_p(:) + ( x(2)%data_p(:) - x(3)%data_p(:) ) / param(1) ) <= param(2) ) 
                  func_result(:) = param(3)
                else where
                  func_result(:) = param(4)
                end where'''),
            dedent('''\
            where ( x(1)%data_p(:) > param(1) .and. x(1)%data_p(:) < param(2) ) 
                  func_result(:) = param(3)
                else where ( x(1)%data_p(:) < param(4) ) 
                  func_result(:) = param(5)
                else where ( x(1)%data_p(:) < param(6) ) 
                  func_result(:) = param(7)
                else where
                  func_result(:) = nodata_dp
                end where'''),
            'func_result(:) = param(1) + asin ( x(1)%data_p(:) )',
            'func_result(:) = param(1) + atan2 ( x(1)%data_p(:) )',
            'func_result(:) = param(1) + tanh ( x(1)%data_p(:) )',
        )
        for raw_string, predictors, translated_name, global_params, math_string in \
                zip(raw_strings, predictors_sets, translated_names, global_params_set, math_strings):
            tf = TF(raw_string, predictors=predictors, global_params=global_params)
            tf.translate_func()
            assert tf.translated_name == translated_name
            assert tf.math_string == math_string

    def test_where_strings(self):
        raw_strings = (
            'where (d1 > p1) d1 else d2',
            'where (d1 > p1) d1 else where (d1 > p2) p2 else d2',
        )
        math_strings = (
            dedent('''\
                where ( x(1)%data_p(:) > param(1) ) 
                      func_result(:) = x(1)%data_p(:)
                    else where
                      func_result(:) = x(2)%data_p(:)
                    end where'''),
            dedent('''\
                where ( x(1)%data_p(:) > param(1) ) 
                      func_result(:) = x(1)%data_p(:)
                    else where ( x(1)%data_p(:) > param(2) ) 
                      func_result(:) = param(2)
                    else where
                      func_result(:) = x(2)%data_p(:)
                    end where'''),
        )

        for raw_string, math_string in zip(raw_strings, math_strings):
            tf = TF(raw_string, predictors=['d1', 'd2'], global_params={'p1': 0, 'p2': 1})
            tf.translate_func()
            assert tf.math_string == math_string

    def test_where_strings_failures(self):
        raw_strings = (
            # missing parenthesis
            'where d1 > y1 d1 else d2',
            'where d1 > y1) d1 else d2',
            # extra parameter or dataarray
            'where (d1 > y3) d1 else d2',
            # typo
            'where (d1 > D2) d1 else d2',
            # forbidden names
            # TODO: this fails and only should if those are parsed from namelist
            #'where (d1 > x1) d1 else d2',
            #'where (d1 > p1) d1 else d2',
        )
        for raw_string in raw_strings:
            tf = TF(raw_string, predictors=['d1', 'd2'], global_params={'y1': 0, 'y2': 1, 'p1': -1, 'x1': -1})
            with pytest.raises(Exception):
                tf.translate_func()

class TestTFSource:
    def test_one(self):
        assert True

class TestDASource:
    def test_one(self):
        assert True

def test_replace_in_string():
    search_strings = (
        '1**a**(',
    )
    string_patterns = (
        '**',
    )
    replacements = (
        ' ** ',
    )
    boundaries = (
        '*=',
    )
    replaceds = (
        '1 ** a ** (',
    )
    for search_string, string_pattern, replacement, boundary, replaced in \
            zip(search_strings, string_patterns, replacements, boundaries, replaceds):
        assert replace_in_string(search_string, string_pattern, replacement, boundary) == replaced

def test_get_index_in_string():
    long_strings = (
         's ss sss',
         'sss ss s',
         '**12*s',
    )
    parts = (
        'ss',
        'ss',
        '*',
    )
    args = (
        [],
        [],
        ['*'],
    )
    results = (
        2,
        4,
        4,
    )
    for long_string, part, arg, result in zip(long_strings, parts, args, results):
        assert get_index_in_string(long_string, part, *arg) == result

