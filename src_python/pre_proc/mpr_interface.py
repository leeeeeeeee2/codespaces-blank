#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : mpr_interface
Project Name: MPR
Description : contains information namelist structure in MPR (mpr.nml)
Author      : Robert Schweppe (robert.schweppe@ufz.de)
Created     : 24.03.20 10:40
"""

# IMPORTS
from copy import deepcopy
import sys
if sys.version_info.major < 3 and sys.version_info.minor < 7:
    raise Exception('Python 3.7 required, dicts need to be ordered')

# FUNCTIONS


# CLASSES
class Option(object):
    def __init__(self, default=None, options=None, dynamic=False, value_type=str, optional=False):
        """
        Parameters
        ----------
        default: object
            default value for Option, also sets type
        options: iterable
            tuple of possible objects
        dynamic:
            this Option is not set initially
        value_type : type or list of types
            types to convert entries to, can be a list for nested arguments,
            e.g. [list, dict, int] for [{'a': 1, 'b': 2}]
        optional : bool
            whether this option is optional
        """
        self.optional = optional
        self.dynamic = dynamic
        self.options = options or set()

        # set value and default
        if default is not None:
            default = self._test_option(default)
        self.default = deepcopy(default)
        self.value = default
        # set default type
        if self.options:
            self._type_hierarchy_default = self._type_hierarchy(next(iter(self.options)))
        elif self.default is not None:
            try:
                self._type_hierarchy_default = self._type_hierarchy(self.default)
            except KeyError:
                if not isinstance(value_type, str):
                    self._type_hierarchy_default = [value_type] if isinstance(value_type, type) else value_type
                else:
                    raise
        else:
            self._type_hierarchy_default = [value_type] if isinstance(value_type, type) else value_type

    def __repr__(self):
        if self.value is None:
            insert = ''
        elif isinstance(self.value, str):
            insert = "'{}'".format(self.value)
        else:
            insert = self.value.__str__()
        return 'Option({})'.format(insert)

    def add_options(self, options):
        """
        adds options to options

        Parameters
        ----------
        options: iterable

        """
        self.options = self.options.union([self._test_type(x) for x in options])

    @property
    def is_multi(self):
        """
        property, whether value is list or not
        """
        if isinstance(self.value, list):
            try:
                return self._type_hierarchy_default != self._type_hierarchy(self.value)
            except KeyError:
                return False
        else:
            return False

    def is_set(self, ignore_dynamic=False):
        """
        checks if a value (other than None) is set

        Parameters
        ----------
        ignore_dynamic: bool
            whether to ignore unset value in case of dynamic option

        Returns
        -------
        bool
        """
        if ignore_dynamic and self.dynamic:
            return True
        return self.value is not None

    @property
    def is_default(self):
        """
        checks if a the value is the same as default

        Returns
        -------
        bool
        """
        return self.value == self.default

    @property
    def is_optional(self):
        """
        checks if self is optional

        Returns
        -------
        bool
        """
        return self.optional

    def set(self, value):
        """
        wrapper for setting values of different types

        Parameters
        ----------
        value: list or object

        """
        if isinstance(value, list):
            if len(self._type_hierarchy_default) == len(self._type_hierarchy(value)):
                self.set_single(value)
            else:
                self.set_multi(value)
        else:
            self.set_single(value)

    def set_single(self, value):
        """
        set single value (no list)

        Parameters
        ----------
        value: object
        """
        self.value = self._test_val(value)

    def set_multi(self, values):
        """
        set multiple values

        Parameters
        ----------
        values: list of objects

        """
        self.value = [self._test_val(x) for x in values]
        # make singularity if list of length 1
        if len(self.value) == 1:
            self.value = self.value[0]

    def _test_val(self, value):
        """
        tests value against all tests

        Parameters
        ----------
        value: object

        Returns
        -------
        tested and optionally transformed value
        """
        return self._test_option(self._test_type(value))

    def _test_type(self, value):
        """
        tests value against correct type

        Parameters
        ----------
        value: object

        Returns
        -------
        tested and optionally transformed value
        """

        return self._type_converter(value, self._type_hierarchy_default)

    def _type_converter(self, value, types):
        """
        checks if a nested object can be converted to specified (nested) types

        Parameters
        ----------
        value: object
        types: list of types
            must have the same length as the nested object is deep

        Returns
        -------
        [converted(value)]
        """
        error_msg = 'The value {} cannot be set because it ' \
                    'cannot be converted to type {}'.format(value, types[0])
        # if it is a list, call method recursively on each item
        if isinstance(value, list) and issubclass(types[0], list):
            return [self._type_converter(i, types[1:]) for i in value]
        elif isinstance(value, list) and len(value) == 1:
            value = value[0]
        # if it is a dict, call method recursively on each value
        if isinstance(value, dict) and issubclass(types[0], dict):
            return {key: self._type_converter(i, types[1:]) for key, i in value.items()}
        # if requested type and value do not match for list or dict, raise
        elif isinstance(value, dict):
            raise TypeError(error_msg)
        # if it needs to be transformed to a boolean, check if in boolean strings
        elif types[0] == bool:
            return value in ['True', '.true.', '.TRUE.', 'TRUE', 'true']
        else:
            # apply type on value
            if types[0] is None:
                return value
            else:
                try:
                    # noinspection PyCallingNonCallable
                    return types[0](value)
                except (TypeError, ValueError):
                    # raise if fails
                    raise Exception(error_msg)

    def _test_option(self, value):
        """
        tests value against options

        Parameters
        ----------
        value: object

        Returns
        -------
        tested and optionally transformed value
        """

        if self.options and value not in self.options:
            raise KeyError('The value {} is not in "{}"'.format(value,
                                                                ', '.join(sorted(map(str, self.options)))))
        return value

    def __contains__(self, item):
        if self.is_multi:
            return item in self.value
        else:
            return item == self.value

    @staticmethod
    def _type_hierarchy(nested_obj):
        types = []
        obj = nested_obj
        while True:
            types.append(type(obj))
            if isinstance(obj, list):
                if not len(obj):
                    raise KeyError('the type of an option cannot be inferred. Please '
                                   'do not provide empty lists.')
                obj = obj[-1]
            elif isinstance(obj, dict):
                obj = obj[list(obj.keys())[-1]]
            else:
                break
        return types

    def modify(self, index):
        """
        modifies self.value
        The idea is to supply an index and this subroutine modifies self.value accordingly at this position
        and optionally creates new lists or None values for intermediate positions

        Parameters
        ----------
        value: object
        index: list of int
            list of indices indicating where to insert/append the object in the current structure

        """
        # self._type_converter(value, self._type_hierarchy_default)
        len_value = 0
        slice_chain = ''
        for iidx, idx in enumerate(index):
            # if the current structure is not a list, indices make no sense
            if not issubclass(self._type_hierarchy_default[iidx], list):
                raise Exception('The {}th index {} from the indices {} is not valid. '
                                'The target structure is not iterable but of '
                                'type {}.'.format(iidx, idx, ', '.join(map(str, index)),
                                                  self._type_hierarchy_default[iidx]))
            if iidx > 0:
                slice_chain = '[' + ']['.join(map(str, [i - 1 for i in index[:iidx]])) + ']'
            # is the index within range of list, then the case represents insertion
            if self.value is not None:
                len_value = eval('len(self.value{slice_chain})'.format(slice_chain=slice_chain))
            # is the index 1 + range of list, then the case represents appending
            if idx >= (len_value + 1):
                # is this the last item in hierarchy, simply use dummy None
                if (iidx + 1) == len(index):
                    eval('self.value{slice_chain}.extend([None]*{multiplier})'.format(slice_chain=slice_chain,
                                                                                      multiplier=idx - len_value))
                # is this not the last item in hierarchy, allow for further nesting
                else:
                    eval('self.value{slice_chain}.extend([[None]]*{multiplier})'.format(slice_chain=slice_chain,
                                                                                        multiplier=idx - len_value))
        slice_chain = '[' + ']['.join(map(str, [i - 1 for i in index])) + ']'
        # TODO: ask someone if this can be overcome
        exec('self.value{slice_chain} = self._type_converter(value, self._type_hierarchy_default[iidx+1:])'.format(
            slice_chain=slice_chain))


# GLOBAL VARIABLES
# noinspection PyTypeChecker,PyTypeChecker
OPTIONS = {
    # 0
    ('main', 'coordinate_aliases'): Option([[]], value_type=[list, list, str]),
    ('main', 'out_filename'): Option(value_type=str),
    ('main', 'write_weights'): Option(value_type=bool, optional=True),
    ('main', 'read_weights'): Option(value_type=bool, optional=True),
    ('main', 'check_for_nodatavalue'): Option(value_type=bool, optional=True),
    ('coordinates', 'coord_name'): Option([], value_type=[list, str]),
    ('coordinates', 'coord_stagger'): Option([], value_type=[list, str], optional=True),
    ('coordinates', 'coord_from_file'): Option([], value_type=[list, str], optional=True),
    ('coordinates', 'coord_from_values'): Option([], value_type=[list, list, float], optional=True),
    ('coordinates', 'coord_from_values_bound'): Option([], value_type=[list, float], optional=True),
    # 10
    ('coordinates', 'coord_from_range_start'): Option([], value_type=[list, float], optional=True),
    ('coordinates', 'coord_from_range_step'): Option([], value_type=[list, float], optional=True),
    ('coordinates', 'coord_from_range_count'): Option([], value_type=[list, int], optional=True),
    ('coordinates', 'coord_attribute_names'): Option([[]], value_type=[list, list, str], optional=True),
    ('coordinates', 'coord_attribute_values'): Option([[]], value_type=[list, list, str], optional=True),
    ('coordinates', 'coord_proj_string'): Option([[]], value_type=[list, str], optional=True),
    ('coordinates', 'coord_sub_dims'): Option([[]], value_type=[list, list, str], optional=True),
    ('coordinates', 'coord_unit'): Option([[]], value_type=[list, str], optional=True),
    ('parameters', 'parameter_names'): Option([], value_type=[list, str], optional=True),
    ('parameters', 'parameter_values'): Option([], value_type=[list, float], optional=True),
    # 20
    ('data_arrays', 'name'): Option([], value_type=[list, str]),
    ('data_arrays', 'from_file'): Option([], value_type=[list, str], optional=True),
    ('data_arrays', 'from_data_arrays'): Option([[]], value_type=[list, list, str], optional=True),
    ('data_arrays', 'transfer_func'): Option([], value_type=[list, str], optional=True),
    ('data_arrays', 'target_coord_names'): Option([[]], value_type=[list, list, str], optional=True),
    ('data_arrays', 'upscale_ops'): Option([[]], value_type=[list, list, str], optional=True),
    ('data_arrays', 'limits'): Option([[]], value_type=[list, list, float], optional=True),
    ('data_arrays', 'to_file'): Option([], value_type=[list, bool], optional=True),
    ('data_arrays', 'from_parameter_values'): Option([], value_type=[list, list, float], optional=True),
    ('data_arrays', 'from_parameter_names'): Option([], value_type=[list, list, str], optional=True),
    # 30
    ('upscalers', 'upscaler_name'): Option([], value_type=[list, str], optional=True),
    ('upscalers', 'from_weight_file'): Option([], value_type=[list, str], optional=True),
    ('upscalers', 'subcell_ids_field_name'): Option([], value_type=[list, str], optional=True),
    ('upscalers', 'weights_field_name'): Option([], value_type=[list, str], optional=True),
    ('upscalers', 'n_subcells_field_name'): Option([], value_type=[list, str], optional=True),
    ('data_arrays', 'transfer_func_label'): Option([], value_type=[list, str], optional=True),
}

# SCRIPT
if __name__ == '__main__':
    pass
