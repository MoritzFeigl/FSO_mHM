import pytest
import xarray as xr
from utils.xarray_combine_patched import get_index_from_ranges, union_non_overlapping_ranges
from functools import partial

"""
File Name   : test_patch
Project Name: 2020_FSO_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 27.08.20 10:01
"""

# IMPORTS

# GLOBAL VARIABLES

# FUNCTIONS

# CLASSES
class TestClass:
    def test_get_index_from_ranges(self):
        ref_ranges = [
            (-1.5, 0),
            (0, 1.5),
            (2.5, 4.0)
        ]
        ranges = [
            (0, 1),
            (1.0, 1.5),
            (-1.0, -0.5),
            (2.5, 4.0),
            (-1.5, 0),
        ]
        order = [1, 1, 0, 2, 0]
        assert get_index_from_ranges(ranges, ref_ranges) == order

        descending_ref_ranges = sorted((sorted(item, reverse=True) for item in ref_ranges), reverse=True)
        descending_ranges = [sorted(item, reverse=True) for item in ranges]
        reversed_order = [(_-1)*(-1)+1 for _ in order]

        assert get_index_from_ranges(descending_ranges, descending_ref_ranges, ascending=False) == reversed_order

        ranges = [(4260050.0, 4299950.0), (4300050.0, 4339950.0), (4340050.0, 4375050.0), (4228050.0, 4259950.0), (4260050.0, 4299950.0), (4300050.0, 4339650.0), (4230950.0, 4259950.0), (4260050.0, 4260750.0), (4281050.0, 4299950.0), (4300050.0, 4339950.0), (4340050.0, 4363550.0), (4252150.0, 4259950.0), (4260050.0, 4299950.0), (4300050.0, 4339950.0), (4340050.0, 4376850.0), (4233950.0, 4259950.0)]
        ref_ranges = [[4228050.0, 4259950.0], [4260050.0, 4299950.0], [4300050.0, 4339950.0], [4340050.0, 4376850.0]]
        order = [1, 2, 3, 0, 1, 2, 0, 1, 1, 2, 3, 0, 1, 2, 3, 0]
        assert get_index_from_ranges(ranges, ref_ranges) == order

    def test_union_non_overlapping_ranges(self):
        # those are taken from https://stackoverflow.com/questions/15273693/union-of-multiple-ranges/15273749#15273749
        ranges = [
            [(7, 10), (11, 13), (11, 15), (14, 20), (23, 39)],
            [(7.0, 10.0), (11.0, 13.0), (11.0, 15.0), (14.0, 20.0), (23.0, 39.0)],
            [(2, 4), (1, 6)],
        ]
        references = [
            [[7, 20], [23, 39]],
            [[7.0, 10.0], [11.0, 20.0], [23.0, 39.0]],
            [[1, 6]],
        ]
        for _range, reference in zip(ranges, references):
            assert union_non_overlapping_ranges(_range) == reference
            reversed_ranges = [sorted(item, reverse=True) for item in _range]
            reversed_reference = sorted([sorted(item, reverse=True) for item in reference], reverse=True)
            assert union_non_overlapping_ranges(reversed_ranges, ascending=False) == reversed_reference

