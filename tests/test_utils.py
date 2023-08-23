from utils import *


def test_resolve_overlaps_coords_list():
    test = [intv(0, 1), intv(1, 3), intv(2, 3), intv(6, 9), intv(7, 20)]
    res_a = [intv(0, 1), intv(2, 3), intv(6, 9), intv(7, 20)]
    assert resolve_overlaps_coords_list(test, 0.9) == res_a
    res_b = [intv(0, 1), intv(2, 3), intv(6, 9)]
    assert resolve_overlaps_coords_list(test, 0) == res_b


def test_get_intervals_overlapping_list():
    test = [intv(0, 1), intv(1, 3), intv(2, 3), intv(6, 9), intv(7, 20)]
    res_a = [(intv(1, 3), intv(2, 3))]
    assert get_intervals_overlapping_list(test, 0.9) == res_a
    res_b = [(intv(1, 3), intv(2, 3)), (intv(6, 9), intv(7, 20))]
    assert get_intervals_overlapping_list(test, 0.0) == res_b
