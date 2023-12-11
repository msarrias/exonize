"""
Module that will perform some profiling tasks, namely wrapping of functions for now with a simple profiling edit.
"""

import pstats
from pathlib import Path

PROFILE_PATH = (Path(__file__).resolve().parent / 'cProfile_dump_stats').as_posix()


def get_run_performance_profile(filename_for_profile_output):
    """
    read the profile data saved in filename_for_profile_output into a pstats.Stats object,
    sort the stats, print it and return the object.
    :param filename_for_profile_output:
    :return:
    """
    profile_stat = pstats.Stats(filename_for_profile_output)
    profile_stat.sort_stats('tottime')
    profile_stat.print_stats()
    return profile_stat
