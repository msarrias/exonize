"""
Module that will perform some profiling tasks, namely wrapping of functions for now with a simple profiling edit.
"""
import pstats


def get_run_performance_profile(
        filename_for_profile_output,
        pr
):
    with open(filename_for_profile_output, 'w') as f:
        pstats.Stats(pr, stream=f).strip_dirs().sort_stats('tottime').print_stats()
