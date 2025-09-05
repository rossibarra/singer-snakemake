"""
Utility functions for converting between formats.

Part of https://github.com/nspope/singer-snakemake.
"""

import numpy as np
import msprime


def bitmask_to_arrays(bitmask: np.ndarray, *, insert_breakpoints: np.ndarray = None) -> msprime.RateMap:
    """
    Convert a bitmask to a binary interval mask, with additional
    `insert_breakpoints` if not `None`. 
    """
    assert np.issubdtype(bitmask.dtype, bool)
    changepoints = np.flatnonzero(bitmask[1:] != bitmask[:-1]) + 1
    changepoints = np.append(np.append(0, changepoints), bitmask.size)
    if insert_breakpoints is not None:
        assert np.issubdtype(insert_breakpoints.dtype, int)
        assert insert_breakpoints.min() >= 0 and insert_breakpoints.max() <= bitmask.size
        changepoints = \
            np.unique(np.append(insert_breakpoints, changepoints))
    values = bitmask[changepoints[:-1]].astype(float)
    return values, changepoints


def ratemap_to_text(ratemap: msprime.RateMap, *, replace_nan_with: float = 0.0) -> str:
    """
    Write a ratemap to a headerless, three-column text file that is
    left coordinate, right coordinate, per-base rate
    """
    text = []
    for left, right, rate in zip(
        ratemap.position[:-1],  # FIXME: use ratemap.left/right
        ratemap.position[1:], 
        ratemap.rate,
    ):
        if np.isnan(rate): rate = replace_nan_with
        text.append(f"{int(left)} {int(right)} {rate:.16f}") # FIXME: precision FIXME: no int
    text = "\n".join(text) + "\n"
    return text 
