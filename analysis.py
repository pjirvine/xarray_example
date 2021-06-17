"""
This python script contains analysis functions used in Irvine et al. 2019 and 2020
"""

import numpy as np
from scipy.stats import ttest_ind_from_stats
import itertools

def ttest_sub(mean_1, std_1, nyears_1, mean_2, std_2, nyears_2, equal_var=True):

    """
    Sub-routine to call ttest_ind_from_stats from scipy
    Checks that shapes match and turns integer years into correct format
    returns pvalue.
    """

    # Convert nobs type
    nyears_1 = int(nyears_1)
    nyears_2 = int(nyears_2)

    # Create arrays like others for nobs
    nobs1_arr = (nyears_1-1) * np.ones_like(mean_1)
    nobs2_arr = (nyears_2-1) * np.ones_like(mean_1)

    """
    # ttest_ind_from_stats
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind_from_stats.html
    """

    ttest_out = ttest_ind_from_stats(mean_1, std_1, nobs1_arr, mean_2, std_2, nobs2_arr)

    # An array of p-values matching the shape of the input arrays
    pvalue_out = ttest_out[1]

    return pvalue_out

