"""
This python script contains analysis functions used in Irvine et al. 2019 and 2020
"""

import numpy as np
from scipy.stats import ttest_ind_from_stats
import itertools
import xarray as xr

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

# A wrapper for the above that works with xarray datasets.
def ttest_func(ds_mean_1, ds_std_1, ds_mean_2, ds_std_2, var, num_years, p_thresh=0.05):
    # p_thresh defaults to 0.05 for a 95% T-test.
    # ttest_sub returns a numpy array of P-values, where P is between 0 and 1. for 95% significance P is below 0.05
    ttest_pvalue = ttest_sub(ds_mean_1[var],ds_std_1[var],num_years,ds_mean_2[var],ds_std_2[var],num_years)

    # Let's put the ttest results into the same format as our xarray datasets
    ds_ttest = xr.full_like(ds_mean_1, 0.0) # copy dataset format from ds_mean_1 and set data values to 0.
    ds_ttest.rename(name_dict={var:'T-test'}) # rename the variable to p_value
    ds_ttest['T-test'] = (['lat','lon'],ttest_pvalue < p_thresh) # Fill in the blank values with our ttest results
    
    # Return ttest data array
    return ds_ttest['T-test']

def ttest_thresh(num_years, p_thresh=0.05, sign_figs=3):
    # returns the number of standard deviations needed to pass a T-test for a given number of years in sample
    # default, 0.05: 95% T-Test, 3 significant figures.
    
    num_steps = 10**sign_figs

    ttest_sds = None
    for IDX in [idx / num_steps for idx in range(num_steps)]:
        if ttest_sds is None:
            if ttest_sub(IDX, 1, num_years, 0, 1, num_years, equal_var=True) < p_thresh:
                ttest_sds = IDX
    return ttest_sds

"""
###
Define functions which determine better off, worse off, don't know and all sub-types
###
"""

def bools_x8_3ttests(test_1,test_2,test_3):
    
    """
    This function produces booleans for all 8 combinations of the three input T-tests:
    test_1 = abs(SRM_anom), abs(CO2_anom), SRM_std, CO2_std
    test_2 = CO2, ctrl
    test_3 = SRM, ctrl
    """
    
    return {'FFF': np.logical_not(test_1) * np.logical_not(test_2) * np.logical_not(test_3),
            'FFT': np.logical_not(test_1) * np.logical_not(test_2) * (test_3),
            'FTF': np.logical_not(test_1) * (test_2) * np.logical_not(test_3),
            'FTT': np.logical_not(test_1) * (test_2) * (test_3),
            'TFF': (test_1) * np.logical_not(test_2) * np.logical_not(test_3),
            'TFT': (test_1) * np.logical_not(test_2) * (test_3),
            'TTF': (test_1) * (test_2) * np.logical_not(test_3),
            'TTT': (test_1) * (test_2) * (test_3),}

# This snippet allows keys + values in dictionaries to be read as variables of form X.key instead of D['key']
class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)

def types_groups_from_bools(bool_dict, ratio):
    
    """
    This function takes the output from bools_x8_3ttests and the ratio and returns them with standard names and in groups
    
    See hypothesis document for details.
    
    ratio = SRM_anom / CO2_anom
    |>1| = exacerbated
    |<1| = moderated
    +ve = same sign
    -ve = reversed sign
    """
    
    """
        1: SRM vs CO2,   2: CO2 vs CTRL,   3: SRM vs CTRL
    """
    type_dict = {'A1': bool_dict['TFF'],
                 'A2': bool_dict['TFT'],
                 'A3': bool_dict['TTF'],
                 'A4': bool_dict['TTT'],
                 'B1': bool_dict['FFF'],
                 'B2': bool_dict['FFT'],
                 'B3': bool_dict['FTF'],
                 'B4': bool_dict['FTT'],}

    td = Bunch(type_dict)
    
    group_dict = {'all': td.A1 + td.A2 + td.A3 + td.A4 + td.B1 + td.B2 + td.B3 + td.B4,
                  'dont_know': td.A1 + td.B1 + td.B2 + td.B3 + td.B4,
                  'better_off': td.A3 + td.A4 * (abs(ratio) < 1),
                  'worse_off': td.A2 + td.A4 * (abs(ratio) > 1),
                  'not_small': td.A2 + td.A3 + td.A4 + td.B4,
                  'certain': td.A2 + td.A3 + td.A4,
                  'dont_know_small': td.A1 + td.B1 + td.B2 + td.B3,
                  'dont_know_big_none': td.B4 * (ratio > 0),
                  'dont_know_big_over': td.B4 * (ratio < 0),
                  'better_off_perfect': td.A3,
                  'better_off_under': td.A4 * (0 < ratio) * (ratio < 1),
                  'better_off_over': td.A4 * (-1 < ratio) * (ratio < 0),
                  'worse_off_novel': td.A2,
                  'worse_off_exacerbate': td.A4 * (ratio > 1),
                  'worse_off_too_much': td.A4 * (ratio < -1),
                  'all_over': (td.A4 + td.B4) * (ratio < 0),
                 }

    return type_dict, group_dict

"""
This function calculates the fraction which are better, worse and don't know
"""

def better_worse_off(SRM_mean, SRM_std, CO2_mean, CO2_std, CTRL_mean, CTRL_std, nyears, ttest_level):
    
    # define anomalies
    CO2_anom = CO2_mean - CTRL_mean
    SRM_anom = SRM_mean - CTRL_mean

    # ratio of anomalies
    try: # check for divide by zero error and create very big number instead
        ratio = SRM_anom / CO2_anom
    except ZeroDivisionError:
        ratio = np.sign(SRM_anom) * 9.999*10**99

    # absolute_double_anom T-Test
    ttest_1 = ttest_sub(abs(SRM_anom), SRM_std, nyears,
                        abs(CO2_anom), CO2_std, nyears) < ttest_level
    # CO2, ctrl T-Test
    ttest_2 = ttest_sub(CO2_mean, CO2_std, nyears,
                        CTRL_mean, CTRL_std, nyears) < ttest_level
    # SRM, ctrl T-Test
    ttest_3 = ttest_sub(SRM_mean, SRM_std, nyears,
                        CTRL_mean, CTRL_std, nyears) < ttest_level
    
    # This geomip_data.py function returns dictionary of combinations of results
    bool_dict = bools_x8_3ttests(ttest_1,ttest_2,ttest_3)
    
    # This geomip_data.py function returns dictionary of types of results
    type_dict, group_dict = types_groups_from_bools(bool_dict, ratio)
    
    return group_dict['better_off'], group_dict['worse_off'], group_dict['dont_know']
# End def

"""
This function returns an array where each point is categorized based on whether it is moderated or exacerbated.
"""

def mod_ex_as_numbers(sg_mean, cc_mean, ctrl_mean, sg_std, cc_std, ctrl_std, years, p_value=0.05, sign=False):
    """
    This function reads numpy arrays for geo, co2, and control and calculates which cells are moderated (absolute magnitude of change reduceded)
    and exacerbated (absolute magnitude of change increased), producing a numpy array with the following values, corresponding to results in the
    category dictionary:
    -2: exacerbated, significant, e.g. +5% to +6% (or -6%)
    -1: exacerbated, insignificant, e.g. +5% to +5.1% (or -5.1%)
    1: moderated, insignificant, e.g. +5% to +4.9% (or -4.9%)
    2: moderated, significant, e.g. +5% to 4% (or -4%)
    [OPTIONAL, only calculated if sign=TRUE, these categories capture the parenthetical results listed above]
    3: moderated, significant, sign-change, e.g. -4%
    4: moderated, insignificant, sign-change, e.g. -4.9%
    5: exacerbated, sign changed, insignificant change, e.g. -5.1%
    6: exacerbated, sign changed, significant increase in magnitude, -6%

    Inputs:
    mean and standard deviation as numpy arrays for the solar geo, climate change and control cases.
    years: number of years for T-Test, e.g., 90 if 3, 30 year runs have been averaged together.
    p_value = 0.05: default 95% T-Test
    sign = False: moderated and exacerbated calculated, True: also include over-effective classification.

    Outputs:
    result: array with same shape as input arrays containing results
    category_dict: dictionary with values used in results and corresponding categories
    """

    # calculate anomalies
    sg_anom = sg_mean - ctrl_mean
    cc_anom = cc_mean - ctrl_mean

    # caclculate where absolute anomaly increased / decreased
    mod = 1 * (abs(sg_anom) < abs(cc_anom)) # < produces a true/false result, multiplying by 1 gives 1s and 0s instead
    ex = 1 * (abs(sg_anom) > abs(cc_anom))

    # calculate where there has been a change in sign:
    cc_pos = 1 * (cc_anom > 0)
    sg_pos = 1 * (sg_anom > 0)

    # sign = sign-change
    # sig = passed significance T-Test.

    sg_sign_same = (cc_pos * sg_pos) + ((1 - cc_pos) * (1 - sg_pos)) # tests whether both pos or both neg. multiplying is same as AND, 1 - X is same as NOT.
    sg_sign_change = 1 - sg_sign_same 

    # call function which applies mod_ex t-test.
    mod_sig, ex_sig, dunno = better_worse_off(sg_mean, sg_std, cc_mean, cc_std, ctrl_mean, ctrl_std, years, 0.05)

    # find places moderated / exacerbated but not significant 
    mod_insig = mod - mod_sig
    ex_insig = ex - ex_sig
    
    category_dict = {-2:'exacerbated, significant',
                     -1:'exacerbated, insignificant',
                     1:'moderated, insignificant',
                     2:'moderated, significant',
    }
    
    # Return standard mod, ex with and without significance test.
    if not sign:
        
        category_dict = {-2:'exacerbated, significant',
                         -1:'exacerbated, insignificant',
                         1:'moderated, insignificant',
                         2:'moderated, significant',
                        }
        result = -2*ex_sig + -1*ex_insig + mod_insig + 2*mod_sig
        
        return result, category_dict

    if sign: # split all categories into two based on whether there is a change in sign.
        # exacerbated, significant, sign changed:
        ex_sig_sign = ex_sig * sg_sign_change
        ex_sig_nosign = ex_sig - ex_sig_sign
        # moderated, significant, sign changed:
        mod_sig_sign = mod_sig * sg_sign_change
        mod_sig_nosign = mod_sig - mod_sig_sign
        # exacerbated, insignificant, sign changed:
        ex_insig_sign = ex_insig * sg_sign_change
        ex_insig_nosign = ex_insig - ex_insig_sign
        # moderated, insignificant, sign changed:
        mod_insig_sign = mod_insig * sg_sign_change
        mod_insig_nosign = mod_insig - mod_insig_sign

        category_dict = {-2:'exacerbated, significant, no-sign', # 6%
                         -1:'exacerbated, insignificant, no-sign', # 5.1%
                         1:'moderated, insignificant, no-sign', # 4.9%
                         2:'moderated, significant, no-sign', # 4%
                         3:'moderated, significant, sign', # -4%
                         4:'moderated, insignificant, sign', # -4.9%
                         5:'exacerbated, insignificant, sign', # -5.1%
                         6:'exacerbated, significant, sign', # -6%
                        }
        result = (-2 * ex_sig_nosign + -1 * ex_insig_nosign + # exacerbated, same sign
                       1 * mod_insig_nosign + 2 * mod_sig_nosign + # moderated, without sign change
                       3 * mod_sig_sign + 4 * mod_insig_nosign + # moderated with sign change
                       5 * ex_insig_sign + 6 * ex_sig_sign) # sign change, insignificant magnitude change and exacerbated, changed sign.

        return result, category_dict

"""
A set of functions for sorting data and finding the quantiles of distributions.
"""

def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.

    http://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
    Examples:

    weighted_quantile([1, 2, 9, 3.2, 4], [0.0, 0.5, 1.])
    array([ 1. , 3.2, 9. ])

    weighted_quantile([1, 2, 9, 3.2, 4], [0.0, 0.5, 1.], sample_weight=[2, 1, 2, 4, 1])
    array([ 1. , 3.2, 9. ])
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

    # end weighted_quantiles

def fraction_distribution(data, values, cumulative=False, sample_weight=None):
    """
    This function takes a set of data and calculates what fraction falls around the values given,
    that is [<1st, >1st and <2nd, >2nd and <3rd, ..., >Nth]

    with input:
    data = [-1,4,5,9,15]
    values = [0,10]
    the function will return:
    Return = [0.2,0.6,0.2]

    with optional sample_weight:
    sample_weight = [0.5,0.5,1.0,1.0,2.0]
    Return = [0.1,0.5,0.4]

    With cumulative option:
    Return = [0.1,0.6,1.0]
    """

    # Copy and Flatten data
    flat_data = data.copy().flatten()

    # sort weighting
    if sample_weight is not None:
        # Copy and Flatten weight
        weight = sample_weight.copy().flatten()
        # Normalize weight
        weight = weight / weight.sum()
    else:
        weight = np.ones_like(flat_data)
        weight = weight / weight.sum()

    # Calculate
    lt_1st = np.sum((flat_data < values[0]) * weight)
    gt_last = np.sum((flat_data > values[-1]) * weight)

    if len(values) == 1:
        lt_gt = [lt_1st,gt_last]
    if len(values) > 1:
        lt_gt_rest = [ np.sum((flat_data < x) * (flat_data > y) * weight) for x, y in zip(values[1:], values[:-1]) ]
        lt_gt = [lt_1st]
        lt_gt.extend(lt_gt_rest)
        lt_gt.append(gt_last)

    if cumulative:
        temp = np.cumsum(lt_gt)
        return [x if x<1.0 else 1.0 for x in temp]
    else:
        return lt_gt

    #end fraction_distribution

def sort_data_distribution(data_in, data_sort, values, distribution=False, sort_weight=None, values_sorted=False):

    """
    Return = data_out{list of ndarrays}, fractions{list of floats}

    Sorts data_in by data_sort according to values. 2 modes of operation:

    Value-sort: data_sort is binned by values and the matching data_in grid_points are returned binned in
    the same way.
    e.g. values = [0,1] will return 3 ndarrays with values binned by data_sort into <0, 0>1, >1

    Distribution sort: data_sort is binned into percentiles (using a weighting) and this binning is applied
    to data_in.
    e.g. values = [0.25,0.5,0.75] will return 4 ndarrays with values binned by data_sort into quartiles

    If distribution = True, interpret values as percentiles, otherwise treat as values.
    
    If values_sorted = True, assume data already sorted, if False sort data and weighting.
    This allows point-wise comparisons and the recovery of a fully sorted 1D array.

    NOTES:

    Values must be in ascending order

    """
    
    if not values_sorted:
        sorter = np.argsort(data_sort)
        data_sort = data_sort[sorter]
        data_in = data_in[sorter]
        if sort_weight is not None:
            sort_weight = sort_weight[sorter]
    
    """
    Setup weighting if included:
    """
    if sort_weight is not None:
        # Copy and Flatten weight
        weight = sort_weight.copy().flatten()
        # Normalize weight
        weight = weight / weight.sum()
    else:
        weight = np.ones_like(data_sort)
        # Normalize weight
        weight = weight / weight.sum()

    """
    Distribution-sort Mode
    """
    if distribution:
        quantiles = weighted_quantile(data_sort, values, sample_weight=weight)
        values = quantiles

    """
    Value-Sort Mode:
    """
    #     if not distribution:

    # Record all sort values less than 1st value and greater than last
    lt_1st = data_sort < values[0]
    gt_last = data_sort > values[-1]

    if len(values) == 1:
        sort_list = [lt_1st,gt_last]
    elif len(values) > 1:
        # less than n, greater than n+1:
        lt_gt_rest = [ (data_sort < x) * (data_sort > y) for x, y in zip(values[1:], values[:-1]) ]
        # Combine all together into a list:
        sort_list = [lt_1st]
        sort_list.extend(lt_gt_rest) # extend as it's already a list
        sort_list.append(gt_last)
    else:
        return 'Values wrong length / type: ', values


    data_in_sorted = [data_in[X] for X in sort_list]
    data_frac = [np.sum(X * weight) for X in sort_list]

    return data_in_sorted, data_frac

    # end sort_data_distribution
