# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 12:13:41 2015

@author: imchugh
"""

import numpy as np
import operator
import copy as cp
import warnings
import pdb

import datetime_functions as dtf
reload(dtf)
#------------------------------------------------------------------------------
# Numpy functions

def count_nans_in_array(arr):
    
    index = ~np.isnan(arr)
    start_len = len(arr)
    end_len = len(arr[index])
    
    return {'Total_obs': start_len,
            'Avail_obs': end_len,
            'Pct_avail_obs': round(end_len / float(start_len) * 100, 1)}

def set_numpy_array_to_nan(data_array, boolean_index):
    
    data_array[~boolean_index] = np.nan
    return data_array

def threshold_numpy_array(data_array, threshold, operation):
    """
    Creates a boolean index indicating whether values of numpy array are 
    greater than, less than, equal to or not equal to a given threshold
    Pass: 1) data_array - a numpy arrays of data
          2) threshold - numerical value of desired threshold
          3) operator - which data to keep (choices are <, >, ==, !=)
    Returns: numpy boolean array
    """
    ops = {">": operator.gt, "<": operator.lt, "==": operator.eq, 
           "!=": operator.ne}
    return ops[operation](data_array, threshold)

#------------------------------------------------------------------------------
# Basic Numpy array dictionary filtering functions 
#   Note that it is assumed that all dictionary
#   keys contain numpy arrays of equal length - indexing will most likely fail 
#   out of range if non-equal length arrays are contained in dict)

def set_numpy_dict_to_nan(data_dict, boolean_index):
    
    for var in data_dict.keys():
        data_dict[var] = set_numpy_array_to_nan(data_dict[var], [boolean_index])
    return data_dict

def subset_numpy_dict(data_dict, boolean_index):
    
    return {var: data_dict[var][boolean_index] for var in data_dict.keys()}

#------------------------------------------------------------------------------
# Higher level Numpy array dictionary filtering functions
#   Note that it is assumed that all dictionary
#   keys contain numpy arrays of equal length - indexing will most likely fail 
#   out of range if non-equal length arrays are contained in dict)

def sort_dict_on_index_variable(data_dict, sort_var):
    """
    Sorts a dictionary of equal-length numpy arrays on the basis of a
    sorting variable
    Pass: 1) data_dict - a dictionary containing arrays
          2) sort_var - the variable whose sorted order will dictate the 
                        ordering of all others (str)
    Returns: sorted dictionary
    """
    index = data_dict[sort_var].argsort()
    for key in data_dict.keys():
        data_dict[key] = data_dict[key][index]
    return data_dict

def subset_arraydict_on_nan(data_dict, var_list = False, condition = 'any', 
                            subset = True):
    """
    Removes all cases where either any or all variables (casewise) in array 
    are nan and return remaining data (ignores non-numeric dtypes)
    """    
    boolean_list = []
    these_vars = var_list if var_list else data_dict.keys()
    for var in these_vars:
        try:
            boolean_list.append(~np.isnan(data_dict[var]))
        except:
            continue
    if condition == 'any':
        all_boolean_index = [all(rec) for rec in zip(*boolean_list)]
    elif condition == 'all':
        all_boolean_index = [any(rec) for rec in zip(*boolean_list)]
    else:
        raise Exception('Valid keywords are "any" or "all"')    
    
    boolean_index = np.array(all_boolean_index)
    if subset:
        return subset_numpy_dict(data_dict, boolean_index)
    else:
        return boolean_index
    
def subset_arraydict_on_threshold(data_dict, threshold_var, threshold, 
                                  keep_cond, drop = False):
    """
    Pass: 1) data_dict - a dictionary containing numpy data arrays
          2) threshold_var - namestring of variable that is used for thresholding
          3) threshold - numerical value of desired threshold
          4) keep_cond - which data to keep (choices are <, >, ==, !=)
          5) drop - optional kwarg (default = False) specifying whether to drop
                    filtered data or set to nan
    Returns: filtered data dictionary
    """
    boolean_index = threshold_numpy_array(data_dict[threshold_var], threshold,
                                          keep_cond)
    if drop:
        return subset_numpy_dict(data_dict, boolean_index)
    else:
        return set_numpy_dict_to_nan(data_dict, boolean_index)
    
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
    
def IQR_filter(data_array, outlier_value = 1.5, minimum_data_avail = 50,
               inplace = True):
    
    valid_data_array = data_array[~np.isnan(data_array)]
    
    if not len(valid_data_array) / float(len(data_array)) * 100 > minimum_data_avail:
        print 'Percentage of valid data below minimum threshold - returning...'
        return

    lo_qtl = np.percentile(valid_data_array, 25)
    hi_qtl = np.percentile(valid_data_array, 75)
    qtl_range = hi_qtl - lo_qtl
    lo_threshold = lo_qtl - outlier_value * qtl_range
    hi_threshold = hi_qtl + outlier_value * qtl_range
    print 'Lower threshold is ' + str(lo_threshold)
    print 'Upper threshold is ' + str(hi_threshold)
    lo_bool_array = data_array < lo_threshold
    hi_bool_array = data_array > hi_threshold
    all_bool_array = lo_bool_array | hi_bool_array
    if not inplace:
        new_array = cp.copy(data_array)
    else:
        new_array = data_array
    new_array[all_bool_array] = np.nan
    
    return new_array

def slide_IQR_filter(data_array, outlier_value = 2, window_size = 11,
                     inplace = True):

    if window_size > len(data_array):
        raise Exception('Window size cannot exceed array size! Quitting...')

    if window_size == len(data_array):
        iter_array = np.array([0])
    else:
        if window_size % 2 == 0:
            window_size = window_size + 1
            if window_size == len(data_array):
                iter_array = np.array([0])
        else:
            iter_array = np.arange(0, len(data_array) - window_size)

    # Create lower and upper threshold arrays
    lo_threshold_array = np.empty(len(data_array))
    lo_threshold_array[:] = np.nan
    hi_threshold_array = np.empty(len(data_array))
    hi_threshold_array[:] = np.nan

    # Calculate outliers for each window
    rslt_index_int = int(window_size / 2)
    for i in iter_array:
        this_array = data_array[i: i + window_size]
        lo_qtl = np.percentile(this_array, 25)
        hi_qtl = np.percentile(this_array, 75)
        qtl_range = hi_qtl - lo_qtl
        lo_threshold_array[i + rslt_index_int] = lo_qtl - outlier_value * qtl_range
        hi_threshold_array[i + rslt_index_int] = hi_qtl + outlier_value * qtl_range

    # Fill gaps using median of low and high thresholds over whole array
    lo_median = np.median(lo_threshold_array[~np.isnan(lo_threshold_array)])
    lo_threshold_array[np.isnan(lo_threshold_array)] = lo_median
    hi_median = np.median(hi_threshold_array[~np.isnan(hi_threshold_array)])
    hi_threshold_array[np.isnan(hi_threshold_array)] = hi_median    

    # Create boolean array (where values outside range will be True)
    bool_array = ((data_array < lo_threshold_array) |
                  (data_array > hi_threshold_array))
    
    
    # Filter original series (or copy if requested)
    if not inplace:
        new_array = cp.copy(data_array)
        new_array[bool_array] = np.nan
        return new_array
    else:
        data_array[bool_array] = np.nan
        return
        
# Remove low ustar values
def screen_low_ustar(data_dict, ustar_threshold, noct_threshold):

    if isinstance(ustar_threshold, dict):
        years_data_dict = dtf.subset_datayear_from_arraydict(data_dict, 'date_time')
        threshold_keys = [int(key) for key in ustar_threshold.keys()]
        miss_list = [year for year in years_data_dict.keys() 
                     if not year in threshold_keys]
        if not len(miss_list) == 0:
            miss_string = ', '.join([str(this_year) for this_year in miss_list])
            raise Exception('Missing years: %s' %miss_string + '; please edit ' \
                            'your configuration file so that years specified ' \
                            'for ustar threshold match those available in ' \
                            'data file; alternatively, specify a single ' \
                            'value (float or int) for ustar. Exiting...')
        extra_list = [year for year in threshold_keys 
                      if not year in years_data_dict.keys()]
        if not len(extra_list) == 0:
            extra_string = ', '.join([str(this_year) for this_year in miss_list])
            warnings.warn('Years %s' %extra_string + ' specified in ustar ' \
                          'dictionary are not present in passed data dictionary' \
                          'and have been ignored! Continuing...')
        data_list = []
        for this_year in years_data_dict.keys():
            this_threshold = ustar_threshold[str(this_year)]
            this_NEE = years_data_dict[this_year]['NEE_series']
            this_NEE[(years_data_dict[this_year]['ustar'] < this_threshold) &
                     (years_data_dict[this_year]['Fsd'] < noct_threshold)] = np.nan
            data_list.append(this_NEE)
        data_dict['NEE_series'] = np.concatenate(data_list)
    elif isinstance(ustar_threshold, (float, int)):
        data_dict['NEE_series'][(data_dict['ustar'] < ustar_threshold) &
                                (data_dict['Fsd'] < noct_threshold)] = np.nan
    else:
        raise TypeError('ustar_threshold object must be dtype float or dict!' \
                        'Please edit your configuration file')

    return

def check_missing_data(data_dict, var_list = False):
    
    """
    Finds missing data (and reports total cross-record percentage) in one or more 
    arrays of data
    Pass: 1) a data dictionary containing the relevant arrays
          2) (optional) a variable list specifying the variables to be checked
             (defaults to all data arrays in the dictionary if var_list is not
              passed)
    Returns: none - output is printed to screen.
    """

    if not var_list: var_list = data_dict.keys()
    miss_list = []
    for var in var_list:
        miss_list.append(subset_arraydict_on_nan(data_dict, 
                                                 var_list = [var],
                                                 condition = 'any', 
                                                 subset = False))
        any_miss_list = subset_arraydict_on_nan(data_dict, 
                                                var_list = var_list,
                                                condition = 'any', 
                                                subset = False)
    if not all(any_miss_list):
        any_miss_array = np.array(any_miss_list)
        missing_drivers_pct = round(len(any_miss_array[~any_miss_array]) / 
                                    float(len(any_miss_array)) * 100, 1)
        print 'Warning: ' + str(missing_drivers_pct) + '% of all available ' \
              'records are missing data for at least one variable; individual ' \
              'variable breakdown follows: '
        for i, var in enumerate(var_list):
            this_array = np.array(miss_list[i])
            count = len(this_array[~this_array])
            if not count == 0:
                print '    - ' + var + ' is missing ' + str(count) + ' records;' 
    
    return
    