# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:37:38 2015

@author: imchugh
"""
import pdb
import numpy as np
from scipy import stats

def estimate_model_error(data_dict, configs_dict):

    # Create data arrays, then count all
    obs_array = data_dict['NEE_series']
    mod_array = data_dict['NEE_model']
    total_records = len(obs_array)

    # Rescale to gC m-2
    for arr in [obs_array, mod_array]:
        arr[...] = arr * configs_dict['measurement_interval'] * 60 * 12 * 10 ** -6

    # Calculate annual sum for obs and model combined
    annual_sum = np.where(np.isnan(obs_array), mod_array, obs_array).sum()

    # Subset arrays to remove nans and calculate proportion remaining
    nan_index = ~np.isnan(obs_array)
    obs_array = obs_array[nan_index]
    mod_array = mod_array[nan_index]
    avail_records = len(obs_array)
    
    # Get the amount of data that will be removed from the sample (based on the 
    # proportion missing from the complete dataset)
    sample_missing = 1000 - int(1000 * (avail_records / float(total_records)))
    
    # Draw a random sample of 1000 data from the timeseries, then calculate the
    # difference between the observed and model-spliced series (appropriately 
    # scaled to gC m-2)
    random_index = np.random.randint(0, len(obs_array), 1000)
    subset_obs_array = obs_array[random_index]
    subset_mod_array = mod_array[random_index]
    subset_splice_array = np.concatenate([subset_mod_array[:sample_missing],
                                          subset_obs_array[sample_missing:]])
    obs_sum = subset_obs_array.sum()
    splice_sum = subset_splice_array.sum()
    return (obs_sum - splice_sum) / obs_sum * annual_sum

def propagate_model_error(data_dict, configs_dict):

    # Calculate critical t-statistic for p = 0.095
    crit_t = stats.t.isf(0.025, configs_dict['num_trials'])

    # Create arrayn to hold results then iterate over num_trials
    error_array = np.empty(configs_dict['num_trials'])        
    for this_trial in xrange(configs_dict['num_trials']):
        error_array[this_trial] = estimate_model_error(data_dict, configs_dict)
        
    return np.round(error_array.std() * crit_t, 2)    