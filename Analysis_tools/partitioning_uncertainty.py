# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 10:41:38 2015

@author: imchugh
"""

# Python standard modules
import os
import numpy as np
import copy as cp
import random
import time

# My modules
import DataIO as io
import respiration as re
import random_error as ra
import datetime_functions as dtf
import data_formatting as dt_fm

reload(re)

# Get the data and format appropriately
def get_data(configs_dict):

    # Get file extension and target
    paths_dict = configs_dict['files']
    ext = os.path.splitext(paths_dict['input_file'])[1]
    data_input_target = os.path.join(paths_dict['input_path'],
                                     paths_dict['input_file'])

    # Initialise name change dictionary with new names via common keys
    oldNames_dict = configs_dict['variables']
    newNames_dict = {'carbon_flux':'Fc',
                     'temperature': 'TempC',
                     'solar_radiation': 'Fsd',
                     'vapour_pressure_deficit': 'VPD',
                     'friction_velocity': 'ustar',
                     'wind_speed': 'ws'}
    names_dict = {oldNames_dict[key]: newNames_dict[key] for key in oldNames_dict}                     

    # get data (screen only the Fc data to obs only)
    if ext == '.nc':
        Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                  var_list = [oldNames_dict
                                                              ['carbon_flux']],
                                                  QC_accept_codes = [0])
        Fc_dict.pop('date_time')
        ancillary_vars = [oldNames_dict[var] for var in oldNames_dict.keys() 
                          if not var == 'carbon_flux']
        ancillary_dict, global_attr = io.OzFluxQCnc_to_data_structure(
                                          data_input_target,
                                          var_list = ancillary_vars,
                                          return_global_attr = True)
        data_dict = dict(Fc_dict, **ancillary_dict)
    elif ext == '.df':
        data_dict, global_attr = io.DINGO_df_to_data_structure(
                                     data_input_target,
                                     var_list = oldNames_dict.values(),
                                     return_global_attr = True)

    # Add model data
    data_dict['Fc_model'] = (dt_fm.get_model_NEE_from_OzFluxQCncL6
                             (data_input_target)['Fc_model'])
    
    # Rename relevant variables    
    data_dict = dt_fm.rename_data_dict_vars(data_dict, names_dict)

    return data_dict, global_attr
    
#------------------------------------------------------------------------------    
# Do the standard respiration fit

# Get configurations
configs_dict = io.config_to_dict(io.file_select_dialog())

# Get data
data_dict, attr = get_data(configs_dict)

# Make a respiration config dict from config file and write measurement 
# interval to it
re_configs_dict = configs_dict['respiration_configs']
re_configs_dict['measurement_interval'] = int(attr['time_step'])

# Make an uncertainty configs dict from config file
uncert_configs_dict = configs_dict['partitioning_uncertainty']

# Local var names for config items
step = re_configs_dict['step_size_days']
window = re_configs_dict['window_size_days']
ustar_threshold = re_configs_dict['ustar_threshold']
num_trials = uncert_configs_dict['num_trials']
gaps = uncert_configs_dict['gaps']
gap_type = uncert_configs_dict['gap_type']
if not gap_type == 'obs':
    if not (isinstance(gap_type, float) or isinstance(gap_type, int)):
        print 'Variable gap_type needs to be a float or int between 1 and 100... ' \
              'reverting to observational gaps!'  
    gap_type = 'obs'

# Assign observational Fc to the 'Fc_series' var
data_dict['Fc_series'] = cp.copy(data_dict['Fc'])

# Get the datetime variable so can construct a new partitioned dataset later
datetime_array = data_dict['date_time']

# Get stats on data availability for each year
years_input_index_dict = dtf.get_year_indices(datetime_array)
obs_year_stats = {'n_cases_total': {},
                  'n_cases_avail': {}}
for year in years_input_index_dict.keys():
    indices = years_input_index_dict[year]
    this_dict = {var: data_dict[var][indices[0]: indices[1] + 1] 
                 for var in ['Fc_series', 'Fsd', 'ustar']}
    obs_year_stats['n_cases_total'][year] = indices[1] - indices[0]
    obs_year_stats['n_cases_avail'][year] = len(this_dict['Fc_series']
                                                [(this_dict['Fsd'] < 5) &
                                                 (this_dict['ustar'] > ustar_threshold) &
                                                 (~np.isnan(this_dict['Fc_series']))])

# Remove low ustar values according to threshold
data_dict['Fc_series'][data_dict['ustar'] < ustar_threshold] = np.nan

# Create a boolean array for data indicating presence of nans
nan_boolean = np.isnan(data_dict['Fc_series'])

# Calculate Re by sending data to main respiration function
re_dict, params_dict = re.main(data_dict, re_configs_dict)
data_dict['Re'] = re_dict['Re']

# Calculate sums for each year
sums_dict = {}
for year in years_input_index_dict.keys():
    indices = years_input_index_dict[year]
    sums_dict[year] = (data_dict['Re'][indices[0]: indices[1] + 1] 
                       * 12 * 0.0018).sum()                                                       
#------------------------------------------------------------------------------
                                                       
# Get the indices of the start and end rows of each unique date in the source 
# data array
dates_input_index_dict = dtf.get_day_indices(datetime_array)

# Get the indices of the start and end rows of each window in the source 
# data array
step_dates_input_index_dict = dtf.get_moving_window_indices(datetime_array, 
                                                            window, step)
                                                            
#------------------------------------------------------------------------------

# Get random error estimate using model data as input
ra_configs_dict = configs_dict['random_error_configs']
ra_configs_dict['measurement_interval'] = int(attr['time_step'])
ra_fig, ra_stats_dict = ra.regress_sigma_delta(data_dict, ra_configs_dict)
sigma_delta = ra.estimate_sigma_delta(data_dict['Re'], ra_stats_dict)

#------------------------------------------------------------------------------

# Initalise parameter dicts with prior estimates
all_noct_dict = re.filtering(data_dict)
params_in_dict = {'Eo_prior': 100,
                  'rb_prior': all_noct_dict['Fc_series'].mean()}

# Make results arrays
trial_sums_dict = {year: np.zeros([num_trials]) for year in 
                   years_input_index_dict.keys()}

#------------------------------------------------------------------------------
a = time.time()
# Loop through number of trials
for i in xrange(num_trials):
    
    # Generate a results dictionary for the parameter values (1 for each day)
    params_out_dict = re.generate_results_array(datetime_array)

    # Generate noise estimate
    noise = ra.estimate_random_error(sigma_delta)

    # Sum noise with model estimate of Re and assign to 'Fc_series' variable
    data_dict['Fc_series'] = data_dict['Re'] + noise

    # If requested, impose either the gaps observed in the observational data 
    # or introduce randomly chosen gaps as a percentage of model data
    if gaps:
        if gap_type == 'obs':
            data_dict['Fc_series'][nan_boolean] = np.nan
        else:
            for year in years_input_index_dict.keys():
                indices = years_input_index_dict[year]
                num = int(round(gap_type / 100.0 * (indices[1] - indices[0])))
                nan_index = random.sample(np.arange(indices[1] - indices[0] + 1), num)
                data_dict['Fc_series'][indices [0]: indices[1] + 1][nan_index] = np.nan
            
    # Partition data into year and step
    years_data_dict = re.segment_data(data_dict, years_input_index_dict)
    step_data_dict = re.segment_data(data_dict, step_dates_input_index_dict)

    # Get Eo for all years
    re.calculate_Eo(years_data_dict, 
                    re_configs_dict,
                    params_in_dict,
                    params_out_dict)
    
    # Get rb for all steps
    re.calculate_rb(step_data_dict,
                    re_configs_dict,
                    params_in_dict,
                    params_out_dict)
    
    # Estimate Re for all data
    this_dict = {'Re': re.estimate_Re(data_dict,
                                      params_out_dict,
                                      dates_input_index_dict)}

    # Calculate sums for each year                     
    for year in years_input_index_dict.keys():
        indices = years_input_index_dict[year]
        trial_sums_dict[year][i] = (this_dict['Re'][indices[0]: 
                                    indices[1] + 1] * 12 * 0.0018).sum()
                                        
b = time.time()

print b - a