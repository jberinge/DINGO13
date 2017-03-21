# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:23:04 2015

@author: imchugh
"""
# Python modules
import os
import numpy as np
import copy as cp
import calendar
import matplotlib.pyplot as plt
import datetime as dt
from scipy.optimize import curve_fit
import pdb

# My modules
import datetime_functions as dtf
import data_filtering as filt
import gap_filling as gf


#------------------------------------------------------------------------------
# Data optimisation algorithm
def TRF(data_dict, Eo, rb):
    
    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_dict['TempC'] + 46.02)))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Write error messages to dictionary with codes as keys
def error_codes():
    
    d = {0:'Optimisation successful',
         1:'Value of Eo failed range check - rejecting all parameters',
         2:'Value of rb has wrong sign - rejecting all parameters',
         3:'Optimisation reached maximum number of iterations ' \
           'without convergence',
         10:'Data did not pass minimum percentage threshold - ' \
            'skipping optimisation'}
    
    return d
#------------------------------------------------------------------------------    

#------------------------------------------------------------------------------
# run optimisation and raise error code for combined rb and Eo
def optimise_all(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    drivers_dict = {driver: data_dict[driver] for driver in ['TempC']}
    response_array = data_dict['NEE_series']

    try:
        params = curve_fit(lambda x, a, b:
                           TRF(x, a, b),
                           drivers_dict, 
                           response_array, 
                           p0 = [params_dict['Eo_prior'], 
                                 params_dict['rb_prior']])[0]
    except RuntimeError:
        params = [np.nan, np.nan]
        error_state = 3

    # If negative rb returned, set to nan
    if params[0] < 50 or params[0] > 400: 
        error_state = 1
        params = [np.nan, np.nan]
    elif params[1] < 0:
        error_state = 2
        params = [np.nan, np.nan]

    return {'Eo': params[0], 'rb': params[1], 'error_code': error_state}
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------    
# run optimisation and raise error code for rb with fixed Eo
def optimise_rb(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    
    # Get drivers and response
    drivers_dict = {driver: data_dict[driver] for driver in ['TempC']}
    response_array = data_dict['NEE_series']        
    
    try:
        params = curve_fit(lambda x, b:
                           TRF(x, params_dict['Eo_default'], b),
                           drivers_dict, 
                           response_array, 
                           p0 = [params_dict['rb_prior']])[0]                           
    except RuntimeError:
        params = [np.nan]

    # If negative rb returned, set to nan
    if params[0] < 0:
        error_state = 2
        params = [np.nan]
       
    return {'rb': params[0], 'error_code': error_state}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Iterate through dates and write optimisation results for rb parameter to
# parameter results dictionary
def calculate_rb(data_dict,
                 configs_dict,
                 params_in_dict,
                 params_out_dict):

    # Create local vars from configs
    meas_int = configs_dict['measurement_interval']    
    min_pct = configs_dict['minimum_pct_window']    
    
    # Calculate rb for windows    
    for date in data_dict.keys():

        # Make a temporary dict with nans and daytime dropped
        bool_filter = data_dict[date]['all_bool']
        temp_dict = {var: data_dict[date][var][bool_filter] 
                     for var in ['TempC', 'NEE_series']}
        
        # Get the index for the current date
        date_index = np.where(params_out_dict['date'] == date)
        
        # Check data availability 
        data_pct = int(len(temp_dict['NEE_series']) / 
                       float((1440 / meas_int) / 2) * 100)
        
        # If enough data, go ahead
        if not data_pct < min_pct:

            # Do fit
            params_in_dict['Eo_default'] = params_in_dict['Eo_years'][date.year]            
            fit_dict = optimise_rb(temp_dict, params_in_dict)
            fit_dict['rb_error_code'] = fit_dict.pop('error_code')

            # Write data to results arrays
            for key in fit_dict.keys():
                params_out_dict[key][date_index] = fit_dict[key]
            
        else:

            # Error code for not enough data
            params_out_dict['rb_error_code'][date_index] = 10
        
    # Interpolate rb
    params_out_dict['rb'] = gf.generic_2d_linear(params_out_dict['rb'])

    return
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Iterate through years and write optimisation results for Eo parameter to
# parameter results dictionary
def calculate_Eo(data_dict,
                 configs_dict,
                 params_in_dict,
                 params_out_dict):

    # Create local vars from configs
    meas_int = configs_dict['measurement_interval']
    min_pct = configs_dict['minimum_pct_annual']

    # Do annual fits for Eo
    Eo_annual_data_dict = {}
    Eo_annual_error_dict = {}
    Eo_pass_keys = []
    Eo_fail_keys = []
    for year in data_dict.keys():

        # Make a temporary dict with nans and daytime dropped
        bool_filter = data_dict[year]['all_bool']
        temp_dict = {var: data_dict[year][var][bool_filter] 
                     for var in ['TempC', 'NEE_series']}

        # Calculate number of nocturnal recs for year
        days = 366 if calendar.isleap(year) else 365
        recs = days * (1440 / meas_int) / 2

        # Input indices
        data_pct = int(len(temp_dict['NEE_series']) / float(recs) * 100)
        if not data_pct < min_pct:
            fit_dict = optimise_all(temp_dict, params_in_dict)
        else:
            fit_dict = {'Eo': np.nan, 'error_code' : 10}
        Eo_annual_data_dict[year] = fit_dict['Eo']
        Eo_annual_error_dict[year] = fit_dict['error_code']
        if fit_dict['error_code'] == 0: 
            Eo_pass_keys.append(year)
        else:
            Eo_fail_keys.append(year)
            
    # Fill gaps 
    if np.all(np.isnan(Eo_annual_data_dict.values())):
        raise Exception('Could not find any values of Eo for selected year(s)! Exiting...')
    if np.any(np.isnan(Eo_annual_data_dict.values())):
        Eo_mean = np.array([Eo_annual_data_dict[year] for year in Eo_pass_keys]).mean()    
        for year in Eo_fail_keys:
            Eo_annual_data_dict[year] = Eo_mean
            Eo_pass_keys.append(year)

    # Attach the yearly Eo to the parameter dictionary
    params_in_dict['Eo_years'] = Eo_annual_data_dict

    # Project Eo to the appropriate indices of the results array
    years_array = np.array([date.year for date in params_out_dict['date']])
    for year in Eo_pass_keys:
        index = np.where(years_array == year)[0]
        out_indices = [index[0], index[-1]]
        params_out_dict['Eo'][out_indices[0]: out_indices[1] + 1] = (
            Eo_annual_data_dict[year])
        params_out_dict['Eo_error_code'][out_indices[0]: out_indices[1] + 1] = (
            Eo_annual_error_dict[year])

    return
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Use observed meteorology and optimised respiration parameters to estimate Re
def estimate_Re(data_dict, 
                all_params_dict, 
                datetime_input_index_dict):
    
    # Create output dicts for estimates
    results_dict = {}
    results_dict['Re'] = np.empty(len(data_dict['TempC']))
    results_dict['Re'][:] = np.nan
    results_dict['date_time'] = data_dict['date_time']
    
    # Estimate time series Re
    for i, date in enumerate(all_params_dict['date']):
        
        indices = datetime_input_index_dict[date]
        this_Eo = all_params_dict['Eo'][i]
        this_rb = all_params_dict['rb'][i]
        this_dict = {'TempC': data_dict['TempC'][indices[0]: indices[1] + 1]}
        Re = TRF(this_dict, this_Eo, this_rb)                                                         
        results_dict['Re'][indices[0]: indices[1] + 1] = Re

    return results_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Build the dictionary for 
def generate_results_dict(datetime_array):
    
    dates_input_index_dict = dtf.get_day_indices(datetime_array)
    date_array = np.array(dates_input_index_dict.keys())
    date_array.sort()
    generic_array = np.empty(len(dates_input_index_dict))
    generic_array[:] = np.nan
    rb_error_code_array = np.ones(len(generic_array)) * 20
    return {'date': date_array,
            'Eo': cp.copy(generic_array),
            'Eo_error_code': cp.copy(generic_array),
            'rb': cp.copy(generic_array),
            'rb_error_code': rb_error_code_array}  
#------------------------------------------------------------------------------
            
#------------------------------------------------------------------------------            
# Nocturnal fits for each window
def plot_windows(step_data_dict, configs_dict, params_dict):

    # Don't send plots to screen
    if plt.isinteractive():
        is_on = True
        plt.ioff()
    else:
        is_on = False

    # Set parameters from dicts
    window = configs_dict['window_size_days']
    
    x_lab = '$Temperature\/(^{o}C$)'
    
    for date in step_data_dict.keys():

        Eo = params_dict['Eo'][params_dict['date'] == date]
        rb = params_dict['rb'][params_dict['date'] == date]

        bool_filter = step_data_dict[date]['night_bool']
        x_var = step_data_dict[date]['TempC'][bool_filter]
        y_var1 = step_data_dict[date]['NEE_series'][bool_filter]
        index = x_var.argsort()
        x_var = x_var[index]
        y_var1 = y_var1[index]
        y_var2 = TRF({'TempC': x_var}, Eo, rb)
          
        # Plot
        date_str = dt.datetime.strftime(date,'%Y-%m-%d')
        fig = plt.figure(figsize = (12,8))
        fig.patch.set_facecolor('white')
        ax = plt.gca()
        ax.plot(x_var, y_var1, 'o' , markerfacecolor = 'none',
                 markeredgecolor = 'black', label = 'NEE_obs', color = 'black')
        ax.plot(x_var, y_var2, linestyle = ':', color = 'black', 
                 label = 'NEE_est')
        ax.set_title('Fit for ' + str(window) + ' day window centred on ' + 
                      date_str + '\n', fontsize = 22)
        ax.set_xlabel(x_lab, fontsize = 18)
        ax.set_ylabel('$NEE\/(\mu mol C\/m^{-2} s^{-1}$)', fontsize = 18)
        ax.axhline(y = 0, color = 'black')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plot_out_name = 'noct' + '_' + date_str + '.jpg'
        plt.tight_layout()
        fig.savefig(os.path.join(configs_dict['output_path'],
                                 plot_out_name))
        plt.close(fig)
    
    if is_on:
        plt.ion()
        
    return
#------------------------------------------------------------------------------
                
#------------------------------------------------------------------------------
#############
# Main code #
#############        
def main(data_dict, configs_dict):

    """
    Calculates Re using Lloyd and Taylor function where Eo is fitted to annual
    data and rb is fitted to specified combination of step and window;
    
    Pass: 1) a data dictionary containing the following key / value pairs (it 
             is assumed that there are no gaps in the time series, but this is 
             not yet enforced!; also, all values in dict must be numpy arrays, 
             and all must be of same length):
                 - 'date_time': numpy array of Python datetimes
                 - 'NEE_series': numpy array of the NEE time series to be used 
                                as the optimisation target (float); note:
                                - missing values must be np.nan
                                - no QC is done on values - this must be done 
                                  prior to passing the data
                 - 'TempC': numpy array of temperatures to be used as 
                            optimisation input (float)
                            
          2) a configs dict containing the following key / value pairs:
                 - 'step_size_days': step size in days between fitting windows
                                     (int; range 0 < x < n days)
                 - 'window_size_days': width of fitting window in days 
                                     (int; range 0 < x < n days)
                 - 'minimum_pct_annual': minimum acceptable percentage of 
                                         available annual data for fitting of 
                                         Eo (int; range 0 <= x <= 100)
                 - 'minimum_pct_window': minimum acceptable percentage of 
                                         available window data for fitting of 
                                         rb (int; range 0 <= x <= 100)
                 - 'measurement_interval': measurement interval (minutes) of 
                                           the input data (int)
                 - 'output_fit_plots': whether to output plots showing the fit
                                       of the parameters to the data (boolean)
                 - 'output_path': path for output of plots (str)
                 
    Returns: 1) a results dictionary containing 2 key / value pairs:
                    - 'date_time': numpy array of Python datetimes for each
                      datum in the original time series
                    - 'Re': numpy array of half-hourly estimates of Re
                    
             2) a results dictionary containing 5 key / value pairs:
                    - 'date_time': numpy array of Python dates for each day in 
                                   the original time series 
                    - 'Eo': numpy array of Eo estimates for each day (note that 
                            each year has a constant value)
                    - 'Eo_error_code': numpy array of Eo diagnostic errors for 
                                       each day (note that each year has a 
                                       constant value)
                    - 'rb': numpy array of rb estimates for each day (note that 
                            linear interpolation is used to gap fill between 
                            steps)
                    - 'rb_error_code': numpy array of rb diagnostic errors for 
                                       each day (including whether the value is 
                                       interpolated or calculated)                    
    """
    
    # Create boolean indices for masking daytime and nan values
    night_mask = data_dict['Fsd'] < 5
    nan_mask = filt.subset_arraydict_on_nan(data_dict,
                                            var_list = ['NEE_series', 'TempC'],
                                            subset = False)
    all_mask = [all(rec) for rec in zip(night_mask, nan_mask)]
    data_dict['night_bool'] = np.array(night_mask)
    data_dict['all_bool'] = np.array(all_mask)

    # Partition the data into year and step pieces
    years_data_dict = dtf.get_year_window(data_dict,
                                          'date_time')

    step_data_dict = dtf.get_moving_window(data_dict, 
                                           'date_time', 
                                           configs_dict['window_size_days'],
                                           configs_dict['step_size_days'])

    # Get the indices of the start and end rows of each unique date in the source 
    # data array - no data dict is built from this, since these indices are used to
    # assign Re estimates to the estimated time series output array only
    dates_input_index_dict = dtf.get_day_indices(data_dict['date_time'])

    # Initalise parameter dicts with prior estimates
    params_in_dict = {'Eo_prior': 100,
                      'rb_prior': data_dict['NEE_series'][data_dict['all_bool']]
                      .mean()}
    
    # Generate a results dictionary for the parameter values (1 for each day)
    params_out_dict = generate_results_dict(data_dict['date_time'])

    # Get Eo for all years
    calculate_Eo(years_data_dict, 
                 configs_dict,
                 params_in_dict,
                 params_out_dict)
    
    # Get rb for all steps
    calculate_rb(step_data_dict,
                 configs_dict,
                 params_in_dict,
                 params_out_dict)
    
    # Estimate Re for all data
    rslt_dict = estimate_Re(data_dict,
                            params_out_dict,
                            dates_input_index_dict)

    # Get error codes
    error_dict = error_codes()
    error_dict[20] = 'Data are linearly interpolated from nearest '\
                     'non-interpolated neighbour'

    # Do plotting if specified
    if configs_dict['output_fit_plots']:
        plot_windows(step_data_dict, configs_dict, params_out_dict)

    return rslt_dict, params_out_dict, error_dict
#------------------------------------------------------------------------------