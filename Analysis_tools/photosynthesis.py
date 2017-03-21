# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:23:04 2015

@author: imchugh
"""
# Python modules
import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from scipy.optimize import curve_fit
import pdb

# My modules
import datetime_functions as dtf
import data_filtering as filt
import gap_filling as gf

reload(filt)

#------------------------------------------------------------------------------
# Data optimisation algorithm

def LRF_part(data_d, Eo, rb, alpha, beta, k):
    beta_VPD = beta * np.exp(-k * (data_d['VPD'] - 1))
    index = data_d['VPD'] <= 1
    beta_VPD[index] = beta
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
           (alpha * data_d['PAR'] / beta_VPD))
#    GPP = (alpha * beta_VPD * data_d['PAR']) / (alpha * data_d['PAR'] + beta_VPD)    
    index = data_d['PAR'] < 5
    GPP[index] = 0
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
    return GPP, Reco

def LRF(data_d, Eo, rb, alpha, beta, k):
    beta_VPD = beta * np.exp(-k * (data_d['VPD'] - 1))
    index = data_d['VPD'] <= 1
    beta_VPD[index] = beta
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
           (alpha * data_d['PAR'] / beta_VPD))    
#    GPP = (alpha * beta_VPD * data_d['PAR']) / (alpha * data_d['PAR'] + beta_VPD)     
    index = data_d['PAR'] < 5
    GPP[index] = 0
    Reco = rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_d['TempC'] + 46.02)))
    return GPP + Reco
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
# Write error messages to dictionary with codes as keys
def error_codes():
    
    d = {0:'Optimisation successful',
         1:'Value of k failed range check - setting to zero and ' \
           'recalculating other parameters',
         2:'Value of alpha failed range check - using previous ' \
           'estimate (zero if unavailable) and recalculating other ' \
           'parameters',
         3:'Value of beta failed range check - rejecting all parameters',
         4:'Value of daytime rb has wrong sign - rejecting all parameters',
         5:'Optimisation reached maximum number of iterations' \
           'without convergence',
         10:'Data did not pass minimum percentage threshold - ' \
            'skipping optimisation'}
    
    return d
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Daytime rb, alpha, beta, k

def optimise_free_rb(data_dict, params_dict):
    """
    This script simultaneously finds optimal parameter values of: 
        i)  a rectangular hyperbolic light response function of the form:
            GPP = (alpha * PAR) / (1 - (PAR / 2000) + (alpha * PAR / beta))
            where alpha is the initial slope of the light response curve and beta is
            the magnitude of GPP at 2000umol photons m^-2 s^-1
        ii) the reference respiration parameter of the lloyd and Taylor function:
            Re = rb * e^(Eo * (1 / (10 + 46.02) - 1 / (T + 46.02)))
    Positional arguments to pass in are: 
        1) a dictionary of numpy arrays (data_dict) with the following keys:
               - 'NEE_series' - carbon flux in umol m^-2 s^-1
               - 'TempC' - temperature in celsius
               - 'PAR' - photosynthetically active radiation in umol photons m^-2 s^-1 
               - 'VPD' - vapour pressure deficit in kPa 
        2) a dictionary containing required configurations (configs_dict), with the 
           following keys:
               'Eo_default' - used to fix value of Eo
               'alpha_default' - used if optimisation with free alpha fails
               'k_default' - used if optimisation with free k fails
               'rb_prior' - initial guess for rb
               'alpha_prior' - initial guess for alpha
               'beta_prior' - initial guess for beta
               'k_prior' - intial guess for k
    Note - no filtering is done here - any missing data values or nans will 
           propagate!!!               
    """

    # Initialise error state variable
    error_state = 0        

    # Get drivers and response
    drivers_d = {driver: data_dict[driver] for driver in ['PAR', 'TempC', 'VPD']}
    response_array = data_dict['NEE_series']      
    
    # Prepare the beta variations for iteration
    beta_dict = {'beta': params_dict['beta_prior'],
                 'beta_half': 0.5 * params_dict['beta_prior'],
                 'beta_double': 2 * params_dict['beta_prior']}
    beta_params_dict = {}
    beta_RMSE_dict = {}
    beta_error_dict = {}

    # Go!
    for key in beta_dict.keys():

        beta_dynamic = beta_dict[key]

        # Do the initial fit with largest number of parameters
        try:
            params = curve_fit(lambda x, b, c, d, e: 
                               LRF(x,
                                   params_dict['Eo_default'], 
                                   b, c, d, e),
                               drivers_d, 
                               response_array, 
                               p0 = [params_dict['rb_prior'], 
                                     params_dict['alpha_prior'], 
                                     beta_dynamic, 
                                     params_dict['k_prior']])[0] 
        except RuntimeError:
            params = [np.nan, np.nan, np.nan, np.nan]
        rb, alpha, beta, k = params[0], params[1], params[2], params[3]

        # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
        if ((np.isnan(k)) or (k < 0) or (k > 2)):
            error_state = 1            
            k = params_dict['k_default']
            try:
                params = curve_fit(lambda x, b, c, d:
                                   LRF(x,
                                       params_dict['Eo_default'],
                                       b, c, d,  
                                       params_dict['k_default']), 
                                   drivers_d, 
                                   response_array, 
                                   p0 = [params_dict['rb_prior'], 
                                         params_dict['alpha_prior'], 
                                         beta_dynamic])[0]
            except RuntimeError:
                params = [np.nan, np.nan, np.nan]
            rb, alpha, beta = params[0], params[1], params[2]
    
        # If a nan, positive or otherwise out of range value of alpha was returned,
        # rerun with previous value of alpha if available, otherwise set to zero
        if ((np.isnan(alpha)) | (alpha > 0) | (alpha < -0.22)):
            error_state = 2   
            k = params_dict['k_default']
            alpha = params_dict['alpha_default']
            try:            
                params = curve_fit(lambda x, b, d:
                                   LRF(x, 
                                       params_dict['Eo_default'],
                                       b,
                                       params_dict['alpha_default'],
                                       d, 
                                       params_dict['k_default']), 
                                   drivers_d, 
                                   response_array, 
                                   p0 = [params_dict['rb_prior'], 
                                         beta_dynamic])[0]
            except RuntimeError:
                params = [np.nan, np.nan]
            rb, beta = params[0], params[1]

        # If an out of range value of beta or rb was returned, then fagedaboutit!
        if ((beta > 0) | (beta < -100)):
            error_state = 3
            rb, alpha, beta, k = np.nan, np.nan, np.nan, np.nan
        elif rb < 0:
            error_state = 4
            rb, alpha, beta, k = np.nan, np.nan, np.nan, np.nan
        elif np.isnan(beta):
            error_state = 5            
            rb, alpha, beta, k = np.nan, np.nan, np.nan, np.nan

        # If valid beta was obtained, calculate RMSE and save error state and params
        beta_error_dict[key] = error_state
        beta_params_dict[key] = [rb, alpha, beta, k]
        if not np.isnan(beta):
            NEE_est = LRF(drivers_d,
                          params_dict['Eo_default'], 
                          rb, alpha, beta, k)
            beta_RMSE_dict[key] = np.sqrt(((response_array - 
                                            NEE_est[0])**2).mean())
        else:
            beta_RMSE_dict[key] = np.nan
    
    # If the beta RMSE dictionary is not filled with nans, find the lowest 
    # value, get its key and then get the corresponding parameters; otherwise
    # just return the parameter set acquired using the original beta estimate
    if not np.all(np.isnan(beta_RMSE_dict.values())):
        min_RMSE_array = np.array(beta_RMSE_dict.values())
        min_RMSE = min_RMSE_array[~np.isnan(min_RMSE_array)].min()
        for key in beta_RMSE_dict:
            if beta_RMSE_dict[key] == min_RMSE:
                params = beta_params_dict[key]
                error_state = beta_error_dict[key]
                break
    else:
        params = beta_params_dict['beta']
        error_state = beta_error_dict['beta']

    params_dict = {'rb': params[0],
                   'alpha': params[1],
                   'beta': params[2],
                   'k': params[3],
                   'error_code': error_state}

    return params_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def optimise_fixed_rb(data_dict, params_dict): 
    """
    This script finds optimal parameter values of a rectangular hyperbolic light 
    response function of the form:
        GPP = (alpha * PAR) / (1 - (PAR / 2000) + (alpha * PAR / beta))
        where alpha is the initial slope of the light response curve and beta is
        the magnitude of GPP at 2000umol photons m^-2 s^-1
    Positional arguments to pass in are: 
        1) a dictionary of numpy arrays (data_dict) with the following keys:
               - 'NEE_series' - carbon flux in umol m^-2 s^-1
               - 'TempC' - temperature in celsius
               - 'PAR' - photosynthetically active radiation in umol photons m^-2 s^-1 
               - 'VPD' - vapour pressure deficit in kPa 
        2) a dictionary containing required configurations (configs_dict), with the 
           following keys:
               'Eo_default' - used to fix value of Eo
               'alpha_default' - used if optimisation with free alpha fails
               'k_default' - used if optimisation with free k fails
               'alpha_prior' - initial guess for alpha
               'beta_prior' - initial guess for beta
               'k_prior' - initial guess for k
    Note - no filtering is done here - any missing data values or nans will 
           propagate!!!
    """
    
    # Initialise error state variable
    error_state = 0        

    # Get drivers and response
    drivers_d = {driver: data_dict[driver] for driver in ['PAR', 'TempC', 'VPD']}
    response_array = data_dict['NEE_series']      

    # Prepare the beta variations for iteration
    beta_dict = {'beta': params_dict['beta_prior'],
                 'beta_half': 0.5 * params_dict['beta_prior'], 
                 'beta_double': 2 * params_dict['beta_prior']}
    beta_params_dict = {}
    beta_RMSE_dict = {}
    beta_error_dict = {}

    for key in beta_dict.keys():

        # Set beta
        beta_dynamic = beta_dict[key]

        # Do the initial fit with largest number of parameters        
        try:
            params = curve_fit(lambda x, c, d, e:
                               LRF(x,
                                   params_dict['Eo_default'], 
                                   params_dict['rb_default'],
                                   c, d, e), 
                               drivers_d, 
                               response_array, 
                               p0 = [params_dict['alpha_prior'], 
                                     beta_dynamic, 
                                     params_dict['k_prior']])[0]
        except RuntimeError:
            params = [np.nan, np.nan, np.nan]
        alpha, beta, k = params[0], params[1], params[2]
               
        # If nan returned or a negative VPD coefficient, rerun with VPD set to zero
        if ((k == np.nan) | (k < 0) | (k > 2)):
            error_state = 1            
            k = params_dict['k_default']
            try:
                params = curve_fit(lambda x, c, d: 
                                   LRF(x,
                                       params_dict['Eo_default'], 
                                       params_dict['rb_default'],
                                       c, d, 
                                       params_dict['k_default']), 
                                   drivers_d, 
                                   response_array, 
                                   p0 = [params_dict['alpha_prior'], 
                                         beta_dynamic])[0]
            except RuntimeError:
                params = [np.nan, np.nan]
            alpha, beta = params[0], params[1]
    
        # If a nan, positive or otherwise out of range value of alpha was returned,
        # rerun with previous value of alpha if available, otherwise set to zero
        if ((alpha == np.nan) | (alpha > 0) | (alpha < -0.22)):
            error_state = 2   
            k = params_dict['k_default']
            alpha = params_dict['alpha_default']
            try:            
                params = curve_fit(lambda x, d:
                                   LRF(x,
                                       params_dict['Eo_default'],
                                       params_dict['rb_default'],
                                       params_dict['alpha_default'],
                                       d, 
                                       params_dict['k_default']),
                                   drivers_d, 
                                   response_array, 
                                   p0 = [beta_dynamic])[0]
            except RuntimeError:
                error_state = 3
                params = [np.nan]
            beta = params[0]

        # If an out of range value of beta was returned, then fagedaboutit!
        if ((beta > 0) | (beta < -100)):
            error_state = 3
            alpha, beta, k = np.nan, np.nan, np.nan
        elif np.isnan(beta):
            error_state = 5            
            alpha, beta, k = np.nan, np.nan, np.nan
       
        # If valid beta was obtained, calculate RMSE
        beta_error_dict[key] = error_state
        beta_params_dict[key] = [alpha, beta, k]
        if not np.isnan(beta):
            beta_params_dict[key] = [alpha, beta, k]
            NEE_est = LRF(drivers_d, 
                          params_dict['Eo_default'], params_dict['rb_default'], 
                          alpha, beta, k)
            beta_RMSE_dict[key] = np.sqrt(((response_array - 
                                            NEE_est[0])**2).mean())

    # If the beta RMSE dictionary is not filled with nans, find the lowest 
    # value, get its key and then get the corresponding parameters; otherwise
    # just return the parameter set acquired using the original beta estimate
    if not np.all(np.isnan(beta_RMSE_dict.values())):
        min_RMSE_array = np.array(beta_RMSE_dict.values())
        min_RMSE = min_RMSE_array[~np.isnan(min_RMSE_array)].min()
        for key in beta_RMSE_dict:
            if beta_RMSE_dict[key] == min_RMSE:
                params = beta_params_dict[key]
                error_state = beta_error_dict[key]
                break
    else:
        params = beta_params_dict['beta']
        error_state = beta_error_dict['beta']

    params_dict = {'alpha': params[0],
                   'beta': params[1],
                   'k': params[2],
                   'error_code': error_state}

    return params_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def append_results_array(params_out_dict, noct_flag):
    
    for var in ['alpha', 'beta', 'k', 'error_code']:
        params_out_dict[var] = np.empty(len(params_out_dict['rb']))
        if not var == 'error_code':
            params_out_dict[var][:] = np.nan
        else:
            params_out_dict[var][:] = 20
    if not noct_flag:
        params_out_dict.pop('rb_error_code')
        params_out_dict['rb'] = np.empty(len(params_out_dict['rb']))
        params_out_dict['rb'][:] = np.nan
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
def calculate_light_response(data_dict,
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
                     for var in ['TempC', 'VPD', 'PAR', 'NEE_series']}

        # Get the index for the current date
        date_index = np.where(params_out_dict['date'] == date)

        # Check data availability 
        data_pct = int(len(temp_dict['NEE_series']) / 
                       float((1440 / meas_int) / 2) * 100)

        # Do (or do not - there is no try... actually there is!) the fit
        if not data_pct < min_pct:

            # Get values of respiration parameters
            params_in_dict['Eo_default'] = params_out_dict['Eo'][date_index]
            
            # If using nocturnal rb...
            if configs_dict['use_nocturnal_rb']:
                
                params_in_dict['rb_default'] = params_out_dict['rb'][date_index]
                fit_dict = optimise_fixed_rb(temp_dict, params_in_dict)
            
            # If using daytime rb...                                                             
            else:
                fit_dict = optimise_free_rb(temp_dict, params_in_dict)
            
            # Write data to results arrays
            for key in fit_dict.keys():
                params_out_dict[key][date_index] = fit_dict[key]

            # Write alpha default values to default dictionary if valid
            if fit_dict < 2:
                params_in_dict['alpha_default'] = fit_dict['alpha']
            else:
                params_in_dict['alpha_default'] = 0

        else:

            # Error code for not enough data
            params_out_dict['error_code'][date_index] = 10

    # Filter and interpolate
    filter_index = np.where(params_out_dict['error_code'] < 3)
    for key in fit_dict.keys():
        if not key == 'error_code':
            valid_array = params_out_dict[key][filter_index]
            filt.slide_IQR_filter(valid_array, outlier_value = 1.5)
            for i, this_index in enumerate(filter_index[0]):
                params_out_dict[key][this_index] = valid_array[i]
            params_out_dict[key] = gf.generic_2d_linear(params_out_dict[key])

    # Rename the error code variable
    params_out_dict['light_response_error_code'] = params_out_dict.pop('error_code')
    
    return
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def estimate_GPP_Re(data_dict, 
                    all_params_dict, 
                    datetime_input_index_dict):
    
    # Create output dicts for estimates
    results_dict = {}
    results_dict['date_time'] = data_dict['date_time']
    for var in ['Re', 'GPP']:
        results_dict[var] = np.empty(len(data_dict['TempC']))
        results_dict[var][:] = np.nan
    
    # Estimate time series GPP and Re and write to results dictionary
    for i, date in enumerate(all_params_dict['date']):
        
        indices = datetime_input_index_dict[date]
        this_Eo = all_params_dict['Eo'][i]
        this_rb = all_params_dict['rb'][i]
        this_alpha = all_params_dict['alpha'][i]
        this_beta = all_params_dict['beta'][i]
        this_k = all_params_dict['k'][i]
        
        this_dict = {var: data_dict[var][indices[0]: indices[1] + 1]
                     for var in ['NEE_series', 'VPD', 'PAR', 'TempC']}
        GPP, Re = LRF_part(this_dict, 
                           this_Eo,
                           this_rb, 
                           this_alpha,
                           this_beta,
                           this_k) 
        results_dict['GPP'][indices[0]: indices[1] + 1] = GPP
        results_dict['Re'][indices[0]: indices[1] + 1] = Re

    return results_dict
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
    
    x_lab = '$PPFD\/(\mu mol\/photons\/m^{-2} s^{-1})$'
    
    for date in step_data_dict.keys():

        Eo = params_dict['Eo'][params_dict['date'] == date]
        rb = params_dict['rb'][params_dict['date'] == date]
        alpha = params_dict['alpha'][params_dict['date'] == date]
        beta = params_dict['beta'][params_dict['date'] == date]
        k = params_dict['k'][params_dict['date'] == date]

        bool_filter = step_data_dict[date]['day_bool']
        x_var = step_data_dict[date]['PAR'][bool_filter]
        y_var1 = step_data_dict[date]['NEE_series'][bool_filter]
        index = x_var.argsort()
        x_var = x_var[index]
        y_var1 = y_var1[index]
        GPP, Re = LRF_part(step_data_dict[date], 
                           Eo, 
                           rb,
                           alpha,
                           beta,
                           k)
        y_var2 = (GPP + Re)[bool_filter]
        y_var2 = y_var2[index]

        # Plot
        date_str = dt.datetime.strftime(date,'%Y-%m-%d')
        fig = plt.figure(figsize = (12,8))
        fig.patch.set_facecolor('white')
        ax = plt.gca()
        ax.plot(x_var, y_var1, 'o' , markerfacecolor = 'none',
                 markeredgecolor = 'black', label = 'Observed', color = 'black')
        ax.plot(x_var, y_var2, '^', color = 'black', label = 'Estimated')
        ax.set_title('Fit for ' + str(window) + ' day window centred on ' + 
                      date_str + '\n', fontsize = 22)
        ax.set_xlabel(x_lab, fontsize = 18)
        ax.set_ylabel('$NEE\/(\mu mol C\/m^{-2} s^{-1}$)', fontsize = 18)
        ax.axhline(y = 0, color = 'black')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(loc = 'upper right', fontsize = 16, frameon = False, 
                  numpoints = 1)
        plot_out_name = 'day' + '_' + date_str + '.jpg'
        plt.tight_layout()
        fig.savefig(os.path.join(configs_dict['output_path'],
                                 plot_out_name))
        plt.close(fig)

    if is_on: plt.ion()
        
    return
#------------------------------------------------------------------------------
                
#------------------------------------------------------------------------------
#############
# Main code #
#############        

def main(data_dict, configs_dict, params_out_dict):
    """
    Calculates GPP and Re using rectangular hyperbolic radiation response 
    function in conjunction with Lloyd and Taylor Arrhenius-style temperature 
    response function (either Eo or both Eo and rb can be fitted using 
    nocturnal data, and the remaining parameters are fitted to specified 
    combination of step and window;
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
                 - 'VPD': vapour pressure deficit in kilopascals
                 - 'PAR': photosynthetically active radiation in umol photons
                          m^2 s^-1
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
                 - 'use_nocturnal_rb': whether to use nocturnal estimates of
                                       rb or add rb as an additional parameter
                                       to the daytime fit (boolean)
                 - 'output_fit_plots': whether to output plots showing the fit
                                       of the parameters to the data (boolean)
                 - 'output_path': path for output of plots
    Returns: 1) a results dictionary containing 2 key / value pairs:
                    - 'date_time': numpy array of Python datetimes for each
                      datum in the original time series
                    - 'Re': numpy array of half-hourly estimates of Re
                    - 'GPP': numpy array of half-hourly estimates of Re
             2) a results dictionary containing 5 key / value pairs:
                    - 'date_time': numpy array of Python dates for each day in 
                                   the original time series 
                    - 'Eo': numpy array of Eo estimates for each day (note that 
                            each year has a constant value)
                    - 'Eo_error_code': numpy array of Eo diagnostic errors for 
                                       each day (note that each year has a 
                                       constant value)
                    - 'rb_noct': numpy array of rb estimates from nocturnal fit
                                 for each day (note that linear interpolation 
                                 is used to gap fill between steps)
                    - 'rb_error_code': numpy array of rb diagnostic errors for 
                                       each day (including whether the value is 
                                       interpolated or calculated)
                    - 'alpha': numpy array of initial slope of radiation / NEE
                               curve
                    - 'beta': numpy array of maximum photosynthetic rate
                    - 'k': numpy array of exponent coefficient of VPD 
                           function used to modify beta
    """

    # Create boolean indices for masking daytime and nan values
    day_mask = data_dict['Fsd'] > 5
    nan_mask = filt.subset_arraydict_on_nan(data_dict,
                                            var_list = ['NEE_series', 'TempC',
                                                        'VPD', 'PAR'],
                                            subset = False)
    all_mask = [all(rec) for rec in zip(day_mask, nan_mask)]
    data_dict['day_bool'] = np.array(day_mask)
    data_dict['all_bool'] = np.array(all_mask)

    # Partition the data into year and step pieces
    step_data_dict = dtf.get_moving_window(data_dict, 
                                           'date_time', 
                                           configs_dict['window_size_days'],
                                           configs_dict['step_size_days'])

    # Get the indices of the start and end rows of each unique date in the source 
    # data array - no data dict is built from this, since these indices are used to
    # assign Re estimates to the estimated time series output array only
    dates_input_index_dict = dtf.get_day_indices(data_dict['date_time'])

    # Initalise parameter dicts with prior estimates    
    params_in_dict = {'k_prior': 0,
                      'alpha_prior': -0.01,
                      'rb_prior': params_out_dict['rb'].mean(),
                      'beta_prior': (np.percentile(data_dict['NEE_series']
                                                   [data_dict['all_bool']], 5) - 
                                     np.percentile(data_dict['NEE_series']
                                                   [data_dict['all_bool']], 95)),
                      'alpha_default': 0,
                      'beta_default': 0,
                      'k_default': 0 }

    # Append the results dictionary for the parameter values (1 for each day)
    append_results_array(params_out_dict, configs_dict['use_nocturnal_rb'])

    # Write the light response parameters to the existing parameters dict
    calculate_light_response(step_data_dict,
                             configs_dict,
                             params_in_dict,
                             params_out_dict)
    
    # Estimate Re and GPP for all data
    rslt_dict = estimate_GPP_Re(data_dict,
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