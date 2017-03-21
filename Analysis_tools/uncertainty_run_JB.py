# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:41 2015

@author: imchugh
"""

# Standard modules
import numpy as np
import os
import copy as cp
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.gridspec as gridspec
import warnings
import pdb

# My modules
import DataIO as io
import data_filtering as filt
import datetime_functions as dtf
import random_error as rand_err
import model_error as mod_err
import data_formatting as dt_fm
import respiration as re
import photosynthesis as ps
import gap_filling as gf

#JB extra
import pandas as pd

#------------------------------------------------------------------------------
# Fetch data from configurations
def get_data(configs_dict):

    # Get data (screen Fc data to obs only - keep gap-filled drivers etc)
    data_input_target = os.path.join(configs_dict['files']['input_path'],
                                     configs_dict['files']['input_file'])

    ext = os.path.splitext(data_input_target)[1]
    if ext == '.nc':
        Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                  var_list = [configs_dict['variables']
                                                                          ['carbon_flux']],
                                                  QC_accept_codes = [0])
        Fc_dict.pop('date_time')
        ancillary_vars = [configs_dict['variables'][var] for var in 
                          configs_dict['variables'] if not var == 'carbon_flux']
        ancillary_dict, attr = io.OzFluxQCnc_to_data_structure(
                                   data_input_target,
                                   var_list = ancillary_vars,
                                   return_global_attr = True)
        data_dict = dict(Fc_dict, **ancillary_dict)

    elif ext == '.df':
        data_dict, attr = io.DINGO_df_to_data_structure(data_input_target,
                              var_list = configs_dict['variables'].values(),
                              return_global_attr = True)

    # Rename to generic names used by scripts
    old_names_dict = configs_dict['variables']
    std_names_dict = dt_fm.standard_names_dictionary()
    map_names_dict = {old_names_dict[key]: std_names_dict[key] 
                      for key in old_names_dict}
    data_dict = dt_fm.rename_data_dict_vars(data_dict, map_names_dict)

        
    return data_dict, attr
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Rebuild the master configuration file for passing to respiration and light 
# response (if requested) functions
def build_config_file(configs_master_dict):

    # Build a specific configuration file
    configs_dict = {'files': configs_master_dict['global_configs']['files'],
                    'global_options': (configs_master_dict['global_configs']
                                                          ['options']),
                    'uncertainty_options': configs_master_dict['NEE_uncertainty_configs']}                                                          
    configs_dict['variables'] = {}
    dict_list = [configs_master_dict['random_error_configs']['variables'],
                 configs_master_dict['model_error_configs']['variables'],
                 configs_master_dict['respiration_configs']['variables'],
                 configs_master_dict['photosynthesis_configs']['variables']]
    for d in dict_list:
        configs_dict['variables'].update(d)
    return configs_dict                                                         
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check whether ustar and Fsd (both used for filtering of the dataset) contain
# missing data when NEE data is available - if so, exclude these cases from 
# analysis
def check_data_consistency(data_dict):
    
    warnings.simplefilter('always')
    for var in ['Fsd', 'ustar']:
        
        flag_index = np.where(~np.isnan(data_dict['NEE_series']) & 
                              np.isnan(data_dict[var]))
        count_str = str(len(flag_index[0]))
        if not len(flag_index[0]) == 0:
            warnings.warn('There are %s' %count_str + ' instances where NEE ' 
                          'element contains valid data and %s' %var + ' element ' 
                          'is not a number - NEE values for these ' 
                          'instances will be excluded from analysis!')
        data_dict['NEE_series'][flag_index] = np.nan

#------------------------------------------------------------------------------
# Check whether all model drivers are complete - if not, warn user (may expand 
# this to force exception, since results will be nan)
def check_driver_consistency(data_dict):
    
    warnings.simplefilter('always')
    arr = np.array([])
    for var in ['Fsd', 'TempC', 'VPD']:
        
        flag_index = np.where(np.isnan(data_dict['NEE_series']) & 
                              np.isnan(data_dict[var]))
        arr = np.concatenate([arr, flag_index[0]])
        count_str = str(len(flag_index[0]))                              
        if not len(flag_index[0]) == 0:
            warnings.warn('There are %s' %count_str + ' instances where neither ' \
                          'the NEE nor model driver %s' %var + ' element ' \
                          'contains valid data - model estimates cannot be ' \
                          'calculated for these instances!')
        data_dict['NEE_series'][flag_index] = np.nan        
    arr = np.unique(arr)
    if not len(arr) == 0:
        print 'Total number of instances with missing driver data is ' + str(len(arr))
#------------------------------------------------------------------------------
        
#------------------------------------------------------------------------------
# Set all sigma delta values to nan where there are no observational data
def filter_sigma_delta(data_dict):
    data_dict['sigma_delta'][np.isnan(data_dict['NEE_series'])] = np.nan
#------------------------------------------------------------------------------    

#------------------------------------------------------------------------------
# Split into day and night
def separate_night_day(data_dict, noct_threshold):
    subset_dict = {}
    subset_dict['day'] = filt.subset_arraydict_on_threshold(
                             data_dict, 'Fsd', noct_threshold, '>', drop = True)
    subset_dict['night'] = filt.subset_arraydict_on_threshold(
                               data_dict, 'Fsd', noct_threshold, '<', drop = True)
    return subset_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Make a results dictionary
def init_interm_rslt_dict(num_trials, do_ustar, do_random, do_model):

    var_list = ['obs_avail_day', 'obs_avail_night']
    if do_ustar:
        var_list = var_list + ['u_star', 'ustar_error']
    if do_random:
        var_list = var_list + ['random_error_day', 'random_error_night']
    if do_model:
        var_list = var_list + ['model_error_day', 'model_error_night']

    nan_array = np.zeros(num_trials)
    nan_array[:] = np.nan
    return {var: cp.copy(nan_array) for var in var_list} 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def run_model(data_dict, NEE_model, re_configs_dict, ps_configs_dict):
    
    if NEE_model == 'LT':    
        
        re_rslt_dict, re_params_dict = re.main(data_dict, re_configs_dict)[0: 2]
        ps_rslt_dict = ps.main(data_dict, ps_configs_dict, re_params_dict)[0]
        data_dict['NEE_model'] = ps_rslt_dict['GPP'] + ps_rslt_dict['Re']
        data_dict['NEE_filled'] = np.where(np.isnan(data_dict['NEE_series']),
                                           data_dict['NEE_model'],
                                           data_dict['NEE_series'])
                                           
    elif NEE_model == 'ANN':
        
        len_int = len(data_dict['NEE_series'])
        input_array = np.empty([len_int, 4])
        for i, var in enumerate(['TempC', 'Sws', 'Fsd', 'VPD']):
            input_array[:, i] = data_dict[var]
        target_array = np.empty([len_int, 1])
        target_array[:, 0] = data_dict['NEE_series']
        
        data_dict['NEE_model'] = gf.train_ANN(input_array, target_array, 
                                              100, 
                                              [4, 24, 16, 1])[:, 0]
        data_dict['NEE_filled'] = np.where(np.isnan(data_dict['NEE_series']),
                                           data_dict['NEE_model'],
                                           data_dict['NEE_series'])                                                     

    else:
        
        raise Exception('\'' + NEE_model + '\' is not a valid model type! ' \
                        'Valid choices are \'ANN\' or \'LT\'')
                                           
    return    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def plot_data(data_d,configs_dict):
    """
    Pass a dictionary containing the following key / value pairs:
    keys: year(int or str) / values: dictionary containing following 
    key / value pairs: 
    keys: 'ustar_error', 'model_error_day', 'model_error_night', 
    'random_error_day' and 'random_error_night' as keys / values: equal length 
    numpy arrays (may contain np.nan - will be filtered)
    """
    
    key_id_list = [isinstance(key, int) for key in data_d.keys()]    
    if not all(key_id_list):
        raise Exception('Expected integer years as outer dictionary key!')
    
    for this_year in data_d.keys():
    
        this_d = data_d[this_year]    
    
        ustar_err = this_d['ustar_error'][~np.isnan(this_d['ustar_error'])]    
        
        rand_err = this_d['random_error_day'] + this_d['random_error_night']
        rand_err = rand_err[~np.isnan(this_d['ustar_error'])]
        
        mod_err = this_d['model_error_day'] + this_d['model_error_night']
        mod_err = mod_err[~np.isnan(this_d['ustar_error'])]
        
        total_err = rand_err + mod_err + ustar_err
    
        labels = ['total', 'model', 'random', 'ustar']
        colors = ['grey', 'magenta', 'cyan', 'blue']
        pos = [0.9, 0.6, 0.3, 0]
    
        # Do the stats
        mu_dict = {}
        sig_dict = {}
        this_list = [total_err, mod_err, rand_err, ustar_err]
        for i, this_one in enumerate(this_list):
            mu_dict[labels[i]] = this_one.mean()
            sig_dict[labels[i]] = this_one.std()
                
        # Create the plot
        fig = plt.figure(figsize = (12, 10))
        fig.patch.set_facecolor('white')
        gs = gridspec.GridSpec(2, 1, height_ratios=[4,1.5])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
    
        # Set up the first subplot
        ax1.set_xlabel('$Uncertainty (g C m^{-2}a^{-1})$',
                      fontsize=18)
        ax1.set_ylabel('$Frequency$', fontsize=18)
        ax1.tick_params(axis = 'y', labelsize = 14)
        ax1.tick_params(axis = 'x', labelsize = 14)
        ax1.axvline(0, color = 'black', linewidth = 2, linestyle = '--')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        
        # Plot the histogram
        for i, var in enumerate(this_list):
            if i == 0:
                ax1.hist(total_err, 50, facecolor = colors[i], edgecolor = 'none',
                         orientation = 'vertical', label = labels[i], normed = True)
            else:
                ax1.hist(var, 50, facecolor = 'none', edgecolor = colors[i], 
                         histtype = 'step', orientation = 'vertical',
                         label = labels[i], normed = True)
        ax1.legend(loc='upper right', frameon = False)

        # Put year in plot
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        ax1.text(xmin + (xmax - xmin) / 10, ymax - (ymax - ymin) / 20, 
                 str(this_year), fontsize = 20)
    
        # Plot the normal distribution
        x = np.linspace(xmin, xmax, 100)
        p = stats.norm.pdf(x, mu_dict['total'], sig_dict['total'])
        ax1.plot(x, p, color = 'black')
    
        # Set up the second plot
        ax2.axes.get_yaxis().set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.xaxis.set_ticks_position('bottom')
        ax2.set_xticklabels([])                    
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_ylim([0, 1])
        the_buffer = 0.12
        
        # Plot the confidence intervals
        for i, this_one in enumerate(labels[:-1]):
            ax2.plot((mu_dict[this_one] - sig_dict[this_one] * 2, 
                      mu_dict[this_one] + sig_dict[this_one] * 2), 
                     (pos[i], pos[i]), color = colors[i], linewidth = 2)
            ax2.plot(mu_dict[this_one] - sig_dict[this_one] * 2, pos[i], 
                     marker = '|', color = colors[i], markersize = 10,
                     mew = 2)
            ax2.plot(mu_dict[this_one] + sig_dict[this_one] * 2, pos[i], 
                     marker = '|', color = colors[i], markersize = 10,
                     mew = 2)
            ax2.plot(mu_dict[this_one], pos[i], 
                     marker = 'o', color = colors[i], markersize = 10,
                     mec = 'none')
            ax2.text(mu_dict[this_one] - sig_dict[this_one] * 2, pos[i] - the_buffer, 
                     str(round(mu_dict[this_one] - sig_dict[this_one] * 2, 1)),
                     verticalalignment = 'center',
                     horizontalalignment = 'right',
                     fontsize = 14)
            ax2.text(mu_dict[this_one] + sig_dict[this_one] * 2, pos[i] - the_buffer, 
                     str(round(mu_dict[this_one] + sig_dict[this_one] * 2, 1)),
                     verticalalignment = 'center',
                     horizontalalignment = 'left',
                     fontsize = 14)
        mypathforResults=configs_dict['files']['output_path']
        fig.savefig(mypathforResults+'/Uncertainty plot for '+str(this_year))
        print "finish plot for " + str(this_year)
        
    return    

def main(output_trial_results = True, output_plot = True, CPDdata = None, DataDF = None, output_path = None):

    # Update
    reload(rand_err)
    reload(mod_err)
    reload(io)
    reload(filt)
    reload(re)
    reload(gf)

    #-----------------------------------
    # General preparation and formatting
    #-----------------------------------

    # Get master config file
    configs_master_dict = io.config_to_dict(io.file_select_dialog())

    # Build custom configuration file for this script
    configs_dict = build_config_file(configs_master_dict)

    # Get the data
    data_dict, attr = get_data(configs_dict)

    # Build required configuration files for imported scripts (random error,
    # model error, respiration, light response)
    rand_err_configs_dict = configs_master_dict['random_error_configs']['options']
    mod_err_configs_dict = configs_master_dict['model_error_configs']['options']
    re_configs_dict = configs_master_dict['respiration_configs']['options']
    ps_configs_dict = configs_master_dict['photosynthesis_configs']['options']

    # Save the time step information into the individual configuration files
    for d in [configs_dict, rand_err_configs_dict, mod_err_configs_dict, 
              re_configs_dict, ps_configs_dict]: 
        d['measurement_interval'] = int(attr['time_step'])
    
    # For respiration and light response, turn off the output option for 
    # window fits of the functions, even if requested in the configuration file
    # - they are WAY too computationally expensive!!!
    if re_configs_dict['output_fit_plots']:
        re_configs_dict['output_fit_plots'] = False 
    if ps_configs_dict['output_fit_plots']:
        ps_configs_dict['output_fit_plots'] = False 
    
    # Write universal config items to local variables
    noct_threshold = configs_dict['global_options']['noct_threshold']
    ustar_threshold = configs_dict['global_options']['ustar_threshold']
    num_trials = configs_dict['uncertainty_options']['num_trials']
    do_ustar_uncertainty = configs_dict['uncertainty_options']['do_ustar_uncertainty']
    do_random_uncertainty = configs_dict['uncertainty_options']['do_random_uncertainty']
    do_model_uncertainty = configs_dict['uncertainty_options']['do_model_uncertainty']
    NEE_model = configs_dict['uncertainty_options']['NEE_model']

    # Print stuff
    print '\nRunning uncertainty analysis for:'
    error_list = ['ustar', 'random', 'model']
    mode_count = 0
    for i, var in enumerate([do_ustar_uncertainty, do_random_uncertainty, 
                             do_model_uncertainty]):
        if var:
            mode_count = mode_count + 1
            print '- ' + error_list[i] + ' error'
    if mode_count == 0:
        raise Exception('Processing flags for all uncertainty sources ' \
                        'currently set to False: set at least one ' \
                        'uncertainty source to True in configuration file ' \
                        'before proceeding!')
    print '---------------------------------'


    #-----------------
    # Data preparation
    #-----------------

    # Sum Fc and Sc if storage is to be included, otherwise if requested, 
    # remove all Fc where Sc is missing
    if configs_dict['global_options']['use_storage']:
        data_dict['NEE_series'] = (data_dict['NEE_series'] + 
                                   data_dict['Sc'])
    elif configs_dict['global_options']['unify_flux_storage_cases']:
        data_dict['NEE_series'][np.isnan(data_dict['Sc'])] = np.nan

    # Convert insolation to PPFD for light response calculations
    data_dict['PAR'] = data_dict['Fsd'] * 0.46 * 4.6       

    # Check no NEE values with missing ustar values
    check_data_consistency(data_dict)

    # Check no drivers missing where NEE is missing
    check_driver_consistency(data_dict)

    #----------------------------------------
    # Random error calculation and statistics
    #----------------------------------------

    if do_random_uncertainty:

        # Generate initial model series for Re and GPP then combine
        # (note: low u* data is left in intentionally)
        run_model(data_dict, NEE_model, re_configs_dict, ps_configs_dict)
    
        # Calculate the linear regression parameters of sigma_delta as a function 
        # of flux magnitude
        fig, stats_dict, rslt_dict = rand_err.regress_sigma_delta(
                                         data_dict, rand_err_configs_dict)
    
        # Calculate estimated sigma_delta for each data point, and remove records 
        # where no observational estimate is available (only has an effect if the 
        # propagation series is a model - which is recommended!!!);                             
        sig_del_array = (rand_err.estimate_sigma_delta
                            (data_dict[rand_err_configs_dict['propagation_series']], 
                             stats_dict))
        data_dict['sigma_delta'] = sig_del_array


    #---------------------
    # Uncertainty analysis
    #---------------------

    # Create dataset separated into years
    years_data_dict = dtf.subset_datayear_from_arraydict(data_dict, 
                                                          'date_time')   

    # Create a results dictionary
    final_rslt_dict = {}
    final_summary_dict = {}
    
    # Get the t-statistic for the 95% CI
    t_stat = stats.t.ppf(0.975, num_trials)
    
    years = years_data_dict.keys()
    years.sort()
    
    # Do the uncertainty analysis for each year        
    for this_year in years:

        print 'Running analysis for year ' + str(this_year) + ':'

        # Make an intermediate results dictionary
        interm_rslt_dict = init_interm_rslt_dict(num_trials,
                                                 do_ustar_uncertainty,
                                                 do_random_uncertainty,
                                                 do_model_uncertainty)

        # Make an intermediate summary dict
        interm_summary_dict = {}
        
        # Write ustar thresholds and uncertainties for years to local variables
        ustar_list = []
        
        #Here check to see if input list has been passed.  If so use it otherwise read from config file
        #if ustar_uncertainty_list:
            #ustar_list = ustar_uncertainty_list
        #else:
        for var in ['ustar_threshold', 'ustar_uncertainty']:
            if isinstance(configs_dict['global_options'][var], dict):
                ustar_list.append(configs_dict['global_options']
                                              [var][str(this_year)])
            else:
                ustar_list.append(configs_dict['global_options'][var])
        ustar_threshold = ustar_list[0]    
        ustar_uncertainty = ustar_list[1]
            
        ustar_threshold = CPDdata['ustar_mean'][CPDdata[0] == this_year].value() 
        ustar_uncertainty = ustar_list[1]            
            
        # Generate a standard data dictionary; this will be overwritten if
        # ustar uncertainty is set to True, but the reference value for NEE
        # will be retained
        this_dict = cp.deepcopy(years_data_dict[this_year])
#        bool_array = np.array([all(rec) for rec in zip(this_dict['Fsd'] < 5, 
#                                                       ~np.isnan(this_dict['NEE_series']))])
#        max_ustar = np.percentile(this_dict['ustar'][bool_array],
#                                  100 - re_configs_dict['minimum_pct_annual'])
#        bool_array = this_dict['Fsd'] < 5
#        miss_array = filt.subset_arraydict_on_nan(this_dict, 
#                                                  var_list = ['ustar',
#                                                              'PAR',
#                                                              'VPD',
#                                                              'TempC',
#                                                              'NEE_series'],
#                                                  condition = 'any',
#                                                  subset = False)
#        twin_array = np.array([all(rec) for rec in zip(bool_array, miss_array)])
#        temp_ustar_array = this_dict['ustar'][twin_array]
#        temp_NEE_array = this_dict['NEE_series'][twin_array]
#        filt_ustar_array = temp_ustar_array[~np.isnan(temp_NEE_array)]
#        max_ustar = np.percentile(filt_ustar_array, 
#                                  100 - re_configs_dict['minimum_pct_annual'])
#        print 'Maximum ustar for year ' + str(this_year) + ' is ' + str(round(max_ustar, 3))                          

        # Screen current dict copy using the best estimate of u_star threshold,
        # run model and calculate annual sum
        filt.screen_low_ustar(this_dict, ustar_threshold, noct_threshold)
        try:
            run_model(this_dict, NEE_model, re_configs_dict, ps_configs_dict)
        except Exception, e:
            print ('    - Excluding the year ' + str(this_year) + 
                   ' - insufficient data!')
            continue # Do the next year
        NEE_sum = (this_dict['NEE_filled'] * 
                   configs_dict['measurement_interval'] *
                   60 * 12 * 10**-6).sum()

        # If including ustar uncertainty:
        #   1) generate an array of ustar values based 
        #      on mu and sigma from change point detection analysis
        #   2) add the resulting array to the intermediate results dict
        #   3) add an empty array to keep NEE error due to ustar 
        if do_ustar_uncertainty:
            ustar_array = np.random.normal(loc = ustar_threshold,
                                           scale = ustar_uncertainty / 2,
                                           size = num_trials)
            interm_rslt_dict['u_star'][:] = ustar_array
           
        # Create a resample disable boolean, which will be set to True if ustar
        # threshold uncertainty is disabled   
        resample_disable_bool = False
            
        # Do trials
        for this_trial in xrange(num_trials):

            # Print progress
            if this_trial == 0:
                print '    - Trial: ' + str(this_trial + 1),
            elif this_trial == num_trials - 1:
                print str(this_trial + 1) + ' ... Done!'
            else:
                print this_trial + 1,

            # If including ustar uncertainty:
            #   1) make a deep copy of the original data so it doesn't get overwritten
            #   2) set ustar threshold 
            #   3) filter for low ustar
            #   4) gap fill the filtered dataset
            #   5) sum, calculate difference relative to best u* and output to dict
            if do_ustar_uncertainty:
                this_dict = cp.deepcopy(years_data_dict[this_year])
                ustar_threshold = ustar_array[this_trial]
                if ustar_threshold > 0: # and ustar_threshold < max_ustar:
                    filt.screen_low_ustar(this_dict, ustar_threshold, noct_threshold)
                    try:
                        run_model(this_dict, NEE_model, re_configs_dict, ps_configs_dict)
                        this_sum = (this_dict['NEE_filled'] * 
                                    configs_dict['measurement_interval'] * 60 *
                                    12 * 10**-6).sum()
                        interm_rslt_dict['ustar_error'][this_trial] = NEE_sum - this_sum                                
                    except Exception, e:
                        print 'Skipped! ustar was: ' + str(ustar_threshold)
                        continue

            # Switch off resampling after first pass if not doing ustar uncertainty                
            if not do_ustar_uncertainty:
                if this_trial == 1:
                    resample_disable_bool = True

            # Before random and model error estimation:
            #   1) screen all sigma_delta values where observations are missing
            #   2) split into day and night
            # * Note this only needs to be done once if ustar uncertainty is disabled
            if not resample_disable_bool:
                if do_random_uncertainty:
                    filter_sigma_delta(this_dict)
                this_dict = separate_night_day(this_dict, noct_threshold)

            # For each of day and night
            for cond in this_dict.keys():

                # If including ustar uncertainty, write the available n to the 
                # intermediate results dictionary for day and night; otherwise,
                # just write it once (since available n won't vary if ustar
                # doesn't)
                if do_ustar_uncertainty:
                    interm_rslt_dict['obs_avail_' + cond][this_trial] = len(
                        this_dict[cond]['NEE_series']
                            [~np.isnan(this_dict[cond]['NEE_series'])])
                else:
                    if not resample_disable_bool:
                        interm_rslt_dict['obs_avail_' + cond][:] = len(
                            this_dict[cond]['NEE_series']
                                [~np.isnan(this_dict[cond]['NEE_series'])])
                    
                # Do the random error and write to correct position in 
                # intermediate results dict
                if do_random_uncertainty:
                    sig_del_array = (this_dict[cond]['sigma_delta']
                                     [~np.isnan(this_dict[cond]['sigma_delta'])])
                    error_array = rand_err.estimate_random_error(sig_del_array)                
                    interm_rslt_dict['random_error_' + cond][this_trial] = (
                        error_array.sum() * configs_dict['measurement_interval'] 
                                          * 60 * 12 * 10 ** -6)

                # Do the model error and write to correct position in 
                # intermediate results dict
                if do_model_uncertainty:
                    sub_dict = cp.deepcopy(this_dict[cond])
                    interm_rslt_dict['model_error_' + cond][this_trial] = (
                        mod_err.estimate_model_error(sub_dict, 
                                                     mod_err_configs_dict))

        # Write the results for the year to the results dictionary
        final_rslt_dict[this_year] = interm_rslt_dict
        
        # Write the results for the year to the intermediate summary dictionary
        vars_list = []
        if do_random_uncertainty:
            vars_list.extend(['random_error_day', 'random_error_night'])
            interm_summary_dict['random_error_day'] = (
                interm_rslt_dict['random_error_day']
                [~np.isnan(interm_rslt_dict['random_error_day'])].std() * t_stat)
            interm_summary_dict['random_error_night'] = (
                interm_rslt_dict['random_error_night']
                [~np.isnan(interm_rslt_dict['random_error_night'])].std() * t_stat)
            interm_summary_dict['random_error_tot'] = (
                interm_rslt_dict['random_error_day'] +
                interm_rslt_dict['random_error_night']).std() * t_stat
        if do_model_uncertainty:
            vars_list.extend(['model_error_day', 'model_error_night'])
            interm_summary_dict['model_error_day'] = (
                interm_rslt_dict['model_error_day']
                [~np.isnan(interm_rslt_dict['model_error_day'])].std() * t_stat)
            interm_summary_dict['model_error_night'] = (
                interm_rslt_dict['model_error_night']
                [~np.isnan(interm_rslt_dict['model_error_night'])].std() * t_stat)
            interm_summary_dict['model_error_tot'] = (
                interm_rslt_dict['model_error_day'] +
                interm_rslt_dict['model_error_night']).std() * t_stat
        if do_ustar_uncertainty:
            vars_list.append('ustar_error')            
            interm_summary_dict['ustar_error'] = (
                interm_rslt_dict['ustar_error']
                [~np.isnan(interm_rslt_dict['ustar_error'])].std() * t_stat)
        l = [interm_rslt_dict[this_var] for this_var in vars_list]
        interm_summary_dict['total_error'] = sum(l).std() * t_stat

        # Write summary results for the year to final summary dictionary
        final_summary_dict[this_year] = interm_summary_dict
        
        
        # Ouput plots
        if output_plot:
            plot_data({this_year: interm_rslt_dict},configs_dict)
    
    import csv
    mypathforResults=configs_dict['files']['output_path']    

    print "Writing out results files"
    for this_year in years: 
        with open(mypathforResults+'Uncertainty_results_' + str(this_year) +'.csv', 'wb') as f:  # Just use 'w' mode in 3.x
            Dicttoprocess = final_rslt_dict[this_year]
            fieldnames= Dicttoprocess.keys()
            w = csv.writer(f)
            w.writerow(fieldnames)
            obs_avail_night=Dicttoprocess['obs_avail_night'].tolist()
            obs_avail_day=Dicttoprocess['obs_avail_day'].tolist()
            u_star= Dicttoprocess['u_star'].tolist()
            model_error_day= Dicttoprocess['model_error_day'].tolist()
            random_error_night= Dicttoprocess['random_error_night'].tolist()
            ustar_error= Dicttoprocess['ustar_error'].tolist()
            random_error_day= Dicttoprocess['random_error_day'].tolist()
            model_error_night = Dicttoprocess['model_error_night'].tolist()

            w.writerows(zip(obs_avail_night, obs_avail_day, u_star, model_error_day, random_error_night, ustar_error, random_error_day, model_error_night))
               
        with open(mypathforResults+'Uncertainty_summary_' + str(this_year) +'.csv', 'wb') as f:  # Just use 'w' mode in 3.x
            Dicttoprocess = final_summary_dict[this_year]
            w = csv.DictWriter(f, Dicttoprocess.keys())
            w.writeheader()
            w.writerow(Dicttoprocess)  
    print "Done writing files.  Returning results" 
    
    if output_trial_results:    
        return final_rslt_dict, final_summary_dict
    else:
        return final_summary_dict
    print "Done" 

CPDfile = "E:/My Dropbox/Dropbox/Data_flux_data/Site data processing/HowardSprings/Advanced_v12a/CPD/Results/annual_statistics.csv"
DataDF = pd.read_pickle("E:/My Dropbox/Dropbox/Data_flux_data/Site data processing/HowardSprings/Advanced_v12a/Advanced_processed_data_HowardSprings_v12a.df")
output_path = "E:/My Dropbox/Dropbox/Data_flux_data/Site data processing/HowardSprings/Advanced_v12a/Diagnostics/Uncertainty/"

CPDdata = pd.read_csv(CPDfile)

main(True, True, CPDdata, DataDF, output_path)