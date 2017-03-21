# -*- coding: utf-8 -*-
import math
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import gridspec
import pdb

#-----------------------------------------------------------------------------#
def regress_sigma_delta(data_dict, configs_dict):

    """
    Pass the following arguments: 1) dict containing non-gap filled QC'd NEE, 
                                     Fsd, Ta, ws, NEE_model
                                  2) dict containing config options as 
                                     specified below
    
    Returns the linear regression statistics for daytime and nocturnal data
    
    Config dictionary should contain the following:
        'measurement_interval' (of fluxes, in minutes)
        'pos_averaging_bins' (number of bins for calculation of random error 
                              variance as function of +ve flux magnitude)
        'neg_averaging_bins' (number of bins for calculation of random error 
                              variance as function of -ve flux magnitude)                              
        'radiation_difference_threshold' (the maximum allowable Fsd difference 
                                          between NEE pairs; see below for 
                                          threshold values)
        'temperature_difference_threshold' (as above)
        'windspeed_difference_threshold' (as above)
        'mean_series' (series to use to calculate the NEE mean that is paired 
                       with the random error estimate - obs can be used but 
                       not recommended since there are only 2 values and both 
                       are affected by random error)
                                     
    Algorithm from reference:
        Hollinger, D.Y., Richardson, A.D., 2005. Uncertainty in eddy covariance measurements 
        and its application to physiological models. Tree Physiol. 25, 873–885.
    
    Uses daily differencing procedure to estimate random error 
    (note this will overestimate by factor of up to 2)
    
    NEE pairs must pass difference constraints as follows 
    (as suggested in ref above):
        Fsd:35W m^-2
        Ta: 3C
        ws: 1m s^-1
    """      
    
    print '\nCalculating random error regression coefficients'
    print '------------------------------------------------\n'

    # Calculate records per day from measurement interval
    recs_per_day = 1440 / configs_dict['measurement_interval']

    # Do the paired differencing (absolute)
    diff_dict = {}
    diff_dict = {'NEE_mean': (data_dict[configs_dict['mean_series']] + 
                              np.roll(data_dict[configs_dict['mean_series']], 
                                     recs_per_day)) / 2,
                 'NEE_diff_abs': abs(data_dict['NEE_series'] -
                                     np.roll(data_dict['NEE_series'], recs_per_day)),
                 'NEE_diff': (data_dict['NEE_series'] -
                             np.roll(data_dict['NEE_series'], recs_per_day)),                
                 'Ta_diff': abs(data_dict['TempC'] -
                                np.roll(data_dict['TempC'], recs_per_day)),
                 'ws_diff': abs(data_dict['ws'] -
                                np.roll(data_dict['ws'], recs_per_day)),
                 'Fsd_diff': abs(data_dict['Fsd'] -
                                 np.roll(data_dict['Fsd'], recs_per_day))}

    # Remove entire record if nan for any variable, then count available pairs
    temp_array = np.empty([len(diff_dict['NEE_diff_abs']), len(diff_dict)])
    for i, var in enumerate(diff_dict.keys()):
        temp_array[:, i] = diff_dict[var]
    temp_array = temp_array[recs_per_day:, :]
    QCdata_index = np.where(np.all(~np.isnan(temp_array), axis=1))    
    temp_array = temp_array[QCdata_index]
    diff_dict = {var: temp_array[:, i] for i, var in enumerate(diff_dict.keys())}
    total_tuples = str(len(diff_dict['NEE_diff_abs']))
    
    # Remove any values that don't pass the difference constraints
    pass_index = np.where((diff_dict['Ta_diff'] < 
                           configs_dict['temperature_difference_threshold']) & 
                          (diff_dict['ws_diff'] < 
                           configs_dict['windspeed_difference_threshold']) & 
                          (diff_dict['Fsd_diff'] < 
                           configs_dict['radiation_difference_threshold']))
    for key in diff_dict.keys():
        diff_dict[key] = diff_dict[key][pass_index]
    passed_tuples = str(len(diff_dict['NEE_diff_abs']))
               
    # Separate out positive and negative values
    neg_dict = {}
    pos_dict = {}
    for var in ['NEE_mean', 'NEE_diff_abs']:
        neg_dict[var] = diff_dict[var][diff_dict['NEE_mean'] < 0]
        pos_dict[var] = diff_dict[var][diff_dict['NEE_mean'] > 0]
    dict_list = [neg_dict, pos_dict]

    # Report stats
    num_pos = str(len(pos_dict['NEE_diff_abs']))
    num_pos_per_bin = str(int(float(num_pos) / configs_dict['pos_averaging_bins']))
    num_neg = str(len(neg_dict['NEE_diff_abs']))
    num_neg_per_bin = str(int(float(num_neg) / configs_dict['neg_averaging_bins']))
    print (passed_tuples +' of ' + total_tuples + 
           ' available tuples passed difference constraints (Fsd < ' +
           str(configs_dict['radiation_difference_threshold']) + 'Wm^-2, Ta < ' + 
           str(configs_dict['temperature_difference_threshold']) + 
           'C, ws < ' + str(configs_dict['windspeed_difference_threshold']) + 
           'ms^-1):')
    print ('    - ' + num_neg + ' records for NEE < 0 (' + num_neg_per_bin + 
           ' records per bin)')
    print ('    - ' + num_pos + ' records for NEE > 0 (' + num_pos_per_bin + 
           ' records per bin)') 

    # Create arrays to takes results of quantile categorisation
    neg_array = np.empty([configs_dict['neg_averaging_bins'], 2])
    pos_array = np.empty([configs_dict['pos_averaging_bins'], 2])
    array_list = [neg_array, pos_array]
    cond_list = ['neg', 'pos']
    stats_list = ['slope','intcpt','r_val','p_val', 'SE']
    stats_dict = {}

    # Now categorise on basis of percentiles and calculate means and sigmas
    for i in range(2):
        this_dict = dict_list[i]
        this_array = array_list[i]
        num_bins = np.shape(this_array)[0]
        for j, num in enumerate(np.linspace(100.0 / num_bins, 100, num_bins)):
            pctl_hi = np.percentile(this_dict['NEE_mean'], num)
            pctl_lo = np.percentile(this_dict['NEE_mean'],
                                    num - 100.0 / num_bins)
            quantile_index = np.where((this_dict['NEE_mean'] > pctl_lo) &
                                      (this_dict['NEE_mean'] <= pctl_hi))                                
            this_array[j, 0] = this_dict['NEE_mean'][quantile_index].mean()
            this_array[j, 1] = ((abs(this_dict['NEE_diff_abs'][quantile_index] - 
                                     this_dict['NEE_diff_abs'][quantile_index].mean()))
                                 .mean() * np.sqrt(2))
                                  
        # Calculate linear fit for positive and negative values...
        linreg_stats_list = stats.linregress(this_array[:, 0],
                                             this_array[:, 1]) 
        stats_dict[cond_list[i]] = {stats_list[stat]: linreg_stats_list[stat] 
                                    for stat in range(5)}

    # Combine pos and neg arrays
    combined_array = np.concatenate([neg_array, pos_array])
    rslt_dict = {'NEE_mean': combined_array[:, 0],
                 'sig_del': combined_array[:, 1]}

    ### Plotting ###

    # Instantiate plot
    fig = plt.figure(figsize = (12, 6))
    fig.patch.set_facecolor('white')
    gs = gridspec.GridSpec(1, 2, width_ratios = [1, 2])
    ax1 = plt.subplot(gs[0])    
    ax2 = plt.subplot(gs[1])

    ### Histogram ###    
    
    # Calculate scaling parameter for Laplace (sigma / sqrt(2)) and Gaussian 
    # (sigma) distributions over entire dataset 
    beta = (abs(diff_dict['NEE_diff'] - diff_dict['NEE_diff'].mean())).mean()
    sig = diff_dict['NEE_diff'].std()

    # Get edge quantiles and range, then calculate Laplace pdf over range
    x_low = myround(np.percentile(diff_dict['NEE_diff'], 0.5))
    x_high = myround(np.percentile(diff_dict['NEE_diff'], 99.5))
    x_range = (x_high - x_low)
    x = np.arange(x_low, x_high, 1 / (x_range * 10.))
    pdf_laplace = np.exp(-abs(x / beta)) / (2. * beta)
	    	
    # Plot normalised histogram with Laplacian and Gaussian pdfs
    ax1.hist(np.array(diff_dict['NEE_diff']), bins = 200, 
             range = [x_low, x_high], normed = True, color = '0.7', 
             edgecolor = 'none')
    ax1.plot(x,pdf_laplace,color='black', label='Laplacian')
    ax1.plot(x,mlab.normpdf(x,0,sig),color='black', linestyle = '--', 
            label='Gaussian')
    ax1.set_xlabel(r'$\delta\/(\mu mol\/m^{-2} s^{-1}$)',fontsize=22)
    ax1.set_ylabel('$Frequency$', fontsize=22)
    ax1.axvline(x=0, color = 'black', linestyle = ':')
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax1.tick_params(axis = 'y', labelsize = 14)

    ### Line plot ###

    # Create series for regression lines
    influx_x = np.append(neg_array[:, 0], np.array(0))
    influx_y = np.polyval([stats_dict['neg']['slope'], 
                           stats_dict['neg']['intcpt']],
                           influx_x)
    efflux_x = np.append(np.array(0), pos_array[:, 0])
    efflux_y = np.polyval([stats_dict['pos']['slope'], 
                           stats_dict['pos']['intcpt']],
                           efflux_x)
    
    # Do formatting
    ax2.set_xlim(round(rslt_dict['NEE_mean'][0]), 
                 math.ceil(rslt_dict['NEE_mean'][-1]))
    ax2.set_xlabel(r'$C\/flux\/(\mu molC\/m^{-2} s^{-1}$)',fontsize=22)
    ax2.set_ylim([0, math.ceil(rslt_dict['sig_del'].max())])    
    ax2.set_ylabel('$\sigma(\delta)\/(\mu molC\/m^{-2} s^{-1})$',fontsize=22)
    ax2.tick_params('x', labelsize = 14)    
    ax2.yaxis.set_ticklabels([])
    ax2.tick_params(axis='y', which='both', left='off', right='off')
    ax3 = ax2.twinx()
    ax3.spines['right'].set_position('zero')
    ax3.set_ylim(ax2.get_ylim())
    ax3.tick_params(axis = 'y', labelsize = 14)
    plt.setp(ax3.get_yticklabels()[0], visible = False)

    # Do plotting
    ax2.plot(rslt_dict['NEE_mean'], rslt_dict['sig_del'], 'o', 
             markerfacecolor='0.8', markeredgecolor='black', markersize=6)
    ax2.plot(influx_x, influx_y, linestyle=':', color='black')
    ax2.plot(efflux_x, efflux_y, linestyle=':', color='black')
    
    fig.tight_layout()
    plt.close()

    return fig, stats_dict, rslt_dict

#-----------------------------------------------------------------------------#
def estimate_sigma_delta(NEE_array, stats_dict):
    
    """
    Pass the following positional arguments:
        1) array containing NEE estimates (a model is recommended)
        2) dict containing the regression statistics for positive and negative 
           values of NEE
           
    Returns a numpy array of sigma_delta estimates
    """    

    # Calculate the estimated sigma_delta for each datum
    return np.where(NEE_array > 0, 
                    NEE_array * stats_dict['pos']['slope'] 
                              + stats_dict['pos']['intcpt'],
                    NEE_array * stats_dict['neg']['slope'] 
                              + stats_dict['neg']['intcpt'])

#-----------------------------------------------------------------------------#
def estimate_random_error(sig_del_array):
    """
    Pass the following arguments: 1) array containing sigma_delta estimate for NEE
    
    Returns a numpy array of random error estimate drawn from the Laplace 
    distribution with sigma_delta as the scaling parameter
    (location parameter is 0)
    """

    return np.random.laplace(0, sig_del_array / np.sqrt(2))

def propagate_random_error(sig_del_array, configs_dict):
    """
    Pass the following arguments: 1) array containing sigma_delta estimate for NEE
                                  2) dict containing config options as specified
                                     below
                                     
    Config dictionary should contain the following:
    'measurement_interval' (of fluxes)
    'num_trials' (number of required Monte Carlo simulations)    
    """

    # Calculate critical t-statistic for p = 0.095
    crit_t = stats.t.isf(0.025, configs_dict['num_trials'])
    
    # Iterate over number of trials
    result_array = np.empty(configs_dict['num_trials'])
    for this_trial in xrange(configs_dict['num_trials']):
        error_array = estimate_random_error(sig_del_array)
        result_array[this_trial] = (error_array.sum() * 
                                    configs_dict['measurement_interval'] * 60 *
                                    12 * 10 ** -6)
                                
    return np.round(result_array.std() * crit_t, 2)

def myround(x,base=10):
    return int(base*round(x/base))