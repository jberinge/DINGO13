# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 11:50:13 2016

@author: imchugh
"""

import numpy as np

import data_filtering as data_filter
import respiration_photosynthesis_run as rp_run

def combine_series(be_filled_array, fill_array):

    return np.where(~np.isnan(be_filled_array), be_filled_array, fill_array)

def calculate_annual(data_dict, use_var,
                     mean_or_sum = 'mean', multiplier = 0, meas_int = 30):

    if len(data_dict[use_var]) != len(data_dict[use_var]
                                               [~np.isnan(data_dict[use_var])]):
        print 'Warning: series ' + use_var + ' contains missing data!'
        
    # Create dataset separated into years
    years_data_dict = data_filter.subset_datayear_from_arraydict(data_dict, 
                                                                'date_time')     
    # 
    rslt_dict = {}
    for yr in years_data_dict.keys():
       if mean_or_sum == 'sum':
           scale_to_sum = len(years_data_dict[yr][use_var]) * 60 * meas_int
       else:
           scale_to_sum = 1
       rslt_dict[yr] = np.round(years_data_dict[yr][use_var].mean() * 
                                multiplier * 
                                scale_to_sum, 1)
    
    return rslt_dict

def NEE_fill_and_sum():

    reload(rp_run)
    
    rslt, params = rp_run.main(do_light_response = True)
    
    rslt['NEE_filled'] = combine_series(rslt['NEE_series'], 
                                        rslt['GPP'] + rslt['Re'])
    
    return calculate_annual(rslt, 'NEE_filled', mean_or_sum = 'sum', 
                            multiplier = 12 * 10**-6)    
    
#    # Generate model estimates
#    rslt_dict['NEE_mod'] = rslt_dict['Re'] + rslt_dict['GPP']
#    
#    n_cases = len(rslt_dict['NEE_series'][~np.isnan(rslt_dict['NEE_series'])])
#    print 'The number of available data = ' + str(n_cases)
#    
#    # Gap fill
#    rslt_dict['NEE_filled'] = rslt_dict['NEE_series']
#    rslt_dict['NEE_filled'][np.isnan(rslt_dict['NEE_filled'])] = \
#        rslt_dict['NEE_mod'][np.isnan(rslt_dict['NEE_filled'])]
#    
#    # Create dataframe
#    df = pd.DataFrame(rslt_dict, index = rslt_dict['date_time'])
#    df.drop('date_time', inplace = True, axis = 1)
#    df['NEE_sum'] = df['NEE_filled'] * 0.0018 * 12
#    
#    for yr in ['2012', '2013', '2014']:
#        
#        print yr + ': ' + str(df.loc[yr, 'NEE_sum'].sum())