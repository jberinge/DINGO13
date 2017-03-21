# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 02:57:37 2015

Routines for cleaning up profile data

@author: imchugh
"""

import pandas as pd
import os
import pdb
import numpy as np
from scipy import stats
import datetime as dt

def truncate_data(data_df):
    """
    This script takes a non-truncated dataframe of profile CO2 concentrations 
    as argument, then truncates according to the specified window and lag to 
    the specified output interval
    """
    output_interval = 30
    bin_size = 8
    lag = 0

    if not bin_size < output_interval: 
        print 'Bin size must be less than output interval'
        return

    df = data_df.copy()
    df_columns = df.columns

    # Find the time window to be used given the lag and bin size settings above
    df['mins'] = df.index.minute % output_interval
    arr = np.unique(df['mins'])
    select_arr = np.arange(bin_size) - int((bin_size - 1) / 2.0) + lag
    valid_arr = arr[select_arr]

    # Create boolean to retrieve central timestamp for each averaged observation
    df['valid_tmstmp'] = df['mins'] == 0
      
    # Create boolean to retrieve appropriate subset of data
    if select_arr[0]<0:
        df['valid_cases'] = (df['mins'] >= valid_arr[0]) | \
                            (df['mins'] <= valid_arr[-1])
    else:
        df['valid_cases'] = (df['mins'] >= valid_arr[0]) & \
                            (df['mins'] <= valid_arr[-1])
    df = df[df['valid_cases']]
    
    # Create a sequentially numbered grouping variable, then group by it, and
    # set to nan any cases where there was less than a full window of data
    df['grp'] = np.arange(len(df)) / bin_size
    count_df = df.groupby('grp').count()[df_columns]
    mean_df = df.groupby('grp').mean()[df_columns]
    for var in df_columns:
        mean_df[var] = np.where(count_df[var] == bin_size, mean_df[var], np.nan)
    
    # Reindex using the center of the window, and drop working variables
    mean_df.index = pd.to_datetime(df[df['valid_tmstmp']].index)
    
    return mean_df

def get_CO2_data():

    """
    This script should be used to make any adjustments to CO2 data; it drops all
    other variables except the timestamp and the CO2 data because it assumes
    that these variables have been used to make corrections and are 
    redundant when the data is returned (as dictionary)
    """
 
    # Set locations   
    path = '/media/Data/Dropbox/Data/Logger data downloads 30minute/Whroo'
    files = ['Whroo_profile_IRGA_avg_1.dat' , 'Whroo_profile_IRGA_avg_2.dat',
             'Whroo_profile_IRGA_avg_3.dat', 'Whroo_profile_IRGA_avg.dat']

    # Set var names
    CO2_names = ['Cc_LI840_1m', 'Cc_LI840_2m', 'Cc_LI840_4m', 'Cc_LI840_8m', 
                 'Cc_LI840_16m', 'Cc_LI840_32m']

    # Set dates for correction
    last_1min_date = '2012-02-28 12:03:00'
    first_2min_date = '2012-02-28 12:10:00'
    baddata_dates = ['2013-08-24', '2013-10-29']
    badcoeff_dates = ['2012-06-28 11:00:00', '2012-10-17 12:50:00']
    instrument_dates = [['2011-12-02 12:00:00', '2012-06-28 10:58:00'],
                        ['2012-06-28 11:00:00', '2012-10-13 12:00:00'],
                        ['2012-10-13 12:02:00', '2013-08-23 23:58:00'],
               		 ['2013-10-29 12:00:00', '2014-06-02 23:58:00']]

    # Set some other stuff    
    coeff_correct = 2.5
    true_heights = [0.5, 2, 4, 8, 16, 36]
    CO2_range = [300, 600]

    # Import and clean up datasets
    df_dict = {}
    for f in files:
        df = pd.read_csv(os.path.join(path, f), skiprows = [0,2,3])
        df.index = pd.to_datetime(df['TIMESTAMP'])
        df.drop(['TIMESTAMP','RECORD'], axis = 1, inplace = True)
        df = df.astype(float)
        df_dict[f] = df
    
    # Change the earliest data to two minute frequency to match later data, then 
    # concatenate everything
    df = pd.concat([df_dict[files[0]].loc[:last_1min_date].resample('2T').mean(),
                    df_dict[files[0]].loc[first_2min_date:],
                    df_dict[files[1]],
                    df_dict[files[2]],
                    df_dict[files[3]]])
    df.drop_duplicates(inplace = True)
    df.sort_index(inplace = True)
    df = df.reindex(pd.date_range(df.index[0], df.index[-1], freq = '2T'))

    # Correct the data for 1) no data; 2) wrong instrument scaling coefficients; 
    # 3) range checks; 4) reversed label assignment of CO2
    for i in CO2_names:
        df.loc[baddata_dates[0]: baddata_dates[1], i] = np.nan
        df.loc[badcoeff_dates[0]: badcoeff_dates[1], i] = \
        df.loc[badcoeff_dates[0]: badcoeff_dates[1], i] * coeff_correct
        df[i][df[i] < CO2_range[0]] = np.nan
        df[i][df[i] > CO2_range[1]] = np.nan

    # 4 above
    true_heights.reverse()
    new_CO2_names = ['Cc_LI840_' + str(i) + 'm' for i in true_heights]
    reverse_dict = {CO2_names[i]: new_CO2_names[i] for i in range(len(CO2_names))}
    df = df.rename(columns = reverse_dict)

    return df[new_CO2_names]
    
def process_data(df, outer_dict):

    # Set parameters
    r = 8.3143 # Universal gas constant
    K_con = 273.15 # Kelvin conversion
    
    # Create dictionaries
    C_names_dict = outer_dict['C_names_dict']
    met_names_dict = outer_dict['met_names_dict']
    C_heights_dict = outer_dict['C_heights_dict']

    # Create a new dataframe containing the CO2 data, and rename to generic cols
    Cc_df = df[C_names_dict.values()].copy()
    Cc_df = Cc_df.rename(columns = dict(zip(C_names_dict.values(),
                                            C_names_dict.keys())))
    
    # Create a new dataframe containing the met data, and rename to generic cols
    met_df = df[met_names_dict.values()].copy()
    Cc_df = Cc_df.rename(columns = dict(zip(C_names_dict.values(),
                                            C_names_dict.keys())))
    
    # Reverse the C heights dictionary to make the height the key, 
    # then get numerical levels and sort
    C_heights_dict = dict(zip(C_heights_dict.values(),
                              C_heights_dict.keys()))
    levels = C_heights_dict.keys()
    levels.sort()
    
    # Calculate the layer averages (mean of upper and lower boundary);
    # lowest layer is just the value observed at the lowest boundary)
    layers_df = pd.DataFrame(index = df.index)
    layers_dict = {}
    for i in range(len(levels)):
        if i == 0:
            level_name = C_heights_dict[levels[i]]
            layer_thickness = levels[i]
            layer_name = 'Cc_layer_' + str(0) + '-' + str(levels[i]) + 'm'
            layers_df[layer_name] = Cc_df[level_name]
        else:
            upper_level_name = C_heights_dict[levels[i]]
            lower_level_name = C_heights_dict[levels[i - 1]]
            layer_thickness = levels[i] - levels[i - 1]
            layer_name = 'Cc_layer_' + str(levels[i - 1]) + '-' + str(levels[i]) + 'm'
            layers_df[layer_name] = (Cc_df[upper_level_name] + 
                                     Cc_df[lower_level_name]) / 2
        layers_dict[layer_name] = layer_thickness
    
    # Do the time differencing for each layer
    Cc_diff_df = layers_df - layers_df.shift()

    # Calculate the storage term: use ideal gas law to get molar density, then
    # scale to CO2 molar density (note that CO2 mixing ratio should be multiplied
    # by 10^-6, but since we want umol we would then multiply by 10^6 so both 
    # drop out), then scale over layer thickness, then scale to seconds;
    # Finally, sum over all thicknesses to get total
    storage_df = pd.DataFrame(index = df.index)
    molar_series = met_df['ps'] * 10**3 / (r * (K_con + met_df['Ta']))
    for layer_name in layers_df.columns:
        layer_thickness = layers_dict[layer_name]
        new_layer_name = 'umolC' + layer_name[2:]
        storage_df[new_layer_name] = (molar_series * Cc_diff_df[layer_name] *
                                      layer_thickness / 1800)
    storage_df['C_stor_total'] = storage_df[storage_df.columns].sum(axis = 1)
    
    return storage_df