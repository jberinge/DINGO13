# -*- coding: utf-8 -*-
"""
Created on Thu May 26 16:09:18 2016

@author: imchugh
"""

import os
import numpy as np

import DataIO as io

def data_check(check_var, insol_var):
    
    """
    This module quantifies data availability in OzFluxQC .nc or DINGO .df files
    
    Keyword arguments are:
        1) check_var (str): the name of the variable to be checked in the relevant file
        2) insol_var (str): the name of the insolation variable for separating
           day and night data
           
    Returns a dictionary separated into years and night and day, with counts
    for total and available observations, and percentage available
    """    
    
    file_in = io.file_select_dialog()

    var_list = [check_var, insol_var]
    
    ext = os.path.splitext(file_in)[1]
    
    if ext == 'nc':
        data_dict = io.OzFluxQCnc_to_data_structure(file_in, var_list = var_list)
    elif ext == '.df':
        data_dict = io.DINGO_df_to_data_structure(file_in, var_list = var_list)
        
    year_series_list = [i.year for i in data_dict['date_time']]
    year_series_array = np.array(year_series_list)
    years_list = list(set(year_series_list))
    
    results_dict = {}
    for yr in years_list:
        results_dict[yr] = {}
        index = np.where(year_series_array == yr)
        this_dict = {var: data_dict[var][index] for var in var_list}        
        for cond in ['day', 'night']:
            results_dict[yr][cond] = {}
            if cond == 'day': 
                data_array = this_dict[check_var][this_dict[insol_var] > 5]
            else:
                data_array = this_dict[check_var][this_dict[insol_var] < 5]
            tot_int = len(data_array)
            avail_int = len(data_array[~np.isnan(data_array)])
            avail_pct = np.round(avail_int / float(tot_int) * 100, 1)
            results_dict[yr][cond]['potential_records'] = tot_int
            results_dict[yr][cond]['observed_records'] = avail_int
            results_dict[yr][cond]['percent_available'] = avail_pct

    return results_dict