# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:21:38 2015

@author: imchugh
"""

import datetime as dt
import numpy as np
import pdb

def get_timestep(datetime_array):
    """
    Checks timedelta for consistency and finds measurement interval
    Takes datetime array as argument
    Returns measurement interval in minutes
    """
    check_timedelta = datetime_array[1: ] - datetime_array[: -1]
    if not all(check_timedelta[0] == rest for rest in check_timedelta):
        raise Exception('Time series is not continuous!')
    else:
        return check_timedelta[0].seconds / 60.0

def get_moving_window(data_dict,
                      datetime_varname, 
                      window, 
                      step, 
                      return_indices_only = False, 
                      retro_stamp = True):
    """
    Finds available windows given window and step
    Pass args: 1) data dictionary containing a numpy array of python datetimes, 
               2) the variable name for the datetime array (str)
               3) window (int)
               4) step (int)
    Optional kwargs: 1) return_indices_only:
                         - if true will return a dictionary with key / value 
                           pairs of:
                               - date of centre of window: start and end indices 
                                                          (list of ints)'
                         - if false will return a dicitonary with key / value 
                           pairs of:
                               - date of centre of window: complete time series
                                                          within the window
                                                          boundaries
                     2) retro_stamp:
                         - if true then selects dates on the basis that the 
                           datetime format in the file is retrospective (i.e.
                           that the timestmap represents the END of the period
                           it references e.g. 28/07/2015 00:00 for a 
                           half-hourly frequency time series represents the 
                           period 27/07/15 23:30 - 28/07/15 00:00)
    """

    # Get datetime array
    datetime_array = data_dict[datetime_varname]

    # Check measurement interval    
    meas_int = get_timestep(datetime_array)
    
    if retro_stamp:
        shift_by = meas_int
    else:
        shift_by = 0

    # Find part days at beginning
    start_date = (dt.datetime.combine(datetime_array[0].date(), 
                                      dt.datetime.min.time()) +
                  dt.timedelta(minutes = shift_by))
    if start_date > datetime_array[0]:
        start_index = np.where(datetime_array == start_date)[0].item()
    elif start_date < datetime_array[0]:
        start_date = start_date + dt.timedelta(1)
        start_index = np.where(datetime_array == start_date)[0].item()
    else:
        start_index = 0              

    # Find part days at end
    end_date = (dt.datetime.combine(datetime_array[-1].date(), 
                                    dt.datetime.min.time()) +
                                    dt.timedelta(1) - 
                                    dt.timedelta(minutes = meas_int) +
                                    dt.timedelta(minutes = shift_by))
    if end_date > datetime_array[-1]:
        end_date = end_date - dt.timedelta(1)
        end_index = np.where(datetime_array == end_date)[0].item()               
    else:
        end_index = len(datetime_array)

    # Slice a new working array
    work_array = datetime_array[start_index: end_index]

    # Generate dates representing begin, centre and end of window
    num_days = (work_array[-1].date() - 
                work_array[0].date()).days + 1 - window
    num_days_range = range(0, num_days + 1, step)
    centre_datetime_array = np.array([(work_array[0] + 
                                       dt.timedelta(i + window / 2.0))
                                      for i in num_days_range])
    begin_datetime_array = np.array(centre_datetime_array - 
                                    dt.timedelta(window / 2.0))
    end_datetime_array = np.array(centre_datetime_array + 
                                  dt.timedelta(window / 2.0) - 
                                  dt.timedelta(minutes = meas_int))

    # Create dictionary with date as key and indices as values
    step_dates_index_dict = {}
    centre_date_array = np.array([date_time.date() 
                                  for date_time in centre_datetime_array])
    for i, date in enumerate(centre_date_array):
        begin_ind = np.where(datetime_array == begin_datetime_array[i])[0].item()
        end_ind = np.where(datetime_array == end_datetime_array[i])[0].item()
        step_dates_index_dict[date] = [begin_ind, end_ind]

    if return_indices_only:
        return step_dates_index_dict
    else:
        return segment_data(data_dict, step_dates_index_dict)

def get_year_window(data_dict,
                    datetime_varname,
                    return_indices_only = False,
                    retro_stamp = True):    
    """
    Finds the array location indices for the years;
    Takes datetime array as arg
    Returns dictionary containing year as key and indices as value
    Note: retro_stamp indicates that the timestamp applies to the data 
          retrospectively i.e. the 00:00 timestamp indicates a data collection
          period of 23:30-00:00; so for example the 2012 data year is correctly
          represented by 2012-01-01 00:30 to 2013-01-01 00:00; if this is set 
          to false, it will interpret the timestamps literally, so the 2012 
          data year will be represented by 2012-01-01 00:00 to 2012-12-31 23:30. 
          This is only correct if timestamps are not retrospective!
    """    

    # Get datetime array
    datetime_array = data_dict[datetime_varname]
 
    # Check measurement interval    
    meas_int = get_timestep(datetime_array)
    
    if retro_stamp:
        shift_by = meas_int
    else:
        shift_by = 0

    datetime_array = datetime_array - dt.timedelta(minutes = shift_by)
    years_index_dict = {}
    year_array = np.array([i.year for i in datetime_array])
    year_list = list(set(year_array))
    for yr in year_list:
        index = np.where(year_array == yr)[0]
        years_index_dict[yr] = [index[0], index[-1]]

    if return_indices_only:
        return years_index_dict
    else:
        return segment_data(data_dict, years_index_dict)
        
    return years_index_dict

def get_day_indices(datetime_array, retro_stamp = True):
    """
    Finds the array location indices for the days;
    Takes datetime array as arg
    Returns dictionary containing day date as key and indices as value
    Note: retro_stamp indicates that the timestamp applies to the data 
          retrospectively i.e. the 00:00 timestamp indicates a data collection
          period of 23:30-00:00; so for example the last day of 2012 is 
          correctly represented by 2012-12-31 00:30 to 2013-01-01 00:00; 
          if this is set to false, it will interpret the timestamps literally, 
          so the above day will be represented by 2012-01-01 00:00 to 
          2012-12-31 23:30. This is only correct if timestamps are not 
          retrospective!
    """    
    
    # Check measurement interval    
    meas_int = get_timestep(datetime_array)    

    if retro_stamp:
        shift_by = meas_int
    else:
        shift_by = 0

    datetime_array = datetime_array - dt.timedelta(minutes = shift_by)
    days_index_dict = {}
    date_array = np.array([i.date() for i in datetime_array])
    date_list = list(set(date_array))
    for date in date_list:
        index = np.where(date_array == date)[0]
        days_index_dict[date] = [index[0], index[-1]]                               
            
    return days_index_dict
    
def get_unique_dates(datetime_array):
    
    date_array = np.array(list(set([i.date() for i in datetime_array])))
    date_array.sort()
    return {date: i for i, date in enumerate(date_array)}
    
def segment_data(data_dict, indices_dict):

    d = {}    
    for key in indices_dict.keys():
        start = indices_dict[key][0]
        end = indices_dict[key][1]
        this_dict = {var: data_dict[var][start: end + 1] 
                     for var in data_dict.keys()}
        d[key] = this_dict

    return d    
    
def get_DOY_first_day_of_month(year):
    
    return [int(dt.datetime.strftime(dt.date(2012,m,1), '%j')) - 1 
            for m in range(1, 13)]
                
def subset_datayear_from_arraydict(data_dict, date_time_var, year = None):
    """
    Pass: 1) data_dict - a dictionary containing arrays, one of which must be a 
             python datetime;
          2) date_time_var - namestring of datetime variable
          3) year to be returned as optional kwarg
    Returns: if year is specified, return the same dictionary structure with 
             only data for that year; if no year is specified, return a 
             dictionary with each data year contained as the value with the 
             year as the key
    """    
    years_array = np.array([date_.year for date_ in data_dict[date_time_var]])
    if not year:
        year_list = set(list(years_array))    
    else:
        if not isinstance(year, list): year_list = [year]
    
    new_dict = {}
    for yr in year_list:
        year_index = years_array == yr            
        new_dict[yr] = {var: data_dict[var][year_index] 
                        for var in data_dict.keys()}
    return new_dict                